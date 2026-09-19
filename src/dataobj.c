#include "wrassp.h"
#include <math.h>               /* ceil, floor */
#include <dataobj.h>
#include <asspfio.h>
#include <asspmess.h>
#include <headers.h>            /* KDTAB */
#include "ssff_convert.hpp"     /* SSFF record buffer -> R matrices */

#ifndef _WIN32
#include <sys/mman.h>           /* mmap: read SSFF windows without an extra copy */
#include <sys/stat.h>
#include <unistd.h>
#endif


/*
 * Mapping costs a handful of syscalls, which only pays off once the window is
 * large; smaller windows are read into the DOBJ buffer as before.
 */
#define SUPERASSP_MMAP_MIN_BYTES ((size_t) 256 * 1024)

/*
 * SSFF record windows.
 *
 * Maps the requested records of an open data object and returns 1 on success,
 * with *records pointing at the first record, mapBase and mapLen describing the
 * mapping (release with munmap()) and *numRecs possibly reduced when the file
 * holds fewer records than the header claims. Returns 0 when mapping is not
 * possible (no mmap on this platform, non-regular file, ...); the caller then
 * falls back to reading the window into the DOBJ buffer.
 */

#ifdef _WIN32
static int
ssffMapWindow(DOBJ * dop, long startRec, long *numRecs, const void **records,
              void **mapBase, size_t *mapLen)
{
    (void) dop; (void) startRec; (void) numRecs;
    (void) records; (void) mapBase; (void) mapLen;
    return 0;
}
#else
static int
ssffMapWindow(DOBJ * dop, long startRec, long *numRecs, const void **records,
              void **mapBase, size_t *mapLen)
{
    struct stat     st;
    long            pageSize, dataOff, mapOff, delta, avail;
    size_t          len;
    char           *p;
    int             fd;

    if (dop == NULL || dop->fp == NULL || dop->recordSize < 1 ||
        dop->headerSize <= 0 || *numRecs < 1)
        return 0;
    fd = fileno(dop->fp);
    if (fd < 0 || fstat(fd, &st) != 0 || !S_ISREG(st.st_mode))
        return 0;
    dataOff = (long) dop->headerSize +
        (startRec - dop->startRecord) * (long) dop->recordSize;
    if (dataOff < 0 || (long) st.st_size <= dataOff)
        return 0;
    avail = ((long) st.st_size - dataOff) / (long) dop->recordSize;
    if (avail < *numRecs)
        *numRecs = avail;           /* truncated file: return what is there */
    if (*numRecs < 1)
        return 0;
    pageSize = sysconf(_SC_PAGESIZE);
    if (pageSize < 1)
        return 0;
    mapOff = (dataOff / pageSize) * pageSize;
    delta = dataOff - mapOff;
    len = (size_t) delta + (size_t) (*numRecs) * dop->recordSize;
    p = (char *) mmap(NULL, len, PROT_READ, MAP_PRIVATE, fd, (off_t) mapOff);
    if (p == MAP_FAILED)
        return 0;
#ifdef MADV_WILLNEED
    madvise((void *) p, len, MADV_WILLNEED);   /* read-ahead: fewer per-page faults */
#endif
    *mapBase = (void *) p;
    *mapLen = len;
    *records = (const void *) (p + delta);
    return 1;
}
#endif

/*
 * Track selection. 'sel' is a character vector of descriptor identifiers or
 * NULL, which selects everything.
 */
static int
trackWanted(DDESC * desc, SEXP sel, int nsel)
{
    int             k;

    if (nsel == 0)
        return 1;
    for (k = 0; k < nsel; k++) {
        if (strcmp(desc->ident, CHAR(STRING_ELT(sel, k))) == 0)
            return 1;
    }
    return 0;
}

/*
 * Report requested tracks that the file does not contain, listing the tracks
 * that it does contain.
 */
static void
validateTrackSelection(DOBJ * dop, SEXP sel)
{
    char            avail[ONEkBYTE];
    const char     *name;
    DDESC          *desc;
    int             k, found;
    size_t          len = 0, nlen;

    avail[0] = '\0';
    for (desc = &(dop->ddl); desc != NULL; desc = desc->next) {
        nlen = strlen(desc->ident);
        if (len + nlen + 3 >= sizeof(avail))
            break;
        if (len > 0) {
            strcpy(avail + len, ", ");
            len += 2;
        }
        strcpy(avail + len, desc->ident);
        len += nlen;
    }
    for (k = 0; k < LENGTH(sel); k++) {
        name = CHAR(STRING_ELT(sel, k));
        found = 0;
        for (desc = &(dop->ddl); desc != NULL; desc = desc->next) {
            if (strcmp(desc->ident, name) == 0) {
                found = 1;
                break;
            }
        }
        if (!found) {
            error("Track '%s' not found in %s. Available tracks: %s",
                  name, dop->filePath, avail);
        }
    }
}

/*
 * This was the original reading function that did not allow for
 * preselecting time. Should be save to remove now. 
 */
SEXP
getDObj(SEXP fname)
{
    SEXP            res;
    DOBJ           *data = NULL;
    long            numRecs;
    /*
     * read the data
     */
    data =
        asspFOpen(strdup(CHAR(STRING_ELT(fname, 0))), AFO_READ,
                  (DOBJ *) NULL);
    if (data == NULL)
        error("%s", getAsspMsg(asspMsgNum));
    /*
     * error(CHAR(STRING_ELT(fname,0)));
     */
    allocDataBuf(data, data->numRecords);
    data->bufStartRec = data->startRecord;
    if ((numRecs = asspFFill(data)) < 0)
        error("%s", getAsspMsg(asspMsgNum));
    asspFClose(data, AFC_KEEP);
    res = PROTECT(dobj2AsspDataObj(data));
    asspFClose(data, AFC_FREE);
    UNPROTECT(1);
    return res;
}


/*
 * This function loads a DOBJ from a file and return its contents as a
 * SEXP. Arguments include the name of the input file, start and end point 
 * for reading and whether these points are sample values (and not times) 
 */
SEXP
getDObj2(SEXP args)
{

    SEXP            el,
                    ans;
    SEXP            trackSel = R_NilValue;
    DOBJ           *data = NULL;
    long            numRecs;
    int             i;
    char           *fName = NULL;
    const char     *name;
    double          begin = 0,
        end = 0;
    int             isSample = 0;
    int             zeroToNa = 0;
    int             numThreads = 1;

    /*
     * parse args
     */
    args = CDR(args);           /* skip name of function */
    el = CAR(args);
    fName = strdup(CHAR(STRING_ELT(el, 0)));

    args = CDR(args);
    for (i = 0; args != R_NilValue; i++, args = CDR(args)) {
        name = isNull(TAG(args)) ? "" : CHAR(PRINTNAME(TAG(args)));
        el = CAR(args);
        if (strcmp(name, "begin") == 0) {
            begin = REAL(el)[0];
            if (begin < 0)
                begin = 0;
        } else if (strcmp(name, "end") == 0) {
            end = REAL(el)[0];
            if (end < 0)
                end = 0;
        } else if (strcmp(name, "samples") == 0) {
            isSample = INTEGER(el)[0];
        } else if (strcmp(name, "zero_to_na") == 0) {
            zeroToNa = asLogical(el);
            if (zeroToNa == NA_LOGICAL)
                zeroToNa = 0;
        } else if (strcmp(name, "tracks") == 0) {
            if (!isNull(el) && TYPEOF(el) != STRSXP)
                error("'tracks' must be a character vector or NULL.");
            trackSel = el;
        } else if (strcmp(name, "threads") == 0) {
            numThreads = INTEGER(el)[0];
            if (numThreads < 1)
                numThreads = 1;
        } else {
            error("Bad option '%s'.", name);
        }
    }

    if (end < begin && end > 0)
        error("End before begin. That's not clever, dude!");

    /*
     * open the file
     */
    data = asspFOpen(fName, AFO_READ, (DOBJ *) NULL);
    if (data == NULL)
        error("%s (%s)", getAsspMsg(asspMsgNum), fName);
    /*
     * figure out timing
     */
    if (isSample) {
        if (end == 0)
            end = data->startRecord + data->numRecords - 1;
    } else {
        begin = ceil(begin * data->dataRate) + data->startRecord;
        if (end == 0)
            end = data->startRecord + data->numRecords - 1;
        else
            end = floor(end * data->dataRate) + data->startRecord;
    }
    if (end > (data->startRecord + data->numRecords))
        end = data->startRecord + data->numRecords - 1;
    if (begin > (data->startRecord + data->numRecords)) {
        asspFClose(data, AFC_FREE);
        error("Begin after end of data. That's not clever, dude!");
    }

    numRecs = (long) (end - begin) + 1;
    /*
     * read the data.
     *
     * Preferred path: map the requested window and convert it straight into the
     * R matrices. Fallback (no mmap, e.g. Windows, or a non-regular file): read
     * into the DOBJ buffer with asspFFill() as before; that path swaps the
     * buffer to host byte order, so conversion then runs without swapping.
     */
    {
        const void     *records = NULL;
        void           *mapBase = NULL;
        size_t          mapLen = 0;
        int             swapped = 0;

        if (trackSel != R_NilValue)
            validateTrackSelection(data, trackSel);

        if ((size_t) numRecs * data->recordSize >= SUPERASSP_MMAP_MIN_BYTES &&
            ssffMapWindow(data, (long) begin, &numRecs, &records, &mapBase, &mapLen)) {
            ENDIAN          sysEndian = { MSB };
            swapped = DIFFENDIAN(data->fileEndian, sysEndian) ? 1 : 0;
            data->bufStartRec = (long) begin;
            data->bufNumRecs = numRecs;
        } else {
            allocDataBuf(data, numRecs);
            data->bufStartRec = (long) begin;
            if ((numRecs = asspFFill(data)) < 0) {
                asspFClose(data, AFC_FREE);
                error("%s", getAsspMsg(asspMsgNum));
            }
            records = data->dataBuffer;
            swapped = 0;
        }
        asspFClose(data, AFC_KEEP);
        ans = PROTECT(dobj2AsspDataObjEx(data, records, zeroToNa, trackSel,
                                         numThreads, swapped));
        if (mapBase != NULL)
            munmap(mapBase, mapLen);
        asspFClose(data, AFC_FREE);
        UNPROTECT(1);
    }
    return ans;
}

/*
 * Originally, we retained the DOBJ and stored a pointer to it in the
 * SEXP. For that reason, garbage collection was an issue and this
 * function was used to clean up the data object when the SEXP was
 * deleted. No longer needed, should be save to remove. 
 */
static void
DObjFinalizer(SEXP dPtr)
{
    DOBJ           *data = R_ExternalPtrAddr(dPtr);
    asspFClose(data, AFC_FREE);
    R_ClearExternalPtr(dPtr);   /* not really needed */
}

/*
 * This function turns a DOBJ and places the contents in a SEXP of class
 * AsspDataObject. (Hopefully) all information is retained in order to
 * safely rewrite without data loss. 
 */
/*
 * This function turns a DOBJ and places the contents in a SEXP of class
 * AsspDataObject. (Hopefully) all information is retained in order to
 * safely rewrite without data loss.
 *
 * Kept as the plain "all tracks, verbatim values" entry point for callers that
 * hand over an in-memory DOBJ (see performAssp.c).
 */
SEXP
dobj2AsspDataObj(DOBJ * data)
{
    return dobj2AsspDataObjEx(data, data->dataBuffer, 0, R_NilValue, 1, 0);
}

/*
 * As dobj2AsspDataObj(), with the options of the file reader:
 *   records    - interleaved record buffer (e.g. a mapped SSFF window), or NULL
 *                to use the DOBJ's own buffer
 *   zeroToNa   - map values that are exactly 0 to NA for descriptors that are
 *                not sampled audio (SSFF has no NA; 0 is its substitute)
 *   trackSel   - character vector of descriptor identifiers, NULL = all
 *   numThreads - OpenMP threads for the block-wise fill (1 = serial)
 *   swapped    - 1 when 'records' is in file byte order, 0 when it already is
 *                in host byte order
 */
SEXP
dobj2AsspDataObjEx(DOBJ * data, const void *records, int zeroToNa, SEXP trackSel,
                   int numThreads, int swapped)
{
    SEXP            ans,        /* dPtr, */
                    class,
                    rate,
                    tracks,
                    startTime,
                    origRate,
                    filePath,
                    startRec,
                    endRec,
                    trackFormats,
                    finfo,
                    genericVars;
    DDESC          *desc = NULL;
    ssff_track_t   *trk = NULL;
    int             i,
                    n,
                    nsel;

    if (records == NULL)
        records = data->dataBuffer;
    if (records == NULL)
        error("No data buffer to read from.");

    nsel = isNull(trackSel) ? 0 : LENGTH(trackSel);

    /*
     * count the tracks to return
     */
    for (n = 0, desc = &(data->ddl); desc != NULL; desc = desc->next) {
        if (trackWanted(desc, trackSel, nsel))
            n++;
    }

    /*
     * create result, a list with a matrix for each track
     */
    PROTECT(ans = allocVector(VECSXP, n));
    /*
     * create list of tracks and formats
     */
    PROTECT(tracks = allocVector(STRSXP, n));
    PROTECT(trackFormats = allocVector(STRSXP, n));
    trk = (ssff_track_t *) R_alloc((size_t) (n > 0 ? n : 1),
                                   sizeof(ssff_track_t));

    for (i = 0, desc = &(data->ddl); desc != NULL; desc = desc->next) {
        SEXP            mat;
        int             isFloat;

        if (!trackWanted(desc, trackSel, nsel))
            continue;
        SET_STRING_ELT(tracks, i, mkChar(desc->ident));
        SET_STRING_ELT(trackFormats, i,
                       mkChar(asspDF2ssffString(desc->format)));
        /*
         * fill tracks with data
         */
        isFloat = ssff_format_is_float(desc->format);
        mat = allocMatrix(isFloat ? REALSXP : INTSXP,
                          (int) data->bufNumRecs, (int) desc->numFields);
        SET_VECTOR_ELT(ans, i, mat);
        trk[i].offset = desc->offset;
        trk[i].numFields = desc->numFields;
        trk[i].format = desc->format;
        trk[i].destType = isFloat ? REALSXP : INTSXP;
        trk[i].dest = isFloat ? (void *) REAL(mat) : (void *) INTEGER(mat);
        trk[i].zeroToNa = (zeroToNa && desc->type != DT_SMP) ? 1 : 0;
        i++;
    }

    /*
     * a single pass over the record buffer fills every track
     */
    if (n > 0) {
        if (ssff_convert_records(records, data->bufNumRecs, data->recordSize,
                                 n, trk, swapped, numThreads) != 0)
            error("Unsupported data format.");
    }
    /*
     * set the names
     */
    setAttrib(ans, R_NamesSymbol, tracks);
    setAttrib(ans, install("trackFormats"), trackFormats);

    /*
     * PROTECT (dPtr = R_MakeExternalPtr (data, install ("DOBJ"), 
     * install ("something"))); 
     * R_RegisterCFinalizerEx (dPtr, DObjFinalizer, TRUE); 
     * setAttrib (ans, install ("data pointer"), dPtr); 
     */

    PROTECT(rate = allocVector(REALSXP, 1));
    REAL(rate)[0] = data->dataRate;
    setAttrib(ans, install("sampleRate"), rate);
    if (data->filePath == NULL || strlen(data->filePath) == 0){
        // Rprintf("at non caps protect call\n");
        PROTECT(filePath = R_NilValue);
    } else {
        PROTECT(filePath = allocVector(STRSXP, 1));
        SET_STRING_ELT(filePath, 0, mkCharCE(data->filePath, CE_UTF8));
    }
    setAttrib(ans, install("filePath"), filePath);
    PROTECT(origRate = allocVector(REALSXP, 1));
    if (data->fileFormat == FF_SSFF) {
        REAL(origRate)[0] = data->sampFreq;
    } else {
        REAL(origRate)[0] = 0;
    }
    setAttrib(ans, install("origFreq"), origRate);
    PROTECT(startTime = allocVector(REALSXP, 1));
    /*
     * REAL (startTime)[0] = data->Start_Time + 
     * (data->bufStartRec / data->dataRate); 
     */
    REAL(startTime)[0] = data->Start_Time;
    setAttrib(ans, install("startTime"), startTime);

    PROTECT(startRec = allocVector(INTSXP, 1));
    INTEGER(startRec)[0] = (int) (data->bufStartRec + 1);
    setAttrib(ans, install("startRecord"), startRec);
    PROTECT(endRec = allocVector(INTSXP, 1));
    INTEGER(endRec)[0] = (int) (data->bufStartRec + data->bufNumRecs);
    setAttrib(ans, install("endRecord"), endRec);

    PROTECT(class = allocVector(STRSXP, 1));
    SET_STRING_ELT(class, 0, mkChar(WRASSP_CLASS));
    classgets(ans, class);

    PROTECT(finfo = allocVector(INTSXP, 2));
    INTEGER(finfo)[0] = (int) data->fileFormat;
    INTEGER(finfo)[1] = (int) data->fileData;
    setAttrib(ans, install("fileInfo"), finfo);


    PROTECT(genericVars = getGenericVars(data));
    setAttrib(ans, install("genericVars"), genericVars);

    UNPROTECT(12);
    return ans;

}

/*
 * This function parses generic variables from a DOBJ (SFFF only) and
 * returns them in a useful format for R 
 */
SEXP
getGenericVars(DOBJ * dop)
{
    SEXP            ans,
                    var,
                    names,
                    varNames,
                    value;
    TSSFF_Generic  *genVar;
    SSFFST         *ssff_types;
    int             i;
    PROTECT(names = allocVector(STRSXP, 2));
    SET_STRING_ELT(names, 1, mkChar("Type"));
    SET_STRING_ELT(names, 0, mkChar("Value"));
    /*
     * count generic variables 
     */
    for (i = 0, genVar = &(dop->meta); genVar != NULL;
         genVar = genVar->next, i++) {
    }
    if (i == 0) {
        UNPROTECT(1);
        return (R_NilValue);
    }

    PROTECT(ans = allocVector(VECSXP, i));
    PROTECT(varNames = allocVector(STRSXP, i));
    for (i = 0, genVar = &(dop->meta); genVar != NULL;
         i++, genVar = genVar->next) {
        if (genVar->ident == NULL) {
            UNPROTECT(3);
            return (R_NilValue);
        }
        PROTECT(var = allocVector(VECSXP, 2));
        for (ssff_types = SSFF_TYPES; ssff_types->type != SSFF_UNDEF;
             ssff_types++) {
            if (ssff_types->type == genVar->type)
                break;
        }
        if (ssff_types->type == SSFF_UNDEF)
            error("Invalid type for SSFF generic variable.");
        PROTECT(value = allocVector(STRSXP, 1));
        SET_STRING_ELT(value, 0, mkChar(ssff_types->ident));
        SET_VECTOR_ELT(var, 1, value);
        switch (genVar->type) {
        case SSFF_CHAR:
        case SSFF_BYTE:
            PROTECT(value = allocVector(STRSXP, 1));
            SET_STRING_ELT(value, 0, mkChar(genVar->data));
            SET_VECTOR_ELT(var, 0, value);
            UNPROTECT(1);
            break;
        case SSFF_SHORT:
        case SSFF_LONG:
        case SSFF_FLOAT:
        case SSFF_DOUBLE:
            PROTECT(value = allocVector(REALSXP, 1));
            REAL(value)[0] = strtod(genVar->data, NULL);
            SET_VECTOR_ELT(var, 0, value);
            UNPROTECT(1);
        case SSFF_UNDEF:
            break;
        }
        setAttrib(var, R_NamesSymbol, names);
        SET_VECTOR_ELT(ans, i, var);
        SET_STRING_ELT(varNames, i, mkChar(genVar->ident));
        UNPROTECT(2);
    }
    setAttrib(ans, R_NamesSymbol, varNames);
    UNPROTECT(3);
    return (ans);
}

/*
 * This function generates a string vector of track names from the data
 * descriptors in a DOBJ and returns it. 
 */
SEXP
getDObjTracks(SEXP dobj)
{
    SEXP            ans,
                    ptr;
    ptr = getAttrib(dobj, install("data pointer"));
    DOBJ           *data = R_ExternalPtrAddr(ptr);
    DDESC          *desc;
    int             i = 0,
        n = 0;
    /*
     * count tracks
     */
    for (desc = &(data->ddl); desc != NULL; desc = desc->next) {
        n++;
    }
    /*
     * Rprintf("Number of descs = %i.", n);
     */
    /*
     * create result
     */
    PROTECT(ans = allocVector(STRSXP, n));
    for (desc = &(data->ddl); desc != NULL; desc = desc->next) {
        SET_STRING_ELT(ans, i, mkChar(desc->ident));
        i++;
    }
    /*
     * for (; i<n; i++) 
     */
    /*
     * SET_STRING_ELT(ans, i, mkChar("")); 
     */
    UNPROTECT(1);
    return (ans);
}

/*
 * switch trough assp data formats and return corresponding string 
 */
char           *
asspDF2ssffString(int df)
{
    switch ((dform_e) df) {
    case DF_BIT:
        return "BIT";
        break;
    case DF_STR:
        return "STR";
        break;
    case DF_CHAR:
        return "CHAR";
        break;
    case DF_UINT8:
        return "UINT8";
        break;
    case DF_INT8:
        return "INT8";
        break;
    case DF_UINT16:
        return "UINT16";
        break;
    case DF_INT16:
        return "INT16";
        break;
    case DF_UINT24:
        return "UINT24";
        break;
    case DF_INT24:
        return "INT24";
        break;
    case DF_UINT32:
        return "UINT32";
        break;
    case DF_INT32:
        return "INT32";
        break;
    case DF_UINT64:
        return "UINT64";
        break;
    case DF_INT64:
        return "INT64";
        break;
    case DF_REAL32:
        return "REAL32";
        break;
    case DF_REAL64:
        return "REAL64";
        break;
    default:
        return NULL;
    }
}

/*
 * This function is the inverse of dobj2AsspDataObj. It takes a SEXP of
 * class AsspDataObj and turns it into a DOBJ 
 */
DOBJ           *
sexp2dobj(SEXP rdobj)
{
    DOBJ           *dop = NULL;
    DDESC          *desc = NULL;
    int             FIRST = 1,
        i = 0,
        myBool = 0;
    size_t          numFields = -1;
    SEXP            attr,
                    tracks,
                    formats,
                    track,
                    var,
                    varNames;
    SSFFST         *ssff_types;
    TSSFF_Generic  *genVar;
    KDTAB          *entry;
    char           *format;

    /*
     * check for right class
     */
    attr = getAttrib(rdobj, R_ClassSymbol);
    for (i = 0; i < LENGTH(attr); i++) {
        if (strcmp(CHAR(STRING_ELT(attr, i)), WRASSP_CLASS) == 0) {
            myBool = 1;
            break;
        }
    }
    if (!myBool) {              /* classname does not match */
        error("Argument must be of class %s", WRASSP_CLASS);
    }
    /*
     * create DObj
     */
    dop = allocDObj();
    desc = &(dop->ddl);
    if (dop == NULL) {
        error("%s", getAsspMsg(asspMsgNum));
    }

    /*
     * Set openMode to AFO_READ to indicate this DOBJ is ready for processing
     * This is required for DSP functions to accept the DOBJ
     */
    dop->openMode = AFO_READ;

    /*
     * assign attributes
     */
    attr = getAttrib(rdobj, install("filePath"));
    if (!isNull(attr) && LENGTH(attr) > 0) {
        dop->filePath = strdup(CHAR(STRING_ELT(attr, 0)));
    }

    attr = getAttrib(rdobj, install("sampleRate"));
    if (isNull(attr)) {
        freeDObj(dop);
        error("Invalid argument: no 'sampleRate' attribute.");
    }
    dop->dataRate = REAL(attr)[0];

    attr = getAttrib(rdobj, install("origFreq"));
    if (!isNull(attr))
        dop->sampFreq = REAL(attr)[0];

    /* If sampFreq is 0 (e.g., for WAVE files), use dataRate as sampFreq */
    if (dop->sampFreq == 0.0)
        dop->sampFreq = dop->dataRate;

    attr = getAttrib(rdobj, install("startTime"));
    if (!isNull(attr))
        dop->Start_Time = REAL(attr)[0];

    attr = getAttrib(rdobj, install("startRecord"));
    if (!isNull(attr))
        /* R uses 1-indexed records, C uses 0-indexed */
        dop->startRecord = INTEGER(attr)[0] - 1;

    attr = getAttrib(rdobj, install("fileInfo"));
    if (LENGTH(attr) != 2) {
        dop->fileFormat = FF_SSFF;
        dop->fileData = FDF_BIN;
        warning("Incomplete 'fileInfo' attribute. Writing to binary"
                "SSFF format (dafault).");
    }
    dop->fileFormat = (fform_e) INTEGER(attr)[0];
    dop->fileData = (fdata_e) INTEGER(attr)[1];

    /*
     * Set fileEndian based on file format
     * For memory-based DOBJs, data is in native endianness
     * Set it to match the file format's expected endianness
     */
    switch(dop->fileFormat) {
    case FF_WAVE:
    case FF_WAVE_X:
    case FF_CSL:
    case FF_CSRE:
        /* These formats are always little-endian (MSB last) */
        SETMSBLAST(dop->fileEndian);
        break;
    case FF_AIFF:
    case FF_AIFC:
    case FF_SND:
    case FF_KTH:
        /* These formats are always big-endian (MSB first) */
        SETMSBFIRST(dop->fileEndian);
        break;
    default:
        /* For other formats including FF_SSFF and FF_RAW, use native */
        /* Native endianness on this system */
        {
            ENDIAN sysEndian={MSB};
            CPYENDIAN(dop->fileEndian, sysEndian);
        }
        break;
    }

    /*
     * Install generic variables
     */
    attr = getAttrib(rdobj, install("genericVars"));
    tracks = getAttrib(attr, R_NamesSymbol);
    if (!isNull(attr)) {
        for (i = 0; i < LENGTH(attr); i++) {
            var = VECTOR_ELT(attr, i);
            /*
             * determine ssff type
             */
            for (ssff_types = SSFF_TYPES; ssff_types->type != SSFF_UNDEF;
                 ssff_types++) {
                format = strdup(CHAR(STRING_ELT(VECTOR_ELT(var, 1), 0)));
                if (strncmp
                    (format, ssff_types->ident,
                     strlen(ssff_types->ident)) == 0)
                    break;
            }
            if (ssff_types->type == SSFF_UNDEF)
                error("Invalid type for SSFF generic variable.");
            if (FIRST) {
                genVar = &(dop->meta);
            } else {
                genVar = addTSSFF_Generic(dop);
                if (genVar == NULL)
                    error("Unable to add Generic Variable (%s).",
                          getAsspMsg(asspMsgNum));
            }
            FIRST = 0;

            genVar->type = ssff_types->type;
            genVar->ident = strdup(CHAR(STRING_ELT(tracks, i)));
            switch (genVar->type) {
            case SSFF_CHAR:
            case SSFF_BYTE:
                genVar->data =
                    strdup(CHAR(STRING_ELT(VECTOR_ELT(var, 0), 0)));
                break;
            case SSFF_SHORT:
            case SSFF_LONG:
            case SSFF_FLOAT:
            case SSFF_DOUBLE:
                // While replacing sprintf with snprintf, I am wondering how much
                // reserved memory `format` is actually pointing at; because that
                // is what I want to pass to snprintf as the second variable. The
                // memory was reserved by strdup() above. I will assume that the
                // strdup() function only reserved as many bytes as it needs. This
                // would mean I can safely cut off snprintf after strlen(format)+1
                // bytes, which is what I am now doing.
                //
                // The only way this change might introduce a regression is if:
                // A. strdup() reserved more memory than is necessary AND
                // B. The value that gets written into *format needs that excess memory.
                //
                snprintf(format, strlen(format) + 1, "%f", REAL(VECTOR_ELT(var, 0))[0]);
                genVar->data = strdup(format);
                break;
            case SSFF_UNDEF:
                break;
            }
            free((void *) format);
        }
    }
    /*
     * prepare ddescs
     * check tracks and formats and are there enough data parts?
     */
    FIRST = 1;
    if (isNull(tracks = getAttrib(rdobj, R_NamesSymbol)) ||
        LENGTH(tracks) == 0 || TYPEOF(tracks) != STRSXP) {
        freeDObj(dop);
        error("There are no data tracks!");
    }

    if (isNull(formats = getAttrib(rdobj, install("trackFormats"))) ||
        TYPEOF(tracks) != STRSXP) {
        freeDObj(dop);
        error("There are no track format specifiers!");
    }
    if (LENGTH(tracks) > LENGTH(formats)) {
        freeDObj(dop);
        error("Not enough format specifiers for the data tracks.");
    }

    for (i = 0; i < LENGTH(tracks); i++) {
        /*
         * get dimensions
         */
        track = VECTOR_ELT(rdobj, i);
        attr = getAttrib(track, R_DimSymbol);
        /*
         * if there is more than one track, add descriptor
         */
        if (FIRST) {
            dop->numRecords = INTEGER(attr)[0];
            FIRST = 0;
        } else {
            desc = addDDesc(dop);
            if (desc == NULL) {
                freeDObj(dop);
                error("%s", getAsspMsg(asspMsgNum));
            }
            if (dop->numRecords != INTEGER(attr)[0]) {
                freeDObj(dop);
                error("Dimensions of tracks do not match."
                      "(%ld rows in first track, but %d rows in track %d).",
                      dop->numRecords, INTEGER(attr)[0], i);
            }
        }
        desc->ident = strdup(CHAR(STRING_ELT(tracks, i)));
        format = strdup(CHAR(STRING_ELT(formats, i)));
        entry = keyword2entry(desc->ident, KDT_SSFF);   /* search SSFF
                                                         * info */
        if (entry != NULL) {
            desc->type = entry->dataType;
            if (entry->factor != NULL)
                strcpy(desc->factor, entry->factor);
            if (entry->unit != NULL)
                strcpy(desc->unit, entry->unit);
        } else {
            /* For audio data without keyword entry, assume sampled data */
            desc->type = DT_SMP;
        }
        if (strcmp(format, "BIT") == 0) {
            desc->format = DF_BIT;
            desc->numBits = 1;
        } else if (strcmp(format, "STR") == 0) {
            desc->format = DF_STR;
            desc->numBits = 1;
        } else if (strcmp(format, "CHAR") == 0) {
            desc->format = DF_CHAR;
            desc->numBits = 8;
        } else if (strcmp(format, "UINT8") == 0) {
            desc->format = DF_UINT8;
            desc->numBits = 8;
        } else if (strcmp(format, "INT8") == 0) {
            desc->format = DF_INT8;
            desc->numBits = 8;
        } else if (strcmp(format, "UINT16") == 0) {
            desc->format = DF_UINT16;
            desc->numBits = 16;
        } else if (strcmp(format, "INT16") == 0) {
            desc->format = DF_INT16;
            desc->numBits = 16;
        } else if (strcmp(format, "UINT32") == 0) {
            desc->format = DF_UINT32;
            desc->numBits = 32;
        } else if (strcmp(format, "INT32") == 0) {
            desc->format = DF_INT32;
            desc->numBits = 32;
        } else if (strcmp(format, "UINT64") == 0) {
            desc->format = DF_UINT64;
            desc->numBits = 64;
        } else if (strcmp(format, "INT64") == 0) {
            desc->format = DF_INT64;
            desc->numBits = 64;
        } else if (strcmp(format, "REAL32") == 0) {
            desc->format = DF_REAL32;
            desc->numBits = 32;
        } else if (strcmp(format, "REAL64") == 0) {
            desc->format = DF_REAL64;
            desc->numBits = 64;
        } else {
            freeDObj(dop);
            error("Cannot handle data format %s.", format);
        }

        desc->coding = DC_LIN;
        desc->numFields = (size_t) INTEGER(attr)[1];
        free((void *) format);
    }
    setRecordSize(dop);
    dop->frameDur = -1;
    checkRates(dop);
    allocDataBuf(dop, dop->numRecords);
    if (dop->dataBuffer == NULL) {
        freeDObj(dop);
        error("%s", getAsspMsg(asspMsgNum));
    }

    for (i = 0, desc = &(dop->ddl); i < LENGTH(tracks);
         i++, desc = desc->next) {
        track = VECTOR_ELT(rdobj, i);
        if (!addTrackData(dop, desc, track)) {
            freeDObj(dop);
            error("Adding Trackdata did not work...");
        }
    }
    dop->bufNumRecs = dop->numRecords;
    /* startRecord is now correctly 0-indexed, use it directly */
    dop->bufStartRec = dop->startRecord;
    return dop;
}


/*
 * This function takes a SEXP of class AsspDataFormat, turns it into a
 * DOBJ and writes it to file. The DOBJ is deleted after wards. 
 */
SEXP writeDObj_(SEXP data, SEXP fname)
{
    DOBJ           *dop = sexp2dobj(data);
    dop = asspFOpen(strdup(CHAR(STRING_ELT(fname, 0))), AFO_WRITE, dop);
    if (dop == NULL) {
        freeDObj(dop);
        error("%s", getAsspMsg(asspMsgNum));
    }
    asspFWrite(dop->dataBuffer, dop->bufNumRecs, dop);
    asspFClose(dop, AFC_FREE);
    return R_NilValue;
}



/*
 * this function takes trackdata in the form of an R Matrix (rdobj) and
 * adds its contents to a DOBJ in correspondence with a given data
 * descriptor. 
 */
int
addTrackData(DOBJ * dop, DDESC * ddl, SEXP rdobj)
{
    void           *bufPtr;
    int             i,
                    m,
                    n,
                    unp = 0;
    /*
     * various pointers for variuos data sizes
     */
    uint8_t        *u8Ptr;
    int8_t         *i8Ptr;
    uint16_t       *u16Ptr;
    int16_t        *i16Ptr;
    uint32_t       *u32Ptr;
    int32_t        *i32Ptr;
    float          *f32Ptr;
    double         *f64Ptr;

    SEXP            numMat;
    double         *numPtr;
    uint8_t        *bPtr;

    if (isReal(rdobj))
        numMat = rdobj;
    else if (isInteger(rdobj)) {
        PROTECT(numMat = coerceVector(rdobj, REALSXP));
        unp++;
    } else
        error("Bad data type, must be INTEGER or REAL.");
    numPtr = REAL(numMat);

    i = 0;                      /* initial index in buffer */

    for (m = 0; m < dop->numRecords; m++) {
        bufPtr = (void *)((char *)dop->dataBuffer + m * dop->recordSize);
        bPtr = (uint8_t *) bufPtr;
        switch (ddl->format) {
        case DF_UINT8:
            {
                u8Ptr = &bPtr[ddl->offset];
                for (n = 0; n < ddl->numFields; n++) {
                    double v = numPtr[m + n * dop->numRecords];
                    if (ISNAN(v))
                        v = 0.0;            /* SSFF cannot store NA/NaN */
                    u8Ptr[n] = (uint8_t) v;
                }
            }
            break;
        case DF_INT8:
            {
                i8Ptr = (int8_t *) & bPtr[ddl->offset];
                for (n = 0; n < ddl->numFields; n++) {
                    double v = numPtr[m + n * dop->numRecords];
                    if (ISNAN(v))
                        v = 0.0;            /* SSFF cannot store NA/NaN */
                    i8Ptr[n] = (int8_t) v;
                }
            }
            break;
        case DF_UINT16:
            {
                u16Ptr = (uint16_t *) & bPtr[ddl->offset];
                for (n = 0; n < ddl->numFields; n++) {
                    double v = numPtr[m + n * dop->numRecords];
                    if (ISNAN(v))
                        v = 0.0;            /* SSFF cannot store NA/NaN */
                    u16Ptr[n] = (uint16_t) v;
                }
            }
            break;
        case DF_INT16:
            {
                i16Ptr = (int16_t *) & bPtr[ddl->offset];
                for (n = 0; n < ddl->numFields; n++) {
                    double v = numPtr[m + n * dop->numRecords];
                    if (ISNAN(v))
                        v = 0.0;            /* SSFF cannot store NA/NaN */
                    i16Ptr[n] = (int16_t) v;
                }
            }
            break;
        case DF_UINT32:
            {
                u32Ptr = (uint32_t *) & bPtr[ddl->offset];
                for (n = 0; n < ddl->numFields; n++) {
                    double v = numPtr[m + n * dop->numRecords];
                    if (ISNAN(v))
                        v = 0.0;            /* SSFF cannot store NA/NaN */
                    u32Ptr[n] = (uint32_t) v;
                }
            }
            break;
        case DF_INT32:
            {
                i32Ptr = (int32_t *) & bPtr[ddl->offset];
                for (n = 0; n < ddl->numFields; n++) {
                    double v = numPtr[m + n * dop->numRecords];
                    if (ISNAN(v))
                        v = 0.0;            /* SSFF cannot store NA/NaN */
                    i32Ptr[n] = (int32_t) v;
                }
            }
            break;
        case DF_REAL32:
            {
                f32Ptr = (float *) &bPtr[ddl->offset];
                for (n = 0; n < ddl->numFields; n++) {
                    double v = numPtr[m + n * dop->numRecords];
                    if (ISNAN(v))
                        v = 0.0;            /* SSFF cannot store NA/NaN */
                    f32Ptr[n] = (float) v;
                }
            }
            break;
        case DF_REAL64:
            {
                f64Ptr = (double *) &bPtr[ddl->offset];
                for (n = 0; n < ddl->numFields; n++) {
                    double v = numPtr[m + n * dop->numRecords];
                    if (ISNAN(v))
                        v = 0.0;            /* SSFF cannot store NA/NaN */
                    f64Ptr[n] = (double) v;
                }
            }
            break;
        default:
            error
                ("Hi, I just landed in the default of a switch in dataobj.c."
                 "I am sorry, I should not be here and I don't know what to do.");
            break;
        }
    }

    UNPROTECT(unp);
    return 1;
}
