#include <R.h>
#include <Rinternals.h>
#include "wrassp.h"
#include "assp/asspfio.h"
#include "assp/asspmess.h"

/*
 * Debug function to print all DOBJ fields for comparison
 */
void debugPrintDObj(DOBJ *dop, const char *label) {
    DDESC *desc;

    Rprintf("\n========== DOBJ Debug: %s ==========\n", label);

    if (dop == NULL) {
        Rprintf("DOBJ is NULL!\n");
        return;
    }

    Rprintf("filePath: %s\n", dop->filePath ? dop->filePath : "(NULL)");
    Rprintf("fp: %p\n", (void*)dop->fp);
    Rprintf("openMode: %d\n", dop->openMode);
    Rprintf("fileFormat: %d\n", (int)dop->fileFormat);
    Rprintf("fileData: %d\n", (int)dop->fileData);
    /* Skip fileEndian - it's not a simple int */
    Rprintf("version: %ld\n", dop->version);
    Rprintf("headerSize: %ld\n", dop->headerSize);
    Rprintf("sampFreq: %.2f\n", dop->sampFreq);
    Rprintf("dataRate: %.2f\n", dop->dataRate);
    Rprintf("frameDur: %ld\n", dop->frameDur);
    Rprintf("recordSize: %zu\n", dop->recordSize);
    Rprintf("startRecord: %ld\n", dop->startRecord);
    Rprintf("numRecords: %ld\n", dop->numRecords);
    Rprintf("Time_Zero: %.6f\n", dop->Time_Zero);
    Rprintf("Start_Time: %.6f\n", dop->Start_Time);
    Rprintf("dataBuffer: %p\n", dop->dataBuffer);
    Rprintf("maxBufRecs: %ld\n", dop->maxBufRecs);
    Rprintf("bufStartRec: %ld\n", dop->bufStartRec);
    Rprintf("bufNumRecs: %ld\n", dop->bufNumRecs);
    Rprintf("bufNeedsSave: %d\n", dop->bufNeedsSave);

    Rprintf("\nData Descriptors:\n");
    int i = 0;
    for (desc = &(dop->ddl); desc != NULL; desc = desc->next, i++) {
        Rprintf("  Descriptor %d:\n", i);
        Rprintf("    ident: %s\n", desc->ident ? desc->ident : "(NULL)");
        Rprintf("    type: %d\n", (int)desc->type);
        Rprintf("    format: %d\n", (int)desc->format);
        Rprintf("    coding: %d\n", (int)desc->coding);
        Rprintf("    numBits: %u\n", desc->numBits);
        Rprintf("    offset: %zu\n", desc->offset);
        Rprintf("    numFields: %zu\n", desc->numFields);
        Rprintf("    unit: %s\n", desc->unit);
        Rprintf("    factor: %s\n", desc->factor);
    }

    Rprintf("========================================\n\n");
}

/*
 * R-callable function to debug DOBJ conversion
 */
SEXP debugDObjConversion(SEXP args) {
    DOBJ *dop;
    SEXP rdobj;

    /* Skip function name in args */
    args = CDR(args);
    rdobj = CAR(args);

    Rprintf("\n=== Starting DOBJ Conversion Debug ===\n");

    /* Convert using sexp2dobj */
    dop = sexp2dobj(rdobj);

    if (dop == NULL) {
        Rprintf("ERROR: sexp2dobj returned NULL!\n");
        return R_NilValue;
    }

    /* Print all fields */
    debugPrintDObj(dop, "From sexp2dobj");

    /* Clean up */
    freeDObj(dop);

    return R_NilValue;
}

/*
 * R-callable function to debug DOBJ from file
 */
SEXP debugDObjFromFile(SEXP args) {
    DOBJ *dop;
    const char *fname_const;
    char *fname;
    SEXP filename;

    /* Skip function name in args */
    args = CDR(args);
    filename = CAR(args);

    fname_const = CHAR(STRING_ELT(filename, 0));
    /* asspFOpen expects non-const char*, so we need to copy it */
    fname = (char*)R_alloc(strlen(fname_const) + 1, sizeof(char));
    strcpy(fname, fname_const);

    Rprintf("\n=== Opening file: %s ===\n", fname);

    /* Open using asspFOpen */
    dop = asspFOpen(fname, AFO_READ, NULL);

    if (dop == NULL) {
        Rprintf("ERROR: asspFOpen returned NULL!\n");
        Rprintf("Error message: %s\n", getAsspMsg(asspMsgNum));
        return R_NilValue;
    }

    /* Print all fields */
    debugPrintDObj(dop, "From asspFOpen");

    /* Clean up */
    asspFClose(dop, AFC_FREE);

    return R_NilValue;
}
