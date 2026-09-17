/*
 * A Snack MP3 format handler using libmpg123.
 *
 * BSD Copyright 2009 - Peter MacDonald
 *
 * Implements mp3 as a loadable module using libmpg123 (which is LGPL).
 * Replaces snacks builtin MP3 driver which:
 *
 *   - has noise major artifacts when used with -file on 48000 sound cards.
 *   - fails on small mp3 files (under 20k?).
 *   - has a restrictive (non-commercial only) licence
 *   - Can't be easily removed from snack (to avoid patent issue)
 * 
 * TODO:
 *   - Check if file changed on multiple opens (for header).
 *   - Check return codes.
 *   - Add encoding support option (ie. use lame library?).
 *
 */
#include <math.h>
#include <tcl.h>
#include "snack.h"
#include <stdlib.h>
#include <time.h>
#include "mpg123.h"
#include <string.h>
#include <ctype.h>

#if defined(__WIN32__)
#  include <io.h>
#  include <fcntl.h>
#  define WIN32_LEAN_AND_MEAN
#  include <windows.h>
#  undef WIN32_LEAN_AND_MEAN
#  define EXPORT(a,b) __declspec(dllexport) a b
BOOL APIENTRY
DllMain(HINSTANCE hInst, DWORD reason, LPVOID reserved)
{
    return TRUE;
}
#else
#  define EXPORT(a,b) a b
#endif

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

#include <stdio.h>

/* #define MPG_NODIRECT_FILES 1 */
/* If above is defined we never let libmpg123 use native files directly. */

#define MPG123_STRING "MPG"
#define SNACK_MPG123_INT 21
#define DECODEBUFSIZE (10*BUFSIZ)
#define READBUFSIZE 8500

typedef struct Mpg123_File {
    mpg123_handle *m;
    int maxbitrate;
    int minbitrate;
    int nombitrate;
    double quality;
    long rate;
    int channels, enc;
    mpg123_id3v1 *v1;
    mpg123_id3v2 *v2;
    Tcl_Obj *fname, *nfname;
    struct mpg123_frameinfo fi;
    int ref;
    size_t savepos[10];
    int lastret;
    Tcl_Channel datasource;
    long ttllen;
    int isFile;
    int noFiles;
    char *chanType;
    int started;
    int opened;
    int gotformat;
    unsigned char *pcmbuf;
    int buffer_size;
    off_t current_frame;
    off_t frames_left;
    double current_seconds;
    double seconds_left;
    int seeksync;
} Mpg123_File;

#ifdef __cplusplus
}
#endif /* __cplusplus */

static int mpgIsInit = 0;


Mpg123_File * AllocMpg(Sound *s) {
    Mpg123_File *of;
    of = (Mpg123_File*) ckalloc(sizeof(Mpg123_File));
    memset(of, 0, sizeof(Mpg123_File));
    s->extHead2 = (char *) of;
    s->extHead2Type = SNACK_MPG123_INT;
    of->nombitrate = 128000;
    of->maxbitrate = -1;
    of->minbitrate = -1;
    of->quality = -1.0;
    of->seeksync = 5000;
#ifdef MPG_NODIRECT_FILES
    of->noFiles = 1;
#endif
    return of;
}

Mpg123_File *MpgObj(Sound *s) {
    Mpg123_File *of = (Mpg123_File *)s->extHead2;
    if (of == NULL) {
        of = AllocMpg(s);
    }
    return of;
}

static int guessByMagic = 1;
/* If above is 1 GuessMpg123File() looks only at header magic alone. */
/* Meaning the header bits start with 0xFFF or the string ID3 or RIFF. */
/* This avoids the more expensive decoding first chunk to see if we have mp3. */

char *
GuessMpg123File(char *buf, int len)
{
    long rate;
    int channels, enc;
    int fnd = 0, ret; size_t done;
    mpg123_handle *m;
    unsigned char *ubuf = buf;
    unsigned char pcmout[4*sizeof(short)*20000];
    int decsiz = 4*sizeof(short)*20000;

    if (len < 4) return(QUE_STRING);
    if ((ubuf[0] == 0xff && (ubuf[1]&0xf0) == 0xf0)) {
        return MPG123_STRING;
    }
    if (buf[0] == 'I' && buf[1] == 'D' && buf[2] == '3') {
        return MPG123_STRING;
    }
    if (len > 20 && buf[20] == 0x55 &&
        toupper(buf[0]) == 'R' && toupper(buf[0]) == 'I' &&
        toupper(buf[0]) == 'F' && toupper(buf[0]) == 'F') {
        return(MP3_STRING);
    }
    if (guessByMagic) {
        return NULL;
    }
    if (!mpgIsInit) {
        mpgIsInit = 1;
        mpg123_init();
    }
    m = mpg123_new(NULL, &ret);
    if(m == NULL)
    {
        fprintf(stderr, "mp3 fail\n" );
        return NULL;
    }
    mpg123_open_feed(m);
    ret = mpg123_decode(m, buf, len, pcmout, decsiz, &done);
    if (ret != MPG123_ERR ) {
        ret = mpg123_getformat(m, &rate, &channels, &enc);
        if (channels<=0) {
            ret = MPG123_ERR;
        }
    }
    mpg123_delete(m);
    if (ret != MPG123_ERR) {
        return(MPG123_STRING);
    }
    return NULL;
}

char *
ExtMpg123File(char *s)
{
    int l1 = strlen(".mp3");
    int l2 = strlen(s);

    if (strncasecmp(".mp3", &s[l2 - l1], l1) == 0) {
        return(MPG123_STRING);
    }
    return(NULL);
}


static int
Mpg123Setup(Sound *s, Tcl_Interp *interp, Tcl_Channel ch)
{
    mpg123_handle *m;
    Mpg123_File *of;
    int ret, fd, rc;
    long mlen;
    const Tcl_ChannelType *cType;

    of = MpgObj(s);
    of->isFile = 0;
    Tcl_SetChannelOption(interp, ch, "-translation", "binary");
#ifdef TCL_81_API
    Tcl_SetChannelOption(interp, ch, "-encoding", "binary");
#endif
    cType = Tcl_GetChannelType(ch);
    if (of->noFiles == 0 && of->opened) {
        of->isFile = !strcmp(cType->typeName, "file");
    }
    if (s->debug)
        fprintf(stderr, "CHANTYPE(%d,%d): %s, BUF=%d\n", of->isFile, of->noFiles, cType->typeName, DECODEBUFSIZE);

    if (!mpgIsInit) {
        mpgIsInit = 1;
        mpg123_init();
    }
    m = of->m;
    /* TODO: check file name didn't change */
    if (m != NULL) {
        /* If used with */
        if (of->ref<10) {
            if (of->isFile){
                of->savepos[of->ref] = mpg123_tell(m);
            } else {
            }
        }
        of->ref++;
    }
    if (of->isFile){
        of->fname = Tcl_NewStringObj(s->fcname, -1);
        Tcl_IncrRefCount(of->fname);
        of->nfname = Tcl_FSGetNormalizedPath(interp, of->fname);
    } else {
        of->lastret = MPG123_NEED_MORE;
    }
    of->datasource = ch;
    m = mpg123_new(NULL, &ret);
    if(m == NULL) {
        Tcl_AppendResult(interp, "Unable to create mpg123 handle: ", mpg123_plain_strerror(ret), 0);
        return TCL_ERROR;
    }
    of->m = m;
    if (of->isFile){
        if (mpg123_open(m, Tcl_GetString(of->nfname)) != MPG123_OK) {
            Tcl_AppendResult(interp, "Open mpg123 failed: ", mpg123_plain_strerror(ret), 0);
            return TCL_ERROR;
        }
        if (s->debug) mpg123_param(m, MPG123_VERBOSE, 2, 0);
        if (s->debug == 0 ) mpg123_param(m, MPG123_ADD_FLAGS, MPG123_QUIET, 0);
        /*mpg123_param(m, MPG123_ADD_FLAGS, MPG123_SEEKBUFFER, 0);
        mpg123_param(m, MPG123_ADD_FLAGS, MPG123_FUZZY, 0);
        mpg123_param(m, MPG123_REMOVE_FLAGS, MPG123_GAPLESS, 0);*/
    } else {
        mpg123_open_feed(m);
    }
    if (of->pcmbuf)  ckfree( of->pcmbuf );
    of->buffer_size = mpg123_outblock( m );
    of->pcmbuf = ckalloc( of->buffer_size );
    mlen = (long)mpg123_length(m);
    if (mlen<=0) {
        return TCL_OK;
    }
    of->gotformat = 1;
    Snack_SetLength(s, mlen);
    mpg123_info(of->m, &of->fi);
    mpg123_getformat(of->m, &of->rate, &of->channels, &of->enc);
    if (s->debug) fprintf(stderr, "MPG FORMAT: channels=%d, rate=%ld enc=0x%x\n", of->channels, of->rate, of->enc);
    Snack_SetSampleRate(s, of->rate);
    Snack_SetNumChannels(s, of->channels);
    Snack_SetSampleEncoding(s, LIN16);
    of->nombitrate = of->rate;

    rc = mpg123_id3(of->m, &of->v1, &of->v2);

    Snack_SetBytesPerSample(s, 2);
    Snack_SetHeaderSize(s, 0);
    return TCL_OK;
}

static int
OpenMpg123File(Sound *s, Tcl_Interp *interp, Tcl_Channel *ch, char *mode)
{
    mpg123_handle *m;
    Mpg123_File *of;
    int ret, fd, rc;
    long mlen;
    Tcl_ChannelType *cType;

    if (s->debug) fprintf(stderr, "MPG Open: %p : %s\n", s, s->fcname);
    *ch = Tcl_OpenFileChannel(interp, s->fcname, mode, 420);
    if (*ch == NULL) {
        Tcl_AppendResult(interp, "Mpg123: unable to open file: ",
            Snack_GetSoundFilename(s), NULL);
            return TCL_ERROR;
    }
    of = MpgObj(s);
    of->opened = 1;
    return Mpg123Setup(s, interp, *ch);
    
}

static int
FreeRes(Mpg123_File *of)
{
    if (of->fname) {
        Tcl_DecrRefCount(of->fname);
    }
    of->fname = NULL;
    of->nfname = NULL;
    of->v1 = NULL;
    of->v2 = NULL;
    if (of->m) {
        mpg123_delete(of->m);
    }
    if (of->pcmbuf) {
        ckfree(of->pcmbuf);
    }
    of->pcmbuf = NULL;
    of->m = NULL;
}

static int
CloseMpg123File(Sound *s, Tcl_Interp *interp, Tcl_Channel *ch)
{
    Mpg123_File *of;

    of = MpgObj(s);

    if (s->debug) fprintf(stderr, "MPG Close: %p\n", s);
    if (of->ref > 0 && of->m) {
        of->ref--;
        if (of->ref<10) {
            if (of->isFile){
                mpg123_seek(of->m, of->savepos[of->ref], SEEK_SET);
            }
        }
        return TCL_OK;
    }

    FreeRes(of);
    if (of->started == 0) {
        *ch = NULL;
    } else {
        of->started = 0;
    }
    if (ch != NULL) {
        Tcl_Close(interp, *ch);
    }
    *ch = NULL;

    return TCL_OK;
}

static int
ReadMpg123Samples(Sound *s, Tcl_Interp *interp, Tcl_Channel ch, char *ibuf,
          float *obuf, int len)
{
    Mpg123_File *of;
    int i, iread, rc, cnt, bigendian = Snack_PlatformIsLittleEndian() ? 0 : 1;
    float *f  = obuf;
    size_t done = 0, rlen, nread = 0;
    short *r;
    long bytes;
    char buffer[READBUFSIZE];

    of = MpgObj(s);
    of->started = 1;
    memset(obuf, 0, len*sizeof(float));
    /*rlen = (len<DECODEBUFSIZE?len:DECODEBUFSIZE); */
    while (len>0) {
        rlen = (len * sizeof(short));
        if (rlen>of->buffer_size) rlen = of->buffer_size;
        if (of->isFile) {
            rc = mpg123_read(of->m, of->pcmbuf, rlen, &done);
        } else {
            if (of->lastret == MPG123_NEED_MORE) {
                bytes=Tcl_Read(of->datasource, buffer, READBUFSIZE); 
                if (bytes <= 0) {
                    if (s->debug) fprintf(stderr, "MPG ERR\n");
                    return 0;
                }
                rc = mpg123_decode(of->m, buffer, bytes, of->pcmbuf, rlen, &done);
            } else {
                rc = mpg123_decode(of->m, NULL, 0, of->pcmbuf, rlen, &done);
            }
        }
        of->lastret = rc;
        if (rc == MPG123_NEW_FORMAT || !of->gotformat) {
            of->gotformat = 1;
            mpg123_getformat(of->m, &of->rate, &of->channels, &of->enc);
            if (s->debug) fprintf(stderr, "MPG FORMAT: channels=%d, rate=%ld enc=0x%x\n", of->channels, of->rate, of->enc);
            Snack_SetSampleRate(s, of->rate);
            Snack_SetNumChannels(s, of->channels);
        }
        if (rc == MPG123_DONE) {
            if (s->debug) fprintf(stderr, "MPG DONE: %zu\n", nread);
            return nread;
        }
        if (rc == MPG123_ERR) {
            if (s->debug) fprintf(stderr, "MPG ERROR: %zu\n", nread);
            return 0;
        }
        r = (short *) of->pcmbuf;
        cnt = (done / sizeof(short));
        for (i = 0; i < cnt; i++) {
            float fv = (float)*r;
            *f++ = fv;
            r++;
        }
        nread += cnt;
        if (s->debug) fprintf(stderr, "MPG READ (%zu of %d): %d\n", nread, len, rc);
        if (cnt >= len) break;
        len -= cnt;
    }

    if (nread>0) {
        of->ttllen += nread;
        if (of->isFile == 0) {
            /* Don't know length for channels so we lie. */
            Snack_SetLength(s, of->ttllen+1);
        }
    }
    if (done < 0) {
        return 0;
    }
    if (of->isFile) {
        mpg123_position(of->m, 0, nread * sizeof(short),
            &of->current_frame, &of->frames_left, &of->current_seconds,
            &of->seconds_left);
    }

    iread = (int)nread;
    if (iread < 0)
    iread = 1;
    if (s->debug) fprintf(stderr, "MPG READ RET: %zu\n", nread);
    return iread;
}

/* SeekMpg123File:
 *
 * Seek to sound-sample position.
 * This happens a lot when global rate does not equal decode rate.
 * We work around quanticize bug (resyncing?) 
 * by seeking back an extra amount, then read forward again by that amount. 
 * TODO: check return codes.
 */
static int
SeekMpg123File(Sound *s, Tcl_Interp *interp, Tcl_Channel ch, int pos)
{
    Mpg123_File *of;
    int opos;

    of = MpgObj(s);
    if (s->debug) fprintf(stderr, "MPG SEEK: %d\n", pos);
    if (of->started  == 0 && pos == 0) {
        if (s->debug) fprintf(stderr, "MPG SEEK SKIPPED\n");
        return pos;
    }
    opos = mpg123_tell(of->m);
    if (pos == opos) {
        if (s->debug) fprintf(stderr, "MPG SEEK NOMOVE: %d->%d\n", opos, pos);
    }
    opos = pos;

    if (of->datasource) {
        int extra = (pos>of->seeksync?of->seeksync:pos);
        size_t done;
        
        if (of->isFile) {
            if (of->seeksync > 0 && extra > 0) {
                mpg123_seek(of->m, pos-extra, SEEK_SET);
                mpg123_read(of->m, of->pcmbuf, extra, &done);
            } else {
                mpg123_seek(of->m, pos, SEEK_SET);
            }
        } else {
            off_t ioffs;
            if (of->seeksync > 0 && extra > 0) {
                mpg123_feedseek(of->m, pos-extra, SEEK_SET, &ioffs);
                Tcl_Seek(of->datasource, ioffs, SEEK_SET);
                Tcl_Read(of->datasource, of->pcmbuf, extra);
                mpg123_decode(of->m, of->pcmbuf, extra, NULL, 0, &done);
                mpg123_decode(of->m, NULL, 0, of->pcmbuf, extra, &done);
            } else {
                mpg123_feedseek(of->m, pos, SEEK_SET, &ioffs);
                Tcl_Seek(of->datasource, ioffs, SEEK_SET);
            }
        }
    }
    pos = mpg123_tell(of->m);
    if (s->debug) fprintf(stderr, "MPG SEEKPOS: %d -> %d\n", opos, pos);
    if (pos<0) {
        return(-1);
    } else {
        return pos;
    }
}

static int
GetMpg123Header(Sound *s, Tcl_Interp *interp, Tcl_Channel ch, Tcl_Obj *obj,
         char *buf)
{
    Mpg123_File *of;
    long mlen;
    size_t done;
    int i, ret, rc;
    mpg123_id3v1 *v1; mpg123_id3v2 *v2;
    
    of = MpgObj(s);
    if (!of->opened) {
        return Mpg123Setup(s, interp, ch);
    }

    if (s->debug) fprintf(stderr, "MPG Header\n");

    /* For the case when Tcl_Open has been done somewhere else */

    if (s->extHead2 != NULL && s->extHead2Type != SNACK_MPG123_INT) {
        Snack_FileFormat *ff;
    
        for (ff = Snack_GetFileFormats(); ff != NULL; ff = ff->nextPtr) {
            if (strcmp(s->fileType, ff->name) == 0) {
                if (ff->freeHeaderProc != NULL) {
                    (ff->freeHeaderProc)(s);
                }
            }
        }
    }
    of = MpgObj(s);
    of->started = 1;
    
    mlen = (long)mpg123_length(of->m);
    if (mlen<=0) {
        return TCL_OK;
        return TCL_ERROR;
    }
    Snack_SetLength(s, mlen);
    mpg123_info(of->m, &of->fi);
    mpg123_getformat(of->m, &of->rate, &of->channels, &of->enc);
    if (s->debug) fprintf(stderr, "MPG FORMAT: channels=%d, rate=%ld enc=0x%x\n", of->channels, of->rate, of->enc);
    Snack_SetSampleRate(s, of->rate);
    Snack_SetNumChannels(s, of->channels);
    Snack_SetSampleEncoding(s, LIN16);
    of->nombitrate = of->rate;

    rc = mpg123_id3(of->m, &of->v1, &of->v2);

    Snack_SetBytesPerSample(s, 2);
    Snack_SetHeaderSize(s, 0);

    return TCL_OK;
}


void
FreeMpg123Header(Sound *s)
{
    Mpg123_File *of = (Mpg123_File *)s->extHead2;

    if (s->extHead2 != NULL) {
        FreeRes(of);
        ckfree((char *)s->extHead2);
        s->extHead2 = NULL;
        s->extHead2Type = 0;
    }
}

int
ConfigMpg123(Sound *s, Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[])
{
    Mpg123_File *of;
    int arg, index;
    static CONST char *optionStrings[] = {
        "-comment", "-album", "-seeksync",
        "-artist", "-year", "-tag", "-title", "-genre",
        "-maxbitrate", "-minbitrate", "-nominalbitrate",
        "-quality", "-nofiles", "-magiconly",  "-played", "-remain", NULL
    };
    enum options {
        COMMENT, ALBUM, SEEKSYNC, ARTIST, YEAR, TAG, TITLE, GENRE, MAX, MIN, NOMINAL, QUALITY, NOFILES, USEMAGIC, SECONDS, REMAIN
    };
  
    of = MpgObj(s);

    if (s->extHead2 != NULL && s->extHead2Type != SNACK_MPG123_INT) {
        Snack_FileFormat *ff;
    
        for (ff = Snack_GetFileFormats(); ff != NULL; ff = ff->nextPtr) {
            if (strcmp(s->fileType, ff->name) == 0) {
                if (ff->freeHeaderProc != NULL) {
                    (ff->freeHeaderProc)(s);
                }
            }
        }
    }
  
    if (objc < 3) return 0;


    if (objc == 3) { 
        /* get option */
        if (Tcl_GetIndexFromObj(interp, objv[2], optionStrings, "option", 0,
            &index) != TCL_OK) {
                Tcl_AppendResult(interp, ", or\n", NULL);
                return 0;
        }

#define NSO(str) Tcl_NewStringObj((of->v1 && of->v1->str)?(of->v1->str):"",-1)

        switch ((enum options) index) {
            case COMMENT:
            {
                Tcl_SetObjResult(interp, NSO(comment));
                break;
            }
            case ALBUM:
            {
                Tcl_SetObjResult(interp, NSO(album));
                break;
            }
            case TITLE:
            {
                Tcl_SetObjResult(interp, NSO(title));
                break;
            }
            case TAG:
            {
                Tcl_SetObjResult(interp, NSO(tag));
                break;
            }
            case YEAR:
            {
                Tcl_SetObjResult(interp, NSO(year));
                break;
            }
            case GENRE:
            {
                if (of->v1)
                Tcl_SetObjResult(interp, Tcl_NewIntObj(of->v1?of->v1->genre:-1));
                break;
            }
            case NOFILES:
            {
                Tcl_SetObjResult(interp, Tcl_NewIntObj(of->noFiles));
                break;
            }
            case MAX:
            {
                Tcl_SetObjResult(interp, Tcl_NewIntObj(of->maxbitrate));
                break;
            }
            case MIN:
            {
                Tcl_SetObjResult(interp, Tcl_NewIntObj(of->minbitrate));
                break;
            }
            case NOMINAL:
            {
                Tcl_SetObjResult(interp, Tcl_NewIntObj(of->nombitrate));
                break;
            }
            case SEEKSYNC:
            {
                Tcl_SetObjResult(interp, Tcl_NewIntObj(of->seeksync));
                break;
            }
            case REMAIN:
            {
                Tcl_SetObjResult(interp, Tcl_NewIntObj(of->seconds_left));
                break;
            }
            case SECONDS:
            {
                Tcl_SetObjResult(interp, Tcl_NewIntObj(of->current_seconds));
                break;
            }
            case QUALITY:
            {
                Tcl_SetObjResult(interp, Tcl_NewDoubleObj(of->quality));
                break;
            }
            case USEMAGIC:
            {
                Tcl_SetObjResult(interp, Tcl_NewIntObj(guessByMagic));
                break;
            }
        }
    } else { 
        /* set option */
        for (arg = 2; arg < objc; arg+=2) {
            int index;
      
            if (Tcl_GetIndexFromObj(interp, objv[arg], optionStrings, "option", 0,
                &index) != TCL_OK) {
                    return 0;
            }
      
            if (arg + 1 == objc) {
                Tcl_AppendResult(interp, "No argument given for ",
                    optionStrings[index], " option\n", (char *) NULL);
                    return 0;
            }
      
            switch ((enum options) index) {
                case NOFILES:
                {
#ifndef MPG_NODIRECT_FILES
                    if (Tcl_GetIntFromObj(interp,objv[arg+1], &of->noFiles) != TCL_OK)
#endif
                    return 0;
                    break;
                }
                case COMMENT:
                {
                    int i, n;
                    break;
                }
                case MAX:
                {
                    if (Tcl_GetIntFromObj(interp,objv[arg+1], &of->maxbitrate) != TCL_OK)
                    return 0;
                    break;
                }
                case MIN:
                {
                    if (Tcl_GetIntFromObj(interp,objv[arg+1], &of->minbitrate) != TCL_OK)
                    return 0;
                    break;
                }
                case NOMINAL:
                {
                    if (Tcl_GetIntFromObj(interp,objv[arg+1], &of->nombitrate) != TCL_OK)
                    return 0;
                    break;
                }
                case SEEKSYNC:
                {
                    if (Tcl_GetIntFromObj(interp,objv[arg+1], &of->seeksync) != TCL_OK)
                    return 0;
                    break;
                }
                case USEMAGIC:
                {
                    if (Tcl_GetIntFromObj(interp,objv[arg+1], &guessByMagic) != TCL_OK)
                    return 0;
                    break;
                }
                case QUALITY:
                {
                    if (Tcl_GetDoubleFromObj(interp, objv[arg+1], &of->quality) !=TCL_OK)
                    return 0;
                    break;
                }
            }
        }
    }

    return 1;
}

#define MPG123FILE_VERSION "1.3"

Snack_FileFormat snackMpg123Format = {
    MPG123_STRING,
        GuessMpg123File,
        GetMpg123Header,
        ExtMpg123File,
        NULL, /* PutMpg123Header, */
        OpenMpg123File,
        CloseMpg123File,
        ReadMpg123Samples,
        NULL, /* WriteMpg123Samples, */
        SeekMpg123File,
        FreeMpg123Header,
        ConfigMpg123,
        (Snack_FileFormat *) NULL
};

/* Called by "load libsnackmpg" */
EXPORT(int, Snackmpg_Init) _ANSI_ARGS_((Tcl_Interp *interp))
{
    int res;
  
#ifdef USE_TCL_STUBS
    if (Tcl_InitStubs(interp, TCL_VERSION, 0) == NULL) {
        return TCL_ERROR;
    }
#endif
  
#ifdef USE_SNACK_STUBS
    if (Snack_InitStubs(interp, "2", 0) == NULL) {
        return TCL_ERROR;
    }
#endif
  
    res = Tcl_PkgProvide(interp, "snackmpg", MPG123FILE_VERSION);
  
    if (res != TCL_OK) return res;

    Tcl_SetVar(interp, "snack::snackmpg", MPG123FILE_VERSION, TCL_GLOBAL_ONLY);

    Snack_CreateFileFormat(&snackMpg123Format);

    return TCL_OK;
}

EXPORT(int, Snackmpg_SafeInit)(Tcl_Interp *interp)
{
    return Snackmpg_Init(interp);
}
