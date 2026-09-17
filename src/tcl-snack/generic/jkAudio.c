/* 
 * Copyright (C) 1997-2004 Kare Sjolander <kare@speech.kth.se>
 *
 * This file is part of the Snack Sound Toolkit.
 * The latest version can be found at http://www.speech.kth.se/snack/
 *
 * The authors hereby grant permission to use, copy, modify, distribute,
 * and license this software and its documentation for any purpose, provided
 * that existing copyright notices are retained in all copies and that this
 * notice is included verbatim in any distributions. No written agreement,
 * license, or royalty fee is required for any of the authorized uses.
 * Modifications to this software may be copyrighted by their authors
 * and need not follow the licensing terms described here, provided that
 * the new terms are clearly indicated on the first page of each file where
 * they apply.
 *
 * IN NO EVENT SHALL THE AUTHORS OR DISTRIBUTORS BE LIABLE TO ANY PARTY
 * FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES
 * ARISING OUT OF THE USE OF THIS SOFTWARE, ITS DOCUMENTATION, OR ANY
 * DERIVATIVES THEREOF, EVEN IF THE AUTHORS HAVE BEEN ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * THE AUTHORS AND DISTRIBUTORS SPECIFICALLY DISCLAIM ANY WARRANTIES,
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE, AND NON-INFRINGEMENT.  THIS SOFTWARE
 * IS PROVIDED ON AN "AS IS" BASIS, AND THE AUTHORS AND DISTRIBUTORS HAVE
 * NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
 * MODIFICATIONS.
 */

#include <stdlib.h>
#include <stdio.h>
#include <signal.h>
#include <math.h>
#include <string.h>
#include "tcl.h"
#include "snack.h"

extern int rop, wop;
extern double startDevTime;
extern struct jkQueuedSound *soundQueue;
extern struct jkQueuedSound *rsoundQueue;

#if defined(HPUX) || defined(MAC)
/* Choosing a good generic value for HP-UX is not easy */
#  define BUFSECS 2.0
#else
#  define BUFSECS 0.25
#endif

double globalLatency = BUFSECS;
float globalScaling = 1.0f;
int failRate = 48000;

char defaultOutDevice[MAX_DEVICE_NAME_LENGTH];
char defaultInDevice[MAX_DEVICE_NAME_LENGTH];

char *
SnackStrDup(const char *str)
{
  size_t len = strlen(str);
  char *new = ckalloc(len + 1);

  if (new) {
    memcpy(new, str, len + 1);
  }

  return new;
}

static int
outDevicesCmd(Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[])
{
  int i, n;
  char *arr[MAX_NUM_DEVICES];
  Tcl_Obj *list = Tcl_NewListObj(0, NULL);

  n = SnackGetOutputDevices(arr, MAX_NUM_DEVICES);

  for (i = 0; i < n; i++) {
    Tcl_ListObjAppendElement(interp, list, Tcl_NewStringObj(arr[i], -1));
    ckfree(arr[i]);
  }

  Tcl_SetObjResult(interp, list);

  return TCL_OK;
}

static int
inDevicesCmd(Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[])
{
  int i, n;
  char *arr[MAX_NUM_DEVICES];
  Tcl_Obj *list = Tcl_NewListObj(0, NULL);

  n = SnackGetInputDevices(arr, MAX_NUM_DEVICES);

  for (i = 0; i < n; i++) {
    Tcl_ListObjAppendElement(interp, list, Tcl_NewStringObj(arr[i], -1));
    ckfree(arr[i]);
  }

  Tcl_SetObjResult(interp, list);

  return TCL_OK;
}

static int
selectOutCmd(Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[])
{
  int i, n, found = 0;
  char *arr[MAX_NUM_DEVICES];
  char *devstr;

  n = SnackGetOutputDevices(arr, MAX_NUM_DEVICES);

  if (objc == 3) {
    devstr = Tcl_GetStringFromObj(objv[2], NULL);
    for (i = 0; i < n; i++) {
      if (strncmp(devstr, arr[i], strlen(devstr)) == 0 && found == 0) {
	strncpy(defaultOutDevice, arr[i], MAX_DEVICE_NAME_LENGTH - 1);
	defaultOutDevice[MAX_DEVICE_NAME_LENGTH - 1] = '\0';
	found = 1;
      }
      ckfree(arr[i]);
    }
    if (found == 0) {
      Tcl_AppendResult(interp, "No such device: ", devstr, (char *) NULL);
      return TCL_ERROR;
    }
  } else {
    Tcl_WrongNumArgs(interp, 1, objv, "selectOutput device");
    return TCL_ERROR;
  }
  
  return TCL_OK;
}

static int
selectInCmd(Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[])
{
  int i, n, found = 0;
  char *arr[MAX_NUM_DEVICES];
  char *devstr;

  n = SnackGetInputDevices(arr, MAX_NUM_DEVICES);

  if (objc == 3) {
    devstr = Tcl_GetStringFromObj(objv[2], NULL);
    for (i = 0; i < n; i++) {
      if (strncmp(devstr, arr[i], strlen(devstr)) == 0 && found == 0) {
	strncpy(defaultInDevice, arr[i], MAX_DEVICE_NAME_LENGTH - 1);
	defaultInDevice[MAX_DEVICE_NAME_LENGTH - 1] = '\0';
	found = 1;
      }
      ckfree(arr[i]);
    }
    if (found == 0) {
      Tcl_AppendResult(interp, "No such device: ", devstr, (char *) NULL);
      return TCL_ERROR;
    }
  } else {
    Tcl_WrongNumArgs(interp, 1, objv, "selectInput device");
    return TCL_ERROR;
  }
  
  return TCL_OK;
}

static int
encodingsCmd(Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[])
{
  char *str = "Lin16 Mulaw Alaw Lin8offset Lin8 Lin24 Lin24packed Lin32 Float";

  Tcl_SetObjResult(interp, Tcl_NewStringObj(str, -1));

  return TCL_OK;
}

static int
ratesCmd(Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[])
{
  char tmpstr[QUERYBUFSIZE];

  SnackAudioGetRates(defaultOutDevice, tmpstr, QUERYBUFSIZE);
  Tcl_SetObjResult(interp, Tcl_NewStringObj(tmpstr, -1));

  return TCL_OK;
}

static int
activeCmd(Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[])
{
  if (wop == IDLE && rop == IDLE) {
    Tcl_SetObjResult(interp, Tcl_NewIntObj(0));
  } else {
    Tcl_SetObjResult(interp, Tcl_NewIntObj(1));
  }
  
  return TCL_OK;
}

static int
play_gainCmd(Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[])
{
  int g;
  
  if (objc == 3) {
    if (Tcl_GetIntFromObj(interp, objv[2], &g) != TCL_OK) return TCL_ERROR;
    ASetPlayGain(g);
  } else {
#ifdef HPUX
    if (wop == IDLE)
#endif
      Tcl_SetObjResult(interp, Tcl_NewIntObj(AGetPlayGain()));
  }

  return TCL_OK;
}

static int
record_gainCmd(Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[])
{
  int g;
  
  if (objc == 3) {
    if (Tcl_GetIntFromObj(interp, objv[2], &g) != TCL_OK) return TCL_ERROR;
	ASetRecGain(g);
  } else {
#ifdef HPUX
    if (rop == IDLE)
#endif
      Tcl_SetObjResult(interp, Tcl_NewIntObj(AGetRecGain()));
  }

  return TCL_OK;
}

static int
elapsedTimeCmd(Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[])
{
  double elapsedTime = SnackCurrentTime() - startDevTime;

  if (wop == IDLE && rop == IDLE) {
    Tcl_SetObjResult(interp, Tcl_NewDoubleObj(0.0));
  }  else if (wop == PAUSED || rop == PAUSED) {
    Tcl_SetObjResult(interp, Tcl_NewDoubleObj(startDevTime));
  } else {
    Tcl_SetObjResult(interp, Tcl_NewDoubleObj(elapsedTime));
  }

  return TCL_OK;
}

static int
currentSoundCmd(Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[])
{
  jkQueuedSound *p;
  char *res;
  Tcl_HashSearch hashSearch;
  Tcl_HashEntry *entryPtr;

  if (soundQueue == NULL) {
    Tcl_SetObjResult(interp, Tcl_NewStringObj("", -1));
    return TCL_OK;
  }
  for (p = soundQueue; p->next != NULL && p->next->status == SNACK_QS_DONE;
       p = p->next);

  entryPtr = Tcl_FirstHashEntry(p->sound->soundTable, &hashSearch);

  if (p->sound != (Sound *) Tcl_GetHashValue(entryPtr)) {
    entryPtr = Tcl_NextHashEntry(&hashSearch);
  }
  res = Tcl_GetHashKey(p->sound->soundTable, entryPtr);
  Tcl_SetObjResult(interp, Tcl_NewStringObj(res, -1));

  return TCL_OK;
}

static int
playLatencyCmd(Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[])
{
  double d;

  if (objc == 2) {
    Tcl_SetObjResult(interp, Tcl_NewDoubleObj(1000.0*globalLatency));
  } else if (objc == 3) {
    if (Tcl_GetDoubleFromObj(interp, objv[2], &d) != TCL_OK) {
      return TCL_ERROR;
    }
    globalLatency = d / 1000.0;
  } else {
    Tcl_WrongNumArgs(interp, 1, objv, "playLatency ?milliseconds?");
    return TCL_ERROR;
  }
  return TCL_OK;
}

static int
scalingCmd(Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[])
{
  double d = 0.0;

  if (objc == 2) {
    Tcl_SetObjResult(interp, Tcl_NewDoubleObj(globalScaling));
  } else if (objc == 3) {
    if (Tcl_GetDoubleFromObj(interp, objv[2], &d) != TCL_OK) {
      return TCL_ERROR;
    }
    globalScaling = (float) d;
  } else {
    Tcl_WrongNumArgs(interp, 1, objv, "scaling ?factor?");
    return TCL_ERROR;
  }
  return TCL_OK;
}

static int
failrateCmd(Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[])
{
  int d = -1;

  if (objc == 2) {
    Tcl_SetObjResult(interp, Tcl_NewIntObj(failRate));
  } else if (objc == 3) {
    if (Tcl_GetIntFromObj(interp, objv[2], &d) != TCL_OK) {
      return TCL_ERROR;
    }
    globalScaling = d;
  } else {
    Tcl_WrongNumArgs(interp, 1, objv, "altrate ?factor?");
    return TCL_ERROR;
  }
  return TCL_OK;
}

static int
audioPlayCmd(Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[])
{
  if (rop == PAUSED || wop == PAUSED) {
    SnackPauseAudio();
  }

  return TCL_OK;
}

static int
audioStopCmd(Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[])
{
  jkQueuedSound *p;

  if (rop == READ || rop == PAUSED) {
    for (p = rsoundQueue; p != NULL; p = p->next) {
      Snack_StopSound(p->sound, interp);
    }
  }
  if (wop == WRITE || wop == PAUSED) {
    for (p = soundQueue; p != NULL; p = p->next) {
      Snack_StopSound(p->sound, interp);
      /*
       * The soundQueue can be remooved during a stop, so check it
       * otherwise p is garbage
       */
      if (soundQueue == NULL)
	break;
    }
  }

  return TCL_OK;
}

static int
audioPauseCmd(Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[])
{
  SnackPauseAudio();

  return TCL_OK;
}

#define NAUDIOCOMMANDS   17
#define MAXAUDIOCOMMANDS 25

int nAudioCommands   = NAUDIOCOMMANDS;
int maxAudioCommands = MAXAUDIOCOMMANDS;

CONST84 char *audioCmdNames[MAXAUDIOCOMMANDS] = {
  "outputDevices",
  "inputDevices",
  "selectOutput",
  "selectInput",
  "formats",
  "frequencies",
  "active",
  "play_gain",
  "record_gain",
  "elapsedTime",
  "currentSound",
  "playLatency",
  "scaling",
  "encodings",
  "rates",
  "play",
  "stop",
  "pause",
  "fallbackrate",
  NULL
};

/* NOTE: NAUDIOCOMMANDS needs updating when new commands are added. */

audioCmd *audioCmdProcs[MAXAUDIOCOMMANDS] = {
  outDevicesCmd,
  inDevicesCmd,
  selectOutCmd,
  selectInCmd,
  encodingsCmd,
  ratesCmd,
  activeCmd,
  play_gainCmd,
  record_gainCmd,
  elapsedTimeCmd,
  currentSoundCmd,
  playLatencyCmd,
  scalingCmd,
  encodingsCmd,
  ratesCmd,
  audioPlayCmd,
  audioStopCmd,
  audioPauseCmd,
  failrateCmd
};

audioDelCmd *audioDelCmdProcs[MAXAUDIOCOMMANDS] = {
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL
};

int
Snack_AudioCmd(ClientData cdata, Tcl_Interp *interp, int objc,
	       Tcl_Obj *CONST objv[])
{
  int index;

  if (objc < 2) {
    Tcl_WrongNumArgs(interp, 1, objv, "option ?arg?");
    return TCL_ERROR;
  }
  
  if (Tcl_GetIndexFromObj(interp, objv[1], audioCmdNames, "option", 0,
			  &index) != TCL_OK) {
    return TCL_ERROR;
  }

  return((audioCmdProcs[index])(interp, objc, objv)); 
}

void
Snack_AudioDeleteCmd(ClientData clientData)
{
  int i;

  for (i = 0; i < nAudioCommands; i++) {
    if (audioDelCmdProcs[i] != NULL) {
      (audioDelCmdProcs[i])();
    }
  }
}
