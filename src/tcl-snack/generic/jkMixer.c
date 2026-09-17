/* 
 * Copyright (C) 1997-2002 Kare Sjolander <kare@speech.kth.se>
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

char defaultMixerDevice[MAX_DEVICE_NAME_LENGTH];

static int
devicesCmd(Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[])
{
  int i, n;
  char *arr[MAX_NUM_DEVICES];
  Tcl_Obj *list = Tcl_NewListObj(0, NULL);

  n = SnackGetMixerDevices(arr, MAX_NUM_DEVICES);
  for (i = 0; i < n; i++) {
    Tcl_ListObjAppendElement(interp, list, Tcl_NewStringObj(arr[i], -1));
    ckfree(arr[i]);
  }

  Tcl_SetObjResult(interp, list);
  
  return TCL_OK;
}

static int
selectCmd(Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[])
{
  int i, n, found = 0;
  char *arr[MAX_NUM_DEVICES];
  char *devstr;

  n = SnackGetMixerDevices(arr, MAX_NUM_DEVICES);

  if (objc == 3) {
    devstr = Tcl_GetStringFromObj(objv[2], NULL);
    for (i = 0; i < n; i++) {
      if (strncmp(devstr, arr[i], strlen(devstr)) == 0 && found == 0) {
	strncpy(defaultMixerDevice, arr[i], MAX_DEVICE_NAME_LENGTH - 1);
	defaultMixerDevice[MAX_DEVICE_NAME_LENGTH - 1] = '\0';
	found = 1;
      }
      ckfree(arr[i]);
    }
    if (found == 0) {
      Tcl_AppendResult(interp, "No such device: ", devstr, (char *) NULL);
      return TCL_ERROR;
    }
  } else {
    Tcl_WrongNumArgs(interp, 1, objv, "select device");
    return TCL_ERROR;
  }

  return TCL_OK;
}

static int
inputCmd(Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[])
{
  char *jack, tmpstr[QUERYBUFSIZE];

  if (objc < 3) {
    SnackMixerGetInputJack(tmpstr, QUERYBUFSIZE);
    Tcl_SetObjResult(interp, Tcl_NewStringObj(tmpstr, -1));
  } else {
    jack = Tcl_GetStringFromObj(objv[2], NULL);
    if (objc == 3) {
      if (SnackMixerSetInputJack(interp, jack, "1")) {
	Tcl_AppendResult(interp, "Error setting input jack", NULL);
	return TCL_ERROR;
      };
    } else {
      SnackMixerLinkJacks(interp, jack, objv[3]);
    }
  }
  
  return TCL_OK;
}

static int
inputsCmd(Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[])
{
  char tmpstr[QUERYBUFSIZE];

  SnackMixerGetInputJackLabels(tmpstr, QUERYBUFSIZE);
  Tcl_SetObjResult(interp, Tcl_NewStringObj(tmpstr, -1));
  
  return TCL_OK;
}

static int
outputCmd(Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[])
{
  char *jack, tmpstr[QUERYBUFSIZE];
  
  if (objc < 3) {
    SnackMixerGetOutputJack(tmpstr, QUERYBUFSIZE);
    Tcl_SetObjResult(interp, Tcl_NewStringObj(tmpstr, -1));
  } else {
    jack = Tcl_GetStringFromObj(objv[2], NULL);
    if (objc == 3) {
      SnackMixerSetOutputJack(jack, "1");
    } else {
      SnackMixerLinkJacks(interp, jack, objv[3]);
    }
  }
  
  return TCL_OK;
}

static int
outputsCmd(Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[])
{
  char tmpstr[QUERYBUFSIZE];

  SnackMixerGetOutputJackLabels(tmpstr, QUERYBUFSIZE);
  Tcl_SetObjResult(interp, Tcl_NewStringObj(tmpstr, -1));
  
  return TCL_OK;
}

static int
channelsCmd(Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[])
{
  char *line, tmpstr[QUERYBUFSIZE];
  
  if (objc == 3) {
    line = Tcl_GetStringFromObj(objv[2], NULL);
    SnackMixerGetChannelLabels(line, tmpstr, QUERYBUFSIZE);
    Tcl_SetObjResult(interp, Tcl_NewStringObj(tmpstr, -1));
  } else {
    Tcl_WrongNumArgs(interp, 1, objv, "channels line");
    return TCL_ERROR;
  }

  return TCL_OK;
}

static int
volumeCmd(Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[])
{
  char *line, tmpstr[QUERYBUFSIZE];

  line = Tcl_GetStringFromObj(objv[2], NULL);
#ifdef HPUX
    if (rop != IDLE || wop != IDLE) return TCL_OK;
#endif
  if (objc == 3) {
    SnackMixerGetVolume(line, -1, tmpstr, QUERYBUFSIZE);
    Tcl_SetObjResult(interp, Tcl_NewStringObj(tmpstr, -1));
  } else if (objc == 4) {
  } else if (objc == 5) {
    SnackMixerGetChannelLabels(line, tmpstr, QUERYBUFSIZE);
    if (strcmp("Mono", tmpstr) == 0) {
      Tcl_AppendResult(interp, "Line is single channel", NULL);
      return TCL_ERROR;
    }
  } else {
    Tcl_WrongNumArgs(interp, 1, objv, "audio volume line [leftVar] [rightVar]");
    return TCL_ERROR;
  }
  SnackMixerLinkVolume(interp, line, objc - 3, objv);

  return TCL_OK;
}

static int
linesCmd(Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[])
{
  char tmpstr[QUERYBUFSIZE];

  SnackMixerGetLineLabels(tmpstr, QUERYBUFSIZE);
  Tcl_SetObjResult(interp, Tcl_NewStringObj(tmpstr, -1));

  return TCL_OK;
}

static int
updateCmd(Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[])
{
  SnackMixerUpdateVars(interp);
  
  return TCL_OK;
}

#define NMIXERCOMMANDS   10
#define MAXMIXERCOMMANDS 20

int nMixerCommands   = NMIXERCOMMANDS;
int maxMixerCommands = MAXMIXERCOMMANDS;

CONST84 char *mixerCmdNames[MAXMIXERCOMMANDS] = {
  "devices",
  "select",
  "input",
  "inputs",
  "output",
  "outputs",
  "channels",
  "volume",
  "lines",
  "update",
  NULL
};

/* NOTE: NMIXERCOMMANDS needs updating when new commands are added. */

mixerCmd *mixerCmdProcs[MAXMIXERCOMMANDS] = {
  devicesCmd,
  selectCmd,
  inputCmd,
  inputsCmd,
  outputCmd,
  outputsCmd,
  channelsCmd,
  volumeCmd,
  linesCmd,
  updateCmd,
};

mixerDelCmd *mixerDelCmdProcs[MAXMIXERCOMMANDS] = {
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
Snack_MixerCmd(ClientData cdata, Tcl_Interp *interp, int objc,
	       Tcl_Obj *CONST objv[])
{
  int index;

  if (objc < 2) {
    Tcl_WrongNumArgs(interp, 1, objv, "option ?arg?");
    return TCL_ERROR;
  }
  
  if (Tcl_GetIndexFromObj(interp, objv[1], mixerCmdNames, "option", 0,
			  &index) != TCL_OK) {
    return TCL_ERROR;
  }

  return((mixerCmdProcs[index])(interp, objc, objv)); 
}

void
Snack_MixerDeleteCmd(ClientData clientData)
{
  int i;

  for (i = 0; i < nMixerCommands; i++) {
    if (mixerDelCmdProcs[i] != NULL) {
      (mixerDelCmdProcs[i])();
    }
  }
}
