/* 
 * Copyright (C) 1997-2003 Kare Sjolander <kare@speech.kth.se>
 *
 * This file contains a sample module which demonstrates how to write
 * extensions to the Snack sound extension for Tcl/Tk.
 * The latest version of Snack can be found at http://www.speech.kth.se/snack/
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

#include "snack.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Windows magic */

#if defined(__WIN32__)
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#undef WIN32_LEAN_AND_MEAN
#define EXPORT(a,b) __declspec(dllexport) a b
BOOL APIENTRY
DllMain(HINSTANCE hInst, DWORD reason, LPVOID reserved)
{
  return TRUE;
}
#else
#define EXPORT(a,b) a b
#endif

int
Square(ClientData cdata, Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[])
{
  Sound *sound;
  int i;

  /* Get the sound structure for this sound */
      
  sound = Snack_GetSound(interp, Tcl_GetStringFromObj(objv[0], NULL));

  /* create a simple square wave */

  for (i = 0; i < Snack_GetLength(sound); i++) {
    if ((i/10)%2) {
      Snack_SetSample(sound, 0, i, 10000.0f);
    } else {
      Snack_SetSample(sound, 0, i, -10000.0f);
    }
  }

  /* update the max/min members of the sound structure */

  Snack_UpdateExtremes(sound, 0, sound->length, SNACK_NEW_SOUND);

  /* execute callbacks for stuff like canvas items */

  Snack_ExecCallbacks(sound, SNACK_NEW_SOUND);

  return TCL_OK;
}

/*
  Initialize the square package and create a new sound command 'square'.
  The syntax is: sndName square
 */

EXPORT(int, Square_Init)(Tcl_Interp *interp)
{
#ifdef USE_TCL_STUBS
  if (Tcl_InitStubs(interp, "8", 0) == NULL) {
    return TCL_ERROR;
  }
#endif

#ifdef USE_TK_STUBS
    if (Tk_InitStubs(interp, "8", 0) == NULL) {
      return TCL_ERROR;
    }
#endif

#ifdef USE_SNACK_STUBS
  if (Snack_InitStubs(interp, "2", 0) == NULL) {
    return TCL_ERROR;
  }
#endif

  if (Tcl_PkgProvide(interp, "square", "1.0") != TCL_OK) {
    return TCL_ERROR;
  }

  Snack_AddSubCmd(SNACK_SOUND_CMD, "square", (Snack_CmdProc *) Square, NULL);

  return TCL_OK;
}

EXPORT(int, Square_SafeInit)(Tcl_Interp *interp)
{
  return Square_Init(interp);
}

#ifdef __cplusplus
}
#endif
