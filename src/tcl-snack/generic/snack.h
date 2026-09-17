/* 
 * Copyright (C) 1997-2005 Kare Sjolander <kare@speech.kth.se>
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

#include "tcl.h"

#ifdef BUILD_snack
# undef TCL_STORAGE_CLASS
# define TCL_STORAGE_CLASS DLLEXPORT
#endif

#ifndef CONST84
#   define CONST84 const
#endif

#ifndef CONST86
#   define CONST86 const
#endif

#ifndef CONST
#   define CONST const
#endif

#define _ANSI_ARGS_(x) x

#ifndef Tk_Offset
#   define Tk_Offset(x,y) offsetof(x,y)
#endif

#ifndef TCL_SIZE_MAX
#   define Tcl_Size int
#endif

#include "jkSound.h"
#include "jkAudIO.h"

#define SNACK_VERSION     "2.2"
#define SNACK_PATCH_LEVEL "2.2.10"
