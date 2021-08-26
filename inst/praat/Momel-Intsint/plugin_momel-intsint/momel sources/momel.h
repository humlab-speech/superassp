/* @(#)momel.h
 */

#ifndef _MOMEL_H
#define _MOMEL_H 1

/*  cible.h */
#define PAS 10.   /* en ms */
#define NVALM 100000   /* 15 mn == 90000 avec pas de 10 ms  */
#define SEUILV 50.
/*  cible.h */

#define PARMAX (NVALM/2)

#define FAUX 0
#define VRAI 1
#define FSIGMA 1.

#define PAS_TRAME 10.
#define ECART_MIN 50.   /*en ms*/
#define RAPP_MIN_FREQ 0.05
#define RAPP_GLITCH 0.05
#define HALO_BORNE_TRAME 4

FILE *outf;

#endif /* _MOMEL_H */

  