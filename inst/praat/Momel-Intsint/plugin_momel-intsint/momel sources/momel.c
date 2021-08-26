/*
 * Copyright (C) 1997-1998 Centre National de la Recherche Scientifique
 * AUTHOR(S)  Robert ESPESSER, Daniel HIRST 
 * (Laboratoire Parole et Langage CNRS-ESA 6057)
 * Combination of cible.c, reduc.c and borne.c into Momel.c: 
 * Stefan Werner 2001 v0.1 22/4
 */
 
/*  Combine cible, reduc and borne into one program: Momel. */
/*  Arguments are: */
/*                 lfen1 hzinf hzsup maxec (for cible)   */
/*                 lfen2 seuildiff_x seuilrapp_y (for reduc) */
/*  stdin:         measured F0 values */
/*  stdout:        target values */

/****************************************************************************/
/*     cible $fen_cible $hzinf $hzsup 1.05 < $F0  \                         */
/*     	| reduc $fen_reduc $ecart_min_trame $rapport_min_freq  > $CIBLES    */
/*                                                                          */
/*     borne $F0 $CIBLES  $ecart_min_trame                                  */
/****************************************************************************/
/*     tst 300 50 600 1.04 200 50 0.05 < tstf0 */
/****************************************************************************/
/* momel 30 50 600  1.04 20 5 0.05
/* GPK 2007
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/*  #include "mt_config.h" */
#include <math.h>
#include <assert.h>
#include "momel.h"

/*  global vars */
/* int   nval=0; */
/*  float hz[NVALM]; */
/* float bhz[100 + NVALM + 100]; */
/* float *hzptr = bhz + 100; */

/* int   lfen1; */
/* float hzinf, hzsup, maxec; */
/* int   lfen2; */
/* float seuildiff_x = (ECART_MIN/PAS_TRAME) ; */
/* float seuilrapp_y = RAPP_MIN_FREQ ; */

/* int nred2; */
FILE *outf;

struct st_cib
{
   float     x;
   float     y;
};

struct st_cibred 
{
   float     x;
   float     y;
   int       p;
};



void eliminer_glitch(int nval, float *hzptr)
{
    int i;
    for(i=0; i<nval; i++)
    {
    	if((hzptr[i] > hzptr[i-1]*(1+RAPP_GLITCH))
		&& (hzptr[i] > hzptr[i+1]*(1+RAPP_GLITCH)))
	{
		hzptr[i] = 0.0;
	}
    }
}

/*****************************************************************************/
/*            cible ****************************** */
/*****************************************************************************/

int   calcrgp (	/* entree: */
		float *pond, int dpx, int fpx,
		float *hzptr,
		/* Retour: */
		float *pa0, float *pa1, float *pa2, float *hzes)
/* valeur calcrgp 0:ras        1:Pb calcul */
{
   double   muet;
   double    pn;
   double   sx,
            sx2,
            sx3,
            sx4,
            sy,
            sxy,
            sx2y;
   double   spdxy,
            spdx2,
            spdx3,
            spdx4,
            spdx2y;
   int   ix;

   pn = 0.;
   sx = sx2 = sx3 = sx4 = sy = sxy = sx2y = 0.;

   for (ix = dpx; ix <= fpx; ix++)
   {
      const double p = pond[ix];
      if (p != 0.)
      {
	 const double val_ix = (double) ix;
	 const double y = hzptr[ix];
	 const double x2 = val_ix * val_ix;
	 const double x3 = x2 * val_ix;
	 const double x4 = x2 * x2;
	 const double xy = val_ix * y;
	 const double x2y = x2 * y;

	 pn += p;
	 sx += (p * val_ix);
	 sx2 += p * x2;
	 sx3 += p * x3;
	 sx4 += p * x4;
	 sy += p * y;
	 sxy += p * xy;
	 sx2y += p * x2y;
      }
   }
   if (pn < 3.)
      return (1);

   spdxy = sxy - (sx * sy) / pn;
   spdx2 = sx2 - (sx * sx) / pn;
   spdx3 = sx3 - (sx * sx2) / pn;
   spdx4 = sx4 - (sx2 * sx2) / pn;
   spdx2y = sx2y - (sx2 * sy) / pn;

   muet = spdx2 * spdx4 - spdx3 * spdx3;
   if (spdx2 == 0.|| muet == 0.)
      return (1);

/*fprintf(outf,"avant calcul pa210 muet %g  spdx2 %g pn %g\n",muet,spdx2,pn);*/

   *pa2 = (spdx2y * spdx2 - spdxy * spdx3) / muet;
   *pa1 = (spdxy - *pa2 * spdx3) / spdx2;
   *pa0 = (sy - *pa1 * sx - *pa2 * sx2) / pn;

   for (ix = dpx; ix <= fpx; ix++)
      hzes[ix] = *pa0 + (*pa1 + *pa2 * (float) ix) * (float) ix;

   return (0);
}



int cible(int nval, float *hzptr, int lfen1,
		float maxec, float hzinf, float hzsup,
		struct st_cib *cib)
{
/*     static float bhz[100 + NVALM + 100]; */
   static float  bpond[100 + NVALM + 100];
   static float bhzes[100 + NVALM + 100];

   float   *hzes;
   float *pond;

   float    a0,
            a1,
            a2;
   float    xc,
            yc;
   float    vxc,
            vyc;

   int   ix,
         x;
   int   dpx,
         fpx;
   int   nsup,
         nsupr;
   int   lfens2;
   int   retour_rgp;
   int verbeux = 0;

   assert(hzsup > hzinf);
   for (ix = 0; ix < NVALM + 200; ix++)
   {
      bpond[ix] = 0.0;
   }

/*     hzptr = bhz + 100; */
   pond = bpond + 100;
   hzes = bhzes + 100;
   

   assert(nval >= 0);
   for (ix = 0; ix < nval; ix++)
   {
      if (hzptr[ix] > SEUILV)
	 pond[ix] = 1.0;
   }


   assert(lfen1 > 0);
   lfens2 = lfen1 / 2;

   for (ix = 0; ix < nval; ix++)
   {
      dpx = ix - lfens2;
      fpx = dpx + lfen1;
      nsup = 0;
      nsupr = -1;

	  if(ix==24)
		  printf("");
      {
/*  SW kludge */
	 float  bpondloc[1000 + NVALM + 1000];
/*  SW kludge */
	 float *pondloc;
	 int   i;

	 pondloc = bpondloc + 1000;
 /* recopie */
	 for (i = dpx; i <= fpx; i++)
	   {
/*  	     printf("++++++++++++++++++++\n\t%d\n", i); */
	     pondloc[i] = pond[i];
/*  	     printf("\t%d\t%d\t%d\n", i, pond[i], dpx); */
	   }

	 while (nsup > nsupr)
	 {
	    nsupr = nsup;
	    nsup = 0;
	    retour_rgp = calcrgp (pondloc, dpx, fpx, hzptr, &a0, &a1, &a2, hzes);

/************************** SW **********************************************/
/*  	    printf("*cible*  %g\n",*hzes); */
/************************** SW **********************************************/
	    if (retour_rgp != 0)
	       break;		       /* quitte le while */

	    for (x = dpx; x <= fpx; x++)
	    {
	       if (hzptr[x] == 0.|| hzes[x] / hzptr[x] > maxec)
	       {
		  pondloc[x] = 0;
		  nsup++;
	       }
	    }
	 }			       /* while */

      }				       /* bloc */

      xc = yc = 0.;
      if (retour_rgp == 0 && a2 != 0.)
      {				       /* tout va bien */
	 vxc = -a1 / (a2 + a2);
	 if ((vxc > ix - lfen1) && (vxc < ix + lfen1))
	 {
	    vyc = a0 + (a1 + a2 * vxc) * vxc;
	    if (vyc > hzinf && vyc < hzsup)
	    {
	       xc = vxc;
	       yc = vyc;
	    }
	 }
      }

	  if(verbeux){
		  fprintf(outf,"%f ", xc);
		  fprintf(outf,"%f\n", yc);
	  }
//      	fprintf (outf, "ret %d a0 %g a1 %g a2 %g\n", retour_rgp, a0, a1, a2);
//      	fprintf (outf, "%g %g %g\n", a0, a1, a2);
 /* ecriture ptcible pour chaque hz */
/*        fwrite (&xc, sizeof (float), 1, stdout); */
/*        printf ("%g\t", xc); */
/*        fwrite (&yc, sizeof (float), 1, stdout); */
/*        printf ("%f\n", yc); */

/************************** SW **********************************************/
      cib[ix].x = xc; 
      cib[ix].y = yc;
/*        printf("%g\t", cib[ix].x); */
/*        printf("%g\n", cib[ix].y); */
/************************** SW **********************************************/
   }				       /* for ix */

   return (0);
}				       /* main */

/*******************
 */




/*****************************************************************************/
/*            reduc ****************************** */
/*****************************************************************************/

   
/****************
 * cle pour qsort :ordre temporel croissant 
 */
int cb_compare(const void *va, const void *vb)
{
/* vu le proto de qsort :
 extern void qsort(void *, size_t, size_t,
	int (*)(const void *, const void *));

 const : interdit de modifier ce qui est pointe par;
 le cast , sinon warning
*/

struct st_cib *pa = (struct st_cib *) va;
struct st_cib *pb = (struct st_cib *) vb;

                  if (pa->x > pb->x)
                          return (1);
                  if (pa->x < pb->x  )
                          return (-1);
                  return (0);
}




/***************************/
int reduc(int nval, int lfen2,
		float seuildiff_x, float seuilrapp_y,
		/* modifie: */
		struct st_cib *cib,
		/* retour: */
		int *pnred2,
		struct st_cibred *cibred2
		)
{
/*     struct st_cib cib[NVALM] ; */
   struct st_cibred cibred[PARMAX] ;
/*     struct st_cibred cibred2[PARMAX] ; */
   
   struct st_cibred  cibred_cour ;

   static float xdist[NVALM];
   static float ydist[NVALM];
   static float dist[NVALM];
   static int xd[PARMAX + 1];

   int       i;
   int       j;
   int       j1, j2;
   int       ng, nd;
   int       ncible;
   int       n, np, xmax;
   int       susseuil;


   float     px, py;
   float     xds, yds;
   float     sxg, syg, sxd, syd;

   float     seuil;

   int       lf;
   int       verbeux = 0;

      
/************************** SW **********************************************/
/*     for (i=0; i<nval; i++) { */
/*        printf("*reduc* cib[%d].x %g\t", i, cib[i].x); */
/*        printf("*reduc* cib[%d].y %g\n", i, cib[i].y); */
/*     } */

/*  for(i=0; i < nval; i++) */
/*     printf ("%g %g\n", cib[i].x * PAS_TRAME, cib[i].y); */

/************************** SW **********************************************/
   

   lf = lfen2 / 2;
   xds = yds = 0.;
   np = 0;
   for (i = 0; i < nval - 1; i++)
   {
      j1 = 0;
      if (i > lf)
	 j1 = i - lf;
      j2 = nval - 1;
      if (i + lf < nval - 1)
	 j2 = i + lf;
      sxg = syg = 0.;
      ng = 0;
      for (j = j1; j < i + 1; j++)			   /* a gauche */
      {
	 if (cib[j].y > SEUILV)
	 {
	    sxg += cib[j].x;
	    syg += cib[j].y;
	    ng++;
	 }
      }
      sxd = syd = 0.;
      nd = 0;
      for (j = i + 1; j < j2; j++)
      {
	 if (cib[j].y > SEUILV)
	 {
	    sxd += cib[j].x;
	    syd += cib[j].y;
	    nd++;
	 }
      }

      xdist[i] = ydist[i] = -1.;
      if (nd * ng > 0)
      {
	 xdist[i] = fabs (sxg / ng - sxd / nd);
	 ydist[i] = fabs (syg / ng - syd / nd);
	 xds += xdist[i];
	 yds += ydist[i];
	 np++;
      }
   }							   /* boucle i */

/* on pondere par 1/ distance moyenne
 */

   px = np / xds;
   py = np / yds;

   for (i = 0; i < nval; i++)
   {
      dist[i] = -1;
      if (xdist[i] > 0.)
      {
	 dist[i] = (xdist[i] * px + ydist[i] * py) / (px + py);
      }
   }

/* seuil = moy des dist ponderees
 */

   seuil = 2. / (px + py);

   xd[0] = 0;
   ncible = 0;
   susseuil = FAUX;

/*cherche les maxs des pics de dist > seuil
 */

   for (i = 0; i < nval; i++)
   {
      if (ncible == PARMAX)
      {
	 fprintf (outf, "reduc: trop de partitions (max: %d)\n", PARMAX);
	 return (1);
      }

      if (susseuil == FAUX)
      {
	 if (dist[i] > seuil)
	 {
	    susseuil = VRAI;
	    xmax = i;
	 }
      }
      if (susseuil == VRAI)
      {
	 if (dist[i] > dist[xmax])
	 {
	    xmax = i;
	 }
	 if (dist[i] < seuil)
	 {
	    ncible++;
	    xd[ncible] = xmax;
	    susseuil = FAUX;
	 }
      }
   }							   /* boucle i */

   if (susseuil == VRAI)
   {
      ncible++;
      xd[ncible] = xmax;
   }

   if (xmax < nval)
   {
      ncible++;
      xd[ncible] = nval;
   }

/* ncibles est bien le nbre de partitions
 * on [0 a ncible] debut de partition xd[i]
 */

   if (verbeux)
   {
      int       k;

      fprintf (outf, " ncible %d partitions\n", ncible);
      for (k = 0; k < ncible; k++)
	 fprintf (outf, "k %d debut %d a %d\n", k, xd[k], xd[k + 1]);
   }


/**********************
 *partition sur les x
 */
   {
      int       ip, j;
      int       k, ncibr;
      float     sx, sx2, sy, sy2;
      double    varx, vary;
      double    xm, ym;
      float     et2x, et2y;
      float     seuilbx, seuilby, seuilhx, seuilhy;
      int       parinf, parsup;

      int ncibr_brut ;
/*        int nred2; */
      
/* ncibr pointe sur la cible reduite precedente
 * on termine avec des items indices [0 a ncibr]
 */
      ncibr = -1;
      for (ip = 0; ip < ncible; ip++)			   /* sur les partitions */
      {
	 parinf = xd[ip];				   /* bornes partition courante */
	 parsup = xd[ip + 1];
	 for (k = 0; k < 1; k++)			   /* 1 (!) passages */
	 {
	    sx = sx2 = sy = sy2 = 0.;
	    n = 0;

	    if (verbeux)
	       fprintf (outf,
			"parti %d k %d   debut %d fin %d\n", ip, k, xd[ip], xd[ip + 1]);
	    /* moyenne sigma */
	    for (j = parinf; j < parsup; j++)		   /* surla pop dune parttion */
	    {
	       if (cib[j].y > 0.)
	       {
/*fprintf(outf,"cible %d x %g y %g\n",j,cib[j].x ,cib[j].y);*/
		  sx += cib[j].x;
		  sx2 += cib[j].x * cib[j].x;
		  sy += cib[j].y;
		  sy2 += cib[j].y * cib[j].y;
		  n++;
	       }
	    }						   /* boucle j */

	    if (n > 1)					   /* pour la variance */
	    {
	       xm = (double) sx / (double) n;
	       ym = (double) sy / (double) n;
	       varx = (double) sx2 / (double) n - xm * xm;
	       vary = (double) sy2 / (double) n - ym * ym;
/*fprintf(outf," X k%d n %d %g %g %g \n", k,n,(double)sx2/(double)n,xm*xm,varx);
 fprintf(outf,"Y k%d n %d  %g %g %g\n", k,n,(double)sy2/(double)n,ym*ym, vary);
 */
	       if (varx <= 0.)				   /* cas ou variance devarit etre == +epsilon */
		  varx = 0.1;
	       if (vary <= 0.)
		  vary = 0.1;

	       et2x = FSIGMA * sqrt (varx);
	       et2y = FSIGMA * sqrt (vary);
	       seuilbx = xm - et2x;
	       seuilhx = xm + et2x;
	       seuilby = ym - et2y;
	       seuilhy = ym + et2y;

	       for (j = parinf; j < parsup; j++)	   /* elimination */
	       {
		  if (cib[j].y > 0. &&
		      (cib[j].x < seuilbx || cib[j].x > seuilhx
		       || cib[j].y < seuilby || cib[j].y > seuilhy))
		  {
		     cib[j].x = cib[j].y = 0.;
		  }
	       }					   /* boucle j */

	    }						   /* if n>1 */
	 }						   /* boucle k */

/*recalcule moyennes */
	 sx = sy = 0.;
	 n = 0;
	 for (j = parinf; j < parsup; j++)
	 {
	    if (cib[j].y > 0.)
	    {
	       sx += cib[j].x;
	       sy += cib[j].y;
	       n++;
	    }
	 }
	 if (n > 0)
	 {
	    cibred_cour.x = sx / n;
	    cibred_cour.y = sy / n;
	    cibred_cour.p = n;

	    if (verbeux)
	       fprintf (outf,
			"ncibr %d  x %g  y %g    p %d\n", ncibr,
			cibred_cour.x, cibred_cour.y, cibred_cour.p
		  );
	    if (ncibr < 0)
	    {
	       ncibr++;
	       cibred[ncibr] = cibred_cour;
	    }
	    else
	    {
 /*  si les cibred[].x ne sont pas strictement croissants 
   on ecrase  la cible ayant le poids le moins fort 
   */

	       if (cibred_cour.x > cibred[ncibr].x)/* 1 cibred en +  car t croissant */
	       {
		  ncibr++;
		  cibred[ncibr] = cibred_cour;
	       }
	       else		       /* t <= precedent */
	       {
	          if(verbeux)
	             fprintf(outf,"reduc:cb arriere %g %d    %g %d\n",
	             cibred_cour.x, cibred_cour.p, 
	             cibred[ncibr].x, cibred[ncibr].p);
	             
		  if (cibred_cour.p > cibred[ncibr].p)
		  {		       /* si p courant >, ecrase la precedente */
		     cibred[ncibr] = cibred_cour;
		  }
	       }
	    }		  

         }  /*if n>0*/

      }							   /* boucle ip */

ncibr_brut= ncibr+1;  /* ( cibred[0,...ncibr_brut[  */

/* 2eme filtrage  des cb trop proches en t [et Hz]*/

/* classe ordre temporel croissant les cibred*/
qsort( (char *) cibred, ncibr_brut, sizeof(struct st_cibred), cb_compare);

*pnred2 = 0;
cibred2[0] = cibred[0];

      assert(seuilrapp_y > 0.0);
      for (i = 1; i < ncibr_brut; i++)
      {
         if (cibred[i].x - cibred2[*pnred2].x < seuildiff_x)
         {
               
            if ( fabs( (double)(cibred[i].y - cibred2[*pnred2].y) )/ cibred2[*pnred2].y  < seuilrapp_y )
            {
               if(verbeux)
                  fprintf(outf,"reduc:fusion_moyenne %g %g %d   %g %g %d\n",
                  cibred2[*pnred2].x,cibred2[*pnred2].y,cibred2[*pnred2].p,
                  cibred[i].x, cibred[i].y, cibred[i].p 
                  );
               cibred2[*pnred2].x = (cibred2[*pnred2].x +cibred[i].x)/2. ;
               cibred2[*pnred2].y = (cibred2[*pnred2].y +cibred[i].y)/2. ;
               cibred2[*pnred2].p += cibred[i].p  ; 
            }
            else
            {
               if(verbeux)
                  fprintf(outf,"reduc:choix_poids  %g %d   %g %d\n",
                  cibred2[*pnred2].x, cibred2[*pnred2].p,
                  cibred[i].x, cibred[i].p
                  ); 
                  
               if (cibred2[*pnred2].p < cibred[i].p)
               {
                  cibred2[*pnred2] = cibred[i] ;
               }
            }
         }
         else
         {            
            (*pnred2)++ ;
            cibred2[*pnred2] = cibred[i] ;
         } 
      }


(*pnred2)++ ; /*le nbre de cible reduite*/
if (verbeux)
   fprintf(outf,"reduc: ncibred_brut: %d  ncibred2 %d\n", ncibr_brut,*pnred2);

/************************** SW **********************************************/
/*  for(i=0; i < *pnred2; i++) */
/*     printf ("*reduc*  %g %g\n", cibred2[i].x * PAS_TRAME, cibred2[i].y); */
/************************** SW **********************************************/

   }							   /* du bloc partition */

return(0);
}							   /* main */


/*****************************************************************************/
/*            borne ****************************** */
/*****************************************************************************/

int borne(int nval, int nred2,
		const struct st_cibred *cibred2,
		const float *hzptr
		)
{
/************************** SW **********************************************/
struct st_cibred ancre, borne;
/*  int ncib; */
/************************** SW **********************************************/

int i, j;
int dernier_voise, premier_voise;
float frontiere;
float sx2y, sx4, x2 ;
float a;
/*  float    seuildiff_x = (ECART_MIN/PAS_TRAME) ; */
int halo = HALO_BORNE_TRAME ;
int verbeux = 0;

/************************** SW **********************************************/
/************************** SW **********************************************/

/******** principes
* calcul borne G (D)  si 1ere (derniere) cible est
*   >   (debut_voisement+halo)
* ( < (fin_voisement -halo)  )
*
* ce pt de debut(fin) voisement  == frontiere
* cible extremite == ancre
* regression quadratique sur Hz de  [frontiere ancre]
*/

/*borne gauche*/
/*recherche 1er voise*/
i=0;
while(i <nval && hzptr[i] < SEUILV)
   i++ ;

premier_voise = i;   

if( (int)cibred2[0].x > premier_voise + halo )
{
    ancre = cibred2[0] ;
 /*origine des t: ancre.x, des y : ancre.y*/
  
  sx2y = 0. ;
  sx4 = 0. ;
  
  j = 0;
  for(i = (int)ancre.x  ; i >= 0; i--)
  {
     if (hzptr[i] > SEUILV)
     {
     x2 = (float)j*(float)j ;
     sx2y += x2* (hzptr[i] - ancre.y) ;
     sx4 += x2*x2 ;
     }
     j++ ;
  }

  frontiere = (float)premier_voise ;

  a = sx2y/sx4 ;
  
  borne.x = frontiere - (ancre.x - frontiere ) ;
  borne.y = ancre.y + 2*a * (ancre.x - frontiere)*(ancre.x - frontiere);

/*    fprintf(fluxcib,"%g %g\n", borne.x*PAS_TRAME, borne.y);     */
/************************** SW **********************************************/
  printf("%g %g\n", borne.x*PAS_TRAME, borne.y);
/************************** SW **********************************************/

  if(verbeux)
     fprintf(outf,"borne:bornegauche %g %g   ancre %g %g  frontiere %g\n",
    borne.x, borne.y, ancre.x, ancre.y ,frontiere );
}


/*reecrit les initiaux*/
for(i=0; i <nred2; i++)

/*borne droite eventuelle*/

/* recherche dernier voisement*/
i = nval-1 ;
while ( i >=0 && hzptr[i] < SEUILV)
   i--;

dernier_voise = i;

/************************** SW **********************************************/
for(i=0; i < nred2; i++)
   printf ("%g %g\n", cibred2[i].x * PAS_TRAME, cibred2[i].y);
/************************** SW **********************************************/

if( (int)cibred2[nred2-1].x < dernier_voise - halo)
{
    ancre = cibred2[nred2-1] ;

 /*origine des t: ancre.x, des y : ancre.y*/
  
  sx2y = 0. ;
  sx4 = 0. ;
  
  j = 0;
  for( i = (int)ancre.x ; i < nval; i++)
  {
     if(hzptr[i] > SEUILV)
     {
     x2 = (float)j*(float)j ;
     sx2y += (x2 * (hzptr[i] - ancre.y)) ;
     sx4 += x2*x2 ;
     }
     j++ ;
  }

  frontiere = (float)dernier_voise ;

  a = sx2y/sx4 ;
  
  borne.x = frontiere + (  frontiere -ancre.x) ;
  borne.y = ancre.y + 2.* a * (ancre.x - frontiere)*(ancre.x - frontiere);

if(verbeux)
   fprintf(outf,"borne: bornedroit %g %g   ancre %g %g   frontiere %g\n",
     borne.x,borne.y, ancre.x, ancre.y ,frontiere);
   
   printf("%g %g\n", borne.x*PAS_TRAME, borne.y);
}

return(0);
}


/*  *********************************************************************** */
/*  *********************************************************************** */

//30 50 600  1.04 20 5 0.05  test2.dat test2.txt
int main(int argc, char **argv)
{

  int i, err=0;
  int nval = 0;
//  float bhz[100 + NVALM + 100];
  float *bhz;
  float *hzptr; 
  int lfen1, lfen2;
  int nred2 = 0;
  float seuilrapp_y = 0.0, seuildiff_x = 0.0;
  float hzinf = 0.0, hzsup = 0.0, maxec = 0.0;
//  struct st_cib cib[NVALM];
//  struct st_cibred cibred2[PARMAX];
  struct st_cib *cib;
  struct st_cibred *cibred2;
  int elim_glitch = VRAI;

  FILE *inf = fopen(argv[8], "r");
  outf = fopen("tmp.txt", "w");

  bhz=(float*)malloc(sizeof(float)*(100 + NVALM + 100));
  cib = (struct st_cib *)malloc(sizeof(struct st_cib)*NVALM);
  cibred2 = (struct st_cibred *)malloc(sizeof(struct st_cibred)*PARMAX);
  hzptr = bhz + 100;


  argc--;
  if(argc < 7 ||strcmp(argv[1],"--usage") ==0 ||strcmp(argv[1],"-h") ==0)
    {
      fprintf(outf,"\nusage: \tMtargets [options] win1 lo hi maxerr win2 mind minr < [F0 values]\n\n");
      fprintf(outf,"win1:\t\tcible window length\n");
      fprintf(outf,"lo:\t\tF0 threshold\n");
      fprintf(outf,"hi:\t\tF0 ceiling\n");
      fprintf(outf,"maxerr:\t\tmaximum error\n");
      fprintf(outf,"win2:\t\treduc window length\n");
      fprintf(outf,"mind:\t\tminimal distance\n");
      fprintf(outf,"minr:\t\tminimal frequency ratio\n");
      fprintf(outf,"\noptions:\t\t--non-elim-glitch\n");
      fprintf(outf,"\nF0 values, sampled at 10 ms, in ASCII file.\n");
      fprintf(outf,"\n");
      return 1;
    }
  if(argc>1 && strcmp(argv[1], "--non-elim-glitch")==0)
  {
  	elim_glitch = FAUX;
	argc--, argv++;
  }

  lfen1       = atoi (*++argv);
  //* lfen1 est en pointes echantillon, pas milliseconds. *
  hzinf       = atof (*++argv);
  hzsup       = atof (*++argv);
  //* En Hirst, Di Cristo et Espesser, hzsup est calcule automatiquement. *
  maxec       = atof (*++argv);
  //* Maxec est 1+Delta en Hirst, Di Cristo et Espesser. *
  lfen2       = atoi (*++argv);
  //* lfen2 est en pointes echantillon, pas milliseconds. *
  seuildiff_x = atof (*++argv);
  //* seuildiff_x est en pointes echantillon, pas milliseconds. 
  seuilrapp_y = atof (*++argv);

//*    printf(" lfen1: %d\n hzinf: %g\n hzsup: %g\n maxec: %g\n lfen2: %d\n seuildiff_x: %g\n seuilrapp_y: %g\n",  lfen1, hzinf, hzsup, maxec, lfen2,  seuildiff_x, seuilrapp_y); 

//*  from cible 
  for (i = 0; i < NVALM + 200; i++)
    {
      bhz[i] = 0.;
    }
  
//  while (! feof(stdin) && nval < NVALM)
	while (! feof(inf) && nval < NVALM)
    {
      //*        scanf ("%f", &(hz[i]));
//      scanf ("%f", hzptr+nval);
      fscanf(inf,"%f\n", hzptr+nval);
      nval++;
    }
	fclose(inf);

  if(nval == NVALM)
    {
      fprintf(outf,"Too many F0 values (max: %d)\n" ,NVALM);
      return 1;
    }

  if(elim_glitch)
  	eliminer_glitch(nval, hzptr);
  
  err = cible(nval, hzptr, lfen1,
		maxec, hzinf, hzsup,
		cib);
 if(err)
 	return err;

  err = reduc(nval, lfen2,
		seuildiff_x, seuilrapp_y,
		//* modifie:
		cib,
		//* retour: 
		&nred2,
		cibred2
		);
 if(err)
 	return err;

  err = borne(nval, nred2,
		cibred2,
		hzptr
		);

  fclose(outf);
  return err;
}
