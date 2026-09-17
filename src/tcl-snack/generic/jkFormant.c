/* formant.c */
/*
 * This software has been licensed to the Centre of Speech Technology, KTH
 * by AT&T Corp. and Microsoft Corp. with the terms in the accompanying
 * file BSD.txt, which is a BSD style license.
 *
 *    "Copyright (c) 1987-1990  AT&T, Inc.
 *    "Copyright (c) 1986-1990  Entropic Speech, Inc. 
 *    "Copyright (c) 1990-1994  Entropic Research Laboratory, Inc. 
 *                   All rights reserved"
 *
 * Written by:  David Talkin
 * Revised by: John Shore
 *
 */

#include "snack.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "jkFormant.h"


int debug = 0;
int w_verbose = 0;

/*	dpform.c       */

/* a formant tracker based on LPC polynomial roots and dynamic programming */
				/***/
/* At each frame, the LPC poles are ordered by increasing frequency.  All
   "reasonable" mappings of the poles to F1, F2, ... are performed.
   The cost of "connecting" each of these mappings with each of the mappings
   in the previous frame is computed.  The lowest cost connection is then
   chosen as the optimum one.  At each frame, each mapping has associated
   with it a cost based on the formant bandwidths and frequencies.  This
   "local" cost is finally added to the cost of the best "connection."  At
   end of utterance (or after a reasonable delay like .5sec) the best
   mappings for the entire utterance may be found by retracing back through
   best candidate mappings, starting at end of utterance (or current frame).
*/

/* Here are the major fudge factors for tweaking the formant tracker. */
#define MAXCAN	300  /* maximum number of candidate mappings allowed */
static double MISSING = 1, /* equivalent delta-Hz cost for missing formant */
	NOBAND = 1000, /* equivalent bandwidth cost of a missing formant */
	DF_FACT =  20.0, /* cost for proportional frequency changes */
	/* with good "stationarity" function:*/
  /*        DF_FACT =  80.0, *//*  cost for proportional frequency changes */
	DFN_FACT = 0.3, /* cost for proportional dev. from nominal freqs. */
	BAND_FACT = .002, /* cost per Hz of bandwidth in the poles */
/*	F_BIAS	  = 0.0004,   bias toward selecting low-freq. poles */
	F_BIAS	  = 0.000, /*  bias toward selecting low-freq. poles */
	F_MERGE = 2000.0; /* cost of mapping f1 and f2 to same frequency */
static double	*fre,
		fnom[]  = {  500, 1500, 2500, 3500, 4500, 5500, 6500},/*  "nominal" freqs.*/
		fmins[] = {   50,  400, 1000, 2000, 2000, 3000, 3000}, /* frequency bounds */
		fmaxs[] = { 1500, 3500, 4500, 5000, 6000, 6000, 8000}; /* for 1st 5 formants */

static int	maxp,	/* number of poles to consider */
		maxf,	/* number of formants to find */
		ncan,  domerge = TRUE;

static short **pc;

static int canbe(int pnumb, int fnumb) /* can this pole be this freq.? */
{
return((fre[pnumb] >= fmins[fnumb])&&(fre[pnumb] <= fmaxs[fnumb]));
}

/* This does the real work of mapping frequencies to formants. */
static void candy(int cand, int pnumb, int fnumb)
     /* cand:  candidate number being considered */
     /* pnumb: pole number under consideration */
     /* fnumb: formant number under consideration */
{
  int i,j;

  if(fnumb < maxf) pc[cand][fnumb] = -1;
  if((pnumb < maxp)&&(fnumb < maxf)){
    /*   printf("\ncan:%3d  pnumb:%3d  fnumb:%3d",cand,pnumb,fnumb); */
    if(canbe(pnumb,fnumb)){
      pc[cand][fnumb] = pnumb;
      if(domerge&&(fnumb==0)&&(canbe(pnumb,fnumb+1))){ /* allow for f1,f2 merger */
	ncan++;
	pc[ncan][0] = pc[cand][0];
	candy(ncan,pnumb,fnumb+1); /* same pole, next formant */
      }
      candy(cand,pnumb+1,fnumb+1); /* next formant; next pole */
      if(((pnumb+1) < maxp) && canbe(pnumb+1,fnumb)){
	/* try other frequencies for this formant */
	ncan++;			/* add one to the candidate index/tally */
	/*		printf("\n%4d  %4d  %4d",ncan,pnumb+1,fnumb); */
	for(i=0; i<fnumb; i++)	/* clone the lower formants */
	  pc[ncan][i] = pc[cand][i];
	candy(ncan,pnumb+1,fnumb);
      }
    } else {
      candy(cand,pnumb+1,fnumb);
    }
  }
  /* If all pole frequencies have been examined without finding one which
     will map onto the current formant, go on to the next formant leaving the
     current formant null. */
  if((pnumb >= maxp) && (fnumb < maxf-1) && (pc[cand][fnumb] < 0)){
    if(fnumb){
      j=fnumb-1;
      while((j>0) && pc[cand][j] < 0) j--;
      i = ((j=pc[cand][j]) >= 0)? j : 0;
    } else i = 0;
    candy(cand,i,fnumb+1);
  }
}

/* Given a set of pole frequencies and allowable formant frequencies
   for nform formants, calculate all possible mappings of pole frequencies
   to formants, including, possibly, mappings with missing formants. */
void get_fcand(int npole, double *freq, double *band, int nform, short **pcan)
     /* poles ordered by increasing FREQUENCY */
{	

  ncan = 0;
  pc = pcan;
  fre = freq;
  maxp = npole;
  maxf = nform;
  candy(ncan, 0, 0);
  ncan++;	/* (converts ncan as an index to ncan as a candidate count) */
}

void set_nominal_freqs(double f1)
{
  int i;
  for(i=0; i < MAXFORMANTS; i++) {
    fnom[i] = ((i * 2) + 1) * f1;
    fmins[i] = fnom[i] - ((i+1) * f1) + 50.0;
    fmaxs[i] = fnom[i] + (i * f1) + 1000.0;
  }
}

/*      ----------------------------------------------------------      */
/* find the maximum in the "stationarity" function (stored in rms) */
double get_stat_max(register POLE **pole, register int nframes)
{
  register int i;
  register double amax, t;

  for(i=1, amax = (*pole++)->rms; i++ < nframes; )
    if((t = (*pole++)->rms) > amax) amax = t;

  return(amax);
}

Sound *dpform(Sound *ps, int nform, double nom_f1)
{
  double pferr, conerr, minerr, dffact, ftemp, berr, ferr, bfact, ffact,
         rmsmax, fbias, **fr, **ba, rmsdffact, merger=0.0, merge_cost,
         FBIAS;
  register int	i, j, k, l, ic, ip, mincan=0;
  short	**pcan;
  FORM	**fl;
  POLE	**pole; /* raw LPC pole data structure array */
  Sound *fbs;
  int dmaxc,dminc,dcountc,dcountf;
  
  if(ps) {
    if(nom_f1 > 0.0)
      set_nominal_freqs(nom_f1);
    pole = (POLE**)ps->extHead;
    rmsmax = get_stat_max(pole, ps->length);
    FBIAS = F_BIAS /(.01 * ps->samprate);
    /* Setup working values of the cost weights. */
    dffact = (DF_FACT * .01) * ps->samprate; /* keep dffact scaled to frame rate */
    bfact = BAND_FACT /(.01 * ps->samprate);
    ffact = DFN_FACT /(.01 * ps->samprate);
    merge_cost = F_MERGE;
    if(merge_cost > 1000.0) domerge = FALSE;

    /* Allocate space for the formant and bandwidth arrays to be passed back. */
    if(debug & DEB_ENTRY){
      printf("Allocating formant and bandwidth arrays in dpform()\n");
    }
    fr = (double**)ckalloc(sizeof(double*) * nform * 2);
    ba = fr + nform;
    for(i=0;i < nform*2; i++){
      fr[i] = (double*)ckalloc(sizeof(double) * ps->length);
    }
    /*    cp = new_ext(ps->name,"fb");*/
    /*    if((fbs=new_signal(cp,SIG_UNKNOWN,dup_header(ps->header),fr,ps->length,		       ps->samprate, nform * 2))) {*/
    if (1) {
      /* Allocate space for the raw candidate array. */
      if(debug & DEB_ENTRY){
	printf("Allocating raw candidate array in dpform()\n");
      }
      pcan = (short**)ckalloc(sizeof(short*) * MAXCAN);
      for(i=0;i<MAXCAN;i++) pcan[i] = (short*)ckalloc(sizeof(short) * nform);

      /* Allocate space for the dp lattice */
      if(debug & DEB_ENTRY){
	printf("Allocating DP lattice structure in dpform()\n");
      }
      fl = (FORM**)ckalloc(sizeof(FORM*) * ps->length);
      for(i=0;i < ps->length; i++)
	fl[i] = (FORM*)ckalloc(sizeof(FORM));

      /*******************************************************************/
      /* main formant tracking loop */
      /*******************************************************************/
      if(debug & DEB_ENTRY){
	printf("Entering main computation loop in dpform()\n");
      }
      for(i=0; i < ps->length; i++){	/* for all analysis frames... */

	ncan = 0;		/* initialize candidate mapping count to 0 */

	/* moderate the cost of frequency jumps by the relative amplitude */
	rmsdffact = pole[i]->rms;
	rmsdffact = rmsdffact/rmsmax;
	rmsdffact = rmsdffact * dffact;

	/* Get all likely mappings of the poles onto formants for this frame. */
	if(pole[i]->npoles){	/* if there ARE pole frequencies available... */
	  get_fcand(pole[i]->npoles,pole[i]->freq,pole[i]->band,nform,pcan);

	  /* Allocate space for this frame's candidates in the dp lattice. */
	  fl[i]->prept =  (short*)ckalloc(sizeof(short) * ncan);
	  fl[i]->cumerr = (double*)ckalloc(sizeof(double) * ncan);
	  fl[i]->cand =   (short**)ckalloc(sizeof(short*) * ncan);
	  for(j=0;j<ncan;j++){	/* allocate cand. slots and install candidates */
	    fl[i]->cand[j] = (short*)ckalloc(sizeof(short) * nform);
	    for(k=0; k<nform; k++)
	      fl[i]->cand[j][k] = pcan[j][k];
	  }
	}
	fl[i]->ncand = ncan;
	/* compute the distance between the current and previous mappings */
	for(j=0;j<ncan;j++){	/* for each CURRENT mapping... */
	  if( i ){		/* past the first frame? */
	    minerr = 0;
	    if(fl[i-1]->ncand) minerr = 2.0e30;
	    mincan = -1;
	    for(k=0; k < fl[i-1]->ncand; k++){ /* for each PREVIOUS map... */
	      for(pferr=0.0, l=0; l<nform; l++){
		ic = fl[i]->cand[j][l];
		ip = fl[i-1]->cand[k][l];
		if((ic >= 0)	&& (ip >= 0)){
		  ftemp = 2.0 * fabs(pole[i]->freq[ic] - pole[i-1]->freq[ip])/
		           (pole[i]->freq[ic] + pole[i-1]->freq[ip]);
    /*		  ftemp = pole[i]->freq[ic] - pole[i-1]->freq[ip];
		  if(ftemp >= 0.0)
		    ftemp = ftemp/pole[i-1]->freq[ip];
		  else
		    ftemp = ftemp/pole[i]->freq[ic]; */
		  /* cost prop. to SQUARE of deviation to discourage large jumps */
		  pferr += ftemp * ftemp;
		}
		else pferr += MISSING;
	      }
	      /* scale delta-frequency cost and add in prev. cum. cost */
	      conerr = (rmsdffact * pferr) + fl[i-1]->cumerr[k]; 
	      if(conerr < minerr){
		minerr = conerr;
		mincan = k;
	      }
	    }			/* end for each PREVIOUS mapping... */
	  }	else {		/* (i.e. if this is the first frame... ) */
	    minerr = 0;
	  }

	  fl[i]->prept[j] = mincan; /* point to best previous mapping */
	  /* (Note that mincan=-1 if there were no candidates in prev. fr.) */
	  /* Compute the local costs for this current mapping. */
	  for(k=0, berr=0, ferr=0, fbias=0; k<nform; k++){
	    ic = fl[i]->cand[j][k];
	    if(ic >= 0){
	      if( !k ){		/* F1 candidate? */
		ftemp = pole[i]->freq[ic];
		merger = (domerge &&
			  (ftemp == pole[i]->freq[fl[i]->cand[j][1]]))?
			  merge_cost: 0.0;
	      }
	      berr += pole[i]->band[ic];
	      ferr += (fabs(pole[i]->freq[ic]-fnom[k])/fnom[k]);
	      fbias += pole[i]->freq[ic];
	    } else {		/* if there was no freq. for this formant */
	      fbias += fnom[k];
	      berr += NOBAND;
	      ferr += MISSING;
	    }
	  }

	  /* Compute the total cost of this mapping and best previous. */
	  fl[i]->cumerr[j] = (FBIAS * fbias) + (bfact * berr) + merger +
	                     (ffact * ferr) + minerr;
	}			/* end for each CURRENT mapping... */

	if(debug & DEB_LPC_PARS){
	  printf("\nFrame %4d  # candidates:%3d stat:%f prms:%f",i,ncan,rmsdffact,pole[i]->rms);
	  for (j=0; j<ncan; j++){
	    printf("\n	");
	    for(k=0; k<nform; k++)
	      if(pcan[j][k] >= 0)
		printf("%6.0f ",pole[i]->freq[fl[i]->cand[j][k]]);
	      else
		printf("  NA   ");
	    printf("  cum:%7.2f pp:%d",fl[i]->cumerr[j], fl[i]->prept[j]);
	  }
	}
      }				/* end for all analysis frames... */	
      /**************************************************************************/

      /* Pick the candidate in the final frame with the lowest cost. */
      /* Starting with that min.-cost cand., work back thru the lattice. */
      if(debug & DEB_ENTRY){
	printf("Entering backtrack loop in dpform()\n");
      }
      dmaxc = 0;
      dminc = 100;
      dcountc = dcountf = 0;
      for(mincan = -1, i=ps->length - 1; i>=0; i--){
	if(debug & DEB_LPC_PARS){
	  printf("\nFrame:%4d mincan:%2d ncand:%2d ",i,mincan,fl[i]->ncand);
	}
	if(mincan < 0)		/* need to find best starting candidate? */
	  if(fl[i]->ncand){	/* have candidates at this frame? */
	    minerr = fl[i]->cumerr[0];
	    mincan = 0;
	    for(j=1; j<fl[i]->ncand; j++)
	      if( fl[i]->cumerr[j] < minerr ){
		minerr = fl[i]->cumerr[j];
		mincan = j;
	      }
	  }
	if(mincan >= 0){	/* if there is a "best" candidate at this frame */
	  if((j = fl[i]->ncand) > dmaxc) dmaxc = j;
	  else
	    if( j < dminc) dminc = j;
	  dcountc += j;
	  dcountf++;
	  for(j=0; j<nform; j++){
	    k = fl[i]->cand[mincan][j];
	    if(k >= 0){
	      fr[j][i] = pole[i]->freq[k];
	      if(debug & DEB_LPC_PARS){
		printf("%6.0f",fr[j][i]);
	      }
	      ba[j][i] = pole[i]->band[k];
	    } else {		/* IF FORMANT IS MISSING... */
	      if(i < ps->length - 1){
		fr[j][i] = fr[j][i+1]; /* replicate backwards */
		ba[j][i] = ba[j][i+1];
	      } else {
		fr[j][i] = fnom[j]; /* or insert neutral values */
		ba[j][i] = NOBAND;
	      }
	      if(debug & DEB_LPC_PARS){
		printf("%6.0f",fr[j][i]);
	      }
	    }
	  }
	  mincan = fl[i]->prept[mincan];
	} else {		/* if no candidates, fake with "nominal" frequencies. */
	  for(j=0; j < nform; j++){
	    fr[j][i] = fnom[j];
	    ba[j][i] = NOBAND;
	    if(debug & DEB_LPC_PARS){
	      printf("%6.0f",fr[j][i]);
	    }
	  } 
	}			/* note that mincan will remain =-1 if no candidates */
      }				/* end unpacking formant tracks from the dp lattice */
      /* Deallocate all the DP lattice work space. */
      /*if(debug & DEB_ENTRY){
	printf("%s complete; max. cands:%d  min. cands.:%d average cands.:%f\n",
	     fbs->name,dmaxc,dminc,((double)dcountc)/dcountf);
	printf("Entering memory deallocation in dpform()\n");
      }*/
      for(i=ps->length - 1; i>=0; i--){
	if(fl[i]->ncand){
	  if(fl[i]->cand) {
	    for(j=0; j<fl[i]->ncand; j++) ckfree((void *)fl[i]->cand[j]);
	    ckfree((void *)fl[i]->cand);
	    ckfree((void *)fl[i]->cumerr);
	    ckfree((void *)fl[i]->prept);
	  }
	}
      }
      for(i=0; i<ps->length; i++)	ckfree((void *)fl[i]);
      ckfree((void *)fl);
      fl = 0;
      
      for(i=0; i<ps->length; i++) {
	ckfree((void *)pole[i]->freq);
	ckfree((void *)pole[i]->band);
	ckfree((void *)pole[i]);
      }
      ckfree((void *)pole);

      /* Deallocate space for the raw candidate aray. */
      for(i=0;i<MAXCAN;i++) ckfree((void *)pcan[i]);
      ckfree((void *)pcan);

      fbs = Snack_NewSound(ps->samprate, SNACK_FLOAT, nform * 2);
      Snack_ResizeSoundStorage(fbs, ps->length);
      for (i = 0; i < ps->length; i++) {
	for (j = 0; j < nform * 2; j++) {
	  Snack_SetSample(fbs, j, i, (float)fr[j][i]);
	}
      }
      fbs->length = ps->length;

      for(i = 0; i < nform*2; i++) ckfree((void *)fr[i]);
      ckfree((void *)fr);

      return(fbs);
    } else
      printf("Can't create a new Signal in dpform()\n");
  } else
    printf("Bad data pointers passed into dpform()\n");
  return(NULL);
}

/* lpc_poles.c */

/* computation and I/O routines for dealing with LPC poles */

#define MAXORDER 30

extern int dlpcwtd(double *s, int *ls, double *p, int *np, double *c, double *phi, double *shi, double *xl, double *w);
extern int formant(int lpc_order, double s_freq, double *lpca, int *n_form, double *freq, double *band, int init);
extern int lpc(int lpc_ord, double lpc_stabl, int wsize, short *data, double *lpca,
        double *ar, double *lpck, double *normerr, double *rms, double preemp, int type);
extern int lpcbsa(int np, double lpc_stabl, int wind, short *data, double *lpc, double *rho, double *nul1, double *nul2, double *energy, double preemp);
extern int w_covar(short *xx, int *m, int n, int istrt, double *y, double *alpha, double *r0, double preemp, int w_type);

/*************************************************************************/
double integerize(register double time, register double freq)
{
  register int i;

  i = (int) (.5 + (freq * time));
  return(((double)i)/freq);
}

/*	Round the argument to the nearest integer.			*/
int eround(register double	flnum)
{
	return((flnum >= 0.0) ? (int)(flnum + 0.5) : (int)(flnum - 0.5));
}

/*************************************************************************/
Sound *lpc_poles(Sound *sp, double wdur, double frame_int, int lpc_ord, double preemp, int lpc_type, int w_type)
{
  int i, j, size, step, nform, init, nfrm;
  POLE **pole;
  double lpc_stabl = 70.0, energy, lpca[MAXORDER], normerr,
         *bap=NULL, *frp=NULL, *rhp=NULL;
  short *datap, *dporg;
  Sound *lp;

  if(lpc_type == 1) { /* force "standard" stabilized covariance (ala bsa) */
    wdur = 0.025;
    preemp = exp(-62.831853 * 90. / sp->samprate); /* exp(-1800*pi*T) */
  }
  if((lpc_ord > MAXORDER) || (lpc_ord < 2)/* || (! ((short**)sp->data)[0])*/)
    return(NULL);
  /*  np = (char*)new_ext(sp->name,"pole");*/
  wdur = integerize(wdur,(double)sp->samprate);
  frame_int = integerize(frame_int,(double)sp->samprate);
  nfrm= 1 + (int) (((((double)sp->length)/sp->samprate) - wdur)/(frame_int));
  if(nfrm >= 1/*lp->buff_size >= 1*/) {
    size = (int) (.5 + (wdur * sp->samprate));
    step = (int) (.5 + (frame_int * sp->samprate));
    pole = (POLE**)ckalloc(nfrm/*lp->buff_size*/ * sizeof(POLE*));
    datap = dporg = (short *) ckalloc(sizeof(short) * sp->length);
    for (i = 0; i < Snack_GetLength(sp); i++) {
      datap[i] = (short) Snack_GetSample(sp, 0, i);
    }
    for(j=0, init=TRUE/*, datap=((short**)sp->data)[0]*/; j < nfrm/*lp->buff_size*/;j++, datap += step){
      pole[j] = (POLE*)ckalloc(sizeof(POLE));
      pole[j]->freq = frp = (double*)ckalloc(sizeof(double)*lpc_ord);
      pole[j]->band = bap = (double*)ckalloc(sizeof(double)*lpc_ord);

      switch(lpc_type) {
      case 0:
	if(! lpc(lpc_ord,lpc_stabl,size,datap,lpca,rhp,NULL,&normerr,
		 &energy, preemp, w_type)){
	  printf("Problems with lpc in lpc_poles()");
	  break;
	}
	break;
      case 1:
	if(! lpcbsa(lpc_ord,lpc_stabl,size,datap,lpca,rhp,NULL,&normerr,
		    &energy, preemp)){
          printf("Problems with lpc in lpc_poles()");
	  break;
	}
	break;
      case 2:
	{
	  int Ord = lpc_ord;
	  double alpha, r0;

	  w_covar(datap, &Ord, size, 0, lpca, &alpha, &r0, preemp, 0);
	  if((Ord != lpc_ord) || (alpha <= 0.0))
	    printf("Problems with covar(); alpha:%f  Ord:%d\n",alpha,Ord);
	  energy = sqrt(r0/(size-Ord));
	}
	break;
      }
      pole[j]->change = 0.0;
       /* don't waste time on low energy frames */
       if((pole[j]->rms = energy) > 1.0){
	 formant(lpc_ord,(double)sp->samprate, lpca, &nform, frp, bap, init);
	 pole[j]->npoles = nform;
	 init=FALSE;		/* use old poles to start next search */
       } else {			/* write out no pole frequencies */
	 pole[j]->npoles = 0;
	 init = TRUE;		/* restart root search in a neutral zone */
       }
/*     if(debug & 4) {
	 printf("\nfr:%4d np:%4d rms:%7.0f  ",j,pole[j]->npoles,pole[j]->rms);
	 for(k=0; k<pole[j]->npoles; k++)
	   printf(" %7.1f",pole[j]->freq[k]);
	 printf("\n                   ");
	 for(k=0; k<pole[j]->npoles; k++)
	   printf(" %7.1f",pole[j]->band[k]);
	 printf("\n");
	 }*/
     } /* end LPC pole computation for all lp->buff_size frames */
    /*     lp->data = (caddr_t)pole;*/
    ckfree((void *)dporg);
    lp = Snack_NewSound((int)(1.0/frame_int), LIN16, lpc_ord);
    Snack_ResizeSoundStorage(lp, nfrm);
    for (i = 0; i < nfrm; i++) {
      for (j = 0; j < lpc_ord; j++) {
      Snack_SetSample(lp, j, i, (float)pole[i]->freq[j]);
      }
    }
    lp->length = nfrm;
    lp->extHead = (char *)pole;
    return(lp);
  } else {
    printf("Bad buffer size in lpc_poles()\n");
  }
  return(NULL);
}

/**********************************************************************/
double frand()
{
  return (((double)rand())/(double)RAND_MAX);
}
    
/**********************************************************************/
/* a quick and dirty interface to bsa's stabilized covariance LPC */
#define NPM	30	/* max lpc order		*/

int lpcbsa(int np, double lpc_stabl, int wind, short *data, double *lpc, double *rho, double *nul1, double *nul2, double *energy, double preemp)
{
  static int i, mm, owind=0, wind1;
  static double w[1000];
  double rc[NPM],phi[NPM*NPM],shi[NPM],sig[1000];
  double xl = .09, fham, amax;
  register double *psp1, *psp3, *pspl;

  if(owind != wind) {		/* need to compute a new window? */
    fham = 6.28318506 / wind;
    for(psp1=w,i=0;i<wind;i++,psp1++)
      *psp1 = .54 - .46 * cos(i * fham);
    owind = wind;
  }
  wind += np + 1;
  wind1 = wind-1;

  for(psp3=sig,pspl=sig+wind; psp3 < pspl; )
    *psp3++ = (double)(*data++) + .016 * frand() - .008;
  for(psp3=sig+1,pspl=sig+wind;psp3<pspl;psp3++)
    *(psp3-1) = *psp3 - preemp * *(psp3-1);
  for(amax = 0.,psp3=sig+np,pspl=sig+wind1;psp3<pspl;psp3++)
    amax += *psp3 * *psp3;
  *energy = sqrt(amax / (double)owind);
  amax = 1.0/(*energy);
	
  for(psp3=sig,pspl=sig+wind1;psp3<pspl;psp3++)
    *psp3 *= amax;
  if((mm=dlpcwtd(sig,&wind1,lpc,&np,rc,phi,shi,&xl,w))!=np) {
    printf("LPCWTD error mm<np %d %d\n",mm,np);
    return(FALSE);
  }
  return(TRUE);
}

/*      ----------------------------------------------------------      */

int ratprx(double a, int *k, int *l, int qlim)
{
    double aa, af, q, em, qq = 0, pp = 0, ps, e;
    int	ai, ip, result = FALSE;
    
    aa = fabs(a);
    ai = (int) aa;
    af = aa - ai;
    q = 0;
    em = 1.0;
    while(++q <= qlim) {
	ps = q * af;
	ip = (int) (ps + 0.5);
	e = fabs((ps - (double)ip)/q);
	if(e < em) {
	    em = e;
	    pp = ip;
	    qq = q;
	    result = TRUE;
	}
    };
    *k = (int) ((ai * qq) + pp);
    *k = (a > 0)? *k : -(*k);
    *l = (int) qq;
    return(result);
}

/* ----------------------------------------------------------------------- */

extern float *downsample(float *input, int samsin, int state_idx, double freq, int *samsout, int decimate, int first_time, int last_time);

Sound *Fdownsample(Sound *s, double freq2, int start, int end)
{
  float	*bufin, *bufout, *bufp;
  int	frame_size = 1024, act_size, first_time, last_time;
  double	ratio, freq1;
  int	ncoeff, insert, decimate, total_samps, out_samps, ndone;
  Sound *so;

  register int i, j;

  freq1 = s->samprate;
  ratio = freq2/freq1;

  if (!ratprx(ratio, &insert, &decimate, 10))
    return(NULL);

  if (decimate <= insert)
    return(NULL);

  freq1 *= (double)insert;
  freq2 = freq1/((double)decimate);

  /* filter length used in downsample(): 5ms */
  ncoeff = ((int)(freq1 * 0.005))/2 + 1;

  total_samps = (end - start + 1) * insert;
  if (total_samps < (ncoeff * decimate * 3))	/* signal too short */
    return(NULL);

  if ((bufin = (float *) ckalloc(sizeof(float) * total_samps))) {
    for (bufp = bufin, i = start; i <= end; i++) {
      *bufp++ = Snack_GetSample(s, 0, i) * ((float)insert);
      for(j = 1; j < insert; j++)
        *bufp++ = 0.0f;	/* insert zeros to boost the sampling frequency */
    }

    if ((frame_size * 2) > total_samps)
      frame_size = (total_samps + 1)/2;

    frame_size -= frame_size % decimate;

    first_time = 1;	/* new filter coefficients need to be computed */

    for (ndone = 0, last_time = 0; !last_time; ndone += act_size) {
      act_size = total_samps - ncoeff - ndone;
      if (act_size > frame_size) {
        act_size = frame_size;
        out_samps = act_size/decimate;
      } else {
        out_samps = act_size/decimate;
        if (!first_time && ((act_size + ncoeff) <= frame_size)) {
          act_size += ncoeff;
          last_time = 1;
        } else
          act_size = out_samps * decimate;
      }

      if ((bufout = downsample(bufin+ndone, total_samps-ndone, act_size, freq1,
                               &out_samps, decimate, first_time, last_time))) {
        if (first_time) {
          first_time = 0;
          so = Snack_NewSound((int)freq2, LIN16, s->nchannels);
          if (!so) {
            printf("Can't create a new Signal in downsample()\n");
            break;
          }
          Snack_ResizeSoundStorage(so, total_samps/decimate);
          so->length = 0;
        }

        Snack_PutSoundData(so, so->length, bufout, out_samps);
        so->length += out_samps;
      } else {
        printf("Problems in downsample()\n");
        break;
      }
    }
    ckfree((void *)bufin);
  } else
       printf("Can't create a new Signal in downsample()\n");
  
  return(so);
}

/*      ----------------------------------------------------------      */

extern void do_ffir(register float *buf, register int in_samps, register float *bufo, register int *out_samps, int idx, register int ncoef, float *fc, register int invert, register int skip, register int init);

Sound *highpass(Sound *s)
{

  float *datain, *dataout;
  static float *lcf;
  static int len = 0;
  double scale, fn;
  register int i;
  int	frame_size = 1024, act_size, total_samps, out_samps, ndone, init;
  Sound *so;

  /*  Header *h, *dup_header();*/
  
#define LCSIZ 101
  /* This assumes the sampling frequency is 10kHz and that the FIR
     is a Hanning function of (LCSIZ/10)ms duration. */

  if(!len) {		/* need to create a Hanning FIR? */
    lcf = (float *) ckalloc(sizeof(float) * LCSIZ);
    len = 1 + (LCSIZ/2);
    fn = M_PI * 2.0 / (LCSIZ - 1);
    scale = 1.0/(.5 * LCSIZ);
    for (i=0; i < len; i++)
      lcf[i] = (float) (scale * (.5 + (.4 * cos(fn * ((double)i)))));
  }

  total_samps = s->length;
  if (total_samps < (len * 3))
    total_samps = len * 3;

  datain = (float *) ckalloc(sizeof(float) * total_samps);
  dataout = (float *) ckalloc(sizeof(float) * total_samps);
  if (!datain || !dataout) {
    printf("Can't create a new Signal in highpass()\n");
    return(NULL);
  }

  Snack_GetSoundData(s, 0, datain, s->length);

  for (i = s->length; i < total_samps; i++)
    datain[i] = 0.0;

  if (frame_size > total_samps)
    frame_size = total_samps;

  for (ndone = 0, init = 1; !(init & 2); ndone += act_size) {
    act_size = total_samps - len - ndone;
    if (act_size > frame_size) {
      out_samps = act_size = frame_size;
    } else {
      out_samps = act_size;
      if (!(init & 1) && ((act_size + len) <= frame_size)) {
        act_size += len;
        init = 2;
      }
    }

    do_ffir(datain+ndone, total_samps-ndone, dataout+ndone,
            &out_samps, act_size, len, lcf, 1, 1, init);

    if (init & 1)
      init = 0;
  }

  so = Snack_NewSound(s->samprate, LIN16, s->nchannels);
  if (!so) {
    printf("Can't create a new Signal in highpass()\n");
  } else {
    Snack_ResizeSoundStorage(so, s->length);
    Snack_PutSoundData(so, 0, dataout, s->length);
    so->length = s->length;
  }

  ckfree((void *)dataout);
  ckfree((void *)datain);
  return(so);
}

int
formantCmd(Sound *s, Tcl_Interp *interp, int objc,
	   Tcl_Obj *CONST objv[])
{
  int nform, i,j, lpc_ord, lpc_type, w_type;
  char *w_type_str = NULL;
  double frame_int, wdur, 
  ds_freq, nom_f1 = -10.0, preemp;
  double cor_wdur;
  Sound *dssnd = NULL, *hpsnd = NULL, *polesnd = NULL;
  Sound *formantsnd = NULL, *hpsrcsnd, *polesrcsnd;
  Tcl_Obj *list;
  int arg, startpos = 0, endpos = -1;
  static CONST84 char *subOptionStrings[] = {
    "-start", "-end", "-progress",
    "-framelength", "-preemphasisfactor", "-numformants",
    "-lpcorder", "-windowlength", "-windowtype", "-lpctype",
    "-ds_freq", "-nom_f1_freq", NULL
  };
  enum subOptions {
    START, END, PROGRESS, FRAME, PREEMP, NUMFORM, ORDER, WINLEN,
    WINTYPE, LPCTYPE, DSFREQ, NOMFREQ
  };

  lpc_ord = 12;
  lpc_type = 0;			/* use bsa's stabilized covariance if != 0 */
  w_type = 2;			/* window type: 0=rectangular; 1=Hamming; 2=cos**4 */
  ds_freq = 10000.0;
  wdur = .049;			/* for LPC analysis */
  cor_wdur = .01;		/* for crosscorrelation F0 estimator */
  frame_int = .01;
  preemp = .7;
  nform = 4;

  for (arg = 2; arg < objc; arg += 2) {
    int index;
	
    if (Tcl_GetIndexFromObj(interp, objv[arg], subOptionStrings,
			    "option", 0, &index) != TCL_OK) {
      return TCL_ERROR;
    }

    if (arg + 1 == objc) {
      Tcl_AppendResult(interp, "No argument given for ",
		       subOptionStrings[index], " option", (char *) NULL);
      return TCL_ERROR;
    }
    
    switch ((enum subOptions) index) {
    case START:
      {
	if (Tcl_GetIntFromObj(interp, objv[arg+1], &startpos) != TCL_OK)
	  return TCL_ERROR;
	break;
      }
    case END:
      {
	if (Tcl_GetIntFromObj(interp, objv[arg+1], &endpos) != TCL_OK)
	  return TCL_ERROR;
	break;
      }
    case PROGRESS:
      {
	char *str = Tcl_GetStringFromObj(objv[arg+1], NULL);
	
	if (strlen(str) > 0) {
	  Tcl_IncrRefCount(objv[arg+1]);
	  s->cmdPtr = objv[arg+1];
	}
	break;
      }
    case FRAME:
      {
	if (Tcl_GetDoubleFromObj(interp, objv[arg+1], &frame_int)
	    != TCL_OK)
	  return TCL_ERROR;
	break;
      }
    case PREEMP:
      {
	if (Tcl_GetDoubleFromObj(interp, objv[arg+1], &preemp)
	    != TCL_OK)
	  return TCL_ERROR;
	break;
      }
    case NUMFORM:
      {
	if (Tcl_GetIntFromObj(interp, objv[arg+1], &nform) != TCL_OK)
	  return TCL_ERROR;
	break;
      }
    case ORDER:
      {
	if (Tcl_GetIntFromObj(interp, objv[arg+1], &lpc_ord) != TCL_OK)
	  return TCL_ERROR;
	break;
      }
    case WINLEN:
      {
	if (Tcl_GetDoubleFromObj(interp, objv[arg+1], &wdur)
	    != TCL_OK)
	  return TCL_ERROR;
	break;
      }
    case WINTYPE:
      {
	w_type_str = Tcl_GetStringFromObj(objv[arg+1], NULL);
	break;
      }
    case LPCTYPE:
      {
	if (Tcl_GetIntFromObj(interp, objv[arg+1], &lpc_type) != TCL_OK)
	  return TCL_ERROR;
	break;
      }
    case DSFREQ:
      {
	if (Tcl_GetDoubleFromObj(interp, objv[arg+1], &ds_freq)
	    != TCL_OK)
	  return TCL_ERROR;
	break;
      }
    case NOMFREQ:
      {
	if (Tcl_GetDoubleFromObj(interp, objv[arg+1], &nom_f1)
	    != TCL_OK)
	  return TCL_ERROR;
	break;
      }
    }
  }
  if (startpos < 0) startpos = 0;
  if (endpos >= (s->length - 1) || endpos == -1)
    endpos = s->length - 1;
  if (startpos > endpos) return TCL_OK;
  
  /*
   * Check for errors in specifying parameters
   */

  if(nform > (lpc_ord-4)/2){
    Tcl_AppendResult(interp, "Number of formants must be <= (lpc order - 4)/2", NULL);
    return TCL_ERROR;
  }
  
  if(nform > MAXFORMANTS){
    Tcl_AppendResult(interp, "A maximum of 7 formants are supported at this time", NULL);
    return TCL_ERROR;
  }

  if (s->storeType != SOUND_IN_MEMORY ) {
    Tcl_AppendResult(interp, "formant only works with in-memory sounds",
		     (char *) NULL);
    return TCL_ERROR;
  }

  if (w_type_str != NULL) {
    int len = strlen(w_type_str);
    if (strncasecmp(w_type_str, "rectangular", len) == 0 ||
	strncasecmp(w_type_str, "0", len) == 0) {
      w_type = 0;
    } else if (strncasecmp(w_type_str, "hamming", len) == 0 ||
	       strncasecmp(w_type_str, "1", len) == 0) {
      w_type = 1;
    } else if (strncasecmp(w_type_str, "cos^4", len) == 0 ||
	       strncasecmp(w_type_str, "2", len) == 0) {
      w_type = 2;
    } else if (strncasecmp(w_type_str, "hanning", len) == 0 ||
	       strncasecmp(w_type_str, "3", len) == 0) {
      w_type = 3;
    } else {
      Tcl_AppendResult(interp, "unknown window type: ", w_type_str, 
		       (char *) NULL);
      return TCL_ERROR;
    }
  }

  Snack_ProgressCallback(s->cmdPtr, interp,"Computing formants",0.05);

  if(ds_freq < s->samprate) {
    dssnd = Fdownsample(s,ds_freq, startpos, endpos);
  }

  Snack_ProgressCallback(s->cmdPtr, interp, "Computing formants",
			 0.5);

  hpsrcsnd = (dssnd ? dssnd : s);
    
  if (preemp < 1.0) { /* be sure DC and rumble are gone! */
    hpsnd = highpass(hpsrcsnd);
  }

  Snack_ProgressCallback(s->cmdPtr, interp, "Computing formants",
			 0.6);

  polesrcsnd = (hpsnd ? hpsnd : s);
  
  if(!(polesnd = lpc_poles(polesrcsnd, wdur, frame_int, lpc_ord,
			   preemp, lpc_type, w_type))) {
    Tcl_AppendResult(interp, "Problems in lpc_poles()", NULL);
    return TCL_ERROR;
  }

  Snack_ProgressCallback(s->cmdPtr, interp, "Computing formants",
			 0.7);

  /* LPC poles are now available for the formant estimator. */
  if (!(formantsnd = dpform(polesnd, nform, nom_f1))) {
    Tcl_AppendResult(interp, "Problems in dpform()", NULL);
    return TCL_ERROR;
  }

  Snack_ProgressCallback(s->cmdPtr, interp, "Computing formants",
			 0.95);

  /*
    SaveSound(formantsnd, interp, "outt.wav", NULL,
    0, NULL, 0, formantsnd->length, WAV_STRING);
  */

  if (dssnd) Snack_DeleteSound(dssnd);
  if (hpsnd) Snack_DeleteSound(hpsnd);
  Snack_DeleteSound(polesnd);

  list = Tcl_NewListObj(0, NULL);
  
  for (i = 0; i < formantsnd->length; i++) {
    Tcl_Obj *frameList;
    frameList = Tcl_NewListObj(0, NULL);
    Tcl_ListObjAppendElement(interp, list, frameList);
    for (j = 0; j < nform * 2; j++) {
      Tcl_ListObjAppendElement(interp, frameList,
	  Tcl_NewDoubleObj((double) Snack_GetSample(formantsnd, j, i)));
    }
  }

  Snack_DeleteSound(formantsnd);

  Snack_ProgressCallback(s->cmdPtr, interp,"Computing formants",1.0);

  Tcl_SetObjResult(interp, list);

  return TCL_OK;
}
