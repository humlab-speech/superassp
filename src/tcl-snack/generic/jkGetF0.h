/*
 * This software has been licensed to the Centre of Speech Technology, KTH
 * by Microsoft Corp. with the terms in the accompanying file BSD.txt,
 * which is a BSD style license.
 *
 *    "Copyright (c) 1990-1996 Entropic Research Laboratory, Inc. 
 *                   All rights reserved"
 *
 * Written by:  David Talkin
 * Checked by:
 * Revised by:
 * @(#)f0.h	1.4 9/9/96 ERL
 * Brief description:
 *
 */



/* f0.h */
/* Some definitions used by the "Pitch Tracker Software". */
       
typedef struct f0_params {
float cand_thresh,	/* only correlation peaks above this are considered */
      lag_weight,	/* degree to which shorter lags are weighted */
      freq_weight,	/* weighting given to F0 trajectory smoothness */
      trans_cost,	/* fixed cost for a voicing-state transition */
      trans_amp,	/* amplitude-change-modulated VUV trans. cost */
      trans_spec,	/* spectral-change-modulated VUV trans. cost */
      voice_bias,	/* fixed bias towards the voiced hypothesis */
      double_cost,	/* cost for octave F0 jumps */
      mean_f0,		/* talker-specific mean F0 (Hz) */
      mean_f0_weight,	/* weight to be given to deviations from mean F0 */
      min_f0,		/* min. F0 to search for (Hz) */
      max_f0,		/* max. F0 to search for (Hz) */
      frame_step,	/* inter-frame-interval (sec) */
      wind_dur;		/* duration of correlation window (sec) */
int   n_cands,		/* max. # of F0 cands. to consider at each frame */
      conditioning;     /* Specify optional signal pre-conditioning. */
} F0_params;

/* Possible values returned by the function f0(). */
#define F0_OK		0
#define F0_NO_RETURNS	1
#define F0_TOO_FEW_SAMPLES	2
#define F0_NO_INPUT	3
#define F0_NO_PAR	4
#define F0_BAD_PAR	5
#define F0_BAD_INPUT	6
#define F0_INTERNAL_ERR	7

/* Bits to specify optional pre-conditioning of speech signals by f0() */
/* These may be OR'ed together to specify all preprocessing. */
#define F0_PC_NONE	0x00		/* no pre-processing */
#define F0_PC_DC	0x01		/* remove DC */
#define F0_PC_LP2000	0x02		/* 2000 Hz lowpass */
#define F0_PC_HP100	0x04		/* 100 Hz highpass */
#define F0_PC_AR	0x08		/* inf_order-order LPC inverse filter */
#define F0_PC_DIFF	0x010		/* 1st-order difference */

#define Fprintf (void)fprintf

/* f0_structs.h */

#define BIGSORD 100

typedef struct cross_rec { /* for storing the crosscorrelation information */
	float	rms;	/* rms energy in the reference window */
	float	maxval;	/* max in the crosscorr. fun. q15 */
	short	maxloc; /* lag # at which max occured	*/
	short	firstlag; /* the first non-zero lag computed */
	float	*correl; /* the normalized corsscor. fun. q15 */
} Cross;

typedef struct dp_rec { /* for storing the DP information */
	short	ncands;	/* # of candidate pitch intervals in the frame */
	short	*locs; /* locations of the candidates */
	float	*pvals; /* peak values of the candidates */
	float	*mpvals; /* modified peak values of the candidates */
	short	*prept; /* pointers to best previous cands. */
	float	*dpvals; /* cumulative error for each candidate */
} Dprec;

typedef struct windstat_rec {  /* for lpc stat measure in a window */
    float rho[BIGSORD+1];
    float err;
    float rms;
} Windstat;

typedef struct sta_rec {  /* for stationarity measure */
  float *stat;
  float *rms;
  float *rms_ratio;
} Stat;

typedef struct frame_rec{
  Cross *cp;
  Dprec *dp;
  float rms;
  struct frame_rec *next;
  struct frame_rec *prev;
} Frame;

extern F0_params *new_f0_params(void);
extern int eround(double x);
extern int window(float *din, float *dout, int n, float preemp, int type);
extern int get_window(float *dout, int n, int type);
extern void get_fast_cands(float *fdata, float *fdsdata, int ind, int step,
                           int size, int dec, int start, int nlags,
                           float *engref, int *maxloc, float *maxval,
                           Cross *cp, float *peaks, int *locs,
                           int *ncand, F0_params *par);
extern void crossf(float *data, int size, int start, int nlags,
                   float *engref, int *maxloc, float *maxval, float *correl);
extern void crossfi(float *data, int size, int start0, int nlags0, int nlags,
                    float *engref, int *maxloc, float *maxval, float *correl,
                    int *locs, int nlocs);
extern int  xlpc(int lpc_ord, float lpc_stabl, int wsize, float *data,
                 float *lpca, float *ar, float *lpck, float *normerr,
                 float *rms, float preemp, int type);
extern float wind_energy(float *data, int size, int w_type);
extern void xa_to_aca(float *a, float *b, float *c, int p);
extern float xitakura(int p, float *b, float *c, float *r, float *gain);
extern Frame *alloc_frame(int nlags, int ncands);
