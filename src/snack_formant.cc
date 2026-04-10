/*
 * snack_formant.cc — Snack formant tracker adapted for superassp
 *
 * Algorithms from the Snack Sound Toolkit (KTH, AT&T, Entropic):
 *   jkFormant.c  — LPC polynomial roots + DP formant tracking (David Talkin)
 *   sigproc2.c   — LPC analysis, root finding, windowing (AT&T / Entropic)
 *
 * Adapted for C++ compilation, Sound* replaced with direct array I/O,
 * wrapped in namespace snack_formant to avoid symbol conflicts.
 *
 * BSD license — see src/tcl-snack/generic/BSD.txt
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cfloat>

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#define M_PI_2 1.57079632679489661923
#endif

namespace snack_formant {

/* ================================================================ */
/*  Memory helpers (match Snack ckalloc/ckfree interface)           */
/* ================================================================ */
static void *ckalloc(size_t sz) { return std::malloc(sz); }
static void  ckfree(void *p)   { std::free(p); }
static void *ckrealloc(void *p, size_t sz) { return std::realloc(p, sz); }

/* ================================================================ */
/*  Struct definitions (from jkFormant.h — no Tcl dependency)       */
/* ================================================================ */

#define DEB_PAUSE    8
#define DEB_LPC_PARS 4
#define DEB_PARAMS   2
#define DEB_ENTRY    1
#define MAXFORMANTS 7

typedef struct form_latt {
  short ncand;
  short **cand;
  short *prept;
  double *cumerr;
} FORM;

typedef struct pole_array {
  double rms;
  double rms2;
  double f0;
  double pv;
  double change;
  short npoles;
  double *freq;
  double *band;
} POLE;

/* ================================================================ */
/*  Forward declarations                                            */
/* ================================================================ */
static int formant_find(int lpc_order, double s_freq, double *lpca,
                        int *n_form, double *freq, double *band, int init);
static int lbpoly(double *a, int order, double *rootr, double *rooti);

/* ================================================================ */
/*  sigproc2.c functions — LPC analysis + windowing                */
/* ================================================================ */

#define MAXORDER 60

/* --- Various windowing functions --- */
static void rwindow(short *din, double *dout, int n, double preemp) {
  if (preemp != 0.0) {
    short *p = din + 1;
    for (int i = 0; i < n; i++)
      *dout++ = (double)(*p++) - preemp * (double)(*din++);
  } else {
    for (int i = 0; i < n; i++)
      *dout++ = (double)(*din++);
  }
}

static void cwindow(short *din, double *dout, int n, double preemp) {
  static int wsize = 0;
  static double *wind = NULL;
  if (wsize != n) {
    if (wind) wind = (double *)ckrealloc(wind, n * sizeof(double));
    else wind = (double *)ckalloc(n * sizeof(double));
    wsize = n;
    double arg = M_PI * 2.0 / wsize;
    for (int i = 0; i < n; i++) {
      double co = 0.5 * (1.0 - std::cos((0.5 + (double)i) * arg));
      wind[i] = co * co * co * co;
    }
  }
  if (preemp != 0.0) {
    short *p = din + 1;
    double *q = wind;
    for (int i = 0; i < n; i++)
      *dout++ = *q++ * ((double)(*p++) - preemp * (double)(*din++));
  } else {
    double *q = wind;
    for (int i = 0; i < n; i++)
      *dout++ = *q++ * (double)(*din++);
  }
}

static void hwindow(short *din, double *dout, int n, double preemp) {
  static int wsize = 0;
  static double *wind = NULL;
  if (wsize != n) {
    if (wind) wind = (double *)ckrealloc(wind, n * sizeof(double));
    else wind = (double *)ckalloc(n * sizeof(double));
    wsize = n;
    double arg = M_PI * 2.0 / wsize;
    for (int i = 0; i < n; i++)
      wind[i] = 0.54 - 0.46 * std::cos((0.5 + (double)i) * arg);
  }
  if (preemp != 0.0) {
    short *p = din + 1;
    double *q = wind;
    for (int i = 0; i < n; i++)
      *dout++ = *q++ * ((double)(*p++) - preemp * (double)(*din++));
  } else {
    double *q = wind;
    for (int i = 0; i < n; i++)
      *dout++ = *q++ * (double)(*din++);
  }
}

static void hnwindow(short *din, double *dout, int n, double preemp) {
  static int wsize = 0;
  static double *wind = NULL;
  if (wsize != n) {
    if (wind) wind = (double *)ckrealloc(wind, n * sizeof(double));
    else wind = (double *)ckalloc(n * sizeof(double));
    wsize = n;
    double arg = M_PI * 2.0 / wsize;
    for (int i = 0; i < n; i++)
      wind[i] = 0.5 - 0.5 * std::cos((0.5 + (double)i) * arg);
  }
  if (preemp != 0.0) {
    short *p = din + 1;
    double *q = wind;
    for (int i = 0; i < n; i++)
      *dout++ = *q++ * ((double)(*p++) - preemp * (double)(*din++));
  } else {
    double *q = wind;
    for (int i = 0; i < n; i++)
      *dout++ = *q++ * (double)(*din++);
  }
}

static void w_window(short *din, double *dout, int n, double preemp, int type) {
  switch (type) {
  case 0: rwindow(din, dout, n, preemp); return;
  case 1: hwindow(din, dout, n, preemp); return;
  case 2: cwindow(din, dout, n, preemp); return;
  case 3: hnwindow(din, dout, n, preemp); return;
  default: rwindow(din, dout, n, preemp); return;
  }
}

/* --- Autocorrelation --- */
static void autoc(int windowsize, double *s, int p, double *r, double *e) {
  double sum0 = 0.0;
  for (int i = 0; i < windowsize; i++)
    sum0 += s[i] * s[i];
  r[0] = 1.0;
  if (sum0 == 0.0) {
    *e = 1.0;
    for (int i = 1; i <= p; i++) r[i] = 0.0;
    return;
  }
  for (int i = 1; i <= p; i++) {
    double sum = 0.0;
    for (int j = 0; j < windowsize - i; j++)
      sum += s[j] * s[j + i];
    r[i] = sum / sum0;
  }
  *e = std::sqrt(sum0 / windowsize);
}

/* --- Durbin recursion --- */
static void durbin(double *r, double *k, double *a, int p, double *ex) {
  double b[MAXORDER];
  double e = *r;
  k[0] = -r[1] / e;
  a[0] = k[0];
  e *= (1.0 - k[0] * k[0]);
  for (int i = 1; i < p; i++) {
    double s = 0.0;
    for (int j = 0; j < i; j++)
      s -= a[j] * r[i - j];
    k[i] = (s - r[i + 1]) / e;
    a[i] = k[i];
    for (int j = 0; j <= i; j++) b[j] = a[j];
    for (int j = 0; j < i; j++)
      a[j] += k[i] * b[i - j - 1];
    e *= (1.0 - k[i] * k[i]);
  }
  *ex = e;
}

/* --- LPC analysis --- */
static int lpc(int lpc_ord, double lpc_stabl, int wsize, short *data,
               double *lpca, double * /*ar_unused*/, double * /*lpck_unused*/,
               double * /*normerr_unused*/, double *rms, double preemp,
               int type) {
  static double *dwind = NULL;
  static int nwind = 0;
  double rho[MAXORDER + 1], k[MAXORDER], a[MAXORDER + 1];
  double en, er;

  if (wsize <= 0 || !data || lpc_ord > MAXORDER) return FALSE;

  if (nwind != wsize) {
    if (dwind) dwind = (double *)ckrealloc(dwind, wsize * sizeof(double));
    else dwind = (double *)ckalloc(wsize * sizeof(double));
    nwind = wsize;
  }

  w_window(data, dwind, wsize, preemp, type);
  autoc(wsize, dwind, lpc_ord, rho, &en);

  if (lpc_stabl > 1.0) {
    double ffact = 1.0 / (1.0 + std::exp((-lpc_stabl / 20.0) * std::log(10.0)));
    double rho_copy[MAXORDER + 1];
    rho_copy[0] = rho[0];
    for (int i = 1; i <= lpc_ord; i++) rho_copy[i] = ffact * rho[i];
    durbin(rho_copy, k, &a[1], lpc_ord, &er);
  } else {
    durbin(rho, k, &a[1], lpc_ord, &er);
  }

  a[0] = 1.0;
  if (lpca) std::memcpy(lpca, a, (lpc_ord + 1) * sizeof(double));
  if (rms) *rms = en;
  return TRUE;
}

/* --- Cholesky, covariance LPC helpers --- */
static void dlwrtrn(double *a, int *n, double *x, double *y) {
  x[0] = y[0] / a[0];
  double *pa = a + *n;
  for (int i = 1; i < *n; i++) {
    double sm = y[i];
    double *pa1 = pa;
    for (int j = 0; j < i; j++)
      sm -= pa1[j] * x[j];
    pa += *n;
    x[i] = sm / pa1[i];
  }
}

static void dreflpc(double *c, double *a, int *n) {
  a[0] = 1.0;
  a[1] = c[0];
  for (int i = 2; i <= *n; i++) {
    a[i] = c[i - 1];
    for (int j = 1; j <= i / 2; j++) {
      double ta1 = a[j] + c[i - 1] * a[i - j];
      a[i - j] += a[j] * c[i - 1];
      a[j] = ta1;
    }
  }
}

static int dchlsky(double *a, int *n, double *t, double *det) {
  *det = 1.0;
  int m = 0;
  int nn = *n;
  for (int i = 0; i < nn; i++) {
    for (int j = 0; j <= i; j++) {
      double sm = a[i * nn + j];
      for (int k = 0; k < j; k++)
        sm -= a[i * nn + k] * a[j * nn + k];
      if (i == j) {
        if (sm <= 0.0) return m;
        t[m] = std::sqrt(sm);
        *det *= t[m];
        a[i * nn + j] = t[m];
        m++;
        t[m - 1] = 1.0 / t[m - 1];
      } else {
        a[i * nn + j] = sm * t[j];
      }
    }
  }
  return m;
}

static int dcovlpc(double *p, double *s, double *a, int *n, double *c) {
  double d;
  int m = dchlsky(p, n, c, &d);
  dlwrtrn(p, n, c, s);
  double thres = 1.0e-31;
  int n1 = *n + 1;
  double ps = a[*n];
  double ps1 = 1.0e-8 * ps;
  int mm = 0;
  for (int i = 0; i < m; i++) {
    if (p[i * n1] < thres) break;
    mm++;
  }
  double ee = ps;
  for (int i = 0; i < mm; i++) {
    ee -= c[i] * c[i];
    if (ee < thres) { mm = i + 1; break; }
    if (ee < ps1) { /* losing accuracy */ }
    a[i] = std::sqrt(ee);
  }
  m = mm;
  c[0] = -c[0] / std::sqrt(ps);
  for (int i = 1; i < m; i++)
    c[i] = -c[i] / a[i - 1];
  dreflpc(c, a, &m);
  for (int i = m + 1; i <= *n; i++) a[i] = 0.0;
  return m;
}

static void dcwmtrx(double *s, int *ni, int *nl, int *np,
                     double *phi, double *shi, double *ps, double *w) {
  *ps = 0.0;
  for (int i = *ni; i < *nl; i++)
    *ps += s[i] * s[i] * w[i - *ni];

  for (int k = 0; k < *np; k++) {
    shi[k] = 0.0;
    int off = *ni - k - 1;
    for (int i = *ni; i < *nl; i++)
      shi[k] += s[i] * s[off + (i - *ni)] * w[i - *ni];
  }

  for (int i = 0; i < *np; i++)
    for (int j = 0; j <= i; j++) {
      double sm = 0.0;
      for (int l = *ni; l < *nl - i; l++)
        sm += s[l - i] * s[l - j] * w[l - *ni];
      phi[(*np) * i + j] = sm;
      phi[(*np) * j + i] = sm;
    }
}

static int dlpcwtd(double *s, int *ls, double *p, int *np, double *c,
                   double *phi, double *shi, double *xl, double *w) {
  int np1 = *np + 1;
  dcwmtrx(s, np, ls, np, phi, shi, &p[*np], w);

  if (*xl >= 1.0e-4) {
    /* Save diagonal + compute Cholesky */
    for (int i = 0; i < *np; i++)
      p[i] = phi[i * np1];
    p[*np] = p[*np]; /* pss */

    int mm = dchlsky(phi, np, c, &p[*np]); /* reuse p[*np] for det temp */
    (void)mm;
    dlwrtrn(phi, np, c, shi);

    double ee = p[*np], pss7 = 0.0000001 * p[*np];
    for (int i = 0; i < *np; i++) {
      if (phi[i * np1] < 0.0) break;
      ee -= c[i] * c[i];
      if (ee < 0.0) break;
      if (ee < pss7) { /* losing accuracy */ }
    }

    double pre = ee * *xl;
    /* Rebuild matrix with ridge regression */
    for (int i = 0; i < *np; i++) {
      for (int j = i + 1; j < *np; j++) {
        phi[j * (*np) + i] = phi[i * (*np) + j]; /* copy upper to lower */
      }
    }
    double pre3 = 0.375 * pre, pre2 = 0.25 * pre, pre0 = 0.0625 * pre;
    for (int i = 0; i < *np; i++) {
      phi[i * np1] = p[i] + pre3;
      if (i >= 1) {
        phi[(i - 1) * (*np) + i] = phi[i * (*np) + (i - 1)] =
            phi[(i - 1) * (*np) + i] - pre2;
      }
      if (i >= 2) {
        phi[(i - 2) * (*np) + i] = phi[i * (*np) + (i - 2)] =
            phi[(i - 2) * (*np) + i] + pre0;
      }
    }
    shi[0] -= pre2;
    if (*np > 1) shi[1] += pre0;
    p[*np] = p[*np] + pre3;
  }
  int m = dcovlpc(phi, shi, p, np, c);
  return m;
}

/* --- lpcbsa: stabilized covariance LPC --- */
static double frand_val() {
  return ((double)rand()) / (double)RAND_MAX;
}

#define NPM 30

static int lpcbsa(int np, double lpc_stabl, int wind, short *data,
                  double *lpc_out, double * /*rho*/, double * /*nul1*/,
                  double * /*nul2*/, double *energy, double preemp) {
  static int owind = 0;
  static double w[1000];
  double rc[NPM], phi[NPM * NPM], shi[NPM], sig[1000];
  double xl = 0.09;

  if (owind != wind) {
    double fham = 6.28318506 / wind;
    for (int i = 0; i < wind; i++)
      w[i] = 0.54 - 0.46 * std::cos(i * fham);
    owind = wind;
  }
  int wind_orig = wind;
  wind += np + 1;
  int wind1 = wind - 1;

  for (int i = 0; i < wind; i++)
    sig[i] = (double)(*data++) + 0.016 * frand_val() - 0.008;
  for (int i = 1; i < wind; i++)
    sig[i - 1] = sig[i] - preemp * sig[i - 1];
  double amax = 0.0;
  for (int i = np; i < wind1; i++)
    amax += sig[i] * sig[i];
  *energy = std::sqrt(amax / (double)wind_orig);
  amax = 1.0 / (*energy);

  for (int i = 0; i < wind1; i++)
    sig[i] *= amax;

  int mm = dlpcwtd(sig, &wind1, lpc_out, &np, rc, phi, shi, &xl, w);
  if (mm != np) return FALSE;
  return TRUE;
}

/* --- w_covar: covariance LPC analysis --- */
static int w_covar(short *xx, int *m, int n, int istrt,
                   double *y, double *alpha, double *r0,
                   double preemp, int w_type) {
  static double *x = NULL;
  static int nold = 0;
  static int mold = 0;
  static double *b = NULL, *beta = NULL, *grc = NULL, *cc = NULL;

  if ((n + 1) > nold) {
    if (x) ckfree(x);
    x = (double *)ckalloc((n + 1) * sizeof(double));
    std::memset(x, 0, (n + 1) * sizeof(double));
    nold = n + 1;
  }

  if (*m > mold) {
    if (b) ckfree(b); if (beta) ckfree(beta);
    if (grc) ckfree(grc); if (cc) ckfree(cc);
    int mnew = *m;
    b    = (double *)ckalloc(sizeof(double) * ((mnew + 1) * (mnew + 1) / 2 + 1));
    beta = (double *)ckalloc(sizeof(double) * (mnew + 3));
    grc  = (double *)ckalloc(sizeof(double) * (mnew + 3));
    cc   = (double *)ckalloc(sizeof(double) * (mnew + 3));
    mold = mnew;
  }

  w_window(xx, x, n, preemp, w_type);

  int ibeg = istrt - 1, ibeg1 = ibeg + 1, mp = *m + 1;
  int ibegm1 = ibeg - 1, ibeg2 = ibeg + 2, ibegmp = ibeg + mp;
  int msq = (*m + (*m) * (*m)) / 2;
  for (int i = 1; i <= msq; i++) b[i] = 0.0;
  *alpha = 0.0; cc[1] = 0.0; cc[2] = 0.0;
  for (int np = mp; np <= n; np++) {
    int np1 = np + ibegm1, np0 = np + ibeg;
    *alpha += x[np0] * x[np0];
    cc[1] += x[np0] * x[np1];
    cc[2] += x[np1] * x[np1];
  }
  *r0 = *alpha;
  b[1] = 1.0; beta[1] = cc[2];
  grc[1] = -cc[1] / cc[2];
  y[0] = 1.0; y[1] = grc[1];
  *alpha += grc[1] * cc[1];
  if (*m <= 1) return FALSE;

  int mf = *m;
  for (int minc = 2; minc <= mf; minc++) {
    for (int j = 1; j <= minc; j++) {
      int jp = minc + 2 - j;
      int n1 = ibeg1 + mp - jp, n2 = ibeg1 + n - minc, n3 = ibeg2 + n - jp;
      cc[jp] = cc[jp - 1] + x[ibegmp - minc] * x[n1] - x[n2] * x[n3];
    }
    cc[1] = 0.0;
    for (int np = mp; np <= n; np++)
      cc[1] += x[np + ibeg - minc] * x[np + ibeg];

    int msub = (minc * minc - minc) / 2, mm1 = minc - 1;
    b[msub + minc] = 1.0;
    for (int ip = 1; ip <= mm1; ip++) {
      int isub = (ip * ip - ip) / 2;
      if (beta[ip] <= 0.0) { *m = minc - 1; return TRUE; }
      double gam = 0.0;
      for (int j = 1; j <= ip; j++)
        gam += cc[j + 1] * b[isub + j];
      gam /= beta[ip];
      for (int jp = 1; jp <= ip; jp++)
        b[msub + jp] -= gam * b[isub + jp];
    }
    beta[minc] = 0.0;
    for (int j = 1; j <= minc; j++)
      beta[minc] += cc[j + 1] * b[msub + j];
    if (beta[minc] <= 0.0) { *m = minc - 1; return TRUE; }
    double s = 0.0;
    for (int ip = 1; ip <= minc; ip++)
      s += cc[ip] * y[ip - 1];
    grc[minc] = -s / beta[minc];
    for (int ip = 1; ip < minc; ip++)
      y[ip] += grc[minc] * b[msub + ip];
    y[minc] = grc[minc];
    s = grc[minc] * grc[minc] * beta[minc];
    *alpha -= s;
    if (*alpha <= 0.0) { if (minc < *m) *m = minc; return TRUE; }
  }
  return TRUE;
}

/* ================================================================ */
/*  Polynomial root finder (lbpoly from sigproc2.c)                */
/* ================================================================ */

#define MAX_ITS   100
#define MAX_TRYS  100
#define MAX_ERR   1.0e-6

static int qquad(double a, double b, double c,
                 double *r1r, double *r1i, double *r2r, double *r2i) {
  if (a == 0.0) {
    if (b == 0.0) return FALSE;
    *r1r = -c / b; *r1i = *r2r = *r2i = 0.0;
    return TRUE;
  }
  double numi = b * b - 4.0 * a * c;
  if (numi >= 0.0) {
    *r1i = *r2i = 0.0;
    if (b < 0.0) {
      double y = -b + std::sqrt(numi);
      *r1r = y / (2.0 * a);
      *r2r = (2.0 * c) / y;
    } else {
      double y = -b - std::sqrt(numi);
      *r1r = (2.0 * c) / y;
      *r2r = y / (2.0 * a);
    }
    return TRUE;
  } else {
    double den = 2.0 * a;
    *r1i = std::sqrt(-numi) / den;
    *r2i = -*r1i;
    *r2r = *r1r = -b / den;
    return TRUE;
  }
}

static int lbpoly(double *a, int order, double *rootr, double *rooti) {
  double b[MAXORDER], c[MAXORDER];
  double lim0 = 0.5 * std::sqrt(DBL_MAX);

  for (int ord = order; ord > 2; ord -= 2) {
    int ordm1 = ord - 1, ordm2 = ord - 2;
    if (std::fabs(rootr[ordm1]) < 1.0e-10) rootr[ordm1] = 0.0;
    if (std::fabs(rooti[ordm1]) < 1.0e-10) rooti[ordm1] = 0.0;
    double p = -2.0 * rootr[ordm1];
    double q = rootr[ordm1] * rootr[ordm1] + rooti[ordm1] * rooti[ordm1];
    int ntrys, itcnt;
    for (ntrys = 0; ntrys < MAX_TRYS; ntrys++) {
      int found = FALSE;
      for (itcnt = 0; itcnt < MAX_ITS; itcnt++) {
        double lim = lim0 / (1.0 + std::fabs(p) + std::fabs(q));
        b[ord] = a[ord]; b[ordm1] = a[ordm1] - p * b[ord];
        c[ord] = b[ord]; c[ordm1] = b[ordm1] - p * c[ord];
        int k;
        for (k = 2; k <= ordm1; k++) {
          int mmk = ord - k;
          b[mmk] = a[mmk] - p * b[mmk + 1] - q * b[mmk + 2];
          c[mmk] = b[mmk] - p * c[mmk + 1] - q * c[mmk + 2];
          if (b[mmk] > lim || c[mmk] > lim) break;
        }
        if (k > ordm1) {
          b[0] = a[0] - p * b[1] - q * b[2];
          if (b[0] <= lim) k++;
        }
        if (k <= ord) break;

        double err = std::fabs(b[0]) + std::fabs(b[1]);
        if (err <= MAX_ERR) { found = TRUE; break; }

        double den = c[2] * c[2] - c[3] * (c[1] - b[1]);
        if (den == 0.0) break;
        double delp = (c[2] * b[1] - c[3] * b[0]) / den;
        double delq = (c[2] * b[0] - b[1] * (c[1] - b[1])) / den;
        p += delp; q += delq;
      }
      if (found) break;
      p = ((double)rand() - 0.5 * RAND_MAX) / (double)RAND_MAX;
      q = ((double)rand() - 0.5 * RAND_MAX) / (double)RAND_MAX;
    }
    if (itcnt >= MAX_ITS && ntrys >= MAX_TRYS) return FALSE;

    if (!qquad(1.0, p, q,
               &rootr[ordm1], &rooti[ordm1], &rootr[ordm2], &rooti[ordm2]))
      return FALSE;
    for (int i = 0; i <= ordm2; i++) a[i] = b[i + 2];
  }

  if (order >= 2) {
    if (order == 2) {
      return qquad(a[2], a[1], a[0], &rootr[1], &rooti[1], &rootr[0], &rooti[0]);
    }
  }
  if (order == 1) {
    if (a[1] != 0.0) rootr[0] = -a[0] / a[1];
    else rootr[0] = 100.0;
    rooti[0] = 0.0;
    return TRUE;
  }
  return TRUE;
}

/* ================================================================ */
/*  formant(): find LPC polynomial roots, convert to freq/band     */
/* ================================================================ */
static int formant_find(int lpc_order, double s_freq, double *lpca,
                        int *n_form, double *freq, double *band, int init) {
  static double rr[MAXORDER], ri[MAXORDER];

  if (init) {
    double x = M_PI / (lpc_order + 1);
    for (int i = 0; i <= lpc_order; i++) {
      double flo = lpc_order - i;
      rr[i] = 2.0 * std::cos((flo + 0.5) * x);
      ri[i] = 2.0 * std::sin((flo + 0.5) * x);
    }
  }

  /* Make a copy of lpca since lbpoly modifies its input */
  double a_copy[MAXORDER + 1];
  std::memcpy(a_copy, lpca, (lpc_order + 1) * sizeof(double));

  if (!lbpoly(a_copy, lpc_order, rr, ri)) {
    *n_form = 0;
    return FALSE;
  }

  double pi2t = M_PI * 2.0 / s_freq;
  int fc = 0;
  for (int ii = 0; ii < lpc_order; ii++) {
    if (rr[ii] != 0.0 || ri[ii] != 0.0) {
      double theta = std::atan2(ri[ii], rr[ii]);
      freq[fc] = std::fabs(theta / pi2t);
      band[fc] = 0.5 * s_freq *
                 std::log(rr[ii] * rr[ii] + ri[ii] * ri[ii]) / M_PI;
      if (band[fc] < 0.0) band[fc] = -band[fc];
      fc++;
      if (ii + 1 < lpc_order &&
          rr[ii] == rr[ii + 1] && ri[ii] == -ri[ii + 1] && ri[ii] != 0.0)
        ii++; /* skip conjugate */
    }
  }

  /* Order poles by frequency, push real poles to end */
  double fold = s_freq / 2.0;
  for (int i = 0; i < fc - 1; i++) {
    for (int ii = 0; ii < fc - 1 - i; ii++) {
      int iscomp1 = (freq[ii] > 1.0) && (freq[ii] < fold);
      int iscomp2 = (freq[ii + 1] > 1.0) && (freq[ii + 1] < fold);
      int swit = (freq[ii] > freq[ii + 1]) && iscomp2;
      if (swit || (iscomp2 && !iscomp1)) {
        double tmp = band[ii + 1]; band[ii + 1] = band[ii]; band[ii] = tmp;
        tmp = freq[ii + 1]; freq[ii + 1] = freq[ii]; freq[ii] = tmp;
      }
    }
  }

  int ii = 0;
  double th = fold - 1.0;
  for (int i = 0; i < fc; i++)
    if (freq[i] > 1.0 && freq[i] < th) ii++;
  *n_form = ii;
  return TRUE;
}

/* ================================================================ */
/*  jkFormant.c functions — DP formant tracking                    */
/* ================================================================ */

#define MAXCAN 300

static double MISSING = 1.0, NOBAND = 1000.0;
static double DF_FACT = 20.0, DFN_FACT = 0.3, BAND_FACT = 0.002;
static double F_BIAS = 0.0, F_MERGE = 2000.0;

static double *fre_ptr;
static double fnom[] = {500, 1500, 2500, 3500, 4500, 5500, 6500};
static double fmins[] = {50, 400, 1000, 2000, 2000, 3000, 3000};
static double fmaxs[] = {1500, 3500, 4500, 5000, 6000, 6000, 8000};
static int maxp, maxf, ncan_var;
static int domerge_var = TRUE;
static short **pc_ptr;

static int canbe(int pnumb, int fnumb) {
  return (fre_ptr[pnumb] >= fmins[fnumb]) && (fre_ptr[pnumb] <= fmaxs[fnumb]);
}

static void candy(int cand, int pnumb, int fnumb) {
  if (fnumb < maxf) pc_ptr[cand][fnumb] = -1;
  if (pnumb < maxp && fnumb < maxf) {
    if (canbe(pnumb, fnumb)) {
      pc_ptr[cand][fnumb] = pnumb;
      if (domerge_var && fnumb == 0 && canbe(pnumb, fnumb + 1)) {
        ncan_var++;
        pc_ptr[ncan_var][0] = pc_ptr[cand][0];
        candy(ncan_var, pnumb, fnumb + 1);
      }
      candy(cand, pnumb + 1, fnumb + 1);
      if ((pnumb + 1) < maxp && canbe(pnumb + 1, fnumb)) {
        ncan_var++;
        for (int i = 0; i < fnumb; i++)
          pc_ptr[ncan_var][i] = pc_ptr[cand][i];
        candy(ncan_var, pnumb + 1, fnumb);
      }
    } else {
      candy(cand, pnumb + 1, fnumb);
    }
  }
  if (pnumb >= maxp && fnumb < maxf - 1 && pc_ptr[cand][fnumb] < 0) {
    int i, j;
    if (fnumb) {
      j = fnumb - 1;
      while (j > 0 && pc_ptr[cand][j] < 0) j--;
      i = ((j = pc_ptr[cand][j]) >= 0) ? j : 0;
    } else i = 0;
    candy(cand, i, fnumb + 1);
  }
}

static void get_fcand(int npole, double *freq, double *band,
                      int nform, short **pcan) {
  (void)band;
  ncan_var = 0;
  pc_ptr = pcan;
  fre_ptr = freq;
  maxp = npole;
  maxf = nform;
  candy(ncan_var, 0, 0);
  ncan_var++;
}

static void set_nominal_freqs(double f1) {
  for (int i = 0; i < MAXFORMANTS; i++) {
    fnom[i] = ((i * 2) + 1) * f1;
    fmins[i] = fnom[i] - ((i + 1) * f1) + 50.0;
    fmaxs[i] = fnom[i] + (i * f1) + 1000.0;
  }
}

static double get_stat_max(POLE **pole, int nframes) {
  double amax = pole[0]->rms;
  for (int i = 1; i < nframes; i++)
    if (pole[i]->rms > amax) amax = pole[i]->rms;
  return amax;
}

static double integerize(double time, double freq) {
  int i = (int)(0.5 + freq * time);
  return ((double)i) / freq;
}

/* ================================================================ */
/*  lpc_poles_array: adapted from lpc_poles(), no Sound*            */
/* ================================================================ */
static POLE **lpc_poles_array(short *data, int length, int samprate,
                              double wdur, double frame_int, int lpc_ord,
                              double preemp, int lpc_type, int w_type,
                              int *out_nframes) {
  if (lpc_type == 1) {
    wdur = 0.025;
    preemp = std::exp(-62.831853 * 90.0 / samprate);
  }
  if (lpc_ord > MAXORDER || lpc_ord < 2) return NULL;

  wdur = integerize(wdur, (double)samprate);
  frame_int = integerize(frame_int, (double)samprate);
  int nfrm = 1 + (int)((((double)length / samprate) - wdur) / frame_int);
  if (nfrm < 1) return NULL;

  int size = (int)(0.5 + wdur * samprate);
  int step = (int)(0.5 + frame_int * samprate);
  POLE **pole = (POLE **)ckalloc(nfrm * sizeof(POLE *));

  short *datap = data;
  int init = TRUE;
  for (int j = 0; j < nfrm; j++, datap += step) {
    pole[j] = (POLE *)ckalloc(sizeof(POLE));
    double *frp = (double *)ckalloc(sizeof(double) * lpc_ord);
    double *bap = (double *)ckalloc(sizeof(double) * lpc_ord);
    pole[j]->freq = frp;
    pole[j]->band = bap;

    double energy = 0.0;
    double lpca[MAXORDER + 1];
    double normerr;

    switch (lpc_type) {
    case 0:
      if (!lpc(lpc_ord, 70.0, size, datap, lpca, NULL, NULL,
               &normerr, &energy, preemp, w_type))
        break;
      break;
    case 1:
      if (!lpcbsa(lpc_ord, 70.0, size, datap, lpca, NULL, NULL,
                  &normerr, &energy, preemp))
        break;
      break;
    case 2: {
      int Ord = lpc_ord;
      double alpha, r0;
      w_covar(datap, &Ord, size, 0, lpca, &alpha, &r0, preemp, 0);
      if (Ord != lpc_ord || alpha <= 0.0) { /* problem */ }
      energy = std::sqrt(r0 / (size - Ord));
    } break;
    }

    pole[j]->change = 0.0;
    if ((pole[j]->rms = energy) > 1.0) {
      int nform;
      formant_find(lpc_ord, (double)samprate, lpca, &nform, frp, bap, init);
      pole[j]->npoles = nform;
      init = FALSE;
    } else {
      pole[j]->npoles = 0;
      init = TRUE;
    }
  }

  *out_nframes = nfrm;
  return pole;
}

/* ================================================================ */
/*  dpform_array: adapted from dpform(), returns arrays directly   */
/* ================================================================ */
static int dpform_array(POLE **pole, int nframes, int samprate,
                        int nform, double nom_f1,
                        double **out_fr, double **out_ba) {
  if (!pole || nframes < 1) return FALSE;

  if (nom_f1 > 0.0) set_nominal_freqs(nom_f1);

  double rmsmax = get_stat_max(pole, nframes);
  double FBIAS_val = F_BIAS / (0.01 * samprate);
  double dffact = (DF_FACT * 0.01) * samprate;
  double bfact = BAND_FACT / (0.01 * samprate);
  double ffact = DFN_FACT / (0.01 * samprate);
  double merge_cost = F_MERGE;
  domerge_var = (merge_cost <= 1000.0) ? TRUE : FALSE;
  /* Note: original code checks merge_cost > 1000 to disable;
     default F_MERGE=2000 means domerge=FALSE. Keep original logic. */
  if (merge_cost > 1000.0) domerge_var = FALSE;

  /* Allocate frequency and bandwidth result arrays */
  double **fr = (double **)ckalloc(sizeof(double *) * nform * 2);
  double **ba = fr + nform;
  for (int i = 0; i < nform * 2; i++)
    fr[i] = (double *)ckalloc(sizeof(double) * nframes);

  /* Candidate array */
  short **pcan = (short **)ckalloc(sizeof(short *) * MAXCAN);
  for (int i = 0; i < MAXCAN; i++)
    pcan[i] = (short *)ckalloc(sizeof(short) * nform);

  /* DP lattice */
  FORM **fl = (FORM **)ckalloc(sizeof(FORM *) * nframes);
  for (int i = 0; i < nframes; i++)
    fl[i] = (FORM *)ckalloc(sizeof(FORM));

  /* Main tracking loop */
  for (int i = 0; i < nframes; i++) {
    ncan_var = 0;
    double rmsdffact = (pole[i]->rms / rmsmax) * dffact;

    if (pole[i]->npoles) {
      get_fcand(pole[i]->npoles, pole[i]->freq, pole[i]->band, nform, pcan);
      fl[i]->prept  = (short *)ckalloc(sizeof(short) * ncan_var);
      fl[i]->cumerr = (double *)ckalloc(sizeof(double) * ncan_var);
      fl[i]->cand   = (short **)ckalloc(sizeof(short *) * ncan_var);
      for (int j = 0; j < ncan_var; j++) {
        fl[i]->cand[j] = (short *)ckalloc(sizeof(short) * nform);
        for (int k = 0; k < nform; k++)
          fl[i]->cand[j][k] = pcan[j][k];
      }
    }
    fl[i]->ncand = ncan_var;

    for (int j = 0; j < ncan_var; j++) {
      double minerr = 0.0;
      int mincan = -1;
      if (i) {
        if (fl[i - 1]->ncand) minerr = 2.0e30;
        for (int k = 0; k < fl[i - 1]->ncand; k++) {
          double pferr = 0.0;
          for (int l = 0; l < nform; l++) {
            int ic = fl[i]->cand[j][l];
            int ip = fl[i - 1]->cand[k][l];
            if (ic >= 0 && ip >= 0) {
              double ftemp = 2.0 * std::fabs(pole[i]->freq[ic] - pole[i - 1]->freq[ip]) /
                             (pole[i]->freq[ic] + pole[i - 1]->freq[ip]);
              pferr += ftemp * ftemp;
            } else
              pferr += MISSING;
          }
          double conerr = rmsdffact * pferr + fl[i - 1]->cumerr[k];
          if (conerr < minerr) { minerr = conerr; mincan = k; }
        }
      }

      fl[i]->prept[j] = mincan;

      double berr = 0.0, ferr = 0.0, fbias = 0.0, merger = 0.0;
      for (int k = 0; k < nform; k++) {
        int ic = fl[i]->cand[j][k];
        if (ic >= 0) {
          if (k == 0) {
            double ftemp = pole[i]->freq[ic];
            merger = (domerge_var && nform > 1 &&
                      ftemp == pole[i]->freq[fl[i]->cand[j][1]])
                         ? merge_cost : 0.0;
          }
          berr += pole[i]->band[ic];
          ferr += std::fabs(pole[i]->freq[ic] - fnom[k]) / fnom[k];
          fbias += pole[i]->freq[ic];
        } else {
          fbias += fnom[k];
          berr += NOBAND;
          ferr += MISSING;
        }
      }
      fl[i]->cumerr[j] = FBIAS_val * fbias + bfact * berr +
                          merger + ffact * ferr + minerr;
    }
  }

  /* Backtrace: pick lowest-cost path */
  int mincan = -1;
  for (int i = nframes - 1; i >= 0; i--) {
    if (mincan < 0 && fl[i]->ncand) {
      double minerr = fl[i]->cumerr[0];
      mincan = 0;
      for (int j = 1; j < fl[i]->ncand; j++)
        if (fl[i]->cumerr[j] < minerr) { minerr = fl[i]->cumerr[j]; mincan = j; }
    }
    if (mincan >= 0) {
      for (int j = 0; j < nform; j++) {
        int k = fl[i]->cand[mincan][j];
        if (k >= 0) {
          fr[j][i] = pole[i]->freq[k];
          ba[j][i] = pole[i]->band[k];
        } else {
          if (i < nframes - 1) {
            fr[j][i] = fr[j][i + 1];
            ba[j][i] = ba[j][i + 1];
          } else {
            fr[j][i] = fnom[j];
            ba[j][i] = NOBAND;
          }
        }
      }
      mincan = fl[i]->prept[mincan];
    } else {
      for (int j = 0; j < nform; j++) {
        fr[j][i] = fnom[j];
        ba[j][i] = NOBAND;
      }
    }
  }

  /* Cleanup DP lattice */
  for (int i = nframes - 1; i >= 0; i--) {
    if (fl[i]->ncand && fl[i]->cand) {
      for (int j = 0; j < fl[i]->ncand; j++) ckfree(fl[i]->cand[j]);
      ckfree(fl[i]->cand);
      ckfree(fl[i]->cumerr);
      ckfree(fl[i]->prept);
    }
  }
  for (int i = 0; i < nframes; i++) ckfree(fl[i]);
  ckfree(fl);

  /* Cleanup pole data */
  for (int i = 0; i < nframes; i++) {
    ckfree(pole[i]->freq);
    ckfree(pole[i]->band);
    ckfree(pole[i]);
  }
  ckfree(pole);

  /* Cleanup candidate array */
  for (int i = 0; i < MAXCAN; i++) ckfree(pcan[i]);
  ckfree(pcan);

  /* Copy results to output arrays (nform freqs then nform bands) */
  for (int j = 0; j < nform; j++) {
    out_fr[j] = fr[j];    /* caller takes ownership */
    out_ba[j] = ba[j];    /* fr+nform = ba pointer */
  }
  /* fr[nform .. 2*nform-1] are the bandwidth arrays = ba[0..nform-1] */
  /* They're already assigned via out_ba[j] = fr[nform+j] */
  /* Free the pointer array but NOT the data arrays */
  ckfree(fr);  /* This only frees the pointer array, not data pointed to */

  return TRUE;
}

/* ================================================================ */
/*  Downsampling + highpass (adapted from jkFormant.c)             */
/* ================================================================ */

static int lc_lin_fir(double fc, int *nf, double *coef) {
  if ((*nf % 2) != 1 || *nf > 127) {
    if (*nf <= 126) *nf += 1; else *nf = 127;
  }
  int n = (*nf + 1) / 2;
  double fn = M_PI * 2.0 * fc;
  coef[0] = 2.0 * fc;
  for (int i = 1; i < n; i++)
    coef[i] = std::sin(i * fn) / (M_PI * i);
  fn = M_PI * 2.0 / ((double)(*nf - 1));
  for (int i = 0; i < n; i++)
    coef[i] *= (0.5 + 0.5 * std::cos(fn * (double)i));
  return TRUE;
}

static void do_fir(short *buf, int in_samps, short *bufo,
                   int ncoef, short *ic, int invert) {
  short co[256], mem[256];
  int k = (ncoef * 2) - 1;
  /* Build full symmetric filter */
  for (int i = ncoef - 1, j = 0, m = k - 1; i > 0; i--, j++, m--) {
    if (!invert) co[j] = co[m] = ic[i];
    else { co[j] = co[m] = -ic[i]; }
  }
  if (!invert) co[ncoef - 1] = ic[0];
  else {
    int integral = 0;
    for (int i = 1; i < ncoef; i++) integral += ic[i];
    integral *= 2;
    integral += ic[0];
    co[ncoef - 1] = integral - ic[0];
  }

  /* Initialize memory */
  for (int i = 0; i < ncoef - 1; i++) mem[i] = 0;
  for (int i = 0; i < ncoef; i++) mem[ncoef - 1 + i] = buf[i];

  int l = 16384, m = 15;
  short *inp = buf + ncoef;
  for (int i = 0; i < in_samps - ncoef; i++) {
    int sum = 0;
    for (int j = 0; j < k; j++)
      sum += ((co[j] * mem[j]) + l) >> m;
    /* shift memory */
    for (int j = 0; j < k - 1; j++) mem[j] = mem[j + 1];
    mem[k - 1] = *inp++;
    *bufo++ = (short)sum;
  }
  for (int i = 0; i < ncoef; i++) {
    int sum = 0;
    for (int j = 0; j < k; j++)
      sum += ((co[j] * mem[j]) + l) >> m;
    for (int j = 0; j < k - 1; j++) mem[j] = mem[j + 1];
    mem[k - 1] = 0;
    *bufo++ = (short)sum;
  }
}

static int get_abs_maximum(short *d, int n) {
  short amax = (d[0] >= 0) ? d[0] : -d[0];
  for (int i = 1; i < n; i++) {
    short t = (d[i] >= 0) ? d[i] : -d[i];
    if (t > amax) amax = t;
  }
  return (int)amax;
}

static int ratprx(double a, int *k, int *l, int qlim) {
  double aa = std::fabs(a);
  int ai = (int)aa;
  double af = aa - ai;
  double em = 1.0;
  double pp = 0, qq = 0;
  for (int q = 1; q <= qlim; q++) {
    double ps = q * af;
    int ip = (int)(ps + 0.5);
    double e = std::fabs((ps - (double)ip) / q);
    if (e < em) { em = e; pp = ip; qq = q; }
  }
  *k = (int)(ai * qq + pp);
  if (a < 0) *k = -(*k);
  *l = (int)qq;
  return TRUE;
}

static int dwnsamp(short *buf, int in_samps, short **buf2,
                   int *out_samps, int insert, int decimate,
                   int ncoef, short *ic, int *smin, int *smax) {
  short *buft = (short *)ckalloc(sizeof(short) * insert * in_samps);
  if (!buft) return FALSE;
  *buf2 = buft;

  int k = get_abs_maximum(buf, in_samps);
  if (k == 0) k = 1;
  if (insert > 1) k = (32767 * 32767) / k;
  else k = (16384 * 32767) / k;
  int l = 16384, m = 15;

  for (int i = 0; i < in_samps; i++) {
    *buft++ = (short)((k * (int)buf[i] + l) >> m);
    for (int j = 1; j < insert; j++) *buft++ = 0;
  }
  buft = *buf2; /* reset to start */

  do_fir(buft, in_samps * insert, buft, ncoef, ic, 0);

  int j = (in_samps * insert) / decimate;
  *out_samps = j;
  int imax, imin;
  short *src = *buf2;
  short *dst = *buf2;
  imax = imin = *src;
  for (int i = 0; i < j; i++, src += decimate) {
    *dst++ = *src;
    if (*src > imax) imax = *src;
    if (*src < imin) imin = *src;
  }
  *smin = imin; *smax = imax;
  *buf2 = (short *)ckrealloc(*buf2, sizeof(short) * (*out_samps));
  return TRUE;
}

static short *fdownsample_array(short *data, int length, int samprate,
                                double freq2, int *out_length, int *out_rate) {
  static double beta = 0.0, b[256];
  static int ncoeff = 127, ncoefft = 0, nbits = 15;
  static short ic[256];
  int insert, decimate, out_samps, smin, smax;

  double freq1 = (double)samprate;
  double ratio = freq2 / freq1;
  ratprx(ratio, &insert, &decimate, 10);
  double ratio_t = (double)insert / (double)decimate;

  if (ratio_t > 0.99) {
    /* No downsampling needed — copy input */
    short *out = (short *)ckalloc(sizeof(short) * length);
    std::memcpy(out, data, sizeof(short) * length);
    *out_length = length;
    *out_rate = samprate;
    return out;
  }

  freq2 = ratio_t * freq1;
  double beta_new = (0.5 * freq2) / (insert * freq1);
  if (beta != beta_new) {
    beta = beta_new;
    lc_lin_fir(beta, &ncoeff, b);
    double maxi = (1 << nbits) - 1;
    int j = (ncoeff / 2) + 1;
    ncoefft = 0;
    for (int i = 0; i < j; i++) {
      ic[i] = (short)(0.5 + maxi * b[i]);
      if (ic[i]) ncoefft = i + 1;
    }
  }

  short *bufout = NULL;
  if (dwnsamp(data, length, &bufout, &out_samps, insert, decimate,
              ncoefft, ic, &smin, &smax)) {
    *out_length = out_samps;
    *out_rate = (int)freq2;
    return bufout;
  }
  *out_length = 0; *out_rate = 0;
  return NULL;
}

#define LCSIZ 101

static short *highpass_array(short *data, int length, int samprate,
                             int *out_length) {
  static short *lcf = NULL;
  static int len = 0;

  short *dataout = (short *)ckalloc(sizeof(short) * length);

  if (!len) {
    lcf = (short *)ckalloc(sizeof(short) * LCSIZ);
    len = 1 + (LCSIZ / 2);
    double fn = M_PI * 2.0 / (LCSIZ - 1);
    double scale = 32767.0 / (0.5 * LCSIZ);
    for (int i = 0; i < len; i++)
      lcf[i] = (short)(scale * (0.5 + 0.4 * std::cos(fn * (double)i)));
  }
  do_fir(data, length, dataout, len, lcf, 1);
  *out_length = length;
  (void)samprate;
  return dataout;
}

/* ================================================================ */
/*  Public API: compute_formants()                                 */
/* ================================================================ */

int compute_formants(const double *audio, int n_samples, int sample_rate,
                     int n_formants, int lpc_order, double window_dur,
                     double frame_interval, double preemphasis,
                     int lpc_type, int window_type, double ds_freq,
                     double nom_f1,
                     double **out_freqs, double **out_bands,
                     int *out_n_frames) {

  if (!audio || n_samples < 1 || sample_rate < 1) return FALSE;
  if (n_formants < 1 || n_formants > MAXFORMANTS) return FALSE;
  if (n_formants > (lpc_order - 4) / 2) return FALSE;

  /* Convert double audio to short (scale by 32767 if needed) */
  short *sdata = (short *)ckalloc(sizeof(short) * n_samples);
  double amax = 0.0;
  for (int i = 0; i < n_samples; i++) {
    double v = std::fabs(audio[i]);
    if (v > amax) amax = v;
  }
  /* If audio is normalized [-1,1], scale up; if INT16 range, keep as-is */
  double scale = (amax <= 1.0) ? 32767.0 : (amax < 32767.0 ? 32767.0 / amax : 1.0);
  for (int i = 0; i < n_samples; i++) {
    double v = audio[i] * scale;
    if (v > 32767.0) v = 32767.0;
    if (v < -32768.0) v = -32768.0;
    sdata[i] = (short)v;
  }

  int cur_rate = sample_rate;
  short *cur_data = sdata;
  int cur_length = n_samples;
  short *ds_data = NULL;
  short *hp_data = NULL;

  /* Downsample if needed */
  if (ds_freq < cur_rate) {
    int ds_len, ds_rate;
    ds_data = fdownsample_array(cur_data, cur_length, cur_rate,
                                ds_freq, &ds_len, &ds_rate);
    if (ds_data) {
      cur_data = ds_data;
      cur_length = ds_len;
      cur_rate = ds_rate;
    }
  }

  /* Highpass filter */
  if (preemphasis < 1.0) {
    int hp_len;
    hp_data = highpass_array(cur_data, cur_length, cur_rate, &hp_len);
    if (hp_data) {
      cur_data = hp_data;
      cur_length = hp_len;
    }
  }

  /* LPC pole analysis */
  int nfrm;
  POLE **poles = lpc_poles_array(cur_data, cur_length, cur_rate,
                                 window_dur, frame_interval, lpc_order,
                                 preemphasis, lpc_type, window_type, &nfrm);
  if (!poles) {
    ckfree(sdata);
    if (ds_data) ckfree(ds_data);
    if (hp_data) ckfree(hp_data);
    return FALSE;
  }

  /* DP formant tracking */
  /* Allocate arrays for results: out_freqs[0..n_formants-1], out_bands[0..n_formants-1] */
  if (!dpform_array(poles, nfrm, cur_rate, n_formants, nom_f1,
                    out_freqs, out_bands)) {
    ckfree(sdata);
    if (ds_data) ckfree(ds_data);
    if (hp_data) ckfree(hp_data);
    return FALSE;
  }

  *out_n_frames = nfrm;

  ckfree(sdata);
  if (ds_data) ckfree(ds_data);
  if (hp_data) ckfree(hp_data);
  return TRUE;
}

} /* namespace snack_formant */
