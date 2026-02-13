// Rcpp wrapper for MOMEL C algorithm
// Wraps cible(), reduc(), borne() from momel.c

#include <Rcpp.h>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <cassert>

using namespace Rcpp;

// --- Constants from momel.h ---
#define PAS 10.
#define NVALM 100000
#define SEUILV 50.
#define PARMAX (NVALM/2)
#define FAUX 0
#define VRAI 1
#define FSIGMA 1.
#define PAS_TRAME 10.
#define ECART_MIN 50.
#define RAPP_MIN_FREQ 0.05
#define RAPP_GLITCH 0.05
#define HALO_BORNE_TRAME 4

// --- Structs ---
struct st_cib {
  float x;
  float y;
};

struct st_cibred {
  float x;
  float y;
  int   p;
};

// --- eliminer_glitch ---
// hzptr is offset into padded buffer, so hzptr[-1] is valid (=0)
static void eliminer_glitch(int nval, float *hzptr) {
  for (int i = 0; i < nval; i++) {
    if ((hzptr[i] > hzptr[i-1] * (1 + RAPP_GLITCH))
        && (hzptr[i] > hzptr[i+1] * (1 + RAPP_GLITCH))) {
      hzptr[i] = 0.0;
    }
  }
}

// --- calcrgp ---
static int calcrgp(float *pond, int dpx, int fpx, float *hzptr,
                   float *pa0, float *pa1, float *pa2, float *hzes) {
  double pn = 0.;
  double sx = 0., sx2 = 0., sx3 = 0., sx4 = 0.;
  double sy = 0., sxy = 0., sx2y = 0.;

  for (int ix = dpx; ix <= fpx; ix++) {
    const double p = pond[ix];
    if (p != 0.) {
      const double val_ix = (double)ix;
      const double y = hzptr[ix];
      const double x2 = val_ix * val_ix;
      const double x3 = x2 * val_ix;
      const double x4_v = x2 * x2;
      const double xy = val_ix * y;
      const double x2y_v = x2 * y;
      pn += p;
      sx += (p * val_ix);
      sx2 += p * x2;
      sx3 += p * x3;
      sx4 += p * x4_v;
      sy += p * y;
      sxy += p * xy;
      sx2y += p * x2y_v;
    }
  }
  if (pn < 3.) return 1;

  double spdxy  = sxy  - (sx * sy) / pn;
  double spdx2  = sx2  - (sx * sx) / pn;
  double spdx3  = sx3  - (sx * sx2) / pn;
  double spdx4  = sx4  - (sx2 * sx2) / pn;
  double spdx2y = sx2y - (sx2 * sy) / pn;
  double muet = spdx2 * spdx4 - spdx3 * spdx3;
  if (spdx2 == 0. || muet == 0.) return 1;

  *pa2 = (float)((spdx2y * spdx2 - spdxy * spdx3) / muet);
  *pa1 = (float)((spdxy - *pa2 * spdx3) / spdx2);
  *pa0 = (float)((sy - *pa1 * sx - *pa2 * sx2) / pn);

  for (int ix = dpx; ix <= fpx; ix++)
    hzes[ix] = *pa0 + (*pa1 + *pa2 * (float)ix) * (float)ix;
  return 0;
}

// --- cible ---
static int cible(int nval, float *hzptr, int lfen1,
                 float maxec, float hzinf, float hzsup,
                 struct st_cib *cib) {
  // Use heap allocation for large arrays
  float *bpond = (float*)calloc(100 + NVALM + 100, sizeof(float));
  float *bhzes = (float*)calloc(100 + NVALM + 100, sizeof(float));
  float *pond = bpond + 100;
  float *hzes = bhzes + 100;

  for (int ix = 0; ix < nval; ix++) {
    if (hzptr[ix] > SEUILV)
      pond[ix] = 1.0;
  }

  int lfens2 = lfen1 / 2;
  float *bpondloc = (float*)calloc(1000 + NVALM + 1000, sizeof(float));
  float *pondloc = bpondloc + 1000;
  for (int ix = 0; ix < nval; ix++) {
    int dpx = ix - lfens2;
    int fpx = dpx + lfen1;
    int nsup = 0, nsupr = -1;

    // Local copy of weights — reset and copy
    for (int i = dpx; i <= fpx; i++)
      pondloc[i] = pond[i];

    float a0, a1, a2;
    int retour_rgp = 1;
    while (nsup > nsupr) {
      nsupr = nsup;
      nsup = 0;
      retour_rgp = calcrgp(pondloc, dpx, fpx, hzptr, &a0, &a1, &a2, hzes);
      if (retour_rgp != 0) break;
      for (int x = dpx; x <= fpx; x++) {
        if (hzptr[x] == 0. || hzes[x] / hzptr[x] > maxec) {
          pondloc[x] = 0;
          nsup++;
        }
      }
    }

    float xc = 0., yc = 0.;
    if (retour_rgp == 0 && a2 != 0.) {
      float vxc = -a1 / (a2 + a2);
      if ((vxc > ix - lfen1) && (vxc < ix + lfen1)) {
        float vyc = a0 + (a1 + a2 * vxc) * vxc;
        if (vyc > hzinf && vyc < hzsup) {
          xc = vxc;
          yc = vyc;
        }
      }
    }
    cib[ix].x = xc;
    cib[ix].y = yc;
  }
  free(bpondloc);
  free(bpond);
  free(bhzes);
  return 0;
}

// --- cb_compare for qsort ---
static int cb_compare(const void *va, const void *vb) {
  struct st_cib *pa = (struct st_cib *)va;
  struct st_cib *pb = (struct st_cib *)vb;
  if (pa->x > pb->x) return 1;
  if (pa->x < pb->x) return -1;
  return 0;
}

// --- reduc ---
static int reduc(int nval, int lfen2,
                 float seuildiff_x, float seuilrapp_y,
                 struct st_cib *cib,
                 int *pnred2,
                 struct st_cibred *cibred2) {
  struct st_cibred *cibred = (struct st_cibred*)malloc(sizeof(struct st_cibred) * PARMAX);
  float *xdist = (float*)calloc(NVALM, sizeof(float));
  float *ydist = (float*)calloc(NVALM, sizeof(float));
  float *dist  = (float*)calloc(NVALM, sizeof(float));
  int *xd = (int*)calloc(PARMAX + 1, sizeof(int));

  struct st_cibred cibred_cour;
  int lf = lfen2 / 2;
  float xds = 0., yds = 0.;
  int np = 0;

  for (int i = 0; i < nval - 1; i++) {
    int j1 = 0;
    if (i > lf) j1 = i - lf;
    int j2 = nval - 1;
    if (i + lf < nval - 1) j2 = i + lf;
    float sxg = 0., syg = 0.;
    int ng = 0;
    for (int j = j1; j < i + 1; j++) {
      if (cib[j].y > SEUILV) {
        sxg += cib[j].x;
        syg += cib[j].y;
        ng++;
      }
    }
    float sxd = 0., syd = 0.;
    int nd = 0;
    for (int j = i + 1; j < j2; j++) {
      if (cib[j].y > SEUILV) {
        sxd += cib[j].x;
        syd += cib[j].y;
        nd++;
      }
    }
    xdist[i] = ydist[i] = -1.;
    if (nd * ng > 0) {
      xdist[i] = fabs(sxg / ng - sxd / nd);
      ydist[i] = fabs(syg / ng - syd / nd);
      xds += xdist[i];
      yds += ydist[i];
      np++;
    }
  }

  float px = np / xds;
  float py = np / yds;

  for (int i = 0; i < nval; i++) {
    dist[i] = -1;
    if (xdist[i] > 0.)
      dist[i] = (xdist[i] * px + ydist[i] * py) / (px + py);
  }

  float seuil = 2. / (px + py);
  xd[0] = 0;
  int ncible = 0;
  int susseuil = FAUX;
  int xmax = 0;

  for (int i = 0; i < nval; i++) {
    if (ncible == PARMAX) {
      free(cibred); free(xdist); free(ydist); free(dist); free(xd);
      return 1;
    }
    if (susseuil == FAUX) {
      if (dist[i] > seuil) {
        susseuil = VRAI;
        xmax = i;
      }
    }
    if (susseuil == VRAI) {
      if (dist[i] > dist[xmax]) xmax = i;
      if (dist[i] < seuil) {
        ncible++;
        xd[ncible] = xmax;
        susseuil = FAUX;
      }
    }
  }
  if (susseuil == VRAI) {
    ncible++;
    xd[ncible] = xmax;
  }
  if (xmax < nval) {
    ncible++;
    xd[ncible] = nval;
  }

  // Partition on x
  int ncibr = -1;
  for (int ip = 0; ip < ncible; ip++) {
    int parinf = xd[ip];
    int parsup = xd[ip + 1];
    for (int k = 0; k < 1; k++) {
      float sx = 0., sx2 = 0., sy = 0., sy2 = 0.;
      int n = 0;
      for (int j = parinf; j < parsup; j++) {
        if (cib[j].y > 0.) {
          sx += cib[j].x;
          sx2 += cib[j].x * cib[j].x;
          sy += cib[j].y;
          sy2 += cib[j].y * cib[j].y;
          n++;
        }
      }
      if (n > 1) {
        double xm = (double)sx / (double)n;
        double ym = (double)sy / (double)n;
        double varx = (double)sx2 / (double)n - xm * xm;
        double vary = (double)sy2 / (double)n - ym * ym;
        if (varx <= 0.) varx = 0.1;
        if (vary <= 0.) vary = 0.1;
        float et2x = FSIGMA * sqrt(varx);
        float et2y = FSIGMA * sqrt(vary);
        float seuilbx = xm - et2x;
        float seuilhx = xm + et2x;
        float seuilby = ym - et2y;
        float seuilhy = ym + et2y;
        for (int j = parinf; j < parsup; j++) {
          if (cib[j].y > 0. &&
              (cib[j].x < seuilbx || cib[j].x > seuilhx ||
               cib[j].y < seuilby || cib[j].y > seuilhy)) {
            cib[j].x = cib[j].y = 0.;
          }
        }
      }
    }

    // Recalculate means
    float sx = 0., sy = 0.;
    int n = 0;
    for (int j = parinf; j < parsup; j++) {
      if (cib[j].y > 0.) {
        sx += cib[j].x;
        sy += cib[j].y;
        n++;
      }
    }
    if (n > 0) {
      cibred_cour.x = sx / n;
      cibred_cour.y = sy / n;
      cibred_cour.p = n;
      if (ncibr < 0) {
        ncibr++;
        cibred[ncibr] = cibred_cour;
      } else {
        if (cibred_cour.x > cibred[ncibr].x) {
          ncibr++;
          cibred[ncibr] = cibred_cour;
        } else {
          if (cibred_cour.p > cibred[ncibr].p)
            cibred[ncibr] = cibred_cour;
        }
      }
    }
  }

  int ncibr_brut = ncibr + 1;
  qsort((char*)cibred, ncibr_brut, sizeof(struct st_cibred), cb_compare);

  *pnred2 = 0;
  cibred2[0] = cibred[0];

  for (int i = 1; i < ncibr_brut; i++) {
    if (cibred[i].x - cibred2[*pnred2].x < seuildiff_x) {
      if (fabs((double)(cibred[i].y - cibred2[*pnred2].y)) / cibred2[*pnred2].y < seuilrapp_y) {
        cibred2[*pnred2].x = (cibred2[*pnred2].x + cibred[i].x) / 2.;
        cibred2[*pnred2].y = (cibred2[*pnred2].y + cibred[i].y) / 2.;
        cibred2[*pnred2].p += cibred[i].p;
      } else {
        if (cibred2[*pnred2].p < cibred[i].p)
          cibred2[*pnred2] = cibred[i];
      }
    } else {
      (*pnred2)++;
      cibred2[*pnred2] = cibred[i];
    }
  }
  (*pnred2)++;

  free(cibred); free(xdist); free(ydist); free(dist); free(xd);
  return 0;
}

// --- borne: returns results via out_t, out_f vectors ---
static int borne_collect(int nval, int nred2,
                         const struct st_cibred *cibred2,
                         const float *hzptr,
                         std::vector<float> &out_t,
                         std::vector<float> &out_f) {
  out_t.clear();
  out_f.clear();

  // Find first voiced
  int premier_voise = 0;
  for (int i = 0; i < nval; i++) {
    if (hzptr[i] >= SEUILV) { premier_voise = i; break; }
  }

  // Left boundary
  if ((int)cibred2[0].x > premier_voise + HALO_BORNE_TRAME) {
    float ancre_x = cibred2[0].x;
    float ancre_y = cibred2[0].y;
    float sx2y = 0., sx4 = 0.;
    int j = 0;
    for (int i = (int)ancre_x; i >= 0; i--) {
      if (hzptr[i] > SEUILV) {
        float x2 = (float)j * (float)j;
        sx2y += x2 * (hzptr[i] - ancre_y);
        sx4 += x2 * x2;
      }
      j++;
    }
    if (sx4 != 0.) {
      float frontiere = (float)premier_voise;
      float a = sx2y / sx4;
      float bx = frontiere - (ancre_x - frontiere);
      float by = ancre_y + 2 * a * (ancre_x - frontiere) * (ancre_x - frontiere);
      out_t.push_back(bx * PAS_TRAME);
      out_f.push_back(by);
    }
  }

  // Original targets
  for (int i = 0; i < nred2; i++) {
    out_t.push_back(cibred2[i].x * PAS_TRAME);
    out_f.push_back(cibred2[i].y);
  }

  // Find last voiced
  int dernier_voise = nval - 1;
  for (int i = nval - 1; i >= 0; i--) {
    if (hzptr[i] >= SEUILV) { dernier_voise = i; break; }
  }

  // Right boundary
  if ((int)cibred2[nred2 - 1].x < dernier_voise - HALO_BORNE_TRAME) {
    float ancre_x = cibred2[nred2 - 1].x;
    float ancre_y = cibred2[nred2 - 1].y;
    float sx2y = 0., sx4 = 0.;
    int j = 0;
    for (int i = (int)ancre_x; i < nval; i++) {
      if (hzptr[i] > SEUILV) {
        float x2 = (float)j * (float)j;
        sx2y += x2 * (hzptr[i] - ancre_y);
        sx4 += x2 * x2;
      }
      j++;
    }
    if (sx4 != 0.) {
      float frontiere = (float)dernier_voise;
      float a = sx2y / sx4;
      float bx = frontiere + (frontiere - ancre_x);
      float by = ancre_y + 2. * a * (ancre_x - frontiere) * (ancre_x - frontiere);
      out_t.push_back(bx * PAS_TRAME);
      out_f.push_back(by);
    }
  }

  return 0;
}

// [[Rcpp::export]]
Rcpp::DataFrame momel_c(Rcpp::NumericVector f0_values,
                        int window_length,
                        double min_f0, double max_f0,
                        double max_error,
                        int reduced_window_length,
                        double minimal_distance,
                        double minimal_frequency_ratio) {
  int nval = f0_values.size();
  if (nval <= 0 || nval >= NVALM) {
    return Rcpp::DataFrame::create(
      Rcpp::Named("time") = Rcpp::NumericVector(0),
      Rcpp::Named("frequency") = Rcpp::NumericVector(0)
    );
  }

  // Allocate padded buffer like original C code
  float *bhz = (float*)calloc(100 + NVALM + 100, sizeof(float));
  float *hzptr = bhz + 100;

  // Copy input as-is (matching C binary behavior — no clamping)
  for (int i = 0; i < nval; i++) {
    hzptr[i] = (float)f0_values[i];
  }

  // Eliminate glitches
  eliminer_glitch(nval, hzptr);

  // Allocate cible output
  struct st_cib *cib = (struct st_cib*)malloc(sizeof(struct st_cib) * NVALM);
  struct st_cibred *cibred2 = (struct st_cibred*)malloc(sizeof(struct st_cibred) * PARMAX);

  int err = cible(nval, hzptr, window_length, (float)max_error,
                  (float)min_f0, (float)max_f0, cib);
  if (err) {
    free(bhz); free(cib); free(cibred2);
    return Rcpp::DataFrame::create(
      Rcpp::Named("time") = Rcpp::NumericVector(0),
      Rcpp::Named("frequency") = Rcpp::NumericVector(0)
    );
  }

  int nred2 = 0;
  err = reduc(nval, reduced_window_length,
              (float)minimal_distance, (float)minimal_frequency_ratio,
              cib, &nred2, cibred2);
  if (err || nred2 == 0) {
    free(bhz); free(cib); free(cibred2);
    return Rcpp::DataFrame::create(
      Rcpp::Named("time") = Rcpp::NumericVector(0),
      Rcpp::Named("frequency") = Rcpp::NumericVector(0)
    );
  }

  std::vector<float> out_t, out_f;
  err = borne_collect(nval, nred2, cibred2, hzptr, out_t, out_f);

  free(bhz); free(cib); free(cibred2);

  if (err || out_t.empty()) {
    return Rcpp::DataFrame::create(
      Rcpp::Named("time") = Rcpp::NumericVector(0),
      Rcpp::Named("frequency") = Rcpp::NumericVector(0)
    );
  }

  // Convert to R — borne output is already in ms (multiplied by PAS_TRAME)
  // but momel() R function expects frame-based time, so divide back
  Rcpp::NumericVector r_time(out_t.size()), r_freq(out_f.size());
  for (size_t i = 0; i < out_t.size(); i++) {
    r_time[i] = out_t[i] / PAS_TRAME;  // back to frame units
    r_freq[i] = out_f[i];
  }

  return Rcpp::DataFrame::create(
    Rcpp::Named("time") = r_time,
    Rcpp::Named("frequency") = r_freq
  );
}
