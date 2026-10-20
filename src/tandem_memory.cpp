/*
 * tandem_memory.cpp
 *
 * Memory-based interface for TANDEM pitch tracking
 * Adds in-memory processing without modifying original TANDEM source
 *
 * Copyright (C) 2025 superassp contributors
 * Based on TANDEM by G. Hu & D. L. Wang, Ohio State University
 */

#include <Rcpp.h>
#include <R.h>
#include <Rinternals.h>

// Include TANDEM headers only once through voicedMask.h
// (voicedMask.h includes pitch.h, which includes feature.h, etc.)
#include "tandem/tandem_64/voicedMask.h"
#include "tandem/tandem_64/gammaTone.h"

#include <vector>
#include <cmath>
#include <cstring>
#include <memory>
#include <new>

#ifdef _OPENMP
#include <omp.h>
#endif

// Memory-based audio processing (replacement for readInput)
double* processAudioBufferTandem(double *samples, int numSamples, double *scale_out) {
    // Calculate RMS for scaling to 60-dB loudness level
    double sumE = 0.0;
    for (int i = 0; i < numSamples; i++) {
        sumE += samples[i] * samples[i];
    }
    sumE = sqrt(sumE / numSamples);

    double scale = 1000.0 / (sumE + 1e-10);  // Avoid division by zero
    *scale_out = scale;

    // Allocate and scale audio
    double *signal = new double[numSamples];
    for (int i = 0; i < numSamples; i++) {
        signal[i] = scale * samples[i];
    }

    return signal;
}

// TANDEM analysis parameters (shared by the constructors below and by the
// feature stage, which must build an identical envelope band-pass filter).
static featurePara tandemPara(void) {
    featurePara x;
    gammaTonePara y;

    y.lCF = 50;      // Lower center frequency
    y.rCF = 8000;    // Upper center frequency
    y.nChan = 64;    // 64-channel filterbank
    y.sf = 20000;    // 20 kHz sample rate

    x.gtP = y;
    x.bP1 = 50;      // Min pitch period
    x.bP2 = 450;     // Max pitch period
    x.bPTs = 20;     // Pitch threshold
    x.theta_p = 0.5; // Pitch probability threshold
    return x;
}

// Initialize TANDEM voicedMask and gammaTone filterbank
void initVoicedMaskTandem(gammaToneFilterBank *&AudiPery, voicedMask *&TGroup) {
    featurePara x = tandemPara();

    AudiPery = new gammaToneFilterBank(x.gtP);
    TGroup = new voicedMask(x);

    // Explicitly set sample rate (may not be set in constructor)
    TGroup->fs = 20000;
    TGroup->numberChannel = 64;

    // corrLgm is allocated by computeFeatureTandem() (not by the vendored
    // feature::newFeature()), so it is owned here and released by
    // releaseTandemFeatures().
    TGroup->corrLgm = NULL;
}

// Free the correlogram allocated by computeFeatureTandem(). Must be called
// before deleting TGroup (the vendored destructor does not own it).
void releaseTandemFeatures(voicedMask *TGroup) {
    if (TGroup && TGroup->corrLgm) {
        delete [] TGroup->corrLgm;
        TGroup->corrLgm = NULL;
    }
}

// ---------------------------------------------------------------------------
// Parallel re-implementation of the vendored feature::computeFeature()
// (src/tandem/tandem_64/feature.cpp). The arithmetic is a verbatim copy of
// feature::fftACF() / feature::computeCross(); only the scratch buffers are
// made per-thread (the vendored version uses the single shared member
// `data`, which serialises the 64 channels). Channels are independent:
// each writes only corrLgm[*].acf/acfEv/zc/zcEv[chan]. computeCross() for
// frame n reads only corrLgm[n]. The output is therefore bit-identical to
// the serial vendored routine for any thread count.
// ---------------------------------------------------------------------------

struct tandemAcfWork {
    std::vector<double> sumS;
    std::vector<double> inputXR, inputXI, inputXWR, inputXWI;
};

// Verbatim copy of feature::fftACF() operating on thread-local buffers.
static void tandemFftACF(const double* const response, int frame, tandemAcfWork &data,
                         int acfLen, int acfOrder, int window, int max_delay, int sigLength)
{
    int step, delay;

    for(step=0; step<acfLen; step++)
        data.inputXR[step] = data.inputXI[step] = data.inputXWR[step] = data.inputXWI[step] = 0;

    for(step=0; step<(window+max_delay); step++)
    {
        int tim = (frame+2) * window/2 - (step+1);
        if ( (tim>= 0)  && (tim<sigLength) )
        {
            data.inputXR[step] = response[tim];

            if(step<window) data.inputXWR[step] = data.inputXR[step];
        }
    }

    data.sumS[0]=0;
    for(step=0; step<window; step++)
    {
        data.sumS[0] += data.inputXR[step]*data.inputXR[step]/double(window);
    }

    for(step=1; step<max_delay; step++)
    {
        double f1 = data.inputXR[step+window-1] - data.inputXR[step-1];
        double f2 = data.inputXR[step+window-1] + data.inputXR[step-1];

        data.sumS[step] = data.sumS[step-1] + f1*f2/double(window);

        if (data.sumS[step]<0) data.sumS[step]=0;
    }

    fft(data.inputXR.data(), data.inputXI.data(), acfOrder, 1);
    fft(data.inputXWR.data(), data.inputXWI.data(), acfOrder, 1);

    for(step=0; step<acfLen; step++)
    {
        double f1 = data.inputXR[step]*data.inputXWR[step] + data.inputXI[step]*data.inputXWI[step];
        double f2 = -data.inputXR[step]*data.inputXWI[step] + data.inputXI[step]*data.inputXWR[step];

        data.inputXR[step]=f1;
        data.inputXI[step]=f2;
    }

    fft(data.inputXR.data(), data.inputXI.data(), acfOrder, -1);

    for(delay=0; delay<max_delay; delay++)
    {
        data.inputXR[delay] /= double(acfLen*window);
        data.inputXR[delay] /= double(sqrt(data.sumS[0]*data.sumS[delay]+1e-5*data.sumS[0])) / data.sumS[0];
    }
}

// Verbatim copy of feature::computeCross().
static void tandemComputeCross(acfFeature *corrLgm, int n, int numberChannel, int max_delay)
{
    for(int chan=0; chan<(numberChannel-1); chan++)
    {
        double m1=0, m2=0, sum=0, sum1=0, sum2=0;

        for(int delay=0; delay<max_delay; delay++)
        {
            m1 += corrLgm[n].acf[chan][delay];
            m2 += corrLgm[n].acf[chan+1][delay];

            sum += corrLgm[n].acf[chan][delay] * corrLgm[n].acf[chan+1][delay];
            sum1 += corrLgm[n].acf[chan][delay] * corrLgm[n].acf[chan][delay];
            sum2 += corrLgm[n].acf[chan+1][delay] * corrLgm[n].acf[chan+1][delay];
        }

        m1 /= double(max_delay);
        m2 /= double(max_delay);
        corrLgm[n].cross[chan] = (sum-m1*m2) /sqrt( (sum1-m1*m1) * (sum2-m2*m2) + 1e-100 );

        m1=0, m2=0, sum=0, sum1=0, sum2=0;

        for(int delay=0; delay<max_delay; delay++)
        {
            m1 += corrLgm[n].acfEv[chan][delay];
            m2 += corrLgm[n].acfEv[chan+1][delay];

            sum += corrLgm[n].acfEv[chan][delay] * corrLgm[n].acfEv[chan+1][delay];
            sum1 += corrLgm[n].acfEv[chan][delay] * corrLgm[n].acfEv[chan][delay];
            sum2 += corrLgm[n].acfEv[chan+1][delay] * corrLgm[n].acfEv[chan+1][delay];
        }

        m1 /= double(max_delay);
        m2 /= double(max_delay);
        corrLgm[n].crossEv[chan] = (sum-m1*m2) /sqrt( (sum1-m1*m1) * (sum2-m2*m2) + 1e-100 );
    }

    corrLgm[n].cross[numberChannel-1]=0;
    corrLgm[n].crossEv[numberChannel-1]=0;
}

static inline int tandemThreadId(void) {
#ifdef _OPENMP
    return omp_get_thread_num();
#else
    return 0;
#endif
}

// Replacement for TGroup->computeFeature(AudiPery, sLen, 0).
static void computeFeatureTandem(const gammaToneFilterBank* const AudiPery,
                                 voicedMask *TGroup, int sLen, int nThreads)
{
    const featurePara x = tandemPara();
    const int fs = TGroup->fs;
    const int numberChannel = TGroup->numberChannel;

    // feature::feature() / feature::newFeature() derived quantities
    const int max_delay = fs*3/200;
    TGroup->sigLength = sLen;
    TGroup->numFrame = sLen*100/(fs);
    const int numFrame = TGroup->numFrame;
    const int window = fs/50;
    int acfLen = max_delay + window;
    int acfOrder = ceil(log(double(acfLen)) / log(2.0));
    acfLen = pow(2.0, acfOrder);

    releaseTandemFeatures(TGroup);
    TGroup->corrLgm = new acfFeature[numFrame];
    acfFeature *corrLgm = TGroup->corrLgm;

    // Per-thread scratch, all allocated outside the parallel region.
    std::vector<std::unique_ptr<filter> > bandPass(nThreads);
    std::vector<std::vector<double> > temp(nThreads), envelope(nThreads);
    std::vector<tandemAcfWork> work(nThreads);
    for (int t = 0; t < nThreads; t++) {
        bandPass[t].reset(new filter( (x.bP1 + x.bP2) / double(fs), (x.bP2 - x.bP1) / double(fs),
                                      double(x.bPTs*2) / double(fs), fs, 0));
        temp[t].assign(sLen, 0.0);
        envelope[t].assign(sLen, 0.0);
        work[t].sumS.assign(max_delay, 0.0);
        work[t].inputXR.assign(acfLen, 0.0);
        work[t].inputXI.assign(acfLen, 0.0);
        work[t].inputXWR.assign(acfLen, 0.0);
        work[t].inputXWI.assign(acfLen, 0.0);
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(nThreads) if (nThreads > 1)
#endif
    for (int chan = 0; chan < numberChannel; chan++)
    {
        const int t = tandemThreadId();
        double *tmp = temp[t].data();
        double *env = envelope[t].data();
        tandemAcfWork &data = work[t];
        const double *resp = AudiPery->response[chan];

        for (int n = 0; n < sLen; n++)
            tmp[n] = (resp[n]>0) ? resp[n]:0;

        bandPass[t]->filtering(tmp, env, sLen);

        for (int n = 0; n < numFrame; n++)
        {
            tandemFftACF(env, n, data, acfLen, acfOrder, window, max_delay, sLen);
            for (int m = 0; m < max_delay; m++)
                corrLgm[n].acfEv[chan][m] = data.inputXR[m];
            corrLgm[n].zcEv[chan] = zeroCross(data.inputXR.data(), max_delay);

            tandemFftACF(resp, n, data, acfLen, acfOrder, window, max_delay, sLen);
            for (int m = 0; m < max_delay; m++)
                corrLgm[n].acf[chan][m] = data.inputXR[m];
            corrLgm[n].zc[chan] = zeroCross(data.inputXR.data(), max_delay);
        }
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads) if (nThreads > 1)
#endif
    for (int n = 0; n < numFrame; n++)
        tandemComputeCross(corrLgm, n, numberChannel, max_delay);
}

// Memory-based pitch estimation (replacement for voicedMaskEst)
void voicedMaskEstMemory(double *audio, int sigLength,
                        gammaToneFilterBank *AudiPery,
                        voicedMask *TGroup,
                        const char *netPath,
                        int nThreads) {
    // Note: Neural networks are loaded automatically in pitchMask constructor
    // from "net/MLP*.64.dat" - R wrapper must ensure this directory exists
    if (nThreads < 1) nThreads = 1;
#ifndef _OPENMP
    nThreads = 1;
#endif

    // Scale audio to proper loudness
    double scale;
    std::unique_ptr<double[]> Input(processAudioBufferTandem(audio, sigLength, &scale));

    // Set number of frames (100 Hz frame rate)
    TGroup->numFrame = sigLength / (TGroup->fs / 100);

    // Allocate pitch parameters
    TGroup->newPitchPara();

    // Process through gammatone filterbank. Each channel owns its filter
    // state (gf[chan]), its memIndi[chan] flag and its output buffer
    // (response[chan]), so channels can run concurrently. An allocation
    // failure is recorded and rethrown outside the parallel region.
    AudiPery->sigLength = sigLength;
    const int nChan = TGroup->numberChannel;
    {
        int failed = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(nThreads) if (nThreads > 1)
#endif
        for (int chan = 0; chan < nChan; chan++) {
            try {
                AudiPery->filtering(Input.get(), sigLength, chan);
            } catch (...) {
#ifdef _OPENMP
#pragma omp atomic write
#endif
                failed = 1;
            }
        }
        if (failed) throw std::bad_alloc();
    }

    // Compute features (parallel equivalent of TGroup->computeFeature)
    computeFeatureTandem(AudiPery, TGroup, sigLength, nThreads);

    // Pre-compute the single-unit probabilities that
    // voicedMask::initPitchEst() requests for every frame >= 10 and every
    // delay in [min_pitch, max_pitch] (sf/400 .. sf/70). The vendored
    // singleUnitProb() is a pure function of corrLgm and the (read-only)
    // networks and writes only Prob[frame].{sProb,hNum,hNumEv}[*][delay], so
    // filling it frame-parallel in advance yields the same values that the
    // lazy (sentinel -2) evaluation inside initPitchEst() would.
    {
        const int numFrame = TGroup->numFrame;
        const int min_pitch = TGroup->fs/400;
        const int max_pitch = TGroup->fs/70;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 4) num_threads(nThreads) if (nThreads > 1)
#endif
        for (int frame = 10; frame < numFrame; frame++)
            for (int delay = min_pitch; delay <= max_pitch; delay++)
                TGroup->pitchProb(frame, delay);
    }

    TGroup->dtmPitchMask();
}

// Structure to hold TANDEM results
struct TandemResults {
    std::vector<double> pitch;       // F0 values (Hz)
    std::vector<double> voicing;     // Voicing probability (0-1)
    std::vector<double> pitch_confidence;  // Pitch confidence scores
    int n_frames;
    int n_contours;
    int sample_rate;
};

// Extract pitch contours from TANDEM voicedMask object
TandemResults extractPitchContoursTandem(voicedMask *TGroup) {
    TandemResults results;
    results.n_frames = TGroup->numFrame;
    results.n_contours = TGroup->numContour;
    results.sample_rate = TGroup->fs;

    // Initialize vectors with NaN/zeros
    results.pitch.resize(TGroup->numFrame, R_NaReal);
    results.voicing.resize(TGroup->numFrame, 0.0);
    results.pitch_confidence.resize(TGroup->numFrame, 0.0);

    // Extract from TANDEM pitch contours
    // TANDEM can track multiple pitch contours simultaneously
    // We'll take the most confident one for each frame
    for (int c = 0; c < TGroup->numContour; c++) {
        pitchContour *pc = &(TGroup->Pitch[c]);

        for (int f = pc->sFrame; f <= pc->eFrame && f < TGroup->numFrame; f++) {
            int idx = f - pc->sFrame;

            // Check if this index is valid
            if (idx >= 0 && pc->value != NULL) {
                int delay = pc->value[idx];

                // Convert delay samples to Hz
                if (delay > 0) {
                    double freq = (double)TGroup->fs / (double)delay;

                    // Sanity check frequency range (20-2000 Hz)
                    if (freq >= 20.0 && freq <= 2000.0) {
                        // Get probability/confidence from mask structure
                        double prob = 0.5;  // Default confidence
                        if (pc->mProb != NULL && idx < (pc->eFrame - pc->sFrame + 1)) {
                            // mProb is array of mask structures
                            // Try to get voicing probability
                            for (int ch = 0; ch < TGroup->numberChannel; ch++) {
                                if (pc->mProb[idx].value[ch] > prob) {
                                    prob = pc->mProb[idx].value[ch];
                                }
                            }
                        }

                        // Take pitch with highest confidence for this frame
                        if (prob > results.pitch_confidence[f]) {
                            results.pitch[f] = freq;
                            results.voicing[f] = prob;
                            results.pitch_confidence[f] = prob;
                        }
                    }
                }
            }
        }
    }

    return results;
}

// Get neural network paths from R package directory
void getTandemNetPaths(const char *pkg_dir,
                      char *net1, char *net2, char *net3) {
    snprintf(net1, 512, "%s/tandem_net/MLP1.64.dat", pkg_dir);
    snprintf(net2, 512, "%s/tandem_net/MLP2.64.dat", pkg_dir);
    snprintf(net3, 512, "%s/tandem_net/MLP3.64.dat", pkg_dir);
}
