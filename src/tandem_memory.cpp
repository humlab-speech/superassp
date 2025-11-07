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

// Initialize TANDEM voicedMask and gammaTone filterbank
void initVoicedMaskTandem(gammaToneFilterBank *&AudiPery, voicedMask *&TGroup) {
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
    
    AudiPery = new gammaToneFilterBank(x.gtP);
    TGroup = new voicedMask(x);
}

// Memory-based pitch estimation (replacement for voicedMaskEst)
void voicedMaskEstMemory(double *audio, int sigLength, 
                        gammaToneFilterBank *AudiPery, 
                        voicedMask *TGroup,
                        const char *netPath) {
    // Scale audio to proper loudness
    double scale;
    double *Input = processAudioBufferTandem(audio, sigLength, &scale);
    
    // Set number of frames (100 Hz frame rate)
    TGroup->numFrame = sigLength / (TGroup->fs / 100);
    
    // Allocate pitch parameters
    TGroup->newPitchPara();
    
    // Process through gammatone filterbank
    AudiPery->sigLength = sigLength;
    for(int chan = 0; chan < TGroup->numberChannel; chan++) {
        AudiPery->filtering(Input, sigLength, chan);
    }
    
    // Compute features and estimate pitch/mask
    TGroup->computeFeature(AudiPery, sigLength, 1);
    TGroup->dtmPitchMask();
    
    delete [] Input;
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
            
            if (pc->indicate[idx] > 0 && idx >= 0) {
                // Convert delay samples to Hz
                int delay = pc->value[idx];
                if (delay > 0) {
                    double freq = (double)TGroup->fs / (double)delay;
                    
                    // Get probability/confidence
                    double prob = 0.0;
                    if (pc->mProb && pc->mProb[idx].value[0] > 0) {
                        prob = pc->mProb[idx].value[0];
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
    
    return results;
}

// Get neural network paths from R package directory
void getTandemNetPaths(const char *pkg_dir, 
                      char *net1, char *net2, char *net3) {
    snprintf(net1, 512, "%s/tandem_net/MLP1.64.dat", pkg_dir);
    snprintf(net2, 512, "%s/tandem_net/MLP2.64.dat", pkg_dir);
    snprintf(net3, 512, "%s/tandem_net/MLP3.64.dat", pkg_dir);
}
