/*
 * tandem_wrapper.cpp
 * 
 * Rcpp wrapper for TANDEM pitch tracking and voiced speech segregation
 * 
 * Based on: Hu, G., & Wang, D. L. (2010). "A tandem algorithm for pitch 
 * estimation and voiced speech segregation." IEEE Trans. Audio, Speech, 
 * Lang. Process., 18(8), 2067-2079.
 *
 * Copyright (C) 2025 superassp contributors
 * Original TANDEM code: G. Hu & D. L. Wang, Ohio State University
 */

#include <Rcpp.h>
#include <vector>
#include <string>
#include <cmath>

// Forward declarations from tandem_memory.cpp
struct TandemResults {
    std::vector<double> pitch;
    std::vector<double> voicing;
    std::vector<double> pitch_confidence;
    int n_frames;
    int n_contours;
    int sample_rate;
};

class gammaToneFilterBank;
class voicedMask;

// C++ linkage (not extern "C") - these are C++ functions
void initVoicedMaskTandem(gammaToneFilterBank *&AudiPery, voicedMask *&TGroup);
void voicedMaskEstMemory(double *audio, int sigLength, 
                        gammaToneFilterBank *AudiPery, 
                        voicedMask *TGroup,
                        const char *netPath);
TandemResults extractPitchContoursTandem(voicedMask *TGroup);
void getTandemNetPaths(const char *pkg_dir, char *net1, char *net2, char *net3);

//' TANDEM Pitch Tracking (C++ Interface)
//'
//' Low-level C++ function for TANDEM algorithm. Users should call trk_tandem() instead.
//'
//' @param audio_signal Numeric vector, audio samples (mono)
//' @param sample_rate Integer, sample rate (must be 20000 Hz for TANDEM)
//' @param min_pitch Numeric, minimum F0 in Hz
//' @param max_pitch Numeric, maximum F0 in Hz
//' @param net_path Character, path to neural network weight files
//' @return List with pitch, voicing_prob, and pitch_confidence
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List tandem_pitch_cpp(
    Rcpp::NumericVector audio_signal,
    int sample_rate = 20000,
    double min_pitch = 50.0,
    double max_pitch = 500.0,
    std::string net_path = ""
) {
    // Validate inputs
    if (sample_rate != 20000) {
        Rcpp::warning("TANDEM requires 20 kHz sample rate. Accuracy may be reduced.");
    }
    
    int n_samples = audio_signal.size();
    
    if (n_samples < 1000) {
        Rcpp::stop("Audio signal too short (minimum 1000 samples required)");
    }
    
    // Initialize TANDEM components
    gammaToneFilterBank *audioPery = NULL;
    voicedMask *tGroup = NULL;
    
    try {
        // Initialize filterbank and voicedMask
        initVoicedMaskTandem(audioPery, tGroup);
        
        // Get network paths (currently unused as networks are embedded in TANDEM)
        // In future, could load from inst/tandem_net/
        const char *net_path_cstr = net_path.c_str();
        
        // Process audio through TANDEM
        double *samples = REAL(audio_signal);
        voicedMaskEstMemory(samples, n_samples, audioPery, tGroup, net_path_cstr);
        
        // Extract pitch contours
        TandemResults results = extractPitchContoursTandem(tGroup);
        
        // Convert to R vectors
        Rcpp::NumericVector pitch(results.pitch.begin(), results.pitch.end());
        Rcpp::NumericVector voicing(results.voicing.begin(), results.voicing.end());
        Rcpp::NumericVector confidence(results.pitch_confidence.begin(), 
                                      results.pitch_confidence.end());
        
        // Cleanup
        delete audioPery;
        delete tGroup;
        
        return Rcpp::List::create(
            Rcpp::Named("pitch") = pitch,
            Rcpp::Named("voicing_prob") = voicing,
            Rcpp::Named("pitch_confidence") = confidence,
            Rcpp::Named("sample_rate") = results.sample_rate,
            Rcpp::Named("n_frames") = results.n_frames,
            Rcpp::Named("n_contours") = results.n_contours,
            Rcpp::Named("status") = "full_tandem"
        );
        
    } catch(std::exception &e) {
        // Cleanup on error
        if (audioPery) delete audioPery;
        if (tGroup) delete tGroup;
        Rcpp::stop("TANDEM processing error: %s", e.what());
    }
}

