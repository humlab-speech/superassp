# Summary table of superassp DSP function output

The summary table produced by this function lists the DSP function
names, default file extension, and a summary of produced SSFF tracks in
the output file, or the number of fields in slice producing functions.

## Usage

``` r
superassp_summary()
```

## Value

A data.frame with function names as row labels, and with "extension" and
"tracks" columns. The output is ordered by file extension in
alphabetical order by default to make it iseasier to make sure DSP data
are not overwritten when multiple functions are applied to the same
recordings.

## Details

If the number of tracks or fields produced by the function is very
large, then the output is truncated to a summary of the number of
tracks.

## Examples

``` r
superassp:::superassp_summary()
#>                           extension                          tracks
#> trk_acf                         acf                             ACF
#> trk_d4c                          ap                    aperiodicity
#> trk_arf                         arf            RMS[dB],gain[dB],ARF
#> lst_avqi                       avqi                      (8 tracks)
#> trk_cepstrum                    cep                           C[dB]
#> trk_formant_cgdzp               cgf                      (5 tracks)
#> lst_ComParE_2016                cmp                   (6373 tracks)
#> trk_cpps                        cps                             cpp
#> trk_covarep_creak               crk              creak_pp,creak_bin
#> trk_pitch_crepe                 crp                  f0,periodicity
#> trk_creak_vat                   crv              creak_pp,creak_bin
#> trk_css_spectrum                css                         CSS[dB]
#> trk_covarep_vad_drugman         cvd                      (4 tracks)
#> lst_covarep_vq                  cvq                                
#> trk_formant_deepformants        dff                              fm
#> trk_dft_spectrum                dft                         DFT[dB]
#> trk_afdiff                      dif                                
#> lst_dsi                         dsi                      (6 tracks)
#> lst_dysprosody                  dyp                                
#> lst_eGeMAPS                     egm                     (88 tracks)
#> lst_emobase                     emb                    (988 tracks)
#> trk_covarep_env_te              ete                   env_te,env_cc
#> trk_pitch_dio                    f0                              f0
#> trk_pitch_harvest                f0                              f0
#> trk_pitch_rapt                   f0                              f0
#> trk_pitch_reaper                 f0                              f0
#> trk_pitch_swipe                  f0                              f0
#> trk_pitch_vat                   f0v                  f0,vad,srh_val
#> trk_affilter                    flt                                
#> trk_formant_forest              fms                     F[Hz],B[Hz]
#> trk_formant_formantnet          fnf                           fm,bw
#> trk_ksvfo                        fo                          fo[Hz]
#> trk_pitch_ksv                    fo                          fo[Hz]
#> trk_gci_vat                    gciv             gci_sample,residual
#> lst_GeMAPS                      gem                     (62 tracks)
#> trk_gfmiaif                     gfm                     (55 tracks)
#> trk_covarep_iaif                glf glottal_flow,glottal_derivative
#> trk_iaif_vat                    glv glottal_flow,glottal_derivative
#> trk_hmpd                        hpd                      ae,pdm,pdd
#> trk_intensity                   int                       intensity
#> trk_lar                         lar            RMS[dB],gain[dB],LAR
#> trk_lpc                         lpc           RMS[dB],gain[dB],LPCi
#> trk_lps_spectrum                lps                         LPS[dB]
#> trk_mdq_vat                     mdq                             mdq
#> trk_mfcc                       mfcc                     (14 tracks)
#> trk_pitch_ac                    pac                              F0
#> trk_pitch_cc                    pcc                              F0
#> trk_pitch_pda                   pda                              F0
#> trk_formant_burg                pfm                     (15 tracks)
#> trk_pitch_mhs                   pit                       pitch[Hz]
#> trk_pitchmark_estk               pm                      pitchmarks
#> trk_praatsauce                  psa                     (37 tracks)
#> trk_pitch_shs                   psh                              F0
#> trk_peakslope                   psl                       peakslope
#> trk_pitch_spinet                psp                              F0
#> trk_peakslope_vat               psv                      peak_slope
#> lst_voice_report                pvr                     (30 tracks)
#> lst_voice_tremor                pvt                     (18 tracks)
#> trk_pitch_pyin                  pyp                         F0,prob
#> trk_rfc                         rfc            RMS[dB],gain[dB],RFC
#> trk_rms                         rms                         RMS[dB]
#> trk_pitchmark_reaper            rpm                              pm
#> trk_pitch_swiftf0               sf0                   f0,confidence
#> trk_formant_snack          snackfmt                           fm,bw
#> trk_pitch_snack          snackpitch                      (4 tracks)
#> trk_cheap_trick                  sp                              sp
#> trk_spectral_moments            spm                      (4 tracks)
#> trk_pitch_srh                   srh                          f0,vad
#> trk_tandem                      tnd              pitch,voicing_prob
#> trk_formant_tvwlp               tvf                           fm,bw
#> trk_covarep_vq_gci              vqg                      (5 tracks)
#> lst_vq_vat                      vqv                                
#> trk_vuv                         vuv                         voicing
#> lst_voxit                       vxt                                
#> trk_pitch_yin                   yip                         F0,prob
#> trk_zcr                         zcr                         ZCR[Hz]
```
