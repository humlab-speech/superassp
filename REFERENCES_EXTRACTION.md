# R Package References Extraction

**Package**: superassp
**Generated**: 2026-03-01
**Total Files**: 58 R files with @references sections

---

## Summary

| Category | Count |
|----------|-------|
| Files with references | 58 |
| Journal articles | 28 |
| Standards/ISO | 8 |
| Technical documents | 12 |
| URLs/Online resources | 8 |
| Books/Theses | 4 |
| Software/Tools | 6 |
| Using \insertAllCited{} | 9 |
| Using \insertRef{} | 6 |
| Empty/TBD | 5 |

---

## References by File

### 1. R/covarep_vq.R

**Function**: `lst_covarep_vq`
**Type**: Journal articles

```
Alku, P. (1992). "Glottal wave analysis with pitch synchronous iterative
adaptive inverse filtering". Speech Communication, 11(2-3), 109-118.

Kane, J., & Gobl, C. (2013). "Evaluation of glottal closure instant
detection in a range of voice qualities". Speech Communication, 55(2), 295-314.
```

**BibTeX Entry**:
```bibtex
@article{Alku1992,
  author = {Alku, P.},
  title = {Glottal wave analysis with pitch synchronous iterative adaptive inverse filtering},
  journal = {Speech Communication},
  year = {1992},
  volume = {11},
  number = {2-3},
  pages = {109--118}
}

@article{Kane2013,
  author = {Kane, J. and Gobl, C.},
  title = {Evaluation of glottal closure instant detection in a range of voice qualities},
  journal = {Speech Communication},
  year = {2013},
  volume = {55},
  number = {2},
  pages = {295--314}
}
```

---

### 2. R/install_phonet.R

**Function**: `install_phonet`
**Type**: Citation reference
**Reference**: `\insertAllCited{}`

**Status**: Uses dynamic bibliography (BibTeX entries in BibTeX file)

---

### 3. R/install_pladdrr.R

**Function**: `install_pladdrr`
**Type**: URLs/Software references

```
pladdrr GitHub: https://github.com/tjmahr/pladdrr
Praat: http://www.praat.org
```

**BibTeX Entry**:
```bibtex
@software{pladdrr,
  title = {pladdrr: R interface to Praat},
  url = {https://github.com/tjmahr/pladdrr}
}

@software{Praat,
  title = {Praat: Doing phonetics by computer},
  url = {http://www.praat.org}
}
```

---

### 4. R/iso226_phon.R

**Functions**:
- `iso226_phon` (data object)
- `ucnv_db_and_hz_to_phon`

**Type**: International Standard

```
ISO 226:2023, Acoustics - Normal equal-loudness-level contours.
International Organization for Standardization, Geneva, Switzerland.

ISO 226:2023, Acoustics - Normal equal-loudness-level contours.
Formula (2), Section 4.2, page 3.

ISO 226:2023, Acoustics - Normal equal-loudness-level contours.
Formula (1), Section 4.1, page 2.
```

**BibTeX Entry**:
```bibtex
@standard{ISO226_2023,
  title = {Acoustics - Normal equal-loudness-level contours},
  organization = {International Organization for Standardization},
  address = {Geneva, Switzerland},
  year = {2023},
  number = {ISO 226:2023}
}
```

---

### 5. R/iso532_sone.R

**Functions**:
- `iso532_sone` (data object)
- `.moore_glasberg_table`

**Type**: International Standards + Journal article

```
ISO 532-1:2017 - Acoustics - Methods for calculating loudness -
Part 1: Zwicker method

ISO 532-2:2017 - Acoustics - Methods for calculating loudness -
Part 2: Moore-Glasberg method

Stevens, S. S. (1936). "A scale for the measurement of a psychological
magnitude: loudness". Psychological Review. 43 (5): 405–416.

MATLAB Audio Toolbox documentation for sone2phon (ISO 532-2 reference table)
```

**BibTeX Entry**:
```bibtex
@standard{ISO532_1_2017,
  title = {Acoustics - Methods for calculating loudness},
  part = {1: Zwicker method},
  organization = {International Organization for Standardization},
  year = {2017},
  number = {ISO 532-1:2017}
}

@standard{ISO532_2_2017,
  title = {Acoustics - Methods for calculating loudness},
  part = {2: Moore-Glasberg method},
  organization = {International Organization for Standardization},
  year = {2017},
  number = {ISO 532-2:2017}
}

@article{Stevens1936,
  author = {Stevens, S. S.},
  title = {A scale for the measurement of a psychological magnitude: loudness},
  journal = {Psychological Review},
  year = {1936},
  volume = {43},
  number = {5},
  pages = {405--416}
}
```

---

### 6. R/list_cpp_opensmile_gemaps.R

**Function**: `lst_GeMAPS`
**Type**: Citation reference
**Reference**: `\insertAllCited{}`

**Status**: Uses dynamic bibliography

---

### 7. R/list_pladdrr_dysprosody.R

**Function**: `lst_dysprosody`
**Type**: Journal article

```
Villarubia et al. (2025). Dysprosody: Comprehensive prosodic feature
extraction for speech analysis. Frontiers in Human Neuroscience, 19, 1566274.
doi:10.3389/fnhum.2025.1566274
```

**BibTeX Entry**:
```bibtex
@article{Villarubia2025,
  author = {Villarubia, et al.},
  title = {Dysprosody: Comprehensive prosodic feature extraction for speech analysis},
  journal = {Frontiers in Human Neuroscience},
  year = {2025},
  volume = {19},
  pages = {1566274},
  doi = {10.3389/fnhum.2025.1566274}
}
```

---

### 8. R/list_pladdrr_pharyngeal.R

**Function**: `lst_pharyngeal`
**Type**: Conference proceeding

```
Iseli, M., & Alwan, A. (2004). An improved correction formula for the
estimation of harmonic magnitudes and its application to open quotient
estimation. IEEE ICASSP, 1, 669-672.
```

**BibTeX Entry**:
```bibtex
@inproceedings{Iseli2004,
  author = {Iseli, M. and Alwan, A.},
  title = {An improved correction formula for the estimation of harmonic magnitudes
           and its application to open quotient estimation},
  booktitle = {IEEE ICASSP},
  year = {2004},
  volume = {1},
  pages = {669--672}
}
```

---

### 9. R/list_pladdrr_vq.R

**Function**: `lst_vq`
**Type**: Software + Software documentation

```
VQ_measurements_V2.praat (Version 2: 20 October 2024)
Boersma, P., & Weenink, D. (2023). Praat: doing phonetics by computer.
```

**BibTeX Entry**:
```bibtex
@software{VQ_measurements_V2,
  title = {VQ\_measurements\_V2.praat},
  version = {2},
  date = {2024-10-20}
}

@software{Boersma2023,
  author = {Boersma, P. and Weenink, D.},
  title = {Praat: doing phonetics by computer},
  year = {2023}
}
```

---

### 10. R/list_python_opensmile_ComParE_2016.R

**Function**: `lst_ComParE_2016`
**Type**: Citation reference
**Reference**: `\insertAllCited{}`

**Status**: Uses dynamic bibliography

---

### 11. R/list_python_opensmile_eGeMAPS.R

**Function**: `lst_eGeMAPS`
**Type**: Citation reference
**Reference**: `\insertAllCited{}` (Empty line)

**Status**: Uses dynamic bibliography

---

### 12. R/list_python_opensmile_emobase.R

**Function**: `lst_emobase`
**Type**: Citation reference
**Reference**: `\insertAllCited{}` (Empty line)

**Status**: Uses dynamic bibliography

---

### 13. R/list_python_pm_pavqi.R

**Function**: `lst_avqip`
**Type**: Journal articles

```
Maryn, Y., Corthals, P., Van Cauwenberge, P., Roy, N., & De Bodt, M. (2010).
Toward improved ecological validity in the acoustic measurement of overall
voice quality: Combining continuous speech and sustained vowels.
Journal of Voice, 24(5), 540-555.

Barsties, B., & Maryn, Y. (2015). The improvement of internal consistency of
the Acoustic Voice Quality Index. American Journal of Otolaryngology, 36(5), 647-656.
```

**BibTeX Entry**:
```bibtex
@article{Maryn2010,
  author = {Maryn, Y. and Corthals, P. and Van Cauwenberge, P. and Roy, N. and De Bodt, M.},
  title = {Toward improved ecological validity in the acoustic measurement of overall
           voice quality: Combining continuous speech and sustained vowels},
  journal = {Journal of Voice},
  year = {2010},
  volume = {24},
  number = {5},
  pages = {540--555}
}

@article{Barsties2015,
  author = {Barsties, B. and Maryn, Y.},
  title = {The improvement of internal consistency of the Acoustic Voice Quality Index},
  journal = {American Journal of Otolaryngology},
  year = {2015},
  volume = {36},
  number = {5},
  pages = {647--656}
}
```

---

### 14. R/list_python_pm_pdsi.R

**Function**: `lst_dsip`
**Type**: Journal article

```
Wuyts, F. L., De Bodt, M. S., Molenberghs, G., Remacle, M., Heylen, L.,
Millet, B., ... & Heyning, P. H. V. D. (2000). The dysphonia severity index:
an objective measure of vocal quality based on a multiparameter approach.
Journal of Speech, Language, and Hearing Research, 43(3), 796-809.
```

**BibTeX Entry**:
```bibtex
@article{Wuyts2000,
  author = {Wuyts, F. L. and De Bodt, M. S. and Molenberghs, G. and Remacle, M. and
            Heylen, L. and Millet, B. and others and Heyning, P. H. V. D.},
  title = {The dysphonia severity index: an objective measure of vocal quality based
           on a multiparameter approach},
  journal = {Journal of Speech, Language, and Hearing Research},
  year = {2000},
  volume = {43},
  number = {3},
  pages = {796--809}
}
```

---

### 15. R/list_python_pm_pvoice_tremor.R

**Function**: `lst_voice_tremorp`
**Type**: Conference proceeding

```
Brückl, M. (2012). Vocal tremor measurement based on autocorrelation of contours.
Proceedings of Interspeech 2012, 2027-2031.
```

**BibTeX Entry**:
```bibtex
@inproceedings{Bruckl2012,
  author = {Br\"{u}ckl, M.},
  title = {Vocal tremor measurement based on autocorrelation of contours},
  booktitle = {Proceedings of Interspeech 2012},
  year = {2012},
  pages = {2027--2031}
}
```

---

### 16. R/list_vat.R

**Function**: `lst_vat`
**Type**: Journal article + Thesis

```
Tsanas, A., Little, M., McSharry, P., & Ramig, L. (2011).
Nonlinear speech analysis algorithms mapped to a standard metric achieve
clinically useful quantification of average Parkinson's disease symptom severity.
Journal of the Royal Society Interface, 8(59), 842-855.

Tsanas, A. (2012). Accurate telemonitoring of Parkinson's disease symptom severity
using nonlinear speech signal processing and statistical machine learning.
D.Phil. thesis, University of Oxford.
```

**BibTeX Entry**:
```bibtex
@article{Tsanas2011,
  author = {Tsanas, A. and Little, M. and McSharry, P. and Ramig, L.},
  title = {Nonlinear speech analysis algorithms mapped to a standard metric achieve
           clinically useful quantification of average Parkinson's disease symptom severity},
  journal = {Journal of the Royal Society Interface},
  year = {2011},
  volume = {8},
  number = {59},
  pages = {842--855}
}

@thesis{Tsanas2012,
  author = {Tsanas, A.},
  title = {Accurate telemonitoring of Parkinson's disease symptom severity using
           nonlinear speech signal processing and statistical machine learning},
  school = {University of Oxford},
  year = {2012},
  type = {D.Phil. thesis}
}
```

---

### 17. R/list_voxit.R

**Function**: `lst_voxit`
**Type**: Conference proceeding + Online resource

```
Voxit toolbox: Voice and articulation complexity measures

SAcC pitch tracker: Ellis, D.P.W., & Weiss, R.J. (2010).
Pitch and voicing estimation via a harmonic model.
```

**BibTeX Entry**:
```bibtex
@software{Voxit,
  title = {Voxit toolbox: Voice and articulation complexity measures}
}

@inproceedings{Ellis2010,
  author = {Ellis, D. P. W. and Weiss, R. J.},
  title = {Pitch and voicing estimation via a harmonic model},
  year = {2010}
}
```

---

### 18. R/lst_voice_sauce.R

**Function**: `lst_voice_sauce`
**Type**: Journal article + Conference proceeding

```
Shue, Y. L., Keating, P., Vicenik, C., & Yu, K. (2011). "VoiceSauce:
A program for voice analysis". In Proceedings of ICPhS XVII (pp. 1846-1849).

Maryn, Y., Corthals, P., Van Cauwenberge, P., Roy, N., & De Bodt, M. (2010).
"Toward improved ecological validity in the acoustic measurement of overall
voice quality: combining microphone and neck-surface accelerometer data".
Journal of Voice, 24(5), 540-554.
```

**BibTeX Entry**:
```bibtex
@inproceedings{Shue2011,
  author = {Shue, Y. L. and Keating, P. and Vicenik, C. and Yu, K.},
  title = {VoiceSauce: A program for voice analysis},
  booktitle = {Proceedings of ICPhS XVII},
  year = {2011},
  pages = {1846--1849}
}

@article{Maryn2010,
  author = {Maryn, Y. and Corthals, P. and Van Cauwenberge, P. and Roy, N. and De Bodt, M.},
  title = {Toward improved ecological validity in the acoustic measurement of overall
           voice quality: combining microphone and neck-surface accelerometer data},
  journal = {Journal of Voice},
  year = {2010},
  volume = {24},
  number = {5},
  pages = {540--554}
}
```

---

### 19. R/prep_recode.R

**Function**: Unknown
**Type**: Citation references
**Reference**: `\insertRef{av2024}{superassp}`, `\insertRef{ffmpeg2024}{superassp}`

**Status**: Uses dynamic bibliography

---

### 20. R/psychoacoustic_scales.R

**Functions**:
- `psychoacoustic_scales` (data object)
- `onLoad_psychoacoustic_units`
- `f` (internal function)
- `ucnv_mel_to_hz`

**Type**: Journal articles + Book

```
Traunmüller, H. (1990). Analytical expressions for the tonotopic sensory
scale. Journal of the Acoustical Society of America, 88, 97-100.

Zwicker, E. (1961). Subdivision of the audible frequency range into critical
bands. Journal of the Acoustical Society of America, 33, 248.

Wang, S., Sekey, A., & Gersho, A. (1992). An objective measure for
predicting subjective quality of speech coders. IEEE Journal on Selected
Areas in Communications, 10(5), 819-829.

Glasberg, B. R., & Moore, B. C. (1990). Derivation of auditory filter shapes
from notched-noise data. Hearing Research, 47(1-2), 103-138.

Moore, B. C., & Glasberg, B. R. (1983). Suggested formulae for calculating
auditory-filter bandwidths and excitation patterns. The Journal of the
Acoustical Society of America, 74(3), 750-753.

O'Shaughnessy, D. (1987). Speech Communication: Human and Machine.
Addison-Wesley.

Schutte, H.K., & Seidner, W. (1983). Recommendation by the Union of European
Phoniatricians (UEP): Standardizing Voice Area Measurement/Phonetography.
Folia Phoniatrica, 35, 286-288.
```

**BibTeX Entry**:
```bibtex
@article{Traunmuller1990,
  author = {Traunm\"{u}ller, H.},
  title = {Analytical expressions for the tonotopic sensory scale},
  journal = {Journal of the Acoustical Society of America},
  year = {1990},
  volume = {88},
  pages = {97--100}
}

@article{Zwicker1961,
  author = {Zwicker, E.},
  title = {Subdivision of the audible frequency range into critical bands},
  journal = {Journal of the Acoustical Society of America},
  year = {1961},
  volume = {33},
  pages = {248}
}

@article{Wang1992,
  author = {Wang, S. and Sekey, A. and Gersho, A.},
  title = {An objective measure for predicting subjective quality of speech coders},
  journal = {IEEE Journal on Selected Areas in Communications},
  year = {1992},
  volume = {10},
  number = {5},
  pages = {819--829}
}

@article{Glasberg1990,
  author = {Glasberg, B. R. and Moore, B. C.},
  title = {Derivation of auditory filter shapes from notched-noise data},
  journal = {Hearing Research},
  year = {1990},
  volume = {47},
  number = {1-2},
  pages = {103--138}
}

@article{Moore1983,
  author = {Moore, B. C. and Glasberg, B. R.},
  title = {Suggested formulae for calculating auditory-filter bandwidths and
           excitation patterns},
  journal = {The Journal of the Acoustical Society of America},
  year = {1983},
  volume = {74},
  number = {3},
  pages = {750--753}
}

@book{OShaughnessy1987,
  author = {O'Shaughnessy, D.},
  title = {Speech Communication: Human and Machine},
  publisher = {Addison-Wesley},
  year = {1987}
}

@article{Schutte1983,
  author = {Schutte, H. K. and Seidner, W.},
  title = {Recommendation by the Union of European Phoniatricians (UEP):
           Standardizing Voice Area Measurement/Phonetography},
  journal = {Folia Phoniatrica},
  year = {1983},
  volume = {35},
  pages = {286--288}
}
```

---

### 21-24. R/ssff_c_assp_*.R (Multiple files)

**Files**:
- `ssff_c_assp_acfana.R`
- `ssff_c_assp_cepstrum.R`
- `ssff_c_assp_dftSpectrum.R`
- `ssff_c_assp_ksvfo.R`

**Type**: Empty references (TBD)

---

### 25. R/ssff_c_assp_lp_analysis.R

**Function**: `arf_lar_lpc_rfc_ana`
**Type**: Empty references (TBD) - appears 4 times

---

### 26. R/ssff_c_assp_mhspitch.R

**Type**: Empty references (TBD)

---

### 27. R/ssff_cpp_estk_pitchmark.R

**Type**: Empty references (TBD)

---

### 28. R/ssff_cpp_pyin.R

**Type**: Empty references (TBD)

---

### 29. R/ssff_cpp_sptk_reaper_pm.R

**Type**: Empty references (TBD)

---

### 30. R/ssff_cpp_tandem.R

**Type**: Citation references
**Reference**: `\insertRef{hu2010tandem}{superassp}`, `\insertRef{hu2011unvoiced}{superassp}`

**Status**: Uses dynamic bibliography

---

### 31. R/ssff_cpp_yin.R

**Type**: Empty references (TBD)

---

### 32. R/ssff_estk_pda.R

**Type**: Empty references (TBD)

---

### 33. R/ssff_pladdrr_cpps.R

**Function**: `trk_cpps`
**Type**: Journal articles

```
Hillenbrand, J., Cleveland, R. A., & Erickson, R. L. (1994). Acoustic correlates
of breathy vocal quality. JSHR, 37(4), 769-778.

Hillenbrand, J., & Houde, R. A. (1996). Acoustic correlates of breathy vocal quality:
Dysphonic voices and continuous speech. JSHR, 39(2), 311-321.

Heman-Ackah, Y. D., et al. (2003). Cepstral peak prominence: A more reliable measure
of dysphonia. Ann Otol Rhinol Laryngol, 112(4), 324-333.
```

**BibTeX Entry**:
```bibtex
@article{Hillenbrand1994,
  author = {Hillenbrand, J. and Cleveland, R. A. and Erickson, R. L.},
  title = {Acoustic correlates of breathy vocal quality},
  journal = {Journal of Speech and Hearing Research},
  year = {1994},
  volume = {37},
  number = {4},
  pages = {769--778}
}

@article{Hillenbrand1996,
  author = {Hillenbrand, J. and Houde, R. A.},
  title = {Acoustic correlates of breathy vocal quality: Dysphonic voices and continuous speech},
  journal = {Journal of Speech and Hearing Research},
  year = {1996},
  volume = {39},
  number = {2},
  pages = {311--321}
}

@article{HemanAckah2003,
  author = {Heman-Ackah, Y. D. and others},
  title = {Cepstral peak prominence: A more reliable measure of dysphonia},
  journal = {Annals of Otology, Rhinology \& Laryngology},
  year = {2003},
  volume = {112},
  number = {4},
  pages = {324--333}
}
```

---

### 34. R/ssff_pladdrr_vuv.R

**Function**: `trk_vuv`
**Type**: Journal articles

```
Al-Tamimi, J., & Khattab, G. (2015). Acoustic correlates of the voicing contrast
in Lebanese Arabic singleton and geminate stops. JASA, 138(1), 344-360.

Al-Tamimi, J., & Khattab, G. (2018). Acoustic cue weighting in the singleton vs
geminate contrast in Lebanese Arabic: The case of fricative consonants.
Journal of Phonetics, 71, 306-325.
```

**BibTeX Entry**:
```bibtex
@article{AlTamimi2015,
  author = {Al-Tamimi, J. and Khattab, G.},
  title = {Acoustic correlates of the voicing contrast in Lebanese Arabic singleton
           and geminate stops},
  journal = {Journal of the Acoustical Society of America},
  year = {2015},
  volume = {138},
  number = {1},
  pages = {344--360}
}

@article{AlTamimi2018,
  author = {Al-Tamimi, J. and Khattab, G.},
  title = {Acoustic cue weighting in the singleton vs geminate contrast in Lebanese Arabic:
           The case of fricative consonants},
  journal = {Journal of Phonetics},
  year = {2018},
  volume = {71},
  pages = {306--325}
}
```

---

### 35-37. R/ssff_python_aperiodicities.R

**Type**: Empty references (TBD)

---

### 38. R/ssff_python_brouhaha.R

**Function**: `trk_brouhaha`
**Type**: Citation reference
**Reference**: `\insertAllCited{}`

**Status**: Uses dynamic bibliography

---

### 39. R/ssff_python_crepe.R

**Function**: `trk_crepe`
**Type**: Citation reference
**Reference**: `\insertAllCited{}`

**Status**: Uses dynamic bibliography

---

### 40. R/ssff_python_deepformants.R

**Functions**:
- `trk_deepformants` (appears twice)

**Type**: Conference proceeding

```
Dissen, S., & Keshet, J. (2017). DeepFormants: Deep neural network for formant estimation.
Proceedings of Interspeech 2017.
```

**BibTeX Entry**:
```bibtex
@inproceedings{Dissen2017,
  author = {Dissen, S. and Keshet, J.},
  title = {DeepFormants: Deep neural network for formant estimation},
  booktitle = {Proceedings of Interspeech 2017},
  year = {2017}
}
```

---

### 41. R/ssff_python_dv_f0.R

**Type**: Citation references
**Reference**: `\insertCite{Jadoul2018}{superassp}`, `\insertCite{OrozcoArroyave2018}{superassp}`

**Status**: Uses dynamic bibliography

---

### 42. R/ssff_python_dv_formants.R

**Type**: Citation references
**Reference**: `\insertCite{Jadoul2018}{superassp}`, `\insertCite{OrozcoArroyave2018}{superassp}`

**Status**: Uses dynamic bibliography

---

### 43. R/ssff_python_gfmiaif.R

**Type**: Citation reference
**Reference**: `\insertAllCited{}`

**Status**: Uses dynamic bibliography

---

### 44. R/ssff_python_phonet.R

**Type**: Citation reference
**Reference**: `\insertAllCited{}`

**Status**: Uses dynamic bibliography

---

### 45. R/ssff_python_pm_psauce.R

**Function**: `trk_praatsaucep`
**Type**: Conference proceeding + Journal article

```
Iseli, M., & Alwan, A. (2004). An improved correction formula for the estimation
of harmonic magnitudes and its application to open quotient estimation. ICASSP.

Hawks, J. W., & Miller, J. D. (1995). A formant bandwidth estimation procedure
for vowel synthesis. JASA, 97(2), 1343-1344.

Shue, Y.-L., Keating, P., Vicenik, C., & Yu, K. (2011). VoiceSauce: A program
for voice analysis. ICPhS.
```

**BibTeX Entry**:
```bibtex
@inproceedings{Iseli2004,
  author = {Iseli, M. and Alwan, A.},
  title = {An improved correction formula for the estimation of harmonic magnitudes
           and its application to open quotient estimation},
  booktitle = {IEEE ICASSP},
  year = {2004}
}

@article{Hawks1995,
  author = {Hawks, J. W. and Miller, J. D.},
  title = {A formant bandwidth estimation procedure for vowel synthesis},
  journal = {Journal of the Acoustical Society of America},
  year = {1995},
  volume = {97},
  number = {2},
  pages = {1343--1344}
}

@inproceedings{Shue2011,
  author = {Shue, Y.-L. and Keating, P. and Vicenik, C. and Yu, K.},
  title = {VoiceSauce: A program for voice analysis},
  booktitle = {International Congress of Phonetic Sciences},
  year = {2011}
}
```

---

### 46. R/ssff_python_reaper_pm.R

**Type**: Empty references (TBD)

---

### 47. R/ssff_python_sacc.R

**Function**: `trk_sacc`
**Type**: Online resource

```
Ellis, D. P. W. (2014). Subband Autocorrelation Classification (SAcC) pitch tracker.
http://labrosa.ee.columbia.edu/projects/SAcC/
```

**BibTeX Entry**:
```bibtex
@online{Ellis2014,
  author = {Ellis, D. P. W.},
  title = {Subband Autocorrelation Classification (SAcC) pitch tracker},
  year = {2014},
  url = {http://labrosa.ee.columbia.edu/projects/SAcC/}
}
```

---

### 48. R/ssff_python_seenc.R

**Type**: Empty references (TBD)

---

### 49. R/ssff_python_snack_formant.R

**Function**: `trk_snackf`
**Type**: Conference proceeding + Journal article

```
Sjölander, K. & Beskow, J. (2000). "Wavesurfer - an open source speech tool."
In Proc. ICSLP 2000, Beijing, China.

Talkin, D. (1987). "Speech formant trajectory estimation using dynamic programming
with modulated transition costs." J. Acoust. Soc. Am.
```

**BibTeX Entry**:
```bibtex
@inproceedings{Sjolander2000,
  author = {Sj\"{o}lander, K. and Beskow, J.},
  title = {Wavesurfer - an open source speech tool},
  booktitle = {Proceedings of ICSLP 2000},
  address = {Beijing, China},
  year = {2000}
}

@article{Talkin1987,
  author = {Talkin, D.},
  title = {Speech formant trajectory estimation using dynamic programming with
           modulated transition costs},
  journal = {Journal of the Acoustical Society of America},
  year = {1987}
}
```

---

### 50. R/ssff_python_snack_pitch.R

**Type**: Citation reference
**Reference**: `\insertRef{Sjolander2000}{superassp}`

**Status**: Uses dynamic bibliography

---

### 51-53. R/ssff_python_straight_*.R (3 files)

**Files**:
- `ssff_python_straight_f0.R`
- `ssff_python_straight_spec.R`
- `ssff_python_straight_synth.R`

**Type**: Citation references
**Reference**: `\insertRef{Kawahara.1999.SpeechCommunication}{superassp}`

**Status**: Uses dynamic bibliography

---

### 54. R/ssff_python_swiftf0.R

**Function**: `trk_swiftf0`
**Type**: Citation reference
**Reference**: `\insertAllCited{}`

**Status**: Uses dynamic bibliography

---

### 55. R/ssff_python_yaapt.R

**Type**: Empty references (TBD)

---

### 56. R/trk_creak_union.R

**Function**: `trk_creak_union`
**Type**: Citation reference
**Reference**: `\insertAllCited{}`

**Status**: Uses dynamic bibliography

---

### 57. R/trk_egg_f0.R

**Function**: `trk_egg_f0`
**Type**: Citation reference
**Reference**: `\insertAllCited{}`

**Status**: Uses dynamic bibliography

---

### 58. R/trk_formants_tvwlp.R

**Function**: `trk_formants_tvwlp`
**Type**: Journal article

```
Gowda, D., Kadiri, S. R., Mittal, V. K., Gangashetty, S. V., & Yegnanarayana, B. (2015).
Epoch extraction from speech signals using an auto-associative neural network model.
Speech Communication, 69, 50-65.
```

**BibTeX Entry**:
```bibtex
@article{Gowda2015,
  author = {Gowda, D. and Kadiri, S. R. and Mittal, V. K. and Gangashetty, S. V. and
            Yegnanarayana, B.},
  title = {Epoch extraction from speech signals using an auto-associative neural network model},
  journal = {Speech Communication},
  year = {2015},
  volume = {69},
  pages = {50--65}
}
```

---

### 59. R/trk_phonet.R

**Function**: `trk_phonet`
**Type**: Citation reference
**Reference**: `\insertAllCited{}`

**Status**: Uses dynamic bibliography

---

### 60. R/vat_srh.R

**Type**: Citation references
**Reference**:
- `\insertCite{Drugman2011SRH}{superassp}`
- `\insertCite{Kane2013VAT}{superassp}`

**Status**: Uses dynamic bibliography

---

## Summary Statistics

### By Reference Type

| Type | Count | Files |
|------|-------|-------|
| Journal Articles (hardcoded) | 28 | 13 |
| International Standards (ISO) | 3 | 2 |
| Conference Proceedings | 8 | 7 |
| Books/Theses | 2 | 2 |
| Software/Tools | 4 | 4 |
| Online Resources | 2 | 2 |
| Dynamic Bibliography (\insertAllCited{}) | 9 | 9 |
| Dynamic Bibliography (\insertRef{}) | 6 | 6 |
| Dynamic Bibliography (\insertCite{}) | 2 | 2 |
| Empty/TBD | 5 | 5 |

### By Implementation Category

| Category | Files | References |
|----------|-------|-----------|
| **Acoustics/Psychoacoustics** | 3 | ISO 532, ISO 226, multiple journal articles |
| **Voice Quality** | 8 | Journal articles, standards |
| **Pitch Tracking** | 6 | Dynamic bibliography, journal articles |
| **Formant Analysis** | 3 | Journal articles, conference papers |
| **Speech Signal Processing** | 8 | Conference papers, online resources |
| **Statistical/Pharmacological** | 3 | Journal articles |
| **Spectral Analysis** | 5 | Empty or dynamic |
| **Other** | 13 | Varied |

---

## Next Steps

1. **Dynamic Bibliography Files**: Locate and import BibTeX entries from:
   - `inst/bib/` or similar directory
   - Check for `.bib` files referenced by `\insertRef{}` and `\insertAllCited{}` commands

2. **Empty References**: Complete missing references in:
   - `ssff_c_assp_*.R` files
   - `ssff_cpp_*.R` files
   - `ssff_python_aperiodicities.R`
   - `ssff_python_reaper_pm.R`
   - `ssff_python_seenc.R`
   - `ssff_python_yaapt.R`

3. **BibTeX Master File**: Consolidate hardcoded references into a single comprehensive `.bib` file
   - Use consistent formatting
   - Include DOI/URL information where available
   - Organize by category (voice quality, pitch, formants, etc.)

