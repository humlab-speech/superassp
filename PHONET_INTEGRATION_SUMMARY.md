# Phonet Integration Summary - superassp

## Overview

Complete integration of Phonet phonological posterior extraction into the superassp R package, providing both list-based (`lst_phonet`) and SSFF track-based (`trk_phonet`) interfaces.

## Files Created/Modified

### New Files in superassp

1. **`R/install_phonet.R`** - Installation helpers
   - `install_phonet()` - Installs Python dependencies including tf-keras
   - `phonet_available()` - Checks if phonet is installed
   - `phonet_info()` - Returns configuration information

2. **`R/ssff_python_phonet.R`** - List-based interface
   - `lst_phonet()` - Returns list/data.frame format
   - For feature extraction and data analysis

3. **`R/trk_phonet.R`** - SSFF track interface ⭐ NEW
   - `trk_phonet()` - Returns SSFF track objects
   - Compatible with emuR framework
   - Supports toFile parameter for batch processing

4. **`tests/test_phonet_integration.R`** - Comprehensive test suite
   - Tests for both lst_ and trk_ functions
   - SSFF file I/O tests
   - Consistency tests between lst and trk outputs

### Modified Files in phonet

5. **`phonet/phonet.py`** - Python 3.12+ compatibility fix
   - Added tf-keras compatibility (lines 20-23)
   - Backward compatible with all Python versions

6. **`CLAUDE.md`** - Updated documentation
   - Python 3.12+ tf-keras requirements
   - R integration guide (200+ lines)
   - Troubleshooting section

7. **`R_INTEGRATION_GUIDE.md`** - Complete R integration documentation

8. **`TESTING_REPORT.md`** - Comprehensive test results

## Function Comparison: lst_phonet vs trk_phonet

### lst_phonet() - List/Data Frame Output

**Purpose**: Feature extraction and data analysis

**Returns**:
```r
list(
  time = numeric vector,
  phoneme = character vector,
  vocalic = numeric vector,
  nasal = numeric vector,
  ...
  file = character
)
```

**Use Cases**:
- Statistical analysis of phonological features
- Feature vectors for machine learning
- Quick exploration of phonological patterns
- Integration with tidyverse workflows

**Example**:
```r
result <- lst_phonet("speech.wav", classes = c("nasal", "stop"))
library(ggplot2)
ggplot(as.data.frame(result), aes(x = time, y = nasal)) +
  geom_line()
```

### trk_phonet() - SSFF Track Output ⭐

**Purpose**: Time-aligned phonological annotation

**Returns**:
- `toFile=TRUE`: Number of files processed
- `toFile=FALSE`: AsspDataObj with SSFF tracks

**Use Cases**:
- Integration with emuR databases
- Time-aligned phonological annotation
- Phonetic segmentation
- Multi-tier annotation workflows
- Publication-quality visualizations

**Example**:
```r
# Return AsspDataObj
result <- trk_phonet("speech.wav",
                     classes = c("vocalic", "consonantal"),
                     toFile = FALSE)

# Access tracks
vocalic_track <- result$vocalic  # matrix
time_vector <- seq(attr(result, "startTime"),
                  by = 1/attr(result, "sampleRate"),
                  length.out = nrow(result$vocalic))

# Write SSFF files for emuR
trk_phonet("speech.wav", classes = "all", toFile = TRUE)
# Creates: speech.phn (18 tracks, 100 Hz)
```

## Technical Specifications

### Phonet Processing
- **Input**: Any audio format (via av package)
- **Internal**: Resampled to 16 kHz mono WAV
- **Feature Extraction**: 33 Mel-filterbanks + energy (34D)
- **Frame Rate**: 100 Hz (10ms intervals)
- **Model**: 2-layer BGRU (128 units each)
- **Output**: 18 phonological posterior probabilities (0-1)

### SSFF Format Details
- **Sample Rate**: 100.0 Hz
- **Track Format**: REAL32
- **File Extension**: `.phn`
- **Tracks**: One per phonological class
- **Compatible with**: emuR, EMU-SDMS, wrassp

## Phonological Classes (18 total)

### Vowel Features
- `vocalic` - Vowel sounds
- `back` - Back vowels (/a/, /o/, /u/)
- `anterior` - Front vowels (/e/, /i/)
- `open` - Open vowels (/a/, /e/, /o/)
- `close` - Close vowels (/i/, /u/)

### Consonant Features
- `consonantal` - Consonant sounds
- `nasal` - Nasal consonants (/m/, /n/)
- `stop` - Stop consonants (/p/, /b/, /t/, /d/, /k/, /g/)
- `continuant` - Continuant sounds
- `lateral` - Lateral sounds (/l/)
- `flap` - Flap/tap sounds (/ɾ/)
- `trill` - Trill sounds (/r/)
- `voice` - Voiced sounds
- `strident` - Strident sounds (/f/, /s/, /ʃ/)

### Place of Articulation
- `labial` - Labial sounds
- `dental` - Dental sounds
- `velar` - Velar sounds

### Other
- `pause` - Silence/pauses

## Installation

### For Users

```r
library(superassp)

# Install Phonet dependencies
install_phonet()

# Check installation
phonet_available()
phonet_info()
```

### For Developers

```r
# In phonet repository
devtools::document()
devtools::test()

# In superassp repository
devtools::document()
devtools::load_all()

# Test integration
source("tests/test_phonet_integration.R")
```

## Python 3.12+ Compatibility Fix

**Issue**: Keras 3.x removed `sample_weight_mode` parameter

**Solution**: Use `tf-keras` for Keras 2.x API
```python
# phonet/phonet.py lines 20-23
try:
    import tf_keras as keras
except ImportError:
    from tensorflow import keras
```

**Installation**:
```bash
pip install tensorflow tf-keras
pip install git+https://github.com/jcvasquezc/phonet.git
```

**R Installation**:
```r
reticulate::py_install(c("tensorflow", "tf-keras"))
reticulate::py_install("git+https://github.com/jcvasquezc/phonet.git")
```

## Usage Examples

### Example 1: Basic Phonological Analysis

```r
library(superassp)

# Track nasality and stops
result <- trk_phonet(
  "speech.wav",
  classes = c("nasal", "stop"),
  toFile = FALSE
)

# Plot posteriors
time <- seq(attr(result, "startTime"),
           by = 1/attr(result, "sampleRate"),
           length.out = nrow(result$nasal))

plot(time, result$nasal, type = "l", col = "blue",
     ylim = c(0, 1), ylab = "Posterior Probability",
     main = "Phonological Posteriors")
lines(time, result$stop, col = "red")
legend("topright", c("Nasal", "Stop"), col = c("blue", "red"), lty = 1)
```

### Example 2: emuR Integration

```r
library(emuR)
library(superassp)

# Process all files in emuR database
db_path <- "/path/to/emuDB"
db <- load_emuDB(db_path)

# Get all audio files
files <- list_files(db)

# Track all phonological classes to SSFF
trk_phonet(
  files$absolute_file_path,
  classes = "all",
  outputDirectory = file.path(db_path, "bundle_dir"),
  toFile = TRUE
)

# Now phonological posteriors can be queried in emuR
query(db, "Phonetic=n")  # Query nasal segments
```

### Example 3: Tidyverse Analysis

```r
library(superassp)
library(tidyverse)

# Extract as list
result <- lst_phonet(
  "speech.wav",
  classes = c("vocalic", "consonantal", "nasal", "stop")
)

# Convert to tibble
df <- as_tibble(result)

# Analysis
df %>%
  mutate(
    vowel = vocalic > 0.5,
    nasal_consonant = consonantal > 0.5 & nasal > 0.5
  ) %>%
  summarize(
    vowel_proportion = mean(vowel),
    nasal_proportion = mean(nasal_consonant)
  )

# Visualization
df %>%
  pivot_longer(c(vocalic, consonantal, nasal, stop),
              names_to = "class", values_to = "posterior") %>%
  ggplot(aes(x = time, y = posterior, color = class)) +
  geom_line() +
  facet_wrap(~class, ncol = 1) +
  theme_minimal() +
  labs(title = "Phonological Posteriors Over Time")
```

### Example 4: Batch Processing

```r
library(superassp)

# Process directory of files
files <- list.files("audio_dir", pattern = "\\.wav$", full.names = TRUE)

# Write SSFF tracks
num_processed <- trk_phonet(
  files,
  classes = "all",
  outputDirectory = "results/",
  toFile = TRUE,
  verbose = TRUE
)

cat(sprintf("Processed %d/%d files\n", num_processed, length(files)))
```

## Performance

- **Model Loading**: ~2-3 seconds (one-time, reused across files)
- **Processing**: ~0.5-1x real-time
  - 3.14s audio in ~1-2 seconds
- **Memory**: ~500 MB (TensorFlow models)
- **Batch Efficiency**: Reuses Phonet instance automatically

## Testing

All tests pass successfully:

✅ Installation and availability checks
✅ Basic audio processing (lst_phonet)
✅ All 18 phonological classes
✅ Time windowing
✅ Multiple file processing
✅ SSFF track output (trk_phonet)
✅ SSFF file I/O
✅ Consistency between lst and trk functions
✅ Input validation

## Known Limitations

1. **Optimizer Warning** (non-critical):
   - Adam optimizer runs slowly on Apple Silicon
   - Can be addressed by using `legacy.Adam`

2. **Model Language**:
   - Pre-trained on Spanish speech
   - Performance may vary on other languages
   - Retraining possible (see phonet/train/)

3. **Audio Requirements**:
   - Best results with 16 kHz mono WAV
   - Other formats auto-converted via av package

## Future Enhancements

Potential additions:
- [ ] Confidence intervals for posteriors
- [ ] Multi-language model support
- [ ] Real-time streaming mode
- [ ] GPU acceleration option
- [ ] Integration with other superassp features (formants, pitch)

## References

Vásquez-Correa, J. C., Klumpp, P., Orozco-Arroyave, J. R., & Nöth, E. (2019).
Phonet: A Tool Based on Gated Recurrent Neural Networks to Extract Phonological
Posteriors from Speech. Proc. Interspeech 2019, 549-553.
doi:10.21437/Interspeech.2019-1405

## Summary

✅ **Complete integration** of Phonet into superassp
✅ **Two interfaces**: lst_phonet() and trk_phonet()
✅ **Python 3.12+ compatible** with tf-keras
✅ **Fully tested** with comprehensive test suite
✅ **emuR compatible** via SSFF format
✅ **Production ready** for research and applications

Both `lst_phonet()` and `trk_phonet()` provide complementary access to phonological posteriors, serving different use cases while sharing the same underlying Phonet BGRU models.
