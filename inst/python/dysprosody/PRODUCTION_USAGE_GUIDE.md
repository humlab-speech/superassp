# Production Usage Guide - Optimized Dysprosody

## Quick Start

### Installation

```bash
# Install dependencies
pip install parselmouth numpy pandas scipy

# No additional dependencies needed!
```

### Basic Usage

```python
from dysprosody_optimized import prosody_measures

# Analyze single file
features = prosody_measures("audio.wav")

# Access features
print(f"Duration: {features['Duration']:.2f}s")
print(f"Pitch Mean: {features['PitchMean']:.0f} Hz")
print(f"INTSINT Labels: {features['IntsIntLabels']:.0f}")

# Save to CSV
features.to_csv("results.csv")
```

## Batch Processing (Recommended)

### For M1 Pro (10 CPU cores)

```python
from dysprosody_optimized import batch_process
import glob
import pandas as pd

# Find all WAV files
files = glob.glob("corpus/**/*.wav", recursive=True)

# Process in parallel (use 8 of 10 cores, leave 2 for OS)
results = batch_process(files, max_workers=8)

# Convert to DataFrame
df = pd.DataFrame(results).T  # Rows = files, Columns = features
df.to_csv("prosody_results.csv")

print(f"Processed {len(results)} files")
print(f"Features per file: {len(df.columns)}")
```

**Expected Performance:**
- ~6 files/second
- 100 files in ~17 seconds
- 1000 files in ~2.8 minutes

### For AMD EPYC VM (4 cores exposed)

```python
from dysprosody_optimized import batch_process
import glob
import pandas as pd

# Find all WAV files
files = glob.glob("corpus/**/*.wav", recursive=True)

# Process in parallel (use all 4 available cores)
results = batch_process(files, max_workers=4)

# Convert to DataFrame
df = pd.DataFrame(results).T
df.to_csv("prosody_results.csv")

print(f"Processed {len(results)} files")
```

**Expected Performance:**
- ~4-5 files/second
- 100 files in ~20-25 seconds
- 1000 files in ~3.3-4.2 minutes

## Advanced Usage

### Progress Tracking

```python
from dysprosody_optimized import prosody_measures
import glob
from tqdm import tqdm

files = glob.glob("*.wav")

results = {}
for file in tqdm(files, desc="Processing"):
    result = prosody_measures(file)
    if result is not None:
        results[file] = result

df = pd.DataFrame(results).T
df.to_csv("results.csv")
```

### Error Handling

```python
from dysprosody_optimized import batch_process
import glob

files = glob.glob("*.wav")

# batch_process handles errors internally
results = batch_process(files, max_workers=4)

# Check for failures
processed = len(results)
total = len(files)
failed = total - processed

print(f"Processed: {processed}/{total}")
if failed > 0:
    print(f"Failed: {failed} files")
```

### Custom Parameters

```python
from dysprosody_optimized import prosody_measures

# Custom F0 range for different speakers
# Male speakers
features = prosody_measures("male_speaker.wav", minF=75, maxF=400)

# Female speakers
features = prosody_measures("female_speaker.wav", minF=100, maxF=600)

# Children
features = prosody_measures("child.wav", minF=150, maxF=800)
```

### Filter by Duration

```python
from dysprosody_optimized import prosody_measures
import glob
import parselmouth

# Only process files longer than 2 seconds
files = glob.glob("*.wav")
long_files = []

for file in files:
    sound = parselmouth.Sound(file)
    duration = sound.get_total_duration()
    if duration >= 2.0:
        long_files.append(file)

results = batch_process(long_files, max_workers=8)
```

## Integration Examples

### With Existing Research Pipeline

```python
import pandas as pd
from dysprosody_optimized import batch_process

# Load metadata
metadata = pd.read_csv("participants.csv", index_col="participant_id")

# Get audio files
files = metadata['audio_path'].tolist()

# Process prosody features
prosody_results = batch_process(files, max_workers=8)

# Merge with metadata
prosody_df = pd.DataFrame(prosody_results).T
prosody_df.index = prosody_df.index.map(lambda x: x.removesuffix('.wav'))

combined = metadata.join(prosody_df, how='inner')
combined.to_csv("combined_results.csv")
```

### With Statistical Analysis

```python
import pandas as pd
import scipy.stats as stats
from dysprosody_optimized import batch_process

# Process two groups
control_files = glob.glob("control/**/*.wav", recursive=True)
patient_files = glob.glob("patients/**/*.wav", recursive=True)

control_results = batch_process(control_files, max_workers=8)
patient_results = batch_process(patient_files, max_workers=8)

control_df = pd.DataFrame(control_results).T
patient_df = pd.DataFrame(patient_results).T

# Compare pitch mean between groups
t_stat, p_value = stats.ttest_ind(
    control_df['PitchMean'],
    patient_df['PitchMean']
)

print(f"Pitch Mean Comparison:")
print(f"  Control: {control_df['PitchMean'].mean():.1f} Hz")
print(f"  Patients: {patient_df['PitchMean'].mean():.1f} Hz")
print(f"  t-statistic: {t_stat:.3f}, p-value: {p_value:.4f}")
```

### With Machine Learning

```python
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from dysprosody_optimized import batch_process

# Process all files
files = glob.glob("data/**/*.wav", recursive=True)
prosody_results = batch_process(files, max_workers=8)

# Prepare data
df = pd.DataFrame(prosody_results).T

# Extract labels from filenames (assuming format: label_id.wav)
df['label'] = df.index.map(lambda x: x.split('_')[0])

# Features (exclude metadata)
feature_cols = [col for col in df.columns if col not in ['label', 'Duration']]
X = df[feature_cols]
y = df['label']

# Train classifier
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)

clf = RandomForestClassifier(n_estimators=100)
clf.fit(X_train, y_train)

score = clf.score(X_test, y_test)
print(f"Accuracy: {score:.2%}")

# Feature importance
importances = pd.Series(clf.feature_importances_, index=feature_cols)
top_features = importances.nlargest(10)
print("\nTop 10 most important features:")
print(top_features)
```

## Performance Optimization Tips

### 1. Use Batch Processing

```python
# ❌ Slow - Sequential
results = {}
for file in files:
    results[file] = prosody_measures(file)

# ✅ Fast - Parallel
results = batch_process(files, max_workers=8)
```

### 2. Optimal Worker Count

```python
import os

# M1 Pro (10 cores)
max_workers = min(8, len(files))  # Don't spawn more workers than files

# AMD EPYC VM (4 cores)
max_workers = min(4, len(files))

# Or automatically detect
max_workers = min(os.cpu_count() - 2, len(files))  # Leave 2 cores for OS
```

### 3. Pre-filter Files

```python
# Filter out very short files before processing
import parselmouth

valid_files = []
for file in all_files:
    try:
        sound = parselmouth.Sound(file)
        if sound.get_total_duration() >= 1.0:
            valid_files.append(file)
    except:
        pass

results = batch_process(valid_files, max_workers=8)
```

### 4. Process in Chunks (Large Corpora)

```python
import numpy as np

# For very large corpora (>10,000 files)
files = glob.glob("huge_corpus/**/*.wav", recursive=True)

# Process in chunks to avoid memory issues
chunk_size = 1000
chunks = [files[i:i+chunk_size] for i in range(0, len(files), chunk_size)]

all_results = {}
for i, chunk in enumerate(chunks):
    print(f"Processing chunk {i+1}/{len(chunks)}")
    chunk_results = batch_process(chunk, max_workers=8)
    all_results.update(chunk_results)

    # Save checkpoint
    pd.DataFrame(chunk_results).T.to_csv(f"results_chunk_{i}.csv")

# Combine all results
df = pd.DataFrame(all_results).T
df.to_csv("complete_results.csv")
```

## Troubleshooting

### Issue: "Skipping utterance with less than 1 second duration"

**Solution:** This is by design. Files < 1 second are too short for meaningful prosody analysis.

```python
# Pre-filter files
import parselmouth

def filter_by_duration(files, min_duration=1.0):
    valid = []
    for file in files:
        sound = parselmouth.Sound(file)
        if sound.get_total_duration() >= min_duration:
            valid.append(file)
    return valid

valid_files = filter_by_duration(all_files)
results = batch_process(valid_files, max_workers=8)
```

### Issue: High memory usage with many files

**Solution:** Process in chunks (see Advanced Usage above)

### Issue: Process slower than expected

**Check:**
1. Too many workers for available cores?
2. I/O bottleneck (slow disk)?
3. Files on network drive?
4. Other processes consuming CPU?

```python
# Reduce workers if CPU usage is high but throughput is low
results = batch_process(files, max_workers=4)  # Try fewer workers
```

## Feature Reference

### Metadata Features (7)
- `Duration` - Audio duration (seconds)
- `PitchKey` - Optimal INTSINT key (Hz)
- `PitchRange` - Optimal pitch range (octaves)
- `PitchMean` - Mean pitch from MOMEL (Hz)
- `IntsIntLabels` - Number of INTSINT labels
- `UniqueIntsInt` - Number of unique tone types
- `IntsIntConcentration` - Labels per second

### Spectral Features (~36)
- `L2L1`, `L2cL1c` - H1-H2 measures
- `L1cLF3c`, `L1LF3` - H1*-A3* measures
- `SLF` - Spectral tilt
- `C1` - First MFCC coefficient
- `Spectral Balance` - Energy ratio
- `SLF6D.1` through `SLF6D.6` - Polynomial coefficients

### Statistical Summaries (~150)
For each feature:
- `_tstd` - Standard deviation
- `_tmean` - Mean
- `_variation` - Coefficient of variation
- `_iqr` - Interquartile range
- `_tmax` - Maximum
- `_tmin` - Minimum

**Total: 193 features**

## Citation

If you use this code in research, please cite:

```bibtex
@article{dysprosody2025,
  title={Dysprosody in Parkinson's Disease},
  journal={Frontiers in Human Neuroscience},
  year={2025},
  doi={10.3389/fnhum.2025.1566274}
}
```

## Support

For issues or questions:
1. Check this guide
2. Review `FINAL_OPTIMIZATION_SUMMARY.md`
3. Check `CLAUDE.md` for architecture details
4. Review example code in `demo_pure.py`

## Version Information

- **dysprosody_pure.py** - v1.0 - Pure Python baseline
- **dysprosody_optimized.py** - v2.0 - ⭐ Recommended for production
- **demo_pure.py** - CLI tool for both versions

**Recommended:** Use `dysprosody_optimized.py` for all production workloads.
