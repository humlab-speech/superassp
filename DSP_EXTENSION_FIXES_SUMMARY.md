# DSP Extension Attribute Fixes - Quick Reference

**Date:** 2025-10-29
**Total Issues:** 5-6 functions need fixes
**Estimated Time:** 2-3 hours

---

## 🔴 HIGH PRIORITY - Fix Immediately (15 minutes)

### 1. Fix `lst_eGeMAPS()` Extension Mismatch

**File:** `R/list_python_opensmile_eGeMAPS.R`

**Line 66 - CURRENT (WRONG):**
```r
lst_eGeMAPS <- function(listOfFiles,
                   beginTime=0,
                   endTime=0,
                   explicitExt="ocp",  # ❌ WRONG
                   use_cpp = TRUE){
```

**Line 66 - FIX TO:**
```r
lst_eGeMAPS <- function(listOfFiles,
                   beginTime=0,
                   endTime=0,
                   explicitExt="ogs",  # ✅ CORRECT
                   use_cpp = TRUE){
```

**Reason:** Attribute correctly uses `"ogs"` (line 125), parameter should match.

---

### 2. Fix `lst_emobase()` Extension Mismatch

**File:** `R/list_python_opensmile_emobase.R`

**Line 26 - CURRENT (WRONG):**
```r
lst_emobase <- function(listOfFiles,
                   beginTime=0,
                   endTime=0,
                   explicitExt="ocp",  # ❌ WRONG
                   use_cpp = TRUE, verbose = FALSE){
```

**Line 26 - FIX TO:**
```r
lst_emobase <- function(listOfFiles,
                   beginTime=0,
                   endTime=0,
                   explicitExt="emo",  # ✅ CORRECT
                   use_cpp = TRUE, verbose = FALSE){
```

**Reason:** Attribute correctly uses `"emo"` (line 87), parameter should match.

---

## 🟡 MEDIUM PRIORITY - Fix Soon (1-2 hours)

### 3. Add Missing Attribute to `trk_creak_union()`

**File:** `R/trk_creak_union.R`

**Current:** Function has `explicitExt = "crk"` parameter but missing attribute.

**Add at end of file:**
```r
# Set function attributes
attr(trk_creak_union, "ext") <- "crk"
attr(trk_creak_union, "tracks") <- c("AM_creak", "CD_creak", "union_creak")
attr(trk_creak_union, "outputType") <- "SSFF"
```

---

### 4. Add Extension Support to `trk_dv_f0()`

**File:** `R/ssff_python_dv_f0.R`

**Current:** No `explicitExt` parameter, no attribute.

**Step 1 - Update function signature:**
```r
trk_dv_f0 <- function(audio_path,
                     frame_shift = 10,
                     min_f0 = 75,
                     max_f0 = 600,
                     include_voicing = TRUE,
                     output_format = "AsspDataObj",
                     explicitExt = "dvf",        # ADD THIS
                     toFile = FALSE,             # ADD THIS
                     outputDirectory = NULL,     # ADD THIS
                     ...)
```

**Step 2 - Add file output logic (if not present):**
```r
# Inside function body, add:
if (toFile) {
  output_file <- .construct_output_path(audio_path, explicitExt, outputDirectory)
  write.AsspDataObj(result, output_file)
}
```

**Step 3 - Add attributes at end of file:**
```r
attr(trk_dv_f0, "ext") <- "dvf"
attr(trk_dv_f0, "tracks") <- c("F0[Hz]", "voicing")
attr(trk_dv_f0, "outputType") <- "SSFF"
```

---

### 5. Add Extension Support to `trk_dv_formants()`

**File:** `R/ssff_python_dv_formants.R`

**Current:** No `explicitExt` parameter, no attribute.

**Step 1 - Update function signature:**
```r
trk_dv_formants <- function(audio_path,
                           frame_shift = 10,
                           max_formant = 5500,
                           n_formants = 4,
                           output_format = "AsspDataObj",
                           explicitExt = "dvfm",      # ADD THIS
                           toFile = FALSE,            # ADD THIS
                           outputDirectory = NULL,    # ADD THIS
                           ...)
```

**Step 2 - Add file output logic (if not present):**
```r
# Inside function body, add:
if (toFile) {
  output_file <- .construct_output_path(audio_path, explicitExt, outputDirectory)
  write.AsspDataObj(result, output_file)
}
```

**Step 3 - Add attributes at end of file:**
```r
attr(trk_dv_formants, "ext") <- "dvfm"
attr(trk_dv_formants, "tracks") <- c("fm", "bw")
attr(trk_dv_formants, "outputType") <- "SSFF"
```

---

### 6. Investigate and Fix `trk_formants_tvwlp()`

**File:** `R/trk_formants_tvwlp.R`

**Action Required:**
1. Read the file to understand implementation
2. Check if it supports file output (`toFile` parameter)
3. If yes, add appropriate `ext` attribute (suggest `"tvw"`)

**Provisional Fix:**
```r
# If function supports file output, add at end of file:
attr(trk_formants_tvwlp, "ext") <- "tvw"
attr(trk_formants_tvwlp, "tracks") <- c("fm", "bw")
attr(trk_formants_tvwlp, "outputType") <- "SSFF"
```

---

## Testing Checklist

### After High Priority Fixes

```r
# Load package
devtools::load_all()

# Test file
test_file <- system.file("samples", "sustained", "a1.wav", package = "superassp")

# Test lst_eGeMAPS
result1 <- lst_eGeMAPS(test_file)
# ✅ Check: No errors
# ✅ Check: If file written, extension is .ogs (NOT .ocp)

# Test lst_emobase
result2 <- lst_emobase(test_file)
# ✅ Check: No errors
# ✅ Check: If file written, extension is .emo (NOT .ocp)

# Verify attributes
stopifnot(attr(lst_eGeMAPS, "ext") == "ogs")
stopifnot(attr(lst_emobase, "ext") == "emo")
```

### After Medium Priority Fixes

```r
# Test trk_creak_union
result3 <- trk_creak_union(test_file, toFile = TRUE)
# ✅ Check: File written with .crk extension
# ✅ Check: attr(trk_creak_union, "ext") == "crk"

# Test trk_dv_f0 (if refactored)
result4 <- trk_dv_f0(test_file, toFile = TRUE)
# ✅ Check: File written with .dvf extension
# ✅ Check: attr(trk_dv_f0, "ext") == "dvf"

# Test trk_dv_formants (if refactored)
result5 <- trk_dv_formants(test_file, toFile = TRUE)
# ✅ Check: File written with .dvfm extension
# ✅ Check: attr(trk_dv_formants, "ext") == "dvfm"
```

---

## Files That Are Correct (No Changes Needed)

### ✅ Summary Functions (Intentionally No ext)
- `list_vat.R` → `lst_vat()`
- `lst_voice_sauce.R` → `lst_voice_sauce()`
- `list_dysprosody.R` → `lst_dysprosody()`
- `list_voxit.R` → `lst_voxit()`
- `covarep_vq.R` → `lst_covarep_vq()`

**Reason:** Return data frames/lists, not SSFF files.

### ✅ Internal C++ Wrappers (Acceptable, Parent Has ext)
- `list_cpp_opensmile_emobase.R` → `lst_emobase_cpp()`
- `list_cpp_opensmile_generic.R` → `lst_eGeMAPS_cpp()`, `lst_ComParE_2016_cpp()`

**Reason:** Internal functions (`@keywords internal`), public wrappers have ext.

---

## Development Workflow

```bash
# 1. Make high priority fixes
# Edit R/list_python_opensmile_eGeMAPS.R line 66
# Edit R/list_python_opensmile_emobase.R line 26

# 2. Regenerate documentation
Rscript -e "devtools::document()"

# 3. Reload package
Rscript -e "devtools::load_all()"

# 4. Test
Rscript -e "devtools::test()"

# 5. Make medium priority fixes (one at a time)
# Edit R/trk_creak_union.R
# Edit R/ssff_python_dv_f0.R
# Edit R/ssff_python_dv_formants.R
# Investigate R/trk_formants_tvwlp.R

# 6. After each change
Rscript -e "devtools::document()"
Rscript -e "devtools::load_all()"
Rscript -e "devtools::test()"

# 7. Final check
Rscript -e "devtools::check()"
```

---

## Summary

| Priority | Issues | Files | Time | Complexity |
|----------|--------|-------|------|------------|
| 🔴 High | 2 | 2 | 15 min | Low |
| 🟡 Medium | 3-4 | 4 | 2 hrs | Medium |
| **TOTAL** | **5-6** | **6** | **2-3 hrs** | **Low-Med** |

---

## Quick Command Reference

```bash
# Check current extensions
cd R
grep -h 'attr.*"ext"' *.R | sort | uniq

# Find functions without ext
for f in *.R; do
  if grep -qE "^(lst_|trk_)" "$f"; then
    if ! grep -q 'attr.*"ext"' "$f"; then
      echo "$f"
    fi
  fi
done

# Verify a specific function
grep -A 2 -B 2 "explicitExt" list_python_opensmile_eGeMAPS.R
grep "attr.*ext" list_python_opensmile_eGeMAPS.R
```

---

**End of Quick Reference**
