# Superassp References Extraction - Complete Index

**Generated**: 2026-03-01
**Source**: All R source files in `/Users/frkkan96/Documents/src/superassp/R/`

---

## Quick Navigation

| File | Purpose | Size | Best For |
|------|---------|------|----------|
| [REFERENCES_EXTRACTION.md](#references_extractionmd) | Detailed reference catalog | 27K | Detailed review, finding specific citations |
| [REFERENCES_BIBTEX.bib](#references_bibtexbib) | BibTeX bibliography | 13K | LaTeX/BibTeX compilation, reference management |
| [REFERENCES_INVENTORY.csv](#references_inventorycsv) | Spreadsheet format | 8.4K | Tracking, filtering, bulk operations |
| [REFERENCES_KEYS.txt](#references_keystxt) | Citation key reference | 6.7K | Quick lookups, citation patterns |
| [REFERENCES_SUMMARY.txt](#references_summarytxt) | Statistical overview | 6.7K | Project status, next steps |
| [REFERENCES_INDEX.md](#) | This file | Navigation | Getting started |

---

## REFERENCES_EXTRACTION.md

**Purpose**: Complete extraction of all @references sections with full BibTeX entries

**Contents**:
- All 58 R files with @references documented
- Function names and line numbers
- Original reference text (formatted)
- Complete BibTeX entries ready to use
- Organized by file

**Use this when**:
- You need the exact reference text as written in the source code
- You want to verify a specific citation
- You're creating a BibTeX master bibliography
- You need to identify where each reference comes from

**Example entry**:
```
### 16. R/list_vat.R

Function: lst_vat
Type: Journal article + Thesis

Tsanas, A., Little, M., McSharry, P., & Ramig, L. (2011).
Nonlinear speech analysis algorithms...
Journal of the Royal Society Interface, 8(59), 842-855.

@article{Tsanas2011,
  author = {...},
  title = {...},
  ...
}
```

---

## REFERENCES_BIBTEX.bib

**Purpose**: Complete BibTeX bibliography ready for use in R package documentation

**Contents**:
- 50+ journal articles, conference papers, standards
- Organized by category (voice quality, phonation, etc.)
- Includes ISO standards (ISO 226, ISO 532)
- Software/tool references
- Notes on dynamic references

**Use this when**:
- Setting up roxygen2 bibliography
- Creating package reference documentation
- Building LaTeX/PDF documentation
- Performing citation management

**Format example**:
```bibtex
@article{Alku1992,
  author = {Alku, P.},
  title = {Glottal wave analysis with pitch synchronous iterative
           adaptive inverse filtering},
  journal = {Speech Communication},
  year = {1992},
  volume = {11},
  number = {2-3},
  pages = {109--118}
}
```

**Note**: Includes placeholder comments for dynamic references that require locating .bib files

---

## REFERENCES_INVENTORY.csv

**Purpose**: Spreadsheet-compatible inventory for tracking and management

**Columns**:
- File: R source file path
- Function: Associated function name
- Reference Type: Journal Article, Conference Paper, Standard, etc.
- Authors/Title: Citation author names or title
- Year: Publication year
- Journal/Conference: Venue
- Volume/Issue: Journal volume and issue
- Pages: Page range
- DOI/URL: Digital Object Identifier or web address
- Status: Hardcoded, Dynamic, or TBD

**Use this when**:
- Importing into spreadsheet applications (Excel, LibreOffice, Google Sheets)
- Tracking completion status of missing references
- Filtering by reference type (e.g., show only standards)
- Performing batch operations on references
- Creating literature review documents

**Example rows**:
```csv
File,Function,Reference Type,Authors/Title,Year,Journal/Conference,Volume/Issue,Pages,DOI/URL,Status
R/covarep_vq.R,lst_covarep_vq,Journal Article,Alku P.,1992,Speech Communication,11(2-3),109-118,,Hardcoded
R/iso226_phon.R,iso226_phon,ISO Standard,ISO 226:2023,2023,International Organization for Standardization,,,ISO 226:2023,Hardcoded
```

---

## REFERENCES_KEYS.txt

**Purpose**: Quick reference guide for citation keys and dynamic bibliography lookups

**Sections**:
1. All citation keys organized by research domain
2. Dynamic references that require .bib file lookup
3. Quality assessment (complete vs incomplete)
4. Citation usage guide and examples
5. Next steps for completing missing references

**Use this when**:
- Looking up a citation key to use in \insertRef{key}{superassp}
- Identifying which .bib file contains a specific reference
- Finding all references in a particular domain (voice quality, pitch, etc.)
- Understanding the citation patterns used in the package
- Tracking incomplete references

**Example**:
```
Voice Quality and Phonation Research:
  Alku1992                 - Glottal wave analysis (IAIF method)
  Kane2013                 - Glottal closure instant detection
  Iseli2004                - Harmonic magnitude estimation

Dynamic References (Requiring .bib file lookup):
  - brouhaha              → Search for brouhaha.bib
  - crepe                 → Search for crepe.bib
```

---

## REFERENCES_SUMMARY.txt

**Purpose**: High-level overview with statistics and recommendations

**Sections**:
1. Overview statistics
2. Reference breakdown by type
3. Files created and purpose
4. Key findings (most referenced authors/topics)
5. Quality assessment
6. Next steps and action items

**Use this when**:
- Getting a bird's eye view of the project
- Understanding what work remains
- Identifying priority areas for completion
- Planning next steps
- Reporting project status

**Statistics provided**:
- Total files analyzed: 196
- Files with @references: 58
- Hardcoded references: ~45
- Dynamic references: 17
- Empty/TBD: 5+

---

## Quick Reference Tables

### Reference Types Distribution

| Type | Count | Percentage |
|------|-------|-----------|
| Journal Articles | 28 | 60% |
| Conference Papers | 8 | 17% |
| Standards (ISO) | 3 | 6% |
| Software/Tools | 4 | 9% |
| Books/Theses | 2 | 4% |
| Online Resources | 2 | 4% |

### Major Research Domains

| Domain | Files | Key References |
|--------|-------|-----------------|
| Voice Quality | 8 | Maryn, Hillenbrand, Kane |
| Pitch Tracking | 6 | Ellis, Villarubia, Gowda |
| Psychoacoustics | 5 | Moore, Glasberg, Zwicker |
| Spectral Analysis | 5 | Talkin, Sjölander, Dissen |
| Formants | 3 | Dissen, Sjölander, Hawks |

### Files Requiring Completion

**Empty/TBD References** (15 files):
- ssff_c_assp_*.R (5 files)
- ssff_cpp_*.R (5 files)
- ssff_python_*.R (5 files)

**Dynamic References Needing .bib Files** (17 files):
- brouhaha, crepe, phonet, gfmiaif, swiftf0, etc.

---

## How to Use These Files

### For Documentation Writers
1. Use **REFERENCES_EXTRACTION.md** to understand existing citations
2. Check **REFERENCES_KEYS.txt** for available citation keys
3. Use **REFERENCES_BIBTEX.bib** for reference entries in documentation

### For Developers Adding New Functions
1. Check **REFERENCES_INVENTORY.csv** for similar functions
2. Look up appropriate references in **REFERENCES_KEYS.txt**
3. Add new @references section with citation keys from **REFERENCES_BIBTEX.bib**

### For Project Managers
1. Review **REFERENCES_SUMMARY.txt** for project status
2. Use **REFERENCES_INVENTORY.csv** to track completion
3. Update status column as references are completed

### For Bibliography Management
1. Import **REFERENCES_BIBTEX.bib** into BibTeX manager
2. Use **REFERENCES_KEYS.txt** to understand key naming conventions
3. Reference entries in R documentation using \insertRef{key}{superassp}

### For Package Maintenance
1. Track **REFERENCES_INVENTORY.csv** for missing references
2. Use **REFERENCES_SUMMARY.txt** to identify priority areas
3. Update **REFERENCES_BIBTEX.bib** as new references are added

---

## Key Findings

### Strengths
- ✓ 28 journal articles fully documented
- ✓ ISO standards properly formatted
- ✓ Voice quality research well-represented
- ✓ Consistent citation patterns in psychoacoustics
- ✓ Software references with URLs

### Gaps to Address
- ✗ 15 functions with empty @references
- ✗ 17 functions with dynamic bibliography (pending .bib file location)
- ✗ Some missing DOI information
- ✗ Some author names need verification (special characters)

### Next Actions (Priority Order)
1. **Locate dynamic .bib files** (inst/bib/ directory)
2. **Complete empty references** (15 files)
3. **Add DOI information** where available
4. **Verify author names** and special characters
5. **Test roxygen2 compilation** with full bibliography

---

## Statistics Summary

| Metric | Value |
|--------|-------|
| Files analyzed | 196 R files |
| Files with references | 58 (29.6%) |
| Hardcoded references | ~45 entries |
| Dynamic references | 17 functions |
| ISO standards | 3 |
| Journal articles | 28 |
| Conference papers | 8 |
| Total citation keys | 50+ |
| Empty @references | 5+ |
| Needs .bib files | 17 |

---

## File Locations

All extracted reference files are located in:
```
/Users/frkkan96/Documents/src/superassp/
```

Files created:
- `REFERENCES_EXTRACTION.md` - Detailed reference catalog (27K)
- `REFERENCES_BIBTEX.bib` - BibTeX bibliography (13K)
- `REFERENCES_INVENTORY.csv` - Spreadsheet inventory (8.4K)
- `REFERENCES_KEYS.txt` - Citation key reference (6.7K)
- `REFERENCES_SUMMARY.txt` - Statistical overview (6.7K)
- `REFERENCES_INDEX.md` - This navigation guide

---

## Related Documentation

For context about the superassp package, see:
- `CLAUDE.md` - Project architecture and conventions
- `README.md` - User documentation
- `NEWS.md` - Version history and changes
- `DESCRIPTION` - Package metadata

For dynamic bibliography files, search for:
- `inst/bib/` directory
- `*.bib` files in project root
- BibTeX entries in roxygen2 comments

---

## Recommendations

### For Immediate Use
1. Use **REFERENCES_BIBTEX.bib** as master bibliography source
2. Import **REFERENCES_INVENTORY.csv** into tracking spreadsheet
3. Reference **REFERENCES_KEYS.txt** for citation key naming

### For Long-term Maintenance
1. Keep **REFERENCES_EXTRACTION.md** updated as functions are added
2. Update **REFERENCES_INVENTORY.csv** status column monthly
3. Consolidate all .bib files into single **REFERENCES_BIBTEX.bib**
4. Test roxygen2 compilation after any bibliography changes

### For Package Documentation
1. Use \insertRef{key}{superassp} for citations in @references
2. Follow citation format examples in REFERENCES_EXTRACTION.md
3. Verify all new references are added to REFERENCES_BIBTEX.bib
4. Include DOI information where available

---

**Generated**: 2026-03-01
**Source**: Automated extraction from R source files
**Status**: Complete - Awaiting dynamic .bib file integration

