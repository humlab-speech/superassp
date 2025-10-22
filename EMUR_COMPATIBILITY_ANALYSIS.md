# emuR Compatibility Analysis for Track Name Migration

## Overview

This document analyzes the impact of migrating track names to Titze 2015 / Nylén 2024 standards on compatibility with the emuR speech database system.

**emuR**: R package for speech database management and analysis
**Website**: https://ips-lmu.github.io/The-EMU-SDMS-Manual/

## Key Integration Points

### 1. SSFF File Format

**Format**: Simple Signal File Format (SSFF) - binary format for time-series data

**Track name storage**:
- Track names are stored in SSFF file headers
- Each track has a name attribute
- emuR reads these names when loading signal files

**Impact**:
- ✅ **SSFF format is flexible** - supports any track name string
- ✅ **No breaking change** - track names are metadata, not structural
- ⚠️ **Database consistency** - existing databases use old names

### 2. AsspDataObj Class

**Usage**: superassp creates `AsspDataObj` objects, emuR can read them

**Integration**:
```r
# superassp creates
result <- trk_forest("audio.wav", toFile = FALSE)
class(result)  # "AsspDataObj"

# emuR can use
library(emuR)
# emuR functions can process AsspDataObj
```

**Impact**:
- ✅ **No breaking change** - `AsspDataObj` structure remains the same
- ✅ **Attribute-based** - track names in `attr(obj, "tracks")`
- ✅ **Backwards compatible** - alias system handles old names

### 3. Automatic Signal File Discovery

**emuR behavior**:
- Scans database directories for SSFF files
- Matches files by extension (`.f0`, `.fms`, `.dft`, etc.)
- Reads track names from file headers

**Impact**:
- ✅ **Extension-based** - track names don't affect file discovery
- ✅ **Header-based** - new names stored in headers
- ⚠️ **Query compatibility** - queries reference track names

### 4. Track Queries

**emuR usage**:
```r
# Query a specific track
data <- get_trackdata(dbHandle, "f0", "Phonetic")
```

**Impact**:
- ⚠️ **Breaking change** - queries use exact track names
- 🔧 **Mitigation**: Provide migration utility to update databases
- 🔧 **Alternative**: emuR could support aliases

## Compatibility Scenarios

### Scenario 1: New Database (Post-Migration)

**Workflow**:
1. User creates SSFF files with superassp (new track names)
2. User imports files into emuR database
3. User queries tracks using new names

**Status**: ✅ **Fully compatible**

**Example**:
```r
# Generate signal with superassp
library(superassp)
trk_forest("audio.wav", toFile = TRUE)  # Creates file with F1[Hz], F2[Hz], ...

# Import to emuR
library(emuR)
import_ssffTrackDefinition(dbHandle, "audio.fms")

# Query using new names
formants <- get_trackdata(dbHandle, "F1[Hz]", "Phonetic")
```

### Scenario 2: Existing Database (Pre-Migration)

**Problem**: Database has SSFF files with old track names

**Status**: ⚠️ **Requires migration**

**Solution Options**:

#### Option A: Database Migration Utility
```r
# Proposed: migrate_track_names() function
library(superassp)

migrate_track_names(
  emuDB_path = "/path/to/database",
  mapping_csv = "TRACK_NAMES_MAPPING.csv",
  backup = TRUE
)

# Process:
# 1. Scan all SSFF files
# 2. Read current track names
# 3. Map to new names via CSV
# 4. Update SSFF headers
# 5. Update database configuration
```

#### Option B: Query Aliasing
```r
# emuR could support old names via aliases
# Requires emuR package update
get_trackdata(dbHandle, "f0", "Phonetic")  # Works
get_trackdata(dbHandle, "f0[Hz]", "Phonetic")  # Also works
```

#### Option C: Keep Both Names
```r
# During transition, SSFF files could have duplicate tracks
# Track 1: f0 (old)
# Track 2: f0[Hz] (new, identical data)
# Allows old queries to work while encouraging new names
```

### Scenario 3: Mixed Environment

**Situation**: Some files use old names, some use new names

**Status**: ⚠️ **Needs handling**

**Solution**:
- Use alias system in superassp
- Provide conversion utilities
- Document clearly which version created files

## Proposed Migration Strategy for emuR Users

### Phase 1: Deprecation Period (v0.7.0 - v0.9.0)

**superassp behavior**:
- Outputs use new track names
- Reading supports both old and new names (aliases)
- Deprecation warnings guide users

**emuR impact**:
- New files use new names
- Old databases still work (via aliases)
- Users can query either name

**User action**:
- No immediate action required
- Start planning database migration
- Test queries with new names

### Phase 2: Migration Tools (v0.8.0+)

**superassp provides**:
```r
# New migration utilities

# 1. Check database for track name compatibility
check_emur_database(
  emuDB_path = "/path/to/db",
  mapping_csv = "TRACK_NAMES_MAPPING.csv"
)
# Returns: List of files and tracks needing migration

# 2. Migrate database
migrate_emur_database(
  emuDB_path = "/path/to/db",
  mapping_csv = "TRACK_NAMES_MAPPING.csv",
  backup_dir = "/path/to/backup",
  dry_run = TRUE  # Preview changes
)

# 3. Verify migration
verify_emur_database(
  emuDB_path = "/path/to/db",
  mapping_csv = "TRACK_NAMES_MAPPING.csv"
)
```

**User action**:
- Run check tool on databases
- Create backups
- Run migration (dry-run first)
- Verify results
- Update query scripts

### Phase 3: Clean Break (v1.0.0)

**superassp behavior**:
- Only new track names supported
- No aliases (or limited to common cases)

**emuR impact**:
- Old databases must be migrated
- Queries must use new names

**User action**:
- Complete migration
- Update all query scripts
- Document database version

## Technical Implementation Details

### SSFF Header Structure

SSFF files store track metadata in ASCII header:

```
SSFF -- (9 bytes)
Machine IBM-PC (13 bytes)
Start-Time 0.000000 (23 bytes)
...
Record-Size 4 (15 bytes)
Start-Record 1 (14 bytes)
Column 1 f0 (10 bytes + track name)
Column 2 voicing (15 bytes + track name)
```

**Migration process**:
1. Read binary SSFF file
2. Parse ASCII header
3. Locate "Column" entries
4. Update track name strings
5. Adjust header size if needed
6. Write updated file

**Complexity**: Moderate - requires binary file manipulation

### AsspDataObj Attribute Update

```r
# Current
attr(obj, "tracks") <- c("f0", "voicing")

# After migration
attr(obj, "tracks") <- c("f0[Hz]", "voicing")

# Alias system (Phase 2)
names(obj) <- c("f0[Hz]", "voicing")
attr(obj, ".track_aliases") <- list(
  "f0" = "f0[Hz]",
  "F0" = "f0[Hz]",
  "fo" = "f0[Hz]"
)
```

**Complexity**: Low - straightforward attribute update

## Coordination with emuR Maintainers

### Recommended Communication

1. **Early notification**
   - Email emuR maintainers about planned changes
   - Share migration plan document
   - Request feedback on compatibility

2. **Collaboration opportunities**
   - emuR could add alias support
   - Joint development of migration tools
   - Coordinated release timeline

3. **Documentation**
   - Update emuR tutorials referencing superassp
   - Cross-reference migration guides
   - Provide example workflows

### Contact Information

**emuR Maintainers**:
- Raphael Winkelmann (IPS, LMU Munich)
- GitHub: https://github.com/IPS-LMU/emuR
- Issues: https://github.com/IPS-LMU/emuR/issues

## Risk Assessment

### High Risk

None identified - migration is backwards compatible

### Medium Risk

1. **User confusion during transition**
   - **Mitigation**: Clear documentation, examples, deprecation warnings
   - **Impact**: Temporary, resolved with updated docs

2. **Database migration complexity**
   - **Mitigation**: Automated tools, dry-run mode, backups
   - **Impact**: One-time effort per database

### Low Risk

1. **SSFF format edge cases**
   - **Mitigation**: Extensive testing, fallback options
   - **Impact**: Rare, fixable

2. **Performance overhead (aliases)**
   - **Mitigation**: O(1) hash lookups, only during deprecation
   - **Impact**: Negligible

## Testing Plan

### Test Cases

1. **New database creation**
   - Create SSFF files with new names
   - Import to emuR
   - Query tracks
   - Verify data integrity

2. **Old database migration**
   - Use test database with old names
   - Run migration tool
   - Verify SSFF headers updated
   - Test queries with new names

3. **Mixed databases**
   - Some files old, some new
   - Verify alias system works
   - Test query consistency

4. **Round-trip compatibility**
   - superassp → SSFF → emuR → export → superassp
   - Verify no data loss
   - Check track name preservation

### Test Databases

**Proposed**:
- Create minimal test database (~10 files)
- Include common track types (f0, formants, voice quality)
- Use for automated testing
- Share with emuR maintainers

## Implementation Checklist

### Preparation
- [ ] Contact emuR maintainers
- [ ] Share migration plan
- [ ] Request feedback
- [ ] Identify test databases

### Development
- [ ] Implement `check_emur_database()`
- [ ] Implement `migrate_emur_database()`
- [ ] Implement `verify_emur_database()`
- [ ] Add SSFF header manipulation functions
- [ ] Create migration examples

### Testing
- [ ] Test new database creation
- [ ] Test old database migration
- [ ] Test mixed scenarios
- [ ] Test round-trip compatibility
- [ ] Performance benchmarks

### Documentation
- [ ] Write migration guide for emuR users
- [ ] Update README with emuR section
- [ ] Create tutorial vignette
- [ ] Add FAQ section

### Release
- [ ] Coordinate with emuR release cycle
- [ ] Announce on R-sig-phonetics mailing list
- [ ] Blog post explaining changes
- [ ] Update emuR documentation (if needed)

## Summary

**Overall Assessment**: ✅ **Migration is compatible with emuR**

**Key Points**:
1. SSFF format is flexible - supports any track names
2. No structural breaking changes
3. Alias system provides backwards compatibility
4. Migration tools will ease transition
5. Coordination with emuR maintainers is important

**Timeline**:
- **v0.7.0**: Start deprecation, announce to emuR community
- **v0.8.0**: Release migration tools
- **v0.9.0**: Finalize based on feedback
- **v1.0.0**: Complete migration

**Recommended Next Steps**:
1. Contact emuR maintainers
2. Implement migration utilities
3. Test with real emuR databases
4. Document thoroughly
5. Coordinate release

## References

- emuR Documentation: https://ips-lmu.github.io/The-EMU-SDMS-Manual/
- SSFF Format Specification: Part of wrassp documentation
- emuR GitHub: https://github.com/IPS-LMU/emuR
- superassp GitHub: https://github.com/humlab-speech/superassp
