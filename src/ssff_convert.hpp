/*
 * ssff_convert.hpp -- conversion of SSFF record buffers into R matrices.
 *
 * The SSFF data section is a fixed-layout array of records; every record holds
 * the values of all descriptors back to back. Turning that into one R matrix per
 * descriptor is therefore pure offset arithmetic plus a typed load, a (possible)
 * byte swap, a widening cast and an optional "0 means missing" mask -- no parser,
 * no intermediate buffers. This module performs that conversion in a single
 * pass over the data: the record buffer is walked once, in cache-sized blocks,
 * and every requested descriptor is written straight into the R matrix it
 * belongs to.
 *
 * Design notes (see planning/2026-09-19-ssff-read-performance.md):
 *   - no per-record temporary copies; the destination is written in place
 *   - one specialised loop per (source format x destination type), the format
 *     switch is hoisted out of the record loop
 *   - endianness swapping and the zero->NA mask are fused into the same loop
 *   - the loops are plain C++ so the compiler's auto-vectoriser emits the SIMD
 *     code (measured to match/beat hand-written intrinsics for these kernels)
 */
#ifndef SUPERASSP_SSFF_CONVERT_HPP
#define SUPERASSP_SSFF_CONVERT_HPP

#include <R.h>
#include <Rinternals.h>
#include <stddef.h>

#include <dataobj.h>

#ifdef __cplusplus
extern "C" {
#endif

/* One descriptor to materialise. */
typedef struct ssff_track_s {
    size_t  offset;       /* byte offset of the descriptor inside a record */
    size_t  numFields;    /* values per record */
    int     format;       /* DF_* storage format of the descriptor */
    int     destType;     /* REALSXP (float formats) or INTSXP (integer formats) */
    void   *dest;         /* REAL() or INTEGER() of the destination matrix */
    int     zeroToNa;     /* map values that are exactly 0 to NA for this track */
} ssff_track_t;

/*
 * Fill the destinations of all tracks from an interleaved record buffer.
 *
 * records      interleaved data, numRecords * recordSize bytes, in *file* byte
 *              order (swapping is done here, per value, when swapBytes is set)
 * numRecords   number of valid records in the buffer
 * recordSize   bytes per record
 * numTracks    number of entries in 'tracks'
 * swapBytes    non-zero if the file byte order differs from the host byte order
 * numThreads   >1 splits the record blocks over OpenMP threads when available,
 *              always 1 without OpenMP support
 *
 * Returns 0 on success, -1 when a descriptor has a format that cannot be
 * converted (the caller reports that as an error, as before).
 */
int ssff_convert_records(const void *records, long numRecords, size_t recordSize,
                         int numTracks, const ssff_track_t *tracks,
                         int swapBytes, int numThreads);

/* Byte size of a DF_* storage format, 0 when unsupported. */
size_t ssff_format_size(int format);

/* Non-zero when a DF_* format is materialised as REALSXP. */
int ssff_format_is_float(int format);

#ifdef __cplusplus
}
#endif

#endif /* SUPERASSP_SSFF_CONVERT_HPP */
