/*
 * ssff_convert.cpp -- see ssff_convert.hpp.
 *
 * The record buffer is walked once, in blocks of records. Inside a block the
 * loop nest is descriptor -> field -> record, which is the order measured to be
 * fastest (planning/2026-09-19-ssff-read-performance.md, section 2.5): the
 * source block stays hot while every destination column is written with
 * unit-stride stores.
 */
#include "ssff_convert.hpp"

#include <stdint.h>
#include <string.h>

namespace {

/* ------------------------------------------------------------------ */
/* loads: alignment safe, byte order aware                             */
/* ------------------------------------------------------------------ */

inline uint16_t bswap16(uint16_t v) { return (uint16_t)((v >> 8) | (v << 8)); }

/* Primary templates; every supported source type below has a definition. */
template <typename T> inline T loadRaw(const char *p);
template <typename T> inline T loadSwapped(const char *p);

template <>
inline float loadRaw<float>(const char *p) { float v; memcpy(&v, p, sizeof(v)); return v; }
template <>
inline double loadRaw<double>(const char *p) { double v; memcpy(&v, p, sizeof(v)); return v; }
template <>
inline uint8_t loadRaw<uint8_t>(const char *p) { return (uint8_t)p[0]; }
template <>
inline int8_t loadRaw<int8_t>(const char *p) { return (int8_t)p[0]; }
template <>
inline uint16_t loadRaw<uint16_t>(const char *p) { uint16_t v; memcpy(&v, p, sizeof(v)); return v; }
template <>
inline int16_t loadRaw<int16_t>(const char *p) { int16_t v; memcpy(&v, p, sizeof(v)); return v; }
template <>
inline uint32_t loadRaw<uint32_t>(const char *p) { uint32_t v; memcpy(&v, p, sizeof(v)); return v; }
template <>
inline int32_t loadRaw<int32_t>(const char *p) { int32_t v; memcpy(&v, p, sizeof(v)); return v; }

template <>
inline uint16_t loadSwapped<uint16_t>(const char *p) { uint16_t v; memcpy(&v, p, sizeof(v)); return bswap16(v); }
template <>
inline int16_t loadSwapped<int16_t>(const char *p) { int16_t v; memcpy(&v, p, sizeof(v)); return (int16_t)bswap16((uint16_t)v); }
template <>
inline uint32_t loadSwapped<uint32_t>(const char *p) {
    uint32_t v; memcpy(&v, p, sizeof(v));
    return __builtin_bswap32(v);
}
template <>
inline int32_t loadSwapped<int32_t>(const char *p) {
    int32_t v; memcpy(&v, p, sizeof(v));
    return (int32_t)__builtin_bswap32((uint32_t)v);
}
template <>
inline float loadSwapped<float>(const char *p) {
    uint32_t v; memcpy(&v, p, sizeof(v));
    v = __builtin_bswap32(v);
    float f; memcpy(&f, &v, sizeof(f));
    return f;
}
template <>
inline double loadSwapped<double>(const char *p) {
    uint64_t v; memcpy(&v, p, sizeof(v));
    v = __builtin_bswap64(v);
    double d; memcpy(&d, &v, sizeof(d));
    return d;
}
template <>
inline uint8_t loadSwapped<uint8_t>(const char *p) { return (uint8_t)p[0]; }
template <>
inline int8_t loadSwapped<int8_t>(const char *p) { return (int8_t)p[0]; }

/*
 * One conversion shape for one (source type, destination type) pair.
 * SWAP and NA are template parameters so that neither branch survives in the
 * inner loop.
 */
template <bool SWAP, bool NA, typename TSrc, typename TDst>
struct Loop {
    static inline TDst convert(const char *p, TDst naValue) {
        const TSrc raw = SWAP ? loadSwapped<TSrc>(p) : loadRaw<TSrc>(p);
        TDst v = (TDst)raw;
        if (NA && raw == (TSrc)0)
            v = naValue;
        return v;
    }

    /* Records are contiguous for this field (single-value, single-track file). */
    static inline void contiguous(const char *src, TDst *dst, long n, TDst naValue) {
        for (long i = 0; i < n; ++i)
            dst[i] = convert(src + i * sizeof(TSrc), naValue);
    }

    /* One field across consecutive records: source stride = recordSize. */
    static inline void strided(const char *src, size_t stride, TDst *dst, long n, TDst naValue) {
        for (long i = 0; i < n; ++i)
            dst[i] = convert(src + i * stride, naValue);
    }
};

template <typename TSrc, typename TDst>
struct Types {
    static inline void contiguous(bool swap, bool na, const char *src, TDst *dst, long n, TDst naValue) {
        if (swap) {
            if (na) Loop<true, true, TSrc, TDst>::contiguous(src, dst, n, naValue);
            else    Loop<true, false, TSrc, TDst>::contiguous(src, dst, n, naValue);
        } else {
            if (na) Loop<false, true, TSrc, TDst>::contiguous(src, dst, n, naValue);
            else    Loop<false, false, TSrc, TDst>::contiguous(src, dst, n, naValue);
        }
    }

    static inline void strided(bool swap, bool na, const char *src, size_t stride, TDst *dst, long n, TDst naValue) {
        if (swap) {
            if (na) Loop<true, true, TSrc, TDst>::strided(src, stride, dst, n, naValue);
            else    Loop<true, false, TSrc, TDst>::strided(src, stride, dst, n, naValue);
        } else {
            if (na) Loop<false, true, TSrc, TDst>::strided(src, stride, dst, n, naValue);
            else    Loop<false, false, TSrc, TDst>::strided(src, stride, dst, n, naValue);
        }
    }
};

/* Per (track, field, block) dispatch: one branch, never inside the loops. */
inline void fillField(int format, int destType, bool contiguous, bool swap, bool na,
                      const char *src, size_t stride, void *dst, long n)
{
    if (destType == INTSXP) {
        int *d = (int *)dst;
        const int naValue = NA_INTEGER;
        switch (format) {
        case DF_UINT8:  contiguous ? Types<uint8_t, int>::contiguous(swap, na, src, d, n, naValue)
                                   : Types<uint8_t, int>::strided(swap, na, src, stride, d, n, naValue); break;
        case DF_INT8:   contiguous ? Types<int8_t, int>::contiguous(swap, na, src, d, n, naValue)
                                   : Types<int8_t, int>::strided(swap, na, src, stride, d, n, naValue); break;
        case DF_UINT16: contiguous ? Types<uint16_t, int>::contiguous(swap, na, src, d, n, naValue)
                                   : Types<uint16_t, int>::strided(swap, na, src, stride, d, n, naValue); break;
        case DF_INT16:  contiguous ? Types<int16_t, int>::contiguous(swap, na, src, d, n, naValue)
                                   : Types<int16_t, int>::strided(swap, na, src, stride, d, n, naValue); break;
        case DF_UINT32: contiguous ? Types<uint32_t, int>::contiguous(swap, na, src, d, n, naValue)
                                   : Types<uint32_t, int>::strided(swap, na, src, stride, d, n, naValue); break;
        case DF_INT32:  contiguous ? Types<int32_t, int>::contiguous(swap, na, src, d, n, naValue)
                                   : Types<int32_t, int>::strided(swap, na, src, stride, d, n, naValue); break;
        default: break;                        /* rejected up front */
        }
    } else {
        double *d = (double *)dst;
        const double naValue = NA_REAL;
        switch (format) {
        case DF_REAL32: contiguous ? Types<float, double>::contiguous(swap, na, src, d, n, naValue)
                                   : Types<float, double>::strided(swap, na, src, stride, d, n, naValue); break;
        case DF_REAL64: contiguous ? Types<double, double>::contiguous(swap, na, src, d, n, naValue)
                                   : Types<double, double>::strided(swap, na, src, stride, d, n, naValue); break;
        default: break;                        /* rejected up front */
        }
    }
}

/*
 * Record-major variant for descriptors with several fields (spectra): the
 * loads stay contiguous inside each record while every destination column is
 * still written sequentially in runs of 'count' entries. Iterating field-major
 * instead would turn every load into a strided one, which measured slower
 * whenever a record holds more than a handful of values.
 */
template <typename TSrc, typename TDst>
struct ByRecord {
    template <bool SWAP, bool NA>
    static inline void loop(const char *records, size_t recordSize, size_t offset,
                            size_t numFields, TDst *dst, long mb, long count,
                            long numRecords, TDst naValue)
    {
        TDst *base = dst + mb;
        for (long m = 0; m < count; ++m) {
            const char *rec = records + (size_t)m * recordSize + offset;
            TDst *out = base + m;
            for (size_t n = 0; n < numFields; ++n)
                out[n * (size_t)numRecords] = Loop<SWAP, NA, TSrc, TDst>::convert(rec + n * sizeof(TSrc), naValue);
        }
    }

    static inline void run(bool swap, bool na, const char *records, size_t recordSize,
                           size_t offset, size_t numFields, TDst *dst, long mb, long count,
                           long numRecords, TDst naValue)
    {
        if (swap) {
            if (na) loop<true, true>(records, recordSize, offset, numFields, dst, mb, count, numRecords, naValue);
            else    loop<true, false>(records, recordSize, offset, numFields, dst, mb, count, numRecords, naValue);
        } else {
            if (na) loop<false, true>(records, recordSize, offset, numFields, dst, mb, count, numRecords, naValue);
            else    loop<false, false>(records, recordSize, offset, numFields, dst, mb, count, numRecords, naValue);
        }
    }
};

/* Same, dispatched on the storage format (called once per track and block). */
inline void fillTrackRecordMajor(int format, int destType, bool swap, bool na,
                                 const char *records, size_t recordSize, size_t offset,
                                 size_t numFields, void *dst, long mb, long count,
                                 long numRecords)
{
    if (destType == INTSXP) {
        int *d = (int *)dst;
        const int naValue = NA_INTEGER;
        switch (format) {
        case DF_UINT8:  ByRecord<uint8_t, int>::run(swap, na, records, recordSize, offset, numFields, d, mb, count, numRecords, naValue); break;
        case DF_INT8:   ByRecord<int8_t, int>::run(swap, na, records, recordSize, offset, numFields, d, mb, count, numRecords, naValue); break;
        case DF_UINT16: ByRecord<uint16_t, int>::run(swap, na, records, recordSize, offset, numFields, d, mb, count, numRecords, naValue); break;
        case DF_INT16:  ByRecord<int16_t, int>::run(swap, na, records, recordSize, offset, numFields, d, mb, count, numRecords, naValue); break;
        case DF_UINT32: ByRecord<uint32_t, int>::run(swap, na, records, recordSize, offset, numFields, d, mb, count, numRecords, naValue); break;
        case DF_INT32:  ByRecord<int32_t, int>::run(swap, na, records, recordSize, offset, numFields, d, mb, count, numRecords, naValue); break;
        default: break;
        }
    } else {
        double *d = (double *)dst;
        const double naValue = NA_REAL;
        (void)naValue;
        switch (format) {
        case DF_REAL32: ByRecord<float, double>::run(swap, na, records, recordSize, offset, numFields, d, mb, count, numRecords, naValue); break;
        case DF_REAL64: ByRecord<double, double>::run(swap, na, records, recordSize, offset, numFields, d, mb, count, numRecords, naValue); break;
        default: break;
        }
    }
}

/* Target ~2 MB of source per block, clamped to keep blocks and thread splits sane. */
inline long blockRecords(size_t recordSize) {
    if (recordSize == 0)
        return 1024;
    long rb = (long)((2u * 1024u * 1024u) / recordSize);
    if (rb < 16)
        rb = 16;
    if (rb > 4096)
        rb = 4096;
    return rb;
}

}  /* namespace */

size_t ssff_format_size(int format)
{
    switch (format) {
    case DF_UINT8:
    case DF_INT8:
        return 1;
    case DF_UINT16:
    case DF_INT16:
        return 2;
    case DF_UINT32:
    case DF_INT32:
    case DF_REAL32:
        return 4;
    case DF_REAL64:
        return 8;
    default:
        return 0;
    }
}

int ssff_format_is_float(int format)
{
    return (format == DF_REAL32 || format == DF_REAL64);
}

extern "C" int ssff_convert_records(const void *records, long numRecords, size_t recordSize,
                                    int numTracks, const ssff_track_t *tracks,
                                    int swapBytes, int numThreads)
{
    if (numRecords < 1 || numTracks < 1 || records == NULL || tracks == NULL)
        return 0;

    /* Reject unsupported formats before anything is written. */
    for (int t = 0; t < numTracks; ++t) {
        if (ssff_format_size(tracks[t].format) == 0)
            return -1;
        if (tracks[t].dest == NULL)
            return -1;
    }

    const unsigned char *base = (const unsigned char *)records;
    const long rb      = blockRecords(recordSize);
    const long nblocks = (numRecords + rb - 1) / rb;
    const bool swap    = swapBytes != 0;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(numThreads) if (numThreads > 1)
#endif
    for (long b = 0; b < nblocks; ++b) {
        const long mb = b * rb;
        const long hi = (mb + rb < numRecords) ? (mb + rb) : numRecords;
        const long count = hi - mb;
        const char *blockBase = (const char *)(base + (size_t)mb * recordSize);

        for (int t = 0; t < numTracks; ++t) {
            const ssff_track_t *tr = &tracks[t];
            const size_t elem = ssff_format_size(tr->format);
            const char *trackBase = blockBase + tr->offset;
            const bool contiguous = (recordSize == elem && tr->offset == 0);

            if (tr->numFields > 1 && !contiguous) {
                fillTrackRecordMajor(tr->format, tr->destType, swap, tr->zeroToNa != 0,
                                     blockBase, recordSize, tr->offset, tr->numFields,
                                     tr->dest, mb, count, numRecords);
                continue;
            }

            for (size_t n = 0; n < tr->numFields; ++n) {
                const char *src = trackBase + n * elem;
                void *dst;
                if (tr->destType == INTSXP)
                    dst = (void *)((int *)tr->dest + (size_t)mb + n * (size_t)numRecords);
                else
                    dst = (void *)((double *)tr->dest + (size_t)mb + n * (size_t)numRecords);
                fillField(tr->format, tr->destType, contiguous, swap, tr->zeroToNa != 0,
                          src, recordSize, dst, count);
            }
        }
    }
    return 0;
}
