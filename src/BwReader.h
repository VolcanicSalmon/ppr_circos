#pragma once
// Minimal RAII C++ wrapper around libBigWig for the openFrameworks circos app.
// Reads BINNED coverage (one bwStats call per chrom/region) — exactly what a ring
// track needs. Values are returned as |value| with NaN->0 (strand-agnostic depth).
//
// libBigWig here is built with -DNOCURL, so LOCAL .bw files only (no http/https).
// Link: -I libBigWig  -L libBigWig -lBigWig -lz   (see config.make)
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>

extern "C" {
#include "bigWig.h"
}

class BwReader {
public:
    struct Chrom { std::string name; uint32_t length; };

    // Scoped mirror of libBigWig's bwStatsType — its unscoped `max`/`min`
    // constants collide with oF / std::max once ofMain.h is also included.
    enum class Stat { Mean = 0, Stdev = 1, Max = 2, Min = 3, Coverage = 4, Sum = 5 };

    BwReader() : fp(nullptr) {}
    explicit BwReader(const std::string& path) : fp(nullptr) { open(path); }
    ~BwReader() { close(); }

    BwReader(const BwReader&) = delete;
    BwReader& operator=(const BwReader&) = delete;

    bool open(const std::string& path) {
        ensureInit();
        close();
        // bwOpen takes non-const char*; callback NULL is fine for local files.
        fp = bwOpen(const_cast<char*>(path.c_str()), nullptr, const_cast<char*>("r"));
        return fp != nullptr;
    }
    void close() { if (fp) { bwClose(fp); fp = nullptr; } }
    bool ok() const { return fp != nullptr; }

    std::vector<Chrom> chroms() const {
        std::vector<Chrom> out;
        if (!fp || !fp->cl) return out;
        out.reserve(fp->cl->nKeys);
        for (int64_t i = 0; i < fp->cl->nKeys; ++i)
            out.push_back({ std::string(fp->cl->chrom[i]), fp->cl->len[i] });
        return out;
    }

    int64_t chromLen(const std::string& name) const {
        if (!fp || !fp->cl) return -1;
        for (int64_t i = 0; i < fp->cl->nKeys; ++i)
            if (name == fp->cl->chrom[i]) return (int64_t)fp->cl->len[i];
        return -1;
    }

    // Binned stat over [start,end) into nBins; returns |value| with NaN->0.
    std::vector<float> binned(const std::string& chrom, uint32_t start, uint32_t end,
                              uint32_t nBins, Stat type = Stat::Max) const {
        std::vector<float> out(nBins, 0.f);
        if (!fp || nBins == 0 || end <= start) return out;
        double* v = bwStats(fp, const_cast<char*>(chrom.c_str()),
                            start, end, nBins, (enum bwStatsType)(int)type);
        if (!v) return out;
        for (uint32_t i = 0; i < nBins; ++i)
            out[i] = std::isnan(v[i]) ? 0.f : (float)std::fabs(v[i]);
        free(v);
        return out;
    }

    // Whole-chromosome convenience.
    std::vector<float> binnedChrom(const std::string& chrom, uint32_t nBins,
                                   Stat type = Stat::Max) const {
        int64_t L = chromLen(chrom);
        if (L <= 0) return std::vector<float>(nBins, 0.f);
        return binned(chrom, 0, (uint32_t)L, nBins, type);
    }

private:
    bigWigFile_t* fp;
    static void ensureInit() {
        static bool inited = false;
        if (!inited) { bwInit(1 << 17); inited = true; }   // 128 KiB read buffer
    }
};
