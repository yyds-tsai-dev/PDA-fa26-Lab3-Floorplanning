// test_shape.cpp -- checks shape.h against a brute-force evaluation of the
// verifier's soft-block rules (design.md §6 checks 5 and 6).
//
// For every tested area A: for every w in [1, wcap] and every h that could
// possibly satisfy the area tolerance, evaluate the verifier predicates
// directly. Asserts that
//   (a) the set of legal widths equals shape.h's list,
//   (b) that list is one contiguous interval [wmin, wmax] (no holes),
//   (c) heights[i] is the smallest legal h for widths[i].
// Tested areas: every A in [10000, 12000], every A in [999000, 1000000],
// and 3000 pseudo-random A in [1e4, 1e6]; R = 3 (released cases), plus a
// smaller sample for R = 1, 2, 4.
#include <cstdint>
#include <cstdio>
#include <vector>

#include "shape.h"

static bool verifier_ok(int64_t A, int64_t R, int64_t w, int64_t h) {
    if (w <= 0 || h <= 0) return false;
    int64_t d = w * h - A;
    if (d < 0) d = -d;
    if (100 * d > A) return false;
    if (R * h < w || R * w < h) return false;
    return true;
}

static int check(int64_t A, int64_t R) {
    fp::ShapeRange s = fp::make_shape_range(A, R);
    // Any legal w satisfies w*w <= R*w*h <= R*1.01A; cap generously above that.
    int64_t wcap = 2 * (int64_t)std::sqrt((double)(R * A)) + 10;
    std::vector<int> legal;
    std::vector<int> hmin;
    for (int64_t w = 1; w <= wcap; ++w) {
        // |w*h - A| <= A/100 needs h in [0.99A/w, 1.01A/w]; scan a superset.
        int64_t hlo = (98 * A) / (100 * w) - 1, hhi = (102 * A) / (100 * w) + 2;
        if (hlo < 1) hlo = 1;
        int64_t best = -1;
        for (int64_t h = hlo; h <= hhi; ++h)
            if (verifier_ok(A, R, w, h)) { best = h; break; }
        if (best > 0) { legal.push_back((int)w); hmin.push_back((int)best); }
    }
    if (legal != s.widths) {
        std::printf("FAIL A=%lld R=%lld: legal width set differs (brute %zu, shape.h %zu)\n",
                    (long long)A, (long long)R, legal.size(), s.widths.size());
        return 1;
    }
    if (legal.empty()) {
        std::printf("FAIL A=%lld R=%lld: no legal width\n", (long long)A, (long long)R);
        return 1;
    }
    if (legal.back() - legal.front() + 1 != (int)legal.size()) {
        std::printf("FAIL A=%lld R=%lld: legal widths are not contiguous\n", (long long)A, (long long)R);
        return 1;
    }
    if (hmin != s.heights) {
        std::printf("FAIL A=%lld R=%lld: hmin differs\n", (long long)A, (long long)R);
        return 1;
    }
    return 0;
}

int main() {
    int fails = 0, tested = 0;
    size_t minw = 1 << 30, maxw = 0;
    auto run = [&](int64_t A, int64_t R) {
        fails += check(A, R);
        ++tested;
        if (R == 3) {
            size_t k = fp::make_shape_range(A, R).widths.size();
            minw = std::min(minw, k);
            maxw = std::max(maxw, k);
        }
    };
    for (int64_t A = 10000; A <= 12000; ++A) run(A, 3);
    for (int64_t A = 999000; A <= 1000000; ++A) run(A, 3);
    uint64_t x = 88172645463325252ULL;
    for (int i = 0; i < 3000; ++i) {
        x ^= x << 13; x ^= x >> 7; x ^= x << 17;
        run(10000 + (int64_t)(x % 990001), 3);
    }
    for (int64_t R = 1; R <= 4; ++R) {
        if (R == 3) continue;
        for (int i = 0; i < 300; ++i) {
            x ^= x << 13; x ^= x >> 7; x ^= x << 17;
            run(10000 + (int64_t)(x % 990001), R);
        }
    }
    std::printf("test_shape: %d areas tested, %d failures; R=3 legal widths per block: %zu..%zu\n",
                tested, fails, minw, maxw);
    return fails ? 1 : 0;
}
