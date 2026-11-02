// shape.h -- legal soft-block shapes, derived exactly from the verifier rules
// (docs/design.md §6 gate checks 5 and 6):
//
//   area tolerance:  100*|w*h - A| <= A   <=>  ceil(99A/100) <= w*h <= floor(101A/100)
//   aspect ratio:    R*h >= w  and  R*w >= h
//
// For a fixed width w the legal heights form the integer interval
//   hmin(w) = max(ceil(whmin/w), ceil(w/R)),  hmax(w) = min(floor(whmax/w), R*w),
// and w is legal iff hmin(w) <= hmax(w). The legal widths are found by a scan
// (no closed-form sqrt endpoints) and stored as an explicit list, so the solver
// never depends on the list being one contiguous interval; test_shape.cpp
// checks against brute force that it is one for A in [1e4, 1e6].
#ifndef FP_SHAPE_H
#define FP_SHAPE_H

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace fp {

inline int64_t ceil_div(int64_t a, int64_t b) { return (a + b - 1) / b; }  // a >= 0, b > 0

struct ShapeRange {
    int64_t A = 0, R = 1;
    int64_t whmin = 0, whmax = 0;   // legal w*h range
    std::vector<int> widths;        // all legal widths, ascending
    std::vector<int> heights;       // heights[i] = hmin(widths[i])
};

inline int64_t shape_hmin(int64_t whmin, int64_t R, int64_t w) {
    return std::max(ceil_div(whmin, w), ceil_div(w, R));
}
inline int64_t shape_hmax(int64_t whmax, int64_t R, int64_t w) {
    return std::min(whmax / w, R * w);
}

inline ShapeRange make_shape_range(int64_t A, int64_t R) {
    ShapeRange s;
    s.A = A;
    s.R = R;
    s.whmin = ceil_div(99 * A, 100);
    s.whmax = (101 * A) / 100;
    // A legal w has h >= w/R, so w*w <= R*w*h <= R*whmax.
    int64_t wlim = (int64_t)std::sqrt((double)(R * s.whmax)) + 2;
    for (int64_t w = 1; w <= wlim; ++w) {
        int64_t lo = shape_hmin(s.whmin, R, w), hi = shape_hmax(s.whmax, R, w);
        if (lo <= hi) {
            s.widths.push_back((int)w);
            s.heights.push_back((int)lo);
        }
    }
    return s;
}

}  // namespace fp

#endif  // FP_SHAPE_H
