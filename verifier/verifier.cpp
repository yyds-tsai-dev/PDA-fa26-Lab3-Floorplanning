// Lab3 verifier -- implements docs/design.md §6.
//
//   ./verifier <input> <output>
//
// Exit codes: 0 feasible output, 1 infeasible output (any gate failure,
// including an unreadable output file), 2 usage error or malformed input.
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "lab3_io.h"

using lab3::Instance;
using lab3::Kind;
using lab3::Rect;

namespace {

const int kNumChecks = 8;
const char *kCheckNames[kNumChecks] = {
    "Output format",
    "Block names",
    "Positive size",
    "Fixed/preplaced shape",
    "Soft area tolerance",
    "Soft aspect ratio",
    "Inside outline",
    "Overlap-free",
};

std::string itos(int64_t v) { return std::to_string(v); }
std::string itos128(__int128 v) {
    if (v == 0) return "0";
    bool neg = v < 0;
    unsigned __int128 a = neg ? (unsigned __int128)(-(v + 1)) + 1 : (unsigned __int128)v;
    std::string s;
    while (a) { s.insert(s.begin(), char('0' + (int)(a % 10))); a /= 10; }
    return neg ? "-" + s : s;
}

std::vector<std::string> split_ws(const std::string &line) {
    std::vector<std::string> t;
    std::istringstream ss(line);
    std::string s;
    while (ss >> s) t.push_back(s);
    return t;
}

// Reported wirelength: [+-]?[0-9]+(\.[0-9]+)? as an exact decimal
// mant / 10^scale. Fraction digits beyond 18 are truncated.
struct Decimal {
    __int128 mant = 0;
    int scale = 0;
};

bool parse_decimal(const std::string &s, Decimal &d) {
    size_t i = 0;
    bool neg = false;
    if (i < s.size() && (s[i] == '+' || s[i] == '-')) { neg = s[i] == '-'; ++i; }
    size_t start = i;
    __int128 m = 0;
    while (i < s.size() && s[i] >= '0' && s[i] <= '9') {
        m = m * 10 + (s[i] - '0');
        if (m > (__int128)1e19) return false;  // absurdly large
        ++i;
    }
    if (i == start) return false;
    int scale = 0;
    if (i < s.size() && s[i] == '.') {
        ++i;
        size_t fs = i;
        while (i < s.size() && s[i] >= '0' && s[i] <= '9') {
            if (scale < 18) { m = m * 10 + (s[i] - '0'); ++scale; }
            ++i;
        }
        if (i == fs) return false;
    }
    if (i != s.size()) return false;
    d.mant = neg ? -m : m;
    d.scale = scale;
    return true;
}

struct OutLine {
    std::string name;
    Rect r;
};

// Gate state.
struct Gate {
    int failed = -1;           // index of first failing check, -1 if none
    std::string msg;
    bool fail(int check, const std::string &m) {
        failed = check;
        msg = m;
        return false;
    }
};

std::string rect_str(const Rect &r) {
    return itos(r.x) + " " + itos(r.y) + " " + itos(r.w) + " " + itos(r.h);
}

// Runs gate checks 1..8 in order (design.md §6). On success fills rects
// (indexed by input block order) and the reported WL.
bool run_gate(const Instance &inst, const std::string &out_path, Gate &g,
              std::vector<Rect> &rects, Decimal &reported, std::string &reported_str) {
    const size_t n = inst.blocks.size();

    // ---- 1. Output parses.
    std::ifstream f(out_path);
    if (!f) return g.fail(0, "Parse error: cannot open output file '" + out_path + "'");
    std::vector<std::pair<int, std::vector<std::string>>> lines;  // (line no, tokens), blank lines skipped
    {
        std::string line;
        int ln = 0;
        while (std::getline(f, line)) {
            ++ln;
            auto t = split_ws(line);
            if (!t.empty()) lines.emplace_back(ln, t);
        }
    }
    auto at = [&](size_t i) { return "Parse error (line " + std::to_string(lines[i].first) + "): "; };
    if (lines.empty()) return g.fail(0, "Parse error: empty output file");
    {
        const auto &t = lines[0].second;
        if (t[0] != "WIRELENGTH") return g.fail(0, at(0) + "expected 'WIRELENGTH <wl>', got '" + t[0] + "'");
        if (t.size() != 2) return g.fail(0, at(0) + "expected 'WIRELENGTH <wl>' (2 tokens), got " + std::to_string(t.size()) + " tokens");
        if (!parse_decimal(t[1], reported)) return g.fail(0, at(0) + "invalid wirelength '" + t[1] + "'");
        reported_str = t[1];
    }
    if (lines.size() < 2) return g.fail(0, "Parse error: missing 'BLOCKS <n>' line");
    int64_t declared;
    {
        const auto &t = lines[1].second;
        if (t[0] != "BLOCKS") return g.fail(0, at(1) + "expected 'BLOCKS <n>', got '" + t[0] + "'");
        if (t.size() != 2) return g.fail(0, at(1) + "expected 'BLOCKS <n>' (2 tokens), got " + std::to_string(t.size()) + " tokens");
        if (!lab3::parse_int64(t[1], declared)) return g.fail(0, at(1) + "invalid block count '" + t[1] + "'");
        if (declared != (int64_t)n) return g.fail(0, "Block count: " + t[1] + " (expected " + itos((int64_t)n) + ")");
    }
    std::vector<OutLine> outs;
    for (size_t i = 2; i < lines.size(); ++i) {
        const auto &t = lines[i].second;
        if (outs.size() == n) return g.fail(0, at(i) + "extra line after " + itos((int64_t)n) + " block lines");
        if (t.size() != 5) return g.fail(0, at(i) + "expected '<name> <x> <y> <w> <h>' (5 tokens), got " + std::to_string(t.size()) + " tokens");
        OutLine o;
        o.name = t[0];
        int64_t *v[4] = {&o.r.x, &o.r.y, &o.r.w, &o.r.h};
        for (int k = 0; k < 4; ++k)
            if (!lab3::parse_int64(t[k + 1], *v[k])) return g.fail(0, at(i) + "non-integer value '" + t[k + 1] + "' for " + t[0]);
        outs.push_back(o);
    }
    if (outs.size() != n)
        return g.fail(0, "Parse error: found " + itos((int64_t)outs.size()) + " block lines, expected " + itos((int64_t)n));

    // ---- 2. Every input block exactly once, no unknown names.
    rects.assign(n, Rect());
    std::vector<char> seen(n, 0);
    for (const OutLine &o : outs) {
        auto it = inst.block_index.find(o.name);
        if (it == inst.block_index.end()) return g.fail(1, "Unknown name: " + o.name);
        if (seen[it->second]) return g.fail(1, "Duplicate block: " + o.name);
        seen[it->second] = 1;
        rects[it->second] = o.r;
    }
    for (size_t i = 0; i < n; ++i)
        if (!seen[i]) return g.fail(1, "Missing block: " + inst.blocks[i].name);

    // Checks 3..8 iterate in input block order.
    // ---- 3. w > 0, h > 0.
    for (size_t i = 0; i < n; ++i)
        if (rects[i].w <= 0 || rects[i].h <= 0)
            return g.fail(2, "Non-positive size: " + inst.blocks[i].name + " (w " + itos(rects[i].w) + ", h " + itos(rects[i].h) + ")");

    // ---- 4. Fixed-shape / preplaced exactness.
    for (size_t i = 0; i < n; ++i) {
        const auto &b = inst.blocks[i];
        const Rect &r = rects[i];
        if (b.kind == Kind::Fixed && (r.w != b.w || r.h != b.h))
            return g.fail(3, "Fixed-shape mismatch: " + b.name + " (w " + itos(r.w) + ", h " + itos(r.h) +
                                 "; expected w " + itos(b.w) + ", h " + itos(b.h) + ")");
        if (b.kind == Kind::Preplaced && (r.x != b.x || r.y != b.y || r.w != b.w || r.h != b.h)) {
            Rect e{b.x, b.y, b.w, b.h};
            return g.fail(3, "Preplaced mismatch: " + b.name + " (got " + rect_str(r) + ", expected " + rect_str(e) + ")");
        }
    }

    // ---- 5. Soft area tolerance: 100*|w*h - area| <= area.
    for (size_t i = 0; i < n; ++i) {
        const auto &b = inst.blocks[i];
        if (b.kind != Kind::Soft) continue;
        __int128 a = (__int128)rects[i].w * rects[i].h;
        __int128 d = a - b.area;
        if (d < 0) d = -d;
        if (100 * d > (__int128)b.area)
            return g.fail(4, "Area tolerance: " + b.name + " (area " + itos128(a) + ", target " + itos(b.area) + ")");
    }

    // ---- 6. Soft aspect ratio: R*h >= w and R*w >= h.
    for (size_t i = 0; i < n; ++i) {
        const auto &b = inst.blocks[i];
        if (b.kind != Kind::Soft) continue;
        const Rect &r = rects[i];
        if ((__int128)inst.R * r.h < r.w || (__int128)inst.R * r.w < r.h)
            return g.fail(5, "Aspect ratio: " + b.name + " (w " + itos(r.w) + ", h " + itos(r.h) + ")");
    }

    // ---- 7. Inside outline.
    for (size_t i = 0; i < n; ++i) {
        const Rect &r = rects[i];
        if (r.x < 0 || r.y < 0 || (__int128)r.x + r.w > inst.W || (__int128)r.y + r.h > inst.H)
            return g.fail(6, "Out of outline: " + inst.blocks[i].name + " (" + rect_str(r) + ")");
    }

    // ---- 8. Overlap-free (touching allowed). O(n^2); n <= 200.
    for (size_t i = 0; i < n; ++i) {
        const Rect &a = rects[i];
        for (size_t j = i + 1; j < n; ++j) {
            const Rect &b = rects[j];
            bool ox = std::max(a.x, b.x) < std::min(a.x + a.w, b.x + b.w);
            bool oy = std::max(a.y, b.y) < std::min(a.y + a.h, b.y + b.h);
            if (ox && oy) return g.fail(7, "Overlap: " + inst.blocks[i].name + " " + inst.blocks[j].name);
        }
    }
    return true;
}

void print_gate(const Gate &g) {
    std::printf("Gate checks:\n");
    for (int c = 0; c < kNumChecks; ++c) {
        const char *st = (g.failed < 0 || c < g.failed) ? "PASS" : (c == g.failed ? "FAIL" : "SKIP");
        std::printf("  %d. %-22s %s", c + 1, kCheckNames[c], st);
        if (c == g.failed) std::printf("  %s", g.msg.c_str());
        std::printf("\n");
    }
}

}  // namespace

int main(int argc, char **argv) {
    if (argc != 3) {
        std::fprintf(stderr, "Usage: %s <input> <output>\n", argv[0]);
        return 2;
    }
    Instance inst;
    std::string err;
    if (!lab3::parse_input_file(argv[1], inst, err)) {
        std::fprintf(stderr, "Input error: %s\n", err.c_str());
        return 2;
    }

    Gate g;
    std::vector<Rect> rects;
    Decimal reported;
    std::string reported_str;
    bool feasible = run_gate(inst, argv[2], g, rects, reported, reported_str);
    print_gate(g);

    if (!feasible) {
        std::printf("Result: INFEASIBLE (%s)\n", g.msg.c_str());
        std::printf("Total 0.00 / 20\n");
        return 1;
    }

    // ---- Scored constraints.
    const int64_t W2 = lab3::wirelength2(inst, rects);   // 2*WL
    const int64_t T2 = inst.T2;                           // 2*T

    // |reported - WL| <= 0.01, exactly: scale both to 10^s with s >= 2.
    int s = reported.scale < 2 ? 2 : reported.scale;
    __int128 p10 = 1;
    for (int k = 0; k < s; ++k) p10 *= 10;
    __int128 rep_scaled = reported.mant;
    for (int k = reported.scale; k < s; ++k) rep_scaled *= 10;
    __int128 wl_scaled = (__int128)W2 * (p10 / 2);       // WL * 10^s (s >= 1, so p10 even)
    __int128 diff = rep_scaled - wl_scaled;
    if (diff < 0) diff = -diff;
    __int128 tol = p10 / 100;                             // 0.01 * 10^s
    bool wl_match = diff <= tol;

    double wl_score;
    if (!wl_match) wl_score = 0.0;
    else if (W2 <= T2) wl_score = 5.0;
    else if (20 * (__int128)W2 <= 21 * (__int128)T2)      // WL <= 1.05T
        wl_score = 5.0 * (double)(21 * (__int128)T2 - 20 * (__int128)W2) / (double)T2;
    else wl_score = 0.0;

    const size_t p = inst.boundaries.size();
    size_t bsat = 0;
    for (const auto &bd : inst.boundaries) {
        const Rect &r = rects[bd.block];
        bool ok = true;
        if ((bd.edges & lab3::EDGE_LEFT) && r.x != 0) ok = false;
        if ((bd.edges & lab3::EDGE_RIGHT) && r.x + r.w != inst.W) ok = false;
        if ((bd.edges & lab3::EDGE_BOTTOM) && r.y != 0) ok = false;
        if ((bd.edges & lab3::EDGE_TOP) && r.y + r.h != inst.H) ok = false;
        if (ok) ++bsat;
    }
    double b_score = p == 0 ? 3.0 : 3.0 * (double)bsat / (double)p;

    const size_t q = inst.mibs.size();
    size_t msat = 0;
    for (const auto &grp : inst.mibs) {
        bool ok = true;
        for (int bi : grp)
            if (rects[bi].w != rects[grp[0]].w || rects[bi].h != rects[grp[0]].h) ok = false;
        if (ok) ++msat;
    }
    double m_score = q == 0 ? 2.0 : 2.0 * (double)msat / (double)q;

    double total = 10.0 + wl_score + b_score + m_score;

    std::printf("Result: FEASIBLE\n");
    std::printf("Reported WL   %s\n", reported_str.c_str());
    std::printf("Computed WL   %s\n", lab3::format_half(W2).c_str());
    std::printf("Threshold T   %s\n", lab3::format_half(T2).c_str());
    if (!wl_match) std::printf("WL mismatch: |reported - computed| > 0.01\n");
    std::printf("Gate          10.00 / 10\n");
    std::printf("WL score      %.2f / 5\n", wl_score);
    std::printf("Boundary      %.2f / 3  (%zu/%zu satisfied)\n", b_score, bsat, p);
    std::printf("MIB           %.2f / 2  (%zu/%zu satisfied)\n", m_score, msat, q);
    std::printf("Total %.2f / 20\n", total);
    return 0;
}
