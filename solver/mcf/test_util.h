// test_util.h -- helpers for the mcf tests: output-file I/O, a Python-
// compatible random.Random (to reproduce tools/lp_polish_proto.py topology
// variants), zlib crc32, a minimal JSON reader, and a verifier runner.
#ifndef MCF_TEST_UTIL_H
#define MCF_TEST_UTIL_H

#include "lab3_io.h"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace tu {

// Reads a Lab3 output file into rects in instance block order.
inline bool read_output(const std::string &path, const lab3::Instance &inst, std::vector<lab3::Rect> &r) {
    std::ifstream f(path);
    if (!f) return false;
    std::string kw, wl;
    long long n;
    if (!(f >> kw >> wl >> kw >> n)) return false;
    r.assign(inst.blocks.size(), lab3::Rect());
    for (long long k = 0; k < n; ++k) {
        std::string name;
        lab3::Rect t;
        if (!(f >> name >> t.x >> t.y >> t.w >> t.h)) return false;
        auto it = inst.block_index.find(name);
        if (it == inst.block_index.end()) return false;
        r[it->second] = t;
    }
    return n == (long long)inst.blocks.size();
}

inline bool write_output(const std::string &path, const lab3::Instance &inst, const std::vector<lab3::Rect> &r) {
    std::ofstream f(path);
    if (!f) return false;
    f << "WIRELENGTH " << lab3::format_half(lab3::wirelength2(inst, r)) << "\nBLOCKS " << inst.blocks.size() << "\n";
    for (size_t i = 0; i < r.size(); ++i)
        f << inst.blocks[i].name << ' ' << r[i].x << ' ' << r[i].y << ' ' << r[i].w << ' ' << r[i].h << '\n';
    return (bool)f;
}

// ------------------------------------------------- Python random.Random

inline uint32_t crc32(const std::string &s) {
    uint32_t c = 0xFFFFFFFFu;
    for (unsigned char ch : s) {
        c ^= ch;
        for (int k = 0; k < 8; ++k) c = (c >> 1) ^ (0xEDB88320u & (0u - (c & 1u)));
    }
    return c ^ 0xFFFFFFFFu;
}

// MT19937 seeded like CPython's random.seed(int) for 0 <= seed < 2^32.
class PyRandom {
public:
    explicit PyRandom(uint32_t seed) {
        init_genrand(19650218u);
        uint32_t key[1] = {seed};
        int i = 1, j = 0, klen = 1;
        for (int k = N > klen ? N : klen; k; --k) {
            mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1664525u)) + key[j] + (uint32_t)j;
            ++i; ++j;
            if (i >= N) { mt[0] = mt[N - 1]; i = 1; }
            if (j >= klen) j = 0;
        }
        for (int k = N - 1; k; --k) {
            mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1566083941u)) - (uint32_t)i;
            ++i;
            if (i >= N) { mt[0] = mt[N - 1]; i = 1; }
        }
        mt[0] = 0x80000000u;
    }
    double random() {
        uint32_t a = next() >> 5, b = next() >> 6;
        return (a * 67108864.0 + b) * (1.0 / 9007199254740992.0);
    }

private:
    static constexpr int N = 624, M = 397;
    uint32_t mt[N];
    int mti = N + 1;
    void init_genrand(uint32_t s) {
        mt[0] = s;
        for (mti = 1; mti < N; ++mti) mt[mti] = 1812433253u * (mt[mti - 1] ^ (mt[mti - 1] >> 30)) + (uint32_t)mti;
    }
    uint32_t next() {
        static const uint32_t mag01[2] = {0u, 0x9908b0dfu};
        if (mti >= N) {
            int kk;
            for (kk = 0; kk < N - M; ++kk) {
                uint32_t y = (mt[kk] & 0x80000000u) | (mt[kk + 1] & 0x7fffffffu);
                mt[kk] = mt[kk + M] ^ (y >> 1) ^ mag01[y & 1u];
            }
            for (; kk < N - 1; ++kk) {
                uint32_t y = (mt[kk] & 0x80000000u) | (mt[kk + 1] & 0x7fffffffu);
                mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ mag01[y & 1u];
            }
            uint32_t y = (mt[N - 1] & 0x80000000u) | (mt[0] & 0x7fffffffu);
            mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ mag01[y & 1u];
            mti = 0;
        }
        uint32_t y = mt[mti++];
        y ^= (y >> 11);
        y ^= (y << 7) & 0x9d2c5680u;
        y ^= (y << 15) & 0xefc60000u;
        y ^= (y >> 18);
        return y;
    }
};

// ------------------------------------------------- minimal JSON

struct Json {
    enum T { NUL, BOOL, NUM, STR, ARR, OBJ } t = NUL;
    double num = 0;
    std::string str;
    std::vector<Json> arr;
    std::map<std::string, Json> obj;
    const Json &operator[](const std::string &k) const {
        static Json none;
        auto it = obj.find(k);
        return it == obj.end() ? none : it->second;
    }
    bool has(const std::string &k) const { return obj.count(k) > 0; }
};

class JsonParser {
public:
    explicit JsonParser(const std::string &s) : s_(s) {}
    bool parse(Json &out) { p_ = 0; return value(out); }

private:
    const std::string &s_;
    size_t p_ = 0;
    void ws() { while (p_ < s_.size() && (s_[p_] == ' ' || s_[p_] == '\n' || s_[p_] == '\t' || s_[p_] == '\r')) ++p_; }
    bool lit(const char *w) {
        size_t n = std::string(w).size();
        if (s_.compare(p_, n, w) != 0) return false;
        p_ += n;
        return true;
    }
    bool string(std::string &o) {
        if (s_[p_] != '"') return false;
        ++p_;
        while (p_ < s_.size() && s_[p_] != '"') {
            if (s_[p_] == '\\') ++p_;
            o += s_[p_++];
        }
        ++p_;
        return true;
    }
    bool value(Json &v) {
        ws();
        if (p_ >= s_.size()) return false;
        char c = s_[p_];
        if (c == '{') {
            v.t = Json::OBJ; ++p_; ws();
            if (s_[p_] == '}') { ++p_; return true; }
            for (;;) {
                ws();
                std::string k;
                if (!string(k)) return false;
                ws();
                if (s_[p_++] != ':') return false;
                if (!value(v.obj[k])) return false;
                ws();
                if (s_[p_] == ',') { ++p_; continue; }
                if (s_[p_] == '}') { ++p_; return true; }
                return false;
            }
        }
        if (c == '[') {
            v.t = Json::ARR; ++p_; ws();
            if (s_[p_] == ']') { ++p_; return true; }
            for (;;) {
                v.arr.emplace_back();
                if (!value(v.arr.back())) return false;
                ws();
                if (s_[p_] == ',') { ++p_; continue; }
                if (s_[p_] == ']') { ++p_; return true; }
                return false;
            }
        }
        if (c == '"') { v.t = Json::STR; return string(v.str); }
        if (lit("true")) { v.t = Json::BOOL; v.num = 1; return true; }
        if (lit("false")) { v.t = Json::BOOL; return true; }
        if (lit("null")) { v.t = Json::NUL; return true; }
        if (lit("NaN")) { v.t = Json::NUM; return true; }
        v.t = Json::NUM;
        char *end;
        v.num = std::strtod(s_.c_str() + p_, &end);
        if (end == s_.c_str() + p_) return false;
        p_ = end - s_.c_str();
        return true;
    }
};

// ------------------------------------------------- verifier

struct VerifierResult {
    bool feasible = false, bnd_full = false, mib_full = false;
    std::string computed_wl;
};

inline VerifierResult run_verifier(const std::string &verifier, const std::string &in, const std::string &out) {
    VerifierResult r;
    std::string cmd = "'" + verifier + "' '" + in + "' '" + out + "' 2>&1";
    FILE *p = popen(cmd.c_str(), "r");
    if (!p) return r;
    std::string txt;
    char buf[4096];
    size_t k;
    while ((k = fread(buf, 1, sizeof buf, p)) > 0) txt.append(buf, k);
    int rc = pclose(p);
    std::istringstream ss(txt);
    std::string line;
    bool res = false;
    while (std::getline(ss, line)) {
        if (line.find("Result: FEASIBLE") != std::string::npos) res = true;
        if (line.rfind("Computed WL", 0) == 0) {
            std::istringstream ls(line.substr(11));
            ls >> r.computed_wl;
        }
        auto full = [&](const char *tag) {
            if (line.rfind(tag, 0) != 0) return false;
            size_t a = line.find('('), sl = line.find('/', a), b = line.find(' ', sl);
            if (a == std::string::npos || sl == std::string::npos) return false;
            return line.substr(a + 1, sl - a - 1) == line.substr(sl + 1, b - sl - 1);
        };
        if (full("Boundary")) r.bnd_full = true;
        if (full("MIB")) r.mib_full = true;
    }
    r.feasible = rc == 0 && res;
    return r;
}

}  // namespace tu

#endif  // MCF_TEST_UTIL_H
