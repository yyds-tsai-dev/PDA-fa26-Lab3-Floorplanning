// lab3_io.h -- Lab3 input parser and wirelength computation (header-only).
//
// Shared by the verifier and the TA reference solver. Follows docs/design.md
// §3 (input format) and §5 (wirelength). No exit(), no global state.
//
// Usage:
//   lab3::Instance inst;
//   std::string err;
//   if (!lab3::parse_input_file(path, inst, err)) { /* report err */ }
//   int64_t wl2 = lab3::wirelength2(inst, rects);   // 2*WL, exact
#ifndef LAB3_IO_H
#define LAB3_IO_H

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace lab3 {

enum class Kind { Soft, Fixed, Preplaced };

inline const char *kind_name(Kind k) {
    switch (k) {
    case Kind::Soft: return "SOFT";
    case Kind::Fixed: return "FIXED";
    default: return "PREPLACED";
    }
}

struct Block {
    std::string name;
    Kind kind = Kind::Soft;
    int64_t area = 0;          // SOFT only
    int64_t w = 0, h = 0;      // FIXED, PREPLACED
    int64_t x = 0, y = 0;      // PREPLACED
};

struct Terminal {
    std::string name;
    int64_t x = 0, y = 0;
};

// Endpoint reference: index >= 0 is a block index; index < 0 is terminal
// number (-index - 1).
struct Conn {
    int a = 0, b = 0;
    int64_t weight = 0;
};

enum Edge { EDGE_LEFT = 1, EDGE_RIGHT = 2, EDGE_TOP = 4, EDGE_BOTTOM = 8 };

struct Boundary {
    int block = 0;
    int edges = 0;             // bitmask of Edge
};

struct Rect {
    int64_t x = 0, y = 0, w = 0, h = 0;
};

struct Instance {
    int64_t W = 0, H = 0;
    int64_t T2 = 0;            // 2 * WLTHRESHOLD (T is a multiple of 0.5)
    int64_t R = 0;
    std::vector<Block> blocks;
    std::vector<Terminal> terminals;
    std::vector<Conn> conns;
    std::vector<std::vector<int>> mibs;   // block indices
    std::vector<Boundary> boundaries;
    std::unordered_map<std::string, int> block_index;     // name -> block idx
    std::unordered_map<std::string, int> terminal_index;  // name -> term idx
};

// ---------------------------------------------------------------- helpers

// Strict signed decimal integer: [+-]?[0-9]+ that fits in int64.
inline bool parse_int64(const std::string &s, int64_t &out) {
    if (s.empty()) return false;
    size_t i = 0;
    bool neg = false;
    if (s[0] == '+' || s[0] == '-') { neg = (s[0] == '-'); i = 1; }
    if (i >= s.size()) return false;
    unsigned __int128 v = 0;
    for (; i < s.size(); ++i) {
        char c = s[i];
        if (c < '0' || c > '9') return false;
        v = v * 10 + (unsigned)(c - '0');
        if (v > (unsigned __int128)INT64_MAX + 1) return false;
    }
    if (!neg && v > (unsigned __int128)INT64_MAX) return false;
    out = neg ? (int64_t)(-(__int128)v) : (int64_t)v;
    return true;
}

// WLTHRESHOLD: non-negative multiple of 0.5 printed with exactly one decimal,
// i.e. [0-9]+\.[05]. Returns 2*T.
inline bool parse_threshold2(const std::string &s, int64_t &t2) {
    size_t dot = s.find('.');
    if (dot == std::string::npos || dot == 0 || dot + 2 != s.size()) return false;
    std::string ip = s.substr(0, dot);
    for (char c : ip) if (c < '0' || c > '9') return false;
    char f = s[dot + 1];
    if (f != '0' && f != '5') return false;
    int64_t v;
    if (!parse_int64(ip, v)) return false;
    if (v > INT64_MAX / 4) return false;
    t2 = v * 2 + (f == '5' ? 1 : 0);
    return true;
}

inline bool parse_edge(const std::string &s, int &e) {
    if (s == "LEFT") e = EDGE_LEFT;
    else if (s == "RIGHT") e = EDGE_RIGHT;
    else if (s == "TOP") e = EDGE_TOP;
    else if (s == "BOTTOM") e = EDGE_BOTTOM;
    else return false;
    return true;
}

inline std::string edge_names(int edges) {
    std::string r;
    const char *nm[4] = {"LEFT", "RIGHT", "TOP", "BOTTOM"};
    for (int i = 0; i < 4; ++i)
        if (edges & (1 << i)) { if (!r.empty()) r += ' '; r += nm[i]; }
    return r;
}

// Token stream over the whole input; records token index for messages.
class Tokens {
public:
    explicit Tokens(std::istream &in) {
        std::string t;
        while (in >> t) toks_.push_back(t);
    }
    bool eof() const { return pos_ >= toks_.size(); }
    size_t pos() const { return pos_; }
    const std::string &peek() const { return toks_[pos_]; }
    bool next(std::string &t) {
        if (eof()) return false;
        t = toks_[pos_++];
        return true;
    }
private:
    std::vector<std::string> toks_;
    size_t pos_ = 0;
};

// ---------------------------------------------------------------- parser

// Parses and validates a Lab3 input (design.md §3). On failure returns false
// and sets err to a human-readable message; inst is then unspecified.
inline bool parse_input(std::istream &in, Instance &inst, std::string &err) {
    inst = Instance();
    Tokens tk(in);
    std::string tok;

    auto fail = [&](const std::string &m) {
        err = m;
        return false;
    };
    auto expect_kw = [&](const char *kw) -> bool {
        if (!tk.next(tok)) { err = std::string("unexpected end of file, expected '") + kw + "'"; return false; }
        if (tok != kw) { err = std::string("expected '") + kw + "', got '" + tok + "'"; return false; }
        return true;
    };
    auto read_int = [&](const char *what, int64_t &v) -> bool {
        if (!tk.next(tok)) { err = std::string("unexpected end of file, expected ") + what; return false; }
        if (!parse_int64(tok, v)) { err = std::string("expected integer ") + what + ", got '" + tok + "'"; return false; }
        return true;
    };
    auto read_name = [&](const char *what, std::string &s) -> bool {
        if (!tk.next(s)) { err = std::string("unexpected end of file, expected ") + what; return false; }
        return true;
    };
    auto read_count = [&](const char *kw, int64_t &c) -> bool {
        if (!expect_kw(kw)) return false;
        if (!read_int((std::string(kw) + " count").c_str(), c)) return false;
        if (c < 0) { err = std::string(kw) + " count is negative"; return false; }
        if (c > 10000000) { err = std::string(kw) + " count is unreasonably large"; return false; }
        return true;
    };
    auto reserved_name = [](const std::string &s) {
        return s == "LEFT" || s == "RIGHT" || s == "TOP" || s == "BOTTOM";
    };
    auto name_used = [&](const std::string &s) {
        return inst.block_index.count(s) || inst.terminal_index.count(s);
    };

    // OUTLINE
    if (!expect_kw("OUTLINE")) return false;
    if (!read_int("outline W", inst.W) || !read_int("outline H", inst.H)) return false;
    if (inst.W <= 0 || inst.H <= 0) return fail("OUTLINE W and H must be positive");

    // WLTHRESHOLD
    if (!expect_kw("WLTHRESHOLD")) return false;
    if (!tk.next(tok)) return fail("unexpected end of file, expected WLTHRESHOLD value");
    if (!parse_threshold2(tok, inst.T2))
        return fail("WLTHRESHOLD must be a non-negative multiple of 0.5 with one decimal, got '" + tok + "'");

    // ASPECTRATIO
    if (!expect_kw("ASPECTRATIO")) return false;
    if (!read_int("aspect ratio R", inst.R)) return false;
    if (inst.R < 1) return fail("ASPECTRATIO must be >= 1");

    // BLOCKS
    int64_t n;
    if (!read_count("BLOCKS", n)) return false;
    for (int64_t i = 0; i < n; ++i) {
        Block b;
        if (!tk.next(tok)) return fail("unexpected end of file in BLOCKS section");
        if (!read_name("block name", b.name)) return false;
        if (reserved_name(b.name)) return fail("block name '" + b.name + "' is reserved (LEFT/RIGHT/TOP/BOTTOM)");
        if (name_used(b.name)) return fail("duplicate name '" + b.name + "'");
        if (tok == "SOFT") {
            b.kind = Kind::Soft;
            if (!read_int("area", b.area)) return false;
            if (b.area < 10000) return fail("SOFT " + b.name + ": area must be >= 10000");
        } else if (tok == "FIXED" || tok == "PREPLACED") {
            b.kind = (tok == "FIXED") ? Kind::Fixed : Kind::Preplaced;
            if (!read_int("width", b.w) || !read_int("height", b.h)) return false;
            if (b.w <= 0 || b.h <= 0) return fail(std::string(kind_name(b.kind)) + " " + b.name + ": width and height must be positive");
            if (b.kind == Kind::Preplaced) {
                if (!read_int("x", b.x) || !read_int("y", b.y)) return false;
                if (b.x < 0 || b.y < 0 || (__int128)b.x + b.w > inst.W || (__int128)b.y + b.h > inst.H)
                    return fail("PREPLACED " + b.name + " is not inside the outline");
                if (b.x == 0 || b.y == 0 || (__int128)b.x + b.w == inst.W || (__int128)b.y + b.h == inst.H)
                    return fail("PREPLACED " + b.name + " touches an outline edge (must be strictly inside)");
            }
        } else {
            return fail("expected SOFT, FIXED or PREPLACED record " + std::to_string(i + 1) + " of " +
                        std::to_string(n) + ", got '" + tok + "'");
        }
        inst.block_index[b.name] = (int)inst.blocks.size();
        inst.blocks.push_back(b);
    }
    // Preplaced pairwise non-overlap.
    for (size_t i = 0; i < inst.blocks.size(); ++i) {
        const Block &a = inst.blocks[i];
        if (a.kind != Kind::Preplaced) continue;
        for (size_t j = i + 1; j < inst.blocks.size(); ++j) {
            const Block &b = inst.blocks[j];
            if (b.kind != Kind::Preplaced) continue;
            bool ox = std::max(a.x, b.x) < std::min(a.x + a.w, b.x + b.w);
            bool oy = std::max(a.y, b.y) < std::min(a.y + a.h, b.y + b.h);
            if (ox && oy) return fail("PREPLACED blocks overlap: " + a.name + " " + b.name);
        }
    }

    // TERMINALS
    int64_t r;
    if (!read_count("TERMINALS", r)) return false;
    for (int64_t i = 0; i < r; ++i) {
        Terminal t;
        if (!expect_kw("TERMINAL")) return false;
        if (!read_name("terminal name", t.name)) return false;
        if (reserved_name(t.name)) return fail("terminal name '" + t.name + "' is reserved (LEFT/RIGHT/TOP/BOTTOM)");
        if (name_used(t.name)) return fail("duplicate name '" + t.name + "'");
        if (!read_int("terminal x", t.x) || !read_int("terminal y", t.y)) return false;
        bool onx = (t.x == 0 || t.x == inst.W) && t.y >= 0 && t.y <= inst.H;
        bool ony = (t.y == 0 || t.y == inst.H) && t.x >= 0 && t.x <= inst.W;
        if (!onx && !ony) return fail("TERMINAL " + t.name + " is not on the outline border");
        inst.terminal_index[t.name] = (int)inst.terminals.size();
        inst.terminals.push_back(t);
    }

    // CONNECTIONS
    int64_t m;
    if (!read_count("CONNECTIONS", m)) return false;
    {
        std::unordered_map<std::string, int> seen_pair;
        for (int64_t i = 0; i < m; ++i) {
            std::string n1, n2;
            int64_t wgt;
            if (!expect_kw("CONN")) return false;
            if (!read_name("connection endpoint", n1) || !read_name("connection endpoint", n2)) return false;
            if (!read_int("connection weight", wgt)) return false;
            auto endpoint = [&](const std::string &s, int &idx) -> bool {
                auto it = inst.block_index.find(s);
                if (it != inst.block_index.end()) { idx = it->second; return true; }
                auto jt = inst.terminal_index.find(s);
                if (jt != inst.terminal_index.end()) { idx = -jt->second - 1; return true; }
                return false;
            };
            Conn c;
            if (!endpoint(n1, c.a)) return fail("CONN " + n1 + " " + n2 + ": unknown name '" + n1 + "'");
            if (!endpoint(n2, c.b)) return fail("CONN " + n1 + " " + n2 + ": unknown name '" + n2 + "'");
            if (c.a < 0 && c.b < 0) return fail("CONN " + n1 + " " + n2 + ": terminal-terminal connection");
            if (c.a == c.b) return fail("CONN " + n1 + " " + n2 + ": block connected to itself");
            if (wgt < 1 || wgt > 10) return fail("CONN " + n1 + " " + n2 + ": weight must be in [1,10]");
            std::string key = (n1 < n2) ? n1 + '\n' + n2 : n2 + '\n' + n1;
            if (seen_pair.count(key)) return fail("CONN " + n1 + " " + n2 + ": pair appears more than once");
            seen_pair[key] = 1;
            c.weight = wgt;
            inst.conns.push_back(c);
        }
    }

    // MIBS
    int64_t q;
    if (!read_count("MIBS", q)) return false;
    {
        std::vector<int> in_group(inst.blocks.size(), -1);
        for (int64_t g = 0; g < q; ++g) {
            int64_t k;
            if (!expect_kw("MIB")) return false;
            if (!read_int("MIB size k", k)) return false;
            if (k < 2) return fail("MIB group " + std::to_string(g + 1) + ": k must be >= 2");
            if (k > (int64_t)inst.blocks.size())
                return fail("MIB group " + std::to_string(g + 1) + ": k exceeds the number of blocks");
            std::vector<int> mem;
            for (int64_t j = 0; j < k; ++j) {
                std::string s;
                if (!read_name("MIB member", s)) return false;
                auto it = inst.block_index.find(s);
                if (it == inst.block_index.end()) return fail("MIB member '" + s + "' is not a block");
                int bi = it->second;
                if (inst.blocks[bi].kind != Kind::Soft) return fail("MIB member '" + s + "' is not a soft block");
                if (in_group[bi] >= 0) return fail("block '" + s + "' belongs to more than one MIB group (or is listed twice)");
                in_group[bi] = (int)g;
                if (!mem.empty() && inst.blocks[bi].area != inst.blocks[mem[0]].area)
                    return fail("MIB member '" + s + "' has a different area target than '" + inst.blocks[mem[0]].name + "'");
                mem.push_back(bi);
            }
            inst.mibs.push_back(mem);
        }
    }

    // BOUNDARIES
    int64_t p;
    if (!read_count("BOUNDARIES", p)) return false;
    {
        std::vector<char> has(inst.blocks.size(), 0);
        for (int64_t i = 0; i < p; ++i) {
            std::string s;
            if (!expect_kw("BOUNDARY")) return false;
            if (!read_name("boundary block name", s)) return false;
            auto it = inst.block_index.find(s);
            if (it == inst.block_index.end()) return fail("BOUNDARY '" + s + "' is not a block");
            int bi = it->second;
            if (inst.blocks[bi].kind == Kind::Preplaced) return fail("BOUNDARY '" + s + "' is a preplaced block");
            if (has[bi]) return fail("block '" + s + "' appears in more than one BOUNDARY record");
            has[bi] = 1;
            int e1, e2;
            if (!tk.next(tok)) return fail("unexpected end of file, expected EDGE for BOUNDARY " + s);
            if (!parse_edge(tok, e1)) return fail("BOUNDARY " + s + ": invalid EDGE '" + tok + "'");
            int edges = e1;
            // Optional second edge: present iff the next token is an EDGE keyword.
            if (!tk.eof() && parse_edge(tk.peek(), e2)) {
                tk.next(tok);
                if (e1 == e2) return fail("BOUNDARY " + s + ": repeated edge " + tok);
                int both = e1 | e2;
                if (both == (EDGE_LEFT | EDGE_RIGHT) || both == (EDGE_TOP | EDGE_BOTTOM))
                    return fail("BOUNDARY " + s + ": edges must be adjacent (a corner)");
                edges = both;
            }
            Boundary bd;
            bd.block = bi;
            bd.edges = edges;
            inst.boundaries.push_back(bd);
        }
    }

    if (!tk.eof()) return fail("trailing token after BOUNDARIES section: '" + tk.peek() + "'");
    return true;
}

inline bool parse_input_file(const std::string &path, Instance &inst, std::string &err) {
    std::ifstream f(path);
    if (!f) { err = "cannot open input file '" + path + "'"; return false; }
    return parse_input(f, inst, err);
}

// ---------------------------------------------------------------- wirelength

// Returns 2*WL exactly (design.md §5). rects[i] is the placement of block i
// (input order). Uses doubled centroids: 2cx = 2x + w, 2cy = 2y + h.
inline int64_t wirelength2(const Instance &inst, const std::vector<Rect> &rects) {
    __int128 sum = 0;
    auto pt = [&](int e, __int128 &x2, __int128 &y2) {
        if (e >= 0) {
            const Rect &r = rects[e];
            x2 = (__int128)2 * r.x + r.w;
            y2 = (__int128)2 * r.y + r.h;
        } else {
            const Terminal &t = inst.terminals[-e - 1];
            x2 = (__int128)2 * t.x;
            y2 = (__int128)2 * t.y;
        }
    };
    for (const Conn &c : inst.conns) {
        __int128 ax, ay, bx, by;
        pt(c.a, ax, ay);
        pt(c.b, bx, by);
        __int128 dx = ax > bx ? ax - bx : bx - ax;
        __int128 dy = ay > by ? ay - by : by - ay;
        sum += (__int128)c.weight * (dx + dy);
    }
    return (int64_t)sum;
}

// Formats a doubled value v2 (= 2*v) as a decimal with one digit: "123.5".
inline std::string format_half(int64_t v2) {
    bool neg = v2 < 0;
    unsigned long long a = neg ? (unsigned long long)(-(v2 + 1)) + 1 : (unsigned long long)v2;
    std::string s = std::to_string(a / 2) + ((a % 2) ? ".5" : ".0");
    return neg ? "-" + s : s;
}

}  // namespace lab3

#endif  // LAB3_IO_H
