// polish -- exact WL polish of a legal Lab3 placement on its own topology.
// Usage: polish <input> <in_output> <out_output>
// Exit 0 when the LP polish succeeded, 1 when the input was copied unchanged
// (illegal input or internal check failure), 2 on usage/IO errors.
#include "wl_lp.h"

#include <chrono>
#include <cstdio>
#include <fstream>
#include <string>

int main(int argc, char **argv) {
    if (argc != 4) {
        std::fprintf(stderr, "usage: %s <input> <in_output> <out_output>\n", argv[0]);
        return 2;
    }
    lab3::Instance inst;
    std::string err;
    if (!lab3::parse_input_file(argv[1], inst, err)) {
        std::fprintf(stderr, "input: %s\n", err.c_str());
        return 2;
    }
    std::ifstream f(argv[2]);
    std::string kw, wl;
    long long n = 0;
    if (!f || !(f >> kw >> wl >> kw >> n) || n != (long long)inst.blocks.size()) {
        std::fprintf(stderr, "cannot read placement '%s'\n", argv[2]);
        return 2;
    }
    std::vector<lab3::Rect> rects(inst.blocks.size());
    std::vector<char> seen(inst.blocks.size(), 0);
    for (long long k = 0; k < n; ++k) {
        std::string name;
        lab3::Rect r;
        if (!(f >> name >> r.x >> r.y >> r.w >> r.h)) { std::fprintf(stderr, "bad placement line\n"); return 2; }
        auto it = inst.block_index.find(name);
        if (it == inst.block_index.end() || seen[it->second]) { std::fprintf(stderr, "bad block '%s'\n", name.c_str()); return 2; }
        seen[it->second] = 1;
        rects[it->second] = r;
    }
    auto t0 = std::chrono::steady_clock::now();
    mcf::PolishInfo info;
    std::vector<lab3::Rect> out = mcf::polish(inst, rects, &info);
    double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    std::ofstream o(argv[3]);
    if (!o) { std::fprintf(stderr, "cannot write '%s'\n", argv[3]); return 2; }
    o << "WIRELENGTH " << lab3::format_half(info.wl2_out) << "\nBLOCKS " << out.size() << "\n";
    for (size_t i = 0; i < out.size(); ++i)
        o << inst.blocks[i].name << ' ' << out[i].x << ' ' << out[i].y << ' ' << out[i].w << ' ' << out[i].h << '\n';
    std::printf("%s WL %s -> %s (%.3f%%), TR edges H %d V %d, %.3f ms\n", info.ok ? "polished" : "UNCHANGED",
                lab3::format_half(info.wl2_in).c_str(), lab3::format_half(info.wl2_out).c_str(),
                info.wl2_in ? 100.0 * (double)(info.wl2_in - info.wl2_out) / (double)info.wl2_in : 0.0,
                info.edges_h, info.edges_v, ms);
    return info.ok ? 0 : 1;
}
