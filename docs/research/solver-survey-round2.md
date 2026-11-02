# Solver survey, round 2: search-side literature for the post-feasible anneal

Status: research note, 2026-10-29 23:10 (Asia/Taipei). Companion to `solver-survey.md` (round 1, papers P1–P27, architecture ranking) and to the measured state in `adhd-solver-round2.md`, `adhd-solver-round3.md`, `legalizer-nogo.md` and `lp-polish-results.md`. Round-1 papers are not re-surveyed; they are cited by their P-numbers. New papers are numbered R1–R57.

Conventions (same as round 1): **[full]** full text read (or the relevant sections), **[abs]** abstract or bibliographic record only, **[sec]** claim taken from another paper's description, **[spec]** derived from `docs/design.md`, **[proposal]** my suggestion, **(uncertain)** not verified. **[meas]** marks a number measured on this project (from the docs above).

Downloaded texts (untrusted data, read only): scratchpad `papers2/` (`dl/` for PDFs and extracted text, `arxiv/` for arXiv texts, `abs/` for OpenAlex/Semantic Scholar records, `searches2/` for the raw search listings).

---

## 0. The problem this round is asked to solve, in measured numbers

The survey is checked against these facts, not against general plausibility.

- **M1 [meas]** Fixed-outline, n ≤ 200, 5–8 % whitespace, ~60 % soft (integer w,h, area ±1 %, aspect ≤ 3), ~28 % fixed-shape, ~12 % preplaced in the interior, boundary constraints on ~20 % of blocks, MIB groups, weighted 2-pin centroid Manhattan WL, 110 s single-thread C++, no libraries.
- **M2 [meas]** Current flow: analytical seed → sequence-pair SA with penalties → min-cost-flow LP polish. SA from its own seeds ends at 1.2–1.6× the known-good (construction) WL; SA started from the known-good topology reaches 0.92–0.95×. So the search, not the representation, is the limit.
- **M3 [meas]** Seed spread (max/median up to 1.3) is created in the post-feasible anneal between first-legal and ~50 s. Continuations forked at first-legal spread as much as independent runs (1.41–1.53 vs 1.33–1.43); forked at 50 s they spread only 1.01–1.14. A 2×55 s or 3×36 s restart portfolio cuts max/median (1.26 → 1.14–1.15) but leaves the median flat (1.41 → 1.40–1.44).
- **M4 [meas]** The exact per-axis LP on a frozen topology gains 2.9 % median (1.8 % at n = 200). Across near-identical topologies, packed WL anti-correlates with LP WL (Spearman median −0.19; the LP-best variant had the worst packed WL in 14/20 cases). HiGHS takes 11 ms per axis; a C++ network simplex is estimated at 0.1–0.5 ms.
- **M5 [meas]** Dead ends with numbers: constraint-graph legalization of GP centres (legal on 0/56 hidden4 seeds, 250–340 ms per repair iteration); guillotine partition around 24 interior preplaced blocks (about 90 regions, ~60 cannot fit their blocks); maximal-empty-rectangle constructive greedy (4–6 % area unplaced, 1.96–3.51× WL); single-block rip-up/re-insert with an exact oracle (4–11 % of blocks have an improving cell, total gain 0.03–1.1 %); elite pool/consensus (variance is late-tail, not early-topology).
- **M6 [meas]** New assets: a coordinate-space overlap-minimization legalizer with guided-local-search weights (legal in 2–10 s at 5–6 % whitespace, hand-off WL ≈ 2.1× the known-good solution), and the round-3 feasibility homotopy (anneal with movable areas scaled to 0.92, GLS absorber at s = 1, then the true anneal): final WL 0.845–0.88× the 110 s control median, p90 1.02–1.05 vs 1.15, legal on 47/48. The displacement pre-test says the 1.2–1.6× gap is region-scale (55–66 % of the WL excess sits on blocks displaced > 2 block-scales), so single-block and small-window moves cannot close it.

Reading guide: a paper is worth stealing from only if it (a) changes the basin the anneal ends in, (b) changes what the late anneal optimises (the geometry the acceptance test sees), or (c) cuts the late-tail variance without cutting the median. Everything else is recorded for completeness.

---

## 1. Findings by direction

### 1.1 Local search over sequence pairs with cost-optimal coordinates (direction 1)

**R1. Imahori, Yagiura, Ibaraki, "Local search algorithms for the rectangle packing problem with general spatial costs"** [abs]. Mathematical Programming 97(3):543–569, 2003. DOI 10.1007/s10107-003-0427-1.
Mechanism: each rectangle has a set of modes (shapes/orientations) and piecewise-linear spatial cost functions of its x and y coordinates; the costs may be non-convex and discontinuous. The search state is a sequence pair plus a mode vector; for a fixed sequence pair the cost-optimal coordinates are computed by dynamic programming over the constraint order (the decoder is a generalisation of the longest-path packing to arbitrary costs). Local search and metaheuristics over several neighbourhoods (shift, swap, mode change; the exact set is (uncertain) from the record). Area minimisation, strip packing and scheduling are special cases.
Steal / why not: this is the published form of "evaluate an SP by its cost-optimal coordinates, not by bottom-left packing" (M4 says our packer ranks neighbouring topologies the wrong way). For our separable L1 objective the "general spatial cost" decoder is exactly the per-axis LP/MCF we already have; the paper shows it is standard practice to run the SP search on top of it. Full text not accessible (Kyoto page gone; the KURENAI copy of the thesis is a 2-page scan), so parameters are not extracted.

**R2. Imahori, Yagiura, Ibaraki, "Improved local search algorithms for the rectangle packing problem with general spatial costs"** [abs]. European Journal of Operational Research 167(1):48–67, 2005. DOI 10.1016/j.ejor.2004.02.020.
Mechanism: speed-up techniques for evaluating all solutions in the neighbourhoods of R1 (per the record: "evaluating solutions in various neighbourhoods" more efficiently). The neighbourhoods are reported in secondary sources to be restricted by the critical paths of the current solution (uncertain).
Steal / why not: the idea that matters is incremental evaluation of an LP-decoded SP neighbourhood; with n = 200 and about 2.5 k arcs per axis, a warm-started network simplex (few pivots after a local SP change) is our equivalent. Companion: Ibaraki, Imahori, Yagiura, "Hybrid metaheuristics for packing problems", Studies in Computational Intelligence 114, 2008, DOI 10.1007/978-3-540-78295-7_7 [abs], a survey of the same group's overlap-minimisation, compaction and separation methods.

**R3. Tang, Tian, Wong, "Minimizing wire length in floorplanning"** [abs]. IEEE TCAD 25(9):1744–1753, 2006, DOI 10.1109/TCAD.2005.858266; conference version "Optimal redistribution of white space for wire length minimization", ASP-DAC 2005, DOI 10.1109/ASPDAC.2005.1466198.
Mechanism: for a fixed floorplan topology (constraint graphs), minimising total wirelength while keeping the compacted area is an LP; its dual is solved as a min-cost flow. Extensions listed in the abstract without loss of optimality: fixed frame, IO pins, preplaced blocks, boundary blocks, range placement, alignment and abutment, rectilinear blocks, soft blocks, 1-D cluster placement, bounded net delay. Runtime < 0.4 s for all MCNC block benchmarks; 4.2 % WL improvement "even on very compact floorplans".
Steal / why not: independent confirmation of M4 (our 1.8–4.2 % post-pass is the same effect). More useful: the extension list says preplaced, boundary and soft blocks all stay inside the min-cost-flow dual, which is what our polish already does. Also the natural place to add proximal anchors (round 3, cluster IV) as extra hinge terms, free of charge.

**R4. Tang, Wong, "Floorplanning with alignment and performance constraints"** [abs]. DAC 2002, DOI 10.1109/DAC.2002.1012740.
Mechanism: SP floorplanner designed for fixed-frame floorplanning; alignment (bus) and bounded-delay constraints; feasibility of an SP with O(n) constrained blocks evaluated in O(n log log n) by weighted LCS instead of O(n³) constraint-graph work.
Steal / why not: the weighted-LCS trick is the fast path for our constraining edges (round-1 P6 uses explicit graphs). Not the bottleneck now (pack ≈ 4 µs [meas]).

**R5. Chen, Yoshimura, "Fixed-outline floorplanning: block-position enumeration and a new method for calculating area costs"** [abs]. IEEE TCAD 27(5):858–871, 2008, DOI 10.1109/TCAD.2008.917968; ISPD 2007 version "A stable fixed-outline floorplanning method", DOI 10.1145/1231996.1232021.
Mechanism: "insertion after remove" (IAR): remove a block from the SP and enumerate all its re-insertion positions, skipping those that violate the fixed outline; a new area cost term. Reported: 100 % fixed-outline success and 20 % lower WL than the fixed-outline floorplanners of the time, with less runtime.
Steal / why not: IAR is our exact insertion oracle (`solver/oracle.h`) used as a perturbation. Measured post hoc on converged outputs it is weak (M5); the in-anneal use (an SA move, not an LNS) is the only untested form and is folded into experiment X4 below. Not ranked on its own.

**R6. Silvanus, Vygen, "Few sequence pairs suffice: representing all rectangle placements"** [abs]. arXiv 1708.09779, 2017.
Mechanism: complete sets of pairwise spatial relations much smaller than the (n!)² sequence pairs, via pattern-avoiding permutations. Theoretical.
Steal / why not: none directly; it quantifies how redundant the SP space is, which is one reason SP moves are often neutral.

**R7. Vorwerk, Kennings, Anjos, "VLSI floorplan repair using dynamic white-space management, constraint graphs, and linear programming"** [abs]. Engineering Optimization 40(7), 2008, DOI 10.1080/03052150801901475.
Mechanism: top-down overlap removal for rough floorplans: geometric shifting plus topological constraint graphs plus LP inside a top-down flow. Reported vs leading methods: 1/5 of the runtime, 4–13 % less WL, up to 43 % less movement; valid placements "in almost all cases".
Steal / why not: same class as our constraint-graph legalizer, which failed on GP centres at 6–8 % whitespace (M5). The reported inputs are rough floorplans with more slack. Not pursued.

**R8. Moffitt, Ng, Markov, Pollack, "Constraint-driven floorplan repair"** (Floorist) [abs]. DAC 2006, DOI 10.1109/DAC.2006.229405 (round-1 P16, then [sec]).
Mechanism: graph-based constraint encoding (outline, fixed locations, regions, proximity, alignment); tracks implications and resolves violations by gradual modifications rather than randomised search.
Steal / why not: see R7; same regime problem.

Also noted [sec] via R28: Murata, Kuh, "Sequence-pair based placement method for hard/soft/pre-placed modules", ISPD 1998, uses convex programming to shape soft modules for each SP and is reported as slow for that reason. Round-1 P7 (Murata 1998 obstacles) remains the SP-adaptation reference for preplaced blocks.

### 1.2 Guided local search and overlap minimisation for packing (direction 2)

**R9. Umetani, Yagiura, Imahori, Imamichi, Nonobe, Ibaraki, "Solving the irregular strip packing problem via guided local search for overlap minimization"** [full, METR 2008-22 tech report]. Int. Trans. in Operational Research 16(6):661–683, 2009, DOI 10.1111/j.1475-3995.2009.00707.x.
Mechanism (details in §2.1): overlap penalty f_ij = min(horizontal, vertical penetration depth) (for rectangles: min of the x- and y-overlap lengths); a block is moved by exact 1-D line searches, alternately horizontal and vertical, each minimising the weighted penalty along the whole line (piecewise linear, breakpoints merged in O(q log n)); first-admissible local search with Fast Local Search activation; GLS weights w_ij ← w_ij + f_ij / max f at each local optimum; outer loop shrinks the strip by 2 % when feasible and extends by 0.5 % otherwise.
Steal / why not: our GLS legalizer is this algorithm family. What to copy exactly: the alternating exact line search (our penalty and WL are both piecewise linear along a coordinate line, so a move can minimise WL + Σ w·overlap exactly), the activation rule (reactivate only blocks that overlapped the moved block before or after), and the weight update. The outer "shrink the container" loop does not map to WL; its analogue for us is the round-3 area-scaling homotopy (see §3).

**R10. Umetani, Murakami, "Coordinate descent heuristics for the irregular strip packing problem of rasterized shapes"** [full, sections 4–6]. EJOR 303(3), 2022, arXiv 2104.04525.
Mechanism: same framework on a raster (integer coordinates); coordinate descent = the same alternating line search with t ∈ Z; GLS weights α_ij initialised to 1.0 and updated as in R9; parameters r_dec = 0.02, r_inc = 0.005, k_max = 200 non-improving GLS iterations; single-threaded C.
Steal / why not: confirms the parameters still hold 13 years later. Integer coordinates match our integer spec.

**R11. Imamichi, Yagiura, Nagamochi, "An iterated local search algorithm based on nonlinear programming for the irregular strip packing problem"** (ILSQN) [full, via Imamichi's doctoral thesis ch. 2]. Discrete Optimization 6(4):345–361, 2009, DOI 10.1016/j.disopt.2009.04.002.
Mechanism: perturb the incumbent by swapping two random pieces (each re-placed at its minimum-overlap position, found by evaluating the penalty at all vertices/intersections of the no-fit polygons; Lemma: finds a zero-overlap position whenever one exists), then "Separate": quasi-Newton descent on the sum of penetration depths with all pieces moving simultaneously; accept if the total overlap improved; stop after N_mo consecutive failures.
Steal / why not: the swap-then-separate step is a two-block coordinated move with a legal re-placement oracle, which is what the region-scale gap (M6) asks for; our oracle is the SP analogue. The quasi-Newton separation is the continuous analogue of our GLS line searches.

**R12. Egeblad, Nielsen, Odgaard, "Fast neighborhood search for two- and three-dimensional nesting problems"** [abs + sec via R9]. EJOR 183(3):1249–1266, 2007, DOI 10.1016/j.ejor.2005.11.063.
Mechanism: GLS on overlap minimisation with the intersection area as penalty; neighbourhood = any horizontal or vertical translation of one polygon, minimised exactly along the line (piecewise quadratic for areas, piecewise linear for our rectangles' 1-D measures). Egeblad's thesis (2008) summarises the strip-packing heuristic as "repeatedly removes overlap and shrinks the container".
Steal / why not: same as R9; the area penalty is smoother than penetration depth but costs more per evaluation. Full text blocked (CORE/KU mirrors behind Cloudflare), so the λ/utility parameters are not extracted.

**R13. Faroe, Pisinger, Zachariasen, "Guided local search for the three-dimensional bin-packing problem"** [abs]. INFORMS J. Computing 15(3):267–283, 2003, DOI 10.1287/ijoc.15.3.267.16080.
Mechanism: start from a greedy upper bound on the number of bins, remove one bin, search a feasible packing with GLS, repeat until the time limit or a lower bound is met; up to 200 boxes.
Steal / why not: the "fix the budget, search feasibility, tighten the budget" outer loop is the same pattern as R9 and as our homotopy; the legality search inside is what our legalizer does.

**R14. Faroe, Pisinger, Zachariasen, "Guided local search for final placement in VLSI design"** [abs]. ICCAD 2001 ("Local search for final placement in VLSI design", DOI 10.1109/ICCAD.2001.968708) and Journal of Heuristics 9(3):269–295, 2003, DOI 10.1023/A:1023721408655; DIKU report 01/01.
Mechanism: GLS with Fast Local Search for general-cell placement, minimising total bounding-box net length with overlap handled by GLS penalties; the FLS sub-neighbourhoods keep the search focused; "could be extended to refining an existing placement". Reported: up to 20 % shorter estimated routing length on large general-cell layouts than existing algorithms.
Steal / why not: the only paper in this group that keeps WL in the objective while removing overlap, i.e. exactly the missing piece between our 2.1× legalizer hand-off and the anneal. The DIKU server was unreachable, so the penalty/λ settings are not extracted; the mechanics can be reconstructed from R9 (same framework) plus a WL term.

**R15. Egeblad, Pisinger, "Heuristic approaches for the two- and three-dimensional knapsack packing problem"** [abs]. Computers & OR 36(4):1026–1049, 2009, DOI 10.1016/j.cor.2007.12.004.
Mechanism: SA over sequence pairs (2-D) and a new "sequence triple" (3-D) for knapsack packing of rectangles.
Steal / why not: nothing beyond what round 1 covers for SP-SA.

How these methods continue past legality: they do not optimise a second objective after a feasible packing is found; they tighten the feasibility budget (strip length, bin count) and re-legalise (R9 §4, R13). Only R14 carries the real objective (WL) inside the GLS cost. For us the "tighten and re-legalise" loop corresponds to the area-scaling homotopy of round 3, and R14 corresponds to running the legalizer with WL + penalties rather than penalties alone.

### 1.3 Analytical and convex fixed-outline floorplanning with preplaced and soft modules (direction 3)

**R16. Lin, Hung, "UFO: unified convex optimization algorithms for fixed-outline floorplanning considering pre-placed modules"** [abs]. IEEE TCAD 30(7):1034–1044, 2011, DOI 10.1109/TCAD.2011.2114531 (round-1 P13).
Mechanism: stage 1 circles with a push-pull convex model spread modules over the outline with wirelength; stage 2 extracts geometric relations into constraint graphs and solves a second-order cone program for locations and shapes with an overlap penalty and the outline as constraints; "it is not necessary to consider wirelength in the legalization phase"; extended to preplaced modules; MATLAB.
Steal / why not: legalisation without WL is exactly what cost us 2.1× (M6). The relation-extraction step is round-1 candidate C's step 2. Needs an SOCP solver; not for our C++.

**R17. Chan, Hsu, Lin, "A flexible fixed-outline floorplanning methodology for mixed-size modules"** [abs]. ASP-DAC 2013, DOI 10.1109/ASPDAC.2013.6509635.
Mechanism: global distribution stage (analytical, WL) then legalisation by shape-curve merging in a slicing tree built by recursive partitioning of circuit and region; 22.5 % less WL than PATOMA and 4.7 % less than DeFer on mixed-size benchmarks.
Steal / why not: partition-based slicing legalisation is what the guillotine experiment tested and lost (M5: ~90 regions, ~60 infeasible). The paper's mixed-size cases do not have interior preplaced blocks. Not pursued.

**R18. Lin, Wu, "F-FM: fixed-outline floorplanning methodology for mixed-size modules considering voltage-island constraint"** [abs]. IEEE TCAD 33(11), 2014, DOI 10.1109/TCAD.2014.2351571.
Mechanism: two-stage analytical + slicing-tree flow (R17 generalised) with voltage islands.
Steal / why not: same as R17.

**R19. Lin, Hung, "SKB-tree: a fixed-outline driven representation for modern floorplanning problems"** [abs]. IEEE TVLSI 20(3), 2012, DOI 10.1109/TVLSI.2011.2104983.
Mechanism: a tree that dynamically allocates regions so every solution fits the outline; SA then optimises WL only; reported better WL than A-FP, Parquet 4.0, ZDS, SAFFOA at zero dead space.
Steal / why not: region allocation needs soft modules to absorb the regions; our 40 % fixed/preplaced blocks and interior pins break it (same failure class as M5's partitions).

**R20. He, Ma, Xu, Wang, Hong, "On handling fixed blocks in incremental fixed-outline floorplanning"** [abs]. ICCCAS 2010, DOI 10.1109/ICCCAS.2010.5581854.
Mechanism: post-process a packing in a TCG: adjust critical paths iteratively to fix outline violations; virtual nodes in the TCG for blocks with fixed positions; ami49 in < 0.3 s; area +1.3 %, WL +19.5 % after repair.
Steal / why not: the +19.5 % WL cost of a topology repair is a warning against post-hoc constraint-graph repair as a WL tool; consistent with M5.

**R21. Lu, Huang, Zhu, "An effective fixed-outline floorplanning algorithm for rectilinear soft modules"** [abs]. ISEDA 2025, DOI 10.1109/ISEDA65950.2025.11100867.
Mechanism: electrostatic global floorplanning with a module-area-growing strategy; negotiation-based legalisation that handles preplaced modules; polygon-clipping refinement for WL.
Steal / why not: the "area-growing" global stage is the analytical cousin of the round-3 homotopy (start with small areas, grow). Rectilinear shapes do not apply. Full text not available.

**R22. Yu, Censor, Luo, "Floorplanning with I/O assignment via feasibility-seeking and superiorization methods"** [full, skimmed]. IEEE TCAD 2024, DOI 10.1109/TCAD.2024.3408106, arXiv 2406.03165.
Mechanism: non-overlap as a union of convex sets; resettable alternating projections (RMAP) for feasibility; superiorisation perturbations steer the iterates toward lower WL (Per-RMAP). MCNC: legal results 166× faster than B&B with +5 % WL; soft modules 15 % faster than PeF; GSRC compared with Parquet-4 and Fast-SA.
Steal / why not: alternating projections are another coordinate-space legalizer; the superiorisation idea (interleave a WL-descent step with feasibility steps) is the same conclusion as R14. Our GLS already does this role; no reason to switch.

**R23. Liu, Dong, Hong, Kajitani, "Fixed-outline floorplanning with constraints through instance augmentation"** [abs]. ISCAS 2005, DOI 10.1109/ISCAS.2005.1464979.
Mechanism: start SA on a sub-instance (fewer modules) and progressively add modules until the full instance is feasible; reported to work "even when tight outline and boundary and pre-placed constraints are imposed".
Steal / why not: a different homotopy (in module count rather than area). The round-3 area homotopy already works (M6); this is a fallback variant if the area scaling stalls on a case.

**R24. Wei, Dong, Hong, Wu, "On handling the fixed-outline constraints of floorplanning using less flexibility first principles"** [abs]. ISCAS 2006, DOI 10.1109/ISCAS.2006.1693903.
Mechanism: constructive LFF packing that keeps the topological relations of a WL reference obtained by partitioning plus multilevel SA.
Steal / why not: constructive packers lost on our instances (M5).

**R25. Huang, Lin, Zhu, Chen, "An improved simulated annealing algorithm with excessive length penalty for fixed-outline floorplanning"** [abs]. IEEE Access 8, 2020, DOI 10.1109/ACCESS.2020.2980135.
Mechanism: SA with a feasible-solution strategy and a two-part outline penalty (violation area plus "excessive violation" that pulls modules toward the outline); 100 % success on all aspect ratios.
Steal / why not: our feasibility is already solved (legal in 2–10 s). Nothing new.

**R26. Zhu, Zheng, Peng, Huang, Wang, Zhu, "Customized diffusion model empowered by heterogeneous graph network for effective floorplanning"** [abs]. IEEE TCAD 2026, DOI 10.1109/TCAD.2026.3681657.
Mechanism: a diffusion model with a heterogeneous GNN generates initial floorplans for a classical analytical floorplanner (PeF lineage); guidance sampling on WL and overlap.
Steal / why not: confirms that the analytical floorplanners of round 1 are "highly sensitive to the quality of initial solutions" (their words), which is what M2 says about our SA. Training pipeline not available to us.

**R27. Sun, Cheng, Wu, Zhu, Chen, "Floorplanning of VLSI by mixed-variable optimization"** [full, skimmed]. arXiv 2401.15317, 2024.
Mechanism: orientations by a distribution evolutionary algorithm, coordinates by a conjugate sub-gradient method on αW + λ√D + μB (D = overlap area, B = outline excess); legalisation by minimising λ₀D + μ₀B̃; claims to beat B*-tree floorplanners on GSRC.
Steal / why not: a population analytical method without a legality guarantee; no preplaced handling. Not pursued.

### 1.4 Deferred decision making and slicing with preplaced or boundary modules (direction 4)

No DeFer (round-1 P19) extension handling preplaced obstacles or boundary constraints was found. The slicing-tree work that does handle them predates DeFer:

**R28. Young, Wong, "Slicing floorplans with pre-placed modules"** [full]. ICCAD 1998, pp. 252–258, DOI 10.1145/288548.288622.
Mechanism (details in §2.3): Wong–Liu Polish-expression SA extended with a reference point per supermodule and four slack numbers (bottom/top/left/right) describing where that point may sit; shape curves of supermodules containing preplaced modules become step functions; two new move classes (swap adjacent operands where one is preplaced; "combine" to move a module out of a congested area); cost A + λW + γD with D the overlap of free modules with preplaced ones. Assumes the preplaced modules themselves form a slicing arrangement (else cut them). MCNC with 3–4 preplaced modules: dead space ≤ 6.6 %, < 6 min on a 300 MHz Pentium II; 50–70-module random cases about 6 min.
Steal / why not: the only published way to keep interior fixed blocks inside a slicing search without partitioning the outline first. Why not for us: 24 interior preplaced blocks (12 % of 200) vs 3–4 in the paper; step-function curves multiply with every preplaced module in a subtree; and the slicing assumption on the preplaced set is not guaranteed by our generator. The measured guillotine failure (M5) is the partition variant of this idea. Recorded as not applicable at our preplaced density.

**R29. Young, Wong, "Slicing floorplans with range constraint"** [abs]. ISPD 1999, DOI 10.1145/299996.300033.
Mechanism: generalises R28 to "module must lie within a rectangle" via a new shape-curve computation; out-performs R28 when specialised to preplaced modules.
Steal / why not: same verdict as R28.

**R30. Young, Wong, "Slicing floorplans with boundary constraint"** [full]. ASP-DAC 1999, DOI 10.1109/ASPDAC.1999.759699.
Mechanism (details in §2.3): a module is on the right (left) boundary iff it is in the right (left) subtree of every ancestor labelled `*`, and on the top (bottom) boundary iff in the right (left) subtree of every `+` ancestor; checked in one right-to-left scan of the Polish expression with a stack of four bits; violations fixed by swapping the module with the nearest module that satisfies the constraint; a distance penalty handles what the fix cannot. ami33 with 12 constrained modules: 1.3–1.9 % dead space in 4–4.6 s; ami49 with 16: 1.5–4.7 % in 36–40 s (300 MHz).
Steal / why not: the "check, fix by swap with the nearest satisfying module, penalise the rest" pattern is what PARSAC (round-1 P2) does in B*-trees and what our SP source/sink invariant does by construction. Confirms the design; nothing to add.

**R31. He, Ji, Li, "An iterative merging placement algorithm for the fixed-outline floorplanning"** [full, skimmed]. arXiv 1401.3172, 2014; journal version Ji, He, Jin, Lan, Li, "An iterative merging algorithm for soft rectangle packing and its extension for application of fixed-outline floorplanning of soft modules", Computers & OR 86, 2017 [abs].
Mechanism: merge the two smallest-area modules repeatedly into a slicing tree, then fix cut directions top-down from the composite aspect ratio; proves zero-dead-space feasibility under a condition weaker than ZDS; extended with WL.
Steal / why not: all-soft, no preplaced, zero whitespace; the "Provocation" of round 2 (dissect around preplaced blocks, then shrink) is this plus obstacles, and the guillotine measurement (M5) already says the obstacle version fragments.

**R32. Hoo, Kanesan, Ramiah, "Enumeration technique in VLSI fixed-outline floorplanning"** [abs]. IET CDS 8(1), 2014, DOI 10.1049/iet-cds.2013.0003.
Mechanism: enumerative floorplanner with dynamic-programming enumerative clustering and a local refinement.
Steal / why not: no constraints of our kind; low relevance.

### 1.5 Multilevel floorplanning (direction 5)

**R33. Lee, Chang, Hsu, Yang, "Multilevel floorplanning/placement for large-scale modules using B*-trees"** (MB*-tree) [full, local `docs/mlbstar-dac2003.pdf`]. DAC 2003, DOI 10.1109/DAC.2003.1219130; journal version IEEE TCAD 26(8), 2007, DOI 10.1109/TCAD.2007.891368 [abs].
Mechanism (details in §2.4): cluster pairs by φ = α·ŝ_ij + β/d̂_ij (ŝ = normalised dead space of the pair, d_ij = nets between the pair divided by the number of primitive modules in it); record the pair's horizontal/vertical relation as a B*-subtree; decluster level by level and refine each level by SA (rotate/move/swap) on area + wirelength; soft modules shaped by Lagrangian relaxation at clustering. Table 1 at 196 modules (area objective): MB*-tree 2.70 % dead space in 2.6 min, B*-tree 6.60 % in 4.7 min, O-tree 9.86 % in 179 min, SP 13.00 % in 309 min; the runtimes are not matched, so only the B*-tree row is a near-equal-time comparison. On a 189-module industrial design, WL-only optimisation: MB*-tree 56 631 vs B*-tree 81 344 vs SP 113 216 (mm), but at 346 min against 12.6 min (B*-tree) and 31 min (SP), so the WL ratios say nothing about equal-time quality.
Steal / why not: the paper's claim is scalability of area utilisation; the one time-matched data point (B*-tree 6.6 % vs 2.7 % at 3–5 min, area objective) is consistent with M2 (flat search at n ≈ 200 unconverged in bounded time) but is not WL evidence. Cost: a second representation and a declustering schedule; interior preplaced blocks must be fixed clusters. Ranked as experiment X6 because the region-scale gap (M6) needs coordinated multi-block moves and clusters are the published way to get them.

**R34. Chen, Wang, "Multilevel large-scale modules floorplanning/placement with improved neighborhood exchange in simulated annealing"** (MBNE) [full]. InTech, Simulated Annealing: Theory with Applications, 2010, DOI 10.5772/10037.
Mechanism (details in §2.4): replaces the MB*-tree refinement SA by Goto's ε-neighbourhood / λ-exchange: pick A, pick B on a common net, pick a random B' within ε B*-tree branches of B, continue to λ modules, try all λ!·2^λ permutations (with rotations) and keep the best; a "null module" lets the exchange act as a move. Results: 70 % of MB*-tree's runtime at equal or lower dead space (2.34/2.11/2.89/2.32 % vs 2.32/2.62/3.18/3.84 % on industry, n100, ami49_20, ami49_60) and lower WL (53 723 vs 55 971 on industry).
Steal / why not: directed, net-driven multi-block exchange as an SA move; λ = 3–4 gives 48–384 candidate permutations per move, affordable at our 4 µs pack. Feeds experiment X4.

**R35. Hu, Chen, Wang, "Fast multilevel floorplanning for large scale modules"** [abs]. ISCAS 2004, DOI 10.1109/ISCAS.2004.1329498.
Mechanism: top-down multilevel genetic floorplanning with clustering only.
Steal / why not: weaker evidence than R33/R34.

**R36. Kahng, Varadarajan, Wang, "Hier-RTLMP: a hierarchical automatic macro placer for large-scale complex IP blocks"** [full, SA sections]. IEEE TCAD 43(5), 2024, DOI 10.1109/TCAD.2023.3346284, arXiv 2304.11761.
Mechanism (details in §2.5): SA over sequence pairs per cluster with four moves at probabilities 0.3/0.3/0.3/0.1 (swap in Γ+, swap in Γ−, swap in both, resize or flip); 500 moves per iteration, 200 iterations, initial acceptance probability 0.9; multi-start with 10 threads; cost = α·area + β·WL + γ·p_outline + penalties; sibling clusters of the parent level act as fixed terminals; preplaced macros are placement blockages.
Steal / why not: a production SP-SA uses pure swaps plus a 10 % resize move and best-of-10 multi-start; the "frozen siblings as terminals" device is the exact-objective window of round 3's windowed re-anneal. Nothing beyond that.

### 1.6 Recent work, 2022–2026 (direction 6)

Contest check: the ICCAD 2025 CAD contest had no floorplanning problem (A: hardware Trojan detection, B: multibit flip-flop optimisation, C: incremental placement optimisation); the FloorSet challenge is ICCAD 2026 Problem C (round-1 P3). Works citing FloorSet in OpenAlex as of today: 7; the only relevant one is R26 (two others are swarm-metaheuristic floorplanners in a low-tier journal, the rest are unrelated). "Flora" returned no floorplanning paper (the hits are an FPGA reconfiguration framework and unrelated items). GoodFloorplan is round-1 P23 and is not re-surveyed. "Lagrangian macro legalization" and "DREAMPlace macro legalization" returned nothing usable beyond R40's description of DREAMPlace 2.0 (R44-style constraint-graph LP plus Puzzle search points) and a 2026 matrix-splitting iteration paper for legalization with macros (Liang, Liao, Zhao, DOI 10.1142/S021812662650297X, not read).

**R37. Xue, Lin, Shi, Kai, Xu, Qian, "Escaping local optima in global placement"** (Hybro) [full]. arXiv 2402.18311, 2024.
Mechanism (details in §2.6): run the analytical placer to convergence, perturb the result (shuffle the locations of 50 % of the macros, or of 50 % of all cells, or re-place all macros greedily with a wire mask), re-run, keep the best over N iterations. Against independent restarts with the same number of placements ("Multiple-DMP"): adaptec4 12.87 vs 19.61 (×10⁷ HPWL), bigblue3 24.76 vs 35.49, adaptec1 6.70 vs 8.01; perturbation costs 12–27 % of the runtime.
Steal / why not: direct evidence that perturbing the incumbent beats independent restarts for a converged placer, with 50 % as the tuned perturbation strength. GPU-scale problem, but the restart logic is representation-free. Feeds experiment X3.

**R38. Liu, Jin, Sun, Lu, Zhang, Yuan, "Coarse-to-fine macro placement via evolutionary search and critical macro tuning"** (C2FPlace) [full, skimmed]. arXiv 2609.34452, 2026.
Mechanism (details in §2.7): population of 20 greedy layouts; tournament of 5; offspring by partial rip-up and re-place of a random macro subset, with the rip-up ratio sampled from a high range early and a low range late; elitist replacement of the worst; final critical-macro tuning (macros on net bounding boxes). HPWL −17.8 % vs EGPlace and −17.9 % vs RollPlace on ISPD 2005; 1000-minute budgets.
Steal / why not: the two-phase rip-up ratio and random (not score-guided) subset selection are the knobs for X3; the budgets are 500× ours.

**R39. Shi et al., "Re²MaP: macro placement by recursively prototyping and packing tree-based relocating"** [skimmed]. arXiv 2511.08054, 2025.
Mechanism: multilevel macro grouping, DREAMPlace prototype, ellipse-based analytical placement, then packing-tree relocating of groups and macros by evolutionary search; only a subset of groups is fixed per iteration, the rest deferred.
Steal / why not: "fix some groups, defer the rest, re-prototype" is a multi-pass version of our seed → SA hand-off; depends on an analytical mixed-size prototype we do not have. Not pursued.

**R40. Jiang et al., "MacroAgent: regularity-aware macro legalization with LLM-agent-designed contour algorithms"** [skimmed]. arXiv 2608.24946, ACM TODAES 2026.
Mechanism: clustering, contour generation, template matching, inter-cluster refinement; the contour heuristics are generated offline by an LLM agent; regularity metric = channels between macros. Reports that DREAMPlace 2.0's macro legalizer is Cong–Xie constraint-graph LP adjustment (R44) plus Puzzle-style search points.
Steal / why not: the objective (regularity, standard-cell channels) is not ours; the LLM-designed-heuristic workflow is an authoring tool, not a solver component.

**R41. Mo, Liu, Xu, Wu, "Order matters: unveiling the hidden impact of macro placement sequences via proxy-guided LLM evolution"** (OrderPlace) [skimmed]. arXiv 2606.08904, ICML 2026.
Mechanism: LLM-evolved placement-order policies for a greedy constructive placer, screened by a cheap deterministic proxy; −34 % WL vs WireMask-EA, −14 % vs EGPlace.
Steal / why not: constructive placement failed here (M5); ordering heuristics do not transfer.

**R42. Li, Kai, Hao, Zhong, Tang et al., "CORE: collaborative optimization with reinforcement learning and evolutionary algorithm for floorplanning"** [abs]. NeurIPS 2025, DOI 10.52202/085713-4662.
Mechanism: clustering-based diversified evolutionary search that perturbs layouts directly, coupled with RL over a B*-tree sequential model; MCNC/GSRC.
Steal / why not: population + direct layout perturbation is again the X3 pattern; the RL part is not available to us.

**R43. Yu, Du, Yang, "A deep reinforcement learning floorplanning algorithm based on sequence pairs"** [abs], Applied Sciences 14(7):2905, 2024; and Yu, Du, "VLSI floorplanning algorithm based on reinforcement learning with obstacles" [abs], 2024, DOI 10.1007/978-3-031-50381-8_110.
Mechanism: an agent explores the SP space; the second paper adds obstacles.
Steal / why not: round 1 already concluded RL over SP matches tuned SA at best (P23); no training budget.

**R44. Cong, Xie, "A robust mixed-size legalization and detailed placement algorithm"** (XDP) [abs]. IEEE TCAD 27(8), 2008, DOI 10.1109/TCAD.2008.925792.
Mechanism: two constraint graphs from the global placement; iterative adjustments until both longest paths fit the chip; LP for macro coordinates minimising displacement (and spreading).
Steal / why not: this is the legalizer class that failed on our GP centres (M5); the mixed-size setting has far more whitespace.

**R45. Chen, Chuang, Chang, Chang, "Constraint graph-based macro placement for modern mixed-size circuit designs"** [abs]. ICCAD 2008, DOI 10.1109/ICCAD.2008.4681577.
Mechanism: constraint-graph macro placer on essential edges only; non-compacted placement by LP and region-cost evaluation.
Steal / why not: same remark as R44; the "essential edges" idea is our transitive reduction (M4 tables).

**R46. Du, Zhong, Kai, Tang, Xu et al., "JigsawPlanner: jigsaw-like floorplanner for eliminating whitespace and overlap among complex rectilinear modules"** [abs]. ICCAD 2024, DOI 10.1145/3676536.3676707.
Mechanism: from analytical centres, partition chip and modules into grids and sub-modules, assign by hierarchical Jonker–Volgenant and cellular automata.
Steal / why not: rectilinear modules; not allowed in our lab.

**R47. Xu, Wang, Deng, "The impact of move schemes on simulated annealing performance"** [skimmed]. arXiv 2504.17949, 2025.
Mechanism: theory plus experiments (Lennard-Jones) showing that moving one randomly chosen coordinate per step beats full-coordinate moves at equal proposal variance.
Steal / why not: supports small moves late in the anneal; generic, continuous setting.

### 1.7 Variance, restarts and acceptance rules (direction 7)

Searched and not found: parallel tempering, population annealing or record-to-record travel with evidence on floorplanning or rectangle packing. Generic sources are listed where the mechanism is portable.

**R48. Cicirello, "Variable annealing length and parallelism in simulated annealing"** [full]. SoCS 2017, DOI 10.1609/socs.v8i1.18424.
Mechanism (details in §2.8): Modified Lam schedule (target acceptance rate 0.44 after the first 15 % of the run, exponential decay over the last 35 %, temperature multiplied or divided by 0.999 per step); restarts with lengths 1000·2^r. Findings on a scheduling benchmark: one long run beats fixed-length restarts at end of run; the doubling schedule matches it at the end and dominates any time before; re-annealing from the best-of-run solution is not better than random restarts in the sequential case (not significant), marginally better in parallel late in the run.
Steal / why not: the adaptive schedule is parameter-free except for the run length, which we know (110 s); the re-annealing finding is counter-evidence to incumbent-perturbation restarts (X3) and is why X3 must be gated by a measurement.

**R49. Lam, Delosme, "Performance of a new annealing schedule"** [abs]. DAC 1988, DOI 10.1109/DAC.1988.14775.
Mechanism: statistics-based schedule; 2–24× speedup over schedules of the time on standard-cell placement and TSP; the 0.44 acceptance target used by R48 comes from here.
Steal / why not: see R48.

**R50. Hu, Kahng, Tsao, "Old bachelor acceptance: a new class of non-monotone threshold accepting methods"** [abs]. ORSA/INFORMS J. Computing 7(4):417–425, 1995, DOI 10.1287/ijoc.7.4.417.
Mechanism: threshold acceptance whose threshold is self-tuning and non-monotone (even negative), designed for a prescribed time bound; Dueck–Scheuer threshold accepting is a special case; TSP experiments.
Steal / why not: designed for exactly our setting (fixed time budget, late tail); no floorplanning evidence. Candidate acceptance rule for X5.

**R51. Burke, Bykov, "The late acceptance hill-climbing heuristic"** [abs]. EJOR 258(1):70–78, 2017, DOI 10.1016/j.ejor.2016.07.012.
Mechanism: accept a candidate if it is no worse than the current solution or than the solution L iterations ago (a length-L history list); one parameter.
Steal / why not: one-parameter alternative to a temperature schedule for the post-feasible phase; candidate for X5. No floorplanning evidence.

**R52. Dueck, "New optimization heuristics: the great deluge algorithm and the record-to-record travel"** [sec]. J. Computational Physics 104(1):86–92, 1993, DOI 10.1006/jcph.1993.1010.
Mechanism: accept if the new cost is below the record plus an allowed deviation; round 3 already lists RRT for the post-feasible phase.
Steal / why not: see X5.

**R53. Vorwerk, Kennings, Greene, "Improving simulated annealing-based FPGA placement with directed moves"** [abs]. IEEE TCAD 28(2):179–192, 2009, DOI 10.1109/TCAD.2008.2009167; thesis "On the use of directed moves for placement in VLSI CAD", Waterloo 2009 (UWSpace, blocked by a bot wall).
Mechanism: intelligent cell selection and target selection (median/WL-driven, timing-driven) interspersed with random moves during the anneal; at equal effort, better wirelength and critical path, faster convergence and "less statistical variability".
Steal / why not: the variability claim is the one that matters for M3. The SP analogue is a move that inserts a block near the SP position of its net-weighted median (our oracle) or swaps it with the block occupying that spot (round 3, cluster VI). Feeds X4.

### 1.8 Boundary, corner and same-size constraints beyond round-1 P6 (direction 8)

**R54. Ma, Dong, Hong, Cai, Cheng, Gu, "VLSI floorplanning with boundary constraints based on corner block list"** [abs]. ASP-DAC 2001, DOI 10.1145/370155.370521.
Mechanism: necessary and sufficient characterisation of boundary modules in a CBL; linear-time check and fix inside SA.
Steal / why not: same pattern as R30; our SP invariant already gives boundary by construction.

**R55. Yang, Dong, Hong, Ma, "A two-stage incremental floorplanning algorithm with boundary constraints"** [abs]. APCCAS 2006, DOI 10.1109/APCCAS.2006.342140.
Mechanism: a slack-push graph (dead space between adjacent modules and the slack push distance per module) and a module-exchange graph (modules with the same size character); stage 1 satisfies boundary constraints with the graphs, stage 2 optimises with exchanges of same-size modules.
Steal / why not: the module-exchange graph is the only published MIB-shaped idea found: identical-size blocks can be swapped at zero geometric cost, so "swap two members of an MIB group" (and more generally two equal-shape blocks) is a WL-only move that never breaks legality. Cheap addition to the move set.

**R56. Lin, Chen, Wang, "Modern floorplanning with boundary and fixed-outline constraints via genetic clustering algorithm"** [abs]. J. Circuits, Systems and Computers 15(1), 2006, DOI 10.1142/S0218126606002940; and **R57. Lin, Chen, Wang, Ho, "Modern floorplanning with abutment and fixed-outline constraints"** [abs], ISCAS 2005, DOI 10.1109/ISCAS.2005.1466060.
Mechanism: slicing floorplans; boundary-constrained modules are clustered into four constrained sub-floorplans and recombined; the floorplan is then gradually fitted into the fixed outline; evolutionary search.
Steal / why not: "cluster the boundary ring, then fit the interior" is round 3's "exact DP for the boundary ring" idea; the slicing representation does not survive our interior preplaced blocks.

Symmetry constraints from analog placement (symmetry islands in B*-trees, symmetric-feasible sequence pairs) were skipped: MIB requires identical dimensions, not mirrored positions, and the identical-dimension part is already solved by the shared shape index (round-1 D2).

---

## 2. Details extracted from the full reads

### 2.1 Umetani et al. 2009 (R9) and Umetani–Murakami 2022 (R10): the GLS overlap-minimisation kernel

- Penalty: f_ij = min(ρ_H, ρ_V), the smaller directional penetration depth; for axis-aligned rectangles this is min(x-overlap length, y-overlap length). Weighted objective F̃ = Σ w_ij f_ij.
- SEARCHNEIGHBOR(P_k): translate P_k along the horizontal line through its position, choosing the t that minimises F̃_k(v_k + t·e_x) over the whole feasible line (inside the container), ties broken by the nearest position; then along the vertical line; alternate until neither direction improves. Each 1-D function is the lower envelope of piecewise-linear pieces; all breakpoints from the O(n) interacting pieces are merged in O(q_k log n). For rectangles the pieces are triangles (ρ_H) and constants (ρ_V) over the overlap interval, so the whole line has O(n) breakpoints.
- IMPROVE: first-admissible local search over blocks with Fast Local Search: every block is "active" at the start; a block whose sub-neighbourhood yields no improvement is deactivated; when P_k moves, every block that overlapped P_k before or after the move is reactivated. Stops when the overlap is zero or no active block remains.
- MINIMIZEOVERLAP (GLS): w_ij = 1.0 initially; at each local optimum, for every overlapping pair w_ij ← w_ij + f_ij / max_kl f_kl; resume from the same solution; stop after max_iter = 200 (R9) / k_max = 200 (R10) consecutive IMPROVE calls without improving the unweighted best.
- Outer loop FITS: if the current placement is feasible, L ← (1 − r_dec)·L with r_dec = 0.02 and the protruding pieces are relocated at random; otherwise L ← (1 + r_inc)·L with r_inc = 0.005; keep L_LB ≤ L < L*. Runs: 10 × 1200 s on a Xeon 2.8 GHz; FITS best on 4 and average-best on 5 of 15 ESICUP instances; its efficiency is within about 1 % of ILSQN (R11) and 2DNest (R12) elsewhere.
- R10 adds integer coordinates (t ∈ Z) and reports the same parameters on a single thread in C (3.2 GHz i7); construction by NFDH plus bottom-left compaction implemented with the same line search (accept any t < 0 without overlap).

What transfers: our legalizer already has the weight scheme; the exact line search and the activation rule are the two mechanics to check for in `solver/gls.h`. With a WL term added to F̃_k (piecewise linear in the moved coordinate with one breakpoint per connection), the same line search minimises WL + Σ w_ij f_ij exactly (R14's objective).

### 2.2 Imamichi et al. 2009 (R11): swap-and-separate

- Perturbation SwapTwoPolygons: remove P_i; place a copy of P_j at P_i's place so P_j cannot stay; move P_j to the minimum-overlap position found by FindBestPosition (evaluate the total penalty at every vertex and edge intersection of the no-fit polygons of all other pieces and the container; Lemma 1: if a zero-overlap position exists it is found); remove the copy; place P_i the same way. A cheaper variant FindBestPosition2 samples K vertices only.
- Local search Separate: quasi-Newton on Σ penetration depths plus protrusion, all pieces moving at once, no rotation.
- ILS: accept the local optimum if it lowers the total overlap; stop after N_mo consecutive failures. Initial solution: pieces in descending area order, each at its best position, bottom-left tie-break.

What transfers: the "two blocks, each re-placed by an exact oracle, then a continuous separation" step is a coordinated two-block move with a legality oracle; our SP oracle plus GLS gives the same step in our representation.

### 2.3 Young–Wong 1998 (R28) and 1999 (R30): slicing with preplaced and boundary modules

- R28: each supermodule X carrying a preplaced module has a reference point ref(X) (the lower-left corner of the inherited preplaced module) and four numbers bottom(X), top(X), left(X), right(X) delimiting where that point may sit inside X. Combining two supermodules that both carry preplaced modules requires y_b > y_a and y_b − y_a ≥ top(A) + bottom(B) (vertical case) and yields closed-form updates of h, w and the four slacks (two cases each for `+` and `*`). Shape curves become step functions for such supermodules; only points with equal x (or y) are combined, which the authors report costs nothing in quality. Moves: M1.1 (swap two adjacent operands, one preplaced), M1.2 (one soft), M2 (complement a chain), M3.1/M3.2 (swap operand and operator, in the reduced or full expression), M4 (combine two operands to move a module out of a congested area). M1.1/M3.1 fire with higher probability while the preplaced relative order is wrong. Cost A + λW + γD, γ large; D → 0 as annealing proceeds. Schedule: T₀ at ~100 % acceptance, constant ratio, 20–30·n iterations per temperature.
- R30: Lemma 1 (boundary membership ⇔ subtree side of all `*` or `+` ancestors); Algorithm Check-Boundary-Constraints is a single right-to-left scan with a stack whose elements carry left/right/above/below bits; fix = exchange the violating module with the closest module in the expression whose position satisfies the constraint; penalty D = distance from the required boundary. Schedule: ratio 0.9, 20·n iterations per temperature.

What transfers: nothing to the solver (see R28's verdict); the boundary check-and-fix pattern confirms the current SP invariant design.

### 2.4 MB*-tree (R33) and MBNE (R34): clustering and λ-exchange

- Clustering metric φ({m_i,m_j}) = α·ŝ_ij + β/d̂_ij with ŝ_ij the normalised dead space of the pair's bounding rectangle and d_ij = c_ij/(n_i + n_j) (nets between the two over the number of primitive modules in them). Pairs are merged until one cluster (or a threshold count) remains; the horizontal/vertical relation of each merge is stored as a B*-subtree.
- Declustering expands one node per step in O(1) and refines the current level by SA with rotate / move / swap, cost γ·ŝ + δ·ŵ. Soft modules: Lagrangian relaxation of the cluster's constraint graphs with closed-form widths w_p = sqrt(Σ_v λ A_p / Σ_h λ) clipped to [L_p, U_p], multipliers projected by a max-flow test.
- MBNE's refinement: ε-neighbourhood = modules within ε B*-tree branches; λ-exchange chain built along nets (A, then B on a common net, then a random module near B, …), all λ!·2^λ permutations evaluated, best kept; a zero-size "null module" inserted as a child turns an exchange into a move. Reported 30 % runtime saving at equal or better area/WL.

What transfers: (i) one time-matched data point that flat B*-tree search at 196 modules is unconverged in a few minutes on an area objective (M2 agrees for WL on our instances); (ii) the net-driven λ-exchange as a directed multi-block move (X4); (iii) connectivity-density clustering as the way to form super-blocks for coordinated moves (X6).

### 2.5 Hier-RTLMP (R36): production SP-SA settings

- Four moves, probabilities 0.3/0.3/0.3/0.1: swap two in Γ+, swap two in Γ−, swap two in both, resize (clusters) or flip (macros). 500 moves per iteration × 200 iterations, initial acceptance probability 0.9, minimum temperature 1e-10 (printed as "1e10" in the arXiv text; (uncertain)). Multi-start with 10 threads; the best run is kept. Cost = α·Area + β·WL + γ·p_outline + ζ·p_bias + η·p_blockage + θ·p_guidance + λ·p_notch, each term normalised by its initial value. Siblings at the parent level are fixed terminals; preplaced macros are placement blockages.

What transfers: nothing we lack; confirms that best-of-k multi-start is the industrial answer to SA variance (our 2×55 result shows the median cost of that answer).

### 2.6 Hybro (R37): perturb-and-reoptimise beats independent restarts

- Loop: run the gradient placer to convergence; perturb (Shuffle: random 50 % of macros exchange locations; Shuffle-all: 50 % of all cells; WireMask: re-place every macro greedily to the grid cell with the least HPWL increase); repeat N times; return the best. 50 % was found to be the right strength ("too high behaves like restarts, too low does not escape").
- ISPD 2005, 3 runs each, HPWL ×10⁷: Multiple-DMP (independent restarts) vs Hybro-WireMask: adaptec1 8.01 → 6.70, adaptec3 14.56 → 13.51, adaptec4 19.61 → 12.87, bigblue3 35.49 → 24.76; adaptec2 and bigblue1 within noise. Perturbation = 12 % (Shuffle) or 27 % (WireMask) of the runtime.

What transfers: the restart policy for X3 (perturb half of the state, re-optimise from there) and its tuned strength.

### 2.7 C2FPlace (R38): rip-up ratio schedule

- Population 20 (greedy WireMask layouts with random tie-breaks); tournament of 5 picks the parent; offspring = rip up a random subset of macros (ratio sampled from a high range in phase 1, a low range in phase 2) and re-insert them greedily in connectivity/area order; offspring replaces the worst only if better; final critical-macro tuning moves bounding-box-defining macros within their swept-free region and then at fine offsets, accepting only improvements.

What transfers: random subset selection and the two-phase ratio are the two knobs to copy into X3; the elitist population is optional (our K = 2–3 trajectories, round 3).

### 2.8 Cicirello 2017 (R48): adaptive schedule and restart evidence

- Modified Lam: track an exponential moving average of the acceptance rate (weight 1/500); target 0.44 + 0.56·560^(−i/(0.15·MaxEvals)) for the first 15 % of the run, 0.44 for the next 50 %, 0.44·440^(−(i/MaxEvals − 0.65)/0.35) for the last 35 %; T ← 0.999·T when the acceptance rate is above target, T ← T/0.999 when below; T₀ = 0.5 on a normalised cost.
- Findings: fixed-length restarts tuned to the full budget win at the end of the run; the doubling schedule (1000·2^r evaluations) matches them at the end (differences not significant from 48 s of 60 s) and dominates earlier; re-annealing from the best-of-run instead of a random start: no significant difference sequentially, marginal gain in parallel during the last 15–20 s.

What transfers: a parameter-free schedule for the post-feasible phase (X5) and the gating condition for X3 (re-annealing is not automatically better than fresh restarts; Hybro says perturbation strength is what decides).

---

## 3. Cross-check against the measured dead ends

| Idea from the literature | Measured on this project | Status |
|---|---|---|
| Constraint-graph / LP repair of an overlapping placement (R7, R8, R16, R20, R44, R45) | GP centres legal on 0/56 hidden4 seeds; 250–340 ms per repair iteration; R20 itself reports +19.5 % WL after repair | closed |
| Partition or slicing legalisation around obstacles (R17, R18, R19, R28, R29, R31, R56) | guillotine around 24 preplaced blocks: ~90 regions, ~60 infeasible; 6 % of free area unusable | closed at our preplaced density |
| Constructive packing with a WL reference (R24, R31, R41) | MER greedy: 4–6 % area unplaced, 1.96–3.51× WL | closed |
| Single-block rip-up / insertion (R5 IAR as LNS) | 4–11 % of blocks improvable, total gain 0.03–1.1 % | closed as LNS; untested as an in-anneal move (X4) |
| Elite pool / consensus | variance is late-tail; continuations from first-legal spread like independent runs | closed |
| Best-of-k independent restarts (R36 multi-start, R48 FAL) | 3×36 s / 2×55 s: max/median 1.26 → 1.15/1.14, median flat | partial; replaced by incumbent perturbation (X3) and go-with-the-winners (round 3) |
| Feasibility homotopy ("fix the budget, legalise, tighten": R9 FITS, R13, R21 area growing, R23 instance augmentation) | area scaling s = 0.92 → GLS absorber → true anneal: 0.845–0.88× control, legal 47/48 | **open and winning**; the literature pattern says: many small steps, each re-legalised (X2) |
| WL inside the legalizer objective (R14, R22 superiorisation) | hand-off 2.1× with penalties only | open (X2) |
| Exact-geometry evaluation of SP neighbours (R1, R2, R3) | packed WL anti-correlates with LP WL (ρ −0.19); LP post-pass 1.8–4.2 % | open (X1) |
| Directed / net-driven moves (R34, R53, R55) | not measured; round 3 lists them as cheap knobs | open (X4) |
| Clustering into super-blocks (R33, R34) | not measured; displacement pre-test says the gap is region-scale | open (X6) |

---

## 4. Access notes and unverified items

- Imahori 2003/2005 (R1, R2): the Kyoto page with the author PDFs is gone (404); the KURENAI record of the doctoral thesis holds a 2-page scanned summary only; Springer/Elsevier are paywalled. Neighbourhood definitions and the critical-path restriction are (uncertain).
- Egeblad 2007/2008 (R12): CORE and KU research-profile mirrors are behind Cloudflare; λ/utility settings not extracted.
- Faroe et al. final placement (R14): the DIKU technical-report server (di.ku.dk) refused connections during this session; the mechanics are reconstructed from R9 plus the abstract.
- Vorwerk thesis (R53): UWSpace bot wall; the TCAD abstract is the only source; the move definitions and probabilities are (uncertain).
- Hier-RTLMP minimum temperature "1e10" in the arXiv text is presumably a typeset 1e-10 (uncertain).
- Chen–Yoshimura (R5), Tang–Tian–Wong (R3), UFO (R16), SKB-tree (R19), Chan–Hsu–Lin (R17), He–Ma (R20), Lu (R21), Ma (R54), Yang (R55), Lin (R56, R57), Lam–Delosme (R49), OBA (R50), LAHC (R51): abstract level; numbers quoted are from the abstracts.
- Semantic Scholar rate-limited most abstract requests; OpenAlex was used instead. Search listings are in scratchpad `searches2/`.

---

## 5. Tables

### 5a. All papers in this round

Access: F = full text, S = sections, A = abstract/record, C = cited description only.

| # | Paper | Venue, year | DOI / URL | Direction | Access | One-line verdict |
|---|---|---|---|---|---|---|
| R1 | Imahori, Yagiura, Ibaraki, rectangle packing with general spatial costs | Math. Prog. 2003 | 10.1007/s10107-003-0427-1 | 1 | A | SP evaluated by cost-optimal coordinates (DP); the published form of X1 |
| R2 | Imahori, Yagiura, Ibaraki, improved local search for the same | EJOR 2005 | 10.1016/j.ejor.2004.02.020 | 1 | A | fast neighbourhood evaluation on top of R1 |
| R3 | Tang, Tian, Wong, minimizing wire length in floorplanning | TCAD 2006 | 10.1109/TCAD.2005.858266 | 1 | A | min-cost-flow WL optimum for a fixed topology; 4.2 % on compact floorplans; handles preplaced/boundary/soft |
| R4 | Tang, Wong, alignment and performance constraints | DAC 2002 | 10.1109/DAC.2002.1012740 | 1, 8 | A | fixed-frame SP with weighted-LCS feasibility |
| R5 | Chen, Yoshimura, block-position enumeration (IAR) | TCAD 2008 / ISPD 2007 | 10.1109/TCAD.2008.917968 | 1 | A | insertion-after-remove = our oracle; weak as LNS, untested as SA move |
| R6 | Silvanus, Vygen, few sequence pairs suffice | arXiv 2017 | 1708.09779 | 1 | A | theory on SP redundancy |
| R7 | Vorwerk, Kennings, Anjos, floorplan repair | Eng. Opt. 2008 | 10.1080/03052150801901475 | 1, 6 | A | constraint-graph + LP repair; class closed by M5 |
| R8 | Moffitt et al., Floorist | DAC 2006 | 10.1109/DAC.2006.229405 | 1, 6 | A | constraint-driven repair; class closed by M5 |
| R9 | Umetani et al., GLS for overlap minimisation (FITS) | ITOR 2009 | 10.1111/j.1475-3995.2009.00707.x | 2 | F | exact 1-D line search, FLS activation, weight update, shrink/extend loop |
| R10 | Umetani, Murakami, coordinate descent heuristics | EJOR 2022 | arXiv 2104.04525 | 2 | S | same kernel on integer coordinates; parameters confirmed |
| R11 | Imamichi, Yagiura, Nagamochi, ILSQN | Discrete Opt. 2009 | 10.1016/j.disopt.2009.04.002 | 2 | S (thesis) | swap two pieces with an exact best-position oracle, then quasi-Newton separation |
| R12 | Egeblad, Nielsen, Odgaard, fast neighborhood search for nesting | EJOR 2007 | 10.1016/j.ejor.2005.11.063 | 2 | A, C | GLS with area penalty and exact 1-D translation |
| R13 | Faroe, Pisinger, Zachariasen, GLS for 3D bin packing | IJOC 2003 | 10.1287/ijoc.15.3.267.16080 | 2 | A | fix bins, legalise, remove a bin, repeat |
| R14 | Faroe, Pisinger, Zachariasen, GLS for final placement | ICCAD 2001 / J. Heuristics 2003 | 10.1023/A:1023721408655 | 2 | A | GLS with WL in the objective; up to 20 % on large layouts |
| R15 | Egeblad, Pisinger, knapsack packing | COR 2009 | 10.1016/j.cor.2007.12.004 | 2 | A | SP/sequence-triple SA |
| R16 | Lin, Hung, UFO | TCAD 2011 | 10.1109/TCAD.2011.2114531 | 3 | A | convex push-pull + SOCP legalisation without WL; preplaced supported |
| R17 | Chan, Hsu, Lin, flexible fixed-outline for mixed-size | ASP-DAC 2013 | 10.1109/ASPDAC.2013.6509635 | 3 | A | analytical + partition slicing legalisation |
| R18 | Lin, Wu, F-FM | TCAD 2014 | 10.1109/TCAD.2014.2351571 | 3 | A | R17 with voltage islands |
| R19 | Lin, Hung, SKB-tree | TVLSI 2012 | 10.1109/TVLSI.2011.2104983 | 3 | A | region-allocating tree, WL-only SA; needs soft-dominated instances |
| R20 | He, Ma et al., fixed blocks in incremental fixed-outline FP | ICCCAS 2010 | 10.1109/ICCCAS.2010.5581854 | 3 | A | TCG repair with virtual nodes; +19.5 % WL |
| R21 | Lu, Huang, Zhu, fixed-outline FP for rectilinear soft modules | ISEDA 2025 | 10.1109/ISEDA65950.2025.11100867 | 3 | A | area-growing global stage; negotiation legalisation with preplaced |
| R22 | Yu, Censor, Luo, feasibility seeking + superiorization | TCAD 2024 | 10.1109/TCAD.2024.3408106 | 3 | S | alternating projections legalizer with WL perturbations |
| R23 | Liu, Dong, Hong, Kajitani, instance augmentation | ISCAS 2005 | 10.1109/ISCAS.2005.1464979 | 3 | A | homotopy in module count |
| R24 | Wei, Dong, Hong, Wu, less-flexibility-first | ISCAS 2006 | 10.1109/ISCAS.2006.1693903 | 3 | A | constructive; class closed |
| R25 | Huang, Lin, Zhu, Chen, excessive length penalty SA | IEEE Access 2020 | 10.1109/ACCESS.2020.2980135 | 3 | A | outline penalty shaping; feasibility already solved here |
| R26 | Zhu et al., diffusion + HGNN initial floorplans | TCAD 2026 | 10.1109/TCAD.2026.3681657 | 3, 6 | A | analytical floorplanners are seed-sensitive; needs training |
| R27 | Sun et al., mixed-variable optimization | arXiv 2024 | 2401.15317 | 3 | S | EDA + conjugate subgradient; no legality guarantee |
| R28 | Young, Wong, slicing with pre-placed modules | ICCAD 1998 | 10.1145/288548.288622 | 4 | F | reference points + step shape curves; 3–4 preplaced modules |
| R29 | Young, Wong, slicing with range constraint | ISPD 1999 | 10.1145/299996.300033 | 4 | A | generalises R28 |
| R30 | Young, Wong, slicing with boundary constraint | ASP-DAC 1999 | 10.1109/ASPDAC.1999.759699 | 4, 8 | F | check–fix–penalise boundary pattern |
| R31 | He, Ji, Li, iterative merging placement (+ Ji et al. 2017) | arXiv 2014 / COR 2017 | 1401.3172 | 4 | S | zero-dead-space slicing for all-soft instances |
| R32 | Hoo, Kanesan, Ramiah, enumeration technique | IET CDS 2014 | 10.1049/iet-cds.2013.0003 | 4 | A | enumerative clustering; low relevance |
| R33 | Lee, Chang, Hsu, Yang, MB*-tree | DAC 2003 / TCAD 2007 | 10.1109/DAC.2003.1219130 | 5 | F | flat SP/B*-tree unconverged at 196 modules; clustering metric |
| R34 | Chen, Wang, MBNE | InTech 2010 | 10.5772/10037 | 5, 7 | F | net-driven λ-exchange; 30 % less runtime |
| R35 | Hu, Chen, Wang, fast multilevel floorplanning | ISCAS 2004 | 10.1109/ISCAS.2004.1329498 | 5 | A | genetic multilevel |
| R36 | Kahng, Varadarajan, Wang, Hier-RTLMP | TCAD 2024 | 10.1109/TCAD.2023.3346284 | 5, 6 | S | production SP-SA settings; frozen siblings as terminals |
| R37 | Xue et al., Hybro | arXiv 2024 | 2402.18311 | 6, 7 | F | perturb 50 % of the incumbent and re-optimise beats restarts |
| R38 | Liu et al., C2FPlace | arXiv 2026 | 2609.34452 | 6, 7 | S | random-subset rip-up with a high→low ratio schedule |
| R39 | Shi et al., Re²MaP | arXiv 2025 | 2511.08054 | 6 | S | recursive prototype + packing-tree relocation |
| R40 | Jiang et al., MacroAgent | arXiv / TODAES 2026 | 2608.24946 | 6 | S | LLM-designed legalisation heuristics for regularity |
| R41 | Mo et al., OrderPlace | ICML 2026 | arXiv 2606.08904 | 6 | S | LLM-evolved placement orders |
| R42 | Li et al., CORE | NeurIPS 2025 | 10.52202/085713-4662 | 6 | A | EA layout perturbation + RL |
| R43 | Yu, Du (+Yang), RL over sequence pairs (+ obstacles) | Appl. Sci. 2024 | 10.3390/app14072905 | 6 | A | RL ≈ SA; no training budget |
| R44 | Cong, Xie, XDP mixed-size legalization | TCAD 2008 | 10.1109/TCAD.2008.925792 | 6 | A | constraint graphs + LP; class closed |
| R45 | Chen, Chuang, Chang, Chang, constraint-graph macro placement | ICCAD 2008 | 10.1109/ICCAD.2008.4681577 | 6 | A | essential edges + LP |
| R46 | Du et al., JigsawPlanner | ICCAD 2024 | 10.1145/3676536.3676707 | 6 | A | rectilinear grid assignment |
| R47 | Xu, Wang, Deng, move schemes in SA | arXiv 2025 | 2504.17949 | 7 | S | single-coordinate moves win |
| R48 | Cicirello, variable annealing length | SoCS 2017 | 10.1609/socs.v8i1.18424 | 7 | F | Modified Lam schedule; long run > restarts; re-annealing ≈ restart |
| R49 | Lam, Delosme, new annealing schedule | DAC 1988 | 10.1109/DAC.1988.14775 | 7 | A | 0.44 acceptance target; 2–24× speedup on placement |
| R50 | Hu, Kahng, Tsao, old bachelor acceptance | ORSA JOC 1995 | 10.1287/ijoc.7.4.417 | 7 | A | self-tuning non-monotone threshold for fixed time budgets |
| R51 | Burke, Bykov, late acceptance hill-climbing | EJOR 2017 | 10.1016/j.ejor.2016.07.012 | 7 | A | one-parameter acceptance rule |
| R52 | Dueck, great deluge and record-to-record travel | J. Comp. Phys. 1993 | 10.1006/jcph.1993.1010 | 7 | C | bounded-drawdown acceptance |
| R53 | Vorwerk, Kennings, Greene, directed moves | TCAD 2009 | 10.1109/TCAD.2008.2009167 | 7 | A | directed moves cut variability and improve QoR at equal effort |
| R54 | Ma et al., boundary constraints with CBL | ASP-DAC 2001 | 10.1145/370155.370521 | 8 | A | check-and-fix in CBL |
| R55 | Yang, Dong, Hong, Ma, two-stage incremental with boundary | APCCAS 2006 | 10.1109/APCCAS.2006.342140 | 8 | A | module-exchange graph of same-size modules |
| R56 | Lin, Chen, Wang, boundary + fixed-outline via genetic clustering | JCSC 2006 | 10.1142/S0218126606002940 | 8 | A | cluster the boundary ring, fit the interior |
| R57 | Lin, Chen, Wang, Ho, abutment + fixed-outline | ISCAS 2005 | 10.1109/ISCAS.2005.1466060 | 8 | A | slicing feasibility properties |

### 5b. Ranked experiments for our solver (≤ 6)

All gains are stated against the numbers in §0; where no number exists the entry says so. Each experiment is phrased to sit on top of the round-3 homotopy flow, which is the current best (0.845–0.88× control).

**X1. Exact-geometry acceptance in the post-feasible anneal (LP/MCF as the evaluator).** Replace packed WL by the per-axis LP WL in the acceptance test once the state is feasible (or for the top-k candidates per temperature step), using a warm-started network simplex so a local SP change costs few pivots. Also allow hinge terms |x_i − anchor_i| so the round-3 proximal anchor and R3's constraint set live in the same LP.
Evidence: M4 (ρ −0.19; the LP-best topology had the worst packed WL in 14/20 cases; post-pass 1.8–4.2 %); R3 (min-cost-flow WL optimum, 4.2 % even on compact floorplans, < 0.4 s on MCNC, preplaced/boundary/soft inside the model); R1/R2 (the SP literature with spatial costs searches on cost-optimal coordinates, not packings).
Expected gain: at least the post-pass 2–4 % moved into the search signal; the real target is part of the 1.2–1.6× gap that M4 attributes to ranking error. Unquantified beyond that; gate on a 30-seed A/B at equal wall time with the solve time measured first (estimate 0.1–0.5 ms per axis).
Cost: network simplex (already planned, ~300–600 lines) plus warm start.

**X2. FITS-style fine-step homotopy with WL inside the legalizer objective.** Replace the single ramp 0.92 → 1 by many small steps (r_inc ≈ 0.5–1 % of area scale per step, as in R9/R10's 0.5 % extension step), each followed by the GLS absorber; run the absorber on WL + Σ w_ij f_ij with the exact alternating 1-D line search (WL is piecewise linear along a coordinate line), FLS activation, and the w_ij ← w_ij + f_ij/max f update, so legality is regained without discarding the basin's WL.
Evidence: M6 (homotopy 0.845–0.88×, but the ramp itself costs 8–40 % WL and stalls at s ≥ 0.98 on 7/14 cases; the GLS absorber lifts reach to 47/48); R9/R10/R13 (every GLS packer tightens the budget in small steps and re-legalises each time; 2 % shrink, 0.5 % extend); R14 (GLS with WL in the objective improved large placements by up to 20 %); R22 (interleaving WL descent with feasibility steps is a published pattern).
Expected gain: round 3 measured that the true anneal recovers the ramp penalty when given +60 s ("ramp then GLS" 0.878× and "GLS straight from s = 0.92" 0.879× of control, same ratio [meas]), so the ramp ratio itself (final/relaxed median 1.24 [meas]) is not the target. The target is budget: a WL-aware absorber hands the true anneal a state closer to the relaxed basin's WL, so more of the fixed 110 s goes to improvement instead of recovery. Measure final WL at 110 s against the Phase F flow, plus the hand-off ratio of the GLS legalizer used cold (2.1× today).
Cost: changes inside `solver/gls.h` (objective and line search) and the homotopy driver; no new data structure.

**X3. Incumbent-perturbation restarts instead of independent restarts.** At the point where the tail stops improving (round 3 estimates 35–50 s), perturb the incumbent by ripping up a random 20–50 % of the movable blocks (SP positions and shapes), re-insert them (oracle or random), re-anneal at low temperature, keep the better of parent and child; ratio high first, low later. Couples with round 3's go-with-the-winners scheduling.
Evidence: M3 (continuations from 50 s spread only 1.01–1.14; 2×55 s cut max/median from 1.26 to 1.14 at a flat median); R37 (perturbing 50 % of the incumbent beats independent restarts on 4/6 ISPD 2005 cases, e.g. 12.87 vs 19.61 on adaptec4); R38 (random-subset rip-up with a high→low ratio beats fixed-ratio and score-guided selection by 17.8 %); counter-evidence R48 (re-annealing the best-of-run gave no significant gain over random restarts in the sequential case), so the experiment is gated: measure Spearman(signal at t, final WL) first, as round 3 already plans.
Expected gain: the tail spread (max/median up to 1.3) down to about 1.1 at a flat or slightly better median; not a median lever.
Cost: small (a perturbation operator and a trajectory loop).

**X4. Directed moves in the SA move set.** With probability 5–10 %: (a) net-median insertion — pick a block, compute its net-weighted median point, insert it at the SP position of the block nearest that point (our exact oracle, 130–230 µs per block [meas], or a cheaper cell guess); (b) λ-exchange along nets (R34: chain of 3–4 blocks connected by nets, try all permutations, keep the best; with MIB/equal-shape pairs the exchange is geometry-free, R55); (c) swap-and-separate (R11) when the legalizer is active.
Evidence: R53 (directed moves interspersed with random ones give better QoR and "less statistical variability" at equal effort); R34 (30 % less runtime at equal or better quality vs MB*-tree SA); R55 (same-size exchange); M5 says post-hoc single-block relocation is weak, which only rules out the LNS form.
Expected gain: unknown; the variability claim is the reason to test it, scored on max/median over 30 seeds.
Cost: small; the oracle exists.

**X5. Acceptance rule for the post-feasible phase.** Try, in order of evidence: Modified Lam adaptive temperature (R48/R49: target acceptance 0.44, parameter-free given the known budget), record-to-record travel with a drawdown bound (R52, already in round 3), late acceptance (R51) and old bachelor acceptance (R50, built for fixed time bounds).
Evidence: R49 (2–24× speedup on standard-cell placement), R48 (parameter-free; one long run > restarts), R50/R51 (general). No floorplanning-specific evidence for R50–R52.
Expected gain: unknown; cheap to test (a few lines each) and the schedule is the one piece of the late tail not yet varied.
Cost: trivial.

**X6. Cluster-then-decluster seeding for region-scale moves.** Cluster tightly connected soft blocks by R33's φ (dead space of the pair plus connectivity density), with preplaced blocks as fixed one-member clusters and boundary blocks kept at the outer level; anneal ~50 super-blocks under the relaxed homotopy, then decluster level by level with the current SA plus the λ-exchange of R34.
Evidence: M6 (the gap is region-scale: 55–66 % of the WL excess sits on blocks displaced > 2 block-scales, so single-block moves cannot fix it); R33 (the only time-matched point: at 196 modules flat B*-tree reaches 6.6 % dead space in 4.7 min where MB*-tree reaches 2.7 % in 2.6 min, area objective; the paper's WL comparisons are not time-matched); R34 (declustering refinement by λ-exchange at 70 % of MB*-tree's runtime).
Expected gain: unknown; the literature evidence is about area and scalability, not WL at a fixed budget, so this stands mainly on M6. It is the only literature-backed answer to region-scale misplacement besides the homotopy. Highest implementation cost in the list, so last.
Cost: a clustering pass, super-block shape handling (sum of member shapes under a fixed internal B*-subtree or SP), a declustering schedule.

Not ranked but recommended as a check: X2's line search doubles as a test of whether `solver/gls.h` already implements R9's exact 1-D minimisation and activation rule; if it uses gradient-like steps instead, that alone may explain part of the 2.1× hand-off.

### 5c. Ideas judged not applicable, and why

- **Slicing-tree SA with preplaced reference points (R28, R29), partition-based slicing legalisation (R17, R18, R56), region-allocating trees (R19), zero-dead-space merging (R31):** all assume few obstacles or none, and soft-dominated instances. 24 interior preplaced blocks fragment the free space (M5); step-function shape curves grow with every preplaced module in a subtree.
- **Constraint-graph / LP repair as a legalizer (R7, R8, R20, R44, R45):** failed on our GP centres (M5); R20 itself reports +19.5 % WL after repair.
- **Constructive packers and order heuristics (R24, R41, R31):** MER greedy left 4–6 % of the area unplaced at 6–7 % whitespace (M5); LLM-evolved orders do not change that.
- **"Shrink the container" continuation (R9 FITS, R13):** the tightened quantity is strip length or bin count, not WL. Its useful form for us is the area-scaling homotopy (X2), not a WL loop.
- **Convex / SOCP / GP legalisation (R16, round-1 P1, P10):** needs a conic solver; and R16 legalises without WL, which is the 2.1× failure mode.
- **RL, diffusion and LLM-in-the-loop methods (R26, R40–R43, round-1 P23–P27):** no training pipeline, no GPU, 110 s single thread; the LLM-designed-heuristic workflow is usable by us as authors but is not a solver component.
- **Rectilinear-module methods (R21, R46, round-1 P5, P14):** shapes not allowed in the lab.
- **Parallel tempering / population annealing:** no evidence on floorplanning or packing was found in any source searched; the single-thread interleaved variant is covered by round 3's go-with-the-winners design, which has the stronger evidence (R37).
- **IAR as an LNS (R5):** measured weak (M5); only the in-anneal move form survives, inside X4.
- **Analog symmetry constraints:** MIB is identical-dimension, not mirrored; already handled by the shared shape index.
- **Independent best-of-k multi-start as the variance fix (R36, R48 FAL):** measured: cuts spread, median flat (M3); replaced by X3.
- **Restricted SP sets (R6), QUBO/SP formulations:** theoretical or hardware-specific.

---
