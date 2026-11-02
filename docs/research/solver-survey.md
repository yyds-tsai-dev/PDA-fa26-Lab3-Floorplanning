# Solver survey: fixed-outline floorplanning for the Lab3 TA reference solver

Status: research note, 2026-10-29 (Asia/Taipei). Scope: algorithms for the problem in `docs/design.md` §1, §3, §5, §6:

- n ≤ 200 blocks: ~60% soft (integer w,h, area within ±1%, `R·h ≥ w` and `R·w ≥ h` with R = 3), ~28% fixed-shape, ~12% preplaced.
- Fixed outline, possibly ≤ 8% whitespace. No rotation.
- Objective: weighted 2-pin centroid Manhattan WL. MIB groups and BOUNDARY (edge/corner) constraints earn extra points.
- Budget: 24 s, single thread, C++.

Conventions used below:

- `[Pn]` cites a paper in §1.
- **[full]** means I read the full text (or the relevant sections). **[abs]** means abstract or bibliographic record only. **[sec]** means the claim comes from another paper's description of the work.
- **[spec]** means the claim is derived from `docs/design.md`, not from a paper.
- **[proposal]** means my engineering suggestion, not a published result.
- **(uncertain)** marks anything not verified.

---

## 1. Most relevant papers

Compute budgets are listed on purpose. Most published quality numbers come from far more compute than 24 s on one thread, or from commercial solvers (MOSEK) or a GPU. Read the quality columns in that light.

### 1a. Constrained / modern fixed-outline floorplanning (closest to our spec)

| # | Paper | Venue / year | Link | Method (1–2 lines) | Relevance to our constraints | Reported results / runtime scale |
|---|---|---|---|---|---|---|
| P1 | Huang, Ai, Chen, Wu, Zhu, Liu, Zhu, "Toward Near-Optimal Floorplanning with Complex Design Constraints for Modern SoCs" **[full]** | ICCAD 2026 | arXiv [2610.34046](https://arxiv.org/abs/2610.34046), DOI 10.1145/3831252.3834205 | Three stages. (1) Two-step SDP init on circle abstraction: anchors = preplaced + boundary blocks, then pruned non-overlap pairs + MST grouping. (2) Incremental augmented-Lagrangian global FP: smooth HPWL, LSE area, bell-shaped pairwise overlap, MIB aspect-variance penalty, boundary slack; projection for outline/area/AR. (3) SOCP legalization on HCG/VCG. | Same constraint family as ours: soft/hard/preplaced/MIB/boundary/fixed outline, soft AR in [1/3, 3]. Soft area as SOC `‖(h−w)/2, √A‖ ≤ (h+w)/2` plus linear AR bounds. Boundary = equality to HCG/VCG source/sink. MIB = shared (w,h) variable in the SOCP. | FloorSet-Lite (21–120 blocks): zero violations on all 11 cases. HPWL and area are 2% and 1% above ground truth; vs PARSAC −90% HPWL, −16% area. **Compute:** C++ + PyTorch + MOSEK (SDP and SOCP), best of 100 parallel trials, RT 11.7–67.6 s per case on an EPYC 7763. Whether RT is wall time for all 100 trials is unclear (uncertain). |
| P2 | Mostafa, Mallappa, Galkin, Phielipp, Majumdar, "PARSAC: Fast, Human-quality Floorplanning for Modern SoCs with Complex Design Constraints" **[full]** + code read | arXiv 2024 | [2405.05495](https://arxiv.org/abs/2405.05495), code [IntelLabs/parsac](https://github.com/IntelLabs/parsac) | B*-tree SA ("CA-SA") with two additions. (a) Constraint-fixing moves, probability ≈ 0.0005 and always accepted: swap a boundary-violating block with an unconstrained block on that edge, or attach it as a child of a constrained edge block. (b) "Anchored" preplaced blocks. MIB handled by perturbing the aspect ratios of all members identically; AR clamped. Outline and grouping are penalty terms. | Directly answers "preplaced in B*-tree" and "boundary in SA". Code check (`src/c_src/ca_sa.cpp`, `BtreeToFloorplan`, `BtreePreorderTraverseU`): preplaced blocks go into `placed_list` before the DFS. Each movable block is then placed at the lowest y-gap that clears every already-placed block overlapping its x-span. This is O(n) per block, not a contour, and gives no overlap with preplaced blocks by construction. Fixed-outline overflow is not prevented (penalty only). | GSRC n100/n200/n300, 10% whitespace, terminals unscaled. HPWL best–worst over 112 workers: PARSAC(soft) 288–301 / 526–553 / 627–663 vs Parquet 344/654/799, B*-FOPP 323/599/757, GoodFloorplan 309/558/691 (paper's "mm" units). **Compute:** 112 cores × 10 M SA steps; 9.3 / 31.2 / 67.8 min per run. P1 reports PARSAC violations on FloorSet cases larger than lite60. |
| P3 | Mallappa et al., "FloorSet — a VLSI Floorplanning Dataset with Design Constraints of Real-World SoCs" **[full, skimmed]** | ICCAD 2024 | arXiv [2405.05480](https://arxiv.org/abs/2405.05480), DOI 10.1145/3676536.3676814; [IntelLabs/FloorSet](https://github.com/IntelLabs/FloorSet) | Synthetic SoC floorplans with near-optimal-by-construction WL. Two sets: Prime (rectilinear, fully abutted) and Lite (rectangular, < 5% whitespace). Constraints: shape, preplaced, MIB, cluster, boundary. Terminals on the outline. | Source of ICCAD 2026 Problem C, which our lab derives from (`design.md` §1). The repo says the contest uses 21–120 blocks and ships a B*-tree SA baseline in `optimizer_template.py`. The repo is archived. | 1 M training samples, 100 validation, 100 hidden test. No leaderboard or winner info on the repo page as of 2026-10-29. |
| P4 | Yeh et al., "Trinity: One Differentiable Physics for Training, Refining and Scoring Generative Floorplanners" **[full, skimmed]** | arXiv 2026 | [2610.04957](https://arxiv.org/abs/2610.04957) | Flow-matching generator trained with six differentiable constraint/objective terms. The same energy drives a post-sampling refiner. Legalizer: extract a **sequence-pair relation from the refined layout, then one LP** positions the blocks, plus a contact snap. | The legalizer (SP extracted from a continuous layout, then LP) is the same handoff we need in candidate B/C. Also states what ML can reach on FloorSet. | FloorSet val: mean ICCAD-2026 hard cost 1.014 (1.0 ≈ reference) at 1.63 s/case on one RTX 5090 GPU with parallel legalization. Finding: "the shared refiner decides more of the final cost than the generator". |
| P5 | Chen, Lin, Hsu, Jiang, Chen, Chen, Chang, "Modern Fixed-Outline Floorplanning with Rectilinear Soft Modules" **[abs]** | ICCAD 2024 | DOI [10.1145/3676536.3676818](https://doi.org/10.1145/3676536.3676818) | Analytical. Modules are shape-adjustable rectangles with a differentiable shape mechanism. Legalization: overlap graph, then shortest-path migration of area from overlaps through chains of modules to whitespace; then expand/shrink refinement. | Handles a fixed outline with **preplaced** modules. The area-migration legalizer is an alternative to constraint-graph LP. Rectilinear shapes are not allowed in our lab. | Abstract claims better WL and runtime than prior work and the 2023 contest winner. RPolyPlan (P14) lists its GSRC 10%-whitespace runtimes: n100 6.83 s, n200 42.92 s, n300 62.34 s. |
| P6 | Young, Chu, Ho, "Placement Constraints in Floorplan Design" **[full, relevant sections]** | IEEE TVLSI 12(7):735–745, 2004 | DOI [10.1109/TVLSI.2004.830915](https://doi.org/10.1109/TVLSI.2004.830915) | Sequence pair + constraint graphs. Preplace, range, boundary, alignment, abutment and clustering become extra "constraining edges" (difference constraints) to the source/sink boundary nodes. A positive cycle means infeasible, which is penalized in SA. Uses incremental updates and a reduced graph for speed. | The cleanest way to put preplaced + boundary (edge and corner) into **SP + longest-path / LP**. It maps 1:1 onto an LP/min-cost-flow WL polish. PARSAC [P2] calls it (its ref. [32], confirmed in P2's reference list) the one prior work that covers all placement constraints, but with soft penalties and sub-optimal whitespace. | Tested with about 1/8 of modules constrained; the authors call results "very promising". Small classic benchmarks only. |
| P7 | Murata, Fujiyoshi, Kaneko, "VLSI/PCB Placement with Obstacles Based on Sequence Pair" **[abs]** | IEEE TCAD 17(1):60–68, 1998 | DOI 10.1109/43.673633 | Preplaced modules make some SPs inconsistent. An "adaptation" procedure minimally modifies the SP to restore consistency, inside SA. | The classic way to handle preplaced blocks in SP. | ami49 with 20% of modules preplaced; one PCB example. |
| P8 | Lin, Yi, Chang, "Module Placement with Boundary Constraints Using B*-trees" **[abs]** | IEE Proc. CDS 149(4):251–256, 2002 | DOI 10.1049/ip-cds:20020433 | Derives feasibility conditions for boundary constraints in a B*-tree and keeps the tree feasible after every SA perturbation. | Boundary handling native to B*-tree. Left edge = leftmost branch, bottom = root's left chain (see the B*-tree definition). | Abstract: smaller area than SP-based boundary work. Area-driven, not fixed-outline WL. |

### 1b. Analytical global floorplanning + legalization

| # | Paper | Venue / year | Link | Method | Relevance | Reported results / runtime |
|---|---|---|---|---|---|---|
| P9 | Li, Peng, Huang, Zhu, "PeF: Poisson's Equation-Based Large-Scale Fixed-Outline Floorplanning" **[full]** | IEEE TCAD 2022 | arXiv [2210.03293](https://arxiv.org/abs/2210.03293), DOI 10.1109/TCAD.2022.3213609 | Flat electrostatic (ePlace-style) global FP. Soft-module **widths are variables** in the potential energy; Nesterov for positions alternates with projected gradient for widths. Legalization treats soft modules as hard. Small cases: Parquet's `pl2sp()` (convert placement to SP, then pack). Large cases: Floorist-style HCG/VCG edge flipping and compression until the bounding box fits. | Shows a published **analytical-to-SP handoff** for small cases. Shows that soft shaping during global FP pays off. | GSRC, 10% whitespace, pads on the boundary (like our terminals), all soft. HPWL n100/n200/n300 = 194 273 / 349 301 / 468 235. That is 6% below DeFer and 29% below Parquet-4 on average. HB+ (~900–1500 modules, 20–30% whitespace): 4.3–8.1 s on the first four cases (i3-9100). About 80% of runtime is global FP. |
| P10 | Yimamu, Pastor-Rué, Cortadella, "Chip Floorplanning Combining Convex and Non-convex Optimization" **[full, relevant sections]** | arXiv 2026 | [2607.21408](https://arxiv.org/abs/2607.21408) | Quadratic placement init, then Adam projected gradient on WL + normalized smooth overlap, then legalization on HCG/VCG. Legalization uses linear-domain coordinates and **log-domain sizes**, which makes soft-area constraints convex (a GP). Solved with IPOPT. | Log-domain sizing is a convex, solver-friendly way to shape soft blocks inside a fixed constraint graph. | GSRC, 10% whitespace, all soft: n100/n200/n300 = 189 635 / 345 246 / 466 365, best published rectangular numbers I found (1–2% below PeF). Python/IPOPT; about 5× PeF runtime (26–107 s on HB+). |
| P11 | Huang, Liu, Li, Yu, Zhu, "Handling Orientation and Aspect Ratio of Modules in Electrostatics-Based Large Scale Fixed-Outline Floorplanning" **[abs]** | ICCAD 2023 | DOI 10.1109/ICCAD57390.2023.10323841 | Electrostatic global FP with sizing. Legalization = constraint-graph adjustment, then module sizing, then an **LP that minimizes wirelength** at the end. | Published use of an LP WL-refinement step after constraint-graph legalization. P1 builds its HCG/VCG with this paper's method. | Abstract: ≥ 5.9% (HB+) and 11% (ami49_x) HPWL reduction vs prior SOTA. |
| P12 | Li, Wang, Moura, Blanton, "Global Floorplanning via Semidefinite Programming" **[abs/sec]** | DAC 2023 | DOI 10.1109/DAC56929.2023.10247967 | Circle abstraction. SDP with a rank-2 constraint replaced by a direction matrix and iterative rank control. | Basis of P1's initialization. A robust topology seed, but it needs an SDP solver. | Appendix (via web summary): −19.1% WL vs AR, −45.4% vs UFO. P1's ablation: SDP init + incremental ALM gives about 4% HPWL over center init. |
| P13 | Lin, Hung, "UFO: Unified Convex Optimization Algorithms for Fixed-Outline Floorplanning Considering Pre-Placed Modules" **[abs]** | IEEE TCAD 30(7):1034–1044, 2011 (ASP-DAC 2010 version) | DOI 10.1109/TCAD.2011.2114531 | Stage 1: circles with a convex push-pull spreading model. Stage 2: SOCP for location + shape with overlap penalty and boundary constraints. Extended to **pre-placed modules**. | Early convex flow with preplaced support. | Abstract: outperforms prior published GSRC results. ASP-DAC version ran in MATLAB. |
| P14 | Huang, Zhu, Liu, Zhu, "RPolyPlan: A Unified Global Floorplanning and Legalization Framework for Rectilinear Soft Modules" **[abs + tables]** | arXiv 2026 | [2610.28276](https://arxiv.org/abs/2610.28276) | Augmented-Lagrangian global FP with differentiable polygon shaping. Legalization: **min-cost flow** decides target areas, HCG/VCG give relations, then a convex program (MOSEK; Lemon for MCF). | Shows MCF + constraint graph + convex shaping as a legalization pattern. Rectilinear shapes do not apply here. | GSRC 10% whitespace (rectilinear, so not comparable with P9/P10): n200 10.4 s, n300 13.7 s. |
| P15 | Meng et al., "CSF: Fixed-outline Floorplanning Based on the Conjugate Subgradient Algorithm Assisted by Q-Learning" **[abs]** | Electronics 2025 / arXiv | [2504.03796](https://arxiv.org/abs/2504.03796), DOI 10.3390/electronics14244893 | Non-smooth analytical model solved by a population-based conjugate subgradient method; Q-learning adapts the step size. Also legalizes via CSA. | Population (multi-start) analytical is a published idea. | Abstract: competitive with SOTA on hard-module-only cases. Claims more efficient legalization than constraint-graph legalization. |
| P16 | Moffitt, Ng, Markov, Pollack, "Constraint-Driven Floorplan Repair" (Floorist) **[sec]** | DAC 2006 (TODAES 2008) | DOI 10.1109/DAC.2006.229405 | Repairs an overlapping floorplan by choosing HCG/VCG relations and compressing (constraint graph + LP). | P9's legalizer is "based on Floorist". | Not read. |

### 1c. Packing-based SA and soft-block shaping

| # | Paper | Venue / year | Link | Method | Relevance | Reported results / runtime |
|---|---|---|---|---|---|---|
| P17 | Chen, Chang, "Modern Floorplanning Based on B*-Tree and Fast Simulated Annealing" **[abs]** | IEEE TCAD 25(4):637–650, 2006 | DOI 10.1109/TCAD.2006.870076 | B*-tree + three-stage Fast-SA. For fixed outline, an **adaptive** Fast-SA changes the cost weights during annealing. | Strong SA baseline; most students' starting point. | Abstract: 100% fixed-outline success at various aspect ratios. |
| P18 | Adya, Markov, "Fixed-Outline Floorplanning: Enabling Hierarchical Design" (Parquet) **[sec]** | IEEE TVLSI 11(6):1120–1135, 2003 | DOI 10.1109/TVLSI.2003.817546 | SP-based SA for fixed outline. Its `pl2sp` conversion is reused by PeF [P9]. | Baseline. The open-source Parquet code is the classic reference. | P9: Parquet-4 is 29% worse in HPWL than PeF on GSRC 10% whitespace. |
| P19 | Yan, Chu, "DeFer: Deferred Decision Making Enabled Fixed-Outline Floorplanning Algorithm" **[full, local `docs/defer-tcad10.pdf`]** | IEEE TCAD 2010 | DOI 10.1109/TCAD.2010.2041850 | Generalized slicing tree + shape-curve deferral; one shape curve stores many slicing solutions; then compaction. | Fast, deterministic soft-block packer; a seed or fallback. | Abstract: best success rate and WL among compared tools then. P9: DeFer runtime is 0.89× PeF on HB+; HPWL 6% worse than PeF on GSRC. |
| P20 | Yan, Chu, "Optimal Slack-Driven Block Shaping Algorithm in Fixed-Outline Floorplanning" (SDS) **[full, local `docs/slackFP-ispd12.pdf`]** | ISPD 2012 | DOI 10.1145/2160916.2160956 | For a fixed non-slicing topology (constraint graphs) and a width bound, iteratively shapes soft blocks by distributing global slack. Layout height decreases monotonically. A GP step guarantees optimality. Each iteration is linear time. | Exactly the "given SP/constraint graph, fit into the outline by reshaping soft blocks" subproblem. | 253× and 33× faster than two prior optimal shapers at similar quality (MCNC/HB). |
| P21 | Young, Chu, Luk, Wong, "Handling Soft Modules in General Nonslicing Floorplan Using Lagrangian Relaxation" **[abs-level title only]** + published comment DOI 10.1109/TCAD.2003.819428 | IEEE TCAD 2001 | DOI 10.1109/43.920707 | Lagrangian relaxation for optimal soft-module shaping in a fixed non-slicing (constraint graph) topology; area objective. | Closest match to the requested "Young–Chu–Shen, *Optimal shapes for soft modules*". I could not find a paper with that exact title or author set. A 2003 comment on this paper exists, so read both before relying on the derivation. | Not read. |
| P22 | Chang, Chang, Wu, Wu, "B*-Trees: A New Representation for Non-Slicing Floorplans" **[full, local `docs/Bstar-chang-dac00.pdf`]** | DAC 2000 | local PDF | B*-tree. Preplaced modules: if a preplaced module cannot reach its fixed position during packing, **exchange** it with another module whose DFS position allows the fixed position (choose the one closest to it). Modules earlier in DFS order stay unchanged. | Primary source for preplaced in B*-tree. PARSAC [P2] says this approach "leads to large white-space". | Area-driven MCNC results only. |

### 1d. ML / learning-based (what is achievable, not what to build)

| # | Paper | Venue / year | Link | Method | Relevance | Results |
|---|---|---|---|---|---|---|
| P23 | Xu, Geng, Chen, Yuan, Zhuo, Kang, Wen, "GoodFloorplan: Graph Convolutional Network and RL-Based Floorplanning" **[sec]** | IEEE TCAD 41(10):3492–3502, 2022 | DOI 10.1109/TCAD.2021.3131550 | GCN + actor-critic searching neighbors in a sequence-pair space [P2's description]. | Shows RL over SP is about as good as tuned SA, not better. | Per P2 Table 2 (GSRC fixed outline): 309/558/691 vs PARSAC(hard) best 291/530/628. |
| P24 | Liu, Ju, et al., "GraphPlanner: Floorplanning with Graph Neural Network" **[abs]** | ACM TODAES 28, 2022 | DOI 10.1145/3555804 | Variational GCN for macro floorplanning with density limits. | Low. Targets the mixed-size placement setting. | Abstract: 25% runtime and 4% WL improvement vs mixed-size placers alone. |
| P25 | Lee, Nguyen, Elzeiny, Deng, Abbeel, Wawrzynek, "Chip Placement with Diffusion Models" **[abs]** | ICML 2025 | arXiv [2407.12282](https://arxiv.org/abs/2407.12282) | Diffusion model with guided sampling; synthetic pre-training data. | One of the four diffusion placers that Trinity [P4] re-implemented on FloorSet. | Zero-shot placement competitive with SOTA (abstract). |
| P26 | Shi, Xue, Song, Qian, "Macro Placement by Wire-Mask-Guided Black-Box Optimization" **[abs]** | NeurIPS 2023 | arXiv [2306.16844](https://arxiv.org/abs/2306.16844), code [lamda-bbo/WireMask-BBO](https://github.com/lamda-bbo/WireMask-BBO) | Greedy wire-mask placement used as the decoder inside black-box optimization (EA/BO). | A sequence-to-placement greedy decoder + EA is a cheap alternative to SA. It also fine-tunes existing placements. | Abstract: up to 50% HPWL improvement when fine-tuning existing placements (macro placement, not fixed-outline FP). |
| P27 | Cheng, Kahng, Kundu, Wang, Wang, "Assessment of Reinforcement Learning for Macro Placement" **[sec]** | ISPD 2023 | DOI 10.1145/3569052.3578926 | Independent re-evaluation of Google's RL macro placer. | P1 cites this work (and its TCAD 2025 update) for the limited generalization of RL. | Not read. |

Also seen but not tabulated: Luo/Anjos/Vannelli ASP-DAC 2008 convex FP (DOI 10.1109/ASPDAC.2008.4483939); Zhan/Feng/Sapatnekar ASP-DAC 2006 analytical fixed-die FP; Lee et al. DAC 2003 MB*-tree (local `docs/mlbstar-dac2003.pdf`); He et al. 2014 iterative merging for soft rectangles (arXiv 1401.3172); JigsawPlanner (ICCAD 2024) and FlexPlanner (NeurIPS 2024), both cited in P1.

### 1e. Open-source code worth reading

| Repo | What it gives us | Notes |
|---|---|---|
| [IntelLabs/parsac](https://github.com/IntelLabs/parsac) (C++ core `src/c_src/ca_sa.cpp`, Python wrapper) | B*-tree SA with boundary-fixing moves, anchored preplaced blocks with gap-aware packing, MIB identical-AR perturbation. | Closest public code to our constraint set. Packing is O(n²) per evaluation (loop over `placed_list`). Uses an HPWL + area cost, not our objective. |
| [IntelLabs/FloorSet](https://github.com/IntelLabs/FloorSet) | Contest data loaders, `cost.py`/`validate.py`, B*-tree SA baseline (`optimizer_template.py`). | Archived. Useful to see which baseline the ICCAD 2026 entrants started from. |
| [Aurora-73/PeF-Poissons_Equation-Based_Floorplanning](https://github.com/Aurora-73/PeF-Poissons_Equation-Based_Floorplanning) | Independent C++ PeF re-implementation (2 stars). | Not the authors' code. Quality unverified (uncertain). |
| [lamda-bbo/WireMask-BBO](https://github.com/lamda-bbo/WireMask-BBO) | Wire-mask greedy decoder + BBO. | Macro-placement setting. |
| [limbo018/DREAMPlace](https://github.com/limbo018/DREAMPlace) | ePlace-style electrostatic density + Nesterov reference implementation (C++/CUDA/PyTorch). | Reference for the density FFT/DCT and Nesterov step size, if candidate B is chosen. |
| OpenROAD `src/mpl` (Hier-RTLMP, [arXiv 2304.11761](https://arxiv.org/abs/2304.11761)) | Production SA on sequence pair for macro clusters; move probabilities 0.3/0.3/0.3/0.1. | C++. Other clusters are treated as fixed terminals during optimization. |
| [hankshyu/Fixed-outline_Floorplanner](https://github.com/hankshyu/Fixed-outline_Floorplanner), [EnHsiangChu/2023_ICCADContest_ProblemD](https://github.com/EnHsiangChu/2023_ICCADContest_ProblemD) | Student contest solvers: fixed outline, rectilinear soft blocks, preplaced. | Example of what strong students build. Not audited. |

---

## 2. Spec-derived facts that drive the ranking

No paper matches our exact spec. These properties decide which ideas transfer.

- **D1 [spec] The objective is exact L1, separable in x and y.** WL is a sum of weighted 2-pin centroid |Δx| + |Δy| (`design.md` §5). There is no multi-pin HPWL, so no max/min auxiliary variables per net are needed.
  - Fix the relative order (an SP, or HCG/VCG) and every block shape. Then minimizing WL subject to non-overlap edges `x_j ≥ x_i + w_i`, outline bounds, preplaced equalities and boundary equalities is **one LP per axis**.
  - Every constraint is a difference constraint, and the objective is a weighted sum of |x_a − x_b| and |x_a − t|. So the LP dual is a **min-cost flow**, solvable exactly by a hand-rolled network simplex.
  - P11 (LP WL-minimization after legalization) and P4 (SP extraction + one LP) are published uses of this step. The MCF-dual observation is standard LP duality, not taken from a specific paper here.
- **D1b [spec + papers] The final stage must optimize WL, not area.**
  - P1's SOCP legalization minimizes W + H (an area proxy) and leaves WL to the global stage [P1 eq. 9]. FloorSet's cost includes area [P3; P1 eq. 1].
  - Our objective has no area term. The outline is given, and WL is the only scored continuous objective (`design.md` §6).
  - So the legalization/polish objective should be WL, as in P11's final LP. That is why the D1 polish is the right last stage.
- **D2 [spec] Soft shapes are a small discrete set.** w and h are integers, with `100·|w·h − A| ≤ A` and `3h ≥ w`, `3w ≥ h`. So each soft block has a precomputable finite list of legal (w, h), and an MIB group shares one index into it.
  - The tolerance permits `w·h` down to 0.99·A. With 60% soft area, shrinking every soft block by up to 1% frees about 0.6% of the outline as extra whitespace. That is a meaningful share when total whitespace is 8%.
  - Consequence: continuous shaping methods (SOCP [P1], GP [P10], SDS [P20]) need a final integer rounding step. A discrete shape index fits SA directly.
- **D3 [spec + papers] Low whitespace with a large non-soft fraction.** Published analytical results on GSRC/MCNC use 10–15% whitespace with all-soft modules [P9, P10]. P1 succeeds on FloorSet-Lite (< 5% whitespace [P3]), but FloorSet-Lite is about 85–90% soft by block count (P1 Table 1), and P1 uses MOSEK SOCP with shape freedom in legalization.
  - Our cases have ~40% fixed-shape or preplaced blocks. That is less shaping freedom, and interior preplaced blocks fragment the free space.
  - A constraint-graph legalizer has no feasibility guarantee in that regime (uncertain how often it fails). A packing representation with an outline penalty plus fit repair is the safer floor.
- **D4 [spec] 24 s, single thread, plain C++.** Several top results depend on resources we do not have:
  - P1 uses MOSEK + PyTorch and 100 parallel trials.
  - P2 uses 112 cores × 9–68 min.
  - P4 uses a GPU and a trained model.
  - P10 uses IPOPT.
  - The single-thread-capable published flows are PeF [P9] (seconds for about 1000 modules), DeFer [P19] and SDS [P20], plus classic SA.
- **D5 [papers] Initialization and refinement both matter.**
  - P1's ablation: the full flow is about 4.0% better in HPWL than "w/o Init & Incre". About 2.0% of that comes from the incremental ALM ("full vs w/o Incre"), so SDP init accounts for roughly the other 2%.
  - P4: the refiner decides more of the final cost than the generator.
  - P26: fine-tuning an existing placement improved HPWL by up to 50% (macro setting).
  - Taken together, a good topology seed plus a strong local polish beats either one alone (moderate confidence; the settings differ).

**Open question for the caller:** may the solver link a third-party LP/MCF library (for example LEMON or OR-Tools), or must it be dependency-free? If it must be dependency-free, the D1 polish needs a hand-rolled network simplex or a dual longest-path method. That is about 300–600 lines (my estimate) and changes the effort ranking below.

---

## 3. Candidate architectures, ranked

Ranking criterion: expected lowest WL that is feasible every time, within 24 s single-thread C++, for n ≤ 200, with soft/fixed/preplaced blocks, MIB and boundary constraints.

### Rank 1 — C: Hybrid (analytical seeds, packing-based SA, exact LP/MCF polish) **[proposal]**

Pipeline:

1. **Analytical seed (≈ 1–3 s total, several seeds).**
   - Start with quadratic placement on the 2-pin net graph, with terminals and preplaced centers as anchors (P10 stage 1). Boundary blocks are anchored to their edge, as in P1's anchor set.
   - Then run a short Nesterov/Adam spread. Use either an electrostatic density (P9) or a pairwise bell-shaped overlap (P1, from Kahng–Wang) with soft aspect ratio as a variable.
   - n ≤ 200, so a pairwise O(n²) overlap term is cheap. No FFT is needed (an estimate: 200² = 40k pairs per gradient).
2. **Convert to a representation.** Extract a sequence pair from the continuous layout: for each pair, the dominant separation direction (P4 does exactly this; P9 uses Parquet `pl2sp`).
3. **SA in SP space (bulk of the budget).** The state is the SP plus a discrete shape index per soft block (D2), with one index per MIB group.
   - Evaluate with longest path on HCG/VCG. Preplaced and boundary constraints become constraining edges to the source/sink (P6). Infeasibility (positive cycle, or outline overflow) is a graded penalty.
   - Add PARSAC-style constraint-fixing moves (P2), adapted to SP [proposal, derived from the SP rules in P6 §III-A]. The rules: b is left of a ⟺ b is before a in both Γ+ and Γ−; b is below a ⟺ b is after a in Γ+ and before a in Γ−.
     - Exact condition for LEFT: no block is before a in both sequences. That leaves a with no HCG predecessor, so `x = 0` is reachable. The other edges are analogous.
     - Sufficient single-position moves:
       - LEFT: first in Γ+ or first in Γ−.
       - RIGHT: last in Γ+ or last in Γ−.
       - BOTTOM: last in Γ+ or first in Γ−.
       - TOP: first in Γ+ or last in Γ−.
     - Corners need one position in one sequence:
       - bottom-left: first in Γ−.
       - top-left: first in Γ+.
       - bottom-right: last in Γ+.
       - top-right: last in Γ−.
     - Each corner position holds one block. The generator guarantees that a fully satisfying solution exists (`design.md` §3).
     - Preplaced blocks still constrain the real coordinates. An edge block can be pushed off its edge only by a preplaced block's fixed-position edge, so check the constraining edges (P6).
   - Use adaptive fixed-outline weights (P17). Start at a low temperature, since the seed already has a good topology (P1's locality/proximal idea plays the same role).
   - Fast SP evaluation: FAST-SP LCS, O(n log log n) (local `docs/fastsp-*.pdf`). With constraining edges you need the explicit graph or an incremental longest path (P6 uses incremental updates).
4. **Exact WL polish (D1).** Fix the SP and shapes, then solve the per-axis LP via MCF for WL-optimal coordinates inside the outline, with boundary equalities.
   - Packing-based SA puts every block at its bottom-left position, which biases WL. The LP removes that bias by spreading whitespace optimally.
   - Run the polish on the incumbent periodically, or on the top-k states.
   - Then do a greedy shape-index sweep with the LP re-solved per accepted change. This is cheap at n = 200, but the per-solve time is not measured yet (uncertain).
5. **Fallback.** If no feasible state exists at the deadline, return the best feasible packing found. A DeFer-style [P19] or B*-tree packer can guarantee a fit early in the run.

Pros:

- Feasibility comes from the packing representation (D3).
- Topology quality comes from analytical seeding (D5).
- The final coordinates are WL-optimal for the chosen topology (D1).
- Every component runs single-threaded in plain C++ (D4).
- It extends candidate A, so A can be built first and kept as a fallback.

Cons:

- Most code of the four candidates: global placer, SP extraction, SA, MCF.
- How much the analytical seed helps at our scale is unmeasured. P1's 4% is on FloorSet with SDP.

Techniques to borrow, with sources:

- Preplaced in SP: constraining edges `x_s + X ≤ x_i ≤ X` (P6). Alternative: P7's SP adaptation.
- Preplaced in B*-tree: PARSAC gap-aware packing with preplaced blocks pre-inserted (P2 code), or the DFS-position exchange (P22).
- Boundary in SP: edge-to-source/sink equality constraints (P6). Boundary fixing moves (P2). Corner = two equalities (P6 Fig. 8).
- MIB: shared shape variable (P1 legalization; P2 identical perturbation).
- Soft shaping inside a fixed graph: SDS slack distribution (P20) or the log-domain GP (P10), then integer rounding (D2).

### Rank 2 — A: Packing SA (SP or B*-tree) with discrete soft shapes plus LP/MCF post-polish **[baseline + proposal]**

The same as C without the analytical seed: random or cluster-based initial SP, then a full Fast-SA schedule (P17), then the D1 polish.

- This is what strong students will build, likely as B*-tree Fast-SA from P17/P22 or the FloorSet baseline [P3]. Our edge over them must come from the D1 exact polish, the discrete-shape SA with the −1% area slack (D2), and engineering speed (more evaluations in 24 s).
- Pros: simplest robust design. Proven at 100% fixed-outline success [P17]. Easy to reason about.
- Cons: WL gap vs analytical. On GSRC at 10% whitespace, Parquet-4 is 29% worse than PeF [P9], although Parquet is an old SA. Convergence within 24 s at n = 200 is the main risk (uncertain).
- Calibration from P2 Table 2: PARSAC's 112 workers serve best-of-N.
  - The best–worst spread across workers is 288–301 (n100), 526–553 (n200) and 627–663 (n300). That is a 4.5–5.7% spread, which bounds what multi-start buys for SA at that step count.
  - One worker runs 10 M steps in 9.3 min (n100), about 18 k steps/s in their B*-tree C++ with an O(n²) packer.
  - These are calibration points only, not predictions for our code.
- SP vs B*-tree here:
  - SP + constraint graphs takes preplaced and boundary as edges directly (P6) and feeds the D1 LP with no conversion.
  - B*-tree needs anchoring (P2) or exchange (P22) for preplaced blocks and feasibility rules for boundaries (P8). Its compacted packings are also bottom-left biased, so the D1 polish is still needed.
  - Recommendation [proposal]: SP, because the polish and the constraints share one data structure.

### Rank 3 — B: Pure analytical (global FP, constraint-graph legalization, LP/convex refinement)

The P1, P9, P10, P11 family, without an SA stage.

- Pros: best published WL quality on fixed-outline benchmarks [P9, P10]. Near-optimal on FloorSet with constraints [P1]. Fast for small n [P9]. Deterministic.
- Cons:
  - Legalization feasibility under ≤ 8% whitespace with 40% non-soft blocks and interior preplaced blocks is unproven (D3).
  - P1 relies on MOSEK for the SOCP; a hand-rolled SOCP/GP in 24 s is a serious engineering risk (D4).
  - Constraint-graph edge-flipping legalizers (P9, Floorist P16) are heuristic and can fail to fit.
  - Any failure needs a fallback, and the natural fallback is a packing SA. That turns B into C.
- Worth borrowing anyway: P1's ALM penalty balancing (ρ0 = ‖∇F‖∞ / max(1, ‖C‖∞)), the proximal locality term, and the MIB variance penalty.

### Rank 4 — D: ML-guided (diffusion / RL / GNN)

- Not viable for a single-thread 24 s C++ solver with no training pipeline. Trinity [P4] needs a GPU and a trained model, even though it reaches hard cost 1.014 at 1.63 s/case.
- RL floorplanners match SA rather than beat it on GSRC [P2 Table 2 vs P23]. Independent assessments question RL generalization [P27, as cited by P1].
- Value for us: FloorSet ground truth and Trinity's numbers show how close to optimal good solvers get on FloorSet-like data (within about 1–2% WL [P1, P4]). That supports setting WL thresholds from a strong optimized solver (ADR 0002).

### Suggested 24 s budget for candidate C **[proposal, unmeasured]**

| Phase | Time | Purpose |
|---|---|---|
| Parse + shape tables (D2) | < 0.1 s | Discrete shape lists; MIB shared shapes |
| Analytical seeds, 3–5 runs | 2–4 s | Different penalty schedules / random perturbations |
| SA from the best seeds (low T) + periodic MCF polish of the incumbent | 16–18 s | Main search |
| Final MCF polish + greedy shape sweep + verification | 1–2 s | WL-optimal coordinates; recheck every gate item from `design.md` §6 |

Build order that de-risks the work:

1. A: SP + constraint graph + constraining edges + SA, with a feasibility-first penalty.
2. The D1 MCF polish, verified against a generic LP on small cases.
3. The analytical seed (C).
4. Measure C vs A on the generated cases. Keep whichever wins; C's runtime always includes A as a fallback.

---

## 4. Claim checklist (checkable, cited)

1. P1 reaches zero violations on all 11 FloorSet-Lite cases, with HPWL and area 2% and 1% above ground truth, using MOSEK + PyTorch and 100 parallel trials. [P1 §4, Table 2]
2. P1's SOCP soft-block constraint is `‖((h_i − w_i)/2, √A_i)‖ ≤ (h_i + w_i)/2` plus linear AR bounds. MIB members share one (w, h) variable. Boundary blocks get equality to the HCG/VCG source/sink. [P1 §3.3]
3. PARSAC's constraint-fixing moves fire with probability ≈ 0.0005 and are always accepted. [P2 §3.2.1, Alg. 1]
4. In the PARSAC code, preplaced blocks are inserted into `placed_list` before the B*-tree DFS. Movable blocks are placed at the lowest y-gap clear of every already-placed block in their x-span. [P2 code: `BtreeToFloorplan`, `BtreePreorderTraverseU`]
5. PARSAC runs 112 workers × 10 M steps: 9.3 / 31.2 / 67.8 min for n100/n200/n300. [P2 Table 2]
6. PeF HPWL on GSRC (10% whitespace, pads on the boundary) is 194 273 / 349 301 / 468 235, 6% below DeFer and 29% below Parquet-4. [P9 Table IV]
7. PeF legalizes small floorplans with Parquet's `pl2sp()`. [P9 §IV-B]
8. Yimamu–Pastor-Rué–Cortadella improve on PeF to 189 635 / 345 246 / 466 365 on the same setting, with a log-domain convex legalization solved by IPOPT. [P10 Table II, §III-B]
9. Young–Chu–Ho encode preplace and boundary constraints as extra constraint-graph edges in an SP floorplanner and detect infeasibility as a positive cycle. [P6 §IV]
10. B*-tree (DAC 2000) handles a preplaced module that cannot reach its fixed position by exchanging it with a module whose DFS position allows it. [P22, pre-placed module discussion and Fig. 5]
11. SDS reshapes soft blocks optimally for a fixed topology and width bound, with linear-time iterations, and is 253× / 33× faster than two prior optimal shapers. [P20]
12. Trinity's legalizer extracts a sequence-pair relation from the refined layout and solves one LP. It reaches hard cost 1.014 at 1.63 s/case on a GPU. [P4 §3.3, §4.4]
13. FloorSet-Lite has < 5% whitespace. ICCAD 2026 Problem C uses 21–120 blocks. [P3 abstract; FloorSet repo README]

Not verified and flagged:

- Whether P1's RT column is wall time for all 100 trials.
- Whether analytical legalization stays feasible at ≤ 8% whitespace with our soft/hard/preplaced mix.
- The exact content of P5, P11, P13, P16, P21 beyond their abstracts.
- The existence of a paper titled "Optimal shapes for soft modules" by Young–Chu–Shen.
- All runtime estimates in §3.
