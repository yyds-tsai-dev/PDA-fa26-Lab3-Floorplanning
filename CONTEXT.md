# Lab3 Floorplanning

PDA Fall 2026 Lab3: a fixed-outline floorplanning assignment derived from ICCAD 2026 CAD Contest Problem C (Intel FloorSet). Each test case is scored by a verifier against fixed thresholds; students are never ranked against each other.

## Language

### Geometry

**Outline**:
The fixed W×H rectangle, lower-left at (0,0), that every block must lie inside. Given in the input file.
_Avoid_: canvas, die, chip boundary, bounding box

**Block**:
A rectangle to be placed; every block is exactly one of soft, fixed-shape, or preplaced.
_Avoid_: module, macro, partition

**Soft block**:
A block whose width and height are chosen by the solver, subject to the area tolerance around its area target.
_Avoid_: flexible block, soft module

**Area tolerance**:
The rule that a soft block's realized area w·h stays within 1% of its area target.

**Aspect ratio bound**:
The rule that every soft block satisfies 1/3 ≤ h/w ≤ 3; part of the gate.

**Fixed-shape block**:
A block whose width and height are given and must be used exactly; its position is free.
_Avoid_: hard macro, hard block

**Preplaced block**:
A block whose position and dimensions are both given and must be used exactly.
_Avoid_: fixed block, pre-placed module

**Whitespace ratio**:
The fraction of outline area not covered by blocks: 1 − (total block area)/(W·H).
_Avoid_: dead space, slack

### Placement constraints

**MIB group**:
A set of blocks (multi-instantiation blocks) that must all have identical width and height.
_Avoid_: same-shape group, multi-instance

**Boundary constraint**:
A requirement that a block touch one or more specified outline edges; requiring two adjacent edges makes it a corner constraint.
_Avoid_: edge constraint, bbox constraint

**Grouping**:
ICCAD Problem C's abutment constraint. Excluded from this lab.

### Grading

**Gate**:
The set of checks every output must pass before anything else is scored: well-formed, all blocks present, overlap-free, inside the outline, area tolerance, aspect ratio bound, fixed-shape and preplaced exactness, and the runtime limit. Failing the gate scores zero for that case.

**Feasible output**:
An output that passes the gate.
_Avoid_: legal solution, valid solution

**Scored constraint**:
A constraint checked only on a feasible output and rewarded with partial points: MIB groups, boundary constraints, and the wirelength threshold.
_Avoid_: soft constraint (that is ICCAD's term, with a different penalty model)

**Wirelength (WL)**:
The weighted sum of centroid-to-centroid Manhattan distances over all block–block and block–terminal connections, as in ICCAD Problem C.
_Avoid_: HPWL (in this lab there are no multi-pin nets)

**WL threshold**:
The per-case wirelength limit given in the input, calibrated from the TA reference solver.

**Terminal**:
A fixed point on the outline that blocks connect to.
_Avoid_: pin, pad, I/O

**Public case** / **Hidden case**:
The one test case released to students, and the four withheld for grading. All five contain every constraint type; difficulty rises through block count, whitespace ratio, and constraint density.
