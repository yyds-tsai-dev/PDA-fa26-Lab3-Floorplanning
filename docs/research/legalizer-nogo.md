# Topology-preserving legalizer: no-go (2026-10-29)

Goal: turn the strong mode's analytical placement (GP) into a legal placement directly, then extract a sequence pair, so SA starts near the GP's wirelength instead of 1.6–1.9× construction WL.

Code (experimental, not used by `solver/fp`): `solver/gp.h`, `solver/legalize.h`, `solver/test_legalize.cpp`, `solver/Makefile.legalize`.

## What works
The legalizer (constraint graphs from centres → longest-path check → repair flips → exact WL LP via `solver/mcf`) legalizes the construction's own centres at 0.95–0.98× construction WL in about 3 ms, and centres jittered by ±20% of block size at about 0.98×.

## Why no-go
- Original GP centres: legal on 2/256 seed runs (public only); 0 on hidden1–4. Initial excess 25k–73k units; repair removes 85–93% then stalls; about 12 s per hidden4 seed.
- Obstacle-aware GP (`GlobalPlacer::obstacle`, `--gp2`: guillotine cuts that avoid preplaced blocks, Stockmeyer shape-curve sizing, boundary blocks forced to their side): halves the initial excess but is still legal on 0/56 hidden4 acceptance seeds. Cuts around 24 interior preplaced blocks leave about 90 small regions, about 60 of which cannot fit their blocks; about 6% of free area stays empty.
- One repair iteration costs 250–340 ms at n = 200, so the approach only pays if the GP output is already nearly legal.

## Not tried
A single sizing tree with preplaced blocks as fixed leaves realized by an LP; a non-guillotine final fit (maximal free rectangles).
