# TA reference solver (Phase F; --mode weak = frozen Phase A). Strong mode
# includes the header-only LP polish from mcf/.
#   make          -> ./fp
#   make test     -> builds and runs the shape unit test
CXX      ?= g++
OPT      ?= -O2
CXXFLAGS ?= $(OPT) -Wall -Wextra
CXXFLAGS += -std=c++17 -I../verifier
# Phase H: optional host-tuned reference build (see docs/research/adhd-solver-round3.md).
# `make fp-ref` builds ./fp_ref with REF_OPT on the build host's own CPU (-march=native:
# do not copy the binary to another machine); -ffp-contract=off keeps the floating-point
# results identical to the default -O2 build (no FMA contraction). The default `make`
# target is unchanged (-O2). `make OPT=... fp` / `make -C mcf test OPT=...` take the same flags.
REF_OPT  ?= -O2 -march=native -ffp-contract=off

all: fp

fp: fp.cpp gls.h shape.h extract.h ../verifier/lab3_io.h mcf/wl_lp.h mcf/network_simplex.h
	$(CXX) $(CXXFLAGS) -o $@ fp.cpp

fp-ref: OPT = $(REF_OPT)
fp-ref: fp.cpp gls.h shape.h extract.h ../verifier/lab3_io.h mcf/wl_lp.h mcf/network_simplex.h
	$(CXX) $(CXXFLAGS) -o fp_ref fp.cpp

test_shape: test_shape.cpp shape.h
	$(CXX) $(CXXFLAGS) -o $@ test_shape.cpp

test: test_shape
	./test_shape

clean:
	rm -f fp fp_ref test_shape

.PHONY: all fp-ref test clean
