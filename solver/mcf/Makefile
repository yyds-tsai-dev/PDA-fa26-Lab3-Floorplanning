# solver/mcf -- exact fixed-topology WL LP (network simplex), header-only.
#   make            -> ./polish (CLI)
#   make test       -> random MCF vs SSP oracle, brute force, dev cases vs the
#                      scipy optimum of the E3 run, verifier checks, timing
#   make E3_DIR=-   -> skip the E3 comparisons (only needs cases/dev)
CXX      ?= g++
OPT      ?= -O2
CXXFLAGS ?= $(OPT) -Wall -Wextra -Werror
CXXFLAGS += -std=c++17 -I../../verifier
ROOT     := ../..
E3_DIR   ?= testdata/e3
HDRS     := network_simplex.h wl_lp.h ../../verifier/lab3_io.h

all: polish

polish: polish.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -o $@ polish.cpp

test_ns: test_ns.cpp network_simplex.h
	$(CXX) $(CXXFLAGS) -o $@ test_ns.cpp

test_wl: test_wl.cpp test_util.h $(HDRS)
	$(CXX) $(CXXFLAGS) -o $@ test_wl.cpp

test_timing: test_timing.cpp test_util.h $(HDRS)
	$(CXX) $(CXXFLAGS) -o $@ test_timing.cpp

test: test_ns test_wl test_timing polish
	@test -x $(ROOT)/verifier/verifier || { echo "build ../../verifier/verifier first (make -C verifier)"; exit 1; }
	./test_ns
	./test_wl $(ROOT) $(E3_DIR) test_out
	./test_timing $(ROOT)
	./polish $(ROOT)/cases/dev/public_s1001.in $(ROOT)/cases/dev/public_s1001.construction.out test_out/cli.out
	$(ROOT)/verifier/verifier $(ROOT)/cases/dev/public_s1001.in test_out/cli.out > /dev/null
	@echo "ALL MCF TESTS PASSED"

clean:
	rm -rf polish test_ns test_wl test_timing test_out

.PHONY: all test clean
