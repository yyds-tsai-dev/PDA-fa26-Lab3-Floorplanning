# Lab3 verifier (design.md §6).
#   make                  -> ./verifier, libstdc++/libgcc linked statically, glibc
#                            dynamic (needs glibc >= 2.25 on the target host)
#   make STATIC=-static   -> fully static build; needs glibc-static (libc.a),
#                            which this TA host does not have
#   make STATIC=          -> plain dynamic build (dev only)
CXX      ?= g++
CXXFLAGS ?= -O2 -Wall -Wextra
CXXFLAGS += -std=c++17
STATIC   ?= -static-libstdc++ -static-libgcc
LDFLAGS  += $(STATIC) -s

verifier: verifier.cpp lab3_io.h
	$(CXX) $(CXXFLAGS) -o $@ verifier.cpp $(LDFLAGS)

clean:
	rm -f verifier

.PHONY: clean
