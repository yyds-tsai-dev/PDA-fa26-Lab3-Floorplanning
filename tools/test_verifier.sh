#!/usr/bin/env bash
# Regression test for verifier/verifier on the toy case.
#   bash tools/test_verifier.sh [path/to/verifier]
# Exits non-zero on any mismatch.
set -u

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VER="${1:-$ROOT/verifier/verifier}"
TOY="$ROOT/cases/toy"
IN="$TOY/toy.in"

if [[ ! -x "$VER" ]]; then
    echo "verifier not found or not executable: $VER (run 'make -C verifier')" >&2
    exit 2
fi

TMP="$(mktemp -d "${TMPDIR:-/tmp}/lab3_test_verifier.XXXXXX")"
trap 'rm -rf "$TMP"' EXIT

pass=0
fail=0

# check <label> <input> <output> <expected exit> <expected substring>
check() {
    local label="$1" in="$2" out="$3" want_rc="$4" want_txt="$5"
    local got rc
    got="$("$VER" "$in" "$out" 2>&1)"
    rc=$?
    if [[ "$rc" == "$want_rc" ]] && grep -qF -- "$want_txt" <<<"$got"; then
        printf 'ok    %-22s exit=%s  %s\n' "$label" "$rc" "$want_txt"
        pass=$((pass + 1))
    else
        printf 'FAIL  %-22s want exit=%s and "%s"; got exit=%s\n' "$label" "$want_rc" "$want_txt" "$rc"
        sed 's/^/        | /' <<<"$got"
        fail=$((fail + 1))
    fi
}

# ---- feasible outputs
check toy.out        "$IN" "$TOY/toy.out"         0 "Total 20.00 / 20"
check toy.out-wl     "$IN" "$TOY/toy.out"         0 "Computed WL   1625.0"
check partial.out    "$IN" "$TOY/partial.out"     0 "Total 16.50 / 20"
check partial.out-wl "$IN" "$TOY/partial.out"     0 "Computed WL   1620.5"
check wrongwl.out    "$IN" "$TOY/bad/wrongwl.out" 0 "Total 15.00 / 20"

# ---- gate failures (first violation)
B="$TOY/bad"
check overlap.out    "$IN" "$B/overlap.out"    1 "Overlap: S3 P1"
check outline.out    "$IN" "$B/outline.out"    1 "Out of outline: S3"
check area.out       "$IN" "$B/area.out"       1 "Area tolerance: S3 (area 21000, target 20000)"
check aspect.out     "$IN" "$B/aspect.out"     1 "Aspect ratio: S1 (w 300, h 50)"
check preplaced.out  "$IN" "$B/preplaced.out"  1 "Preplaced mismatch: P1"
check fixed.out      "$IN" "$B/fixed.out"      1 "Fixed-shape mismatch: F1"
check missing.out    "$IN" "$B/missing.out"    1 "Block count: 5 (expected 6)"
check duplicate.out  "$IN" "$B/duplicate.out"  1 "Duplicate block: S1"
check unknown.out    "$IN" "$B/unknown.out"    1 "Unknown name: F9"
check terminal.out   "$IN" "$B/terminal.out"   1 "Unknown name: T1"
check nonint.out     "$IN" "$B/nonint.out"     1 "non-integer value '150.0' for S1"
check trailing.out   "$IN" "$B/trailing.out"   1 "(5 tokens), got 6 tokens"
check count.out      "$IN" "$B/count.out"      1 "Block count: 7 (expected 6)"
check shortfile.out  "$IN" "$B/shortfile.out"  1 "found 5 block lines, expected 6"
check extraline.out  "$IN" "$B/extraline.out"  1 "extra line after 6 block lines"
check header.out     "$IN" "$B/header.out"     1 "expected 'WIRELENGTH <wl>'"
check badwl.out      "$IN" "$B/badwl.out"      1 "invalid wirelength '1.625e3'"
check zero.out       "$IN" "$B/zero.out"       1 "Non-positive size: S2 (w 0, h 100)"
check no-such-output "$IN" "$TMP/does-not-exist.out" 1 "cannot open output file"

# ---- WL tolerance boundary (|reported - WL| <= 0.01) and accepted number forms
for v in 1625.01 1624.99 1625 1625.0100000; do
    sed "s/^WIRELENGTH .*/WIRELENGTH $v/" "$TOY/toy.out" > "$TMP/wl_$v.out"
    check "wl=$v" "$IN" "$TMP/wl_$v.out" 0 "Total 20.00 / 20"
done
for v in 1625.011 1624.989; do
    sed "s/^WIRELENGTH .*/WIRELENGTH $v/" "$TOY/toy.out" > "$TMP/wl_$v.out"
    check "wl=$v" "$IN" "$TMP/wl_$v.out" 0 "Total 15.00 / 20"
done

# ---- WL ramp: T = 1600.0, WL = 1625.0 -> 5*(1.05T-WL)/(0.05T) = 3.4375
sed 's/^WLTHRESHOLD .*/WLTHRESHOLD 1600.0/' "$IN" > "$TMP/ramp.in"
check ramp           "$TMP/ramp.in" "$TOY/toy.out" 0 "WL score      3.44 / 5"
check ramp-total     "$TMP/ramp.in" "$TOY/toy.out" 0 "Total 18.44 / 20"
# WL > 1.05T -> 0 (T = 1500.0, 1.05T = 1575 < 1625)
sed 's/^WLTHRESHOLD .*/WLTHRESHOLD 1500.0/' "$IN" > "$TMP/over.in"
check over-1.05T     "$TMP/over.in" "$TOY/toy.out" 0 "Total 15.00 / 20"

# ---- malformed input -> exit 2
badin() {  # <label> <sed expr> <expected substring>
    sed "$2" "$IN" > "$TMP/$1.in"
    check "in:$1" "$TMP/$1.in" "$TOY/toy.out" 2 "$3"
}
badin threshold  's/^WLTHRESHOLD .*/WLTHRESHOLD 1700/'          "WLTHRESHOLD must be"
badin threshold2 's/^WLTHRESHOLD .*/WLTHRESHOLD 1700.3/'        "WLTHRESHOLD must be"
badin smallarea  's/^SOFT S3 20000$/SOFT S3 9999/'              "area must be >= 10000"
badin dupname    's/^SOFT S2 15000$/SOFT S1 15000/'             "duplicate name 'S1'"
badin termdup    's/^TERMINAL T2 /TERMINAL S1 /'                "duplicate name 'S1'"
badin termpos    's/^TERMINAL T2 400 50$/TERMINAL T2 399 50/'   "not on the outline border"
badin preout     's/^PREPLACED P1 .*/PREPLACED P1 100 100 350 100/' "not inside the outline"
badin presedge   's/^PREPLACED P1 .*/PREPLACED P1 100 100 0 100/' "touches an outline edge"
badin resname    's/^SOFT S2 15000$/SOFT LEFT 15000/'           "is reserved"
badin selfconn   's/^CONN S1 S2 2$/CONN S1 S1 2/'               "connected to itself"
badin ttconn     's/^CONN S1 S2 2$/CONN T1 T2 2/'               "terminal-terminal"
badin weight     's/^CONN S1 S2 2$/CONN S1 S2 11/'              "weight must be in [1,10]"
badin duppair    's/^CONN S1 P1 1$/CONN S2 S1 1/'               "appears more than once"
badin mibarea    's/^SOFT S2 15000$/SOFT S2 15001/'             "different area target"
badin mibfixed   's/^MIB 2 S1 S2$/MIB 2 S1 F1/'                 "not a soft block"
badin mibk       's/^MIB 2 S1 S2$/MIB 1 S1 S2/'                 "k must be >= 2"
badin opposite   's/^BOUNDARY F1 LEFT BOTTOM$/BOUNDARY F1 LEFT RIGHT/' "must be adjacent"
badin bpre       's/^BOUNDARY S3 RIGHT$/BOUNDARY P1 RIGHT/'     "is a preplaced block"
badin bedge      's/^BOUNDARY S3 RIGHT$/BOUNDARY S3 NORTH/'     "invalid EDGE 'NORTH'"
badin count      's/^CONNECTIONS 7$/CONNECTIONS 8/'             "expected 'CONN', got 'MIBS'"
badin trailing   '$a EXTRA'                                     "trailing token"
badin fixedw     's/^FIXED F1 150 100$/FIXED F1 0 100/'           "FIXED F1: width and height must be positive"
badin nonint     's/^OUTLINE 400 300$/OUTLINE 400 3e2/'         "expected integer outline H"
check in:missing "$TMP/nope.in" "$TOY/toy.out" 2 "cannot open input file"

echo
echo "passed $pass, failed $fail"
[[ $fail -eq 0 ]]
