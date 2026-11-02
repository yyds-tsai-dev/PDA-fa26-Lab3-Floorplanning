#!/usr/bin/env bash
# Grade a Lab3 binary on one or more cases.
#   bash tools/grade.sh <Lab3_binary> <verifier> <case.in>...
#
# For each case: runs `timeout 120 <Lab3> <case.in> <out>`, measures wall and
# user+sys time, flags MULTITHREAD? when (user+sys) > 1.2*wall and wall > 2 s,
# then runs the verifier and takes its "Total x.xx / 20". A timeout or a wall
# time above 120 s scores 0. Outputs and logs are kept in a temp directory.
# Always exits 0 after printing the table (a grading run, not a test).
set -u

LIMIT="${LAB3_TIME_LIMIT:-120}"   # override only for testing this script

if [[ $# -lt 3 ]]; then
    echo "Usage: $0 <Lab3_binary> <verifier> <case.in>..." >&2
    exit 2
fi
BIN="$1"; VER="$2"; shift 2
for f in "$BIN" "$VER"; do
    if [[ ! -x "$f" ]]; then echo "not executable: $f" >&2; exit 2; fi
done
BIN="$(cd "$(dirname "$BIN")" && pwd)/$(basename "$BIN")"
VER="$(cd "$(dirname "$VER")" && pwd)/$(basename "$VER")"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/lab3_grade.XXXXXX")"
echo "work dir: $WORK"

USE_GNU_TIME=0
if [[ -x /usr/bin/time ]] && /usr/bin/time -f '%e' -o /dev/null true 2>/dev/null; then
    USE_GNU_TIME=1
else
    echo "note: /usr/bin/time (GNU) not available; using bash 'time'"
fi

total=0
rows=()
for IN in "$@"; do
    idx=$((${idx:-0} + 1))
    name="$(basename "$IN")"; name="${name%.in}"
    OUT="$WORK/$idx.$name.out"
    TF="$WORK/$idx.$name.time"
    LOG="$WORK/$idx.$name.run.log"
    VLOG="$WORK/$idx.$name.verify.log"
    rm -f "$OUT"
    if [[ ! -r "$IN" ]]; then
        rows+=("$(printf '%-14s %8s %8s %-12s %-10s %6s  %s' "$name" - - - ERROR 0.00 "cannot read $IN")")
        continue
    fi

    ABS_IN="$(cd "$(dirname "$IN")" && pwd)/$(basename "$IN")"
    # Run the solver from the work dir so stray files do not pollute the repo.
    if [[ $USE_GNU_TIME -eq 1 ]]; then
        ( cd "$WORK" && /usr/bin/time -f '%e %U %S' -o "$TF" \
              timeout "$LIMIT" "$BIN" "$ABS_IN" "$OUT" >"$LOG" 2>&1 )
        rc=$?
    else
        ( cd "$WORK"; TIMEFORMAT='%R %U %S'
          { time timeout "$LIMIT" "$BIN" "$ABS_IN" "$OUT" >"$LOG" 2>&1 ; } 2>"$TF" )
        rc=$?
    fi
    # GNU time may prefix "Command exited with non-zero status N"; take the last line.
    read -r wall usr sys < <(tail -n 1 "$TF" 2>/dev/null)
    wall="${wall:-0}"; usr="${usr:-0}"; sys="${sys:-0}"
    cpu="$(awk -v u="$usr" -v s="$sys" 'BEGIN{printf "%.2f", u+s}')"

    flag="-"
    if awk -v c="$cpu" -v w="$wall" 'BEGIN{exit !(c > 1.2*w && w > 2)}'; then
        flag="MULTITHREAD?"
    fi

    status=""; score="0.00"; note=""
    if [[ $rc -eq 124 || $rc -eq 137 ]] || awk -v w="$wall" -v l="$LIMIT" 'BEGIN{exit !(w > l)}'; then
        status="TIMEOUT"; note="runtime > ${LIMIT}s"
    else
        [[ $rc -ne 0 ]] && note="Lab3 exit $rc; "
        "$VER" "$IN" "$OUT" >"$VLOG" 2>&1
        vrc=$?
        score="$(sed -n 's/^Total \([0-9.]*\) \/ 20$/\1/p' "$VLOG" | tail -n 1)"
        score="${score:-0.00}"
        case $vrc in
            0) status="FEASIBLE" ;;
            1) status="INFEASIBLE"
               note+="$(grep -m1 '^Result:' "$VLOG" | sed 's/^Result: INFEASIBLE //')" ;;
            *) status="VERIFY-ERR"; score="0.00"; note+="$(head -n 1 "$VLOG")" ;;
        esac
    fi
    total="$(awk -v t="$total" -v s="$score" 'BEGIN{printf "%.2f", t+s}')"
    rows+=("$(printf '%-14s %8.2f %8.2f %-12s %-10s %6s  %s' "$name" "$wall" "$cpu" "$flag" "$status" "$score" "$note")")
done

echo
printf '%-14s %8s %8s %-12s %-10s %6s  %s\n' CASE 'WALL(s)' 'CPU(s)' THREADS STATUS SCORE NOTE
printf '%s\n' "${rows[@]}"
echo
printf 'TOTAL %.2f / %d\n' "$total" $((20 * $#))
exit 0
