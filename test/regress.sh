#!/bin/sh
# Correctness regression tests.  Each case runs the simulator and checks a
# property that a refactor could plausibly break.  Exits non-zero on failure.

set -u
BIN="${CARNAVAL:-$(dirname "$0")/../bin/carnaval}"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
fail=0
pass=0

check() {  # check <description> <condition-exit-code>
  if [ "$2" -eq 0 ]; then
    pass=$((pass+1)); printf '  ok    %s\n' "$1"
  else
    fail=$((fail+1)); printf '  FAIL  %s\n' "$1"
  fi
}

TRNA=GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA
HP=GGGGAAAACCCC

echo "carnaval regression tests ($BIN)"

# --- the simulator runs and conserves board validity -----------------------
# assertValid() runs inside toJson(), so a successful --save means the board,
# its cell occupancy and all its chain links are mutually consistent.
for args in \
  "-x 64 -y 64 --init $HP" \
  "-x 128 -y 128 --init $TRNA" \
  "-x 64 -y 64 --init AACCUUGG --density 0.2 --bond 1" \
  "-x 64 -y 64 --init AACCUUGG --density 0.2 --bond 0" \
  "-x 32 -y 32 -z 16 --surface z --tether 4 --init GGCCAATTGC --dna" \
  "-x 64 -y 64 --init $HP --pin-first" \
  "-x 64 -y 64 --init $HP --sort-units" \
  "-x 64 -y 64 --init $HP --legal-moves-only" \
  "-x 48 -y 48 --init $HP" ; do
  # shellcheck disable=SC2086
  $BIN $args --unit-moves 20000 --rnd 5 --period 100000000 --save "$TMP/v.json" >/dev/null 2>&1
  check "board stays valid: $args" $?
done

# --- determinism ------------------------------------------------------------
$BIN -x 64 -y 64 --init "$HP" --unit-moves 50000 --rnd 11 --folds --period 1000 --monochrome >"$TMP/a" 2>/dev/null
$BIN -x 64 -y 64 --init "$HP" --unit-moves 50000 --rnd 11 --folds --period 1000 --monochrome >"$TMP/b" 2>/dev/null
cmp -s "$TMP/a" "$TMP/b"; check "same seed reproduces the same trajectory" $?
$BIN -x 64 -y 64 --init "$HP" --unit-moves 50000 --rnd 12 --folds --period 1000 --monochrome >"$TMP/c" 2>/dev/null
if cmp -s "$TMP/a" "$TMP/c"; then check "different seed gives a different trajectory" 1
else check "different seed gives a different trajectory" 0; fi

# --- JSON round trip --------------------------------------------------------
$BIN -x 32 -y 32 -z 16 --surface z --tether 4 --init GGCCAATTGC --dna \
     --unit-moves 2000 --rnd 1 --period 100000000 --save "$TMP/t1.json" >/dev/null 2>&1
$BIN --load "$TMP/t1.json" --unit-moves 0 --period 100000000 --save "$TMP/t2.json" >/dev/null 2>&1
check "tethered board reloads" $?
python3 -I - "$TMP/t1.json" "$TMP/t2.json" <<'PY'
import json, sys
a = json.load(open(sys.argv[1])); b = json.load(open(sys.argv[2]))
def summary(j):
    u = j["unit"]
    return (len(u),
            sum(1 for x in u if x.get("inert")),
            sum(1 for x in u if x.get("pinned")),
            sum(1 for x in u if "next" in x),
            [x["pos"] for x in u])
sys.exit(0 if summary(a) == summary(b) else 1)
PY
check "round trip preserves units, surface sites, pinning, links, positions" $?

# --- bondProb actually gates ligation --------------------------------------
links() {
  $BIN -x 64 -y 64 --init AACCUUGG --density 0.1 --bond "$1" --unit-moves 50000 \
       --rnd 3 --period 100000000 --save "$TMP/b$1.json" >/dev/null 2>&1
  python3 -I -c "import json,sys; print(sum(1 for x in json.load(open(sys.argv[1]))['unit'] if 'next' in x))" "$TMP/b$1.json"
}
L0=$(links 0); L1=$(links 1)
# with bond 0 the only covalent links are the 7 in the 8-base seed template
[ "$L0" -eq 7 ]; check "--bond 0 forms no new covalent bonds (got $L0, template has 7)" $?
[ "$L1" -gt "$L0" ]; check "--bond 1 forms new covalent bonds ($L1 > $L0)" $?

# --- surface sites never move ----------------------------------------------
$BIN -x 32 -y 32 -z 16 --surface z --tether 4 --init GGCCAATTGC --unit-moves 5000 \
     --rnd 2 --period 100000000 --save "$TMP/s.json" >/dev/null 2>&1
python3 -I -c "
import json,sys
u=json.load(open(sys.argv[1]))['unit']
bad=[x for x in u if x.get('inert') and x['pos'][2]!=0]
sys.exit(0 if not bad else 1)" "$TMP/s.json"
check "inert surface sites stay in their plane" $?

# --- pinning ----------------------------------------------------------------
$BIN -x 64 -y 64 --init "$HP" --pin-first --unit-moves 100000 --rnd 1 \
     --period 100000000 --save "$TMP/p.json" >/dev/null 2>&1
python3 -I -c "
import json,sys
u=json.load(open(sys.argv[1]))['unit']
sys.exit(0 if u[0]['pos']==[0,0,0] and u[0].get('pinned') else 1)" "$TMP/p.json"
check "--pin-first holds unit 0 at its lattice site" $?

# --- DNA mode ---------------------------------------------------------------
$BIN -x 64 -y 64 --dna --init GGCCAATTGC --unit-moves 100 --rnd 1 \
     --period 100000000 --save "$TMP/d.json" >/dev/null 2>&1
python3 -I -c "
import json,sys
j=json.load(open(sys.argv[1]))
bases={x['base'] for x in j['unit']}
sys.exit(0 if 't' in bases and 'u' not in bases and j['params']['wobble'] is False else 1)" "$TMP/d.json"
check "--dna emits t, not u, and disables wobble pairs" $?

# --- constructs used by the folding-time test are self-consistent ----------
python3 -I "$(dirname "$0")/constructs.py" >/dev/null 2>&1
check "reference constructs have complementary, well-formed target pairs" $?

echo "$pass passed, $fail failed"
[ "$fail" -eq 0 ]
