#!/bin/bash
# Builds and runs the full-evolution checks of the Coulomb afterburner (coulomb_evolution_checks.cxx).
# Needs ROOT: root-config in PATH (e.g. an activated conda environment), or ROOTSYS set.
# Build products go to $WORK (default: a new temporary directory).
#
# The afterburner is compiled from a copy of build/src/CoulombAfterburner.cxx with USE_SPECTATORS = false,
# so that the test events contain only the particles defined in the checks. Nothing in the source tree is changed.
#
# Works both when executed (./run.sh) and when sourced (. run.sh): the body runs in a subshell, so that
# `exit`, `set -e` and `cd` do not affect the calling shell (sourcing used to close the terminal).
(
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
TOP=$(cd "$HERE/../.." && pwd)
WORK=${WORK:-$(mktemp -d)}
mkdir -p "$WORK"

if ! command -v root-config >/dev/null; then
    if [ -n "${ROOTSYS:-}" ]; then export PATH="$ROOTSYS/bin:$PATH"; else echo "root-config not found, set ROOTSYS"; exit 1; fi
fi
export LD_LIBRARY_PATH="$(root-config --libdir):${LD_LIBRARY_PATH:-}"

CXXFLAGS="-std=c++14 -O2 -g -Wno-deprecated -I $TOP/build/include -D_DEBUG_LEVEL_=0 $(root-config --cflags)"

sed 's/^static const bool USE_SPECTATORS = true;/static const bool USE_SPECTATORS = false;/' \
    "$TOP/build/src/CoulombAfterburner.cxx" > "$WORK/CoulombAfterburner_nospec.cxx"
grep -q '^static const bool USE_SPECTATORS = false;' "$WORK/CoulombAfterburner_nospec.cxx" \
    || { echo "Could not switch off spectators in the copy of CoulombAfterburner.cxx"; exit 1; }

echo "Building in $WORK"
OBJS=()
for SRC in "$WORK/CoulombAfterburner_nospec.cxx" "$HERE/coulomb_evolution_checks.cxx" \
           $(for f in AbstractAfterburner Particle ParticleCoor ParticleType ParticleDB DecayTable DecayChannel Event Crc32 TGlauberMC; do
                 echo "$TOP/build/src/$f.cxx"; done); do
    OBJ="$WORK/$(basename "${SRC%.cxx}").o"
    g++ -c "$SRC" -o "$OBJ" $CXXFLAGS
    OBJS+=("$OBJ")
done
g++ -o "$WORK/coulomb_evolution_checks" "${OBJS[@]}" $(root-config --libs)

cd "$WORK"
# The afterburner prints progress messages; show only the results of the checks
./coulomb_evolution_checks | grep -E '^(PASS|FAIL|INFO|ALL PASSED|FAILED)'
exit "${PIPESTATUS[0]}"
)
