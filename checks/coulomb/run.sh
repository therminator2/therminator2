#!/bin/bash
# Builds and runs the full-evolution checks of the Coulomb afterburner (coulomb_evolution_checks.cxx).
# Needs ROOT: root-config in PATH (e.g. an activated conda environment), or ROOTSYS set.
# Build products go to $WORK (default: a new temporary directory), graphs of the trajectories to $OUT
# (default: coulomb_checks.root in the current directory).
#
# Spectators are switched off in the checks program through the configuration (CoulombUseSpectators = 0),
# so that the test events contain only the particles defined in the checks.
#
# Works both when executed (./run.sh) and when sourced (. run.sh): the body runs in a subshell, so that
# `exit`, `set -e` and `cd` do not affect the calling shell (sourcing used to close the terminal).
(
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
TOP=$(cd "$HERE/../.." && pwd)
WORK=${WORK:-$(mktemp -d)}
mkdir -p "$WORK"
# Graphs of the trajectories, code vs exact (default: in the directory the script is called from)
OUT=${OUT:-$PWD/coulomb_checks.root}
case "$OUT" in /*) ;; *) OUT="$PWD/$OUT" ;; esac

if ! command -v root-config >/dev/null; then
    if [ -n "${ROOTSYS:-}" ]; then export PATH="$ROOTSYS/bin:$PATH"; else echo "root-config not found, set ROOTSYS"; exit 1; fi
fi
export LD_LIBRARY_PATH="$(root-config --libdir):${LD_LIBRARY_PATH:-}"

CXXFLAGS="-std=c++14 -O2 -g -Wno-deprecated -I $TOP/build/include -D_DEBUG_LEVEL_=0 $(root-config --cflags)"

echo "Building in $WORK"
OBJS=()
for SRC in "$HERE/coulomb_evolution_checks.cxx" \
           $(for f in CoulombAfterburner AbstractAfterburner Particle ParticleCoor ParticleType ParticleDB DecayTable DecayChannel \
                      Event Crc32 TGlauberMC Configurator Parser Hypersurface_Library Vector3D; do
                 echo "$TOP/build/src/$f.cxx"; done); do
    OBJ="$WORK/$(basename "${SRC%.cxx}").o"
    g++ -c "$SRC" -o "$OBJ" $CXXFLAGS
    OBJS+=("$OBJ")
done
g++ -o "$WORK/coulomb_evolution_checks" "${OBJS[@]}" $(root-config --libs)

cd "$WORK"
# The afterburner prints progress messages; show only the results of the checks
./coulomb_evolution_checks "$OUT" "$TOP" | grep -E '^(PASS|FAIL|INFO|ALL PASSED|FAILED)'
exit "${PIPESTATUS[0]}"
)
