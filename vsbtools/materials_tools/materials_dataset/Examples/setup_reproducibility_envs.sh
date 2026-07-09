#!/usr/bin/env bash
set -euo pipefail

VSBTOOLS_REPO_URL="${VSBTOOLS_REPO_URL:-https://github.com/v-baturin/vsbtools.git}"
SCOUT_MATTER_REPO_URL="${SCOUT_MATTER_REPO_URL:-https://github.com/link-lab3629/scout-matter.git}"
VSBTOOLS_REF="${VSBTOOLS_REF:-main}"
SCOUT_MATTER_REF="${SCOUT_MATTER_REF:-main}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
ROOT="${VSBTOOLS_REPRO_ROOT:-$PWD/vsbtools_reproducibility_env}"
LAUNCH=1
FORCE=0
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
    cat <<'USAGE'
Usage: setup_reproducibility_envs.sh [options]

Creates a contained reproducibility workspace with three virtual environments:
  1. vsbtools notebook/kernel environment from https://github.com/v-baturin/vsbtools
  2. scout-matter/MatterGen environment from https://github.com/link-lab3629/scout-matter
  3. GRACE/tensorpotential environment

Options:
  --root PATH               Workspace root. Default: ./vsbtools_reproducibility_env
  --python PATH             Python used to create venvs. Default: python3
  --vsbtools-ref REF        Git ref for vsbtools. Default: main
  --scout-matter-ref REF    Git ref for scout-matter. Default: main
  --no-launch               Install/configure only; do not launch JupyterLab
  --force                   Recreate src/, venvs/, state/, and work/ under --root
  -h, --help                Show this help

Environment overrides:
  VSBTOOLS_REPO_URL         Default: https://github.com/v-baturin/vsbtools.git
  SCOUT_MATTER_REPO_URL     Default: https://github.com/link-lab3629/scout-matter.git
  VSBTOOLS_REF              Default: main
  SCOUT_MATTER_REF          Default: main
  PYTHON_BIN                Default: python3
  VSBTOOLS_REPRO_ROOT       Default: ./vsbtools_reproducibility_env

All Jupyter, IPython, matplotlib, pip, and vsbtools external-path state is kept
inside --root/state. The script does not write to ~/.config/vsbtools or install
a user/global Jupyter kernel.
USAGE
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --root)
            ROOT="$2"
            shift 2
            ;;
        --python)
            PYTHON_BIN="$2"
            shift 2
            ;;
        --vsbtools-ref)
            VSBTOOLS_REF="$2"
            shift 2
            ;;
        --scout-matter-ref)
            SCOUT_MATTER_REF="$2"
            shift 2
            ;;
        --no-launch)
            LAUNCH=0
            shift
            ;;
        --force)
            FORCE=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown argument: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

log() {
    printf '\n[%s] %s\n' "$(date '+%H:%M:%S')" "$*"
}

require_cmd() {
    if ! command -v "$1" >/dev/null 2>&1; then
        echo "Required command not found: $1" >&2
        exit 1
    fi
}

clone_or_checkout() {
    local url="$1"
    local ref="$2"
    local dest="$3"
    if [[ -d "$dest/.git" ]]; then
        log "Updating $dest"
        git -C "$dest" fetch --tags --prune origin
    elif [[ -e "$dest" ]]; then
        echo "Path exists but is not a git checkout: $dest" >&2
        exit 1
    else
        log "Cloning $url into $dest"
        git clone "$url" "$dest"
    fi
    git -C "$dest" checkout "$ref"
}

make_venv() {
    local venv="$1"
    if [[ ! -x "$venv/bin/python" ]]; then
        log "Creating venv $venv"
        "$PYTHON_BIN" -m venv "$venv"
    fi
    "$venv/bin/python" -m pip install --upgrade pip setuptools wheel
}

python_site_packages() {
    "$1" - <<'PY'
import sysconfig
print(sysconfig.get_paths()["purelib"])
PY
}

ROOT="$(mkdir -p "$ROOT" && cd "$ROOT" && pwd)"

if [[ "$FORCE" -eq 1 ]]; then
    log "Removing contained workspace directories under $ROOT"
    rm -rf "$ROOT/src" "$ROOT/venvs" "$ROOT/state" "$ROOT/work"
fi

require_cmd git
require_cmd "$PYTHON_BIN"

SRC_DIR="$ROOT/src"
VENVS_DIR="$ROOT/venvs"
STATE_DIR="$ROOT/state"
WORK_DIR="$ROOT/work"
mkdir -p "$SRC_DIR" "$VENVS_DIR" "$STATE_DIR" "$WORK_DIR"

export XDG_CONFIG_HOME="$STATE_DIR/xdg_config"
export XDG_CACHE_HOME="$STATE_DIR/xdg_cache"
export JUPYTER_CONFIG_DIR="$STATE_DIR/jupyter_config"
export JUPYTER_DATA_DIR="$STATE_DIR/jupyter_data"
export JUPYTER_RUNTIME_DIR="$STATE_DIR/jupyter_runtime"
export IPYTHONDIR="$STATE_DIR/ipython"
export MPLCONFIGDIR="$STATE_DIR/matplotlib"
export PIP_CACHE_DIR="$STATE_DIR/pip_cache"
export VSBTOOLS_EXTERNAL_PATHS_CONFIG="$STATE_DIR/vsbtools_external_paths.json"
mkdir -p "$XDG_CONFIG_HOME" "$XDG_CACHE_HOME" "$JUPYTER_CONFIG_DIR" \
    "$JUPYTER_DATA_DIR" "$JUPYTER_RUNTIME_DIR" "$IPYTHONDIR" "$MPLCONFIGDIR" "$PIP_CACHE_DIR"

VSBTOOLS_SRC="$SRC_DIR/vsbtools"
SCOUT_SRC="$SRC_DIR/scout-matter"
VSBTOOLS_VENV="$VENVS_DIR/vsbtools"
SCOUT_VENV="$VENVS_DIR/scout-matter"
GRACE_VENV="$VENVS_DIR/grace"

clone_or_checkout "$VSBTOOLS_REPO_URL" "$VSBTOOLS_REF" "$VSBTOOLS_SRC"
clone_or_checkout "$SCOUT_MATTER_REPO_URL" "$SCOUT_MATTER_REF" "$SCOUT_SRC"
VSBTOOLS_COMMIT="$(git -C "$VSBTOOLS_SRC" rev-parse HEAD)"
SCOUT_MATTER_COMMIT="$(git -C "$SCOUT_SRC" rev-parse HEAD)"

make_venv "$VSBTOOLS_VENV"
make_venv "$SCOUT_VENV"
make_venv "$GRACE_VENV"

log "Installing scout-matter into its contained venv"
if ! "$SCOUT_VENV/bin/python" -m pip install "$SCOUT_SRC"; then
    log "Non-editable scout-matter install failed; trying editable install"
    "$SCOUT_VENV/bin/python" -m pip install -e "$SCOUT_SRC"
fi
"$SCOUT_VENV/bin/python" - <<'PY'
import mattergen
print("mattergen:", mattergen.__file__)
PY
SCOUT_SITE_PACKAGES="$(python_site_packages "$SCOUT_VENV/bin/python")"

log "Installing GRACE/tensorpotential into its contained venv"
"$GRACE_VENV/bin/python" -m pip install "ase<3.26" tensorpotential
"$GRACE_VENV/bin/python" - <<'PY'
import tensorpotential.calculator
print("tensorpotential.calculator import OK")
PY

log "Installing vsbtools notebook/kernel environment"
"$VSBTOOLS_VENV/bin/python" -m pip install -e "$VSBTOOLS_SRC"
"$VSBTOOLS_VENV/bin/python" -m pip install \
    "cclib>=1.8" \
    dscribe \
    ijson \
    jupyterlab \
    ipykernel \
    matplotlib \
    mp-api \
    networkx \
    pandas \
    phonopy \
    prettytable \
    pymatgen \
    PyYAML \
    scipy \
    torch
if "$VSBTOOLS_VENV/bin/python" -m pip install mysqlclient; then
    log "Optional mysqlclient dependency installed"
else
    log "Optional mysqlclient dependency failed to install; continuing because this notebook does not use local OQMD/MySQL"
fi

export MATTERGEN_PYTHON_PATH="$SCOUT_SITE_PACKAGES"
export GRACE_PYTHON="$GRACE_VENV/bin/python"

log "Validating bridge imports from the vsbtools venv"
"$VSBTOOLS_VENV/bin/python" - <<'PY'
import os
import sys
from pathlib import Path

sys.path.append(os.environ["MATTERGEN_PYTHON_PATH"])
import mattergen.common.data.chemgraph
import mattergen.diffusion.diffusion_loss

import vsbtools
print("vsbtools:", Path(vsbtools.__file__).resolve())
print("MATTERGEN_PYTHON_PATH:", os.environ["MATTERGEN_PYTHON_PATH"])
print("GRACE_PYTHON:", os.environ["GRACE_PYTHON"])
PY

NOTEBOOK_SRC="$VSBTOOLS_SRC/vsbtools/materials_tools/materials_dataset/Examples/mg_generation_postprocessing_pipeline.ipynb"
if [[ ! -f "$NOTEBOOK_SRC" && -f "$SCRIPT_DIR/mg_generation_postprocessing_pipeline.ipynb" ]]; then
    NOTEBOOK_SRC="$SCRIPT_DIR/mg_generation_postprocessing_pipeline.ipynb"
fi
if [[ ! -f "$NOTEBOOK_SRC" ]]; then
    echo "Could not find mg_generation_postprocessing_pipeline.ipynb in cloned vsbtools or beside this script." >&2
    exit 1
fi
NOTEBOOK_DST="$WORK_DIR/mg_generation_postprocessing_pipeline.ipynb"
cp "$NOTEBOOK_SRC" "$NOTEBOOK_DST"

KERNEL_PREFIX="$STATE_DIR/jupyter_kernel_prefix"
"$VSBTOOLS_VENV/bin/python" -m ipykernel install \
    --prefix "$KERNEL_PREFIX" \
    --name vsbtools-repro \
    --display-name "vsbtools reproducibility"

ENV_FILE="$ROOT/reproducibility_env.sh"
cat > "$ENV_FILE" <<EOF
#!/usr/bin/env bash
export VSBTOOLS_REPRO_ROOT="$ROOT"
export VSBTOOLS_SRC="$VSBTOOLS_SRC"
export SCOUT_MATTER_SRC="$SCOUT_SRC"
export VSBTOOLS_COMMIT="$VSBTOOLS_COMMIT"
export SCOUT_MATTER_COMMIT="$SCOUT_MATTER_COMMIT"
export VSBTOOLS_VENV="$VSBTOOLS_VENV"
export SCOUT_MATTER_VENV="$SCOUT_VENV"
export GRACE_VENV="$GRACE_VENV"
export MATTERGEN_PYTHON_PATH="$MATTERGEN_PYTHON_PATH"
export GRACE_PYTHON="$GRACE_PYTHON"
export XDG_CONFIG_HOME="$XDG_CONFIG_HOME"
export XDG_CACHE_HOME="$XDG_CACHE_HOME"
export JUPYTER_CONFIG_DIR="$JUPYTER_CONFIG_DIR"
export JUPYTER_DATA_DIR="$JUPYTER_DATA_DIR"
export JUPYTER_RUNTIME_DIR="$JUPYTER_RUNTIME_DIR"
export JUPYTER_PATH="$KERNEL_PREFIX/share/jupyter"
export IPYTHONDIR="$IPYTHONDIR"
export MPLCONFIGDIR="$MPLCONFIGDIR"
export PIP_CACHE_DIR="$PIP_CACHE_DIR"
export VSBTOOLS_EXTERNAL_PATHS_CONFIG="$VSBTOOLS_EXTERNAL_PATHS_CONFIG"
EOF
chmod +x "$ENV_FILE"

SETUP_MANIFEST="$ROOT/setup_manifest.json"
cat > "$SETUP_MANIFEST" <<EOF
{
  "vsbtools_repo_url": "$VSBTOOLS_REPO_URL",
  "vsbtools_ref": "$VSBTOOLS_REF",
  "vsbtools_commit": "$VSBTOOLS_COMMIT",
  "scout_matter_repo_url": "$SCOUT_MATTER_REPO_URL",
  "scout_matter_ref": "$SCOUT_MATTER_REF",
  "scout_matter_commit": "$SCOUT_MATTER_COMMIT",
  "root": "$ROOT",
  "vsbtools_venv": "$VSBTOOLS_VENV",
  "scout_matter_venv": "$SCOUT_VENV",
  "grace_venv": "$GRACE_VENV",
  "mattergen_python_path": "$MATTERGEN_PYTHON_PATH",
  "grace_python": "$GRACE_PYTHON"
}
EOF

LAUNCHER="$ROOT/run_reproducibility_notebook.sh"
cat > "$LAUNCHER" <<EOF
#!/usr/bin/env bash
set -euo pipefail
source "$ENV_FILE"
cd "$WORK_DIR"
exec "$VSBTOOLS_VENV/bin/python" -m jupyter lab \\
    --notebook-dir "$WORK_DIR" \\
    "$NOTEBOOK_DST"
EOF
chmod +x "$LAUNCHER"

log "Contained reproducibility environment is ready"
cat <<EOF

Workspace:
  $ROOT

Virtual environments:
  vsbtools:      $VSBTOOLS_VENV
  scout-matter:  $SCOUT_VENV
  GRACE:         $GRACE_VENV

Notebook copy:
  $NOTEBOOK_DST

Environment file:
  $ENV_FILE

Setup manifest:
  $SETUP_MANIFEST

Launcher:
  $LAUNCHER

EOF

if [[ "$LAUNCH" -eq 1 ]]; then
    log "Launching JupyterLab"
    exec "$LAUNCHER"
else
    log "Launch skipped. Run this later:"
    echo "  $LAUNCHER"
fi
