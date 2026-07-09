# Reproducibility Notebook Setup

This folder contains a self-contained reproducibility notebook for processing
packaged MatterGen outputs and generating the downstream analysis artifacts.

Notebook:

```text
mg_generation_postprocessing_pipeline.ipynb
```

Setup script:

```text
setup_reproducibility_envs.sh
```

The notebook needs three Python environments:

- `vsbtools`: the Jupyter/kernel environment.
- `scout-matter`: the MatterGen fork used for guidance descriptors and losses.
- `grace`: a `tensorpotential` environment used by the GRACE bridge.

## Recommended Contained Setup

From a machine with `git`, Python, and internet access:

```bash
bash setup_reproducibility_envs.sh --root ./vsbtools_reproducibility_env
```

The script creates everything under `./vsbtools_reproducibility_env`:

```text
src/vsbtools
src/scout-matter
venvs/vsbtools
venvs/scout-matter
venvs/grace
state/
work/mg_generation_postprocessing_pipeline.ipynb
```

It also keeps local runtime state inside `state/`, including Jupyter, IPython,
matplotlib, pip cache, and `vsbtools` external-path configuration. It does not
install a user/global Jupyter kernel and does not write to `~/.config/vsbtools`.

After installation, the script launches JupyterLab automatically. To install
without launching:

```bash
bash setup_reproducibility_envs.sh --root ./vsbtools_reproducibility_env --no-launch
```

Launch later with:

```bash
./vsbtools_reproducibility_env/run_reproducibility_notebook.sh
```

For a fixed review run, pin repository refs:

```bash
bash setup_reproducibility_envs.sh \
  --root ./vsbtools_reproducibility_env \
  --vsbtools-ref <commit-or-tag> \
  --scout-matter-ref <commit-or-tag>
```

Use `--force` to recreate the contained workspace:

```bash
bash setup_reproducibility_envs.sh --root ./vsbtools_reproducibility_env --force
```

## Manual Configuration

If you already have the three environments, open the notebook in a `vsbtools`
kernel and set the first code cell:

```python
MATTER_SCOUT_PATH = Path("/path/to/scout-matter-or-venv")
GRACE_ENV_OR_PYTHON = Path("/path/to/grace-venv-or-bin-python")
```

`MATTER_SCOUT_PATH` may point to:

- a `scout-matter` source tree,
- a `scout-matter` virtual environment root,
- that virtual environment's Python executable,
- or an import root/site-packages directory containing `mattergen`.

`GRACE_ENV_OR_PYTHON` may point to:

- a GRACE/tensorpotential virtual environment root,
- or that virtual environment's `bin/python`.

If these notebook variables are `None`, `vsbtools` falls back to:

```text
MATTERGEN_PYTHON_PATH
GRACE_PYTHON
~/.config/vsbtools/external_paths.json
```

The contained setup script avoids those global fallbacks by exporting explicit
paths only inside its launcher environment.

## Expected Outputs

The notebook writes derived files under:

```text
reproducibility_run/
```

inside the directory where the notebook is run. The contained setup copies the
notebook to:

```text
vsbtools_reproducibility_env/work/
```

so generated outputs stay inside the contained workspace.

