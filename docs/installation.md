# Installation

Requires [conda](https://docs.conda.io/) (miniconda or anaconda) on `PATH`.

This project uses **separate, independent conda environments** - none is
created "inside" another, and they serve different purposes:

| Environment | Created from | Do you activate it? | Required? |
|---|---|---|---|
| `pypkatool` | `environment.yml` | **Yes - every time** you run the tool: `conda activate pypkatool` | Always |
| `py27` | `environment-py27.yml` | **No, never.** `pypkatool` finds it on disk automatically at runtime. | Always |
| `pdbfixer` | `environment-pdbfixer.yml` | **No, never.** `pypkatool` finds it on disk automatically at runtime. | Only for the `fixstructure` command |

Why `py27` exists: PyPKA's compiled Poisson-Boltzmann backend (DelPhi4py)
needs `numpy<2` and `libgfortran4=7.5.0` in the main environment (the
compiled binary was built against NumPy's 1.x C API and links `GFORTRAN_7`,
which `libgfortran.so.5` / GFortran 10+ doesn't export) - that part lives in
`pypkatool`. Separately, PyPKA also shells out to a bare `python2.7`
interpreter internally (via `pdbmender`'s vendored `pdb2pqr.py`), which
cannot coexist with Python 3 in the same environment, hence the second,
minimal, interpreter-only environment.

Why `pdbfixer` exists: the optional `fixstructure` command (see
{doc}`fixstructure`) uses PDBFixer/OpenMM, which require `numpy>=2` - the
opposite pin from `pypkatool`'s own `numpy<2` requirement. The two cannot
share an environment, so this is a third, separate one. It is only needed
if you plan to use `fixstructure`.

## One-time setup

```bash
# 1. Create the required environments (pdbfixer is optional - see below).
#    Order doesn't matter, and it doesn't matter whether any environment is
#    currently active: `conda env create` always builds a new, independent
#    environment from scratch.
conda env create -f environment.yml
conda env create -f environment-py27.yml
conda env create -f environment-pdbfixer.yml   # optional, only for fixstructure

# 2. Activate the MAIN environment (py27 is never activated) and install
#    pypkatool into it - pick ONE of the two:
conda activate pypkatool

# 2a. Editable install (recommended if you cloned this repo to modify or
#     update it): the command reads pypkatool/ from this checkout directly,
#     so `git pull` and local edits take effect immediately, no reinstall.
pip install -e .

# 2b. Regular install (recommended if you just want to use the CLI and
#     don't plan to touch the source): copies the package into the
#     environment's site-packages, same as any other pip package.
pip install .

# 3. Verify
pypkatool --help
```

`pypkatool` locates the Python 2.7 interpreter automatically by looking for a
conda environment literally named `py27` under your home directory
(`~/miniconda3/envs/py27/bin/python2.7` or `~/anaconda3/envs/py27/bin/python2.7`).
If you named it something else, either rename it to `py27` or prepend its
`bin/` to `PATH` yourself before running `pypkatool run`.

The `pdbfixer` environment (if created) is located the same way, by name,
under `~/miniconda3/envs/pdbfixer/bin/python` or
`~/anaconda3/envs/pdbfixer/bin/python`.

## Every time you want to use it

Only `pypkatool` needs activating - `py27` and `pdbfixer` are "install once
and forget":

```bash
conda activate pypkatool
pypkatool run my_protein.pdb --ph 7.0
```
