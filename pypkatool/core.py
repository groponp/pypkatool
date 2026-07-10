#!/usr/bin/env python3
"""pypkatool core module: PyPKA pKa results to CHARMM-GUI protonation states.

This module drives `PyPKA <https://doi.org/10.1021/acs.jcim.0c00718>`_
(Poisson-Boltzmann + Monte Carlo pKa calculation) on an input structure,
cross-validates every predicted pKa against the independent machine-learning
predictor `pKAI+ <https://doi.org/10.1021/acs.jctc.2c00308>`_, maps each
titratable site's most probable tautomer onto a CHARMM36 residue/patch label
(``HSD``/``HSE``/``HSP``, ``ASPP``, ``GLUP``, ``LSN``, ``CNEU``, ``CYSD``,
``TYRD``, ``SERD``, ...), and writes a CHARMM-GUI-ready protonation table.

**Tautomer reference-state convention** (verified against the installed
PyPKA source, ``titsite.py::Titsite.getRefProtState()``, and against the
signed atomic partial charges in the CHARMM36 ``.st`` tautomer files):

* Cationic sites (``HIS``, ``LYS``, ``NTR``) - the *reference* tautomer
  (index ``N + 1``, where ``N`` is the number of regular tautomers) is the
  **protonated**, positively charged state.
* Anionic sites (``ASP``, ``GLU``, ``CTR``, ``CYS``, ``TYR``, ``SER``) - the
  reference tautomer is the **deprotonated**, negatively charged state; the
  regular tautomers (1..N) are the protonated, neutral forms (e.g. COOH for
  ASP/GLU, SH for CYS, OH for TYR).

This is the opposite assignment of what a naive reading of "reference =
default" might suggest for anionic sites, and it is easy to get backwards -
see :func:`_label` for the full mapping table.

Usage:
    ``pypkatool run <pdb> --ph <float> [--outdir <dir>] [--ncpus N] [--epsin F]``

    ``pypkatool reprocess <outdir> --ph <float> [--pdb <pdb>]``

:seealso: PyPKA DOI 10.1021/acs.jcim.0c00718; pKAI/pKAI+ DOI
    10.1021/acs.jctc.2c00308.
"""
from __future__ import annotations

import argparse
import importlib.metadata
import json
import os
import re
import sys
import dataclasses
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

__version__ = "1.0.0"

# ── pKAI+ prerequisite (hard requirement - must be installed) ─────────────────

def _require_pkai() -> None:
    """Abort import with an actionable message if the ``pkai`` package is missing.

    Called once at module import time. pKAI+ cross-validation is not
    optional in this pipeline: every run's output includes a
    ``crossvalidation_report.dat``, so the dependency is enforced up front
    rather than failing deep inside a multi-minute PB+MC run.

    :raises RuntimeError: if the ``pkai`` package cannot be found on
        ``sys.path``.
    """
    import importlib.util
    if importlib.util.find_spec("pkai") is None:
        raise RuntimeError(
            "pKAI is required but not installed.\n"
            "  pip install pKAI\n"
            "pKAI+ cross-validation is mandatory - every output includes it."
        )

_require_pkai()

# ── RTF ───────────────────────────────────────────────────────────────────────

#: Path to the bundled CHARMM36 protein topology file (``top_all36_prot.rtf``).
#: This is the single authority used throughout the module to decide which
#: CHARMM labels are valid RESI/PRES entries and what H atoms each implies.
RTF_PATH = Path(__file__).parent / "data" / "top_all36_prot.rtf"

@dataclass
class RtfBlock:
    """One verbatim ``RESI`` or ``PRES`` block parsed from the CHARMM36 RTF.

    :ivar kind: ``"RESI"`` for a full residue definition or ``"PRES"`` for a
        patch.
    :ivar name: Block name in upper case, e.g. ``"HSD"``, ``"ASPP"``.
    :ivar charge: Net formal charge declared on the block's header line.
    :ivar verbatim: The block's full original text, used both for display in
        ``detail.json`` and for atom-presence checks elsewhere in this module.
    """
    kind: str; name: str; charge: float; verbatim: str

def _load_rtf(path: Path = RTF_PATH) -> dict[str, RtfBlock]:
    """Parse a CHARMM ``.rtf`` topology file into ``RESI``/``PRES`` blocks.

    :param path: Path to the RTF file. Defaults to the bundled CHARMM36
        protein topology (:data:`RTF_PATH`).
    :type path: pathlib.Path
    :returns: Mapping of upper-cased block name to :class:`RtfBlock`.
    :rtype: dict[str, RtfBlock]
    :raises FileNotFoundError: if ``path`` does not exist.
    """
    if not path.exists():
        raise FileNotFoundError(f"RTF not found: {path}")
    text = path.read_text(errors="replace")
    blocks: dict[str, RtfBlock] = {}
    pat = re.compile(r"^(RESI|PRES)\s+(\S+)\s+([-\d.]+)", re.MULTILINE)
    ms = list(pat.finditer(text))
    for i, m in enumerate(ms):
        end = ms[i + 1].start() if i + 1 < len(ms) else len(text)
        name = m.group(2).upper()
        blocks[name] = RtfBlock(m.group(1), name, float(m.group(3)), text[m.start():end].rstrip())
    return blocks

def _rtf_has(name: str, rtf: dict) -> bool:
    """Return whether ``name`` has a RESI/PRES block in the parsed RTF.

    :param name: CHARMM label to look up (case-insensitive).
    :type name: str
    :param rtf: RTF blocks as returned by :func:`_load_rtf`.
    :type rtf: dict[str, RtfBlock]
    :rtype: bool
    """
    return name.upper() in rtf

def _rtf_get(name: str, rtf: dict) -> RtfBlock | None:
    """Look up a single RTF block by label.

    :param name: CHARMM label to look up (case-insensitive).
    :type name: str
    :param rtf: RTF blocks as returned by :func:`_load_rtf`.
    :type rtf: dict[str, RtfBlock]
    :returns: The matching block, or ``None`` if absent from the RTF.
    :rtype: RtfBlock | None
    """
    return rtf.get(name.upper())

# ── PyPKA runner ──────────────────────────────────────────────────────────────

DEFAULT_PARAMS: dict[str, Any] = {
    "epsin": 15, "ionicstr": 0.1, "pbc_dimensions": 0,
    "convergence": 0.1, "pH": "0,14",
    "ncpus": os.cpu_count() or 4,
    # PyPKA defaults to ffID="G54A7" (GROMOS) / ffinput="GROMOS" when these are
    # not set explicitly (see pypka/config.py Config.__init__). This tool's
    # entire premise - the bundled top_all36_prot.rtf, the CHARMM RESI/PRES
    # labels produced by _label(), and structure_output's "charmm" naming
    # scheme - assumes the underlying PB+MC calculation itself used CHARMM36m
    # charges/radii/tautomer .st files, not GROMOS ones. Must be set explicitly.
    "ffID": "CHARMM36m",
}

@dataclass
class SiteResult:
    """One titratable site as returned by PyPKA, before CHARMM mapping.

    :ivar resname: Residue name as reported by PyPKA (``"HIS"``, ``"ASP"``, ...).
    :ivar resid: Residue number (icode-corrected string if applicable).
    :ivar chain: Chain identifier, normalized to a single upper-case letter.
    :ivar pka: Numeric pKa, or ``None`` if outside ``[0, 14)`` (see :func:`_fmt_pka`).
    :ivar pka_str: Display string for ``pka`` (``"<0.0"``, ``">14.0"``, or ``"%.2f"``).
    :ivar site_type: ``"c"`` (cationic) or ``"a"`` (anionic), from PyPKA ``getType()``.
    :ivar most_prob_taut: 1-indexed most probable tautomer at the target pH;
        equals ``n_regular_tautomers + 1`` when the *reference* tautomer wins.
    :ivar n_regular_tautomers: Number of non-reference tautomers for this site.
    :ivar populations: Named tautomer population fractions at the target pH,
        for display only (see :func:`_build_pops`).
    :ivar pct_protonated: Percent protonated at the target pH (see :func:`_charge_split`).
    :ivar pct_deprotonated: Percent deprotonated at the target pH.
    """
    resname: str; resid: int | str; chain: str
    pka: float | None; pka_str: str; site_type: str
    most_prob_taut: int; n_regular_tautomers: int
    populations: dict[str, float] = field(default_factory=dict)
    pct_protonated: float = 0.0; pct_deprotonated: float = 0.0

    def __post_init__(self):
        self.chain = (str(self.chain).strip() or "A").upper()

#: Number of *regular* (non-reference) tautomers per residue type, as defined
#: by the CHARMM36m ``.st`` tautomer files bundled with PyPKA.
_N_REG = {"HIS":2,"ASP":4,"GLU":4,"LYS":3,"NTR":3,"CTR":4,"CYS":3,"TYR":2,"SER":3,"THR":3}

#: Site polarity per residue type: ``"c"`` (cationic, reference = protonated)
#: or ``"a"`` (anionic, reference = deprotonated). See the module docstring.
_STYPE = {"HIS":"c","LYS":"c","NTR":"c","ASP":"a","GLU":"a","CTR":"a","CYS":"a","TYR":"a","SER":"a","THR":"a"}

def _fmt_pka(v: Any) -> tuple[float | None, str]:
    """Normalize a raw PyPKA pKa value to a numeric value plus display string.

    PyPKA scans pH 0-14; a site whose curve never crosses 50% protonation in
    that window has no defined pKa and is reported as always-deprotonated
    (< 0) or always-protonated (> 14).

    :param v: Raw value from ``site.pK`` (may be ``None`` or non-numeric).
    :type v: Any
    :returns: ``(value, text)``; ``value`` is ``None`` when outside ``[0, 14)``.
    :rtype: tuple[float | None, str]
    """
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None, "N/A"
    if f <= 0.0: return None, "<0.0"
    if f >= 14.0: return None, ">14.0"
    return f, f"{f:.2f}"

def _build_pops(resname: str, site_type: str, probs: list[float]) -> dict[str, float]:
    """Attach human-readable tautomer names to a raw probability vector.

    Display-only helper (used by :func:`_taut_detail` and ``detail.json``);
    it does not affect the CHARMM label decision, which is driven solely by
    ``most_prob_taut`` in :func:`_label`.

    :param resname: Residue name (``"HIS"``, ``"ASP"``, ...).
    :type resname: str
    :param site_type: ``"c"`` or ``"a"`` (unused directly, kept for symmetry
        with :func:`_charge_split`).
    :type site_type: str
    :param probs: Tautomer probabilities, ``[regular_1..regular_N, reference]``.
    :type probs: list[float]
    :returns: Mapping of tautomer name to probability.
    :rtype: dict[str, float]
    """
    r = resname.upper()
    names_map = {
        "HIS": ["HSD","HSE","HSP"],
        "ASP": ["ASPH_OD1_syn","ASPH_OD2_syn","ASPH_OD1_anti","ASPH_OD2_anti","ASP"],
        "GLU": ["GLUH_OE1_syn","GLUH_OE2_syn","GLUH_OE1_anti","GLUH_OE2_anti","GLU"],
        "LYS": ["LSN_tau1","LSN_tau2","LSN_tau3","LYS"],
        "NTR": ["NNEU_tau1","NNEU_tau2","NNEU_tau3","NTER"],
        "CTR": ["CNEU_tau1","CNEU_tau2","CNEU_tau3","CNEU_tau4","CTER"],
    }
    names = names_map.get(r, [f"state_{i+1}" for i in range(len(probs)-1)] + ["ref"])
    while len(names) < len(probs): names.insert(-1, f"tau{len(names)}")
    return {n: float(p) for n, p in zip(names[:len(probs)], probs)}

def _charge_split(site_type: str, probs: list[float]) -> tuple[float, float]:
    """Split a tautomer probability vector into percent protonated/deprotonated.

    Uses the reference-tautomer convention documented at module level: for
    cationic sites the reference slot (last element of ``probs``) is the
    protonated population; for anionic sites it is the deprotonated one.

    :param site_type: ``"c"`` (cationic) or ``"a"`` (anionic).
    :type site_type: str
    :param probs: Tautomer probabilities, ``[regular_1..regular_N, reference]``.
    :type probs: list[float]
    :returns: ``(pct_protonated, pct_deprotonated)``, each in ``[0, 100]``.
    :rtype: tuple[float, float]
    """
    if not probs: return 50.0, 50.0
    if site_type == "c":
        return probs[-1]*100, sum(probs[:-1])*100
    return sum(probs[:-1])*100, probs[-1]*100

def _hh_probs(site_type: str, pka: float | None, ph: float, n_reg: int) -> list[float]:
    """Approximate tautomer populations from a scalar pKa via Henderson-Hasselbalch.

    Fallback used only when PyPKA's per-pH titration curve is unavailable
    (see the ``ph_key is None`` branch in :func:`run_pypka`); the exact MC
    populations from :func:`SiteResult.getTautsProb` are always preferred
    when present. The protonated fraction ``fp = 1 / (1 + 10**(pH - pKa))``
    is generic to any weak acid/base pair and independent of site polarity;
    only its placement into the reference vs. regular slots differs by
    ``site_type``.

    :param site_type: ``"c"`` (cationic) or ``"a"`` (anionic).
    :type site_type: str
    :param pka: Site pKa, or ``None`` if outside ``[0, 14)`` - in that case a
        uniform distribution is returned (protonation state is then decided
        by the caller from the boundary, not from this vector).
    :type pka: float | None
    :param ph: Target pH.
    :type ph: float
    :param n_reg: Number of regular (non-reference) tautomers.
    :type n_reg: int
    :returns: Probability vector ``[regular_1..regular_N, reference]``.
    :rtype: list[float]
    """
    total = n_reg + 1
    if pka is None: return [1.0/total]*total
    fp = 1.0 / (1.0 + 10**(ph - pka))
    probs = [0.0]*total
    if site_type == "c":
        probs[-1] = fp
        for i in range(n_reg): probs[i] = (1.0-fp)/n_reg if n_reg else 0
    else:
        probs[-1] = 1.0-fp
        for i in range(n_reg): probs[i] = fp/n_reg if n_reg else 0
    return probs

_STANDARD_AA = frozenset({
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
    "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
    "HSD","HSE","HSP","CYX","CYM",
})

def validate_pdb(pdb_path: Path) -> None:
    """Sanity-check a PDB file before handing it to PyPKA.

    Verifies the file exists, is readable, has well-formed ``ATOM``/``HETATM``
    coordinate columns, and contains at least one standard protein residue.
    Non-fatal issues (insertion codes, non-standard residues) are printed as
    notes rather than raising, since PyPKA can often still process them.

    :param pdb_path: Path to the input PDB file.
    :type pdb_path: pathlib.Path
    :raises SystemExit: on a missing/unreadable file, malformed ATOM/HETATM
        records, or a PDB with no recognizable standard protein residues.
    """
    if not pdb_path.exists():
        sys.exit(f"PDB not found: {pdb_path}")
    try:
        text = pdb_path.read_text(errors="replace")
    except OSError as e:
        sys.exit(f"Cannot read PDB: {e}")

    atom_lines = []
    for idx, ln in enumerate(text.splitlines(), 1):
        if ln.startswith(("ATOM  ", "ATOM ", "HETATM")):
            if len(ln) < 54:
                sys.exit(
                    f"ERROR: PDB format error at line {idx}: ATOM/HETATM record is too short\n"
                    f"  (got {len(ln)} characters, expected at least 54 for coordinates).\n"
                    f"  Line: {ln!r}"
                )
            try:
                float(ln[30:38])
                float(ln[38:46])
                float(ln[46:54])
            except ValueError:
                sys.exit(
                    f"ERROR: PDB format error at line {idx}: coordinates in columns 31-54 are not valid floats.\n"
                    f"  Line: {ln!r}"
                )
            atom_lines.append(ln)

    if not atom_lines:
        sys.exit(f"PDB has no ATOM/HETATM records: {pdb_path}")

    residues: set[str] = set()
    chains:   set[str] = set()
    non_std:  list[str] = []
    has_icode = False

    for ln in atom_lines:
        resname = ln[17:20].strip()
        chain   = (ln[21:22].strip() or "A").upper()
        icode   = ln[26:27].strip()
        residues.add(resname)
        chains.add(chain)
        if icode:
            has_icode = True
        if resname not in _STANDARD_AA and resname not in non_std:
            non_std.append(resname)

    std = residues & _STANDARD_AA
    if not std:
        sys.exit(
            f"PDB appears to contain no standard protein residues "
            f"(found: {sorted(residues)[:10]}). PyPKA requires a protein PDB."
        )
    if has_icode:
        print(f"  Note: PDB has insertion codes - PyPKA will renumber residues internally.")
    if non_std:
        print(f"  Note: non-standard residues present (will be ignored by PyPKA): {non_std[:5]}")
    print(f"  PDB OK - {len(atom_lines)} atoms, {len(chains)} chain(s): {sorted(chains)}")

def _inject_py27() -> None:
    """Put a Python 2.7 interpreter on ``PATH`` if one isn't already there.

    PyPKA's bundled ``pdbmender``/``pdb2pqr`` internally shells out to a
    ``python2.7`` binary. If none is visible on ``PATH``, this looks for the
    conda environment named ``py27`` (created via ``environment-py27.yml``,
    see the project README) under the current user's home directory, or
    ``$HOME``, and prepends its ``bin/`` to ``PATH`` for the current process.
    A no-op if ``python2.7`` is already resolvable.
    """
    import shutil
    if shutil.which("python2.7"): return
    for d in [Path.home() / "miniconda3/envs/py27/bin",
              Path.home() / "anaconda3/envs/py27/bin",
              Path(os.environ.get("HOME","/root")) / "miniconda3/envs/py27/bin"]:
        if (d / "python2.7").exists():
            os.environ["PATH"] = str(d) + ":" + os.environ.get("PATH","")
            return

# ── fixstructure (PDBFixer, separate env) ─────────────────────────────────────

#: Bundled worker script run inside the ``pdbfixer`` conda environment (see
#: :func:`_find_pdbfixer_python`). Kept as a standalone script rather than an
#: importable module because pdbfixer/OpenMM need numpy>=2, incompatible with
#: this package's own numpy<2 (delphi4py) environment.
_FIXSTRUCTURE_WORKER = Path(__file__).parent / "data" / "fixstructure_worker.py"

def _find_pdbfixer_python() -> str:
    """Locate the Python interpreter of a conda environment with pdbfixer installed.

    Looks for a conda environment named ``pdbfixer`` (see
    ``environment-pdbfixer.yml`` / the project README "Repairing fragmented
    structures") under common conda roots. pdbfixer/OpenMM require numpy>=2,
    which conflicts with delphi4py's numpy<2 pin, so this always runs as a
    separate interpreter rather than an in-process import.

    :returns: Path to the ``python`` executable inside the ``pdbfixer`` env.
    :rtype: str
    :raises SystemExit: if no such environment can be found.
    """
    for d in [Path.home() / "miniconda3/envs/pdbfixer/bin",
              Path.home() / "anaconda3/envs/pdbfixer/bin",
              Path(os.environ.get("HOME","/root")) / "miniconda3/envs/pdbfixer/bin"]:
        if (d / "python").exists():
            return str(d / "python")
    sys.exit(
        "fixstructure requires a separate conda environment with pdbfixer + openmm.\n"
        "  conda create -n pdbfixer -c conda-forge python=3.11 pdbfixer openmm -y\n"
        "See README.md 'Repairing fragmented structures'."
    )

def fix_structure(pdb_path: Path, outdir: Path | None = None) -> Path:
    """Repair a structurally incomplete PDB using PDBFixer, before feeding it to PyPKA.

    Delegates to :data:`_FIXSTRUCTURE_WORKER`, run under
    :func:`_find_pdbfixer_python`, which:

    * fills in missing heavy atoms for every residue that is present
      (including at chain termini, e.g. a missing ``OXT``);
    * rebuilds missing residues only when the gap is strictly internal to a
      chain (a present residue on both sides) - a gap at the very start or
      end of a chain is left untouched, since that reflects where the
      deposited construct actually starts/ends, not damage to repair;
    * never adds hydrogens, so the output stays heavy-atom-only like a
      standard deposited PDB.

    :param pdb_path: Input PDB to repair.
    :type pdb_path: pathlib.Path
    :param outdir: Directory for the output file. Defaults to ``pdb_path``'s
        parent directory.
    :type outdir: pathlib.Path | None
    :returns: Path to the repaired PDB (``<stem>_fixed.pdb``).
    :rtype: pathlib.Path
    :raises SystemExit: if the ``pdbfixer`` environment is missing, or the
        worker subprocess fails (e.g. on a PDB with no recognizable sequence).
    """
    import subprocess
    outdir = outdir or pdb_path.parent
    outdir.mkdir(parents=True, exist_ok=True)
    out_pdb = outdir / f"{pdb_path.stem}_fixed.pdb"
    python = _find_pdbfixer_python()
    result = subprocess.run(
        [python, str(_FIXSTRUCTURE_WORKER), "--pdb", str(pdb_path), "--output", str(out_pdb)],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        sys.exit(f"fixstructure failed:\n{result.stderr}")
    print(f"  {result.stdout.strip()}")
    return out_pdb


def run_pypka(pdb_path: Path, target_ph: float, outdir: Path,
              ncpus: int = 4, epsin: float = 15,
              charmm_input: bool = False) -> tuple[list[SiteResult], Path]:
    """Run PyPKA (Poisson-Boltzmann + Monte Carlo) and collect per-site results.

    Writes PyPKA's raw output (protonated PDB, titration curve) under
    ``outdir/output_pypka/`` and a ``protocol.json`` capturing the exact
    parameters used, then reads back every titratable site's most probable
    tautomer at ``target_ph`` via the PyPKA ``Titration`` API.

    :param pdb_path: Input protein structure (PDB format).
    :type pdb_path: pathlib.Path
    :param target_ph: pH at which to evaluate the most probable protonation state.
    :type target_ph: float
    :param outdir: Output directory (created if missing).
    :type outdir: pathlib.Path
    :param ncpus: Number of CPUs for PyPKA's parallel PB solve. Defaults to 4.
    :type ncpus: int
    :param epsin: Protein interior dielectric constant. Defaults to 15
        (PyPKA's literature-optimized value; RMSE 0.82 / MAE 0.57 on the
        PKAD benchmark).
    :type epsin: float
    :param charmm_input: Set ``ffinput="CHARMM"`` (PyPKA/pdbmender default is
        ``"GROMOS"``). Only matters if ``pdb_path`` already carries CHARMM
        protonation-state-specific residue names (``HSD``/``HSE``/``HSP``/
        ``ASPP``/``GLUP``/``CYSM``/...) - e.g. a PDB re-exported from a
        CHARMM-GUI PDB Reader step. Verified against every call site that
        reads ``ffinput`` (``pdbmender/utils.py::identify_tit_sites``,
        ``pypka/clean/cleaning.py`` x2): all three only fire on those exact
        residue names, via ``pdbmender.ffconverter.CHARMM_protomers``. A
        standard PDB (RCSB, AlphaFold, ...) never contains them, so this flag
        is a no-op for the common case - it exists for the CHARMM-GUI
        round-trip case only. Defaults to ``False``.
    :type charmm_input: bool
    :returns: ``(sites, protonated_pdb_path)``.
    :rtype: tuple[list[SiteResult], pathlib.Path]
    :raises SystemExit: if the ``pypka`` package is not installed, or if
        ``Titration()`` fails (most commonly a structurally incomplete PDB -
        see the ``except Exception`` branch below).
    """
    try:
        from pypka import Titration
    except ImportError:
        sys.exit("pypka not installed. Run: pip install pypka")

    outdir.mkdir(parents=True, exist_ok=True)
    raw_dir = outdir / "output_pypka"   # all PyPKA raw files go here
    raw_dir.mkdir(exist_ok=True)
    prot_pdb = raw_dir / f"{pdb_path.stem}_protonated_pH{target_ph}.pdb"
    titration_out = raw_dir / f"{pdb_path.stem}_titration.dat"

    params: dict[str, Any] = {**DEFAULT_PARAMS, "epsin": epsin, "ncpus": ncpus,
        "structure": str(pdb_path.resolve()),
        "structure_output": f"{prot_pdb}, {target_ph}, charmm",
        "titration_output": str(titration_out),
    }
    if charmm_input:
        params["ffinput"] = "CHARMM"
    _save_protocol(params, pdb_path, target_ph, outdir)
    _inject_py27()
    orig = os.getcwd(); os.chdir(raw_dir)   # PyPKA temp files land in output_pypka/
    try:
        try:
            tit = Titration(params, sites="all")
        except Exception as e:
            # PyPKA raises bare Exceptions from deep inside pdbmender's cleanPDB()
            # for structurally incomplete input (e.g. "CTR N has missing atoms" -
            # a truncated C-terminus with no OXT/O2) with no indication of the
            # fix. validate_pdb() cannot catch this upfront: it only checks that
            # ATOM/HETATM records are well-formed, not that each residue has its
            # expected atom set.
            sys.exit(
                f"PyPKA failed while preprocessing the structure:\n  {e}\n\n"
                f"This usually means the PDB is structurally incomplete for at least\n"
                f"one residue (most often a truncated terminus missing its O/OXT/OT2\n"
                f"atoms, or a residue missing backbone/side-chain atoms). PyPKA does\n"
                f"not repair missing heavy atoms - re-run on a structure that has been\n"
                f"completed first (e.g. with PDBFixer, MODELLER, or your PDB's\n"
                f"'*_clean.pdb' / repaired variant if one exists)."
            )
    finally:
        os.chdir(orig)

    results: list[SiteResult] = []
    for site in tit:
        try:
            resname = site.getName(); resid = site.getResNumber(correct_icode=True)
            chain = (site.molecule.chain.strip() or "A").upper()
            pka, pka_str = _fmt_pka(site.pK)
            site_type = site.getType(); n_reg = site.getNTautomers()
            try:
                curve = site.getTitrationCurve()
                ph_key = min(curve, key=lambda k: abs(float(k)-target_ph)) if curve else None
            except Exception: ph_key = None
            if ph_key is not None:
                mpt, _ = site.getMostProbTaut(ph_key)
                probs = site.getTautsProb(ph_key)
            else:
                mpt = (n_reg+1) if ((site_type=="c") == (pka is None or pka >= target_ph)) else 1
                probs = _hh_probs(site_type, pka, target_ph, n_reg)
            pct_p, pct_d = _charge_split(site_type, probs)
            results.append(SiteResult(resname=resname, resid=resid, chain=chain,
                pka=pka, pka_str=pka_str, site_type=site_type,
                most_prob_taut=mpt, n_regular_tautomers=n_reg,
                populations=_build_pops(resname, site_type, probs),
                pct_protonated=pct_p, pct_deprotonated=pct_d))
        except Exception as e:
            print(f"  WARNING: site {site}: {e}")
    return results, prot_pdb

# ── reprocess (from existing titration.dat + protonated PDB) ─────────────────

def _parse_titration_dat(p: Path):
    """Parse a PyPKA ``*_titration.dat`` file into per-site titration curves.

    :param p: Path to the titration data file (``#pH <site>_<chain>_<resname> ...``
        header followed by one row per pH step).
    :type p: pathlib.Path
    :returns: ``(sites, curves)`` where ``sites`` is a list of
        ``(resid, chain, resname)`` tuples in column order and ``curves`` maps
        each original column key to a list of ``(pH, fraction_protonated)``.
    :rtype: tuple[list[tuple], dict[str, list[tuple[float, float]]]]
    :raises ValueError: if the file has no ``#pH`` header line.
    """
    lines = p.read_text().splitlines()
    col_keys, rows = [], []
    for ln in lines:
        s = ln.strip()
        if s.startswith("#pH"):
            col_keys = s.lstrip("#").split()[1:]
        elif s and not s.startswith("#"):
            rows.append(s.split())
    if not col_keys: raise ValueError(f"No #pH header in {p}")
    sites = []
    for k in col_keys:
        parts = k.split("_")
        try: resid = int(parts[0])
        except ValueError: resid = parts[0]
        sites.append((resid, parts[1], "_".join(parts[2:]) if len(parts)>3 else parts[2]))
    curves = {k: [] for k in col_keys}
    for row in rows:
        if len(row) < len(col_keys)+1: continue
        ph = float(row[0])
        for i, k in enumerate(col_keys):
            curves[k].append((ph, float(row[i+1])))
    return sites, curves

def _pka_from_curve(curve) -> tuple[float|None, str]:
    """Recover a pKa by linear interpolation of a titration curve across 50%.

    :param curve: List of ``(pH, fraction_protonated)`` samples, pH-ordered.
    :type curve: list[tuple[float, float]]
    :returns: ``(pka, display_string)`` at the first crossing of 0.5
        protonation found; ``(None, "N/A")`` if the curve never crosses 0.5.
    :rtype: tuple[float | None, str]
    """
    for i in range(1, len(curve)):
        ph0, p0 = curve[i-1]; ph1, p1 = curve[i]
        if (p0-0.5)*(p1-0.5) <= 0:
            dp = p1-p0
            pka = ph0 + (ph1-ph0)*(0.5-p0)/dp if abs(dp)>1e-12 else (ph0+ph1)/2
            return _fmt_pka(pka)
    return None, "N/A"

def _read_pdb_names(pdb: Path) -> dict[tuple[str,int|str], str]:
    """Read the first residue name seen per ``(chain, resid)`` from a PDB.

    :param pdb: Path to a PDB file.
    :type pdb: pathlib.Path
    :returns: Mapping of ``(chain, resid)`` to residue name (upper case).
    :rtype: dict[tuple[str, int | str], str]
    """
    names: dict[tuple[str,int|str], str] = {}
    for line in pdb.read_text().splitlines():
        if len(line) < 26: continue
        if not line.startswith(("ATOM","HETATM")): continue
        rn = line[17:21].strip().upper(); ch = (line[21].strip() or "A").upper()
        try: rid = int(line[22:26])
        except ValueError: rid = line[22:27].strip()
        if (ch, rid) not in names: names[(ch,rid)] = rn
    return names

def reprocess_from_files(titration_dat: Path, protonated_pdb: Path,
                         target_ph: float) -> list[SiteResult]:
    sites, curves = _parse_titration_dat(titration_dat)
    pdb_names = _read_pdb_names(protonated_pdb)
    col_keys = list(curves)
    results: list[SiteResult] = []
    for i, (resid, chain, resname) in enumerate(sites):
        curve = curves[col_keys[i]]
        pka, pka_str = _pka_from_curve(curve)
        prot_frac = min(curve, key=lambda x: abs(x[0]-target_ph))[1] if curve else 0.5
        pct_p, pct_d = prot_frac*100, (1-prot_frac)*100
        st = _STYPE.get(resname.upper(), "a"); nr = _N_REG.get(resname.upper(), 1)
        chain_up = (str(chain).strip() or "A").upper()
        if resname.upper() == "HIS":
            cn = pdb_names.get((chain_up, resid), "HSP")
            mpt = {"HSD":1,"HSE":2}.get(cn, 3)
            nf = 1.0-prot_frac
            probs = [0.0, 0.0, prot_frac]
            if mpt < 3: probs[mpt-1] = nf
            else: probs[0] = probs[1] = nf/2
        else:
            mpt = (nr+1 if pct_p>50 else 1) if st=="c" else (1 if pct_p>50 else nr+1)
            probs = _hh_probs(st, pka, target_ph, nr)
            # Use exact fraction instead of HH approximation
            total = nr+1; probs = [0.0]*total
            if st=="c":
                probs[-1]=prot_frac
                for j in range(nr): probs[j]=(1-prot_frac)/nr if nr else 0
            else:
                probs[-1]=1-prot_frac
                for j in range(nr): probs[j]=prot_frac/nr if nr else 0
        results.append(SiteResult(resname=resname, resid=resid, chain=chain_up,
            pka=pka, pka_str=pka_str, site_type=st, most_prob_taut=mpt,
            n_regular_tautomers=nr, populations=_build_pops(resname, st, probs),
            pct_protonated=pct_p, pct_deprotonated=pct_d))
    return results

# ── State mapper ──────────────────────────────────────────────────────────────

NO_ACTION = frozenset({"ASP","GLU","LYS","ARG","TYR","SER","THR","CYS","NTER","CTER"})
SKIP_TABLE = frozenset({"NTR","CTR"})  # terminal groups: included in JSON, excluded from dat

@dataclass
class MappedResidue:
    """A titratable site after mapping to a CHARMM36 label.

    :ivar resname: Original PyPKA residue name.
    :ivar resid: Residue number.
    :ivar chain: Chain identifier.
    :ivar pka_str: Display pKa string, propagated from :class:`SiteResult`.
    :ivar pct_protonated: Percent protonated at the target pH.
    :ivar pct_deprotonated: Percent deprotonated at the target pH.
    :ivar charmm_label: Resolved CHARMM36 RESI/PRES label (e.g. ``"HSE"``, ``"ASPP"``).
    :ivar rtf_available: Whether ``charmm_label`` has a block in the bundled RTF.
    :ivar is_tie: True when protonated/deprotonated percentages are within
        :data:`TIE_MARGIN` of each other (ambiguous state).
    :ivar tautomer_detail: Human-readable tautomer breakdown for display.
    :ivar needs_action: Whether CHARMM-GUI needs a RESI change or PRES patch.
    :ivar site_type: PyPKA site polarity (``"c"``/``"a"``), kept for cross-validation.
    :ivar mpt: Most probable tautomer index, kept for cross-validation.
    :ivar n_reg: Number of regular tautomers, kept for cross-validation.
    """
    resname: str; resid: int|str; chain: str; pka_str: str
    pct_protonated: float; pct_deprotonated: float
    charmm_label: str; rtf_available: bool; is_tie: bool
    tautomer_detail: str; needs_action: bool
    # PyPKA source fields - kept for label-chain cross-validation
    site_type: str = ""; mpt: int = 0; n_reg: int = 0

    def __post_init__(self):
        self.chain = (str(self.chain).strip() or "A").upper()

#: Percentage-point margin within which protonated/deprotonated populations
#: are considered a tie (ambiguous state, flagged ``ACTION=TIE`` in the report).
TIE_MARGIN = 2.0

def _label(resname: str, mpt: int, n_reg: int, site_type: str) -> str:
    """Map a PyPKA tautomer index to its CHARMM36 RESI/PRES label.

    ``ref = (mpt == n_reg + 1)`` identifies the *reference* tautomer. Per the
    module-level convention (verified against ``getRefProtState()`` in the
    PyPKA source and the CHARMM36 ``.st`` tautomer charge sets):

    * Cationic (``site_type == "c"``): reference = protonated. For ``LYS``
      that is the CHARMM *default* RESI - no patch needed. For ``NTR`` the
      reference maps to ``NTER``, which - unlike ``LYS`` - is a ``PRES`` in
      the RTF, not a ``RESI`` (a chain's N-terminus is always applied as a
      patch on its first residue); it is, however, the default terminal
      patch CHARMM-GUI's PDB Reader applies automatically, so no explicit
      action is needed there either. Non-reference = neutral = needs a
      patch (``LSN``, ``NNEU``).
    * Anionic (``site_type == "a"``): reference = deprotonated. For
      ``ASP``/``GLU`` the deprotonated, charged form *is* the CHARMM
      default RESI - no patch needed; the protonated form needs a patch
      (``ASPP``, ``GLUP``). For ``CTR`` the reference maps to ``CTER``,
      which - like ``NTER`` - is a ``PRES``, not a ``RESI``, but is the
      default C-terminal patch applied automatically. For
      ``CYS``/``TYR``/``SER`` it is the opposite: the CHARMM default RESI is
      the *protonated* neutral form (``CYS``, ``TYR``, ``SER``), so the
      deprotonated reference state is the one that needs a patch
      (``CYSD``, ``TYRD``, ``SERD``).

    ``HIS`` is special-cased on the raw tautomer index rather than the
    reference flag: ``mpt == 1`` -> ``HSD`` (ND1-H), ``mpt == 2`` -> ``HSE``
    (NE2-H), anything else (the reference, doubly-protonated state) -> ``HSP``.

    ``ARG`` is never titrated by PyPKA (pKa ~12.5, always protonated at
    biological pH) and never reaches this function. ``THR`` is titratable
    (``n_reg == 3``) but CHARMM36 has no ``PRES THRD``, so it is always
    returned as plain ``"THR"`` regardless of ``ref``; callers surface a
    warning instead (see :func:`reconcile_non_his_from_pdb`).

    :param resname: Residue name (``"HIS"``, ``"ASP"``, ...).
    :type resname: str
    :param mpt: Most probable tautomer index (1-indexed).
    :type mpt: int
    :param n_reg: Number of regular (non-reference) tautomers for this residue.
    :type n_reg: int
    :param site_type: ``"c"`` (cationic) or ``"a"`` (anionic).
    :type site_type: str
    :returns: CHARMM36 label, e.g. ``"HSE"``, ``"ASPP"``, ``"LSN"``.
    :rtype: str
    """
    r = resname.upper()
    ref = (mpt == n_reg+1)
    if r == "HIS": return {1:"HSD",2:"HSE"}.get(mpt,"HSP")
    if site_type == "c":
        return {"LYS": "LYS" if ref else "LSN",
                "NTR": "NTER" if ref else "NNEU"}.get(r, r)
    return {"ASP": "ASP"  if ref else "ASPP",
            "GLU": "GLU"  if ref else "GLUP",
            "CTR": "CTER" if ref else "CNEU",
            "CYS": "CYSD" if ref else "CYS",
            "TYR": "TYRD" if ref else "TYR",
            "SER": "SERD" if ref else "SER"}.get(r, r)

def _taut_detail(resname: str, pops: dict, pct_p: float, pct_d: float) -> str:
    """Format a short human-readable tautomer breakdown for display.

    Display-only; does not influence :func:`_label`. For HIS, always shows
    the HSP fraction plus the HSD/HSE split within the neutral population.
    For ASP/GLU with a non-negligible protonated fraction, shows which
    carboxylate oxygen (OD1/OD2 or OE1/OE2) carries the proton.

    :param resname: Residue name.
    :type resname: str
    :param pops: Named tautomer populations, as built by :func:`_build_pops`.
    :type pops: dict[str, float]
    :param pct_p: Percent protonated.
    :type pct_p: float
    :param pct_d: Percent deprotonated.
    :type pct_d: float
    :returns: A short descriptive string, or ``"-"`` if not applicable.
    :rtype: str
    """
    r = resname.upper()
    if r == "HIS":
        hsd=pops.get("HSD",0); hse=pops.get("HSE",0); hsp=pops.get("HSP",0)
        nt = hsd+hse
        s = f"HSP:{hsp*100:.0f}%"
        if nt>0: s += f" | neutral→ HSD:{100*hsd/nt:.0f}% / HSE:{100*hse/nt:.0f}%"
        return s
    if r == "ASP" and pct_p >= 5:
        od1=pops.get("ASPH_OD1_syn",0)+pops.get("ASPH_OD1_anti",0)
        od2=pops.get("ASPH_OD2_syn",0)+pops.get("ASPH_OD2_anti",0)
        t=od1+od2
        if t>0: return f"OD2-H(ASPP):{100*od2/t:.0f}% / OD1-H:{100*od1/t:.0f}%"
    if r == "GLU" and pct_p >= 5:
        oe1=pops.get("GLUH_OE1_syn",0)+pops.get("GLUH_OE1_anti",0)
        oe2=pops.get("GLUH_OE2_syn",0)+pops.get("GLUH_OE2_anti",0)
        t=oe1+oe2
        if t>0: return f"OE2-H(GLUP):{100*oe2/t:.0f}% / OE1-H:{100*oe1/t:.0f}%"
    return "-"

def map_residue(site: SiteResult, rtf: dict) -> MappedResidue:
    """Map one PyPKA :class:`SiteResult` to a :class:`MappedResidue`.

    Applies :func:`_label` for the CHARMM label, checks RTF availability,
    flags ties, and (for HIS) always requires CHARMM-GUI action since a RESI
    choice (HSD/HSE/HSP) must be made explicitly regardless of pKa.

    :param site: Raw PyPKA site result.
    :type site: SiteResult
    :param rtf: Parsed CHARMM36 RTF blocks, as returned by :func:`_load_rtf`.
    :type rtf: dict[str, RtfBlock]
    :returns: The mapped residue, prior to PDB-based reconciliation.
    :rtype: MappedResidue
    """
    lbl = _label(site.resname, site.most_prob_taut, site.n_regular_tautomers, site.site_type)
    tie = abs(site.pct_protonated - site.pct_deprotonated) <= TIE_MARGIN
    action = (site.resname.upper() == "HIS") or (lbl not in NO_ACTION)
    return MappedResidue(
        resname=site.resname, resid=site.resid, chain=site.chain,
        pka_str=site.pka_str, pct_protonated=site.pct_protonated,
        pct_deprotonated=site.pct_deprotonated,
        charmm_label=lbl, rtf_available=_rtf_has(lbl, rtf),
        is_tie=tie,
        tautomer_detail=_taut_detail(site.resname, site.populations,
                                     site.pct_protonated, site.pct_deprotonated),
        needs_action=action,
        site_type=site.site_type, mpt=site.most_prob_taut, n_reg=site.n_regular_tautomers,
    )

# RTF-primary protonation profile for HIS RESI blocks.
# RTF is the AUTHORITY: defines which imidazole H atoms MUST and MUST NOT appear
# in the protonated PDB for each HIS label to be consistent with CHARMM36 topology.
#
# Parsed from top_all36_prot.rtf RESI blocks:
#   RESI HSD: ATOM HD1 H  0.32  (only HD1 on ND1 → NE2 is neutral)
#   RESI HSE: ATOM HE2 H  0.32  (only HE2 on NE2 → ND1 is neutral)
#   RESI HSP: ATOM HD1 H  0.44 + ATOM HE2 H  0.44  (both nitrogens protonated)
_RESI_H_PROFILE: dict[str, dict[str, frozenset[str]]] = {
    "HSD": {"required": frozenset({"HD1"}),        "forbidden": frozenset({"HE2"})},
    "HSE": {"required": frozenset({"HE2"}),        "forbidden": frozenset({"HD1"})},
    "HSP": {"required": frozenset({"HD1", "HE2"}), "forbidden": frozenset()},
}
_HIS_RTF_EXPECTED = {lbl: p["required"] for lbl, p in _RESI_H_PROFILE.items()}

def _his_rtf_h_atoms(rtf: dict) -> dict[str, frozenset[str]]:
    """Parse imidazole H atoms (HD1 / HE2) from HSD / HSE / HSP RESI blocks in the RTF.

    Standard CHARMM36 result:
        HSD → frozenset({"HD1"})   proton on ND1
        HSE → frozenset({"HE2"})   proton on NE2
        HSP → frozenset({"HD1","HE2"})  both nitrogens protonated

    :param rtf: Parsed CHARMM36 RTF blocks, as returned by :func:`_load_rtf`.
    :type rtf: dict[str, RtfBlock]
    :returns: Mapping of ``"HSD"``/``"HSE"``/``"HSP"`` to their required
        imidazole hydrogen names per the RTF.
    :rtype: dict[str, frozenset[str]]
    """
    out: dict[str, frozenset[str]] = {}
    for lbl in ("HSD", "HSE", "HSP"):
        blk = _rtf_get(lbl, rtf)
        if blk is None:
            out[lbl] = frozenset(); continue
        h: set[str] = set()
        for line in blk.verbatim.splitlines():
            parts = line.split()
            if parts and parts[0].upper() == "ATOM" and len(parts) > 1:
                name = parts[1].upper()
                if name in ("HD1", "HE2"):
                    h.add(name)
        out[lbl] = frozenset(h)
    return out


def _his_pdb_signals(prot_pdb: Path) -> dict[tuple[str, int], tuple[str, str]]:
    """Read protonated PDB → {(chain, resid): (label_from_name, label_from_atoms)}.

    label_from_name : residue name field (HSD / HSE / HSP / HIS / unknown)
    label_from_atoms: inferred from which imidazole N has a hydrogen:
                        HD1 present + HE2 absent  → HSD
                        HE2 present + HD1 absent  → HSE
                        HD1 + HE2 both present    → HSP
                        neither present           → UNKNOWN
    Both signals must agree; if they disagree the atom-based label is used
    (atoms are ground truth - the name field could be a PyPKA formatting bug).

    :param prot_pdb: PyPKA's output protonated PDB.
    :type prot_pdb: pathlib.Path
    :returns: Mapping of ``(chain, resid)`` to ``(label_from_name, label_from_atoms)``.
    :rtype: dict[tuple[str, int], tuple[str, str]]
    """
    name_first: dict[tuple[str,int], str] = {}
    atom_set:   dict[tuple[str,int], set[str]] = {}

    for line in prot_pdb.read_text().splitlines():
        if len(line) < 26: continue
        if not line.startswith(("ATOM", "HETATM")):
            continue
        resname = line[17:21].strip().upper()
        if resname not in ("HIS", "HSD", "HSE", "HSP"):
            continue
        chain = (line[21].strip() or "A").upper()
        try:
            resid = int(line[22:26])
        except ValueError:
            continue
        key = (chain, resid)
        if key not in name_first:
            name_first[key] = resname
        atom_set.setdefault(key, set()).add(line[12:16].strip().upper())

    result: dict[tuple[str,int], tuple[str,str]] = {}
    for key, atoms in atom_set.items():
        hd1 = "HD1" in atoms
        he2 = "HE2" in atoms
        if hd1 and he2:
            from_atoms = "HSP"
        elif hd1:
            from_atoms = "HSD"
        elif he2:
            from_atoms = "HSE"
        else:
            from_atoms = "UNKNOWN"
        result[key] = (name_first.get(key, "UNKNOWN"), from_atoms)

    return result


def reconcile_his_from_pdb(mapped: list[MappedResidue], prot_pdb: Path,
                            rtf: dict) -> list[MappedResidue]:
    """Override HIS labels using three independent signals.

    Signal 1 - residue name field in protonated PDB (set by PyPKA's formatter).

    Signal 2 - atom inventory: which imidazole N carries H::

        HD1 present, HE2 absent  -> HSD
        HE2 present, HD1 absent  -> HSE
        HD1 + HE2 both present   -> HSP
        neither                  -> UNKNOWN

    Signal 3 - CHARMM36 RTF ATOM definitions for each HIS label: HSD defines
    HD1 only; HSE defines HE2 only; HSP defines both.

    Decision priority: atoms (2) > name (1) > API-derived label.
    RTF (3) validates the final label against the topology; mismatch → ERROR.

    :param mapped: Residues already mapped by :func:`map_residue`.
    :type mapped: list[MappedResidue]
    :param prot_pdb: PyPKA's output protonated PDB (ground truth for atoms).
    :type prot_pdb: pathlib.Path
    :param rtf: Parsed CHARMM36 RTF blocks.
    :type rtf: dict[str, RtfBlock]
    :returns: ``mapped`` with every HIS entry's ``charmm_label`` reconciled
        against the protonated PDB; non-HIS entries are passed through unchanged.
    :rtype: list[MappedResidue]
    """
    signals = _his_pdb_signals(prot_pdb)
    rtf_h   = _his_rtf_h_atoms(rtf)

    # One-time RTF sanity: verify each HIS label block matches the required H atoms.
    # RTF is the AUTHORITY - if this fires, the top_all36_prot.rtf may be non-standard.
    for lbl, profile in _RESI_H_PROFILE.items():
        got = rtf_h.get(lbl, frozenset())
        exp = profile["required"]
        if got != exp:
            print(f"  RTF WARNING: {lbl} imidazole H block defines H={set(got)} "
                  f"(expected required={set(exp)}) - RTF file may be non-standard")

    # Lookup table: atom_lbl → which H atoms were seen in the PDB
    _pdb_h: dict[str, frozenset[str]] = {
        "HSD": frozenset({"HD1"}),
        "HSE": frozenset({"HE2"}),
        "HSP": frozenset({"HD1", "HE2"}),
    }

    out: list[MappedResidue] = []
    for r in mapped:
        if r.resname.upper() != "HIS":
            out.append(r); continue

        key = (r.chain.upper(), r.resid)
        info = signals.get(key)

        if info is None:
            print(f"  WARNING: HIS{r.resid} {r.chain}: not found in protonated PDB "
                  f"- keeping API label {r.charmm_label}")
            out.append(r); continue

        name_lbl, atom_lbl = info

        # ── Resolve final label (signals 1 + 2) ───────────────────────────────
        if atom_lbl == "UNKNOWN":
            final = name_lbl if name_lbl in ("HSD","HSE","HSP") else r.charmm_label
            rtf_h_for_final = set(rtf_h.get(final, frozenset()))
            print(f"  WARNING: HIS{r.resid} {r.chain}: no imidazole H in PDB atoms "
                  f"- using residue name ({final}); "
                  f"RTF[{final}] expects H={rtf_h_for_final}")
        elif name_lbl in ("HSD","HSE","HSP") and name_lbl != atom_lbl:
            final = atom_lbl  # atoms win
            print(f"  ERROR: HIS{r.resid} {r.chain}: name={name_lbl} but "
                  f"atoms(HD1/HE2)→{atom_lbl}; RTF[{atom_lbl}] expects "
                  f"H={set(rtf_h.get(atom_lbl, frozenset()))} "
                  f"- USING ATOM-BASED LABEL (PyPKA naming inconsistency)")
        else:
            final = atom_lbl

        # ── Signal 3: RTF-primary required/forbidden check ────────────────────
        # RTF is the AUTHORITY. For the resolved label, verify:
        #   (a) required H atoms are present in the PDB imidazole ring
        #   (b) forbidden H atoms are absent from the PDB imidazole ring
        # Normal path (atom_lbl → final): these checks trivially pass, confirming
        # internal consistency. They catch logic bugs if the label resolution above
        # produces an inconsistent final label.
        if atom_lbl not in (None, "UNKNOWN") and final in _RESI_H_PROFILE:
            profile    = _RESI_H_PROFILE[final]
            imid_h     = _pdb_h.get(atom_lbl, frozenset())
            missing    = profile["required"] - imid_h
            forbidden  = profile["forbidden"] & imid_h
            if missing:
                print(f"  RTF ERROR: HIS{r.resid} {r.chain}: assigned label={final} "
                      f"but RTF-required H atoms {set(missing)} absent from PDB imidazole")
            if forbidden:
                print(f"  RTF ERROR: HIS{r.resid} {r.chain}: assigned label={final} "
                      f"but RTF-forbidden H atoms {set(forbidden)} present in PDB imidazole")
            # Verify the RTF file itself also defines the expected atoms for this label
            rtf_h_for_final = rtf_h.get(final, frozenset())
            if rtf_h_for_final != profile["required"]:
                print(f"  RTF WARNING: HIS{r.resid} {r.chain}: RTF[{final}] defines "
                      f"H={set(rtf_h_for_final)} but profile expects {set(profile['required'])}"
                      f" - RTF/PDB topology inconsistency")

        if final != r.charmm_label:
            print(f"  NOTE: HIS{r.resid} {r.chain}: API→{r.charmm_label} PDB→{final}")

        out.append(MappedResidue(
            resname=r.resname, resid=r.resid, chain=r.chain, pka_str=r.pka_str,
            pct_protonated=r.pct_protonated, pct_deprotonated=r.pct_deprotonated,
            charmm_label=final, rtf_available=_rtf_has(final, rtf),
            is_tie=r.is_tie, tautomer_detail=r.tautomer_detail, needs_action=True))

    return out

# ── Non-HIS residue validation (RTF patches + PDB atom inventory) ─────────────

# H atom that signals protonation or deprotonation for each CHARMM patch.
# Tuple: (action, atom)  where action="ADD" means the patch introduces that atom,
# "DELETE" means it removes it.  Used to verify RTF and, where PyPKA reflects the
# state in the protonated PDB, to cross-check via atom inventory.
#
# H atom that signals each non-HIS protonation state (parsed from CHARMM36 RTF PRES blocks).
# RTF is the AUTHORITY: these define what atom is added (ADD) or deleted (DELETE)
# when the patch is applied in CHARMM-GUI relative to the default protonated residue.
#
# PyPKA protonated PDB observability - ALL types have PDB atom signals:
#   HIS  : residue renamed HSD/HSE/HSP + HD1/HE2 presence → 3-signal check ✓
#   ASP  : HD2 present if ASPP (confirmed: PyPKA writes HD2 via ffconverter/pdbmender) ✓
#   GLU  : HE2 present if GLUP (confirmed: PyPKA writes HE2 via ffconverter/pdbmender) ✓
#   LYS  : HZ count (3→LYS default; 2→LSN neutral) ✓
#   CYS  : HG1 present→CYS; absent→CYSD (or CYX if SS-bonded) ✓
#   TYR  : HH present→TYR; absent→TYRD ✓
#   SER  : HG1 present→SER; absent→SERD ✓
_PATCH_H_INDICATORS: dict[str, tuple[str, str]] = {
    "ASPP": ("ADD",    "HD2"),   # PRES ASPP (line 1578): ADD ATOM HD2 to OD2
    "GLUP": ("ADD",    "HE2"),   # PRES GLUP (line 1593): ADD ATOM HE2 to OE2
    "LSN":  ("DELETE", "HZ3"),   # PRES LSN  (line 1608): DELETE ATOM HZ3 from NZ
    "CYSD": ("DELETE", "HG1"),   # PRES CYSD (line 1659): DELETE ATOM HG1 from SG
    "SERD": ("DELETE", "HG1"),   # PRES SERD:              DELETE ATOM HG1 from OG
}


def _validate_patches(rtf: dict) -> None:
    """One-time RTF patch verification: each CHARMM non-HIS patch defines the expected H.

    Prints RTF WARNING if a patch is missing or lacks its expected H indicator atom.
    Called once per run from _post(); silent when everything is correct.

    :param rtf: Parsed CHARMM36 RTF blocks.
    :type rtf: dict[str, RtfBlock]
    :rtype: None
    """
    for label, (action, h_atom) in _PATCH_H_INDICATORS.items():
        blk = _rtf_get(label, rtf)
        if blk is None:
            print(f"  RTF WARNING: PRES {label} not found - cannot verify {action} ATOM {h_atom}")
            continue
        verb = blk.verbatim
        if action == "ADD":
            found = any(
                p[0].upper() == "ATOM" and len(p) > 1 and p[1].upper() == h_atom
                for p in (l.split() for l in verb.splitlines()) if p
            )
        else:  # DELETE
            found = any(
                f"DELETE ATOM {h_atom}" in l.upper()
                for l in verb.splitlines()
            )
        if not found:
            print(f"  RTF WARNING: PRES {label}: expected {action} ATOM {h_atom} "
                  f"not found - RTF may be non-standard")


def _lys_tyr_cys_pdb_signals(prot_pdb: Path) -> dict[tuple[str, str, int], str]:
    """Read protonated PDB → state inferred from atom inventory for LYS, TYR, CYS, SER.

    RTF-derived atom indicators (RTF is the AUTHORITY for what each state looks like):
        LYS  : 3 HZ atoms → "LYS" (protonated, all HZ1/HZ2/HZ3 present)
               2 HZ atoms → "LSN" (PRES LSN deletes HZ3; neutral amine)
               other      → "UNKNOWN"
        TYR  : HH present → "TYR" (phenol OH intact)
               HH absent  → "TYRD" (phenol deprotonated; pKa ~10, rare at physiol. pH)
        CYS  : HG1 present → "CYS" (thiol SH intact)
               HG1 absent  → "CYSD" (PRES CYSD deletes HG1; thiolate; excludes SS-bonds)
        SER  : HG1 present → "SER" (hydroxyl OH intact)
               HG1 absent  → "SERD" (PRES SERD deletes HG1; pKa ~13–14, very rare)

    ASP  : HD2 present → "ASPP" (protonated ASP; PyPKA writes HD2 in the protonated PDB)
           HD2 absent  → "ASP" (deprotonated default)
    GLU  : HE2 present → "GLUP" (protonated GLU; PyPKA writes HE2 in the protonated PDB)
           HE2 absent  → "GLU" (deprotonated default)

    Note: PyPKA starts from the fully-protonated model (HD11/HD12/HD21/HD22 for ASP,
    HE11/HE12/HE21/HE22 for GLU) and removes tautomer-specific H atoms via CHARMM_protomers,
    renaming the remaining one to HD2 (ASP) or HE2 (GLU) via ffconverter atom-name maps.
    The residue name stays "ASP"/"GLU" (not "ASPP"/"GLUP") in the protonated PDB;
    the H atom presence is the signal.

    Additional CYS variants handled:
        CYX  : SS-bonded cysteine (PyPKA writes CYX, no HG1 by disulfide geometry)
               → returns "CYX" (not a pKa-driven deprotonation; skip warning)
        CYM  : Anionic cysteine in AMBER nomenclature / RESI CYM in CHARMM36 RTF
               → maps to "CYSD" (equivalent deprotonated state)
        CYSN : PyPKA internal CHARMM tautomer name for protonated CYS
               → maps to "CYS"

    Note: RTF has both RESI CYM (anionic CYS standalone residue) and PRES CYSD (patch).
    Our pipeline uses CYSD (patch approach for CHARMM-GUI). CYM in an input PDB
    (from AMBER or manual editing) is treated as the same deprotonated state.

    :param prot_pdb: PyPKA's output protonated PDB.
    :type prot_pdb: pathlib.Path
    :returns: Mapping of ``(resname, chain, resid)`` to the inferred state label.
    :rtype: dict[tuple[str, str, int], str]
    """
    atom_sets: dict[tuple[str, str, int], set[str]] = {}
    for line in prot_pdb.read_text().splitlines():
        if len(line) < 26: continue
        if not line.startswith(("ATOM", "HETATM")):
            continue
        resname = line[17:20].strip().upper()
        if resname not in ("ASP", "GLU", "LYS", "TYR", "CYS", "SER", "THR", "CYX", "CYM", "CYN"):
            continue
        chain = (line[21].strip() or "A").upper()
        try:
            resid = int(line[22:26])
        except ValueError:
            continue
        atom_sets.setdefault((resname, chain, resid), set()).add(line[12:16].strip().upper())

    result: dict[tuple[str, str, int], str] = {}
    for (rn, ch, ri), atoms in atom_sets.items():
        if rn == "ASP":
            result[(rn, ch, ri)] = "ASPP" if "HD2" in atoms else "ASP"
        elif rn == "GLU":
            result[(rn, ch, ri)] = "GLUP" if "HE2" in atoms else "GLU"
        elif rn == "LYS":
            hz = sum(1 for a in atoms if a.startswith("HZ"))
            result[(rn, ch, ri)] = "LYS" if hz >= 3 else ("LSN" if hz == 2 else "UNKNOWN")
        elif rn == "TYR":
            result[(rn, ch, ri)] = "TYR" if "HH" in atoms else "TYRD"
        elif rn in ("CYS", "CYN"):  # CYN = neutral CYS in some naming conventions
            result[(rn, ch, ri)] = "CYS" if "HG1" in atoms else "CYSD"
        elif rn == "CYX":
            result[(rn, ch, ri)] = "CYX"   # SS-bonded: not a protonation state change
        elif rn == "CYM":
            result[(rn, ch, ri)] = "CYSD"  # anionic CYS = equivalent to CYSD
        elif rn == "SER":
            result[(rn, ch, ri)] = "SER" if "HG1" in atoms else "SERD"
        elif rn == "THR":
            # RTF THR has HG1 (hydroxyl H on OG1). No PRES THRD in CHARMM36 RTF.
            # Deprotonated THR (pKa~14, very rare) = no HG1; kept as "THR" (unsupported patch).
            result[(rn, ch, ri)] = "THR" if "HG1" in atoms else "THR_nohg1"
    return result


def reconcile_non_his_from_pdb(mapped: list[MappedResidue], prot_pdb: Path,
                               rtf: dict) -> list[MappedResidue]:
    """Correct non-HIS labels using the protonated PDB atom inventory as authority.

    For each non-HIS titratable residue, the PDB written by PyPKA is the ground truth.
    If the assigned label disagrees with the PDB, the label is CORRECTED (not just warned),
    preventing wrong outputs in the reprocess path where mpt is reconstructed from a 50%
    cutoff and may differ from the original PyPKA MC result.

    Signals used per type:
      ASP : HD2 present → ASPP (protonated); absent → ASP (default)
      GLU : HE2 present → GLUP (protonated); absent → GLU (default)
      LYS : HZ≥3 → LYS (default); HZ=2 → LSN (neutral)
      TYR : HH present → TYR (default); absent → TYRD (deprotonated, no RTF patch)
      CYS : HG1 present → CYS (default); absent → CYSD
      SER : HG1 present → SER (default); absent → SERD
      THR : HG1 absent → warn only (PRES THRD absent from CHARMM36 RTF)
      CYX : skip (SS-bonded, not a pKa-driven state)
      CYM : treated as CYSD (anionic CYS in AMBER naming)

    :param mapped: Residues already mapped by :func:`map_residue` (and, if
        HIS is present, already passed through :func:`reconcile_his_from_pdb`).
    :type mapped: list[MappedResidue]
    :param prot_pdb: PyPKA's output protonated PDB (ground truth for atoms).
    :type prot_pdb: pathlib.Path
    :param rtf: Parsed CHARMM36 RTF blocks.
    :type rtf: dict[str, RtfBlock]
    :returns: ``mapped`` with non-HIS titratable labels corrected to match
        the protonated PDB's atom inventory where they disagree.
    :rtype: list[MappedResidue]
    """
    pdb_states = _lys_tyr_cys_pdb_signals(prot_pdb)
    result: list[MappedResidue] = []
    for r in mapped:
        rn = r.resname.upper()
        if rn not in ("ASP", "GLU", "LYS", "TYR", "CYS", "SER", "THR", "CYX", "CYM", "CYN"):
            result.append(r)
            continue

        pdb_lbl = pdb_states.get((rn, r.chain.upper(), r.resid))
        effective_label = "CYSD" if rn == "CYM" else r.charmm_label

        if pdb_lbl is None or pdb_lbl == "CYX":
            result.append(r)  # no signal / SS-bonded: keep label as-is
            continue

        if pdb_lbl == "THR_nohg1":
            print(f"  WARNING: THR{r.resid} {r.chain}: HG1 absent in PDB "
                  f"(deprotonated THR pKa~14) - PRES THRD not in CHARMM36 RTF, "
                  f"manual handling required")
            result.append(r)
            continue

        if pdb_lbl == "UNKNOWN":
            print(f"  WARNING: {rn}{r.resid} {r.chain}: ambiguous HZ count in PDB "
                  f"(assigned label={effective_label}); verify manually")
            result.append(r)
            continue

        if rn == "CYS" and pdb_lbl == "CYSD" and effective_label == "CYS":
            result.append(r)  # HG1 absent but label=CYS: SS-bonded, not a pKa event
            continue

        if pdb_lbl != effective_label:
            # PDB is authoritative: correct the label to match PDB
            corrected = pdb_lbl
            print(f"  NOTE: {rn}{r.resid} {r.chain}: label corrected "
                  f"{effective_label}→{corrected} (PDB atom inventory authority)")
            action = (rn == "HIS") or (corrected not in NO_ACTION)
            result.append(dataclasses.replace(
                r,
                charmm_label=corrected,
                rtf_available=_rtf_has(corrected, rtf),
                needs_action=action,
            ))
        else:
            result.append(r)
    return result


def _validate_label_chain(mapped: list[MappedResidue], rtf: dict) -> None:
    """Cross-check each action residue through: PyPKA source → _label() → RTF.

    For every titratable site, validates three things in sequence:

    1. RTF availability: the label must have a RESI or PRES entry in CHARMM36 RTF.
       If not, CHARMM-GUI cannot apply the action (flagged as NO_PATCH in the table).
    2. Patch indicator consistency: for PRES patches (ASPP, GLUP, LSN, CYSD, SERD), confirm
       the patch atom indicator (_PATCH_H_INDICATORS) is coherent with the RTF block.
    3. pct direction: ADD-patch labels require pct_protonated≥50%; DELETE-patch labels
       require pct_deprotonated≥50%. TIE cases are exempt.

    Note: mpt-vs-label source check is NOT done here for non-HIS because
    reconcile_non_his_from_pdb() already ran and may have corrected the label to match the
    PDB (which is authoritative). Checking mpt vs corrected-label would fire false alarms.

    HIS is excluded - it is handled by reconcile_his_from_pdb() (3-signal).
    Prints nothing when everything is consistent (normal case).

    :param mapped: Fully reconciled residues (post HIS and non-HIS reconciliation).
    :type mapped: list[MappedResidue]
    :param rtf: Parsed CHARMM36 RTF blocks.
    :type rtf: dict[str, RtfBlock]
    :rtype: None
    """
    for r in mapped:
        if not r.needs_action or r.resname.upper() in SKIP_TABLE:
            continue
        rn = r.resname.upper()
        if rn == "HIS":
            continue  # covered by reconcile_his_from_pdb()

        # ── Signal 2: RTF must have the label ──────────────────────────────────
        if not r.rtf_available:
            # NO_PATCH already in table; log details here
            blk = _rtf_get(r.charmm_label, rtf)
            if blk is None:
                print(f"  RTF MISSING: {rn}{r.resid} {r.chain}: label={r.charmm_label!r} "
                      f"has no RESI/PRES in CHARMM36 RTF - CHARMM-GUI cannot apply this")
            continue  # nothing more to check without RTF entry

        # ── Signal 2: _PATCH_H_INDICATORS coherence ────────────────────────────
        lbl = r.charmm_label
        if lbl in _PATCH_H_INDICATORS:
            patch_action, h_atom = _PATCH_H_INDICATORS[lbl]
            blk = _rtf_get(lbl, rtf)
            if blk:
                rtf_text = blk.verbatim.upper()
                # RTF PRES adds atoms via "ATOM <name>" (no "ADD" keyword); deletes via "DELETE ATOM <name>"
                h_in_rtf = h_atom.upper() in rtf_text
                del_consistent = (patch_action == "DELETE") == ("DELETE" in rtf_text)
                if not h_in_rtf or not del_consistent:
                    print(f"  RTF MISMATCH: {rn}{r.resid} {r.chain}: "
                          f"_PATCH_H_INDICATORS[{lbl!r}]=({patch_action},{h_atom}) "
                          f"not confirmed in RTF PRES block")

        # ── Signal 3: pct_protonated/pct_deprotonated ↔ label direction ─────────
        # ADD patch (e.g. ASPP, GLUP) = protonated state → pct_protonated must be ≥ 50%
        # DELETE patch (e.g. CYSD, SERD, LSN) = deprotonated state → pct_deprotonated ≥ 50%
        # Skip TIE cases (both ~50%) - the label is ambiguous by design.
        if not r.is_tie and lbl in _PATCH_H_INDICATORS:
            patch_action, _ = _PATCH_H_INDICATORS[lbl]
            if patch_action == "ADD" and r.pct_protonated < 50.0:
                print(f"  PCT MISMATCH: {rn}{r.resid} {r.chain}: label={lbl!r} (protonated) "
                      f"but pct_protonated={r.pct_protonated:.1f}% < 50% - check mpt/titration curve")
            elif patch_action == "DELETE" and r.pct_deprotonated < 50.0:
                print(f"  PCT MISMATCH: {rn}{r.resid} {r.chain}: label={lbl!r} (deprotonated) "
                      f"but pct_deprotonated={r.pct_deprotonated:.1f}% < 50% - check mpt/titration curve")


# ── Cross-validation (pKAI+) ─────────────────────────────────────────────────

@dataclass
class CVRow:
    """One PyPKA-vs-pKAI+ cross-validation comparison row for a titratable site.

    :ivar sign_agree: True if both methods place the site on the same side
        of 50% protonation at the target pH.
    :ivar near_ph: True if the PyPKA pKa is within 2 pH units of the target
        pH (used to compute the "near" agreement statistic, where cross-
        validation is most meaningful).
    """
    resname:str; resid:int|str; chain:str
    pka_pypka:str; pct_pypka:float
    pka_pkai:str; pct_pkai:float
    sign_agree:bool; near_ph:bool

def run_pkai(pdb_path: Path) -> list[dict]:
    """Run pKAI+ (the improved-weights model) on a PDB and collect predictions.

    Works around pKAI 1.2.0's bare (non-package-relative) internal imports
    by inserting the ``pkai`` package directory at the front of ``sys.path``
    before importing ``pkai.pKAI``.

    :param pdb_path: Path to a PDB file (any protonation state; pKAI infers
        the local heavy-atom environment, not existing H positions).
    :type pdb_path: pathlib.Path
    :returns: One dict per titratable residue: ``{"chain", "resid",
        "resname", "pka"}``.
    :rtype: list[dict]
    :raises RuntimeError: if the ``pkai`` package is not installed.
    """
    import importlib.util
    spec = importlib.util.find_spec("pkai")
    if spec is None: raise RuntimeError("pKAI not installed")
    pkai_dir = str(Path(spec.origin).parent)
    if pkai_dir not in sys.path: sys.path.insert(0, pkai_dir)
    from pkai.pKAI import pKAI as _run
    raw = _run(str(pdb_path), model_name="pKAI+")
    out: list[dict] = []
    for item in raw:
        if isinstance(item, (tuple,list)) and len(item)>=4:
            out.append({"chain":str(item[0]).strip(),"resid":int(item[1]),
                        "resname":str(item[2]).strip().upper(),"pka":float(item[3])})
        elif isinstance(item, dict): out.append(item)
    return out

def compare_cv(mapped: list[MappedResidue], pkai: list[dict], ph: float) -> list[CVRow]:
    """Join PyPKA-derived and pKAI+-derived pKa predictions per site.

    :param mapped: Fully reconciled residues from the PyPKA pipeline.
    :type mapped: list[MappedResidue]
    :param pkai: Raw pKAI+ predictions, as returned by :func:`run_pkai`.
    :type pkai: list[dict]
    :param ph: Target pH.
    :type ph: float
    :returns: One :class:`CVRow` per site present in both PyPKA and pKAI+
        output (sites reported by only one method are silently dropped).
    :rtype: list[CVRow]
    """
    idx = {}
    for rec in pkai:
        ch = str(rec.get("chain","A")).strip().upper() or "A"
        try: rid = int(rec.get("resid",-1))
        except (TypeError,ValueError): continue
        rn = str(rec.get("resname","UNK")).strip().upper()
        try: idx[(ch,rid,rn)] = float(rec.get("pka",rec.get("pKa",0)))
        except (TypeError,ValueError): pass
    rows: list[CVRow] = []
    for r in mapped:
        pkai_pka = idx.get((r.chain.upper(), r.resid, r.resname.upper()))
        if pkai_pka is None: continue
        pct_pkai = 100.0/(1.0+10**(ph-pkai_pka))
        try: near = abs(float(r.pka_str)-ph) <= 2.0
        except ValueError: near = True
        rows.append(CVRow(r.resname, r.resid, r.chain, r.pka_str, r.pct_protonated,
            f"{pkai_pka:.2f}", pct_pkai, (r.pct_protonated>=50)==(pct_pkai>=50), near))
    return rows

def write_cv_report(rows: list[CVRow], outdir: Path, ph: float) -> Path:
    """Write ``crossvalidation_report.dat`` summarizing PyPKA vs. pKAI+ agreement.

    :param rows: Comparison rows, as returned by :func:`compare_cv`.
    :type rows: list[CVRow]
    :param outdir: Output directory.
    :type outdir: pathlib.Path
    :param ph: Target pH (for the report header).
    :type ph: float
    :returns: Path to the written report.
    :rtype: pathlib.Path
    """
    near = [r for r in rows if r.near_ph]
    a_all = sum(1 for r in rows if r.sign_agree)
    a_near = sum(1 for r in near if r.sign_agree)
    col = (f"{'RESNAME':<8} {'RESID':>6} {'CHAIN':>6}  {'PKA_PYPKA':>10} "
           f"{'%PROT_PYPKA':>12}  {'PKA_PKAI+':>10} {'%PROT_PKAI+':>12}  {'AGREE':>6}")
    lines = [
        f"# pypkatool cross-validation: PyPKA vs pKAI+ | pH={ph}",
        f"# Agreement (all)       : {a_all}/{len(rows)} ({100*a_all/len(rows):.1f}%)" if rows else "# N/A",
        f"# Agreement (|pKa-pH|≤2): {a_near}/{len(near)} ({100*a_near/len(near):.1f}%)" if near else "# near: N/A",
        "-"*90, col, "-"*90,
    ]
    for r in rows:
        lines.append(f"{r.resname:<8} {str(r.resid):>6} {r.chain:>6}  "
                     f"{r.pka_pypka:>10} {r.pct_pypka:>12.1f}  "
                     f"{r.pka_pkai:>9} {r.pct_pkai:>11.1f}  {'YES' if r.sign_agree else 'NO ':>6}")
    lines.append("-"*90)
    dest = outdir / "crossvalidation_report.dat"
    dest.write_text("\n".join(lines)+"\n", encoding="utf-8")
    return dest

# ── Report generator ──────────────────────────────────────────────────────────

_DAT_HDR = """\
# pypkatool v{v} | engine=PyPKA | epsin={e} | ionicstr={i} | pH={ph}
# protein={prot} | date={date} | input={pdb}
# Reference: PyPKA DOI 10.1021/acs.jcim.0c00718 | Cross-validation: pKAI+ (ML)
#
# Only residues requiring CHARMM-GUI action are listed.
# Terminal groups (NTR/CTR) and default-state residues are omitted.
# Full tautomer detail and RTF blocks are in detail.json.
#
# PKA_PKAI+: N/P = residue not paired with pKAI+ output
# ACTION   : RESI = set in PDB Reader | PATCH = Patches panel | TIE = ambiguous
"""
_DIV = "-"*100
_COL = (f"{'RESNAME':<8} {'RESID':>6} {'CHAIN':>6}  {'PKA_PYPKA':>10}  "
        f"{'%PROT':>7} {'%DEPROT':>8}  {'PKA_PKAI+':>10}  "
        f"{'CHARMM_LABEL':<14} {'RTF':>4}  ACTION")

def _action(r: MappedResidue) -> str:
    """Classify the CHARMM-GUI action needed for one mapped residue.

    :param r: A fully reconciled residue.
    :type r: MappedResidue
    :returns: ``"TIE"`` (ambiguous), ``"NO_PATCH"`` (label has no RTF entry),
        ``"RESI"`` (HIS: pick the residue name in PDB Reader), or ``"PATCH"``
        (apply a PRES patch in the Patches panel).
    :rtype: str
    """
    if r.is_tie: return "TIE"
    if not r.rtf_available: return "NO_PATCH"
    return "RESI" if r.charmm_label in ("HSD","HSE","HSP") else "PATCH"

def write_dat(mapped: list[MappedResidue], outdir: Path, protein: str,
              ph: float, pdb_path: Path, params: dict,
              pkai_map: dict) -> Path:
    """Write ``protonation_inputs.dat``, the CHARMM-GUI-facing action table.

    Only residues that need a CHARMM-GUI action are listed (see
    :data:`NO_ACTION`/:data:`SKIP_TABLE`); default-state and terminal-group
    residues are omitted to keep the table short. Full detail for every
    titratable site is in ``protonation_inputs.json``/``detail.json``.

    :param mapped: Fully reconciled residues.
    :type mapped: list[MappedResidue]
    :param outdir: Output directory.
    :type outdir: pathlib.Path
    :param protein: Protein name/stem, used in the header.
    :type protein: str
    :param ph: Target pH.
    :type ph: float
    :param pdb_path: Input PDB path, used in the header.
    :type pdb_path: pathlib.Path
    :param params: Run parameters (``epsin``, ``ionicstr``, ...).
    :type params: dict
    :param pkai_map: ``(chain, resid, resname) -> pKAI+ pKa`` lookup.
    :type pkai_map: dict
    :returns: Path to the written file.
    :rtype: pathlib.Path
    """
    date = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    epsin = params.get("epsin", 15)
    lines = [_DAT_HDR.format(v=__version__, e=epsin, i=params.get("ionicstr",0.1),
                              ph=ph, prot=protein, date=date, pdb=pdb_path.name),
             _DIV, _COL, _DIV]
    for r in mapped:
        if not r.needs_action or r.resname.upper() in SKIP_TABLE: continue
        key = (r.chain.upper(), r.resid, r.resname.upper())
        pkai_str = f"{pkai_map[key]:.2f}" if key in pkai_map else "N/P"
        lines.append(
            f"{r.resname:<8} {str(r.resid):>6} {r.chain:>6}  {r.pka_str:>10}  "
            f"{r.pct_protonated:>7.1f} {r.pct_deprotonated:>8.1f}  {pkai_str:>10}  "
            f"{r.charmm_label:<14} {'YES' if r.rtf_available else 'NO':>4}  {_action(r)}")
    lines.append(_DIV)
    dest = outdir / "protonation_inputs.dat"
    dest.write_text("\n".join(lines)+"\n", encoding="utf-8")
    return dest

def write_json(mapped: list[MappedResidue], outdir: Path, protein: str,
               ph: float, pdb_path: Path, params: dict,
               pypka_ver: str = "unknown") -> Path:
    """Write ``protonation_inputs.json``, the machine-readable full site table.

    Unlike :func:`write_dat`, this includes every titratable site
    (including default-state and terminal residues), for downstream scripting.

    :param mapped: Fully reconciled residues.
    :type mapped: list[MappedResidue]
    :param outdir: Output directory.
    :type outdir: pathlib.Path
    :param protein: Protein name/stem.
    :type protein: str
    :param ph: Target pH.
    :type ph: float
    :param pdb_path: Input PDB path.
    :type pdb_path: pathlib.Path
    :param params: Run parameters (``epsin``, ``ionicstr``, ...).
    :type params: dict
    :param pypka_ver: Installed PyPKA version string, for provenance.
    :type pypka_ver: str
    :returns: Path to the written file.
    :rtype: pathlib.Path
    """
    dest = outdir / "protonation_inputs.json"
    data = {
        "protein": protein, "target_ph": ph, "tool": "pypkatool",
        "engine": "PyPKA", "pypka_version": pypka_ver,
        "epsin": params.get("epsin",15), "ionicstr": params.get("ionicstr",0.1),
        "rtf_source": "top_all36_prot.rtf",
        "generated": datetime.now(timezone.utc).isoformat(),
        "total_titratable": len(mapped),
        "needs_action_count": sum(1 for r in mapped if r.needs_action),
        "summary": [{"resname":r.resname,"resid":r.resid,"chain":r.chain,
                     "pka":r.pka_str,"pct_protonated":round(r.pct_protonated,2),
                     "pct_deprotonated":round(r.pct_deprotonated,2),
                     "final_label":r.charmm_label,"rtf_available":r.rtf_available,
                     "needs_action":r.needs_action,"tautomer_detail":r.tautomer_detail}
                    for r in mapped],
    }
    dest.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")
    return dest

def write_detail_json(mapped: list[MappedResidue], outdir: Path) -> Path:
    """Write ``detail.json``: per-site tautomer breakdowns and RTF blocks used.

    :param mapped: Fully reconciled residues.
    :type mapped: list[MappedResidue]
    :param outdir: Output directory.
    :type outdir: pathlib.Path
    :returns: Path to the written file.
    :rtype: pathlib.Path
    """
    rtf = _load_rtf()
    rtf_blocks: dict[str, str] = {}
    for r in mapped:
        lbl = r.charmm_label
        if lbl not in rtf_blocks:
            blk = _rtf_get(lbl, rtf)
            if blk: rtf_blocks[lbl] = blk.verbatim
    data = {
        "tautomer_detail": {
            f"{r.resname}_{r.resid}_{r.chain}": r.tautomer_detail
            for r in mapped if r.tautomer_detail and r.tautomer_detail != "-"
        },
        "rtf_blocks_used": rtf_blocks,
    }
    dest = outdir / "detail.json"
    dest.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")
    return dest

def _save_pypka_out(sites: list[SiteResult], outdir: Path, ph: float,
                    pkai_map: dict) -> None:
    """Save a raw pKa-per-residue table (PyPKA + pKAI+ columns) to ``output_pypka/pka.out``.

    :param sites: Raw (pre-mapping) PyPKA site results.
    :type sites: list[SiteResult]
    :param outdir: Output directory.
    :type outdir: pathlib.Path
    :param ph: Target pH (for the header).
    :type ph: float
    :param pkai_map: ``(chain, resid, resname) -> pKAI+ pKa`` lookup.
    :type pkai_map: dict
    :rtype: None
    """
    out = outdir / "output_pypka"
    out.mkdir(exist_ok=True)

    lines = [
        f"# pKa per residue | pH={ph} | PyPKA + pKAI+ cross-validation",
        f"# {'RESNAME':<8} {'RESID':>6} {'CHAIN':>6}  {'pKA_PyPKA':>12}  {'pKA_pKAI+':>12}",
        "# pKA_pKAI+: N/P = residue not paired with pKAI+ output",
        "-" * 56,
    ]
    for s in sites:
        pkai_pka = pkai_map.get((s.chain.upper(), s.resid, s.resname.upper()))
        pkai_str = f"{pkai_pka:.2f}" if pkai_pka is not None else "N/P"
        lines.append(
            f"  {s.resname:<8} {str(s.resid):>6} {s.chain:>6}  "
            f"{s.pka_str:>12}  {pkai_str:>12}"
        )
    (out / "pka.out").write_text("\n".join(lines) + "\n", encoding="utf-8")


def _save_protocol(params: dict, pdb_path: Path, ph: float, outdir: Path) -> None:
    """Write ``protocol.json``: the exact parameters and versions used for a run.

    Written before PyPKA is invoked (from :func:`run_pypka`), so it also
    serves as a record of run intent if PyPKA crashes mid-execution.

    :param params: Run parameters passed to ``pypka.Titration``.
    :type params: dict
    :param pdb_path: Input PDB path.
    :type pdb_path: pathlib.Path
    :param ph: Target pH.
    :type ph: float
    :param outdir: Output directory.
    :type outdir: pathlib.Path
    :rtype: None
    """
    try: pypka_ver = importlib.metadata.version("pypka")
    except Exception: pypka_ver = "unknown"
    (outdir/"protocol.json").write_text(json.dumps({
        "tool":"pypkatool","engine":"PyPKA","pypka_version":pypka_ver,
        "reference":"DOI: 10.1021/acs.jcim.0c00718",
        "parameters":{k:v for k,v in params.items() if k!="structure"},
        "input_pdb":str(pdb_path),"target_ph":ph,
        "generated":datetime.now(timezone.utc).isoformat(),
    }, indent=2), encoding="utf-8")

# ── Shared post-processing ────────────────────────────────────────────────────

def _post(mapped_raw: list[SiteResult], prot_pdb: Path, pdb_path: Path,
          outdir: Path, ph: float, params: dict) -> None:
    """Run the full mapping/reconciliation/cross-validation/report pipeline.

    Shared by both ``run`` and ``reprocess`` commands. Order matters:
    :func:`map_residue` first assigns labels from the PyPKA API result, then
    :func:`reconcile_his_from_pdb` and :func:`reconcile_non_his_from_pdb`
    correct any labels using the protonated PDB's atom inventory as ground
    truth, and only then do the validation and report-writing steps run
    against the final, reconciled labels.

    :param mapped_raw: Raw PyPKA site results (:func:`run_pypka` or
        :func:`reprocess_from_files`).
    :type mapped_raw: list[SiteResult]
    :param prot_pdb: PyPKA's protonated output PDB.
    :type prot_pdb: pathlib.Path
    :param pdb_path: Original input PDB (used for pKAI+ cross-validation).
    :type pdb_path: pathlib.Path
    :param outdir: Output directory.
    :type outdir: pathlib.Path
    :param ph: Target pH.
    :type ph: float
    :param params: Run parameters, propagated into ``protonation_inputs.dat``'s header.
    :type params: dict
    :raises FileNotFoundError: if ``pdb_path`` does not exist (pKAI+
        cross-validation is mandatory and needs the original PDB).
    :rtype: None
    """
    rtf = _load_rtf()
    mapped = [map_residue(s, rtf) for s in mapped_raw]
    mapped = reconcile_his_from_pdb(mapped, prot_pdb, rtf)          # HIS: 3-signal correction
    mapped = reconcile_non_his_from_pdb(mapped, prot_pdb, rtf)      # ASP/GLU/LYS/CYS/SER/TYR correction

    # One-time RTF patch existence check + post-correction audit
    _validate_patches(rtf)
    _validate_label_chain(mapped, rtf)  # PyPKA source → _label() → RTF → pct coherence

    n_action = sum(1 for r in mapped if r.needs_action and r.resname.upper() not in SKIP_TABLE)
    print(f"      {len(mapped)} sites mapped. {n_action} need CHARMM-GUI action.")

    if not pdb_path.exists():
        raise FileNotFoundError(
            f"PDB not found for pKAI+ cross-validation: {pdb_path}\n"
            "pKAI+ is mandatory - provide --pdb or ensure the PDB is in the outdir."
        )
    print("\n[+] Cross-validation (PyPKA vs pKAI+)...")
    pkai_recs = run_pkai(pdb_path)   # hard failure if pKAI not installed or PDB missing
    pkai_map: dict = {(r["chain"].upper(), r["resid"], r["resname"].upper()): r["pka"]
                      for r in pkai_recs}
    cv = compare_cv(mapped, pkai_recs, ph)
    cv_path = write_cv_report(cv, outdir, ph)
    agree = sum(1 for r in cv if r.sign_agree)
    print(f"      {agree}/{len(cv)} sign agreements with pKAI+.")
    print(f"      {cv_path.name}")

    try: pypka_ver = importlib.metadata.version("pypka")
    except Exception: pypka_ver = "unknown"

    # Strip whichever default-outdir prefix produced this directory, to recover
    # a clean protein name for the report header: "pypkatool_" (current),
    # "pypkatools_" / "autopypka_" (older name generations of this same tool).
    protein_name = outdir.name
    for prefix in ("pypkatool_", "pypkatools_", "autopypka_"):
        if protein_name.startswith(prefix):
            protein_name = protein_name[len(prefix):]
            break
    protein_name = protein_name.split("_pH")[0]
    write_dat(mapped, outdir, protein_name, ph, pdb_path, params, pkai_map)
    write_json(mapped, outdir, pdb_path.stem, ph, pdb_path, params, pypka_ver)
    write_detail_json(mapped, outdir)
    _save_pypka_out(mapped_raw, outdir, ph, pkai_map=pkai_map)

# ── Commands ──────────────────────────────────────────────────────────────────

def cmd_run(args: argparse.Namespace) -> None:
    """Entry point for ``pypkatool run``: full PyPKA + mapping + reports pipeline.

    :param args: Parsed CLI arguments (``pdb``, ``ph``, ``outdir``, ``ncpus``, ``epsin``).
    :type args: argparse.Namespace
    :rtype: None
    """
    pdb_path = Path(args.pdb).resolve()
    print("\n[0/2] Validating PDB...")
    validate_pdb(pdb_path)
    ph: float = args.ph
    outdir = Path(args.outdir) if args.outdir else pdb_path.parent / f"pypkatool_{pdb_path.stem}_pH{ph}"
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*60}\npypkatool | {pdb_path.stem} | pH {ph}\n{'='*60}")
    print(f"  Input: {pdb_path}  epsin={args.epsin}  ncpus={args.ncpus}")

    print("\n[1/2] Running PyPKA (PB+MC)...")
    charmm_input = getattr(args, "charmm_input", False)
    params = {**DEFAULT_PARAMS, "epsin": args.epsin, "ncpus": args.ncpus}
    if charmm_input:
        params["ffinput"] = "CHARMM"
    site_results, prot_pdb = run_pypka(pdb_path, ph, outdir, args.ncpus, args.epsin,
                                        charmm_input=charmm_input)
    print(f"      {len(site_results)} sites found. Protonated PDB → {prot_pdb.name}")

    print("\n[2/2] Mapping + reports...")
    _post(site_results, prot_pdb, pdb_path, outdir, ph, params)

    print(f"\n{'='*60}\nDone. Output: {outdir}")
    for f in sorted(outdir.iterdir()):
        print(f"  {f.name}/") if f.is_dir() else print(f"  {f.name}")
    raw_dir = outdir / "output_pypka"
    if raw_dir.exists():
        print(f"\n  output_pypka/ contents:")
        for f in sorted(raw_dir.iterdir()): print(f"    {f.name}")


def cmd_reprocess(args: argparse.Namespace) -> None:
    """Entry point for ``pypkatool reprocess``: regenerate reports without rerunning PyPKA.

    Reconstructs :class:`SiteResult` objects from a previous run's
    ``*_titration.dat`` and protonated PDB (see :func:`reprocess_from_files`),
    then runs the same mapping/reconciliation/cross-validation/report
    pipeline as ``run``. Useful when PyPKA completed the expensive PB+MC step
    but a later stage failed, or to regenerate reports at a different pH
    without rerunning PyPKA (the titration curve already covers pH 0-14).

    :param args: Parsed CLI arguments (``outdir``, ``ph``, ``pdb``, ``epsin``).
    :type args: argparse.Namespace
    :raises SystemExit: if ``outdir`` or its expected PyPKA output files are missing.
    :rtype: None
    """
    outdir = Path(args.outdir).resolve()
    if not outdir.exists(): sys.exit(f"ERROR: outdir not found: {outdir}")
    ph: float = args.ph

    # Search output_pypka/ first (new layout), then outdir root (legacy)
    raw_dir = outdir / "output_pypka"
    search = raw_dir if raw_dir.exists() else outdir
    tit_dats  = sorted(search.glob("*_titration.dat"))  or sorted(outdir.glob("*_titration.dat"))
    prot_pdbs = sorted(search.glob("*_protonated_pH*.pdb")) or sorted(outdir.glob("*_protonated_pH*.pdb"))
    if not tit_dats:  sys.exit(f"ERROR: no *_titration.dat in {outdir}")
    if not prot_pdbs: sys.exit(f"ERROR: no *_protonated_pH*.pdb in {outdir}")

    titration_dat = tit_dats[0]; prot_pdb = prot_pdbs[0]
    stem = prot_pdb.name.split("_protonated_")[0]
    pdb_path = Path(args.pdb).resolve() if args.pdb else outdir / f"{stem}.pdb"

    print(f"\n{'='*60}\npypkatool reprocess | {stem} | pH {ph}\n{'='*60}")
    print(f"  titration : {titration_dat.name}\n  prot_pdb  : {prot_pdb.name}")

    print("\n[1/2] Reconstructing sites...")
    site_results = reprocess_from_files(titration_dat, prot_pdb, ph)
    print(f"      {len(site_results)} sites reconstructed.")

    print("\n[2/2] Mapping + reports...")
    params = {**DEFAULT_PARAMS, "epsin": getattr(args,"epsin",15)}
    _post(site_results, prot_pdb, pdb_path, outdir, ph, params)

    print(f"\n{'='*60}\nDone.")
    for f in sorted(outdir.iterdir()): print(f"  {f.name}")


def cmd_fixstructure(args: argparse.Namespace) -> None:
    """Entry point for ``pypkatool fixstructure``: repair a PDB with PDBFixer.

    :param args: Parsed CLI arguments (``pdb``, ``outdir``).
    :type args: argparse.Namespace
    :rtype: None
    """
    pdb_path = Path(args.pdb).resolve()
    outdir = Path(args.outdir) if args.outdir else pdb_path.parent
    print(f"\n{'='*60}\npypkatool fixstructure | {pdb_path.name}\n{'='*60}")
    out_pdb = fix_structure(pdb_path, outdir)
    print(f"\nDone. Repaired PDB: {out_pdb}")


def main() -> None:
    """CLI entry point (console script ``pypkatool``).

    Defines the ``run`` and ``reprocess`` subcommands and dispatches to
    :func:`cmd_run` / :func:`cmd_reprocess`. Prints help and returns if no
    subcommand is given.

    :rtype: None
    """
    p = argparse.ArgumentParser(prog="pypkatool",
        description="PyPKA pKa → CHARMM-GUI protonation state pipeline")
    sub = p.add_subparsers(dest="command")

    r = sub.add_parser("run", help="Run full pipeline on one PDB")
    r.add_argument("pdb"); r.add_argument("--ph", type=float, required=True)
    r.add_argument("--outdir", default=None)
    r.add_argument("--ncpus", type=int, default=os.cpu_count() or 4)
    r.add_argument("--epsin", type=float, default=15)
    r.add_argument("--charmm-input", action="store_true", dest="charmm_input",
        help="Set ffinput=CHARMM for PDBs that already carry CHARMM "
             "protonation-state residue names (HSD/HSE/HSP/ASPP/GLUP/CYSM/...), "
             "e.g. a structure re-exported from a CHARMM-GUI PDB Reader step. "
             "No effect on a standard PDB (RCSB/AlphaFold/...); default is off.")

    rp = sub.add_parser("reprocess", help="Regenerate outputs from existing partial run")
    rp.add_argument("outdir"); rp.add_argument("--ph", type=float, required=True)
    rp.add_argument("--pdb", default=None); rp.add_argument("--epsin", type=float, default=15)

    fs = sub.add_parser("fixstructure",
        help="Repair missing atoms/internal-gap residues in a PDB with PDBFixer "
             "(separate 'pdbfixer' conda env) before running it through pypkatool")
    fs.add_argument("pdb"); fs.add_argument("--outdir", default=None)

    args = p.parse_args()
    if args.command == "run": cmd_run(args)
    elif args.command == "reprocess": cmd_reprocess(args)
    elif args.command == "fixstructure": cmd_fixstructure(args)
    else: p.print_help()

if __name__ == "__main__":
    main()
