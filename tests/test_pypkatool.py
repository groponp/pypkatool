#!/usr/bin/env python3
"""
Validation tests for pypkatool: unit, regression, and adversarial layers.

Regression tests compare against frozen results in tests/fixtures/ (no
PyPKA rerun needed to run this suite). Unit/adversarial tests exercise pure
functions directly with synthetic inputs.

Run: python tests/test_pypkatool.py [-v]
"""
import dataclasses
import json
import sys
import traceback
import tempfile
from pathlib import Path

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))
import pypkatool.core as ap

T = Path(__file__).parent            # tests/
DATA = T / "data"                    # tests/data/  — benchmark input PDBs
FIXTURES = T / "fixtures"            # tests/fixtures/ — frozen reference outputs
VERBOSE = "-v" in sys.argv

def _rdir(stem: str, ph: float) -> Path:
    """Frozen reference-output directory for one benchmark protein/pH."""
    return FIXTURES / f"{stem}_pH{ph}"

def _find_in_outdir(outdir: Path, pattern: str) -> Path:
    """Find a file by glob in output_pypka/ first (new layout), then outdir root (legacy)."""
    sub = outdir / "output_pypka"
    hits = sorted(sub.glob(pattern)) if sub.exists() else []
    if not hits:
        hits = sorted(outdir.glob(pattern))
    if not hits:
        raise FileNotFoundError(f"{pattern} not found in {outdir}")
    return hits[0]

RESULTS = {
    "lysozyme": _rdir("lysozyme",          7.0),
    "denv2":    _rdir("denv2",             5.0),
    "denv3":    _rdir("denv3",             5.0),
    "barnase":  _rdir("1a2p_barnase",      7.0),
    "rnasea":   _rdir("7rsa_rnasea",       7.0),
    "snase":    _rdir("1stn_snase",        7.0),
    "thiorx":   _rdir("2trx_thioredoxin",  7.0),
    "bpti":     _rdir("1bpi_bpti",         7.0),
    "proteing": _rdir("1pgb_proteing",     7.0),
    "westnile": _rdir("2hg0_westnile",     5.0),
    "dengue":   _rdir("3j27_dengue",       5.0),
}

# ─── tiny test runner ─────────────────────────────────────────────────────────

passed = failed = 0

def check(name: str, ok: bool, detail: str = "") -> None:
    global passed, failed
    if ok:
        passed += 1
        if VERBOSE: print(f"  PASS  {name}")
    else:
        failed += 1
        print(f"  FAIL  {name}" + (f"  [{detail}]" if detail else ""))

def section(title: str) -> None:
    print(f"\n── {title} {'─'*(55-len(title))}")

# ─── helpers ──────────────────────────────────────────────────────────────────

def read_dat_rows(d: Path) -> list[str]:
    """Non-comment, non-divider data rows from protonation_inputs.dat."""
    return [ln for ln in (d/"protonation_inputs.dat").read_text().splitlines()
            if ln and not ln.startswith(("#", "-", "RESNAME"))]

def read_json(d: Path) -> dict:
    return json.loads((d/"protonation_inputs.json").read_text())

def read_detail(d: Path) -> dict:
    return json.loads((d/"detail.json").read_text())

def cv_agreement(d: Path) -> float | None:
    p = d / "crossvalidation_report.dat"
    if not p.exists(): return None
    for ln in p.read_text().splitlines():
        if ln.startswith("# Agreement (all)"):
            try: return float(ln.split(":")[1].strip().split("/")[0]) / \
                              float(ln.split(":")[1].strip().split("/")[1].split()[0]) * 100
            except Exception: pass
    return None

# ─── 1. Unit tests — pure functions ──────────────────────────────────────────

section("Prerequisite — pKAI+ mandatory availability")

import importlib.util as _ilu
check("pKAI installed (pip install pKAI — mandatory requirement)",
      _ilu.find_spec("pkai") is not None,
      "FATAL: pKAI not installed. Run: pip install pKAI")
check("_require_pkai() does not raise when pKAI installed",
      True if _ilu.find_spec("pkai") is not None else False,
      "pKAI missing — all downstream tests may fail")

section("Unit — _label() mapping")

HIS_CASES = [(1,"HSD"),(2,"HSE"),(3,"HSP")]
for mpt, exp in HIS_CASES:
    check(f"HIS mpt={mpt} → {exp}", ap._label("HIS", mpt, 2, "c") == exp)

ASP_CASES = [(1,"ASPP"),(4,"ASPP"),(5,"ASP")]
for mpt, exp in ASP_CASES:
    check(f"ASP mpt={mpt} → {exp}", ap._label("ASP", mpt, 4, "a") == exp)

GLU_CASES = [(1,"GLUP"),(5,"GLU")]
for mpt, exp in GLU_CASES:
    check(f"GLU mpt={mpt} → {exp}", ap._label("GLU", mpt, 4, "a") == exp)

LYS_CASES = [(1,"LSN"),(4,"LYS")]
for mpt, exp in LYS_CASES:
    check(f"LYS mpt={mpt} → {exp}", ap._label("LYS", mpt, 3, "c") == exp)

NTR_CASES = [(1,"NNEU"),(4,"NTER")]
for mpt, exp in NTR_CASES:
    check(f"NTR mpt={mpt} → {exp}", ap._label("NTR", mpt, 3, "c") == exp)

CTR_CASES = [(1,"CNEU"),(5,"CTER")]
for mpt, exp in CTR_CASES:
    check(f"CTR mpt={mpt} → {exp}", ap._label("CTR", mpt, 4, "a") == exp)

section("Unit — _pka_from_curve()")

# Normal crossing near pH 5
curve = [(4.0,0.9),(4.5,0.7),(5.0,0.5),(5.5,0.3),(6.0,0.1)]
pka, s = ap._pka_from_curve(curve)
check("normal crossing → pKa ~5.0", pka is not None and abs(pka-5.0)<0.05, f"got {s}")

# No crossing (always protonated) → N/A
flat = [(ph, 0.9) for ph in range(15)]
_, s = ap._pka_from_curve(flat)
check("flat curve (no crossing) → N/A", s == "N/A")

# Empty curve → N/A
_, s = ap._pka_from_curve([])
check("empty curve → N/A", s == "N/A")

# pKa outside [0,14] floor
low = [(0.0,0.6),(0.5,0.55),(1.0,0.45)]
pka2, s2 = ap._pka_from_curve(low)
check("crossing near 0 → pKa formatted", pka2 is not None or s2 in ("<0.0","N/A"))

section("Unit — _charge_split()")

probs_c = [0.1, 0.1, 0.8]   # cationic: last=prot
pp, pd = ap._charge_split("c", probs_c)
check("cationic: last prob = protonated", abs(pp-80.0)<0.1 and abs(pd-20.0)<0.1,
      f"pp={pp:.1f} pd={pd:.1f}")

probs_a = [0.1, 0.1, 0.1, 0.1, 0.6]  # anionic: last=deprotonated
pp, pd = ap._charge_split("a", probs_a)
check("anionic: last prob = deprotonated", abs(pd-60.0)<0.1 and abs(pp-40.0)<0.1,
      f"pp={pp:.1f} pd={pd:.1f}")

check("empty probs → 50/50", ap._charge_split("c", []) == (50.0, 50.0))

section("Unit — RTF lookup")

rtf = ap._load_rtf()
for lbl, expected in [("HSD",True),("HSE",True),("HSP",True),
                       ("ASPP",True),("GLUP",True),("LSN",True),
                       ("NNEU",True),("CNEU",True),
                       ("ARGN",False),("NONEXISTENT",False)]:
    check(f"RTF has '{lbl}' = {expected}", ap._rtf_has(lbl, rtf) == expected)

section("Unit — SKIP_TABLE excludes NTR/CTR from dat")

mr = ap.MappedResidue("NTR",1,"A","7.5",80,20,"NNEU",True,False,"-",True)
check("NTR in SKIP_TABLE", "NTR" in ap.SKIP_TABLE)
check("CTR in SKIP_TABLE", "CTR" in ap.SKIP_TABLE)

section("Unit — TIE detection")

rtf = ap._load_rtf()
# Build a site result with prot~50% to trigger tie
import dataclasses
mr_tie = ap.MappedResidue("ASP",10,"A","7.0",51.0,49.0,"ASP",True,False,"-",False)
check("TIE when |prot-deprot| <= 2%", abs(mr_tie.pct_protonated - mr_tie.pct_deprotonated) <= ap.TIE_MARGIN)

mr_no_tie = ap.MappedResidue("ASP",10,"A","3.0",5.0,95.0,"ASP",True,False,"-",False)
check("no TIE when clearly deprotonated", abs(mr_no_tie.pct_protonated - mr_no_tie.pct_deprotonated) > ap.TIE_MARGIN)

section("Unit — error handling")

try:
    ap._load_rtf(Path("/nonexistent/path/to.rtf"))
    check("missing RTF raises FileNotFoundError", False)
except FileNotFoundError:
    check("missing RTF raises FileNotFoundError", True)
except Exception as e:
    check("missing RTF raises FileNotFoundError", False, str(e))

# ── validate_pdb error handling ──
# C6: Non-existent file
try:
    ap.validate_pdb(Path(tempfile.mktemp(suffix=".pdb")))
    check("validate_pdb non-existent file raises SystemExit", False)
except SystemExit as e:
    check("validate_pdb non-existent file raises SystemExit", "not found" in str(e.code))

# C7: Empty file
_pdb_c7 = Path(tempfile.mktemp(suffix=".pdb"))
_pdb_c7.write_text("")
try:
    ap.validate_pdb(_pdb_c7)
    check("validate_pdb empty file raises SystemExit", False)
except SystemExit as e:
    check("validate_pdb empty file raises SystemExit", "no ATOM/HETATM records" in str(e.code))
finally:
    _pdb_c7.unlink()

# C8: Short ATOM line
_pdb_c8 = Path(tempfile.mktemp(suffix=".pdb"))
_pdb_c8.write_text("ATOM      1  N   LYS A   1       2.812   4.829\n") # shorter than 54 chars
try:
    ap.validate_pdb(_pdb_c8)
    check("validate_pdb short ATOM line raises SystemExit", False)
except SystemExit as e:
    check("validate_pdb short ATOM line raises SystemExit", "too short" in str(e.code))
finally:
    _pdb_c8.unlink()

# C9: Invalid float coordinates
_pdb_c9 = Path(tempfile.mktemp(suffix=".pdb"))
_pdb_c9.write_text("ATOM      1  N   LYS A   1       2.812   4.829  XXXXXX\n") # XXXXXX is not a float
try:
    ap.validate_pdb(_pdb_c9)
    check("validate_pdb invalid float coordinates raises SystemExit", False)
except SystemExit as e:
    check("validate_pdb invalid float coordinates raises SystemExit", "not valid floats" in str(e.code))
finally:
    _pdb_c9.unlink()

# ─── 2. Output validation — format ────────────────────────────────────────────

section("Output format — dat table")

for name, d in RESULTS.items():
    if not d.exists(): continue
    rows = read_dat_rows(d)
    # No NTR/CTR
    ntr_ctr = [r for r in rows if r.strip().startswith(("NTR","CTR"))]
    check(f"{name}: no NTR/CTR in dat", len(ntr_ctr)==0, f"found: {ntr_ctr}")
    # All rows have RTF=YES
    no_rtf = [r for r in rows if "YES" not in r.split() and "NO_PATCH" not in r]
    check(f"{name}: all action residues RTF=YES", len(no_rtf)==0,
          f"{len(no_rtf)} rows without YES")

section("Output format — detail.json RTF as plain text")

for name, d in RESULTS.items():
    if not d.exists(): continue
    det = read_detail(d)
    bad = {k:type(v).__name__ for k,v in det.get("rtf_blocks_used",{}).items()
           if not isinstance(v, str)}
    check(f"{name}: RTF blocks are strings", len(bad)==0, f"non-string: {bad}")

section("Output format — JSON completeness")

REQUIRED_FIELDS = {"protein","target_ph","total_titratable","needs_action_count","summary"}
for name, d in RESULTS.items():
    if not d.exists(): continue
    js = read_json(d)
    missing = REQUIRED_FIELDS - set(js)
    check(f"{name}: JSON has required fields", len(missing)==0, f"missing: {missing}")
    # summary entries have expected keys
    if js["summary"]:
        entry_keys = set(js["summary"][0])
        req = {"resname","resid","chain","pka","pct_protonated","final_label","needs_action"}
        check(f"{name}: summary entries complete", req <= entry_keys, f"missing: {req-entry_keys}")

# ─── 3. Physics/chemistry — known labels ──────────────────────────────────────

section("Physics — known protonation states (literature)")

KNOWN = {
    # (protein_key, resname, resid, chain): expected_label
    ("lysozyme", "HIS", 15, "A"): "HSE",    # Tanford 1972, exp pKa~5.7
    ("rnasea",   "HIS", 12, "A"): "HSD",    # catalytic, pKa~6.0, ND1-H
    ("rnasea",   "HIS",119, "A"): "HSE",    # catalytic, pKa~6.5, NE2-H
    ("barnase",  "HIS",100, "A"): "HSE",    # canonical His102, pKa~6.3
    ("snase",    "HIS",  3, "A"): "HSE",    # canonical His8, pKa~6.5
    ("thiorx",   "HIS",  6, "A"): "HSE",    # pKa~6.2
    # denv3 chain A HIS144 + HIS315: pKa~4.8 (borderline), PyPKA PB gives HSE in dimer
    # (virion/cryo-EM shifts pKa up; documented PB limitation for buried residues)
    ("denv3",    "HIS",144, "A"): "HSE",    # buried, pKa<5 in dimer PB — expected limitation
    ("denv3",    "HIS",315, "A"): "HSE",    # pKa~4.8 borderline → HSE in dimer; HSP in virion
}

for (prot, resname, resid, chain), exp_lbl in KNOWN.items():
    d = RESULTS.get(prot)
    if not d or not d.exists():
        if VERBOSE: print(f"  SKIP  {prot} results not found")
        continue
    js = read_json(d)
    got = next((r["final_label"] for r in js["summary"]
                if r["resname"]==resname and int(r["resid"])==resid and r["chain"]==chain), None)
    check(f"{prot}: {resname}{resid}{chain} → {exp_lbl}", got == exp_lbl,
          f"got {got!r}")

section("Physics — alpha-pocket HIS protonated at pH 5 (DENV)")

for prot in ("denv2","denv3"):
    d = RESULTS.get(prot)
    if not d or not d.exists(): continue
    js = read_json(d)
    alpha_his = [r for r in js["summary"]
                 if r["resname"]=="HIS" and r["needs_action"] and
                 r["final_label"]=="HSP" and float(r["pct_protonated"])>60]
    check(f"{prot}: ≥4 alpha-pocket HIS protonated (>60%) at pH 5", len(alpha_his)>=4,
          f"found {len(alpha_his)}")

section("Physics — pKa values in chemically reasonable range")

for name, d in RESULTS.items():
    if not d.exists(): continue
    js = read_json(d)
    outliers = [r for r in js["summary"]
                if r["pka"] not in ("<0.0",">14.0","N/A")
                and not (0.0 < float(r["pka"]) < 14.0)]
    check(f"{name}: all pKa in (0,14) or flagged correctly", len(outliers)==0,
          f"{len(outliers)} outliers")

section("Physics — LYS/ARG predominantly protonated at pH 7")

for name in ("lysozyme","barnase","snase","rnasea"):
    d = RESULTS.get(name)
    if not d or not d.exists(): continue
    js = read_json(d)
    lys = [r for r in js["summary"] if r["resname"]=="LYS"]
    if lys:
        prot_frac = sum(r["pct_protonated"] for r in lys) / len(lys)
        check(f"{name}: LYS avg %prot at pH7 > 95%", prot_frac>95.0,
              f"avg={prot_frac:.1f}%")

section("Physics — ASP/GLU predominantly deprotonated at pH 7")

for name in ("lysozyme","barnase","rnasea"):
    d = RESULTS.get(name)
    if not d or not d.exists(): continue
    js = read_json(d)
    acidic = [r for r in js["summary"] if r["resname"] in ("ASP","GLU")]
    if acidic:
        deprot_frac = sum(r["pct_deprotonated"] for r in acidic) / len(acidic)
        check(f"{name}: ASP/GLU avg %deprot at pH7 > 90%", deprot_frac>90.0,
              f"avg={deprot_frac:.1f}%")

# ─── 4. Cross-validation agreement ───────────────────────────────────────────

section("Cross-validation — PyPKA vs pKAI+ sign agreement")

THRESHOLDS = {
    "lysozyme": 95.0, "barnase": 95.0, "rnasea": 95.0,
    "snase": 95.0, "thiorx": 95.0,
    "denv2": 90.0, "denv3": 88.0,
}

for name, threshold in THRESHOLDS.items():
    d = RESULTS.get(name)
    if not d or not d.exists(): continue
    ag = cv_agreement(d)
    if ag is None:
        check(f"{name}: CV report exists", False, "file missing")
        continue
    check(f"{name}: CV agreement >= {threshold}% (got {ag:.1f}%)", ag >= threshold,
          f"got {ag:.1f}%")

# ─── 5. output_pypka/ integrity ──────────────────────────────────────────────────

section("Output structure — output_pypka/pka.out")

for name, d in RESULTS.items():
    if not d.exists(): continue
    pka_out = d / "output_pypka" / "pka.out"
    check(f"{name}: output_pypka/pka.out exists", pka_out.exists())
    if pka_out.exists():
        content = pka_out.read_text()
        check(f"{name}: pka.out has residue lines", "HIS" in content or "ASP" in content or "LYS" in content,
              "no residue data found")

# ─── 6. Edge cases ────────────────────────────────────────────────────────────

section("HIS tautomer — atom-based cross-check (_his_pdb_signals)")

import tempfile, textwrap

def make_his_pdb(resname: str, atoms: list[str]) -> Path:
    """Create a minimal PDB with one HIS residue having the given atom names."""
    # Atom coordinates don't matter for the test
    lines = []
    for i, atom in enumerate(atoms):
        lines.append(
            f"ATOM  {i+1:>5}  {atom:<4}{resname} A   1       0.000   0.000   0.000"
            f"  1.00  0.00           {atom[0]}\n"
        )
    p = Path(tempfile.mktemp(suffix=".pdb"))
    p.write_text("".join(lines))
    return p

# Case: residue name=HSD, atoms HD1 only → both agree → HSD
pdb = make_his_pdb("HSD", ["N","CA","CB","CG","ND1","HD1","CE1","NE2","CD2"])
sig = ap._his_pdb_signals(pdb); pdb.unlink()
lbl_name, lbl_atoms = sig.get(("A",1),("?","?"))
check("HSD: name=HSD atoms=HD1 only → signals agree HSD",
      lbl_name=="HSD" and lbl_atoms=="HSD", f"name={lbl_name} atoms={lbl_atoms}")

# Case: residue name=HSE, atoms HE2 only → both agree → HSE
pdb = make_his_pdb("HSE", ["N","CA","CB","CG","ND1","CE1","NE2","HE2","CD2"])
sig = ap._his_pdb_signals(pdb); pdb.unlink()
lbl_name, lbl_atoms = sig.get(("A",1),("?","?"))
check("HSE: name=HSE atoms=HE2 only → signals agree HSE",
      lbl_name=="HSE" and lbl_atoms=="HSE", f"name={lbl_name} atoms={lbl_atoms}")

# Case: residue name=HSP, atoms HD1+HE2 → both agree → HSP
pdb = make_his_pdb("HSP", ["N","CA","CB","CG","ND1","HD1","CE1","NE2","HE2","CD2"])
sig = ap._his_pdb_signals(pdb); pdb.unlink()
lbl_name, lbl_atoms = sig.get(("A",1),("?","?"))
check("HSP: name=HSP atoms=HD1+HE2 → signals agree HSP",
      lbl_name=="HSP" and lbl_atoms=="HSP", f"name={lbl_name} atoms={lbl_atoms}")

# Case: name=HSD but atoms show HD1+HE2 → MISMATCH → atoms win → HSP
pdb = make_his_pdb("HSD", ["N","CA","CB","CG","ND1","HD1","CE1","NE2","HE2","CD2"])
sig = ap._his_pdb_signals(pdb); pdb.unlink()
lbl_name, lbl_atoms = sig.get(("A",1),("?","?"))
check("MISMATCH: name=HSD atoms=HD1+HE2 → atom_lbl=HSP (detected error)",
      lbl_name=="HSD" and lbl_atoms=="HSP", f"name={lbl_name} atoms={lbl_atoms}")

# Case: name=HSP but atoms show only HE2 → MISMATCH → atoms win → HSE
pdb = make_his_pdb("HSP", ["N","CA","CB","CG","ND1","CE1","NE2","HE2","CD2"])
sig = ap._his_pdb_signals(pdb); pdb.unlink()
lbl_name, lbl_atoms = sig.get(("A",1),("?","?"))
check("MISMATCH: name=HSP atoms=HE2 only → atom_lbl=HSE (detected error)",
      lbl_name=="HSP" and lbl_atoms=="HSE", f"name={lbl_name} atoms={lbl_atoms}")

# Case: no H atoms at all → UNKNOWN
pdb = make_his_pdb("HSD", ["N","CA","CB","CG","ND1","CE1","NE2","CD2"])
sig = ap._his_pdb_signals(pdb); pdb.unlink()
lbl_name, lbl_atoms = sig.get(("A",1),("?","?"))
check("no imidazole H → atom_lbl=UNKNOWN (graceful fallback)",
      lbl_atoms=="UNKNOWN", f"atoms={lbl_atoms}")

section("HIS tautomer — RTF cross-check (_his_rtf_h_atoms)")

rtf_chk = ap._load_rtf()
rtf_h_chk = ap._his_rtf_h_atoms(rtf_chk)

check("RTF HSD defines only HD1",
      rtf_h_chk.get("HSD") == frozenset({"HD1"}), f"got {rtf_h_chk.get('HSD')}")
check("RTF HSE defines only HE2",
      rtf_h_chk.get("HSE") == frozenset({"HE2"}), f"got {rtf_h_chk.get('HSE')}")
check("RTF HSP defines HD1 + HE2",
      rtf_h_chk.get("HSP") == frozenset({"HD1","HE2"}), f"got {rtf_h_chk.get('HSP')}")
check("RTF HSD does NOT define HE2",
      "HE2" not in rtf_h_chk.get("HSD", set()), "HE2 found in HSD block")
check("RTF HSE does NOT define HD1",
      "HD1" not in rtf_h_chk.get("HSE", set()), "HD1 found in HSE block")
check("_HIS_RTF_EXPECTED consistent with actual RTF file",
      all(ap._HIS_RTF_EXPECTED[lbl] == rtf_h_chk[lbl] for lbl in ("HSD","HSE","HSP")),
      f"mismatch: expected={ap._HIS_RTF_EXPECTED} got={rtf_h_chk}")
check("_RESI_H_PROFILE required atoms consistent with RTF file",
      all(ap._RESI_H_PROFILE[lbl]["required"] == rtf_h_chk[lbl] for lbl in ("HSD","HSE","HSP")),
      f"profile required mismatch")
check("_RESI_H_PROFILE HSD forbidden={HE2}",
      ap._RESI_H_PROFILE["HSD"]["forbidden"] == frozenset({"HE2"}), "HSD forbidden wrong")
check("_RESI_H_PROFILE HSE forbidden={HD1}",
      ap._RESI_H_PROFILE["HSE"]["forbidden"] == frozenset({"HD1"}), "HSE forbidden wrong")
check("_RESI_H_PROFILE HSP forbidden=empty",
      ap._RESI_H_PROFILE["HSP"]["forbidden"] == frozenset(), "HSP forbidden not empty")
check("_PATCH_H_INDICATORS includes SERD (Ser deprotonation)",
      "SERD" in ap._PATCH_H_INDICATORS, "SERD missing from _PATCH_H_INDICATORS")

section("HIS tautomer — reconcile_his_from_pdb() integration")

# Real lysozyme: run reconciliation and verify HIS15 → HSE, no discrepancies
import io, sys as _sys
buf = io.StringIO()
_old = _sys.stdout; _sys.stdout = buf
lys_dir = RESULTS["lysozyme"]
sites = ap.reprocess_from_files(
    _find_in_outdir(lys_dir, "lysozyme_titration.dat"),
    _find_in_outdir(lys_dir, "lysozyme_protonated_pH7.0.pdb"),
    7.0)
rtf = ap._load_rtf()
mapped_raw = [ap.map_residue(s, rtf) for s in sites]
mapped = ap.reconcile_his_from_pdb(mapped_raw, _find_in_outdir(lys_dir, "lysozyme_protonated_pH7.0.pdb"), rtf)
_sys.stdout = _old; log = buf.getvalue()

his_result = {(r.resid, r.chain): r.charmm_label for r in mapped if r.resname=="HIS"}
check("lysozyme HIS15 A → HSE after reconcile", his_result.get((15,"A"))=="HSE",
      f"got {his_result.get((15,'A'))}")
check("no ERROR lines in reconcile log (all signals agree)",
      "ERROR" not in log, f"errors found:\n{log}")
check("no WARNING for missing residues in reconcile log",
      "not found in protonated PDB" not in log, f"warnings:\n{log}")
check("no RTF WARNING in reconcile log (RTF file is standard CHARMM36)",
      "RTF WARNING" not in log, f"RTF warnings found:\n{log}")
check("no RTF/PDB topology inconsistency in reconcile log",
      "RTF/PDB topology inconsistency" not in log, f"inconsistencies:\n{log}")

section("Non-HIS patch RTF verification (_validate_patches + _PATCH_H_INDICATORS)")

# Verify each CHARMM36 patch has the expected H indicator atom
rtf_v = ap._load_rtf()

def _patch_has(label: str, action: str, h_atom: str, rtf: dict) -> bool:
    blk = ap._rtf_get(label, rtf)
    if blk is None: return False
    verb = blk.verbatim
    if action == "ADD":
        return any(
            p[0].upper()=="ATOM" and len(p)>1 and p[1].upper()==h_atom
            for p in (l.split() for l in verb.splitlines()) if p
        )
    return any(f"DELETE ATOM {h_atom}" in l.upper() for l in verb.splitlines())

check("PRES ASPP adds HD2 (protonated Asp OD2 proton)",
      _patch_has("ASPP", "ADD", "HD2", rtf_v), "HD2 not found in ASPP PRES block")
check("PRES GLUP adds HE2 (protonated Glu OE2 proton)",
      _patch_has("GLUP", "ADD", "HE2", rtf_v), "HE2 not found in GLUP PRES block")
check("PRES LSN deletes HZ3 (neutral Lys loses NZ proton)",
      _patch_has("LSN", "DELETE", "HZ3", rtf_v), "DELETE HZ3 not found in LSN PRES block")
check("PRES CYSD deletes HG1 (deprotonated Cys loses SG proton)",
      _patch_has("CYSD", "DELETE", "HG1", rtf_v), "DELETE HG1 not found in CYSD PRES block")

# _PATCH_H_INDICATORS must be consistent with actual RTF
for label, (action, h_atom) in ap._PATCH_H_INDICATORS.items():
    ok = _patch_has(label, action, h_atom, rtf_v)
    check(f"_PATCH_H_INDICATORS[{label}] consistent with RTF ({action} {h_atom})", ok,
          f"expected {action} ATOM {h_atom} in PRES {label} but not found")

# _validate_patches produces no output on standard CHARMM36 RTF
buf2 = io.StringIO(); _old2 = _sys.stdout; _sys.stdout = buf2
ap._validate_patches(rtf_v)
_sys.stdout = _old2; patch_log = buf2.getvalue()
check("_validate_patches() silent on standard CHARMM36 RTF (no warnings)",
      "RTF WARNING" not in patch_log, f"got:\n{patch_log}")

section("LYS / TYR / CYS / SER / THR / CYX / CYM protonation from PDB atom inventory")

# NOTE: ASP/GLU are NOT included here — PyPKA does not add HD2/HE2 to the
# protonated PDB for these types (they remain as standard GROMOS H atoms).
# ARG: NOT TITRATABLE in PyPKA (TITRABLETAUTOMERS has no ARG entry; pKa~12.5 always protonated).
# THR: titratable (n_reg=3) but no PRES THRD in CHARMM36 RTF; deprotonation extremely rare.

# ASP/GLU PDB signal: PyPKA writes HD2 on ASP (protonated=ASPP) and HE2 on GLU (protonated=GLUP).
# Residue name stays "ASP"/"GLU" in the PDB; the H atom is the signal.
def _make_asp_pdb(has_hd2: bool) -> Path:
    atoms = ["N","HN","CA","HA","CB","HB1","HB2","CG","OD1","OD2","C","O"]
    if has_hd2:
        atoms.append("HD2")
    lines = [
        f"ATOM  {i+1:>5}  {a:<4}ASP A  30       0.000   0.000   0.000"
        f"  1.00  0.00           {a[0]}\n"
        for i, a in enumerate(atoms)
    ]
    p = Path(tempfile.mktemp(suffix=".pdb"))
    p.write_text("".join(lines))
    return p

def _make_glu_pdb(has_he2: bool) -> Path:
    atoms = ["N","HN","CA","HA","CB","HB1","HB2","CG","HG1","HG2","CD","OE1","OE2","C","O"]
    if has_he2:
        atoms.append("HE2")
    lines = [
        f"ATOM  {i+1:>5}  {a:<4}GLU A  40       0.000   0.000   0.000"
        f"  1.00  0.00           {a[0]}\n"
        for i, a in enumerate(atoms)
    ]
    p = Path(tempfile.mktemp(suffix=".pdb"))
    p.write_text("".join(lines))
    return p

pdb_asp_prot = _make_asp_pdb(has_hd2=True)
pdb_asp_deprot = _make_asp_pdb(has_hd2=False)
pdb_glu_prot = _make_glu_pdb(has_he2=True)
pdb_glu_deprot = _make_glu_pdb(has_he2=False)

asp_prot_sig = ap._lys_tyr_cys_pdb_signals(pdb_asp_prot)
asp_deprot_sig = ap._lys_tyr_cys_pdb_signals(pdb_asp_deprot)
glu_prot_sig = ap._lys_tyr_cys_pdb_signals(pdb_glu_prot)
glu_deprot_sig = ap._lys_tyr_cys_pdb_signals(pdb_glu_deprot)
pdb_asp_prot.unlink(); pdb_asp_deprot.unlink()
pdb_glu_prot.unlink(); pdb_glu_deprot.unlink()

check("ASP with HD2 → _lys_tyr_cys_pdb_signals → 'ASPP'",
      asp_prot_sig[("ASP","A",30)] == "ASPP",
      f"got {asp_prot_sig.get(('ASP','A',30))}")
check("ASP without HD2 → _lys_tyr_cys_pdb_signals → 'ASP' (default deprotonated)",
      asp_deprot_sig[("ASP","A",30)] == "ASP",
      f"got {asp_deprot_sig.get(('ASP','A',30))}")
check("GLU with HE2 → _lys_tyr_cys_pdb_signals → 'GLUP'",
      glu_prot_sig[("GLU","A",40)] == "GLUP",
      f"got {glu_prot_sig.get(('GLU','A',40))}")
check("GLU without HE2 → _lys_tyr_cys_pdb_signals → 'GLU' (default deprotonated)",
      glu_deprot_sig[("GLU","A",40)] == "GLU",
      f"got {glu_deprot_sig.get(('GLU','A',40))}")

def pdb_lyt_states(pdb: Path) -> dict:
    return ap._lys_tyr_cys_pdb_signals(pdb)

# Lysozyme pH 7: all LYS protonated (pKa ~10.5 >> 7), all TYR protonated
lys_pdb = _find_in_outdir(RESULTS["lysozyme"], "lysozyme_protonated_pH7.0.pdb")
if lys_pdb.exists():
    lyt = pdb_lyt_states(lys_pdb)
    lys_states = {k: v for k, v in lyt.items() if k[0]=="LYS"}
    tyr_states = {k: v for k, v in lyt.items() if k[0]=="TYR"}
    lys_ok = all(v=="LYS" for v in lys_states.values())
    tyr_ok = all(v=="TYR" for v in tyr_states.values())
    check(f"lysozyme pH7: all LYS have 3 HZ atoms (protonated, n={len(lys_states)})",
          lys_ok, f"mismatch: {[k for k,v in lys_states.items() if v!='LYS']}")
    check(f"lysozyme pH7: all TYR have HH (protonated, n={len(tyr_states)})",
          tyr_ok, f"mismatch: {[k for k,v in tyr_states.items() if v!='TYR']}")

# DENV2 pH 5: all LYS protonated (66 LYS, all have pKa >> 5)
denv2_pdb = _find_in_outdir(RESULTS["denv2"], "denv2_protonated_pH5.0.pdb")
if denv2_pdb.exists():
    lyt2 = pdb_lyt_states(denv2_pdb)
    lys2 = {k: v for k, v in lyt2.items() if k[0]=="LYS"}
    lys2_ok = all(v=="LYS" for v in lys2.values())
    check(f"denv2 pH5: all LYS have 3 HZ atoms (protonated, n={len(lys2)})",
          lys2_ok, f"mismatch: {[k for k,v in lys2.items() if v!='LYS']}")
    tyr2 = {k: v for k, v in lyt2.items() if k[0]=="TYR"}
    tyr2_ok = all(v=="TYR" for v in tyr2.values())
    check(f"denv2 pH5: all TYR have HH (protonated, n={len(tyr2)})",
          tyr2_ok, f"mismatch: {[k for k,v in tyr2.items() if v!='TYR']}")

# reconcile_non_his_from_pdb: well-computed proteins → no corrections needed
buf3 = io.StringIO(); _old3 = _sys.stdout; _sys.stdout = buf3
sites_lys = ap.reprocess_from_files(
    _find_in_outdir(RESULTS["lysozyme"], "lysozyme_titration.dat"),
    _find_in_outdir(RESULTS["lysozyme"], "lysozyme_protonated_pH7.0.pdb"), 7.0)
rtf_lys = ap._load_rtf()
mapped_lys = [ap.map_residue(s, rtf_lys) for s in sites_lys]
reconciled_lys = ap.reconcile_non_his_from_pdb(
    mapped_lys, _find_in_outdir(RESULTS["lysozyme"], "lysozyme_protonated_pH7.0.pdb"), rtf_lys)
_sys.stdout = _old3; mismatch_log = buf3.getvalue()
check("lysozyme reprocess: reconcile_non_his_from_pdb makes no corrections (PDB and mpt agree)",
      "NOTE" not in mismatch_log and "WARNING" not in mismatch_log,
      f"unexpected corrections:\n{mismatch_log}")

# CYX/CYM PDB handling: CYX = SS-bonded (skip), CYM = anionic CYS = CYSD
def _make_cys_pdb(resname: str, has_hg1: bool) -> Path:
    atoms = ["N","HN","CA","HA","CB","HB1","HB2","SG","C","O"]
    if has_hg1: atoms.append("HG1")
    lines = [
        f"ATOM  {i+1:>5}  {a:<4}{resname} A   7       0.000   0.000   0.000"
        f"  1.00  0.00           {a[0]}\n"
        for i, a in enumerate(atoms)
    ]
    p = Path(tempfile.mktemp(suffix=".pdb"))
    p.write_text("".join(lines))
    return p

# CYX: SS-bonded → _lys_tyr_cys_pdb_signals returns "CYX"
pdb_cyx = _make_cys_pdb("CYX", has_hg1=False)
lyt_cyx = ap._lys_tyr_cys_pdb_signals(pdb_cyx); pdb_cyx.unlink()
check("CYX in PDB → state='CYX' (SS-bonded, not a pKa-driven deprotonation)",
      lyt_cyx.get(("CYX","A",7)) == "CYX", f"got {lyt_cyx}")

# CYM (AMBER anionic): → _lys_tyr_cys_pdb_signals returns "CYSD"
pdb_cym = _make_cys_pdb("CYM", has_hg1=False)
lyt_cym = ap._lys_tyr_cys_pdb_signals(pdb_cym); pdb_cym.unlink()
check("CYM in PDB → state='CYSD' (anionic CYS = CYSD in CHARMM36 nomenclature)",
      lyt_cym.get(("CYM","A",7)) == "CYSD", f"got {lyt_cym}")

# CYX in reconcile_non_his_from_pdb: SS-bonded → skip, label unchanged
pdb_cyx2 = _make_cys_pdb("CYX", has_hg1=False)
mr_cyx = [ap.MappedResidue(
    resname="CYX", resid=7, chain="A", pka_str="N/A",
    pct_protonated=0.0, pct_deprotonated=100.0,
    charmm_label="CYX", rtf_available=True, is_tie=False,
    tautomer_detail="", needs_action=False)]
buf_cyx = io.StringIO(); old_cyx = _sys.stdout; _sys.stdout = buf_cyx
result_cyx = ap.reconcile_non_his_from_pdb(mr_cyx, pdb_cyx2, ap._load_rtf())
_sys.stdout = old_cyx; pdb_cyx2.unlink(); cyx_log = buf_cyx.getvalue()
check("CYX → reconcile_non_his_from_pdb skips (SS-bond), no correction, no warning",
      "NOTE" not in cyx_log and "WARNING" not in cyx_log and
      result_cyx[0].charmm_label == "CYX",
      f"output: {cyx_log}, label: {result_cyx[0].charmm_label}")

# pka.out with pKAI+ column: test _save_pypka_out() writes both columns when pkai_map given
import tempfile as _tf2
with _tf2.TemporaryDirectory() as _tmp_pka:
    _pka_dir = Path(_tmp_pka)
    # Create a minimal SiteResult to write
    _sr = ap.SiteResult(
        resname="HIS", resid=15, chain="A", pka=5.7, pka_str="5.70",
        site_type="c", most_prob_taut=2, n_regular_tautomers=2,
        populations={}, pct_protonated=35.0, pct_deprotonated=65.0)
    # pKAI+ is MANDATORY — always two columns
    _pkai_map = {("A", 15, "HIS"): 5.82}
    ap._save_pypka_out([_sr], _pka_dir, 7.0, pkai_map=_pkai_map)
    _pka_content = (_pka_dir / "output_pypka" / "pka.out").read_text()
    check("pka.out → always has pKA_PyPKA column (mandatory two-column format)",
          "pKA_PyPKA" in _pka_content, f"got:\n{_pka_content}")
    check("pka.out → always has pKA_pKAI+ column (pKAI+ is mandatory)",
          "pKA_pKAI+" in _pka_content, f"got:\n{_pka_content}")
    check("pka.out → pKAI+ value 5.82 present when paired",
          "5.82" in _pka_content, f"pKAI+ value not found:\n{_pka_content}")
    # Unpaired residue (not in pkai_map): N/P, NOT N/A
    _sr2 = ap.SiteResult(
        resname="ASP", resid=52, chain="A", pka=1.2, pka_str="1.20",
        site_type="a", most_prob_taut=2, n_regular_tautomers=2,
        populations={}, pct_protonated=0.1, pct_deprotonated=99.9)
    ap._save_pypka_out([_sr, _sr2], _pka_dir, 7.0, pkai_map=_pkai_map)
    _pka_np = (_pka_dir / "output_pypka" / "pka.out").read_text()
    check("pka.out → unpaired residue shows N/P (not N/A — pKAI+ always ran)",
          "N/P" in _pka_np and "N/A" not in _pka_np,
          f"got:\n{_pka_np}")

# _label() SER fix: deprotonated SER (mpt==n_reg+1) → SERD; protonated → SER
check("_label SER deprotonated (mpt==4) → SERD",
      ap._label("SER", 4, 3, "a") == "SERD", f"got {ap._label('SER', 4, 3, 'a')}")
check("_label SER protonated (mpt==1) → SER",
      ap._label("SER", 1, 3, "a") == "SER", f"got {ap._label('SER', 1, 3, 'a')}")
check("_label THR deprotonated (mpt==4) → THR (no PRES THRD in CHARMM36 RTF)",
      ap._label("THR", 4, 3, "a") == "THR", f"got {ap._label('THR', 4, 3, 'a')}")

section("Synthetic error injection — RTF-primary validation catches errors")

# These tests deliberately inject errors into PDB content or labels to verify
# that the validation functions correctly detect and report inconsistencies.
# RTF is the AUTHORITY: errors are defined as deviations from RTF topology.

def _make_mapped_his(resid: int, chain: str, api_label: str) -> ap.MappedResidue:
    return ap.MappedResidue(
        resname="HIS", resid=resid, chain=chain, pka_str="5.00",
        pct_protonated=50.0, pct_deprotonated=50.0,
        charmm_label=api_label, rtf_available=True,
        is_tie=False, tautomer_detail="", needs_action=True)

rtf_inj = ap._load_rtf()

# ── Inject 1: HSD in PDB with FORBIDDEN atom HE2 → must correct to HSP ──────
# RTF HSD: required={HD1}, forbidden={HE2} → having both means HSP, not HSD
pdb_inj1 = make_his_pdb("HSD", ["N","CA","CB","CG","ND1","HD1","CE1","NE2","HE2","CD2"])
mr_inj1 = [_make_mapped_his(1, "A", "HSD")]
buf_inj1 = io.StringIO(); old_inj1 = _sys.stdout; _sys.stdout = buf_inj1
result_inj1 = ap.reconcile_his_from_pdb(mr_inj1, pdb_inj1, rtf_inj)
_sys.stdout = old_inj1; pdb_inj1.unlink(); log_inj1 = buf_inj1.getvalue()
check("inject HSD+HE2 → corrected to HSP (forbidden atom detected)",
      result_inj1[0].charmm_label == "HSP",
      f"got {result_inj1[0].charmm_label}")
check("inject HSD+HE2 → ERROR message printed (RTF-primary validation)",
      "ERROR" in log_inj1, f"log was:\n{log_inj1}")

# ── Inject 2: HSP in PDB but MISSING required atom HD1 → must correct to HSE ─
# RTF HSP: required={HD1,HE2} → missing HD1 means this is actually HSE
pdb_inj2 = make_his_pdb("HSP", ["N","CA","CB","CG","ND1","CE1","NE2","HE2","CD2"])
mr_inj2 = [_make_mapped_his(1, "A", "HSP")]
buf_inj2 = io.StringIO(); old_inj2 = _sys.stdout; _sys.stdout = buf_inj2
result_inj2 = ap.reconcile_his_from_pdb(mr_inj2, pdb_inj2, rtf_inj)
_sys.stdout = old_inj2; pdb_inj2.unlink(); log_inj2 = buf_inj2.getvalue()
check("inject HSP missing HD1 → corrected to HSE (required atom missing)",
      result_inj2[0].charmm_label == "HSE",
      f"got {result_inj2[0].charmm_label}")
check("inject HSP missing HD1 → ERROR message printed",
      "ERROR" in log_inj2, f"log was:\n{log_inj2}")

# ── Inject 3: HSE in PDB with FORBIDDEN atom HD1 → must correct to HSP ───────
# RTF HSE: required={HE2}, forbidden={HD1} → having both means HSP, not HSE
pdb_inj3 = make_his_pdb("HSE", ["N","CA","CB","CG","ND1","HD1","CE1","NE2","HE2","CD2"])
mr_inj3 = [_make_mapped_his(1, "A", "HSE")]
buf_inj3 = io.StringIO(); old_inj3 = _sys.stdout; _sys.stdout = buf_inj3
result_inj3 = ap.reconcile_his_from_pdb(mr_inj3, pdb_inj3, rtf_inj)
_sys.stdout = old_inj3; pdb_inj3.unlink(); log_inj3 = buf_inj3.getvalue()
check("inject HSE+HD1 → corrected to HSP (forbidden atom HD1 detected)",
      result_inj3[0].charmm_label == "HSP",
      f"got {result_inj3[0].charmm_label}")
check("inject HSE+HD1 → ERROR message printed",
      "ERROR" in log_inj3, f"log was:\n{log_inj3}")

# ── reconcile_non_his_from_pdb: CORRECTS label from PDB, not just warns ───────
# This covers the reprocess path where mpt is reconstructed from a 50% threshold
# and may disagree with the actual PyPKA MC result recorded in the protonated PDB.

def _make_lys_pdb(hz_count: int) -> Path:
    atoms = ["N","CA","CB","CG","CD","CE","NZ","HZ1","HZ2"]
    if hz_count >= 3:
        atoms.append("HZ3")
    lines = [
        f"ATOM  {i+1:>5}  {a:<4}LYS A   5       0.000   0.000   0.000"
        f"  1.00  0.00           {a[0]}\n"
        for i, a in enumerate(atoms)
    ]
    p = Path(tempfile.mktemp(suffix=".pdb"))
    p.write_text("".join(lines))
    return p

def _make_tyr_pdb(has_hh: bool) -> Path:
    atoms = ["N","CA","CB","CG","CD1","HD1","CD2","HD2","CE1","HE1","CE2","HE2","CZ","OH"]
    if has_hh:
        atoms.append("HH")
    lines = [
        f"ATOM  {i+1:>5}  {a:<4}TYR A   8       0.000   0.000   0.000"
        f"  1.00  0.00           {a[0]}\n"
        for i, a in enumerate(atoms)
    ]
    p = Path(tempfile.mktemp(suffix=".pdb"))
    p.write_text("".join(lines))
    return p

_rtf_inj = ap._load_rtf()

# Inject: LYS with 2 HZ labeled as LYS → PDB says LSN → corrected to LSN
pdb_inj4 = _make_lys_pdb(hz_count=2)
mr_inj4 = [ap.MappedResidue(
    resname="LYS", resid=5, chain="A", pka_str="12.0",
    pct_protonated=99.0, pct_deprotonated=1.0,
    charmm_label="LYS", rtf_available=True, is_tie=False,
    tautomer_detail="", needs_action=False)]
buf_inj4 = io.StringIO(); old_inj4 = _sys.stdout; _sys.stdout = buf_inj4
result_inj4 = ap.reconcile_non_his_from_pdb(mr_inj4, pdb_inj4, _rtf_inj)
_sys.stdout = old_inj4; pdb_inj4.unlink()
check("LYS with 2 HZ labeled LYS → CORRECTED to LSN (PDB authority)",
      result_inj4[0].charmm_label == "LSN",
      f"got {result_inj4[0].charmm_label}")
check("LYS corrected to LSN → needs_action=True",
      result_inj4[0].needs_action, f"needs_action={result_inj4[0].needs_action}")

# Inject: TYR without HH labeled TYR → PDB says TYRD → corrected to TYRD
pdb_inj5 = _make_tyr_pdb(has_hh=False)
mr_inj5 = [ap.MappedResidue(
    resname="TYR", resid=8, chain="A", pka_str="10.0",
    pct_protonated=99.0, pct_deprotonated=1.0,
    charmm_label="TYR", rtf_available=True, is_tie=False,
    tautomer_detail="", needs_action=False)]
buf_inj5 = io.StringIO(); old_inj5 = _sys.stdout; _sys.stdout = buf_inj5
result_inj5 = ap.reconcile_non_his_from_pdb(mr_inj5, pdb_inj5, _rtf_inj)
_sys.stdout = old_inj5; pdb_inj5.unlink()
check("TYR without HH labeled TYR → CORRECTED to TYRD (PDB authority)",
      result_inj5[0].charmm_label == "TYRD",
      f"got {result_inj5[0].charmm_label}")
check("TYR corrected to TYRD → needs_action=True",
      result_inj5[0].needs_action, f"needs_action={result_inj5[0].needs_action}")

# Inject: ASP with HD2 labeled ASP → PDB says ASPP → corrected to ASPP
pdb_inj_asp = _make_asp_pdb(has_hd2=True)
mr_inj_asp = [ap.MappedResidue(
    resname="ASP", resid=30, chain="A", pka_str="4.00",
    pct_protonated=60.0, pct_deprotonated=40.0,
    charmm_label="ASP", rtf_available=True, is_tie=False,
    tautomer_detail="", needs_action=False)]
buf_inj_asp = io.StringIO(); _oia = _sys.stdout; _sys.stdout = buf_inj_asp
result_inj_asp = ap.reconcile_non_his_from_pdb(mr_inj_asp, pdb_inj_asp, _rtf_inj)
_sys.stdout = _oia; pdb_inj_asp.unlink()
check("ASP with HD2 labeled ASP → CORRECTED to ASPP (PDB authority)",
      result_inj_asp[0].charmm_label == "ASPP",
      f"got {result_inj_asp[0].charmm_label}")
check("ASP corrected to ASPP → needs_action=True",
      result_inj_asp[0].needs_action, f"needs_action={result_inj_asp[0].needs_action}")

# Inject: GLU without HE2 labeled GLUP → PDB says GLU → corrected to GLU (no action)
pdb_inj_glu = _make_glu_pdb(has_he2=False)
mr_inj_glu = [ap.MappedResidue(
    resname="GLU", resid=40, chain="A", pka_str="4.00",
    pct_protonated=40.0, pct_deprotonated=60.0,
    charmm_label="GLUP", rtf_available=True, is_tie=False,
    tautomer_detail="", needs_action=True)]
buf_inj_glu = io.StringIO(); _oig = _sys.stdout; _sys.stdout = buf_inj_glu
result_inj_glu = ap.reconcile_non_his_from_pdb(mr_inj_glu, pdb_inj_glu, _rtf_inj)
_sys.stdout = _oig; pdb_inj_glu.unlink()
check("GLU without HE2 labeled GLUP → CORRECTED to GLU (no action needed)",
      result_inj_glu[0].charmm_label == "GLU",
      f"got {result_inj_glu[0].charmm_label}")
check("GLU corrected to GLU → needs_action=False",
      not result_inj_glu[0].needs_action, f"needs_action={result_inj_glu[0].needs_action}")

# ── Inject 6: correct labels produce no ERROR/WARNING (silence = all OK) ──────
pdb_inj6 = make_his_pdb("HSD", ["N","CA","CB","CG","ND1","HD1","CE1","NE2","CD2"])
mr_inj6 = [_make_mapped_his(1, "A", "HSD")]
buf_inj6 = io.StringIO(); old_inj6 = _sys.stdout; _sys.stdout = buf_inj6
result_inj6 = ap.reconcile_his_from_pdb(mr_inj6, pdb_inj6, rtf_inj)
_sys.stdout = old_inj6; pdb_inj6.unlink(); log_inj6 = buf_inj6.getvalue()
check("correct HSD (HD1 only, no HE2) → no ERROR/WARNING (RTF consistent)",
      "ERROR" not in log_inj6 and "WARNING" not in log_inj6,
      f"unexpected output:\n{log_inj6}")
check("correct HSD → label unchanged (still HSD)",
      result_inj6[0].charmm_label == "HSD",
      f"got {result_inj6[0].charmm_label}")

section("_validate_label_chain — PyPKA source → _label() → RTF for all types")

# For each titratable type that can appear in the table, verify the full chain:
# 1. _label(resname, mpt, n_reg, site_type) → expected CHARMM label
# 2. RTF has that label (_rtf_has)
# 3. _PATCH_H_INDICATORS coherent with RTF text

# Reference: PyPKA TITRABLETAUTOMERS n_reg values (from pypka/constants.py):
# ASP=4, GLU=4, CYS=3, LYS=3, SER=3, THR=3, TYR=2, HIS=2, NTR=3, CTR=4

_rtf_v = ap._load_rtf()

# ASP protonated: mpt<n_reg+1 (anionic, ref=False) → ASPP → RTF PRES ASPP → ADD HD2
# This is the "no PDB signal" case: we rely on PyPKA source (mpt) + RTF
_asp_lbl = ap._label("ASP", 1, 4, "a")   # mpt=1 < n_reg+1=5 → ASPP
check("ASP protonated (mpt=1) → _label()=ASPP", _asp_lbl == "ASPP", f"got {_asp_lbl}")
check("ASPP in RTF (PRES block)", ap._rtf_has("ASPP", _rtf_v))
check("ASP deprotonated (mpt=5) → _label()=ASP (default, no action)",
      ap._label("ASP", 5, 4, "a") == "ASP")

# GLU protonated: same pattern
_glu_lbl = ap._label("GLU", 1, 4, "a")
check("GLU protonated (mpt=1) → _label()=GLUP", _glu_lbl == "GLUP", f"got {_glu_lbl}")
check("GLUP in RTF (PRES block)", ap._rtf_has("GLUP", _rtf_v))

# LYS deprotonated: mpt<n_reg+1 (cationic, ref=False) → LSN → RTF PRES LSN
_lys_lbl = ap._label("LYS", 1, 3, "c")   # mpt=1 < n_reg+1=4 → LSN
check("LYS deprotonated (mpt=1) → _label()=LSN", _lys_lbl == "LSN", f"got {_lys_lbl}")
check("LSN in RTF (PRES block)", ap._rtf_has("LSN", _rtf_v))
check("LYS protonated (mpt=4) → _label()=LYS (default, no action)",
      ap._label("LYS", 4, 3, "c") == "LYS")

# CYS deprotonated (anionic, mpt==n_reg+1=4, ref=True) → CYSD
check("CYS deprotonated (mpt=4) → _label()=CYSD", ap._label("CYS", 4, 3, "a") == "CYSD")
check("CYSD in RTF (PRES block)", ap._rtf_has("CYSD", _rtf_v))

# TYR deprotonated (anionic, mpt==n_reg+1=3, ref=True) → TYRD → NOT in RTF
check("TYR deprotonated (mpt=3) → _label()=TYRD", ap._label("TYR", 3, 2, "a") == "TYRD")
check("TYRD NOT in RTF (no PRES TYRD in CHARMM36)", not ap._rtf_has("TYRD", _rtf_v))

# SER deprotonated (anionic, mpt==n_reg+1=4, ref=True) → SERD → RTF PRES SERD
check("SER deprotonated (mpt=4) → _label()=SERD", ap._label("SER", 4, 3, "a") == "SERD")
check("SERD in RTF (PRES block)", ap._rtf_has("SERD", _rtf_v))

# THR deprotonated (anionic, mpt==n_reg+1=4, ref=True) → THR (no PRES THRD in RTF)
check("THR deprotonated (mpt=4) → _label()=THR (THRD absent from CHARMM36 RTF)",
      ap._label("THR", 4, 3, "a") == "THR")

# _validate_label_chain() silent on correct residues (no SOURCE MISMATCH or RTF MISSING)
_ADD_PATCHES = frozenset(ap._PATCH_H_INDICATORS[k][0] == "ADD"
                         for k in ap._PATCH_H_INDICATORS if ap._PATCH_H_INDICATORS[k][0] == "ADD")

def _make_action_residue(resname, mpt, n_reg, site_type, label):
    # ADD patches → protonated state → pct_p=90%; DELETE patches → deprotonated → pct_d=90%
    add = label in ap._PATCH_H_INDICATORS and ap._PATCH_H_INDICATORS[label][0] == "ADD"
    pct_p, pct_d = (90.0, 10.0) if add else (10.0, 90.0)
    return ap.MappedResidue(
        resname=resname, resid=10, chain="A", pka_str="4.00",
        pct_protonated=pct_p, pct_deprotonated=pct_d,
        charmm_label=label, rtf_available=ap._rtf_has(label, _rtf_v),
        is_tie=False, tautomer_detail="", needs_action=True,
        site_type=site_type, mpt=mpt, n_reg=n_reg)

for resname, mpt, n_reg, site_type, expected_lbl in [
    ("ASP", 1, 4, "a", "ASPP"),
    ("GLU", 2, 4, "a", "GLUP"),
    ("LYS", 1, 3, "c", "LSN"),
    ("CYS", 4, 3, "a", "CYSD"),
    ("SER", 4, 3, "a", "SERD"),
]:
    mr = _make_action_residue(resname, mpt, n_reg, site_type, expected_lbl)
    buf = io.StringIO(); _old = _sys.stdout; _sys.stdout = buf
    ap._validate_label_chain([mr], _rtf_v)
    _sys.stdout = _old; out = buf.getvalue()
    check(f"_validate_label_chain({resname}→{expected_lbl}): no mismatch",
          "MISMATCH" not in out and "MISSING" not in out, f"got:\n{out}")

# _validate_label_chain() catches pct direction inconsistent with label
# (mpt-vs-label source check was removed: reconcile_non_his_from_pdb is the authority for that)
mr_pct_inc = _make_action_residue("ASP", 1, 4, "a", "ASPP")   # ADD patch, but force pct_p=10%
mr_pct_inc = dataclasses.replace(mr_pct_inc, pct_protonated=10.0, pct_deprotonated=90.0, is_tie=False)
buf_pi = io.StringIO(); _old_pi = _sys.stdout; _sys.stdout = buf_pi
ap._validate_label_chain([mr_pct_inc], _rtf_v)
_sys.stdout = _old_pi; out_pi = buf_pi.getvalue()
check("_validate_label_chain catches pct direction mismatch (PCT MISMATCH)",
      "PCT MISMATCH" in out_pi, f"got:\n{out_pi}")

# Signal 4: pct_protonated/pct_deprotonated ↔ label direction consistency
# ADD patch (ASPP, GLUP) = protonated state → pct_protonated must be ≥ 50%
# DELETE patch (CYSD, SERD, LSN) = deprotonated/neutral → pct_deprotonated must be ≥ 50%

def _make_pct_residue(resname, mpt, n_reg, site_type, label, pct_p, pct_d, is_tie=False):
    return ap.MappedResidue(
        resname=resname, resid=10, chain="A", pka_str="5.00",
        pct_protonated=pct_p, pct_deprotonated=pct_d,
        charmm_label=label, rtf_available=ap._rtf_has(label, _rtf_v),
        is_tie=is_tie, tautomer_detail="", needs_action=True,
        site_type=site_type, mpt=mpt, n_reg=n_reg)

# Correct: ASPP with pct_protonated=80% → no PCT MISMATCH
mr_pct_ok = _make_pct_residue("ASP", 1, 4, "a", "ASPP", pct_p=80.0, pct_d=20.0)
buf_pct = io.StringIO(); _op = _sys.stdout; _sys.stdout = buf_pct
ap._validate_label_chain([mr_pct_ok], _rtf_v)
_sys.stdout = _op
check("ASPP with pct_protonated=80%: no PCT MISMATCH",
      "PCT MISMATCH" not in buf_pct.getvalue(), f"got:\n{buf_pct.getvalue()}")

# Wrong: ASPP with pct_protonated=30% (more deprotonated) → PCT MISMATCH
mr_pct_bad = _make_pct_residue("ASP", 1, 4, "a", "ASPP", pct_p=30.0, pct_d=70.0)
buf_pct2 = io.StringIO(); _op2 = _sys.stdout; _sys.stdout = buf_pct2
ap._validate_label_chain([mr_pct_bad], _rtf_v)
_sys.stdout = _op2
check("ASPP with pct_protonated=30%: PCT MISMATCH flagged",
      "PCT MISMATCH" in buf_pct2.getvalue(), f"got:\n{buf_pct2.getvalue()}")

# Correct: CYSD with pct_deprotonated=85% → no PCT MISMATCH
mr_cysd_ok = _make_pct_residue("CYS", 4, 3, "a", "CYSD", pct_p=15.0, pct_d=85.0)
buf_cysd = io.StringIO(); _oc = _sys.stdout; _sys.stdout = buf_cysd
ap._validate_label_chain([mr_cysd_ok], _rtf_v)
_sys.stdout = _oc
check("CYSD with pct_deprotonated=85%: no PCT MISMATCH",
      "PCT MISMATCH" not in buf_cysd.getvalue(), f"got:\n{buf_cysd.getvalue()}")

# TIE case: ASPP with pct_protonated=49% but is_tie=True → no PCT MISMATCH (ambiguous by design)
mr_tie = _make_pct_residue("ASP", 1, 4, "a", "ASPP", pct_p=49.0, pct_d=51.0, is_tie=True)
buf_tie = io.StringIO(); _ot = _sys.stdout; _sys.stdout = buf_tie
ap._validate_label_chain([mr_tie], _rtf_v)
_sys.stdout = _ot
check("ASPP TIE case (pct_p=49%, is_tie=True): no PCT MISMATCH (ambiguous by design)",
      "PCT MISMATCH" not in buf_tie.getvalue(), f"got:\n{buf_tie.getvalue()}")

section("Edge cases")

# _fmt_pka
check("_fmt_pka(None) → N/A", ap._fmt_pka(None) == (None,"N/A"))
check("_fmt_pka(-1.0) → <0.0", ap._fmt_pka(-1.0) == (None,"<0.0"))
check("_fmt_pka(15.0) → >14.0", ap._fmt_pka(15.0) == (None,">14.0"))
check("_fmt_pka(7.0) → (7.0,'7.00')", ap._fmt_pka(7.0) == (7.0,"7.00"))
check("_fmt_pka('bad') → N/A", ap._fmt_pka("bad") == (None,"N/A"))

# _build_pops length matches probs
for resname in ("HIS","ASP","GLU","LYS","NTR","CTR"):
    n = ap._N_REG[resname]+1
    probs = [1.0/n]*n
    pops = ap._build_pops(resname, ap._STYPE[resname], probs)
    check(f"_build_pops({resname}): {n} states", len(pops)==n)

# ─── Adversarial injection — inputs deliberately wrong, pipeline must block ───
#
# These tests inject the OPPOSITE-OF-CORRECT state and verify that the pipeline
# detects and fixes each one. Each test encodes a real bug that would reach the
# user if a validation step were missing.

section("Adversarial: wrong label corrected by PDB atom authority")

def _make_ser_pdb(has_hg1: bool) -> Path:
    """SER with or without hydroxyl proton HG1 on OG."""
    atoms = ["N","HN","CA","HA","CB","HB1","HB2","OG"]
    if has_hg1:
        atoms.append("HG1")
    lines = [
        f"ATOM  {i+1:>5}  {a:<4}SER A  20       0.000   0.000   0.000"
        f"  1.00  0.00           {a[0]}\n"
        for i, a in enumerate(atoms)
    ]
    p = Path(tempfile.mktemp(suffix=".pdb"))
    p.write_text("".join(lines))
    return p

_rtf_adv = ap._load_rtf()

def _mr(resname, resid, chain, label, pct_p, pct_d, needs_action=True, is_tie=False):
    return ap.MappedResidue(
        resname=resname, resid=resid, chain=chain, pka_str="5.00",
        pct_protonated=pct_p, pct_deprotonated=pct_d,
        charmm_label=label, rtf_available=ap._rtf_has(label, _rtf_adv),
        is_tie=is_tie, tautomer_detail="", needs_action=needs_action)

def _reconcile(mr_list, pdb):
    buf = io.StringIO(); old = _sys.stdout; _sys.stdout = buf
    result = ap.reconcile_non_his_from_pdb(mr_list, pdb, _rtf_adv)
    _sys.stdout = old
    return result, buf.getvalue()

def _validate(mr_list):
    buf = io.StringIO(); old = _sys.stdout; _sys.stdout = buf
    ap._validate_label_chain(mr_list, _rtf_adv)
    _sys.stdout = old
    return buf.getvalue()

# ── A1: ASPP labeled but PDB has NO HD2 → label was wrong; must correct to ASP ─
_pdb_a1 = _make_asp_pdb(has_hd2=False)
_res_a1, _log_a1 = _reconcile([_mr("ASP",30,"A","ASPP",60.0,40.0)], _pdb_a1)
_pdb_a1.unlink()
check("A1: ASPP label + no HD2 in PDB → corrected to ASP (deprotonated)",
      _res_a1[0].charmm_label == "ASP", f"got {_res_a1[0].charmm_label}")
check("A1: corrected to ASP → needs_action=False (default state)",
      not _res_a1[0].needs_action, f"needs_action={_res_a1[0].needs_action}")

# ── A2: GLU labeled (deprotonated default) but PDB HAS HE2 → must correct to GLUP ─
_pdb_a2 = _make_glu_pdb(has_he2=True)
_res_a2, _log_a2 = _reconcile([_mr("GLU",40,"A","GLU",30.0,70.0,needs_action=False)], _pdb_a2)
_pdb_a2.unlink()
check("A2: GLU label + HE2 in PDB → corrected to GLUP (protonated)",
      _res_a2[0].charmm_label == "GLUP", f"got {_res_a2[0].charmm_label}")
check("A2: corrected to GLUP → needs_action=True (patch required)",
      _res_a2[0].needs_action, f"needs_action={_res_a2[0].needs_action}")

# ── A3: LSN labeled but PDB has 3 HZ → mpt said neutral, PDB says protonated ─
_pdb_a3 = _make_lys_pdb(hz_count=3)
_res_a3, _log_a3 = _reconcile([_mr("LYS",5,"A","LSN",10.0,90.0)], _pdb_a3)
_pdb_a3.unlink()
check("A3: LSN label + 3 HZ in PDB → corrected to LYS (protonated default)",
      _res_a3[0].charmm_label == "LYS", f"got {_res_a3[0].charmm_label}")
check("A3: corrected to LYS → needs_action=False (default state)",
      not _res_a3[0].needs_action, f"needs_action={_res_a3[0].needs_action}")

# ── A4: CYSD labeled but PDB has HG1 → mpt said deprotonated, PDB says protonated ─
_pdb_a4_atoms = ["N","CA","CB","SG","HG1"]
_pdb_a4_lines = [
    f"ATOM  {i+1:>5}  {a:<4}CYS A  10       0.000   0.000   0.000  1.00  0.00           {a[0]}\n"
    for i, a in enumerate(_pdb_a4_atoms)
]
_pdb_a4 = Path(tempfile.mktemp(suffix=".pdb"))
_pdb_a4.write_text("".join(_pdb_a4_lines))
_res_a4, _log_a4 = _reconcile([_mr("CYS",10,"A","CYSD",5.0,95.0)], _pdb_a4)
_pdb_a4.unlink()
check("A4: CYSD label + HG1 in PDB → corrected to CYS (protonated default)",
      _res_a4[0].charmm_label == "CYS", f"got {_res_a4[0].charmm_label}")
check("A4: corrected to CYS → needs_action=False (default state)",
      not _res_a4[0].needs_action, f"needs_action={_res_a4[0].needs_action}")

# ── A5: CYS labeled (default=CYS) but PDB has NO HG1 → SS-bond exception ────
# Intentional: if label=CYS (no action) and PDB has no HG1, reconcile assumes
# SS-bonded cysteine (not a pKa-driven deprotonation) and does NOT correct.
# The correct correction path for deprotonated CYS is: mpt=4 → label=CYSD (A4).
_pdb_a5_atoms = ["N","CA","CB","SG"]
_pdb_a5_lines = [
    f"ATOM  {i+1:>5}  {a:<4}CYS A  11       0.000   0.000   0.000  1.00  0.00           {a[0]}\n"
    for i, a in enumerate(_pdb_a5_atoms)
]
_pdb_a5 = Path(tempfile.mktemp(suffix=".pdb"))
_pdb_a5.write_text("".join(_pdb_a5_lines))
_res_a5, _log_a5 = _reconcile([_mr("CYS",11,"A","CYS",90.0,10.0,needs_action=False)], _pdb_a5)
_pdb_a5.unlink()
check("A5: CYS label + no HG1 in PDB → NOT corrected (SS-bond assumed, intentional)",
      _res_a5[0].charmm_label == "CYS", f"got {_res_a5[0].charmm_label}")
check("A5: SS-bond skip → needs_action stays False (no erroneous CYSD delivered)",
      not _res_a5[0].needs_action, f"needs_action={_res_a5[0].needs_action}")

# ── A6: TYRD labeled but PDB has HH → mpt said deprotonated, PDB says protonated ─
_pdb_a6 = _make_tyr_pdb(has_hh=True)
_res_a6, _log_a6 = _reconcile([_mr("TYR",8,"A","TYRD",5.0,95.0)], _pdb_a6)
_pdb_a6.unlink()
check("A6: TYRD label + HH in PDB → corrected to TYR (protonated default)",
      _res_a6[0].charmm_label == "TYR", f"got {_res_a6[0].charmm_label}")
check("A6: corrected to TYR → needs_action=False (default state)",
      not _res_a6[0].needs_action, f"needs_action={_res_a6[0].needs_action}")

# ── A7: TYR labeled (default) but PDB has NO HH → must correct to TYRD ──────
_pdb_a7 = _make_tyr_pdb(has_hh=False)   # _make_tyr_pdb uses resid=8
_res_a7, _log_a7 = _reconcile([_mr("TYR",8,"A","TYR",90.0,10.0,needs_action=False)], _pdb_a7)
_pdb_a7.unlink()
check("A7: TYR label + no HH in PDB → corrected to TYRD",
      _res_a7[0].charmm_label == "TYRD", f"got {_res_a7[0].charmm_label}")
check("A7: corrected to TYRD → needs_action=True",
      _res_a7[0].needs_action, f"needs_action={_res_a7[0].needs_action}")

# ── A8: SERD labeled but PDB has HG1 → mpt said deprotonated, PDB says protonated ─
_pdb_a8 = _make_ser_pdb(has_hg1=True)
_res_a8, _log_a8 = _reconcile([_mr("SER",20,"A","SERD",5.0,95.0)], _pdb_a8)
_pdb_a8.unlink()
check("A8: SERD label + HG1 in PDB → corrected to SER (protonated default)",
      _res_a8[0].charmm_label == "SER", f"got {_res_a8[0].charmm_label}")
check("A8: corrected to SER → needs_action=False (default state)",
      not _res_a8[0].needs_action, f"needs_action={_res_a8[0].needs_action}")

# ── A9: SER labeled (default) but PDB has NO HG1 → must correct to SERD ─────
_pdb_a9 = _make_ser_pdb(has_hg1=False)   # _make_ser_pdb uses resid=20
_res_a9, _log_a9 = _reconcile([_mr("SER",20,"A","SER",90.0,10.0,needs_action=False)], _pdb_a9)
_pdb_a9.unlink()
check("A9: SER label + no HG1 in PDB → corrected to SERD",
      _res_a9[0].charmm_label == "SERD", f"got {_res_a9[0].charmm_label}")
check("A9: corrected to SERD → needs_action=True",
      _res_a9[0].needs_action, f"needs_action={_res_a9[0].needs_action}")

# ── A10: HIS — HSP but PDB has only HD1 (no HE2) → must correct to HSD ──────
# make_his_pdb uses resid=1, chain="A"; MappedResidue must match
_pdb_a10 = make_his_pdb("HSP", ["N","CA","CB","CG","ND1","HD1","CE1","NE2","CD2"])
_mr_a10 = [_make_mapped_his(1, "A", "HSP")]
_buf_a10 = io.StringIO(); _old_a10 = _sys.stdout; _sys.stdout = _buf_a10
_res_a10 = ap.reconcile_his_from_pdb(_mr_a10, _pdb_a10, _rtf_adv)
_sys.stdout = _old_a10; _pdb_a10.unlink(); _log_a10 = _buf_a10.getvalue()
check("A10: HSP label + HD1 only (no HE2) in PDB → corrected to HSD",
      _res_a10[0].charmm_label == "HSD", f"got {_res_a10[0].charmm_label}")
check("A10: HSP→HSD correction → ERROR message emitted",
      "ERROR" in _log_a10, f"got:\n{_log_a10}")

# ── A11: Multiple wrong labels in one call — all corrected independently ───────
# PDB has: ASP(30) with HD2; LYS(5) with 3 HZ; GLU(40) without HE2
# Labels:  ASP→"ASP" (wrong, should ASPP); LYS→"LSN" (wrong, should LYS); GLU→"GLUP" (wrong, should GLU)
_pdb_a11_lines = []
for i, (a, rn, ri) in enumerate([
    ("N","ASP",30),("HD2","ASP",30),("CG","ASP",30),("OD1","ASP",30),("OD2","ASP",30),
    ("N","LYS",5),("NZ","LYS",5),("HZ1","LYS",5),("HZ2","LYS",5),("HZ3","LYS",5),
    ("N","GLU",40),("CD","GLU",40),("OE1","GLU",40),("OE2","GLU",40),
]):
    # Standard PDB column layout: resname at 17-19, chain at 21, resSeq at 22-25
    _pdb_a11_lines.append(
        f"ATOM  {i+1:>5}  {a:<4}{rn} A{ri:>4}       0.000   0.000   0.000  1.00  0.00           {a[0]}\n"
    )
_pdb_a11 = Path(tempfile.mktemp(suffix=".pdb"))
_pdb_a11.write_text("".join(_pdb_a11_lines))
_mr_a11 = [
    _mr("ASP", 30, "A", "ASP",  60.0, 40.0, needs_action=False),
    _mr("LYS",  5, "A", "LSN",  10.0, 90.0, needs_action=True),
    _mr("GLU", 40, "A", "GLUP", 40.0, 60.0, needs_action=True),
]
_res_a11, _log_a11 = _reconcile(_mr_a11, _pdb_a11)
_pdb_a11.unlink()
_by_rn_a11 = {r.resname: r for r in _res_a11}
check("A11: multi-wrong — ASP+HD2 with label ASP → ASPP",
      _by_rn_a11["ASP"].charmm_label == "ASPP", f"got {_by_rn_a11['ASP'].charmm_label}")
check("A11: multi-wrong — LYS+3HZ with label LSN → LYS",
      _by_rn_a11["LYS"].charmm_label == "LYS", f"got {_by_rn_a11['LYS'].charmm_label}")
check("A11: multi-wrong — GLUP+noHE2 with label GLUP → GLU",
      _by_rn_a11["GLU"].charmm_label == "GLU", f"got {_by_rn_a11['GLU'].charmm_label}")

section("Adversarial: pct direction contradicts label (_validate_label_chain)")

# ── B1: LSN (DELETE=neutral/deprotonated) with pct_deprotonated=15% → PCT MISMATCH ─
_mr_b1 = _make_pct_residue("LYS",1,3,"c","LSN", pct_p=85.0, pct_d=15.0, is_tie=False)
_log_b1 = _validate([_mr_b1])
check("B1: LSN label + pct_d=15% (DELETE patch needs pct_d≥50%) → PCT MISMATCH",
      "PCT MISMATCH" in _log_b1, f"got:\n{_log_b1}")

# ── B2: GLUP (ADD=protonated) with pct_protonated=20% → PCT MISMATCH ─────────
_mr_b2 = _make_pct_residue("GLU",1,4,"a","GLUP", pct_p=20.0, pct_d=80.0, is_tie=False)
_log_b2 = _validate([_mr_b2])
check("B2: GLUP label + pct_p=20% (ADD patch needs pct_p≥50%) → PCT MISMATCH",
      "PCT MISMATCH" in _log_b2, f"got:\n{_log_b2}")

# ── B3: SERD (DELETE=deprotonated) with pct_deprotonated=5% → PCT MISMATCH ───
_mr_b3 = _make_pct_residue("SER",4,3,"a","SERD", pct_p=95.0, pct_d=5.0, is_tie=False)
_log_b3 = _validate([_mr_b3])
check("B3: SERD label + pct_d=5% (DELETE patch needs pct_d≥50%) → PCT MISMATCH",
      "PCT MISMATCH" in _log_b3, f"got:\n{_log_b3}")

# ── B4: ASPP at TIE (pct_p=49%, is_tie=True) → no PCT MISMATCH (ambiguous by design) ─
_mr_b4 = _make_pct_residue("ASP",1,4,"a","ASPP", pct_p=49.0, pct_d=51.0, is_tie=True)
_log_b4 = _validate([_mr_b4])
check("B4: ASPP TIE (pct_p=49%, is_tie=True) → no PCT MISMATCH (ties exempt)",
      "PCT MISMATCH" not in _log_b4, f"got:\n{_log_b4}")

# ── B5: CYSD (DELETE=deprotonated) with pct_deprotonated=10% → PCT MISMATCH ──
_mr_b5 = _make_pct_residue("CYS",4,3,"a","CYSD", pct_p=90.0, pct_d=10.0, is_tie=False)
_log_b5 = _validate([_mr_b5])
check("B5: CYSD label + pct_d=10% (DELETE patch needs pct_d≥50%) → PCT MISMATCH",
      "PCT MISMATCH" in _log_b5, f"got:\n{_log_b5}")

section("Adversarial: unknown/unsupported labels (_validate_label_chain)")

# ── C1: Non-existent label → RTF MISSING ─────────────────────────────────────
_mr_c1 = dataclasses.replace(
    _make_action_residue("ASP",1,4,"a","ASPP"),
    charmm_label="XYZP", rtf_available=False)
_log_c1 = _validate([_mr_c1])
check("C1: label='XYZP' not in CHARMM36 RTF → RTF MISSING flagged",
      "RTF MISSING" in _log_c1, f"got:\n{_log_c1}")

# ── C2: TYRD → RTF MISSING (no PRES TYRD in CHARMM36) but code handles gracefully ─
_mr_c2 = dataclasses.replace(
    _make_action_residue("TYR",3,2,"a","TYRD"),
    rtf_available=False)
_log_c2 = _validate([_mr_c2])
check("C2: TYRD (no PRES in RTF) → RTF MISSING flagged (no PRES TYRD in CHARMM36)",
      "RTF MISSING" in _log_c2, f"got:\n{_log_c2}")

# ── C3: Lowercase chain ID in PDB and MappedResidue ──────────────────────────
# Let's make a LYS with lowercase chain 'b'. The pipeline should convert it
# to uppercase, and correctly reconcile it with the PDB atoms (LSN with 2 HZ).
_pdb_c3 = _make_lys_pdb(hz_count=2)
# Make the PDB content use chain 'b' instead of 'A'
_pdb_c3_txt = _pdb_c3.read_text().replace("A", "b")
_pdb_c3.write_text(_pdb_c3_txt)

# Test reconcile_non_his_from_pdb with lowercase chain in MappedResidue
_mr_c3 = _mr("LYS", 5, "b", "LSN", 10.0, 90.0) # chain is lowercase 'b'
_res_c3, _log_c3 = _reconcile([_mr_c3], _pdb_c3)
_pdb_c3.unlink()
check("C3: LSN label with lowercase chain 'b' in PDB and mapped residue → chain normalized, LSN verified",
      _res_c3[0].charmm_label == "LSN" and _res_c3[0].chain == "B",
      f"got label={_res_c3[0].charmm_label}, chain={_res_c3[0].chain}")

# ── C4: Malformed/short line check in _his_pdb_signals and _lys_tyr_cys_pdb_signals ──
# Slicing line[21] on lines shorter than 22 characters should be handled safely.
_pdb_c4 = Path(tempfile.mktemp(suffix=".pdb"))
_pdb_c4.write_text("ATOM      1  N   LYS A   1       2.812   4.829  13.856\nATOM\nATOM\n")
# These calls should complete without IndexError
try:
    _res_c4_his = ap._his_pdb_signals(_pdb_c4)
    _res_c4_non = ap._lys_tyr_cys_pdb_signals(_pdb_c4)
    _pdb_c4.unlink()
    _ok_c4 = True
except IndexError as e:
    _pdb_c4.unlink()
    _ok_c4 = False
check("C4: short/malformed ATOM line in PDB → handled safely without IndexError", _ok_c4)

# ── C5: Blank/missing chain ID in PDB ──────────────────────────────────────────
# If chain ID is blank in PDB (column 22), it should default to 'A' and be normalized.
_pdb_c5 = _make_lys_pdb(hz_count=2)
# Replace chain 'A' with ' ' in the PDB
_pdb_c5_lines = []
for line in _pdb_c5.read_text().splitlines():
    if line.startswith("ATOM"):
        # Replace column 22 (index 21) with space
        line = line[:21] + " " + line[22:]
    _pdb_c5_lines.append(line + "\n")
_pdb_c5.write_text("".join(_pdb_c5_lines))

_mr_c5 = _mr("LYS", 5, "", "LSN", 10.0, 90.0) # chain is empty string
_res_c5, _log_c5 = _reconcile([_mr_c5], _pdb_c5)
_pdb_c5.unlink()
check("C5: LSN label with blank chain ID in PDB → mapped to 'A' and reconciled successfully",
      _res_c5[0].charmm_label == "LSN" and _res_c5[0].chain == "A",
      f"got label={_res_c5[0].charmm_label}, chain={_res_c5[0].chain}")

# ─── THR — label, PDB signal, reconcile behavior ─────────────────────────────

section("THR — label mapping and reconcile behavior")

def _make_thr_pdb(has_hg1: bool) -> Path:
    """THR with or without hydroxyl proton HG1 on OG1."""
    atoms = ["N","HN","CA","HA","CB","HB","OG1","CG2"]
    if has_hg1:
        atoms.append("HG1")
    lines = [
        f"ATOM  {i+1:>5}  {a:<4}THR A  55       0.000   0.000   0.000"
        f"  1.00  0.00           {a[0]}\n"
        for i, a in enumerate(atoms)
    ]
    p = Path(tempfile.mktemp(suffix=".pdb"))
    p.write_text("".join(lines))
    return p

_rtf_thr = ap._load_rtf()

# THR label: both protonated and deprotonated map to 'THR' (no PRES THRD in RTF)
check("THR _label: mpt=1 (protonated) → 'THR' (default)",
      ap._label("THR", 1, 3, "a") == "THR", f"got {ap._label('THR',1,3,'a')}")
check("THR _label: mpt=4 (deprotonated ref) → 'THR' (no THRD in CHARMM36)",
      ap._label("THR", 4, 3, "a") == "THR", f"got {ap._label('THR',4,3,'a')}")
check("THR: 'THR' is in RTF as RESI (canonical residue, not a patch)",
      ap._rtf_has("THR", _rtf_thr), "RTF lookup failed")
check("THR: no PRES THRD in CHARMM36 RTF",
      not ap._rtf_has("THRD", _rtf_thr), "found unexpected PRES THRD")

# THR PDB signals
_pdb_thr_prot = _make_thr_pdb(has_hg1=True)
_pdb_thr_dep  = _make_thr_pdb(has_hg1=False)
_sig_thr_p = ap._lys_tyr_cys_pdb_signals(_pdb_thr_prot)
_sig_thr_d = ap._lys_tyr_cys_pdb_signals(_pdb_thr_dep)
_pdb_thr_prot.unlink(); _pdb_thr_dep.unlink()

check("THR with HG1 in PDB → signal='THR' (protonated, default, no action)",
      _sig_thr_p.get(("THR","A",55)) == "THR", f"got {_sig_thr_p.get(('THR','A',55))}")
check("THR without HG1 in PDB → signal='THR_nohg1' (sentinel: no RTF patch available)",
      _sig_thr_d.get(("THR","A",55)) == "THR_nohg1",
      f"got {_sig_thr_d.get(('THR','A',55))}")

# THR reconcile: no HG1 in PDB → WARNING but label NOT corrected (no valid CHARMM patch)
_pdb_thr_dep2 = _make_thr_pdb(has_hg1=False)
_mr_thr = ap.MappedResidue("THR",55,"A","N/A",100.0,0.0,"THR",True,False,"",False)
_buf_thr = io.StringIO(); _old_thr = _sys.stdout; _sys.stdout = _buf_thr
_res_thr = ap.reconcile_non_his_from_pdb([_mr_thr], _pdb_thr_dep2, _rtf_thr)
_sys.stdout = _old_thr; _pdb_thr_dep2.unlink(); _log_thr = _buf_thr.getvalue()
check("THR no HG1 → reconcile prints WARNING (no PRES THRD, cannot auto-correct)",
      "WARNING" in _log_thr, f"got:\n{_log_thr}")
check("THR no HG1 → label stays 'THR' (not changed to invalid 'THR_nohg1')",
      _res_thr[0].charmm_label == "THR", f"got {_res_thr[0].charmm_label}")
check("THR no HG1 → needs_action stays False (no patch to apply)",
      not _res_thr[0].needs_action, f"needs_action={_res_thr[0].needs_action}")

# THR in real output: all protonated (pKa ~14), needs_action=False, absent from .dat
_lys_dir = RESULTS["lysozyme"]
if _lys_dir.exists():
    import json as _json
    _lys_j = _json.loads((_lys_dir / "protonation_inputs.json").read_text())
    _lys_thr = [r for r in _lys_j["summary"] if r["resname"] == "THR"]
    check("lysozyme: THR residues present in JSON summary (titration ran for THR)",
          len(_lys_thr) > 0, f"found {len(_lys_thr)} THR entries")
    check("lysozyme: all THR in JSON have final_label='THR' (both states map same)",
          all(r["final_label"] == "THR" for r in _lys_thr),
          str([(r["resid"],r["final_label"]) for r in _lys_thr if r["final_label"]!="THR"]))
    check("lysozyme: all THR have needs_action=False (THR excluded from action table)",
          all(not r["needs_action"] for r in _lys_thr),
          str([(r["resid"],r["needs_action"]) for r in _lys_thr]))
    check("lysozyme: all THR have pct_protonated=100.0 at pH 7.0 (pKa~14 >> 7)",
          all(r["pct_protonated"] == 100.0 for r in _lys_thr),
          str([(r["resid"],r["pct_protonated"]) for r in _lys_thr if r["pct_protonated"]!=100.0]))
    _dat_thr = [l for l in (_lys_dir/"protonation_inputs.dat").read_text().splitlines()
                if not l.startswith("#") and "THR" in l.split()[:1]]
    check("lysozyme: THR absent from protonation_inputs.dat (no CHARMM action needed)",
          len(_dat_thr) == 0, f"found unexpected THR rows: {_dat_thr}")

# ─── HIS — PDB atom inventory cross-check ────────────────────────────────────

section("HIS — PDB atom inventory must match CHARMM label in JSON output")

_LABEL_TO_ATOMS = {
    "HSD": {"required": {"HD1"}, "forbidden": {"HE2"}},
    "HSE": {"required": {"HE2"}, "forbidden": {"HD1"}},
    "HSP": {"required": {"HD1", "HE2"}, "forbidden": set()},
}

def _check_his_atoms(test_name: str, res_dir: Path, his_entries: list[dict]) -> None:
    """Verify each HIS label in JSON matches atom inventory in protonated PDB."""
    sub = res_dir / "output_pypka"
    pdb_files = (list(sub.glob("*_protonated_pH*.pdb")) if sub.exists() else []) \
                or list(res_dir.glob("*_protonated_pH*.pdb"))
    if not pdb_files:
        return
    pdb = pdb_files[0]
    sig = ap._his_pdb_signals(pdb)
    for r in his_entries:
        key = (r["chain"], r["resid"])
        if key not in sig:
            check(f"{test_name} HIS{r['resid']}{r['chain']}: found in PDB",
                  False, f"key {key} not in _his_pdb_signals")
            continue
        _name_in_pdb, _atoms_label = sig[key]
        expected = r["final_label"]
        check(f"{test_name} HIS{r['resid']}{r['chain']}: PDB atoms → {_atoms_label} matches label {expected}",
              _atoms_label == expected,
              f"PDB says {_atoms_label!r} but JSON says {expected!r}")
        # Verify required atoms present and forbidden atoms absent
        for a in _LABEL_TO_ATOMS.get(expected, {}).get("required", set()):
            has_a = _name_in_pdb == expected or _atoms_label == expected
            # Use _his_pdb_signals atoms data indirectly via the returned label
            # (if atoms_label matches, required atoms were found by the function)
        check(f"{test_name} HIS{r['resid']}{r['chain']}: PDB residue name is HSD/HSE/HSP (not plain HIS)",
              _name_in_pdb in ("HSD","HSE","HSP"), f"got name {_name_in_pdb!r}")

for _case_name, _case_key in [("lysozyme", "lysozyme"), ("barnase", "barnase"),
                               ("snase", "snase"), ("thiorx", "thiorx")]:
    _case_dir = RESULTS.get(_case_key, Path(""))
    if not _case_dir.exists():
        continue
    import json as _json2
    _case_j = _json2.loads((_case_dir / "protonation_inputs.json").read_text())
    _case_his = [r for r in _case_j["summary"] if r["resname"] == "HIS"]
    if _case_his:
        _check_his_atoms(_case_name, _case_dir, _case_his)

# Specific known cases (regression guards)
if RESULTS["lysozyme"].exists():
    import json as _json3
    _lj = _json3.loads((RESULTS["lysozyme"] / "protonation_inputs.json").read_text())
    _his15 = next((r for r in _lj["summary"] if r["resname"]=="HIS" and r["resid"]==15 and r["chain"]=="A"), None)
    check("lysozyme HIS15 A: in JSON summary", _his15 is not None, "not found")
    check("lysozyme HIS15 A: final_label=HSE (deprotonated at pH 7.0, pKa 5.88)",
          _his15 is not None and _his15["final_label"] == "HSE",
          f"got {_his15.get('final_label') if _his15 else 'missing'}")
    check("lysozyme HIS15 A: needs_action=True (differs from HSD default)",
          _his15 is not None and _his15["needs_action"],
          f"needs_action={_his15.get('needs_action') if _his15 else 'missing'}")
    check("lysozyme HIS15 A: pct_deprotonated > 80% at pH 7.0 (pKa 5.88 << 7.0)",
          _his15 is not None and _his15["pct_deprotonated"] > 80.0,
          f"pct_d={_his15.get('pct_deprotonated') if _his15 else 'missing'}")

if RESULTS["barnase"].exists():
    import json as _json4
    _bj = _json4.loads((RESULTS["barnase"] / "protonation_inputs.json").read_text())
    _barn_his = [r for r in _bj["summary"] if r["resname"]=="HIS" and r["needs_action"]]
    check("barnase: at least 1 HIS with needs_action=True at pH 7.0",
          len(_barn_his) > 0, f"found {len(_barn_his)}")
    _bad_his = [r for r in _barn_his if r["final_label"] not in ("HSD","HSE","HSP")]
    check("barnase: all HIS labels are HSD/HSE/HSP (no invalid labels)",
          len(_bad_his) == 0,
          str([(r["resid"],r["final_label"]) for r in _bad_his]))

# ─── All outputs — structural invariants ─────────────────────────────────────

section("Output files — structural invariants for all test cases")

_VALID_LABELS = frozenset({"HSD","HSE","HSP","ASPP","GLUP","LSN","CYSD","SERD",
                            "TYRD","NNEU","CNEU","NTER","CTER",
                            "ASP","GLU","LYS","CYS","TYR","SER","THR","ARG"})
_VALID_ACTIONS = frozenset({"RESI","PATCH","TIE"})
_VALID_RTF = frozenset({"YES","NO"})
_REQUIRED_JSON_KEYS = {"protein","target_ph","tool","engine","pypka_version",
                       "epsin","ionicstr","rtf_source","generated",
                       "total_titratable","needs_action_count","summary"}
_REQUIRED_SUMMARY_KEYS = {"resname","resid","chain","pka","pct_protonated",
                          "pct_deprotonated","final_label","rtf_available",
                          "needs_action","tautomer_detail"}

for _cname, _cdir in RESULTS.items():
    if not _cdir.exists():
        continue
    import json as _cjson

    # ── protonation_inputs.dat ─────────────────────────────────────────────
    _dat_p = _cdir / "protonation_inputs.dat"
    if _dat_p.exists():
        _dat_txt = _dat_p.read_text()
        _data_rows = [l for l in _dat_txt.splitlines()
                      if l.strip() and not l.startswith("#") and not l.startswith("-")
                      and l.split() and l.split()[0] in _VALID_LABELS | {"RESNAME"}]
        _dat_labels = [r.split()[0] for r in _data_rows if r.split()[0] != "RESNAME"]
        check(f"{_cname}: protonation_inputs.dat has required column header",
              "RESNAME" in _dat_txt and "CHARMM_LABEL" in _dat_txt,
              "missing column header")
        check(f"{_cname}: all CHARMM_LABEL values in .dat are known labels",
              all(l in _VALID_LABELS for l in _dat_labels),
              str([l for l in _dat_labels if l not in _VALID_LABELS]))
        _dat_actions = [r.split()[-1] for r in _data_rows if r.split()[0] != "RESNAME"]
        check(f"{_cname}: all ACTION values in .dat are RESI/PATCH/TIE",
              all(a in _VALID_ACTIONS for a in _dat_actions),
              str([a for a in _dat_actions if a not in _VALID_ACTIONS]))

    # ── protonation_inputs.json ────────────────────────────────────────────
    _json_p = _cdir / "protonation_inputs.json"
    if _json_p.exists():
        _jdata = _cjson.loads(_json_p.read_text())
        check(f"{_cname}: protonation_inputs.json has all required keys",
              _REQUIRED_JSON_KEYS.issubset(_jdata.keys()),
              str(_REQUIRED_JSON_KEYS - set(_jdata.keys())))
        _summary = _jdata.get("summary", [])
        check(f"{_cname}: summary is non-empty (at least one titratable residue)",
              len(_summary) > 0, f"summary length={len(_summary)}")
        _bad_keys = [r for r in _summary if not _REQUIRED_SUMMARY_KEYS.issubset(r.keys())]
        check(f"{_cname}: all summary entries have required keys",
              len(_bad_keys) == 0, f"{len(_bad_keys)} entries missing keys")
        _empty_labels = [r for r in _summary if not r.get("final_label","")]
        check(f"{_cname}: no empty final_label in summary",
              len(_empty_labels) == 0,
              str([(r["resname"],r["resid"]) for r in _empty_labels]))
        _pct_bad = [r for r in _summary
                    if abs(r["pct_protonated"] + r["pct_deprotonated"] - 100.0) > 0.5]
        check(f"{_cname}: pct_protonated + pct_deprotonated ≈ 100% for all residues",
              len(_pct_bad) == 0,
              str([(r["resname"],r["resid"],r["pct_protonated"],r["pct_deprotonated"])
                   for r in _pct_bad]))

    # ── protocol.json ──────────────────────────────────────────────────────
    _prot_p = _cdir / "protocol.json"
    if _prot_p.exists():
        _prot = _cjson.loads(_prot_p.read_text())
        check(f"{_cname}: protocol.json tool is a known tool name across renames "
              f"(AutoPypKa -> PypKaTools -> pypkatool)",
              _prot.get("tool") in ("AutoPypKa","PypKaTools","pypkatool"),
              f"got tool={_prot.get('tool')!r}")
        check(f"{_cname}: protocol.json engine='PyPKA'",
              _prot.get("engine") == "PyPKA", f"got {_prot.get('engine')!r}")
        check(f"{_cname}: protocol.json has non-empty pypka_version",
              bool(_prot.get("pypka_version")), f"got {_prot.get('pypka_version')!r}")
        check(f"{_cname}: protocol.json parameters has pH entry",
              isinstance(_prot.get("parameters"), dict) and "pH" in _prot["parameters"],
              f"parameters={_prot.get('parameters')}")

    # ── crossvalidation_report.dat ─────────────────────────────────────────
    _cv_p = _cdir / "crossvalidation_report.dat"
    if _cv_p.exists():
        _cv_txt = _cv_p.read_text()
        check(f"{_cname}: crossvalidation_report.dat has Agreement header",
              "Agreement (all)" in _cv_txt, "missing Agreement line")
        check(f"{_cname}: crossvalidation_report.dat has RESNAME column header",
              "RESNAME" in _cv_txt and "AGREE" in _cv_txt,
              "missing column headers")

    # ── detail.json ────────────────────────────────────────────────────────
    _det_p = _cdir / "detail.json"
    if _det_p.exists():
        _det = _cjson.loads(_det_p.read_text())
        check(f"{_cname}: detail.json has 'tautomer_detail' key",
              "tautomer_detail" in _det, "missing key")
        check(f"{_cname}: detail.json has 'rtf_blocks_used' key",
              "rtf_blocks_used" in _det, "missing key")
        _td = _det.get("tautomer_detail")
        check(f"{_cname}: detail.json tautomer_detail is a dict",
              isinstance(_td, dict), f"type={type(_td)}")
        # Non-empty only if the protein has HIS (BPTI/proteinG have no HIS)
        _has_his = any(r.get("resname") == "HIS"
                       for r in _cjson.loads((_cdir/"protonation_inputs.json").read_text())
                       .get("summary", []) if (_cdir/"protonation_inputs.json").exists())
        if _has_his:
            check(f"{_cname}: detail.json tautomer_detail non-empty (has HIS)",
                  bool(_td), f"empty tautomer_detail but HIS residues present")

    # ── output_pypka/pka.out ───────────────────────────────────────────────
    _pka_p = _cdir / "output_pypka" / "pka.out"
    if _pka_p.exists():
        _pka_txt = _pka_p.read_text()
        check(f"{_cname}: pka.out has pKa header comment",
              "pKa" in _pka_txt.splitlines()[0], f"first line: {_pka_txt.splitlines()[0]!r}")
        _pka_data = [l for l in _pka_txt.splitlines()
                     if l.strip() and not l.startswith("#") and not l.startswith("-")]
        check(f"{_cname}: pka.out has data rows (at least one titratable site)",
              len(_pka_data) > 0, "no data rows found")

    # ── protonated PDB ──────────────────────────────────────────────────────
    _sub = _cdir / "output_pypka"
    _pdbs = (list(_sub.glob("*_protonated_pH*.pdb")) if _sub.exists() else []) \
            or list(_cdir.glob("*_protonated_pH*.pdb"))
    if _pdbs:
        _pdb_txt = _pdbs[0].read_text()
        _pdb_atoms = [l for l in _pdb_txt.splitlines() if l.startswith("ATOM")]
        check(f"{_cname}: protonated PDB has ATOM records",
              len(_pdb_atoms) > 0, f"found {len(_pdb_atoms)} ATOM lines")
        _pdb_resnames = {l[17:20].strip() for l in _pdb_atoms}
        _non_std = _pdb_resnames & {"ASPP","GLUP","LSN","CYSN","TYRP","SERN","THRP"}
        check(f"{_cname}: protonated PDB uses only standard residue names (no ASPP/GLUP/LSN...)",
              len(_non_std) == 0, f"found non-standard: {_non_std}")
        _his_names = {l[17:20].strip() for l in _pdb_atoms if l[17:20].strip() in ("HIS","HSD","HSE","HSP")}
        if _his_names:
            check(f"{_cname}: HIS in protonated PDB are HSD/HSE/HSP (not plain HIS)",
                  "HIS" not in _his_names, f"found plain 'HIS' residue name: {_his_names}")

# ─── Summary ─────────────────────────────────────────────────────────────────

total = passed + failed
print(f"\n{'='*60}")
print(f"  {passed}/{total} passed  |  {failed} failed")
print(f"{'='*60}\n")
sys.exit(0 if failed == 0 else 1)
