#!/usr/bin/env python3
"""PDBFixer worker process for pypkatool's ``fixstructure`` command.

Runs in a separate conda environment (see
``pypkatool/core.py::_find_pdbfixer_python``) because pdbfixer/OpenMM require
numpy>=2, which conflicts with delphi4py's numpy<2 pin in the main
pypkatool/PyPKA environment - the two cannot share one environment.

Policy (see README.md "Repairing fragmented structures"):
  * Missing heavy atoms are filled in for every residue that IS present in
    the structure, including residues at the very start/end of a chain
    (e.g. a missing OXT on the last residue).
  * Missing residues (gaps in the chain) are only rebuilt when the gap is
    strictly internal - i.e. there is at least one present residue both
    before and after the gap. A gap at the very start or end of a chain
    (the construct simply starts/ends there) is left untouched - this
    command never *extends* a chain, only repairs holes inside it.
  * No hydrogens are added, so the output stays heavy-atom-only, matching
    the convention of a standard deposited PDB.
  * By default, only the protein (and DNA/RNA, if present) is kept in the
    output - waters, ions, ligands, and other hetatoms (HETATM records not
    part of the polymer) are dropped via PDBFixer's own
    `removeHeterogens(keepWater=False)`. Pass --keep-hetatoms to keep them.
    This runs after chain selection, so kept hetatoms are scoped to
    whatever chains --select-chains kept (removeChains() already drops any
    hetatom belonging to a removed chain) - not the whole original
    structure. This runs before gap/atom repair and has no effect on it.

Two mutually exclusive input modes:
  * --pdb-file: repair a local PDB file as-is.
  * --pdb-id: download the structure directly from RCSB (always carries the
    official SEQRES, so internal-gap detection is reliable by construction).

Detecting internal gaps requires a reference sequence (PDBFixer compares
the chain's actual residues against it via findMissingResidues()), and that
reference comes only from SEQRES records. Confirmed directly: on a PDB with
no SEQRES, `PDBFixer(filename=...).sequences` is `[]` and
`findMissingResidues()` returns `{}` regardless of how many residues are
actually missing - i.e. it would silently report "0 gaps" on a file that
may have real ones. Rather than let that pass silently (easy to miss for
someone new to the tool, and it produces a PDB that looks repaired but
isn't), this worker treats a missing reference sequence as a hard error and
refuses to write any output.
"""
import argparse
import json
import sys

from pdbfixer import PDBFixer
from openmm.app import PDBFile


def fix(pdb_out: str, pdb_file: str | None = None, pdb_id: str | None = None,
        select_chains: list[str] | None = None, keep_hetatoms: bool = False) -> dict:
    """Repair one structure with PDBFixer and write the result to ``pdb_out``.

    Exactly one of `pdb_file`/`pdb_id` must be given (enforced by the
    mutually exclusive CLI group in :func:`main`, not by this function
    itself). See the module docstring for the full repair policy (which
    gaps get rebuilt, why a missing SEQRES is a hard error, hetatom
    handling, and why no hydrogens are added).

    Parameters
    ----------
    pdb_out : str
        Path to write the repaired PDB to.
    pdb_file : str, optional
        Local PDB file to repair as-is.
    pdb_id : str, optional
        4-character RCSB PDB code to download and repair instead of a
        local file.
    select_chains : list of str, optional
        Chain IDs to keep; every other chain is removed before repair.
        If None, all chains are kept.
    keep_hetatoms : bool, optional
        If True, keep waters/ions/ligands/other non-polymer ``HETATM``
        records (scoped to whatever `select_chains` kept). If False
        (default), drop them via ``removeHeterogens(keepWater=False)``.

    Returns
    -------
    dict
        Summary with keys ``internal_residues_added``,
        ``residues_with_missing_atoms_filled``, ``hetatoms_removed``, and
        ``terminal_gaps_left_untouched`` (a ``{"chain_<id>_pos_<index>":
        [residue names]}`` map of terminal gaps intentionally left
        unrepaired).

    Raises
    ------
    SystemExit
        If `select_chains` requests a chain not present in the structure,
        or if the structure has no ``SEQRES``-derived reference sequence
        (internal gaps cannot be reliably detected without one).
    """
    fixer = PDBFixer(filename=pdb_file) if pdb_file else PDBFixer(pdbid=pdb_id)

    if select_chains:
        present = {c.id for c in fixer.topology.chains()}
        requested = set(select_chains)
        missing = requested - present
        if missing:
            sys.exit(
                f"--select-chains requested chain(s) not found in the structure: "
                f"{sorted(missing)}. Chains present: {sorted(present)}"
            )
        to_remove = present - requested
        if to_remove:
            fixer.removeChains(chainIds=list(to_remove))

    n_hetatoms_removed = 0
    if not keep_hetatoms:
        n_hetatoms_removed = len(fixer.removeHeterogens(keepWater=False))

    if not fixer.sequences:
        source = pdb_file if pdb_file else f"RCSB entry {pdb_id}"
        sys.exit(
            f"No SEQRES records found in {source}: internal chain gaps cannot "
            f"be reliably detected (PDBFixer has no reference sequence to "
            f"compare the chain against), so this structure cannot be safely "
            f"repaired - refusing to write a PDB that would look fixed but "
            f"might still be missing internal residues. "
            + ("Re-run with --pdb-id <4-char RCSB code> instead of --pdb-file "
               "to download the official deposited structure (always includes "
               "SEQRES)." if pdb_file else
               "This RCSB entry unexpectedly has no SEQRES records; fixstructure "
               "cannot proceed for it.")
        )

    fixer.findMissingResidues()
    chain_lengths = {i: len(list(c.residues())) for i, c in enumerate(fixer.topology.chains())}
    dropped_terminal = {
        k: v for k, v in fixer.missingResidues.items()
        if not (0 < k[1] < chain_lengths[k[0]])
    }
    fixer.missingResidues = {
        k: v for k, v in fixer.missingResidues.items()
        if 0 < k[1] < chain_lengths[k[0]]
    }

    fixer.findMissingAtoms()
    n_missing_atom_residues = len(fixer.missingAtoms) + len(fixer.missingTerminals)
    n_internal_residues_added = sum(len(v) for v in fixer.missingResidues.values())

    fixer.addMissingAtoms()
    # No addMissingHydrogens() call: output stays heavy-atom-only.

    with open(pdb_out, "w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)

    return {
        "internal_residues_added": n_internal_residues_added,
        "residues_with_missing_atoms_filled": n_missing_atom_residues,
        "hetatoms_removed": n_hetatoms_removed,
        "terminal_gaps_left_untouched": {
            f"chain_{k[0]}_pos_{k[1]}": v for k, v in dropped_terminal.items()
        },
    }


def main() -> None:
    """Parse CLI arguments, run :func:`fix`, and print its summary as JSON.

    Invoked as a subprocess by ``pypkatool.core.fix_structure()`` in the
    separate ``pdbfixer`` conda environment (see
    ``pypkatool/core.py::_find_pdbfixer_python``); not intended to be run
    interactively.

    Returns
    -------
    None
    """
    p = argparse.ArgumentParser()
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--pdb-file", default=None)
    src.add_argument("--pdb-id", default=None)
    p.add_argument("--output", required=True)
    p.add_argument("--select-chains", default=None,
        help="Comma-separated chain IDs to keep, e.g. A,B,C. All other chains are dropped.")
    p.add_argument("--keep-hetatoms", action="store_true",
        help="Keep waters, ions, ligands, and other non-polymer HETATM records in the "
             "output (scoped to whatever --select-chains kept, if given). By default "
             "they are dropped (protein/DNA/RNA only).")
    args = p.parse_args()
    select_chains = args.select_chains.split(",") if args.select_chains else None
    summary = fix(args.output, pdb_file=args.pdb_file, pdb_id=args.pdb_id,
                  select_chains=select_chains, keep_hetatoms=args.keep_hetatoms)
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
