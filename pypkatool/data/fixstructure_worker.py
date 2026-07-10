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
    output - waters, ions, ligands, and other heterogens (HETATM records
    not part of the polymer) are dropped via PDBFixer's own
    `removeHeterogens(keepWater=False)`. Pass --keep-heterogens to keep
    them. This runs before gap/atom repair and has no effect on it.

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
        select_chains: list[str] | None = None, keep_heterogens: bool = False) -> dict:
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

    n_heterogens_removed = 0
    if not keep_heterogens:
        n_heterogens_removed = len(fixer.removeHeterogens(keepWater=False))

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
        "heterogens_removed": n_heterogens_removed,
        "terminal_gaps_left_untouched": {
            f"chain_{k[0]}_pos_{k[1]}": v for k, v in dropped_terminal.items()
        },
    }


def main() -> None:
    p = argparse.ArgumentParser()
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--pdb-file", default=None)
    src.add_argument("--pdb-id", default=None)
    p.add_argument("--output", required=True)
    p.add_argument("--select-chains", default=None,
        help="Comma-separated chain IDs to keep, e.g. A,B,C. All other chains are dropped.")
    p.add_argument("--keep-heterogens", action="store_true",
        help="Keep waters, ions, ligands, and other non-polymer HETATM records in the "
             "output. By default they are dropped (protein/DNA/RNA only).")
    args = p.parse_args()
    select_chains = args.select_chains.split(",") if args.select_chains else None
    summary = fix(args.output, pdb_file=args.pdb_file, pdb_id=args.pdb_id,
                  select_chains=select_chains, keep_heterogens=args.keep_heterogens)
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
