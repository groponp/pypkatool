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

Detecting internal gaps requires a reference sequence (PDBFixer compares
the chain's actual residues against it via findMissingResidues()) - by
default that comes only from SEQRES records in the input PDB. A PDB with
no SEQRES (e.g. hand-edited/stripped test files) gives PDBFixer nothing to
compare against, so internal gaps are silently invisible to it even when
real - confirmed directly: `PDBFixer(filename=...).sequences` is `[]` and
`findMissingResidues()` returns `{}` on such a file regardless of how many
residues are actually missing. --pdbid fetches the deposited SEQRES from
RCSB as a substitute reference sequence (atoms/coordinates still come from
the local file - only the reference sequence used for gap detection comes
from RCSB), which is the fix for that case.
"""
import argparse
import json
import sys

from pdbfixer import PDBFixer
from openmm.app import PDBFile


def fix(pdb_in: str, pdb_out: str, pdbid: str | None = None) -> dict:
    fixer = PDBFixer(filename=pdb_in)
    no_sequence_source = not fixer.sequences and not pdbid

    if pdbid:
        ref = PDBFixer(pdbid=pdbid)
        fixer.sequences = ref.sequences

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

    summary = {
        "internal_residues_added": n_internal_residues_added,
        "residues_with_missing_atoms_filled": n_missing_atom_residues,
        "terminal_gaps_left_untouched": {
            f"chain_{k[0]}_pos_{k[1]}": v for k, v in dropped_terminal.items()
        },
    }
    if no_sequence_source:
        summary["warning"] = (
            "No SEQRES records in the input PDB and no --pdbid given: internal "
            "chain gaps cannot be detected (PDBFixer has no reference sequence "
            "to compare against), so none were repaired even if present. Pass "
            "--pdbid <4-char RCSB code> if this structure has a deposited entry."
        )
    return summary


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--pdb", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--pdbid", default=None,
        help="4-char RCSB PDB code to fetch the reference sequence from, for "
             "detecting internal gaps when the input PDB has no SEQRES records.")
    args = p.parse_args()
    summary = fix(args.pdb, args.output, pdbid=args.pdbid)
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
