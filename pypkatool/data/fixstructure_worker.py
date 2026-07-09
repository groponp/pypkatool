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
"""
import argparse
import json
import sys

from pdbfixer import PDBFixer
from openmm.app import PDBFile


def fix(pdb_in: str, pdb_out: str) -> dict:
    fixer = PDBFixer(filename=pdb_in)

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
        "terminal_gaps_left_untouched": {
            f"chain_{k[0]}_pos_{k[1]}": v for k, v in dropped_terminal.items()
        },
    }


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--pdb", required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()
    summary = fix(args.pdb, args.output)
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
