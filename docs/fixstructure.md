# Repairing fragmented structures

PyPKA does not repair a structurally incomplete PDB - it requires every
present residue to have its expected atom set. A truncated terminus missing
its carboxyl oxygens, or a residue missing backbone/side-chain atoms, will
make the underlying preprocessing step fail (see {doc}`failure_modes`). For
inputs like this, run `pypkatool fixstructure` first to produce a repaired
PDB, then feed that into `pypkatool run`.

`fixstructure` wraps [PDBFixer](https://github.com/openmm/pdbfixer) with a
fixed policy, chosen to repair damage without inventing biology that was
never part of the deposited structure:

1. **Missing heavy atoms are always filled in**, for any residue that is
   present in the structure - including residues at the very start or end of
   a chain (e.g. a C-terminal residue missing its `OXT`).
2. **Missing residues (a gap in the chain) are rebuilt only when the gap is
   internal** - i.e. there is at least one present residue on both sides of
   it. A gap at the very start or end of a chain is left untouched:
   `fixstructure` repairs holes inside a chain, it does not extend a chain
   to match a SEQRES record or a reference sequence.
3. **No hydrogens are added.** The output stays heavy-atom-only, the same
   convention as a standard deposited PDB (and the form PyPKA's own cleaning
   step expects as input).
4. **Only the protein (and DNA/RNA, if present) is kept by default.**
   Waters, ions, ligands, and any other heterogen are dropped - pass
   `--keep-heterogens` to keep them. This is independent of the atom/gap
   repair above; it just controls what ends up in the output file.

| Input residue state | Chain position | Result |
|---|---|---|
| Present, missing some atoms (e.g. `OXT`, side-chain atoms) | Anywhere, including chain ends | Atoms filled in |
| Completely absent (gap in numbering) | Internal (residues present on both sides) | Residue rebuilt |
| Completely absent (gap in numbering) | At the start or end of a chain | Left untouched - chain is not extended |

## Two input modes

Exactly one of `--pdb-file` or `--pdb-id` is required:

```bash
# Repair a local file as-is
pypkatool fixstructure --pdb-file my_protein.pdb --outdir results/
# -> results/my_protein_fixed.pdb

# Download directly from RCSB and repair that instead
pypkatool fixstructure --pdb-id 7A3S --outdir results/
# -> results/7A3S_fixed.pdb

pypkatool run results/my_protein_fixed.pdb --pH 7.0
```

`--select-chains A,B,C` keeps only the listed chains (everything else is
dropped) before repair, for either input mode - useful for a large
multi-copy deposition where only some chains/subunits are wanted.

## Non-protein atoms are dropped by default

A deposited structure commonly ships with crystallographic waters,
cryoprotectant/buffer ions, and bound ligands or glycans (`HETATM` records
that are not part of the protein/DNA/RNA polymer). `fixstructure` drops all
of these by default, keeping only the repaired protein - the intent is a
clean structure ready for `pypkatool run`, not a faithful copy of everything
in the deposition. Pass `--keep-heterogens` to keep them instead:

```bash
pypkatool fixstructure --pdb-id 7A3S --keep-heterogens --outdir results/
```

AlphaFold2/ColabFold predictions are heavy-atom-only (no hydrogens) and, for
the residue range you gave the model, have no internal gaps or missing
atoms by construction - `fixstructure` is normally unnecessary for them.
It is intended for experimentally determined structures with genuine
crystallographic disorder or truncated regions.

## Internal gaps require a reference sequence - why `--pdb-file` can refuse to run

Detecting an internal gap means comparing the chain's actual residues
against what *should* be there, and PDBFixer's `findMissingResidues()` gets
that reference sequence only from `SEQRES` records. A PDB with no `SEQRES`
(common for hand-edited or programmatically stripped files) gives PDBFixer
nothing to compare against - confirmed directly: `PDBFixer(filename=...).sequences`
is `[]` on such a file, and `findMissingResidues()` returns `{}` regardless
of how many residues are actually missing.

Rather than silently report "0 gaps repaired" in that situation - easy to
mistake for "there was nothing to fix", especially for someone new to the
tool - **`fixstructure` refuses to write any output at all** when it has no
reference sequence to check against, with an explicit error naming the
problem. This is a hard stop, not a warning: a PDB that *looks* repaired
but might still be silently missing internal residues is worse than an
error message telling you it can't be verified.

The fix is `--pdb-id` instead of `--pdb-file`: downloading directly from
RCSB always carries the official `SEQRES`, so gap detection is reliable by
construction - there is no local-file/no-`SEQRES` case to worry about.

```bash
pypkatool fixstructure --pdb-id 7A3S --outdir results/
```

Requires network access, and can take on the order of a minute for a large
entry (PDBFixer fetches the full deposited structure to derive the
sequence, not just a lightweight SEQRES query).

## Implementation

`fixstructure` runs `pypkatool/data/fixstructure_worker.py` as a subprocess
inside the separate `pdbfixer` conda environment (see
{func}`pypkatool.core._find_pdbfixer_python` and
{func}`pypkatool.core.fix_structure`). The internal/terminal-gap distinction
is implemented by inspecting the `(chain_index, position)` keys PDBFixer's
`findMissingResidues()` returns: a key with `position == 0` is before the
first present residue in its chain (N-terminal gap), and
`position == len(present_residues_in_chain)` is after the last one
(C-terminal gap) - both are filtered out before `addMissingAtoms()` is
called, so only strictly internal gaps (`0 < position < len(...)`) get
rebuilt. `--select-chains` is implemented with PDBFixer's own
`removeChains(chainIds=...)`, applied before gap/atom detection. Dropping
heterogens (the default) uses PDBFixer's own
`removeHeterogens(keepWater=False)`.
