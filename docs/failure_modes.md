# Failure modes

If PyPKA's preprocessing step fails - most commonly because of a
structurally incomplete PDB - `pypkatool run` reports the underlying error
together with a suggestion to repair the structure first (see
{doc}`fixstructure`), rather than a raw Python traceback.

A gap that is silently *ignored* rather than causing a failure is also
possible: if an entire residue is absent from the input and `fixstructure`
was not run first, PyPKA simply omits that position - it will not appear
anywhere in the output tables, and no warning is issued. Compare the number
of sites reported against the residue count you expect for the protein if
this is a concern.

## What actually raises inside PyPKA

The real crash points are external subprocess tools PyPKA shells out to
during structure cleaning (`pdb2pqr`, `addHtaut`), not a simple Python-level
"missing atom" check - `pdbmender` is more tolerant of incomplete input than
that name might suggest (a single missing terminal atom, or even an entire
missing internal residue, is either repaired automatically or silently
skipped, not a fatal error by itself). When those external tools do fail,
they write their error to stderr, and `pdbmender` re-raises it as a bare
`Exception` with no further explanation (`pdbmender/utils.py:107,157,319`).

`pypkatool` wraps the `Titration()` call in {func}`pypkatool.core.run_pypka`
in a `try`/`except` that turns any such failure into an actionable message
naming the likely cause (an incomplete structure) and the fix
(`fixstructure`), instead of letting the raw exception propagate.
