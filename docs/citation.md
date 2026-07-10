# Citation

## References

- PyPKA: Reis, P. B. P. S.; Vila-Viçosa, D.; Rocchia, W.; Machuqueiro, M.
  *J. Chem. Inf. Model.* 2020, 60 (10), 4442-4448.
  [DOI: 10.1021/acs.jcim.0c00718](https://doi.org/10.1021/acs.jcim.0c00718)
- pKAI / pKAI+: Reis, P. B. P. S.; Bertolini, M.; Montanari, F.; Rocchia, W.;
  Machuqueiro, M.; Clevert, D.-A. *J. Chem. Theory Comput.* 2022, 18 (8), 5068-5078.
  [DOI: 10.1021/acs.jctc.2c00308](https://doi.org/10.1021/acs.jctc.2c00308)
- CHARMM36: Best, R. B.; Zhu, X.; Shim, J.; Lopes, P. E. M.; Mittal, J.;
  Feig, M.; MacKerell, A. D. Jr. *J. Chem. Theory Comput.* 2012, 8 (9),
  3257-3273. [DOI: 10.1021/ct300400x](https://doi.org/10.1021/ct300400x) -
  source of `top_all36_prot.rtf`'s RESI/PRES blocks for all protonation
  states used here.
- CHARMM36m: Huang, J.; Rauscher, S.; Nawrocki, G.; Ran, T.; Feig, M.;
  de Groot, B. L.; Grubmüller, H.; MacKerell, A. D. Jr. *Nat. Methods* 2017,
  14 (1), 71-73. [DOI: 10.1038/nmeth.4067](https://doi.org/10.1038/nmeth.4067) -
  the specific, more current CHARMM36 revision PyPKA's tautomer library
  (`CHARMM36m/sts/`) is built on.

## How to cite this repository

If you use `pypkatool` in published work, please cite it alongside the
methods it wraps (PyPKA, pKAI+, CHARMM36/CHARMM36m - see "References"
above), since those are what actually compute the pKa values and define
the protonation-state topology.

Citation metadata for `pypkatool` itself is kept in
[`CITATION.cff`](https://github.com/groponp/PyPkaTool/blob/main/CITATION.cff)
(GitHub reads this automatically and adds a "Cite this repository" button -
APA/BibTeX export - to the repo page). Manually:

```bibtex
@software{ropon_palacios_pypkatool,
  author  = {Ropón-Palacios, G.},
  title   = {pypkatool},
  version = {1.0.0},
  date    = {2026-07-08},
  url     = {https://github.com/groponp/PyPkaTool}
}
```

## Author

**Ropón-Palacios G.**
Department of Physics, UNESP.
[georcki.ropon@unesp.br](mailto:georcki.ropon@unesp.br)

## Disclaimer

This software is provided "as is", without warranty of any kind, express or
implied, including but not limited to the warranties of merchantability,
fitness for a particular purpose, and noninfringement - see the full text in
the [LICENSE](https://github.com/groponp/PyPkaTool/blob/main/LICENSE) file.
pKa predictions and CHARMM protonation-state assignments produced by this
pipeline are computational estimates and must be checked against
experimental data and domain expertise before use in downstream modeling;
the authors assume no liability for outcomes derived from its use.

## License

MIT - see [LICENSE](https://github.com/groponp/PyPkaTool/blob/main/LICENSE).
