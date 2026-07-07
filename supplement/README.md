# Supplement (source)

LaTeX source for the online supplemental material.

```
supplement/
├── supplementary.tex   ← the document
├── figures/            ← all figures it includes (37 files)
└── references.bib      ← bibliography
```

## Compiling

The document uses the `apa7` class with `biblatex` (biber backend). Compile with
`latexmk` (recommended — it runs the passes and biber automatically):

```bash
cd supplement
latexmk -pdf supplementary.tex
```

or manually:

```bash
pdflatex supplementary
biber supplementary
pdflatex supplementary
pdflatex supplementary   # a second pass populates the table of contents + cross-references
```

Requires a TeX distribution with `apa7`, `biblatex`, `biber`, and the usual
packages (`booktabs`, `tabularx`, `threeparttable`, `hyperref`, `pdflscape`, …).

The compiled `supplementary.pdf` is not bundled (it would go stale); build it from
this source. See ../README.md for the analysis pipeline that produces the figures
and results reported here.
