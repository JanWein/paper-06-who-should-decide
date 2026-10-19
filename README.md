# Who Should Decide? Information Access, Strategic Proposals, and Scarce Review

**Dr. Jan-Hendrik Weinert** · Working paper · Editorial edition F1, 19 October 2026

[Read the paper](paper.pdf) · [LaTeX source](manuscript.tex) · [Editable Markdown](manuscript.md) · [Reproduction instructions](REPRODUCIBILITY.md)

## Abstract

Submission-contingent review can improve access to scarce audits while creating adverse filing equilibria. We solve a population proposal game and a two-manager extension with one non-storable audit slot. The finite model characterizes every filing equilibrium under arbitrary randomized rules and derives exact strict-deterrence target frontiers for preassignment and contingent redeployment. Higher-risk-first recourse can support complete deterrence alongside a both-file equilibrium and an interior mixed equilibrium. At an uninformative witness, every recourse rule supporting strict complete deterrence also admits strict both-filing. Target resource costs have a finite-candidate solution, but lower target costs do not establish worst-equilibrium dominance. Complete proofs and synthetic computations distinguish attainable targets from robust implementation within the stated institution.

## Contents

- `paper.pdf`: complete paper, including proofs, bibliography and appendices.
- `manuscript.tex`: authoritative mathematical typesetting source; `manuscript.md`: editable alternative.
- Calculation scripts, recorded synthetic outputs, and the paper-specific source ledgers and prospective research protocols.
- `reproduce.py`, `requirements.txt`: one entry point and dependencies for the supplied analysis.
- `build_paper.sh`: rebuilds the current PDF from the LaTeX source without replacing the supplied PDF.
- `CITATION.cff`: citation metadata; `SHA256SUMS.txt`: published payload checksums.

## Reproduce the calculations

Use Python 3.12. From the repository root, create a virtual environment and run:

```sh
python3 -m venv .venv
. .venv/bin/activate
python3 -m pip install -r requirements.txt
python3 reproduce.py
```

The analysis uses synthetic inputs or exact model primitives supplied with the code. It requires no private organizational data, API credentials or data download. Installing dependencies may require access to a package index. Full analytical proofs remain in the paper; finite computations supplement those proofs.

## Rebuild the paper

With TeX Live and the packages declared in `manuscript.tex` installed:

```sh
bash build_paper.sh
```

The result is `build/manuscript.pdf`. Compile from the repository root so figure paths resolve. The supplied `paper.pdf` is the visually checked publication edition. PDF byte hashes can vary between rebuilds because of timestamps and toolchain metadata. Do not regenerate this edition's LaTeX through Pandoc when a matching layout is required.

## Research status

This is a theoretical working paper with synthetic computational illustrations. Any empirical or field protocol included here is prospective, and no reported synthetic output is an observed firm result. Publication of this repository does not imply journal acceptance or empirical validation. The paper states its model-specific assumptions and limitations.

## Citation

Weinert, Jan-Hendrik (2026). *Who Should Decide? Information Access, Strategic Proposals, and Scarce Review*. Working paper, editorial edition F1.
Use the repository's citation metadata when citing or reusing the work.
