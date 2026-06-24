# Contributing to phytree

Thanks for your interest. phytree is a pure-Python phylogenetics +
visualization library; contributions that keep it dependency-light and well
tested are very welcome.

## Setup

```bash
pip install -e .[dev]      # installs pytest
pytest -q                  # run the test suite
```

## Architecture (where things live)

- `phytree/core/`        — `Tree`/`Node` model and I/O
- `phytree/layout/`      — topology → display coordinates (subclass `Layout`)
- `phytree/infer/`       — alignment, distances, ML, parsimony, bootstrap
- `phytree/comparative/` — ancestral states, stochastic mapping
- `phytree/plot/`        — TreeFigure builder, elements, matplotlib/plotly backends
- `phytree/scene.py`     — backend-agnostic drawing primitives

The key invariant: **layout computes coordinates and emits scene primitives;
backends only render them.** A new element subclasses `plot.figure._Element`
and appends primitives; a new layout subclasses `layout.base.Layout`.

## Guidelines

- Add a test in `tests/` for any new feature or bug fix.
- Match the surrounding style; keep heavy deps (scipy/biopython/plotly)
  import-local where reasonable.
- Run `pytest -q` before opening a PR; CI must stay green.
- For phylogenetic methods, prefer adding a correctness check (e.g. recovering
  a known topology, or agreement with an independent reference as in
  `validation/`).
