# AI Contributor Guide for RNAvigate

This file gives AI agents context and behavioral expectations to contribute effectively
to RNAvigate.

---

## 1. Always read the following files for context:

1. `README.md`: one-paragraph project description and audience
3. `TODO.md`: prioritized backlog with explainations and reasons for each item
2. `ARCHITECTURE.md`: map of modules, classes, patterns, and dependencies
4. `docs/source/resources/changelog.rst`: release history and breaking changes

---

## 2. Project Identity

**Who is it for?** RNA structure biologists, many of whom are non-programmers working
in Jupyter Notebooks.

**What is it?** A Jupyter-compatible Python library for loading, aligning, analyzing,
and visualizing RNA per-nucleotide and inter-nucleotide data alongside secondary and
tertiary structure models and annotations.

**Design priorities:** usability is more important than performance.

**Primary entry points:** `rnavigate.Sample` and `rnavigate.plot_*()`.

---

## 3. Python Standards

Follow these conventions on every file you write or edit.

### Style

- Comply with PEP 8.
- Line length limit: 88 characters.
- Use `ruff`-compatible formatting and import order (stdlib → third-party → local).
- Prefer `pathlib.Path` over `os.path`.
- Use f-strings over `.format()` or `%` formatting.

### Naming

- Classes: `UpperCamelCase`
- Functions and methods: `snake_case`
- Constants: `UPPER_SNAKE_CASE`
- Private helpers: `_leading_underscore`
- Never use single-letter variable names outside of short loop indices or
  mathematical notation where the variable name is standard (e.g. `i`, `j`
  for nucleotide positions is fine here).

### Type hints

- Add PEP 484 type hints to all new public functions and methods.
- For existing functions you modify substantially, add hints at the same time.
- Use `from __future__ import annotations` at the top of the file to allow
  forward references without quotes.

### Docstrings

- Use **NumPy-style** docstrings on public classes, methods, and functions.
- One-line summary on the first line, blank line, then extended description.
- Required sections: `Parameters`, `Returns`. Add `Raises`, `Notes`,`Examples`, and
  `References` where appropriate.

### Error handling

- Raise specific exceptions (`ValueError`, `TypeError`, `FileNotFoundError`)
  with informative messages that name the offending value and explain what was
  expected.
- Never use bare `except:`. Catch the narrowest exception that makes sense.
- Avoid using `print()` for warnings in library code; use `warnings.warn()`
  with an appropriate category (`UserWarning`, `DeprecationWarning`, etc.).

---

## 4. Architecture Rules

These rules encode design decisions that must not be reversed without explicit
discussion. See `ARCHITECTURE.md` for the rationale behind each.

1. **`Sample` does not import from `plots/`**. The data layer and the
   visualization layer must remain decoupled.
2. **`data/` classes do not import from `analysis/` or `plots/`**. Data flows
   downward: `analysis → Sample → data`, `plots → data`. Never upward.
3. **New data types register via `data_loading.py`**. All keyword-to-class
   mappings live in `data_keyword_defaults`. Do not bypass this by
   constructing data objects directly in `Sample.__init__`.
4. **`get_aligned_data(alignment)` must be implemented** on every new `Data`
   subclass. This is what makes multi-sample plotting work.
5. **Do not add new circular imports**. The existing `data.py` ↔
   `alignments.py` circularity is a known issue to be fixed, not a pattern
   to copy.
6. **`plots/functions/`** contains only stateless, axes-level drawing
   primitives. No data loading, no alignment, no `Sample` access.
7. **Analysis entry points should subclass `Sample`** so their results are
   directly usable with all `plot_*()` functions.

---

## 5. Testing Expectations

- Every new public function or class must be accompanied by at least one
  `pytest` test (unit or integration).
- Use fixtures from `rnavigate.examples` where real data is needed — do not
  create synthetic fixtures that don't resemble real file formats.
- Tests live in `tests/` mirroring the package structure:
  `tests/data/`, `tests/plots/`, `tests/analysis/`, etc.
- Smoke tests for `plot_*()` functions should call the function and assert
  the return value is a `plots.Plot` instance. They do not need to validate
  pixel output.
- Run the full test suite with `pytest` before reporting any task complete.
  If tests do not yet exist, note this explicitly.

---

## 6. Workflow

### Before starting any task

1. Read the four context files listed in §1.
2. Check `TODO.md` to understand where the task fits in the broader plan.
3. If the task involves a structural change (new module, moved class, new
   dependency), update `ARCHITECTURE.md` first and confirm with the maintainer
   before writing code.

### While working

- Make the smallest change that accomplishes the goal. Avoid unrelated edits.
- If you discover a bug or issue unrelated to your current task, add it to
  `TODO.md` under the appropriate section rather than fixing it silently.
- Prefer editing existing files over creating new ones.
- When creating a new file, add it to the directory map in `ARCHITECTURE.md`.

### Before reporting a task complete

- Confirm the code runs without import errors.
- Confirm `pytest` passes (or note explicitly if tests don't yet exist).
- Confirm `ARCHITECTURE.md` and `TODO.md` reflect any changes made.
- Provide a one-paragraph summary of what changed and why.

---

## 7. What to Avoid

- **Do not silently refactor** code that isn't part of the requested task.
- **Do not add dependencies** without adding them to both `pyproject.toml`
  and `environment.yml`, and noting them in `ARCHITECTURE.md` §9.
- **Do not remove public API names** (`Sample`, `plot_*`, exported data
  classes) without explicit maintainer approval — this is a breaking change.
- **Do not use `print()` for diagnostic output** in library code.
- **Do not leave `TODO` or `FIXME` comments in code** without adding a
  corresponding entry to `TODO.md`.
- **Do not guess at file formats.** If a new file format is needed, ask for
  a sample file or a format specification before implementing the parser.
