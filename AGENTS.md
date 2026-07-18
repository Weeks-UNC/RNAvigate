# AI Contributor Guide for RNAvigate

---

## 1. Always read for context

1. `README.md`: project description and audience
2. `TODO.md`: prioritized backlog with explanations
3. `ARCHITECTURE.md`: modules, classes, patterns, dependencies
4. `docs/source/resources/changelog.rst`: release history and breaking changes

---

## 2. Project Identity

**Audience:** RNA structure biologists, many non-programmers using Jupyter Notebooks.

**What it is:** A Jupyter-compatible Python library for loading, aligning, analyzing, and visualizing RNA per-nucleotide and inter-nucleotide data alongside secondary/tertiary structure models and annotations.

**Design priority:** Usability over performance.

**Primary entry points:** `rnavigate.Sample` and `rnavigate.plot_*()`.

---

## 3. Python Standards

### Style

- PEP 8. Line limit: 88 characters.
- `ruff`-compatible formatting; import order: stdlib → third-party → local.
- `pathlib.Path` over `os.path`. F-strings over `.format()` or `%`.

### Naming

- Classes: `UpperCamelCase`; functions/methods: `snake_case`; constants: `UPPER_SNAKE_CASE`; private helpers: `_leading_underscore`.
- No single-letter names except standard loop indices (`i`, `j`) or math notation.

### Type hints

- PEP 484 hints on all new public functions/methods.
- Add hints to existing functions when making substantial edits.
- `from __future__ import annotations` at file top for forward references.

### Docstrings

- NumPy-style on all public classes, methods, functions.
- One-line summary, blank line, extended description.
- Required sections: `Parameters`, `Returns`. Add `Raises`, `Notes`, `Examples`, `References` as needed.

### Error handling

- Raise specific exceptions (`ValueError`, `TypeError`, `FileNotFoundError`) with messages naming the offending value and expected input.
- No bare `except:`. Catch the narrowest applicable exception.
- Use `warnings.warn()` with an appropriate category; never `print()` for warnings.

---

## 4. Architecture Rules

Do not reverse these without explicit discussion. See `ARCHITECTURE.md` for rationale.

1. **`Sample` does not import from `plots/`**: data and visualization layers stay decoupled.
2. **`data/` classes do not import from `analysis/` or `plots/`**: data flows downward only: `analysis → Sample → data`, `plots → data`.
3. **New data types register via `data_loading.py`**: keyword-to-class mappings live in `data_keyword_defaults`; don't bypass this in `Sample.__init__`.
4. **`get_aligned_data(alignment)` must be implemented** on every new `Data` subclass (required for multi-sample plotting).
5. **No new circular imports**: the `data.py` ↔ `alignments.py` circularity is a known issue to fix, not a pattern to copy.
6. **`plots/functions/`**: stateless axes-level drawing primitives only; no data loading, alignment, or `Sample` access.
7. **Analysis entry points should subclass `Sample`** so results work with all `plot_*()` functions.

---

## 5. Testing

- Every new public function or class needs at least one `pytest` test (unit or integration).
- Use fixtures from `rnavigate.examples` for real data; don't create synthetic fixtures that don't resemble real file formats.
- Tests mirror package structure: `tests/data/`, `tests/plots/`, `tests/analysis/`, etc.
- `plot_*()` smoke tests: call the function, assert return value is a `plots.Plot` instance. No pixel validation needed.
- Run `pytest` before reporting a task complete. Note explicitly if tests don't yet exist.

---

## 6. Workflow

### Before starting

1. Read the four context files in §1.
2. Check `TODO.md` to understand where the task fits.
3. For structural changes (new module, moved class, new dependency), update `ARCHITECTURE.md` first and confirm with the maintainer before writing code.

### While working

- Make the smallest change that accomplishes the goal.
- Unrelated bugs discovered → add to `TODO.md`; don't fix silently.
- Prefer editing existing files over creating new ones.
- New files → add to the directory map in `ARCHITECTURE.md`.

### Before reporting complete

- Confirm no import errors.
- Confirm `pytest` passes (or note if tests don't yet exist).
- Confirm `ARCHITECTURE.md` and `TODO.md` reflect all changes.
- Provide a one-paragraph summary of what changed and why.

---

## 7. What to Avoid

- **Silent refactoring** of code outside the requested task.
- **Adding dependencies** without updating `pyproject.toml`, `environment.yml`, and `ARCHITECTURE.md` §9.
- **Removing public API names** (`Sample`, `plot_*`, exported data classes) without explicit maintainer approval.
- **`print()` for diagnostic output** in library code.
- **`TODO`/`FIXME` comments in code** without a corresponding `TODO.md` entry.
- **Guessing file formats**: ask for a sample file or spec before implementing a parser.
