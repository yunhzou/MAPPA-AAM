# Self-contained notebook validation

Historical record: the notebook `docs/AAM_SIMPLE.ipynb` validated here used the
September 2026 interface and was later removed from the repository. Its last
copy is in revision `2689e5e`.

The executed example was `docs/AAM_SIMPLE.ipynb`.
It contains both endpoints, display helpers, public imports, and saved outputs.
No benchmark archive, notebook helper module, path injection, xTB calculation,
or downloaded molecule is required. Install the package with its notebook extras first.

All 16 code cells ran from a fresh temporary working directory outside the
checkout, with one worker. The run took about two seconds including kernel
startup; this is an example smoke test, not a runtime benchmark.
The JSON round-trip was checked both structurally and by decoding the restored
AAM again. Both selectable candidates were checked in the symmetry-query cells.
The default R2/R3 swap is allowed for candidate 0 and forbidden for candidate 1
within its retained families.

Chrome loaded the saved py3Dmol output and the actual notebook animation iframe
with HTTP(S) blocked. Both rendered without JavaScript errors or network requests.
Every frame of all eight paths across six sweep contexts was visited; final
mapping previews matched saved terminal mappings, and timed playback advanced.
Static drawings and 3D views were visually inspected.

Focused regression tests:

```sh
python -m pytest -q tests/test_search_trajectory.py tests/test_python_workflow.py
```

All 11 passed. Added checks cover in-memory/archive replay agreement, preservation
of the input AAM, selection validation, inclusive event boundaries, metal
thresholds, and target-index correspondence in the event overlay.

For a fresh execution after installation:

```python
import tempfile
import nbformat
from nbclient import NotebookClient

notebook = nbformat.read("docs/AAM_SIMPLE.ipynb", as_version=4)
with tempfile.TemporaryDirectory() as folder:
    NotebookClient(notebook, timeout=120, kernel_name="python3",
                   resources={"metadata": {"path": folder}}).execute()
assert not any(output.output_type == "error"
               for cell in notebook.cells for output in cell.get("outputs", []))
```

Use the Python environment in which the repository was installed as the kernel.
`validation.json` records the notebook hash, package versions, and browser checks.
GitHub displays saved text and PNGs; open and trust the notebook in Jupyter for
interactive 3D views and playback.
