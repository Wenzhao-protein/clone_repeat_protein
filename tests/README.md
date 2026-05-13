# tests

Smoke and regression tests for the HURDLER toolkit.

These tests assume the HURDLER pipeline has been run at least once and
that the canonical outputs exist under `../output/` (see
[`../docs/workflows/hurdler_site_combinations.md`](../docs/workflows/hurdler_site_combinations.md)).

| Test | Scope |
|------|-------|
| `test_hurdler_system.py` | End-to-end smoke test of the HURDLER pipeline and query layer. |
| `test_hurdler_success_rate.py` | Sanity-checks the success-rate analysis output. |

Run from the repo root:

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
```

Ad-hoc / historical tests are preserved under
[`../archive/tests/`](../archive/tests/) and are **not** part of the
maintained test suite.
