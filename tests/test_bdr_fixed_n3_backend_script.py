from __future__ import annotations

import importlib.util
import json
import time
from pathlib import Path


ROOT = Path(__file__).parents[1]
SCRIPT = ROOT / "scripts" / "bdr_fixed_n3_backend.py"


def _module():
    spec = importlib.util.spec_from_file_location("bdr_fixed_n3_backend", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_installed_revision_reads_pep610_vcs_record(monkeypatch):
    module = _module()

    class FakeDistribution:
        def read_text(self, filename):
            assert filename == "direct_url.json"
            return json.dumps(
                {
                    "url": "https://example.invalid/moment_cone.git",
                    "vcs_info": {"commit_id": module.EXPECTED_REVISION},
                }
            )

    monkeypatch.setattr(module, "distribution", lambda _name: FakeDistribution())
    assert module._installed_revision() == module.EXPECTED_REVISION


def test_installed_checkout_dirty_reads_local_git_state(monkeypatch, tmp_path: Path):
    module = _module()

    class FakeDistribution:
        def read_text(self, filename):
            assert filename == "direct_url.json"
            return json.dumps({"url": tmp_path.as_uri()})

    class Result:
        returncode = 0
        stdout = " M src/example.py\n"

    monkeypatch.setattr(module, "distribution", lambda _name: FakeDistribution())
    monkeypatch.setattr(module.subprocess, "run", lambda *args, **kwargs: Result())
    assert module._installed_checkout_dirty() is True


def test_passthrough_preserves_pending_and_validated():
    module = _module()

    class Dataset:
        def pending(self):
            return (1, 2)

        def validated(self):
            return (3,)

    class Output:
        @staticmethod
        def from_separate(**kwargs):
            return kwargs

    class Step:
        TDataset = Output

    assert module._passthrough_stabilizer(Step(), Dataset()) == {
        "pending": (1, 2),
        "validated": (3,),
    }


def test_sandbox_process_time_is_a_monotone_integer():
    module = _module()
    before = module._sandbox_process_time(object)
    time.process_time_ns()
    after = module._sandbox_process_time(object)
    assert type(before) is int
    assert after >= before


def test_wrapper_selects_a_writable_sage_state_directory(
    monkeypatch, tmp_path: Path
):
    module = _module()
    monkeypatch.delenv("DOT_SAGE", raising=False)
    monkeypatch.setattr(module.tempfile, "gettempdir", lambda: str(tmp_path))

    selected = module._ensure_writable_dot_sage()

    expected = tmp_path / "gpc-census-sage"
    assert selected == str(expected)
    assert expected.is_dir()
    assert module.os.environ["DOT_SAGE"] == str(expected)


def test_wrapper_path_option_is_removed_and_sidecar_is_stable(tmp_path: Path):
    module = _module()
    sidecar = tmp_path / "run.json"
    args, selected = module._pop_single_path_option(
        ["Fermion", "13", "--provenance-output", str(sidecar), "--quiet"],
        "--provenance-output",
    )
    assert args == ["Fermion", "13", "--quiet"]
    assert selected == sidecar

    payload = {"run_status": "completed", "revision": module.EXPECTED_REVISION}
    module._write_provenance(sidecar, payload)
    assert json.loads(sidecar.read_text()) == payload


def test_candidate_capture_preserves_dataset_and_writes_literal(tmp_path: Path):
    module = _module()
    output = tmp_path / "candidates.py"

    class Tau:
        def __init__(self, values):
            self.flattened = values

    class Inequality:
        def __init__(self, values):
            self.wtau = Tau(values)

    class Dataset:
        def __init__(self, pending, validated):
            self._pending = pending
            self._validated = validated

        def pending(self):
            return iter(self._pending)

        def validated(self):
            return iter(self._validated)

    class Output:
        @staticmethod
        def from_separate(**kwargs):
            return kwargs

    class Step:
        TDataset = Output

    source = Dataset([Inequality((1, 0, -1))], [Inequality((0, 0, -1))])
    wrapped = module._capturing_candidate_apply(lambda _step, _data: source, output)
    preserved = wrapped(Step(), object())

    assert tuple(row.wtau.flattened for row in preserved["pending"]) == ((1, 0, -1),)
    assert tuple(row.wtau.flattened for row in preserved["validated"]) == ((0, 0, -1),)
    assert "brut_inequations" in output.read_text()
    assert "(1, 0, -1)" in output.read_text()
