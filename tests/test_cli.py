import json

import pytest

from gpc_census import __version__
from gpc_census.cli import main


def test_export_classification_is_machine_readable(capsys):
    assert main(["export", "-n", "3", "-d", "9", "--kind", "classification"]) == 0
    rows = json.loads(capsys.readouterr().out)
    assert len(rows) == 58
    assert all({"index", "integer_form", "verdict"} <= set(r) for r in rows)
    assert {r["verdict"] for r in rows} <= {"DESIGN-INT", "DESIGN-REAL", "INTERFERENCE"}


def test_export_vertices_kind(capsys):
    assert main(["export", "-n", "3", "-d", "9", "--kind", "vertices"]) == 0
    rows = json.loads(capsys.readouterr().out)
    assert len(rows) == 58
    assert all("spectrum" in r for r in rows)


def test_export_all_has_sections(capsys):
    assert main(["export", "-n", "4", "-d", "9", "--kind", "all"]) == 0
    obj = json.loads(capsys.readouterr().out)
    assert set(obj) == {"system", "constraints", "vertices", "classification", "states"}
    assert obj["system"] == {"n": 4, "d": 9}


def test_states_precomputed_is_json_list(capsys):
    assert main(["states", "-n", "3", "-d", "9"]) == 0
    recs = json.loads(capsys.readouterr().out)
    assert isinstance(recs, list)


def test_vertices_output(capsys):
    assert main(["vertices", "-d", "4", "-n", "1"]) == 0
    out = capsys.readouterr().out.strip().splitlines()
    assert out == ["1 0 0 0", "0 1 0 0", "0 0 1 0", "0 0 0 1"]


def test_version(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["--version"])
    assert exc.value.code == 0
    assert __version__ in capsys.readouterr().out


def test_invalid_arguments_exit_nonzero(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["vertices", "-d", "2", "-n", "3"])
    assert exc.value.code != 0
