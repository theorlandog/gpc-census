import pytest

from gpc_census import __version__
from gpc_census.cli import main


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
