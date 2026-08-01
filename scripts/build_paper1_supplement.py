#!/usr/bin/env python3
"""Build and extract-test the minimal Paper 1 supplementary archive."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
TEMPLATE = REPO / "results" / "report" / "supplement"
FIXED_ZIP_TIME = (2026, 1, 1, 0, 0, 0)
VERTEX_FILENAMES = (
    "vertices_3_6.json",
    "vertices_3_7.json",
    "vertices_3_8.json",
    "vertices_4_8.json",
    "vertices_3_9.json",
    "vertices_4_9.json",
    "vertices_3_10.json",
    "vertices_4_10.json",
    "vertices_5_10.json",
)


FILES = {
    REPO / "LICENSE": Path("LICENSE"),
    TEMPLATE / "README.md": Path("README.md"),
    TEMPLATE / ".python-version": Path(".python-version"),
    TEMPLATE / "pyproject.toml": Path("pyproject.toml"),
    TEMPLATE / "uv.lock": Path("uv.lock"),
    REPO / "results/data/interference_certificates.json": Path(
        "data/interference_certificates.json"
    ),
    REPO / "results/data/constraints/altunbulak_appendix_a_map.json": Path(
        "data/constraint_source_map.json"
    ),
    REPO / "results/data/constraints/ak_thesis_rank9_10_constraints.json": Path(
        "data/source_rows.json"
    ),
    REPO / "src/gpc_census/data/constraints.json": Path("data/constraints.json"),
    REPO / "scripts/verify_states_standalone.py": Path(
        "scripts/verify_states_standalone.py"
    ),
    REPO / "scripts/verify_interference_certificates_standalone.py": Path(
        "scripts/verify_interference_certificates_standalone.py"
    ),
    REPO / "scripts/verify_vertex_exhaustion_standalone.py": Path(
        "scripts/verify_vertex_exhaustion_standalone.py"
    ),
    REPO / "scripts/verify_constraint_source_standalone.py": Path(
        "scripts/verify_constraint_source_standalone.py"
    ),
    REPO / "scripts/verify_paper1_bundle_standalone.py": Path(
        "scripts/verify_paper1_bundle_standalone.py"
    ),
}


def _copy_files(root: Path) -> None:
    for source, relative in FILES.items():
        if not source.is_file():
            raise FileNotFoundError(source)
        destination = root / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
    _write_state_projection(root / "data/states.jsonl")
    for filename in VERTEX_FILENAMES:
        source = REPO / "results/data/vertices" / filename
        if not source.is_file():
            raise FileNotFoundError(source)
        destination = root / "data/vertices" / source.name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
    _write_source_identity(root / "SOURCE.json")


def _write_state_projection(destination: Path) -> None:
    source = REPO / "results/data/states.jsonl"
    fields = (
        "system",
        "index",
        "classified",
        "integer_form",
        "denominator",
        "closed_form",
    )
    rows = [
        {field: record[field] for field in fields}
        for record in (
            json.loads(line)
            for line in source.read_text(encoding="utf-8").splitlines()
            if line.strip()
        )
    ]
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        "".join(json.dumps(row, separators=(",", ":")) + "\n" for row in rows),
        encoding="utf-8",
    )


def _write_source_identity(destination: Path) -> None:
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=REPO,
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()
    dirty = bool(
        subprocess.run(
            ["git", "status", "--porcelain", "--untracked-files=normal"],
            cwd=REPO,
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
    )
    source_map = json.loads(
        (REPO / "results/data/constraints/altunbulak_appendix_a_map.json").read_text(
            encoding="utf-8"
        )
    )
    identity = {
        "schema_version": 1,
        "artifact": "gpc-census-paper1-supplement",
        "artifact_version": "1.0.0",
        "source_repository": "https://github.com/theorlandog/gpc-census",
        "source_commit": commit,
        "source_tree_dirty": dirty,
        # sys.version carries a pre-release suffix on rc interpreters, which
        # the bundle verifier rejects; the triple is what identifies the build.
        "python_version": "{}.{}.{}".format(*sys.version_info[:3]),
        "sympy_version": __import__("sympy").__version__,
        "pdftotext_version": source_map["source"]["pdftotext_version"],
    }
    destination.write_text(
        json.dumps(identity, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def _write_checksums(root: Path) -> None:
    lines = []
    for path in sorted(root.rglob("*")):
        if not path.is_file() or path.name == "SHA256SUMS":
            continue
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        lines.append(f"{digest}  {path.relative_to(root)}")
    (root / "SHA256SUMS").write_text("\n".join(lines) + "\n", encoding="utf-8")


def _run_checks(root: Path) -> None:
    environment = os.environ.copy()
    environment.pop("PYTHONPATH", None)
    environment["PYTHONNOUSERSITE"] = "1"
    commands = [
        [sys.executable, "scripts/verify_paper1_bundle_standalone.py", "."],
        [sys.executable, "scripts/verify_states_standalone.py", "data/states.jsonl"],
        [
            sys.executable,
            "scripts/verify_interference_certificates_standalone.py",
            "data/interference_certificates.json",
        ],
        [
            sys.executable,
            "scripts/verify_vertex_exhaustion_standalone.py",
            "--constraints",
            "data/constraints.json",
            "--vertices",
            "data/vertices",
        ],
    ]
    for command in commands:
        subprocess.run(command, cwd=root, env=environment, check=True)


def _run_locked_environment(root: Path) -> None:
    with tempfile.TemporaryDirectory(prefix="paper1-locked-env-") as directory:
        environment = os.environ.copy()
        environment.pop("PYTHONPATH", None)
        environment["PYTHONNOUSERSITE"] = "1"
        environment["UV_PROJECT_ENVIRONMENT"] = str(Path(directory) / "venv")
        subprocess.run(
            ["uv", "sync", "--locked", "--no-install-project"],
            cwd=root,
            env=environment,
            check=True,
        )
        for arguments in (
            ["scripts/verify_paper1_bundle_standalone.py", "."],
            ["scripts/verify_states_standalone.py", "data/states.jsonl"],
            [
                "scripts/verify_interference_certificates_standalone.py",
                "data/interference_certificates.json",
            ],
            [
                "scripts/verify_vertex_exhaustion_standalone.py",
                "--constraints",
                "data/constraints.json",
                "--vertices",
                "data/vertices",
            ],
        ):
            subprocess.run(
                ["uv", "run", "--no-sync", "python", *arguments],
                cwd=root,
                env=environment,
                check=True,
            )


def _write_zip(root: Path, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(
        output, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9
    ) as archive:
        for path in sorted(root.rglob("*")):
            if not path.is_file():
                continue
            info = zipfile.ZipInfo(str(path.relative_to(root)), FIXED_ZIP_TIME)
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = (0o100644 & 0xFFFF) << 16
            archive.writestr(info, path.read_bytes(), compresslevel=9)


def build(output: Path, test_locked_environment: bool = False) -> None:
    with tempfile.TemporaryDirectory(prefix="paper1-supplement-") as directory:
        root = Path(directory) / "gpc-census-paper1-supplement"
        root.mkdir()
        _copy_files(root)
        _write_checksums(root)
        _run_checks(root)
        if test_locked_environment:
            _run_locked_environment(root)
        _write_zip(root, output)
    with tempfile.TemporaryDirectory(prefix="paper1-extracted-") as directory:
        root = Path(directory) / "gpc-census-paper1-supplement"
        with zipfile.ZipFile(output) as archive:
            archive.extractall(root)
        _run_checks(root)
        if test_locked_environment:
            _run_locked_environment(root)
    digest = hashlib.sha256(output.read_bytes()).hexdigest()
    print(f"PASS: extract-tested Paper 1 supplement written to {output}")
    print(f"SHA-256: {digest}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO / "build/gpc-census-paper1-supplement.zip",
    )
    parser.add_argument(
        "--test-locked-environment",
        action="store_true",
        help="create the documented external uv environment and rerun its checks",
    )
    args = parser.parse_args()
    build(args.output.resolve(), args.test_locked_environment)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
