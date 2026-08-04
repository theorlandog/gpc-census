#!/usr/bin/env python3
"""Run the pinned BDR backend while retaining its prefilter candidates.

This is an optional discovery backend.  It requires the ``moment-cone``
package from the pinned upstream revision and delegates its normal command
line after two deliberate changes:

1. the randomized stabilizer filter is replaced by a pass-through; and
2. ``--candidate-superset`` replaces the default probabilistic inequality
   filters by an empty filter list.

The resulting export is not itself a constraint system.  Feed it to the
local exact trace, Ressayre determinant, Schubert, and redundancy stages.
Even then, a finite failed determinant-evaluation search is unresolved, not
an invalidity proof.  The flag name describes the intended prefilter mode;
this wrapper does not prove that the upstream candidate step is exhaustive.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
import time
from collections.abc import Callable, Iterable
from importlib.metadata import distribution
from pathlib import Path
from typing import Any


EXPECTED_REVISION = "7ab347d9a7837d68d92bde9fac74606913f395ca"


def _pop_single_path_option(args: list[str], option: str) -> tuple[list[str], Path | None]:
    """Remove one wrapper-owned ``OPTION PATH`` pair from upstream arguments."""
    positions = [index for index, value in enumerate(args) if value == option]
    if not positions:
        return args, None
    if len(positions) != 1:
        raise SystemExit(f"{option} may be supplied only once")
    index = positions[0]
    if index + 1 >= len(args) or args[index + 1].startswith("--"):
        raise SystemExit(f"{option} requires a path")
    path = Path(args[index + 1])
    return args[:index] + args[index + 2 :], path


def _write_provenance(path: Path | None, payload: dict[str, object]) -> None:
    """Write an optional deterministic backend sidecar."""
    if path is not None:
        path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def _write_candidate_export(path: Path, inequalities: Iterable[Any]) -> None:
    """Write prefilter homogeneous rows in the safe parser's input format."""
    rows = [tuple(int(value) for value in row.wtau.flattened) for row in inequalities]
    lines = [
        "# Unfiltered inequality candidates captured before probabilistic filters",
        "brut_inequations = [",
        *(f"    {row!r}," for row in rows),
        "]",
        "",
    ]
    path.write_text("\n".join(lines))


def _capturing_candidate_apply(
    original_apply: Callable[[Any, Any], Any], output_path: Path
) -> Callable[[Any, Any], Any]:
    """Wrap the upstream transformer and preserve its prefilter output."""

    def apply(step: Any, tau_dataset: Any) -> Any:
        result = original_apply(step, tau_dataset)
        pending = tuple(result.pending())
        validated = tuple(result.validated())
        _write_candidate_export(output_path, (*pending, *validated))
        return step.TDataset.from_separate(pending=pending, validated=validated)

    return apply


def _installed_revision() -> str | None:
    """Recover a PEP 610 VCS revision or a local checkout revision."""
    direct_url = distribution("moment-cone").read_text("direct_url.json")
    if direct_url is None:
        return None
    record = json.loads(direct_url)
    vcs = record.get("vcs_info", {})
    if vcs.get("commit_id"):
        return str(vcs["commit_id"])
    url = record.get("url", "")
    if not url.startswith("file://"):
        return None
    checkout = Path(url.removeprefix("file://"))
    result = subprocess.run(
        ["git", "-C", str(checkout), "rev-parse", "HEAD"],
        check=False,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip() if result.returncode == 0 else None


def _installed_checkout_dirty() -> bool | None:
    """Report tracked changes for a PEP 610 local checkout, if applicable."""
    direct_url = distribution("moment-cone").read_text("direct_url.json")
    if direct_url is None:
        return None
    url = json.loads(direct_url).get("url", "")
    if not url.startswith("file://"):
        return None
    checkout = Path(url.removeprefix("file://"))
    result = subprocess.run(
        ["git", "-C", str(checkout), "status", "--porcelain", "--untracked-files=no"],
        check=False,
        capture_output=True,
        text=True,
    )
    return bool(result.stdout.strip()) if result.returncode == 0 else None


def _passthrough_stabilizer(self: Any, tau_dataset: Any) -> Any:
    """Preserve every exact upstream tau candidate for local certification."""
    return self.TDataset.from_separate(
        pending=tau_dataset.pending(),
        validated=tau_dataset.validated(),
    )


def _sandbox_process_time(_cls: type[Any]) -> int:
    """Return local CPU time when psutil cannot inspect sandbox children."""
    return time.process_time_ns()


def _ensure_writable_dot_sage() -> str:
    """Select and create a writable Sage state directory before worker startup."""
    configured = os.environ.get("DOT_SAGE")
    path = (
        Path(configured)
        if configured
        else Path(tempfile.gettempdir()) / "gpc-census-sage"
    )
    path.mkdir(parents=True, exist_ok=True)
    os.environ["DOT_SAGE"] = str(path)
    return str(path)


def main() -> None:
    args = sys.argv[1:]
    args, provenance_output = _pop_single_path_option(args, "--provenance-output")
    args, candidate_output = _pop_single_path_option(args, "--candidate-output")
    allow_unpinned = "--allow-unpinned-backend" in args
    allow_dirty = "--allow-dirty-backend" in args
    candidate_superset = "--candidate-superset" in args
    args = [
        arg
        for arg in args
        if arg
        not in {
            "--allow-dirty-backend",
            "--allow-unpinned-backend",
            "--candidate-superset",
        }
    ]

    revision = _installed_revision()
    if revision != EXPECTED_REVISION and not allow_unpinned:
        raise SystemExit(
            "moment-cone revision mismatch: "
            f"expected {EXPECTED_REVISION}, found {revision!r}; "
            "use --allow-unpinned-backend only for an explicitly provisional run"
        )
    checkout_dirty = _installed_checkout_dirty()
    if checkout_dirty and not allow_dirty:
        raise SystemExit(
            "moment-cone local checkout has tracked changes; use a clean pinned "
            "checkout or --allow-dirty-backend only for an explicitly provisional run"
        )

    dot_sage = _ensure_writable_dot_sage()

    from moment_cone.main import moment_cone_from_cmd
    from moment_cone.main_steps import InequalityCandidatesStep, StabilizerConditionStep
    from moment_cone.task import Task

    StabilizerConditionStep.apply = _passthrough_stabilizer
    Task.current_process_time = classmethod(_sandbox_process_time)
    if candidate_output is not None:
        InequalityCandidatesStep.apply = _capturing_candidate_apply(
            InequalityCandidatesStep.apply,
            candidate_output,
        )
    if candidate_superset:
        args.append("--filters")
    sys.argv = [sys.argv[0], *args]
    provenance: dict[str, object] = {
        "backend": "ea-icj/moment_cone",
        "expected_revision": EXPECTED_REVISION,
        "observed_revision": revision,
        "checkout_dirty": checkout_dirty,
        "randomized_stabilizer": "bypassed",
        "task_cpu_accounting": "current_process_only_fallback",
        "dot_sage": dot_sage,
        "inequality_filters": (
            "none_candidate_superset" if candidate_superset else "upstream_cli"
        ),
        "delegated_arguments": args,
        "prefilter_candidate_output": (
            None if candidate_output is None else str(candidate_output)
        ),
        "run_status": "started",
        "proof_status": "discovery_only_requires_local_exact_certification",
    }
    _write_provenance(provenance_output, provenance)
    print(json.dumps(provenance, sort_keys=True), file=sys.stderr)
    moment_cone_from_cmd()
    provenance["run_status"] = "completed"
    _write_provenance(provenance_output, provenance)


if __name__ == "__main__":
    main()
