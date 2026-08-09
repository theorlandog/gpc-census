# Continuation handoff: the five-track campaign

Written for whoever picks this up next, human or agent. Read this before
running anything. It records where the campaign actually stands, what is
verified (almost nothing), and the operational traps that cost three runs.

## Current state in one paragraph

Five research tracks were specified, refreshed against HEAD, committed, and
dispatched as parallel worktree-isolated agents. They were launched three
times. All three launches were stopped by external causes, never by a defect
in the work: once by a user interrupt, twice by an account weekly usage limit.
Each launch got further than the last. Real findings exist and are committed,
but NONE of them has passed `make lint` or `make test`, because no run ever
reached the gates. The next person's first job is auditing, not generating.

## Repository coordinates

- Designated working branch: `claude/new-session-e66aw9`.
- Campaign base commit: `e80eae6`, "Land the five parallel agent tracks,
  refreshed against Theorem S". This is on origin.
- `main` at the time of writing: `476bf0b`, "Prove Conjecture S: the pattern
  set determines the values (#127)".
- Track specs live in `docs/agent_tracks/`. `COORDINATION.md` there holds the
  shared handoff protocol every track agent must satisfy.

## The five branches

All five branch from `e80eae6`. They were LOCAL ONLY at the time of writing
and were never pushed, because the session's instructions restrict pushes to
the designated branch. They were exported as `git format-patch` text so they
could survive outside the container. If you are reading this from a fresh
checkout, the branches may no longer exist anywhere except that export.

| Branch | Tip | Commits | State |
| --- | --- | ---: | --- |
| `claude/track-qvalidate-fermionic` | `1c46de3` | 4 | Phase 2 result committed, Phase 3 in flight |
| `claude/track-qec-magic` | `5b5dc64` | 5 | Prior-art finding committed, benchmark in flight |
| `claude/track-entanglement-marginal-spin` | `c8a5b1f` | 4 | Phase 1 reproduction committed, benchmark in flight |
| `claude/track-anyon-newton-scout` | `627a6da` | 4 | Artifact verified, deliverable document in flight |
| `claude/track-generator-frontier` | `b466f52` | 1 | One script, never run to completion |

Restore from the export with, per track:

```sh
git checkout -b claude/track-qec-magic e80eae6
git am < track-qec-magic.txt
```

## What was claimed, and what that is worth

Every item below is the authoring agent's own claim, taken from its commit
subject. Treat the whole table as unaudited. Two of these would change the
campaign's direction if they hold, so audit those first.

| Commit | Claim | Priority |
| --- | --- | --- |
| `c90ede9` | The QEC/magic track's core idea is PRIOR ART | Audit first. If true, narrow or stop that track. |
| `12787d1` | The 28 uncertified carriers are two-body non-injective | Audit first. If true, it is a real result and the basis of Phase 3. |
| `3874eda` | Exact reproduction of the published qubit entanglement polytopes | Audit second. It is the gate licensing the rest of that track. |
| `87a9e45` | Fixes a sign bug that made the certified interval infeasible on real codes | Audit with the QEC track. |
| `533c106` | Scout artifact verified against its generator, guards added | Lowest risk, that track was closest to done. |

Unbanked but reported in the final moments before a stop, so recorded here as
leads and nothing more: the fermionic Phase 3 "gives a decisive negative", and
the entanglement track found "the threshold is exactly GHZ fidelity 4/5".
Neither is committed. Do not cite either without deriving it yourself.

Commits whose subject begins "Salvage" are uncommitted work that was rescued
when a run died. They are preserved so effort was not lost, not because they
are correct. The QEC salvage had a live traceback in `coset_weight_enumerator`
at one point, so that module is known to have been broken at least once.

## The Theorem S correction, do not lose this

The original generator track was written before `476bf0b` landed and asked its
agent to prove that level spread is bounded, or to find an unbounded-spread
counterexample. Both are settled and re-attacking them is pure waste:

- Theorem S is PROVED and guarded by `tests/test_profiles.py`. The pattern
  differences span `V = 1^perp n v^perp`, so `V^perp = span{1, v}`, the pattern
  set determines the values, spread is bounded at every rank, and candidate
  coverage is DECIDABLE rather than open.
- The explicit bound is `spread <= 2 sqrt(r-1) (N sqrt 2)^(r-2)`.
- A counterexample family of unbounded spread at fixed `d` is therefore ruled
  out. Do not go looking for one.

The surviving problem is quantitative: the bound is about `2.8e6` at
`N = 3, r = 11` against a realized maximum of 20, so the open question is
whether the true maximum spread at rank `d` is polynomial in `r`. The named
route is to make the multiplicity-budget argument `sum_a min(m_a, N + 1)
>= f(W)` quantitative. `docs/agent_tracks/AGENT_GENERATOR_FRONTIER.md` Phase 2
already carries this correction and the supporting evidence; the version in
any older copy of the bundle does not.

## Operational traps that actually cost runs

1. **Worktree isolation refuses compound shell commands.** A guard rejects any
   command it cannot verify stays inside the worktree, which killed the first
   generator run outright. Use plain, separate, single commands. Use
   `git -C <path>` instead of `cd`. Never chain with `&&`, `||` or `;`.
2. **Five parallel Opus agents exhaust a weekly budget within minutes.** This
   happened twice in a row, identically. Run tracks serially, or at most two
   at a time. Fan-out is what killed this campaign, not the difficulty of the
   work.
3. **Commit increments as you go.** The first two runs banked nothing and lost
   everything to scratch. The third run was told to commit incrementally and
   is the only reason any result above exists.
4. **Write artifacts to disk the moment they are computed.** The fermionic
   agent completed Phase 2 in run one and lost the artifact because it was
   still in memory when the stop landed.
5. **Sending a message interrupts in-flight subagents.** Asking for a status
   update is what ended run one.

## Recommended order of work

1. Audit `c90ede9` and `12787d1`. They are cheap to check relative to their
   consequences, and they decide whether two tracks continue at all.
2. Run `make lint` and `make test` against each restored branch and record
   what actually passes. No track should advance a phase until its existing
   work is green.
3. Resume the fermionic and generator tracks, serially, one at a time. These
   are the two with genuine in-repo grounding.
4. Treat the QEC track as provisionally stopped pending the prior-art audit.
5. Finish the anyon/Newton scout, which is closest to a deliverable and is
   explicitly meant to be cheap. A well-argued "no" is its success condition.

## House rules that bite

- Never use em dashes in docs, comments, or commit messages.
- Every evidence tag in `docs/RESEARCH.md` must point at an in-repo certifying
  artifact. A claim with no artifact cannot carry a tag.
- Never ship a result that has not passed the structural invariants in
  `gpc_census.validate`, and never change state-solving code without the v_B
  preflight.
- Any fiber-dimension or rigidity claim must name fixed-spectrum vs fixed-rho
  and reduced vs unreduced.
- If you touch `results/data/`, run `make checksums` then `make checksums-check`
  as separate commands. Never hand-edit `results/data/SHA256SUMS`.
- New dependencies must be at least 14 days old, and are discouraged.

Read `AGENTS.md` and `docs/RESEARCH.md` in full before touching the science.
