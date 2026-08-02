# Observable Manifolds Program Packet

This packet converts the current research direction into four coordinated
tracks:

1. theorem and novelty audit;
2. referee-grade paper plan;
3. polished software architecture;
4. staged implementation and release roadmap.

It is deliberately conservative. Standard geometric-quantum-mechanics facts
are separated from the new support-combinatorial and reduced-observable
claims, and benchmark claims are restricted to their tested algorithm classes.

The packet is a plan, not evidence. Every status in `CLAIMS_LEDGER.csv` names
the artifact under `results/data/` that backs it, and the ledger is the file
to update when a status changes. Two entries already differ from what the
tracks were originally drafted against, because the work landed first:
determinacy closure on the unresolved 28 has run and resolved none, and the
benchmark suite has been rerun on a second machine, where every runtime factor
moved while the method ordering, variable counts and accuracies held. The
prior-art track builds on `docs/prior_art_roadmap.md` rather than starting a
second survey, and the adversarial review in
`docs/adversarial_observable_audit.md` has already run, refuting two readings
of T1 that the tracks were drafted against.
