# AK2008 extracted data — census validation corpus
Source: Altunbulak & Klyachko, arXiv:0802.0918v1 (= CMP 282, 287 (2008)), text-extracted 2026-07-16.
Vertices given as projective ratios exactly as printed (normalize to trace N).
CAUTION: extremal-state surd coefficients in Tables 5-6 suffered line-break mangling in
extraction; vertex ratios are clean, but re-verify any state coefficients against the PDF
before use as ground truth. Ranks 9-10 are NOT in the paper text (zip / Springer ESM only).

## (3,6) inequalities [Table 1]
l1+l6=1, l2+l5=1, l3+l4=1 (Borland-Dennis equalities), l4 <= l5+l6

## (3,7) inequalities [Table 2]
l2+l3+l4+l5 <= 2
l1+l3+l4+l6 <= 2
l1+l2+l4+l7 <= 2
l1+l2+l5+l6 <= 2

## (3,7) vertices [Table 6, first ten rows]
(1:1:1:0:0:0:0)
(2:1:1:1:1:0:0)
(1:1:1:1:1:1:0)
(3:3:2:2:1:1:0)
(1:1:1:1:1:1:1)
(2:2:1:1:1:1:1)
(2:2:2:2:2:1:1)
(3:1:1:1:1:1:1)
(3:3:3:3:1:1:1)
(5:5:3:3:3:1:1)
(rank-7 states: [123]; [123]+[145]; [123]+[145]+[246]+[356]; sqrt2[123]+[145]+[246]; ...)

## (4,8) inequalities [Table 3]
l1 <= 1
l5-l6-l7-l8 <= 0
l1-l2-l7-l8 <= 0
l1-l3-l6-l8 <= 0
l1-l4-l6-l7 <= 0
l1-l4-l5-l8 <= 0
l3-l4-l7-l8 <= 0
l2-l4-l6-l8 <= 0
l2+l3+l5-l8 <= 2
l1+l3+l6-l8 <= 2
l1+l2+l7-l8 <= 2
l1+l2+l3-l4 <= 2
l1+l4+l5-l8 <= 2
l1+l2+l5-l6 <= 2
l1+l3+l5-l7 <= 2

## (4,8) vertices [Table 5] (22)
(1:1:1:1:0:0:0:0) (1:1:1:1:1:1:0:0) (2:2:1:1:1:1:0:0) (1:1:1:1:1:1:1:0)
(2:1:1:1:1:1:1:0) (2:2:2:1:1:1:1:0) (2:2:2:2:2:1:1:0) (3:3:2:2:2:1:1:0)
(3:3:2:2:2:2:2:0) (4:3:3:2:2:1:1:0) (1:1:1:1:1:1:1:1) (3:2:2:1:1:1:1:1)
(3:3:1:1:1:1:1:1) (3:3:3:3:3:1:1:1) (4:2:2:2:2:2:1:1) (4:4:3:3:3:1:1:1)
(4:4:4:2:2:2:1:1) (5:3:3:3:3:1:1:1) (5:5:5:3:3:1:1:1) (7:3:3:3:3:3:3:3)
(7:5:5:3:3:3:1:1) (7:7:7:3:3:3:3:3)

## (3,8) inequalities [Table 4] (31)
l2+l3+l4+l5 <= 2
l1+l2+l4+l7 <= 2
l1+l3+l4+l6 <= 2
l1+l2+l5+l6 <= 2
l1+l2-l3 <= 1
l2+l5-l7 <= 1
l1+l6-l7 <= 1
l2+l4-l6 <= 1
l1+l4-l5 <= 1
l3+l4-l7 <= 1
l1+l8 <= 1
l2-l3-l6-l7 <= 0
l4-l5-l6-l7 <= 0
l1-l3-l5-l7 <= 0
l2+l3+2l4-l5-l7+l8 <= 2
l1+l3+2l4-l5-l6+l8 <= 2
l1+2l2-l3+l4-l5+l8 <= 2
l1+2l2-l3+l5-l6+l8 <= 2
l1+l2-2l3-l4-l5 <= 0
l1-l2-l3+l6-2l7 <= 0
l1-l3-l4-l5+l8 <= 0
l1-l2-l3-l7+l8 <= 0
2l1-l2+l4-2l5-l6+l8 <= 1
l3+2l4-2l5-l6-l7+l8 <= 1
2l1-l2-l4+l6-2l7+l8 <= 1
2l1+l2-2l3-l4-l6+l8 <= 1
l1+2l2-2l3-l5-l6+l8 <= 1
2l1-2l2-l3-l4+l6-3l7+l8 <= 0
-l1+l3+2l4-3l5-2l6-l7+l8 <= 0
2l1+l2-3l3-2l4-l5-l6+l8 <= 0
l1+2l2-3l3-l4-2l5-l6+l8 <= 0

## (3,8) vertices [Table 6] (38, incl. the numerically-found (15:15:15:15:6:6:6:6)/28)
(1:1:1:0:0:0:0:0) (2:1:1:1:1:0:0:0) (1:1:1:1:1:1:0:0) (3:3:2:2:1:1:0:0)
(1:1:1:1:1:1:1:0) (2:2:1:1:1:1:1:0) (2:2:2:2:2:1:1:0) (3:1:1:1:1:1:1:0)
(3:3:3:3:1:1:1:0) (5:5:3:3:3:1:1:0) (1:1:1:1:1:1:1:1) (2:1:1:1:1:1:1:1)
(2:2:1:1:1:1:1:1) (3:3:3:3:3:1:1:1) (4:2:2:2:2:2:1:1) (4:3:2:2:1:1:1:1)
(4:4:2:2:2:2:1:1) (4:4:3:3:1:1:1:1) (4:4:4:4:2:1:1:1) (5:3:2:2:2:2:1:1)
(5:5:3:3:2:1:1:1) (5:5:5:5:2:2:2:2) (6:3:3:3:2:2:1:1) (6:5:5:5:2:2:1:1)
(6:6:3:3:3:2:2:2) (6:6:4:4:4:1:1:1) (7:5:5:5:2:2:2:2) (7:7:4:4:4:2:1:1)
(9:5:5:5:3:3:3:3) (9:6:4:4:4:3:3:3) (9:8:8:8:3:3:3:3) (9:9:5:5:3:3:3:2)
(9:9:9:9:4:4:2:2) (10:10:10:10:4:4:3:3) (11:6:6:5:5:5:2:2) (11:11:6:6:4:4:3:3)
(12:6:6:5:5:5:3:3) (12:12:7:7:4:4:4:4)

## Ranks 9-10 status [Sec 6.2.2, data NOT in text]
(3,9): 52 inequalities, rigorous. (4,9): 60 inequalities, 103 vertices; two vertices
numerical-only: (16:16:16:6:6:6:6:6:6)/21 and (20:14:14:14:14:4:4:4:4)/23 [= our A and B,
now closed]. (4,10): 125 constraints (completeness hinges on same two vertices).
(3,10): 93 inequalities. (5,10): 161 inequalities.
Data source named in paper: fen.bilkent.edu.tr/~murata/N-Representability.zip (dead).
Archival candidate: Springer ESM at doi:10.1007/s00220-008-0552-z.
