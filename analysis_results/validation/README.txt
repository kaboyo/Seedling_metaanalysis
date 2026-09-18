MANUSCRIPT AUDIT — Meta-analysis of the Seedling Microbiome
Recomputation of every headline statistic in Manuscript.docx, 2026-08-04.

READ FIRST
  revision_memo.md ............. every change, with justification, and every
                                 claim that must be softened or dropped
  discrepancy_audit.csv ........ 22 claims: manuscript value vs recomputed
                                 value vs verdict vs mechanism

DELIVERABLES
  Manuscript_revised.docx ...... corrected manuscript. 17 paragraphs rewritten,
                                 13 inserted, 3 tables corrected, new Table 4b
                                 (marginal PERMANOVA), 4 captions rewritten.
                                 The ORIGINAL Manuscript.docx is unmodified.
  claim_evidence_register_verified.csv
                                 the 5 original claims re-graded against
                                 recomputed evidence, expanded to 7
  verify_manuscript.py ......... the full recomputation, runnable end to end

  figures/ ..................... fig1 composition, fig2 alpha, fig3 beta,
                                 fig4 indicator taxa (300 dpi)
  supplementary/ ............... Tables S1-S6 plus the three PERMANOVA /
                                 PERMDISP result files

THE FOUR FINDINGS THAT MATTER
  1. Table 4 was computed on ASV-level top-100 RAW COUNTS at n=1,553, not
     genus-level relative abundances at n=1,883 as the caption states. All
     five F, all five R2 and all five PERMDISP F reproduce to 4 s.f. on that
     basis. The n column came from the genus object; the statistics did not.
  2. Shannon Kruskal-Wallis is H=36.357 (p=1.64e-9, n=1,883), not H=38.2 at
     n=1,886. The published value requires scoring Shannon=0 for three
     samples where the index is undefined.
  3. 98.79% of 119,135 ASVs occur in exactly one study; median between-study
     Jaccard is 0.000. ASV-level cross-study pooling does not hold, and that
     is the basis Table 4 used.
  4. The stage effect on beta diversity falls from R2=0.012 (univariate) to
     R2=0.003 in a marginal model with study and host species.

  The paper's central ordering - host species and cultivar above
  developmental stage - HOLDS and strengthens under correct analysis.

OPEN DECISION
  PERMDISP p is reported as 0.003, the floor of the 299-permutation scheme
  used. At 999 permutations it becomes 0.001 and Table 4 should be updated.
