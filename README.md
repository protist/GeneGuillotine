Split_transcripts
=================
This script reads in a gtf file (e.g. from cuffmerge), and splits it into
non-overlapping genes, with one transcript per gene according to a second
(reference) gff. This might be useful if a downstream program (e.g. DEXSeq)
requires each gene to be separate and not overlapping with its neighbours.
