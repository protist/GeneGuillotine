GeneGuillotine
==============
This script reads in a GTF file (e.g. from Cuffmerge), and prevents transcripts
from overlapping multiple genes, according to a second (reference) GFF. This
might be useful if a downstream program (e.g. DEXSeq) requires each transcript
to be separate and not overlapping with its neighbours.

The position of the split is determined by the genes in the reference GFF. The
default is to constrain each transcript to the limits of the CDS. Transcripts
that lie wholly within intergenic regions will be kept.

Arguments and input files
-------------------------
**Help** (`-h`, `--help`). Display help.

**Verbosity** (`-v`, `--verbose [2]`). Run verbosely. The default level of verbosity
reports at a chromosomal level. Optionally, specify level 2 for gene-level
reporting.

**RNA-seq data** (`-i`, `--input GTF_FILE`). The GTF parser is designed for use with
Cuffmerge GTF files. Hence, the canonical usage is to use the Tuxedo pipeline
(map reads with [Tophat](http://tophat.cbcb.umd.edu/), create gene models with
[Cufflinks](http://cufflinks.cbcb.umd.edu/), and merge samples with
[Cuffmerge](http://cufflinks.cbcb.umd.edu/)).

**Reference gene models** (`-g`, `--ref_gff GFF_FILE`). The GFF parser is designed
for use with GFF files from [EuPathDB](http://eupathdb.org). It only parses
features marked as `CDS` (and `tRNA` and `rRNA`), since UTR information is not
available for all genes.

**Minimal split** (`-m`, `--minimal_split`). This optional flag tells the script to
only split transcripts that overlap multiple genes (i.e. the first part of the
script). This prevents the second part of the script from running, which would
otherwise truncate transcripts that lie on adjacent genes, but overlap with each
other.

Output
------
**Output file** (`-o`, `--output OUTPUT_FILE`). Gene IDs from the reference GFF file
are written to the transcripts, if they cover a single gene. In other cases, the
nearest gene is recorded. e.g. for intergenic transcripts, "after_GENE_ID" or
"before_GENE_ID"; and for transcripts that lie before the first gene or after
the last gene, "before_first_GENE_ID" or "after_last_GENE_ID". If there are no
genes on the reference contig, the gene ID is "No_genes_on_ref_contig".

If transcripts cover multiple genes, then the transcripts will be renamed to
"TRANSCRIPT_ID", "TRANSCRIPT_ID:2", "TRANSCRIPT_ID:3", etc. These strings can be
easily modified from the code.
