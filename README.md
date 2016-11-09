# SRST2_toolkit
Supplementary Python and R scripts for processing results from [SRST2](https://github.com/katholt/srst2).

## A list of scripts
* [clustering_allele_variants.py](#clustering_allele_variants)
* [rm_seqID.py](#rm_seqID)
* [get_alleleSeq.py](#get_alleleSeq)
* [evaluate_argannot_seqlen.py](#evaluate_argannot_seqlen)
* [find_argannot_redundancy.R](#find_argannot_redundancy)

## Manual
### <a name="clustering_allele_variants"></a>clustering_allele_variants.py
This script assigns identifiers to allele variants in the compiled result of SRST2 in accordance with consensus sequences clustered by CD-HIT-EST. A number from cluster IDs will be appended to the end of the corresponding allele name according to results of the sequence clustering. Renamed alleles are printed to STDOUT.

### <a name="rm_seqID"></a>rm_seqID.py
This script removes seqIDs and following characters (\* or ?) from the compiled result of SRST2.

### <a name="get_alleleSeq"></a>get_alleleSeq.py
This script retrieves allele sequences from the ARG-Annot database according to a list of allele names.

### <a name="evaluate_argannot_seqlen"></a>evaluate_argannot_seqlen.py
This script finds discrepancies of sequence lengths in the SRST2-formatted ARG-Annot database.

### <a name="find_argannot_redundancy"></a>find_argannot_redundancy.R
This script finds out redundant sequences in the SRST2-formatted ARG-Annot database.