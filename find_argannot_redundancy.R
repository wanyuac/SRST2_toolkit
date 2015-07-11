# Find out redundant sequences in SRST2's ARG-Annot database
# Yu Wan (wanyuac@gmail.com, https://github.com/wanyuac), 10 July 2015
# License: GNU GPL 2.1

kFILE <- "blast_argannot.txt"  # Alignment summaries in format 6 to be imported to R from megaBLAST
kOUTPUT <- "redundant_seqs.txt"  # the file name of the output
kSEQLEN <- "seqlen.txt"  # a list of sequence lengths described in the original ARG-Annot database

removeReplicates <- function(t) {
  # This function removes the second member of the duples: (a, b) and (b, a)
  n <- nrow(t)
  marker <- rep(TRUE, times = n)  # initialise the marker for selecting rows
  for (i in 1 : n) {
    if (marker[i]) {  # if this marker has not been marked as FALSE
      qseqid <- t[i, "qseqid"]  # get the current qseqid
      selected <- which(t[, "sseqid"] == qseqid)  # selected row indices
      j <- selected[t[selected, "qseqid"] == t[i, "sseqid"]]
      marker[j] <- FALSE  # unselect the row whose sseqid matches qseqid[i]
    }
  }
  return(t[marker,])
}

tab <- read.csv(kFILE, header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
sel <- (tab[,"length"] == tab[, "slen"]) & (tab[,"length"] == tab[, "qlen"]) & 
  (tab[,"pident"] == 100) & (tab[, "evalue"] <= 0.001)  # the definition of an exact match, which must have the highest bit scores
exact.matches <- tab[sel, ]  # remove imperfect matches
redundancy <- exact.matches[exact.matches[, "qseqid"] != exact.matches[, "sseqid"],]  # remove self-matches
result <- removeReplicates(redundancy)
# This relationship must be shown here: nrow(result) = 1/2 * nrow(redundancy)
write.table(result, file = kOUTPUT, quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)
