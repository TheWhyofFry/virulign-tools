# virulign-tools

These scripts are designed to try and augment the performance of virulign. Sometimes virulign cannot resolve a sequence and exists with a FrameShift error. While it can be argued that perhaps the assembly/input sequence is incorrect, it can be shown that virulign fails for \~200 env sequences tested out of ~4000 curated HIV-1 alignments.  The preliminary procedure to circumvent this is:

1. Get all the sequences that successfully align to the respective genes. We will allow zero frameshifts
2. Use this set of aligned sequences to build a blast database (GENEDB)
3. Run virulign on a sequence set
 - For every sequence that can't be aligned
   * Do a blast against GENEDB
   * Take the top hit and create a virulign XML from it
   * Run virulign with an added HXB2 reference sequence
   * Repeat until success/limit number of failures
4. Should we ultimately use HXB2.gb and AGA?



