SLM_RSCRIPT = """\
#!/usr/bin/env Rscript

# Calculate copy number segmentation by SLM.
# Input: log2 coverage data in CNVkit's tabular format
# Output: the SLM output table

%(rlibpath)s
library('SLMSeg')

write("Loading probe coverages into a data frame", stderr())
tbl = read.delim("%(probes_fname)s")

# Drop any 0-depth bins
tbl = tbl[tbl$depth > 0,]

write("Segmenting the probe data", stderr())
segments <- HSLM(tbl$log2, tbl$start, 0.1, 0.00001, 1000000, 1)

write("Setting segment endpoints to original bin start/end positions", stderr())
out_tbl <- data.frame(id               = '%(sample_id)s',
                      chromosome       = tbl$chromosome,
                      start            = tbl$start,
                      end              = tbl$end,
                      nprobes          = 1,
                      log2             = segments[1, ],
                      stringsAsFactors = FALSE
                      )

write("Printing the CBS table to standard output", stderr())
write.table(out_tbl, '', sep='\t', row.names=FALSE)
"""
