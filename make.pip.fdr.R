#!/usr/bin/env Rscript
argv <- commandArgs(trailingOnly = TRUE)
alt.file <- argv[1]
null.file <- argv[2]
out.file <- argv[3]

library(qvalue)
options(stringsAsFactors = FALSE)

null.pip <- scan(null.file)
alt.tab <- read.table(alt.file, sep = '\t', col.names = c('ensg', 'max.pip'))
p.value <- empPvals(alt.tab$max.pip, null.pip)
q.value <- qvalue(p.value)

out.hdr <- gsub(pattern = '.gz', replacement = '', out.file)
pdf(file = paste(out.hdr, '-pvalue.pdf', sep = ''))
hist(p.value, 30)
dev.off()

pdf(file = paste(out.hdr, '-qvalue.pdf', sep = ''))
plot(q.value)
dev.off()

out <- data.frame(alt.tab, p.val = p.value, q.val = q.value$qvalues)
write.table(out, file = gzfile(out.file), quote = FALSE, row.names = FALSE,
            col.names = TRUE)
