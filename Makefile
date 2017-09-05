## Take ROS/MAP and Brain Cloud data

RAW_FILE := /broad/dejagerlab/cogdec/BrainExpn/pipeline/RSEM/Expn/RSEM_Phase1_Phase2_gene_expected_count_annt.txt

GENCODE := gencode.v19.genes.patched_contigs.gtf.gz

PHENOTYPE := phenotypes/pheno_cov_n3033_032315.csv

SAMP_INFO := sample.info.txt

CHR := $(shell seq 1 22)

## make a list of coding genes
step1: data/coding.genes.txt.gz data/qc.samples.txt.gz data/pheno.txt.gz

data/coding.genes.txt.gz: ./make.coding.genes.R $(GENCODE) data/rnaseq.rows.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	./run.sh "Rscript --vanilla $^ $@"

data/qc.samples.txt.gz: ./make.qc.samples.R data/rnaseq.samples.txt.gz $(SAMP_INFO)
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	./run.sh "Rscript --vanilla $^ $@"

data/rnaseq.rows.txt.gz:
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	cat $(RAW_FILE) | awk -F'\t' '{ print $$1 FS NR }' | gzip > $@

data/rnaseq.samples.txt.gz:
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	cat $(RAW_FILE) | head -n1 | tr '\t' '\n' | awk -F'\t' '{ print $$1 FS NR }' | tail -n+4 | gzip > $@

data/pheno.txt.gz: make.matched.phenotype.R $(PHENOTYPE)  data/qc.samples.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	./run.sh "Rscript --vanilla $^ $@"

## separate chr by chr
step2: $(foreach chr, $(CHR), data/rnaseq/chr$(chr)-probes.txt.gz data/rnaseq/chr$(chr)-samples.txt.gz data/rnaseq/chr$(chr)-count.txt.gz)

data/rnaseq/chr%-probes.txt.gz: data/coding.genes.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	zcat $^ | awk '$$1 == $*' | gzip > $@

data/rnaseq/chr%-samples.txt.gz: data/pheno.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	zcat $^ | cut -f2- | gzip > $@

data/rnaseq/chr%-count.txt.gz: data/pheno.txt.gz data/rnaseq/chr%-probes.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	cat $(RAW_FILE) | awk -F'\t' -vROWS=$$(zcat data/rnaseq/chr$*-probes.txt.gz | awk 'NR > 1 { printf "," } { printf $$NF }') -f util_subset_rows.awk | awk -F'\t' -vCOLS=$$(zcat data/pheno.txt.gz | tail -n+2 | cut -f1 | awk 'NR > 1 { printf "," } { printf $$1 }') -f util_subset_cols.awk | awk -f util_transpose.awk | gzip > $@

