## Take ROS/MAP and Brain Cloud data

RAW_FILE := /broad/dejagerlab/cogdec/BrainExpn/pipeline/RSEM/Expn/RSEM_Phase1_Phase2_gene_expected_count_annt.txt

GENCODE := gencode.v19.genes.patched_contigs.gtf.gz

PHENOTYPE := phenotypes/pheno_cov_n3033_032315.csv

SAMP_INFO := sample.info.txt

PLINK_HDR := rosmap-geno/gen/impute/rosmap1709-chr

CHR := $(shell seq 1 22)

################################################################
## make a list of coding genes and sample Q/C
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


################################################################
## separate chr by chr
step2: $(foreach chr, $(CHR), data/rnaseq/chr$(chr)-genes.txt.gz data/rnaseq/chr$(chr)-samples.txt.gz data/rnaseq/chr$(chr)-count.txt.gz) data/rnaseq/size_factors.txt.gz data/rnaseq/peer_factors.txt.gz

data/rnaseq/chr%-genes.txt.gz: data/coding.genes.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	zcat $^ | awk '$$1 == $*' | gzip > $@

data/rnaseq/chr%-samples.txt.gz: data/pheno.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	zcat $^ | cut -f2- | gzip > $@

data/rnaseq/chr%-count.txt.gz: data/pheno.txt.gz data/rnaseq/chr%-genes.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	cat $(RAW_FILE) | awk -F'\t' -vROWS=$$(zcat data/rnaseq/chr$*-genes.txt.gz | awk 'NR > 1 { printf "," } { printf $$NF }') -f util_subset_rows.awk | awk -F'\t' -vCOLS=$$(zcat data/pheno.txt.gz | tail -n+2 | cut -f1 | awk 'NR > 1 { printf "," } { printf $$1 }') -f util_subset_cols.awk | awk -F'\t' -f util_transpose.awk | gzip > $@

## estimate sequencing depth
data/rnaseq/size_factors.txt.gz: data/pheno.txt.gz
	./run.sh ./make.size.factor.R data/rnaseq/chr -count.txt.gz $^ $@

## estimate expression PEER factors
data/rnaseq/peer_factors.txt.gz: data/pheno.txt.gz data/rnaseq/size_factors.txt.gz
	./run.sh ./make.expr.peer.R $^ $@

################################################################
## pre-generate temporary data
TEMPDIR := /broad/hptmp/ypp/AD/twas/qtl/
CHUNK := 33
NCTRL := 3

step3: $(TEMPDIR) data/rnaseq/size_factors.txt.gz jobs/qtl-data-jobs.txt.gz

jobs/qtl-data-jobs.txt.gz: $(foreach chr, $(CHR), jobs/temp-qtl-data-$(chr)-jobs.txt.gz)
	zcat $^ | awk 'system("[ ! -f " $$NF ".x.ft ]") == 0' | gzip > $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N qtl.data -binding "linear:1" -q short -l h_vmem=4g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@
	rm $^

jobs/temp-qtl-data-%-jobs.txt.gz: data/rnaseq/chr%-genes.txt.gz data/rnaseq/chr%-count.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	mkdir -p $(TEMPDIR)/$*/
	./make_job_segments.awk -vNTOT=$$(zcat $< | wc -l) -vCHUNK=$(CHUNK) | awk '{ print "./make.data.qtl.R" FS "$^" FS $$1 FS $(NCTRL) FS ("$(PLINK_HDR)" $*) FS ("$(TEMPDIR)/$*/data-" NR) }' | gzip > $@

$(TEMPDIR):
	[ -d $@ ] || mkdir -p $@


################################################################
step4: $(TEMPDIR) jobs/qtl-run-jobs.txt.gz

jobs/qtl-run-jobs.txt.gz: $(foreach chr, $(CHR), jobs/temp-qtl-run-$(chr)-jobs.txt.gz)
	zcat $^ | awk 'system("[ ! -f " $$NF ".x.ft ]") == 0' | gzip > $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N TWAS.qtl -binding "linear:1" -q short -l h_vmem=4g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@
	rm $^

jobs/temp-qtl-run-%-jobs.txt.gz: data/rnaseq/chr%-genes.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	./make_job_segments.awk -vNTOT=$$(zcat $< | wc -l) -vCHUNK=$(CHUNK) | awk '{ print "./make.summary.qtl.R" FS ("$(TEMPDIR)/$*/data-" NR) FS ("result/qtl/$*/b" NR) }' | gzip > $@

step4-long: jobs/qtl-run-jobs-long.txt.gz

jobs/qtl-run-jobs-long.txt.gz:
	zcat jobs/qtl-run-jobs.txt.gz | awk 'system("[ ! -f " $$NF ".genes.gz ]") == 0' | gzip > $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N TWAS.qtl -binding "linear:1" -q long -l h_vmem=4g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@



################################################################
step5: jobs/poly-jobs.txt.gz

step5-resubmit: jobs/poly-jobs-resubmit.txt.gz

RSEED := 1991 1331 1771

step5_jobs := $(foreach chr, $(CHR), $(shell ls -1 result/qtl/$(chr)/*.resid.gz | awk -F'/' '{ gsub(/.resid.gz/,"",$$NF); print "jobs/temp_poly-" $$(NF -1 ) "/" $$NF "-jobs.txt" }'))

jobs/poly-jobs.txt.gz: $(step5_jobs)
	cat $^ | gzip > $@
	[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N TWAS.poly -binding "linear:1" -q short -l h_vmem=4g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@
	rm -r jobs/temp_poly*

jobs/poly-jobs-resubmit.txt.gz:
	zcat jobs/poly-jobs.txt.gz | awk 'system("[ ! -f " $$NF ".effect.gz ] && [ ! -f " $$(NF-1) ".effect.gz ]") == 0' | gzip > $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N TWAS.poly -binding "linear:1" -q short -l h_vmem=8g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@


# % = $(chr)/$(resid) e.g., 21/b1-peer-sqtl
jobs/temp_poly-%-jobs.txt: result/qtl/%.resid.gz
	@printf "" > $@
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@[ -d $(dir result/poly/$*) ] || mkdir -p $(dir result/poly/$*)
	@[ -f result/poly/$*.effect.gz ] || echo ./make.polygenic.R $< $(TEMPDIR)/$(shell echo $* | awk -F'/' '{ print $$1 }')/data-$(shell basename $* | awk -F'-' '{ gsub(/b/, "", $$1); print $$1 }') result/poly/$* >> $@
	@for rseed in $(RSEED); do [ -d $(dir result/poly/perm.$${rseed}/$*) ] || mkdir -p $(dir result/poly/perm.$${rseed}/$*); done
	@for rseed in $(RSEED); do [ -f result/poly/perm.$${rseed}/$*.effect.gz ] || echo ./make.polygenic.R $< $(TEMPDIR)/$(shell echo $* | awk -F'/' '{ print $$1 }')/data-$(shell basename $* | awk -F'-' '{ gsub(/b/, "", $$1); print $$1 }') result/poly/perm.$${rseed}/$* $${rseed} >> $@; done





################################################################
## Utilities
PLINKZIP := https://www.cog-genomics.org/static/bin/plink170815/plink_linux_x86_64.zip

bin/plink:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	curl $(PLINKZIP) -o bin/plink.zip
	unzip bin/plink.zip -d bin/

