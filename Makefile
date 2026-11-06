# OSCC dendritic-cell atlas -- reproducible workflow
# -------------------------------------------------
# R is invoked with an explicit interpreter so the workflow also runs on
# Windows where `Rscript` is not always on PATH.  Override with
#   make RSCRIPT=/path/to/Rscript all
RSCRIPT ?= Rscript

STEPS := 01_load_qc 02_integrate_cluster 03_annotate_major 04_dc_subset \
         05_dc_annotate_de 06_enrichment 07_trajectory 08_tf_activity 09_composition

.PHONY: all clean dirs params session $(STEPS)

all: $(STEPS)

dirs:
	@mkdir -p results/objects results/tables results/figures results/logs

params: 01_load_qc
	@echo "parameters are written by R/01_load_qc.R"

01_load_qc: | dirs
	$(RSCRIPT) R/01_load_qc.R 2>&1 | tee results/logs/01_load_qc.log

02_integrate_cluster: 01_load_qc
	$(RSCRIPT) R/02_integrate_cluster.R 2>&1 | tee results/logs/02_integrate_cluster.log

03_annotate_major: 02_integrate_cluster
	$(RSCRIPT) R/03_annotate_major.R 2>&1 | tee results/logs/03_annotate_major.log

04_dc_subset: 03_annotate_major
	$(RSCRIPT) R/04_dc_subset.R 2>&1 | tee results/logs/04_dc_subset.log

05_dc_annotate_de: 04_dc_subset
	$(RSCRIPT) R/05_dc_annotate_de.R 2>&1 | tee results/logs/05_dc_annotate_de.log

06_enrichment: 05_dc_annotate_de
	$(RSCRIPT) R/06_enrichment.R 2>&1 | tee results/logs/06_enrichment.log

07_trajectory: 05_dc_annotate_de
	$(RSCRIPT) R/07_trajectory.R 2>&1 | tee results/logs/07_trajectory.log

08_tf_activity: 05_dc_annotate_de
	$(RSCRIPT) R/08_tf_activity.R 2>&1 | tee results/logs/08_tf_activity.log

09_composition: 05_dc_annotate_de
	$(RSCRIPT) R/09_composition.R 2>&1 | tee results/logs/09_composition.log

session:
	$(RSCRIPT) -e "writeLines(capture.output(sessionInfo()), 'docs/sessionInfo.txt')"

clean:
	rm -rf results/objects results/figures results/tables results/logs
