# ADJUST THESE VARIABLES AS NEEDED TO SUIT YOUR COMPUTING ENVIRONMENT
# -------------------------------------------------------------------
# This variable specifies the number of threads to use for the
# parallelization. This could also be specified automatically using
# environment variables. For example, in SLURM, SLURM_CPUS_PER_TASK
# specifies the number of CPUs allocated for each task.
nc = 6

# R scripting front-end. Note that makeCluster sometimes fails to
# connect to a socker when using Rscript, so we are using the "R CMD
# BATCH" interface instead.
rexec = R CMD BATCH --no-save --no-restore

# AVOID EDITING ANYTHING BELOW THIS LINE
# --------------------------------------

# Files in preparation of comparing ngsLD
ngsprep = ./output/ngs_out/updog_format.RDS \
          ./output/ngs_out/llike.tsv.gz \
          ./output/ngs_out/pos.tsv.gz

# Plots looking at bias, se, and mse of various estimators of LD
mlesimplots = ./output/mle_plots/D_bias_nind_100_pA_50_pB_50.pdf \
	      ./output/mle_plots/D_bias_nind_100_pA_50_pB_75.pdf \
	      ./output/mle_plots/D_bias_nind_100_pA_90_pB_90.pdf \
	      ./output/mle_plots/D_bias_nind_1000_pA_50_pB_50.pdf \
	      ./output/mle_plots/D_bias_nind_1000_pA_50_pB_75.pdf \
	      ./output/mle_plots/D_bias_nind_1000_pA_90_pB_90.pdf \
	      ./output/mle_plots/D_mse_nind_100_pA_50_pB_50.pdf \
	      ./output/mle_plots/D_mse_nind_100_pA_50_pB_75.pdf \
	      ./output/mle_plots/D_mse_nind_100_pA_90_pB_90.pdf \
	      ./output/mle_plots/D_mse_nind_1000_pA_50_pB_50.pdf \
	      ./output/mle_plots/D_mse_nind_1000_pA_50_pB_75.pdf \
	      ./output/mle_plots/D_mse_nind_1000_pA_90_pB_90.pdf \
	      ./output/mle_plots/D_se_nind_100_pA_50_pB_50.pdf \
	      ./output/mle_plots/D_se_nind_100_pA_50_pB_75.pdf \
	      ./output/mle_plots/D_se_nind_100_pA_90_pB_90.pdf \
	      ./output/mle_plots/D_se_nind_1000_pA_50_pB_50.pdf \
	      ./output/mle_plots/D_se_nind_1000_pA_50_pB_75.pdf \
	      ./output/mle_plots/D_se_nind_1000_pA_90_pB_90.pdf \
	      ./output/mle_plots/r2_bias_nind_100_pA_50_pB_50.pdf \
	      ./output/mle_plots/r2_bias_nind_100_pA_50_pB_75.pdf \
	      ./output/mle_plots/r2_bias_nind_100_pA_90_pB_90.pdf \
	      ./output/mle_plots/r2_bias_nind_1000_pA_50_pB_50.pdf \
	      ./output/mle_plots/r2_bias_nind_1000_pA_50_pB_75.pdf \
	      ./output/mle_plots/r2_bias_nind_1000_pA_90_pB_90.pdf \
	      ./output/mle_plots/r2_mse_nind_100_pA_50_pB_50.pdf \
	      ./output/mle_plots/r2_mse_nind_100_pA_50_pB_75.pdf \
	      ./output/mle_plots/r2_mse_nind_100_pA_90_pB_90.pdf \
	      ./output/mle_plots/r2_mse_nind_1000_pA_50_pB_50.pdf \
	      ./output/mle_plots/r2_mse_nind_1000_pA_50_pB_75.pdf \
	      ./output/mle_plots/r2_mse_nind_1000_pA_90_pB_90.pdf \
	      ./output/mle_plots/r2_se_nind_100_pA_50_pB_50.pdf \
	      ./output/mle_plots/r2_se_nind_100_pA_50_pB_75.pdf \
	      ./output/mle_plots/r2_se_nind_100_pA_90_pB_90.pdf \
	      ./output/mle_plots/r2_se_nind_1000_pA_50_pB_50.pdf \
	      ./output/mle_plots/r2_se_nind_1000_pA_50_pB_75.pdf \
	      ./output/mle_plots/r2_se_nind_1000_pA_90_pB_90.pdf

# Plots looking at accuracy of standard errors and qq-plots for normality
mleqqplots = ./output/mle_se_plots/qq_nind100_pA50_pB50_r0.pdf \
	     ./output/mle_se_plots/qq_nind100_pA90_pB90_r0.pdf \
             ./output/mle_se_plots/comnorm_se_est.pdf \
             ./output/mle_se_plots/mle_se_est.pdf

# uitdewilligen data
uitdat = ./data/NewPlusOldCalls.headed.vcf \
         CSV-file\ S1\ -\ Sequence\ variants\ filtered\ DP15.csv

all : mle ngsLD

# Pairwise LD estimation simulations ---------------
./output/mle/mle_sims_out.csv : ./code/mle_sims.R
	mkdir -p ./output/mle
	mkdir -p ./output/rout
	$(rexec) '--args nc=$(nc)' $< ./output/rout/$(basename $(notdir $<)).Rout

$(mlesimplots) : ./code/mle_sim_plots.R ./output/mle/mle_sims_out.csv
	mkdir -p ./output/mle_plots
	mkdir -p ./output/rout
	$(rexec) $< ./output/rout/$(basename $(notdir $<)).Rout

$(mleqqplots) : ./code/mle_sim_se_plots.R ./output/mle/mle_sims_out.csv
	mkdir -p ./output/mle_se_plots
	mkdir -p ./output/rout
	$(rexec) $< ./output/rout/$(basename $(notdir $<)).Rout

.PHONY : mle
mle : $(mlesimplots) $(mleqqplots)


# Comparing ngsLD to ldsep --------------------
$(ngsprep) : ./code/ngs_sim_data.R
	mkdir -p ./output/ngs_out
	mkdir -p ./output/rout
	rm ./output/ngs_out/llike.tsv.gz -f
	rm ./output/ngs_out/pos.tsv.gz -f
	$(rexec) $< ./output/rout/$(basename $(notdir $<)).Rout

./output/ngs_out/lsep_out.tsv ./output/ngs_out/ngs_fit.tsv : ./code/ngs_sims.R $(ngsprep)
	mkdir -p ./output/ngs_out
	mkdir -p ./output/rout
	$(rexec) $< ./output/rout/$(basename $(notdir $<)).Rout

./output/fig/D_ngsld_ldsep.pdf : ./code/ngsld_vs_ldsep.R ./output/ngs_out/lsep_out.tsv ./output/ngs_out/ngs_fit.tsv
	mkdir -p ./output/fig
	mkdir -p ./output/rout
	$(rexec) $< ./output/rout/$(basename $(notdir $<)).Rout

.PHONY : ngsLD
ngsLD : ./output/fig/D_ngsld_ldsep.pdf

# Uitdewilligen analysis
$(uitdat) :
	wget --directory-prefix=data -nc https://doi.org/10.1371/journal.pone.0062355.s004
	mv ./data/journal.pone.0062355.s004 ./data/journal.pone.0062355.s004.gz
	wget --directory-prefix=data -nc https://doi.org/10.1371/journal.pone.0062355.s007
	mv ./data/journal.pone.0062355.s007 ./data/journal.pone.0062355.s007.gz
	7z e ./data/journal.pone.0062355.s004.gz -o./data/
	7z e ./data/journal.pone.0062355.s007.gz -o./data/

.PHONY : uit
uit : $(uitdat)


