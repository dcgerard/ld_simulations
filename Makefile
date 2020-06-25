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

# Plots showing flexibility of proportional bivariate normal distribution
normplots = ./output/compare_norm/normdist.pdf \
            ./output/compare_norm/randnorm.pdf

# uitdewilligen data
uitdat = ./data/NewPlusOldCalls.headed.vcf \
         ./data/CSV-file\ S1\ -\ Sequence\ variants\ filtered\ DP15.csv \
         ./data/journal.pone.0062355.s009.xls \
         ./data/journal.pone.0062355.s010.xls

# Subset of uitdewilligen data
uitsnps = ./output/uit/refmat_suc.RDS \
          ./output/uit/sizemat_suc.RDS \
          ./output/uit/uit_suc.csv

# Uitdewilligen LD estimates on subset
uitld = ./output/uit/ldest_hap_genolike.RDS \
        ./output/uit/ldest_hap_geno.RDS \
        ./output/uit/ldest_comp_genolike.RDS \
        ./output/uit/ldest_comp_genolike_flex.RDS \
        ./output/uit/ldest_comp_geno.RDS

# Plots from Uitdewilligen LD estimates on subset
uitfig = ./output/uit/uit_fig/heat_comp_geno.pdf \
         ./output/uit/uit_fig/heat_comp_genolike.pdf \
         ./output/uit/uit_fig/heat_comp_genolike_flex.pdf \
         ./output/uit/uit_fig/heat_hap_geno.pdf \
         ./output/uit/uit_fig/heat_hap_genolike.pdf \
         ./output/uit/uit_fig/uit_pairs.pdf \
         ./output/uit/uit_fig/diff11.pdf \
         ./output/uit/uit_fig/diff22.pdf \
         ./output/uit/uit_fig/box12.pdf

# Properties of Uitdewilligen SNPs
uitprop = ./output/uit/uit_fig/maf.pdf \
          ./output/uit/uit_fig/readdepth.pdf

# Raw data from McAllister et al (2017)
mcadat = ./data/gerardii/McAllister.Miller.all.mergedRefGuidedSNPs.vcf.gz \
         ./data/gerardii/McAllister_Miller_Locality_Ploidy_Info.csv

# Subset of McAllister data
mcasnps = ./output/mca/refmat_hex.RDS \
          ./output/mca/sizemat_hex.RDS \
          ./output/mca/refmat_non.RDS \
          ./output/mca/sizemat_non.RDS

# Updog fits of McAllister data
mcaupdog = ./output/mca/updog_fits_hex.RDS \
           ./output/mca/updog_fits_non.RDS

# McAllister LD estimates from hexaploids
mcaldhex = ./output/mca/ldest_hap_genolike_hex.RDS \
           ./output/mca/ldest_hap_geno_hex.RDS \
           ./output/mca/ldest_comp_genolike_hex.RDS \
           ./output/mca/ldest_comp_genolike_flex_hex.RDS \
           ./output/mca/ldest_comp_geno_hex.RDS

# McAllister LD estimates from nonaploids
mcaldnon = ./output/mca/ldest_hap_genolike_non.RDS \
           ./output/mca/ldest_hap_geno_non.RDS \
           ./output/mca/ldest_comp_genolike_non.RDS \
           ./output/mca/ldest_comp_genolike_flex_non.RDS \
           ./output/mca/ldest_comp_geno_non.RDS

all : mle ngsLD uit mca norm

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

./output/ngs_out/D_ngsld_ldsep.pdf : ./code/ngsld_vs_ldsep.R ./output/ngs_out/lsep_out.tsv ./output/ngs_out/ngs_fit.tsv
	mkdir -p ./output/ngs_out
	mkdir -p ./output/rout
	$(rexec) $< ./output/rout/$(basename $(notdir $<)).Rout

.PHONY : ngsLD
ngsLD : ./output/ngs_out/D_ngsld_ldsep.pdf

# Uitdewilligen analysis
$(uitdat) :
	wget --directory-prefix=data -nc https://doi.org/10.1371/journal.pone.0062355.s004
	mv ./data/journal.pone.0062355.s004 ./data/journal.pone.0062355.s004.gz ## variant annotations
	wget --directory-prefix=data -nc https://doi.org/10.1371/journal.pone.0062355.s007
	mv ./data/journal.pone.0062355.s007 ./data/journal.pone.0062355.s007.gz ## read-depth data
	wget --directory-prefix=data -nc https://doi.org/10.1371/journal.pone.0062355.s009
	mv ./data/journal.pone.0062355.s009 ./data/journal.pone.0062355.s009.xls ## super scaffold annotations
	wget --directory-prefix=data -nc https://doi.org/10.1371/journal.pone.0062355.s010
	mv ./data/journal.pone.0062355.s010 ./data/journal.pone.0062355.s010.xls ## contig annotations
	7z e ./data/journal.pone.0062355.s004.gz -o./data/
	7z e ./data/journal.pone.0062355.s007.gz -o./data/

$(uitsnps) : ./code/uit_extract.R $(uitdat)
	mkdir -p ./output/uit
	mkdir -p ./output/rout
	$(rexec) $< ./output/rout/$(basename $(notdir $<)).Rout

./output/uit/uit_updog_fit.RDS : ./code/uit_fit_updog.R $(uitsnps)
	mkdir -p ./output/uit
	mkdir -p ./output/rout
	$(rexec) '--args nc=$(nc)' $< ./output/rout/$(basename $(notdir $<)).Rout

$(uitld) : ./code/uit_est_ld.R ./output/uit/uit_updog_fit.RDS
	mkdir -p ./output/uit
	mkdir -p ./output/rout
	$(rexec) '--args nc=$(nc)' $< ./output/rout/$(basename $(notdir $<)).Rout

$(uitfig) : ./code/uit_compare_ld.R $(uitld)
	mkdir -p ./output/uit
	mkdir -p ./output/uit/uit_fig
	mkdir -p ./output/rout
	$(rexec) $< ./output/rout/$(basename $(notdir $<)).Rout

$(uitprop) : ./code/uit_prop.R ./output/uit/uit_updog_fit.RDS
	mkdir -p ./output/uit
	mkdir -p ./output/uit/uit_fig
	mkdir -p ./output/rout
	$(rexec) $< ./output/rout/$(basename $(notdir $<)).Rout

.PHONY : uit
uit : $(uitfig) $(uitprop)


# Analysis of Mcallister et al (2017) Data
$(mcasnps) : ./code/mca_extract.R $(mcadat)
	mkdir -p ./output/mca
	mkdir -p ./output/rout
	$(rexec) $< ./output/rout/$(basename $(notdir $<)).Rout

$(mcaupdog) : ./code/mca_fit_updog.R $(mcasnps)
	mkdir -p ./output/mca
	mkdir -p ./output/rout
	$(rexec) '--args nc=$(nc)' $< ./output/rout/$(basename $(notdir $<)).Rout

$(mcaldhex) : ./code/mca_est_ld_hex.R $(mcaupdog)
	mkdir -p ./output/mca
	mkdir -p ./output/rout
	$(rexec) '--args nc=$(nc)' $< ./output/rout/$(basename $(notdir $<)).Rout

$(mcaldnon) : ./code/mca_est_ld_non.R $(mcaupdog)
	mkdir -p ./output/mca
	mkdir -p ./output/rout
	$(rexec) '--args nc=$(nc)' $< ./output/rout/$(basename $(notdir $<)).Rout

.PHONY : mca
mca : $(mcaldhex) $(mcaldnon)


# Proportional normal distribution plots
$(normplots) : ./code/pbnorm_flex.R
	mkdir -p ./output/compare_norm
	mkdir -p ./output/rout
	$(rexec) $< ./output/rout/$(basename $(notdir $<)).Rout

.PHONY : norm
norm: $(normplots)
