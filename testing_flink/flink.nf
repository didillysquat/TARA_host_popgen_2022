#!/usr/bin/env nextflow

// The aim of this script is to run the flink analysis on the POC and POR genomes.
// We started it running as a single analysis across the entire genome but this takes a very long time
// So rather we will try to split it up similar to how they do it in their paper.
// They say that 10^4 variants is enough to estimate parameters so we will at first
// run an flink analysis on subset vcfs that are only a single scaffold for those scaffolds
// that hold more than 10^4 varaiants.

// We have output a list for POC and POR of the scaffolds that have >10^4 variants
// /home/humebc/projects/paper_4/testing_flink/POC.scaff.10000.list.txt
// /home/humebc/projects/paper_4/testing_flink/POR.scaff.10000.list.txt

// The vcfs that we will work with are here:
// /home/humebc/projects/paper_4/testing_flink/POC/Pocillopora_meandrina_v3_11Islands_New_maf05_minQ30_biallelic_nomiss.samples.removed.vcf.gz
// /home/humebc/projects/paper_4/testing_flink/POR/Porites_lobata_v3_11Islands_New_maf05_minQ30_biallelic_nomiss.samples.removed.vcf.gz

// For the flink analysis, we will need to prepare the allele counts input file using Atlas. For this
// we need the input file that defines which population each of the samples belongs to.
// These are here:

poc_sample_populations = file("/home/humebc/projects/paper_4/testing_flink/POC/POC.samplesPopulations.txt")
por_sample_populations = file("/home/humebc/projects/paper_4/testing_flink/POR/POR.samplesPopulations.txt")
poc_sample_populations_lineage_island = file("/home/humebc/projects/paper_4/testing_flink/POC/POC.samplesPopulations.lineage.island.txt")
por_sample_populations_lineage_island = file("/home/humebc/projects/paper_4/testing_flink/POR/POR.samplesPopulations.lineage.island.txt")

calculate_param_averages_lineage_as_populations = file("/home/humebc/projects/paper_4/testing_flink/bin/calc_param_averages_lineage_as_populations.r")

awk_create_flink_input_lineage_island = file("/home/humebc/projects/paper_4/testing_flink/bin/create_flink_input_lineage_island.awk")

// For each of the scaffolds to process we will subset the vcf using bcftools
// We will then create the input files for each of these using Atlas and start a flink analysis for each

poc_vcf = file("/home/humebc/projects/paper_4/testing_flink/POC/Pocillopora_meandrina_v3_11Islands_New_maf05_minQ30_biallelic_nomiss.samples.removed.vcf.gz")
poc_vcf_tbi = file("/home/humebc/projects/paper_4/testing_flink/POC/Pocillopora_meandrina_v3_11Islands_New_maf05_minQ30_biallelic_nomiss.samples.removed.vcf.gz.tbi")
por_vcf = file("/home/humebc/projects/paper_4/testing_flink/POR/Porites_lobata_v3_11Islands_New_maf05_minQ30_biallelic_nomiss.samples.removed.vcf.gz")
por_vcf_tbi = file("/home/humebc/projects/paper_4/testing_flink/POR/Porites_lobata_v3_11Islands_New_maf05_minQ30_biallelic_nomiss.samples.removed.vcf.gz.tbi")

// Make a channel that is the scaffold values
// Then add a species denotion to this
poc_scaffs_list_path = file('/home/humebc/projects/paper_4/testing_flink/POC.scaff.10000.list.txt')
poc_scaffs_list  = poc_scaffs_list_path.readLines()
poc_scaffs_ch = Channel.fromList(poc_scaffs_list)
poc_spec_ch = Channel.from( ["POC"] )
poc_scaffs_spec_ch = poc_spec_ch.combine(poc_scaffs_ch)

por_scaffs_list_path = file('/home/humebc/projects/paper_4/testing_flink/POR.scaff.10000.list.txt')
por_scaffs_list  = por_scaffs_list_path.readLines()
por_scaffs_ch = Channel.fromList(por_scaffs_list)
por_spec_ch = Channel.from( ["POR"] )
por_scaffs_spec_ch = por_spec_ch.combine(por_scaffs_ch)


process subset_vcf {
    tag "subset_vcf: ${scaffold}"
    container 'halllab/bcftools:v1.9'

    input:
    tuple val(species), val(scaffold) from poc_scaffs_spec_ch.concat(por_scaffs_spec_ch)
    file poc_vcf
    file poc_vcf_tbi
    file por_vcf
    file por_vcf_tbi

    output:
    tuple val(species), val(scaffold), file("*.${scaffold}.vcf.gz") into flink_atlas_make_input_ch,flink_atlas_make_input_lineage_as_groups_ch

    script:
    if (species == "POC")
    """
    bcftools view -O z -o POC.${scaffold}.vcf.gz ${poc_vcf} ${scaffold}
    """
    else if (species == "POR")
    """
    bcftools view -O z -o POR.${scaffold}.vcf.gz ${por_vcf} ${scaffold}
    """

}

// We will run two sets of flink
// One which has one group with multiple populations
// where each population is one of our detected resolved lineages
// We will call this: flink_atlas_lineage_as_populations
// Then we will do one where we have multiple groups each containing multiple populations
// In this scenario, the groups will be the resolved lineages and the populations will be the islands.
// We will call this: flink_atlas_lineage_as_groups
process flink_atlas {
    tag "flink_atlas: ${scaffold}"
    container "didillysquat/flink_atlas:latest"
    publishDir "/home/humebc/projects/paper_4/testing_flink/flink_results/lineages_as_populations/${species}/${scaffold}"
    cpus 4

    input:
    tuple val(species), val(scaffold), file(subset_vcf) from flink_atlas_make_input_ch
    file poc_sample_populations
    file por_sample_populations

    output:
    file("*") into flink_atlas_ch_out
    tuple val(species), val(scaffold), file("*.in_MCMCM_hierarchicalParameters.txt") into calc_av_params_lineage_as_populations_ch

    shell:
    if (species == "POC")
    '''
    atlas task=alleleCounts vcf=!{subset_vcf} samples=!{poc_sample_populations}
    zcat *_alleleCounts.txt.gz | awk 'BEGIN{print "-\t-\tGroup1\tGroup1\tGroup1\tGroup1\tGroup1"} {if (match($0, /^chr/)){print "-\t-\tSVD" $3 "\tSVD" $4 "\tSVD" $5 "\tSVD" $6 "\tSVD" $7}else{print $0}}' - > !{species}.!{scaffold}.flink.in.txt
    flink task=estimate numThreads=!{task.cpus} lnMu_g="(-4.0,-0.0)" lnNu_g="(-5.0,-0.0)" s_max=14 beta="(-2.0,1.8)" alpha_max=4.0 numIterations=500000 burnin=300000 thinning=100 lnkappa="(-10.0,-0.1)" logFile=POC.logfile sigmaProp_mu=0.005 sigmaProp_nu=0.05 sigmaProp_kappa=0.05 data=!{species}.!{scaffold}.flink.in.txt
    '''
    else if (species == "POR")
    '''
    atlas task=alleleCounts vcf=!{subset_vcf} samples=!{por_sample_populations}
    zcat *_alleleCounts.txt.gz | awk 'BEGIN{print "-\t-\tGroup1\tGroup1\tGroup1"} {if (match($0, /^chr/)){print "-\t-\tK" $3 "\tK" $4 "\tK" $5}else{print $0}}' - > !{species}.!{scaffold}.flink.in.txt
    flink task=estimate numThreads=!{task.cpus} lnMu_g="(-4.0,-0.0)" lnNu_g="(-5.0,-0.0)" s_max=14 beta="(-2.0,1.8)" alpha_max=4.0 numIterations=500000 burnin=300000 thinning=100 lnkappa="(-10.0,-0.1)" logFile=POC.logfile sigmaProp_mu=0.005 sigmaProp_nu=0.05 sigmaProp_kappa=0.05 data=!{species}.!{scaffold}.flink.in.txt
    '''
}

process flink_atlas_lineage_as_groups {
    tag "flink_atlas: ${scaffold}"
    container "didillysquat/flink_atlas:latest"
    publishDir "/home/humebc/projects/paper_4/testing_flink/flink_results/lineages_as_groups/${species}/${scaffold}"
    cpus 4

    input:
    tuple val(species), val(scaffold), file(subset_vcf) from flink_atlas_make_input_lineage_as_groups_ch
    file poc_sample_populations_lineage_island
    file por_sample_populations_lineage_island
    file awk_create_flink_input_lineage_island

    output:
    file("*") into flink_atlas_ch_lineage_as_groups_out
    tuple val(species), val(scaffold), file("*.in_MCMCM_hierarchicalParameters.txt") into calc_av_params_lineage_as_groups_ch

    shell:
    if (species == "POC")
    '''
    atlas task=alleleCounts vcf=!{subset_vcf} samples=!{poc_sample_populations_lineage_island}
    zcat *_alleleCounts.txt.gz | awk -f !{awk_create_flink_input_lineage_island} > !{species}.!{scaffold}.flink.in.txt
    flink task=estimate numThreads=!{task.cpus} A_max=4.0 B="(-2.0,1.8)" lnK=-10.0,-0.1 lnMu=-4.0,0.0 lnNu=-5.0,0.0 \
    lnMu_g="(-4.0,-0.0)" lnNu_g="(-5.0,-0.0)" s_max=14 beta="(-2.0,1.8)" alpha_max=4.0 numIterations=500000 burnin=300000 \
    thinning=100 lnkappa="(-10.0,-0.1)" logFile=POC.logfile sigmaProp_mu=0.005 sigmaProp_nu=0.05 sigmaProp_kappa=0.05 \
    data=!{species}.!{scaffold}.flink.in.txt
    '''
    else if (species == "POR")
    '''
    atlas task=alleleCounts vcf=!{subset_vcf} samples=!{por_sample_populations_lineage_island}
    zcat *_alleleCounts.txt.gz | awk -f !{awk_create_flink_input_lineage_island} > !{species}.!{scaffold}.flink.in.txt
    flink task=estimate numThreads=!{task.cpus} A_max=4.0 B="(-2.0,1.8)" lnK=-10.0,-0.1 lnMu=-4.0,0.0 lnNu=-5.0,0.0 \
    lnMu_g="(-4.0,-0.0)" lnNu_g="(-5.0,-0.0)" s_max=14 beta="(-2.0,1.8)" alpha_max=4.0 numIterations=500000 burnin=300000 \
    thinning=100 lnkappa="(-10.0,-0.1)" logFile=POR.logfile sigmaProp_mu=0.005 sigmaProp_nu=0.05 sigmaProp_kappa=0.05 \
    data=!{species}.!{scaffold}.flink.in.txt
    '''
}

// We also want to work out whether the same values are being converged to so that we can assess the possibility
// of fixing the parameter values so that we can predict the run selction profiles for the shorter
// scaffolds in addition to the longer scaffolds.
// To do this we will do an additional process where we calculate the parameter averages
process calc_param_averages_lineage_as_population{
    tag "calc_param_averages_lineage_as_population: ${scaffold}"
    container "didillysquat/r_multipurpose:latest"
    publishDir "/home/humebc/projects/paper_4/testing_flink/flink_results/param_averages_lineages_as_populations"
    cpus 1

    input:
    tuple val(species), val(scaffold), file(mcmc) from calc_av_params_lineage_as_populations_ch
    file calculate_param_averages_lineage_as_populations

    output:
    file "*.averages.tsv"  into calc_param_averages_lineage_as_population_out_ch

    script:
    """
    Rscript $calculate_param_averages_lineage_as_populations $species $scaffold $mcmc
    """
}

process calc_param_averages_lineage_as_group{
    tag "calc_param_averages_lineage_as_population: ${scaffold}"
    container "didillysquat/r_multipurpose:latest"
    publishDir "/home/humebc/projects/paper_4/testing_flink/flink_results/param_averages_lineages_as_groups"
    cpus 1

    input:
    tuple val(species), val(scaffold), file(mcmc) from calc_av_params_lineage_as_groups_ch
    file calculate_param_averages_lineage_as_populations

    output:
    file "*.averages.tsv"  into calc_param_averages_lineage_as_group_out_ch

    script:
    """
    Rscript $calculate_param_averages_lineage_as_populations $species $scaffold $mcmc
    """
}