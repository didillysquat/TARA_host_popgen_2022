#!/usr/bin/env nextflow

/*
The first part of the flink analysis was done in flink.nf
In that analysis we were estimating all parameters for those scaffholds that had more than 10000 variants in them
We did not estimate parameters for those scaffolds with less than that number of variants because the estimation is
unreliable when the number of variable sites is low.
From these results we have worked out a set of parameters that we can use for all remaining scaffolds by fixing the parameters
to those values.

There is one special case which is scaffold Pocillopora_meandrina_contig_1. This scaffold did not complete computation for the
lineages as groups structure
before we had to restart the server. However we did get through 700 000 iterations out of 800 000 so we will consider it complete
and work with its output to calculate the selection profiles just like the other scaffolds.

So for this script we will be working with the remaining scaffolds and running the same analysis as before but with fixed parameters.
The complication for this is that we will need to rearrange the columns of the flink input so that we can match them up with the
fixed parameters. This is probably best done with a python script.

NB, we have excluded those scaffolds with <=100 variants due to lack of meaningful information. in POR this represents about 3% of the genome
unassessed. in POC it is << 1%.
*/

poc_sample_populations = file("/home/humebc/projects/paper_4/testing_flink/POC/POC.samplesPopulations.txt")
por_sample_populations = file("/home/humebc/projects/paper_4/testing_flink/POR/POR.samplesPopulations.txt")
poc_sample_populations_lineage_island = file("/home/humebc/projects/paper_4/testing_flink/POC/POC.samplesPopulations.lineage.island.txt")
por_sample_populations_lineage_island = file("/home/humebc/projects/paper_4/testing_flink/POR/POR.samplesPopulations.lineage.island.txt")

calculate_param_averages_lineage_as_populations = file("/home/humebc/projects/paper_4/testing_flink/bin/calc_param_averages_lineage_as_populations.r")

awk_create_flink_input_lineage_island = file("/home/humebc/projects/paper_4/testing_flink/bin/create_flink_input_lineage_island.awk")
python_create_flink_input_from_atlas_output = file("/home/humebc/projects/paper_4/testing_flink/bin/modify_atlas_output.py")

// For each of the scaffolds to process we will subset the vcf using bcftools
// We will then create the input files for each of these using Atlas and start a flink analysis for each

poc_vcf = file("/home/humebc/projects/paper_4/testing_flink/POC/Pocillopora_meandrina_v3_11Islands_New_maf05_minQ30_biallelic_nomiss.samples.removed.vcf.gz")
poc_vcf_tbi = file("/home/humebc/projects/paper_4/testing_flink/POC/Pocillopora_meandrina_v3_11Islands_New_maf05_minQ30_biallelic_nomiss.samples.removed.vcf.gz.tbi")
por_vcf = file("/home/humebc/projects/paper_4/testing_flink/POR/Porites_lobata_v3_11Islands_New_maf05_minQ30_biallelic_nomiss.samples.removed.vcf.gz")
por_vcf_tbi = file("/home/humebc/projects/paper_4/testing_flink/POR/Porites_lobata_v3_11Islands_New_maf05_minQ30_biallelic_nomiss.samples.removed.vcf.gz.tbi")

// Make a channel that is the scaffold values
// Then add a species denotion to this
poc_scaffs_list_path = file('/home/humebc/projects/paper_4/testing_flink/POC.scaff.less.than.10000.grt.thn.100.list.txt')
poc_scaffs_list  = poc_scaffs_list_path.readLines()
poc_scaffs_ch = Channel.fromList(poc_scaffs_list)
poc_spec_ch = Channel.from( ["POC"] )
poc_scaffs_spec_ch = poc_spec_ch.combine(poc_scaffs_ch)

por_scaffs_list_path = file('/home/humebc/projects/paper_4/testing_flink/POR.scaff.less.than.10000.grt.thn.100.list.txt')
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
process atlas_lineage_as_population {
    tag "atlas_lineage_as_population: ${species} ${scaffold}"
    container "didillysquat/flink_atlas:latest"
    cpus 1

    input:
    tuple val(species), val(scaffold), file(subset_vcf) from flink_atlas_make_input_ch
    file poc_sample_populations
    file por_sample_populations

    output:
    tuple val(species), val(scaffold), file("*alleleCounts.txt.gz") into flink_in_as_population_ch
    

    script:
    if (species == "POC")
    """
    atlas task=alleleCounts vcf=${subset_vcf} samples=${poc_sample_populations}
    """
    else if (species == "POR")
    """
    atlas task=alleleCounts vcf=${subset_vcf} samples=${por_sample_populations}
    """
}

process flink_in_as_population{
    tag "atlas_lineage_as_population: ${species} ${scaffold}"
    conda "/home/humebc/projects/paper_4/testing_flink/bin/env/general_python3.yml"
    publishDir "/home/humebc/projects/paper_4/testing_flink/flink_fixed_parameters/flink_results/lineages_as_populations/${species}/${scaffold}"
    cpus 1

    input:
    tuple val(species), val(scaffold), file(alleleCounts) from flink_in_as_population_ch
    file python_create_flink_input_from_atlas_output

    output:
    tuple val(species), val(scaffold), file("${species}.${scaffold}.flink.in.as.populations.txt") into flink_lineage_as_populations_ch

    script:
    """
    python3 $python_create_flink_input_from_atlas_output $species $scaffold population $alleleCounts
    """
}

process atlas_lineage_as_group {
    tag "atlas_lineage_as_groups: ${species} ${scaffold}"
    container "didillysquat/flink_atlas:latest"
    cpus 1

    input:
    tuple val(species), val(scaffold), file(subset_vcf) from flink_atlas_make_input_lineage_as_groups_ch
    file poc_sample_populations_lineage_island
    file por_sample_populations_lineage_island

    output:
    tuple val(species), val(scaffold), file("*alleleCounts.txt.gz") into flink_in_as_group_ch

    script:
    if (species == "POC")
    """
    atlas task=alleleCounts vcf=${subset_vcf} samples=${poc_sample_populations_lineage_island}
    """
    else if (species == "POR")
    """
    atlas task=alleleCounts vcf=${subset_vcf} samples=${por_sample_populations_lineage_island}
    """
}

process flink_in_as_group{
    tag "flink_in_as_group: ${species} ${scaffold}"
    conda "/home/humebc/projects/paper_4/testing_flink/bin/env/general_python3.yml"
    publishDir "/home/humebc/projects/paper_4/testing_flink/flink_fixed_parameters/flink_results/lineages_as_groups/${species}/${scaffold}"
    cpus 1

    input:
    tuple val(species), val(scaffold), file(alleleCounts) from flink_in_as_group_ch
    file python_create_flink_input_from_atlas_output

    output:
    tuple val(species), val(scaffold), file("${species}.${scaffold}.flink.in.as.groups.txt") into flink_lineage_as_groups_ch

    script:
    """
    python3 $python_create_flink_input_from_atlas_output $species $scaffold group $alleleCounts
    """
}

process flink_lineage_as_populations {
    tag "flink_lineage_as_population: ${species} ${scaffold}"
    container "didillysquat/flink_atlas:latest"
    publishDir "/home/humebc/projects/paper_4/testing_flink/flink_fixed_parameters/flink_results/lineages_as_populations/${species}/${scaffold}"
    cpus 4

    input:
    tuple val(species), val(scaffold), file(flink_in) from flink_lineage_as_populations_ch

    output:
    file("*") into flink_lineage_as_populations_ch_out

    script:
    if (species == "POC")
    """
    flink task=estimate numThreads=${task.cpus} lnkappa=-6.11365784 lnMu_g=-0.12166096 lnNu_g=-2.58021025 beta=[2.29458681,0.01010877,0.06741716,-0.05363849,0.11839443] s_max=14 alpha_max=4.0 numIterations=500000 burnin=300000 thinning=100 logFile=POC.as.population.logfile sigmaProp_mu=0.005 sigmaProp_nu=0.05 sigmaProp_kappa=0.05 data=${flink_in}
    """
    else if (species == "POR")
    """
    flink task=estimate numThreads=${task.cpus} lnkappa=-8.0292037 lnMu_g=-0.5415207 lnNu_g=-0.7368341 beta=[0.8124666,0.1223843,-0.5396381] s_max=14 alpha_max=4.0 numIterations=500000 burnin=300000 thinning=100 logFile=POC.as.population.logfile sigmaProp_mu=0.005 sigmaProp_nu=0.05 sigmaProp_kappa=0.05 data=${flink_in}
    """
}

process flink_lineage_as_groups {
    tag "flink_lineage_as_population: ${species} ${scaffold}"
    container "didillysquat/flink_atlas:latest"
    publishDir "/home/humebc/projects/paper_4/testing_flink/flink_fixed_parameters/flink_results/lineages_as_groups/${species}/${scaffold}"
    cpus 4

    input:
    tuple val(species), val(scaffold), file(flink_in) from flink_lineage_as_groups_ch
    
    output:
    file("*") into flink_lineage_as_groups_ch_out

    script:
    if (species == "POC")
    """
    # TODO fix parameters
    flink task=estimate numThreads=${task.cpus} lnK=-7.693636 lnMu=-0.4685393 lnNu=-0.01341808 \
    B=0.1657752,-0.5836575,-0.4376080,-0.5607286,-0.4706092 lnkappa=-3.307948,-3.435418,-4.020781,-4.501050,-5.108612 lnMu_g=-0.0005084079,-0.04975536,-0.005003625,-0.07767124,-0.1316791 lnNu_g=-2.021511,-0.1651979,-1.953806,-0.02662667,-0.02286499 \
    beta=[-1.468845,-1.497420,-1.597289],[-3.501289,-11.68201,-12.17739,-4.984479,-8.927332],[-11.06781,-10.84990,-11.96515,-4.862533],[-4.661715,-5.275214,-3.944286,-4.621625,-6.787444],[-3.026023,-10.99860,-4.225887] \
    A_max=4.0 alpha_max=4.0 s_max=14 numIterations=500000 burnin=300000 \
    thinning=100  logFile=POC.as.group.logfile sigmaProp_mu=0.005 sigmaProp_nu=0.05 sigmaProp_kappa=0.05 \
    data=${flink_in}
    """
    else if (species == "POR")
    """
    flink task=estimate numThreads=${task.cpus} lnK=-9.69016550 lnMu=-1.52804563 lnNu=-0.96370943 \
    B=0.24192145,-0.10296206,-0.70763183, lnkappa=-5.27544072,-5.87896448,-0.64476588 lnMu_g=-0.16881391,-0.29761149,-0.03492038 lnNu_g=-0.11705645,-0.10152420,-2.55296037 \
    beta=[-3.67682280,-9.74113019,-2.23321272],[-2.23966769,-1.93963595,-1.92797906,-1.98811786,-1.82580770,-0.07664835,-1.28976588,-3.30993384],[-10.26140002,-1.67004772,-1.67371249,-2.66060952,-2.51812963,-2.76303691,-1.69630430,-1.69652767] \
    A_max=4.0 alpha_max=4.0 s_max=14 numIterations=500000 burnin=300000 \
    thinning=100  logFile=POR.as.group.logfile sigmaProp_mu=0.005 sigmaProp_nu=0.05 sigmaProp_kappa=0.05 \
    data=${flink_in}
    """
}