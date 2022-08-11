"""
Here we compare the selection profiles that we have created with flink to the SNPs that have been identified by Didier as 
linked to historical climate paramters and genomic outliers.

We will aim to simply append the q-values for the alphas for divergent and balancing selection
from each of the groups results onto the dataframe from Didier and then we can work with this.

One issue is that the following 2 scaffolds in POR have < 100 snps on them and so weren't included in the selection analysis
> scaff.df[scaff.df$scaff == "Porites_lobata_contig_848",]
                        scaff count length
859 Porites_lobata_contig_848    28  40076

> scaff.df[scaff.df$scaff == "Porites_lobata_contig_918",]
                        scaff count length
928 Porites_lobata_contig_918    98  27596

We will have to manually account for this.
"""

import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

class ScreenSelectionProfiles:
    def __init__(self):
        self.original_poc_df_path = "/home/humebc/projects/paper_4/testing_flink/selection_table_results/PocG111_100_outlier_HistPCA_pval_qval.tsv"
        self.new_poc_out_df_path = "/home/humebc/projects/paper_4/testing_flink/selection_table_results/PocG111_100_outlier_HistPCA_pval_qval.w_selection_q_vals.tsv"
        self.original_por_df_path = "/home/humebc/projects/paper_4/testing_flink/selection_table_results/PorG111_100_outlier_HistPCA_pval_qval.tsv"
        self.new_por_out_df_path = "/home/humebc/projects/paper_4/testing_flink/selection_table_results/PorG111_100_outlier_HistPCA_pval_qval.w_selection_q_vals.tsv"
        
        self.alpha = 0.05
        self.svd_headers = [
            "lin_as_grp_divergent_qval_SVD1", "lin_as_grp_balancing_qval_SVD1", "lin_as_grp_divergent_qval_SVD2", "lin_as_grp_balancing_qval_SVD2",
            "lin_as_grp_divergent_qval_SVD3", "lin_as_grp_balancing_qval_SVD3", "lin_as_grp_divergent_qval_SVD4", "lin_as_grp_balancing_qval_SVD4",
            "lin_as_grp_divergent_qval_SVD5", "lin_as_grp_balancing_qval_SVD5"
            ]

        self.svd_headers_grp_divergent_only = [
        "lin_as_grp_divergent_qval_SVD1", "lin_as_grp_divergent_qval_SVD2",
        "lin_as_grp_divergent_qval_SVD3", "lin_as_grp_divergent_qval_SVD4",
        "lin_as_grp_divergent_qval_SVD5"
        ]

        self.svd_headers_grp_hh_divergent_only = [
        "lin_as_grp_divergent_qval_SVD1", "lin_as_grp_divergent_qval_SVD2",
        "lin_as_grp_divergent_qval_SVD3", "lin_as_grp_divergent_qval_SVD4",
        "lin_as_grp_divergent_qval_SVD5", "lin_as_grp_divergent_qval_higher_hierarchy"
        ]

        self.k_headers = [
            "lin_as_grp_divergent_qval_K1", "lin_as_grp_balancing_qval_K1", "lin_as_grp_divergent_qval_K2", "lin_as_grp_balancing_qval_K2",
            "lin_as_grp_divergent_qval_K3", "lin_as_grp_balancing_qval_K3"
            ]

        self.k_headers_grp_divergent_only = [
            "lin_as_grp_divergent_qval_K1", "lin_as_grp_divergent_qval_K2",
            "lin_as_grp_divergent_qval_K3",
            ]

        self.k_headers_grp_hh_divergent_only = [
            "lin_as_grp_divergent_qval_K1", "lin_as_grp_divergent_qval_K2",
            "lin_as_grp_divergent_qval_K3", "lin_as_grp_divergent_qval_higher_hierarchy"
            ]

        try:
            # Try to read in the df to which we have added the selection profile q values
            self.poc_df = pd.read_table(self.new_poc_out_df_path)
            self.por_df = pd.read_table(self.new_por_out_df_path)
        except:
            # Else make the new table
            self.poc_df = pd.read_table(self.original_poc_df_path)
            self.poc_df.index = [f"{chrom}:{pos}" for chrom, pos in zip(self.poc_df["CHROM"],self.poc_df["POS"])]
            self.por_df = pd.read_table(self.original_por_df_path)
            self.por_df.index = [f"{chrom}:{pos}" for chrom, pos in zip(self.por_df["CHROM"],self.por_df["POS"])]
            self.poc_df, self.por_df = self._make_didier_snp_tables()

        self.selection_table_base_dir = "/home/humebc/projects/paper_4/testing_flink/selection_table_results/"

    def _questions(self):
        """ Method for general questions that we want to answer."""
        # Update 20220414
        # After a chat with dider we decided it would be good to leave out the balancing selection
        # and rather focus on the divergent seleciton at the group level and possibly at the higher hierarchy as well
        
        # The first question we want to answer is how the proportion of the number of SNPs under seleciton 
        # across the whole genome compares to the selection of SNPs that didier identified as being linked to
        # the climate variables.
        # We can test statistically be bootstrapping 100 random SNP sites.
        self._output_stats_for_species(species="POC", species_df=self.poc_df)
        self._output_stats_for_species(species="POR", species_df=self.por_df)
        
        foo = "bar"

    def _output_stats_for_species(self, species_df, species="POC"):
        underselection_grp_only = []
        underselection_grp_only_genom_island = []
        underselection_frequencies_grp_only = []
        underselection_frequencies_grp_only_genom_island = []
        underselection_grp_or_hh = []
        underselection_grp_or_hh_genom_island = []
        underselection_grp_and_hh = []
        underselection_grp_and_hh_genom_island = []
        underselection_frequencies_grp_or_hh = []
        underselection_frequencies_grp_or_hh_genom_island = []
        underselection_hh_only = []
        underselection_hh_only_genom_island = []
        genom_island_snps = species_df.loc[species_df["in genomic island"] == 1,].index
        for snp in species_df.index:
            # screen the values to see if any are less than alpha
            if species == "POC":
                selection_values_grp_only = species_df.loc[snp, self.svd_headers_grp_divergent_only]
                selection_values_grp_or_hh = species_df.loc[snp, self.svd_headers_grp_hh_divergent_only]
            elif species == "POR":
                selection_values_grp_only = species_df.loc[snp, self.k_headers_grp_divergent_only]
                selection_values_grp_or_hh = species_df.loc[snp, self.k_headers_grp_hh_divergent_only]
            
            
            selection_count_grp_only = sum(selection_values_grp_only <= self.alpha)
            underselection_frequencies_grp_only.append(selection_count_grp_only)
            if snp in genom_island_snps:
                underselection_frequencies_grp_only_genom_island.append(selection_count_grp_only)
                if selection_count_grp_only > 0:
                    underselection_grp_only_genom_island.append(True)
                else:
                    underselection_grp_only_genom_island.append(False)
            if selection_count_grp_only > 0:
                underselection_grp_only.append(True)
            else:
                underselection_grp_only.append(False)


            selection_count_grp_and_hh = sum(selection_values_grp_or_hh <= self.alpha)
            underselection_frequencies_grp_or_hh.append(selection_count_grp_and_hh)
            if snp in genom_island_snps:
                underselection_frequencies_grp_or_hh_genom_island.append(selection_count_grp_and_hh)
                if selection_count_grp_and_hh > 0:
                    underselection_grp_or_hh_genom_island.append(True)
                else:
                    underselection_grp_or_hh_genom_island.append(False)
            if selection_count_grp_and_hh > 0:
                underselection_grp_or_hh.append(True)
            else:
                underselection_grp_or_hh.append(False)


            selection_value_hh_only = species_df.at[snp, "lin_as_grp_divergent_qval_higher_hierarchy"]
            if snp in genom_island_snps:
                if selection_value_hh_only < self.alpha:
                    underselection_hh_only_genom_island.append(True)
                    if sum(selection_values_grp_only <= self.alpha) > 1:
                        underselection_grp_and_hh_genom_island.append(True)
                    else:
                        underselection_grp_and_hh_genom_island.append(False)
                else:
                    underselection_hh_only_genom_island.append(False)
                    underselection_grp_and_hh_genom_island.append(False)
            if selection_value_hh_only < self.alpha:
                underselection_hh_only.append(True)
                if sum(selection_values_grp_only <= self.alpha) > 1:
                    underselection_grp_and_hh.append(True)
                else:
                    underselection_grp_and_hh.append(False)
            else:
                underselection_hh_only.append(False)
                underselection_grp_and_hh.append(False)


        species_df["under_selection_grp_only"] = underselection_grp_only
        species_df["under_selection_grp_and_hh"] = underselection_grp_or_hh
        species_df["under_selection_hh_only"] = underselection_hh_only
        
        climate_linked_prop = sum(underselection_grp_only)/len(underselection_grp_only)
        
        # We decided to work with the higher hierarchy only in the end so I have commented out the others.
        # print(f"\n\nproportion of {species} SNPs linked to climate history under selection (groups only): {sum(underselection_grp_only)/len(underselection_grp_only)}")
        # print(f"proportion of {species} SNPs linked to climate history under selection (at least 1 group OR higher hierarchy): {sum(underselection_grp_or_hh)/len(underselection_grp_or_hh)}")
        # print(f"proportion of {species} SNPs linked to climate history under selection (at least 1 group AND higher hierarchy): {sum(underselection_grp_and_hh)/len(underselection_grp_and_hh)}")
        print(f"proportion of {species} SNPs linked to climate history under selection (higher hierarchy only): {sum(underselection_hh_only)}/{len(underselection_hh_only)}={sum(underselection_hh_only)/len(underselection_hh_only)}")
        
        # print(f"\nproportion of {species} SNPs linked to climate history and in genomic island under selection (groups only): {sum(underselection_grp_only_genom_island)}/{len(underselection_grp_only_genom_island)}={sum(underselection_grp_only_genom_island)/len(underselection_grp_only_genom_island)}")
        # print(f"proportion of {species} SNPs linked to climate history and in genomic island under selection (at least 1 group OR higher hierarchy): {sum(underselection_grp_or_hh_genom_island)}/{len(underselection_grp_or_hh_genom_island)}={sum(underselection_grp_or_hh_genom_island)/len(underselection_grp_or_hh_genom_island)}")
        # print(f"proportion of {species} SNPs linked to climate history and in genomic island under selection (at least 1 group AND higher hierarchy): {sum(underselection_grp_and_hh_genom_island)}/{len(underselection_grp_and_hh_genom_island)}={sum(underselection_grp_and_hh_genom_island)/len(underselection_grp_and_hh_genom_island)}")
        print(f"proportion of {species} SNPs linked to climate history and in genomic island under selection (higher hierarchy only): {sum(underselection_hh_only_genom_island)}/{len(underselection_hh_only_genom_island)}={sum(underselection_hh_only_genom_island)/len(underselection_hh_only_genom_island)}")

        fig, ax = plt.subplots(ncols=2, nrows=1)
        if species == "POC":
            bins = [-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5]
            x_ticks = [0,1,2,3,4,5]
        elif species == "POR":
            bins = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5]
            x_ticks = [0,1,2,3]
        
        ax[0].hist([_ - 0.5 for _ in underselection_frequencies_grp_only], bins=bins)
        ax[0].hist([_ - 0.5 for _ in underselection_frequencies_grp_only_genom_island], bins=bins)
        ax[0].set_xlabel("number of lineages in which SNP is under\ndivergent selection (groups only)", fontsize=8)
        ax[0].set_xticks(x_ticks)
        ax[0].set_xticklabels(x_ticks)
        ax[0].set_ylabel("number of SNPs")
        # ax[0].set_title("Number of POC lineages in which the climate linked loci are under selection")

        ax[1].hist([_ - 0.5 for _ in underselection_frequencies_grp_or_hh], bins=bins)
        ax[1].hist([_ - 0.5 for _ in underselection_frequencies_grp_or_hh_genom_island], bins=bins)
        ax[1].set_xlabel("number of lineages in which SNP is under\ndivergent selection (groups and HigherHier)", fontsize=8)
        ax[1].set_xticks(x_ticks)
        ax[1].set_xticklabels(x_ticks)
        ax[1].set_ylabel("number of SNPs")
        fig.suptitle(f"Number of {species} lineages in which the climate linked loci are under selection")
        plt.tight_layout()
        plt.savefig(f"/home/humebc/projects/paper_4/testing_flink/selection_table_results/{species}_RDA_SNPs_hist_divergent_only.png", dpi=600)

        # Now we want to know how this compares to the proportion of the SNPs over the whole genoe that are under selection
        # We will have to get the average over all of the SVD groups. THis measn reading in each of the dfs
        fig, ax = plt.subplots()
        tot = 0
        snp_under_selection_set_grp_only = set()
        snp_under_selection_set_grp_or_hh = set()
        snp_under_selection_set_grp_and_hh = set()
        snp_under_selection_set_hh_only = set()
        if species == "POC":
            lineages = [1,2,3,4,5]
            file_base = "poc_lineage_as_group_SVD"
            hh_spec = "poc"
        elif species == "POR":
            lineages = [1,2,3]
            file_base = "por_lineage_as_group_K"
            hh_spec = "por"
        for lin in lineages:
            df = pd.read_table(os.path.join(self.selection_table_base_dir, f"{file_base}{lin}_df.tsv"))
            snps_under_div = df.loc[df["divergent"] < self.alpha,].index.values
            snp_under_selection_set_grp_only.update(snps_under_div)
            snp_under_selection_set_grp_or_hh.update(snps_under_div)
            tot = df.shape[0]
        
        # and the hh
        df = pd.read_table(f"/home/humebc/projects/paper_4/testing_flink/selection_table_results/{hh_spec}_lineage_as_group_higher_hierarchy_df.tsv")
        snps_under_div = df.loc[df["divergent"] < self.alpha,].index.values
        snp_under_selection_set_grp_and_hh.update(set(snps_under_div).intersection(snp_under_selection_set_grp_only))
        snp_under_selection_set_grp_or_hh.update(snps_under_div)
        snp_under_selection_set_hh_only.update(snps_under_div)
        
        # plt.savefig("/home/humebc/projects/paper_4/testing_flink/selection_table_results/POC_selection_proportion.png", dpi=600)
        
        # change this so only div
        # print(f"\nProportion of SNPs under div selection in at least one group across the whole genome in {species}: {len(snp_under_selection_set_grp_only)/tot}")
        # print(f"Proportion of SNPs under div selection in at least one group OR HH across the whole genome in {species}: {len(snp_under_selection_set_grp_or_hh)/tot}")
        # print(f"Proportion of SNPs under div selection in at least one group AND HH across the whole genome in {species}: {len(snp_under_selection_set_grp_and_hh)/tot}")
        print(f"Proportion of SNPs under div selection in the HH across the whole genome in {species}: {len(snp_under_selection_set_hh_only)/tot}")

    def _make_didier_snp_tables(self):

        # We will append the qvalues in order of lineage as pop and then lineage as group
        ## POC
        # lineage as population
        print("poc lin as pop")
        poc_lin_as_pop_df = pd.read_table("/home/humebc/projects/paper_4/testing_flink/selection_table_results/poc_lineage_as_population_0.tsv")
        poc_lin_as_pop_df.index = [f"Pocillopora_meandrina_contig_{chrom}:{pos}" for chrom, pos in zip(poc_lin_as_pop_df["scaff_coord"],poc_lin_as_pop_df["bp_coord"])]
        self.poc_df["lin_as_pop_divergent_qval"] = poc_lin_as_pop_df.loc[self.poc_df.index,"divergent"]
        self.poc_df["lin_as_pop_balancing_qval"] = poc_lin_as_pop_df.loc[self.poc_df.index,"balancing"]

        # lineage as group SVD1
        print("poc lin as group svd1")
        poc_lin_as_grp_svd1_df = pd.read_table("/home/humebc/projects/paper_4/testing_flink/selection_table_results/poc_lineage_as_group_SVD1_df.tsv")
        poc_lin_as_grp_svd1_df.index = [f"Pocillopora_meandrina_contig_{chrom}:{pos}" for chrom, pos in zip(poc_lin_as_grp_svd1_df["scaff_coord"],poc_lin_as_grp_svd1_df["bp_coord"])]
        self.poc_df["lin_as_grp_divergent_qval_SVD1"] = poc_lin_as_grp_svd1_df.loc[self.poc_df.index,"divergent"]
        self.poc_df["lin_as_grp_balancing_qval_SVD1"] = poc_lin_as_grp_svd1_df.loc[self.poc_df.index,"balancing"]

        # lineage as group SVD2
        print("poc lin as group svd2")
        poc_lin_as_grp_svd2_df = pd.read_table("/home/humebc/projects/paper_4/testing_flink/selection_table_results/poc_lineage_as_group_SVD2_df.tsv")
        poc_lin_as_grp_svd2_df.index = [f"Pocillopora_meandrina_contig_{chrom}:{pos}" for chrom, pos in zip(poc_lin_as_grp_svd2_df["scaff_coord"],poc_lin_as_grp_svd2_df["bp_coord"])]
        self.poc_df["lin_as_grp_divergent_qval_SVD2"] = poc_lin_as_grp_svd2_df.loc[self.poc_df.index,"divergent"]
        self.poc_df["lin_as_grp_balancing_qval_SVD2"] = poc_lin_as_grp_svd2_df.loc[self.poc_df.index,"balancing"]

        # lineage as group SVD3
        print("poc lin as group svd3")
        poc_lin_as_grp_svd3_df = pd.read_table("/home/humebc/projects/paper_4/testing_flink/selection_table_results/poc_lineage_as_group_SVD3_df.tsv")
        poc_lin_as_grp_svd3_df.index = [f"Pocillopora_meandrina_contig_{chrom}:{pos}" for chrom, pos in zip(poc_lin_as_grp_svd3_df["scaff_coord"],poc_lin_as_grp_svd3_df["bp_coord"])]
        self.poc_df["lin_as_grp_divergent_qval_SVD3"] = poc_lin_as_grp_svd3_df.loc[self.poc_df.index,"divergent"]
        self.poc_df["lin_as_grp_balancing_qval_SVD3"] = poc_lin_as_grp_svd3_df.loc[self.poc_df.index,"balancing"]

        # lineage as group SVD4
        print("poc lin as group svd4")
        poc_lin_as_grp_svd4_df = pd.read_table("/home/humebc/projects/paper_4/testing_flink/selection_table_results/poc_lineage_as_group_SVD4_df.tsv")
        poc_lin_as_grp_svd4_df.index = [f"Pocillopora_meandrina_contig_{chrom}:{pos}" for chrom, pos in zip(poc_lin_as_grp_svd4_df["scaff_coord"],poc_lin_as_grp_svd4_df["bp_coord"])]
        self.poc_df["lin_as_grp_divergent_qval_SVD4"] = poc_lin_as_grp_svd4_df.loc[self.poc_df.index,"divergent"]
        self.poc_df["lin_as_grp_balancing_qval_SVD4"] = poc_lin_as_grp_svd4_df.loc[self.poc_df.index,"balancing"]

        # lineage as group SVD5
        print("poc lin as group svd5")
        poc_lin_as_grp_svd5_df = pd.read_table("/home/humebc/projects/paper_4/testing_flink/selection_table_results/poc_lineage_as_group_SVD5_df.tsv")
        poc_lin_as_grp_svd5_df.index = [f"Pocillopora_meandrina_contig_{chrom}:{pos}" for chrom, pos in zip(poc_lin_as_grp_svd5_df["scaff_coord"],poc_lin_as_grp_svd5_df["bp_coord"])]
        self.poc_df["lin_as_grp_divergent_qval_SVD5"] = poc_lin_as_grp_svd5_df.loc[self.poc_df.index,"divergent"]
        self.poc_df["lin_as_grp_balancing_qval_SVD5"] = poc_lin_as_grp_svd5_df.loc[self.poc_df.index,"balancing"]

        # lineage as group higher hierarchy
        print("poc lin as group hh")
        poc_lin_as_grp_hh_df = pd.read_table("/home/humebc/projects/paper_4/testing_flink/selection_table_results/poc_lineage_as_group_higher_hierarchy_df.tsv")
        poc_lin_as_grp_hh_df.index = [f"Pocillopora_meandrina_contig_{chrom}:{pos}" for chrom, pos in zip(poc_lin_as_grp_hh_df["scaff_coord"],poc_lin_as_grp_hh_df["bp_coord"])]
        self.poc_df["lin_as_grp_divergent_qval_higher_hierarchy"] = poc_lin_as_grp_hh_df.loc[self.poc_df.index,"divergent"]
        self.poc_df["lin_as_grp_balancing_qval_higher_hierarchy"] = poc_lin_as_grp_hh_df.loc[self.poc_df.index,"balancing"]

        # write out poc_df
        self.poc_df.to_csv(self.new_poc_out_df_path, sep="\t")

        ## POR
        # lineage as population
        print("por lin as pop")
        por_lin_as_pop_df = pd.read_table("/home/humebc/projects/paper_4/testing_flink/selection_table_results/por_lineage_as_population_0.tsv")
        por_lin_as_pop_df.index = [f"Porites_lobata_contig_{chrom}:{pos}" for chrom, pos in zip(por_lin_as_pop_df["scaff_coord"],por_lin_as_pop_df["bp_coord"])]
        self.por_df["lin_as_pop_divergent_qval"] = [por_lin_as_pop_df.at[_,"divergent"] if _ in por_lin_as_pop_df.index else np.nan for _ in self.por_df.index]
        self.por_df["lin_as_pop_balancing_qval"] = [por_lin_as_pop_df.at[_,"balancing"] if _ in por_lin_as_pop_df.index else np.nan for _ in self.por_df.index]

        # lineage as group K1
        print("por lin as group k1")
        por_lin_as_grp_k1_df = pd.read_table("/home/humebc/projects/paper_4/testing_flink/selection_table_results/por_lineage_as_group_K1_df.tsv")
        por_lin_as_grp_k1_df.index = [f"Porites_lobata_contig_{chrom}:{pos}" for chrom, pos in zip(por_lin_as_grp_k1_df["scaff_coord"],por_lin_as_grp_k1_df["bp_coord"])]
        self.por_df["lin_as_grp_divergent_qval_K1"] = [por_lin_as_grp_k1_df.at[_,"divergent"] if _ in por_lin_as_grp_k1_df.index else np.nan for _ in self.por_df.index]
        self.por_df["lin_as_grp_balancing_qval_K1"] = [por_lin_as_grp_k1_df.at[_,"balancing"] if _ in por_lin_as_grp_k1_df.index else np.nan for _ in self.por_df.index]

        # lineage as group K2
        print("por lin as group K2")
        por_lin_as_grp_K2_df = pd.read_table("/home/humebc/projects/paper_4/testing_flink/selection_table_results/por_lineage_as_group_K2_df.tsv")
        por_lin_as_grp_K2_df.index = [f"Porites_lobata_contig_{chrom}:{pos}" for chrom, pos in zip(por_lin_as_grp_K2_df["scaff_coord"],por_lin_as_grp_K2_df["bp_coord"])]
        self.por_df["lin_as_grp_divergent_qval_K2"] = [por_lin_as_grp_K2_df.at[_,"divergent"] if _ in por_lin_as_grp_K2_df.index else np.nan for _ in self.por_df.index]
        self.por_df["lin_as_grp_balancing_qval_K2"] = [por_lin_as_grp_K2_df.at[_,"balancing"] if _ in por_lin_as_grp_K2_df.index else np.nan for _ in self.por_df.index]

        # lineage as group K3
        print("por lin as group K3")
        por_lin_as_grp_K3_df = pd.read_table("/home/humebc/projects/paper_4/testing_flink/selection_table_results/por_lineage_as_group_K3_df.tsv")
        por_lin_as_grp_K3_df.index = [f"Porites_lobata_contig_{chrom}:{pos}" for chrom, pos in zip(por_lin_as_grp_K3_df["scaff_coord"],por_lin_as_grp_K3_df["bp_coord"])]
        self.por_df["lin_as_grp_divergent_qval_K3"] = [por_lin_as_grp_K3_df.at[_,"divergent"] if _ in por_lin_as_grp_K3_df.index else np.nan for _ in self.por_df.index]
        self.por_df["lin_as_grp_balancing_qval_K3"] = [por_lin_as_grp_K3_df.at[_,"balancing"] if _ in por_lin_as_grp_K3_df.index else np.nan for _ in self.por_df.index]

        # lineage as group higher hierarchy
        print("por lin as group hh")
        por_lin_as_grp_hh_df = pd.read_table("/home/humebc/projects/paper_4/testing_flink/selection_table_results/por_lineage_as_group_higher_hierarchy_df.tsv")
        por_lin_as_grp_hh_df.index = [f"Porites_lobata_contig_{chrom}:{pos}" for chrom, pos in zip(por_lin_as_grp_hh_df["scaff_coord"],por_lin_as_grp_hh_df["bp_coord"])]
        self.por_df["lin_as_grp_divergent_qval_higher_hierarchy"] = [por_lin_as_grp_hh_df.at[_,"divergent"] if _ in por_lin_as_grp_hh_df.index else np.nan for _ in self.por_df.index]
        self.por_df["lin_as_grp_balancing_qval_higher_hierarchy"] = [por_lin_as_grp_hh_df.at[_,"balancing"] if _ in por_lin_as_grp_hh_df.index else np.nan for _ in self.por_df.index]

        # write out por_df
        self.por_df.to_csv(self.new_por_out_df_path, sep="\t")
        
        return self.poc_df, self.por_df


ScreenSelectionProfiles()._questions()