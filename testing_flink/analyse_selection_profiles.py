"""
General script for working with the selction profiles that we created.

"""

import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np

class PlotSelectionProfiles:
    def __init__(self):
        """ class to plot up the selction profiles."""
        self.table_in_dir = "/home/humebc/projects/paper_4/testing_flink/selection_table_results/"
        self.poc_lineage_as_population_0_df = pd.read_table(os.path.join(self.table_in_dir, "poc_lineage_as_population_0.tsv"))
        self.por_lineage_as_population_0_df = pd.read_table(os.path.join(self.table_in_dir, "por_lineage_as_population_0.tsv"))
        self.poc_lineage_as_group_SVD1_df = pd.read_table(os.path.join(self.table_in_dir, "poc_lineage_as_group_SVD1_df.tsv"))
        self.poc_lineage_as_group_SVD2_df = pd.read_table(os.path.join(self.table_in_dir, "poc_lineage_as_group_SVD2_df.tsv"))
        self.poc_lineage_as_group_SVD3_df = pd.read_table(os.path.join(self.table_in_dir, "poc_lineage_as_group_SVD3_df.tsv"))
        self.poc_lineage_as_group_SVD4_df = pd.read_table(os.path.join(self.table_in_dir, "poc_lineage_as_group_SVD4_df.tsv"))
        self.poc_lineage_as_group_SVD5_df = pd.read_table(os.path.join(self.table_in_dir, "poc_lineage_as_group_SVD5_df.tsv"))
        self.poc_lineage_as_group_higher_hierarchy_df = pd.read_table(os.path.join(self.table_in_dir, "poc_lineage_as_group_higher_hierarchy_df.tsv"))
        self.por_lineage_as_group_K1_df = pd.read_table(os.path.join(self.table_in_dir, "por_lineage_as_group_K1_df.tsv"))
        self.por_lineage_as_group_K2_df = pd.read_table(os.path.join(self.table_in_dir, "por_lineage_as_group_K2_df.tsv"))
        self.por_lineage_as_group_K3_df = pd.read_table(os.path.join(self.table_in_dir, "por_lineage_as_group_K3_df.tsv"))
        self.por_lineage_as_group_higher_hierarchy_df = pd.read_table(os.path.join(self.table_in_dir, "por_lineage_as_group_higher_hierarchy_df.tsv"))
        # list_of_dfs = [
        #     self.poc_lineage_as_population_0_df, self.por_lineage_as_population_0_df, self.poc_lineage_as_group_0_df, self.poc_lineage_as_group_1_df, self.poc_lineage_as_group_2_df, self.poc_lineage_as_group_3_df, self.poc_lineage_as_group_4_df, 
        #     self.poc_lineage_as_group_H_df, self.por_lineage_as_group_0_df, self.por_lineage_as_group_1_df, self.por_lineage_as_group_2_df, self.por_lineage_as_group_H_df
        # ]

        # We want to create a plot for POC and for POR
        # We will go with the blue gold line as it is displayed in the literature.
        # We will write this up as a method that will take a list of dfs to make the plot from.
        # self._make_plot_from_df_list(
        #     [
        #         self.poc_lineage_as_group_0_df,
        #         self.poc_lineage_as_group_1_df,
        #         self.poc_lineage_as_group_2_df,
        #         self.poc_lineage_as_group_3_df,
        #         self.poc_lineage_as_group_4_df,
        #         self.poc_lineage_as_group_H_df,
        #             ]
        #             )

        # The plotting is not really so useful as the lines merge into one.
        # self._make_plot_from_df_list([self.poc_lineage_as_group_H_df], title="foo")
        # Rather let's write out the percentage of the SNPs that are under selection
        print("proportion of SNPs under selection")
        print("comparison\tdivergent_selection\tbalancing_selection")
        self._write_out_metrics(self.poc_lineage_as_population_0_df, "poc_lineage_as_population_0_df")
        self._write_out_metrics(self.por_lineage_as_population_0_df, "por_lineage_as_population_0_df")
        self._write_out_metrics(self.poc_lineage_as_group_SVD1_df, "poc_lineage_as_group_SVD1_df")
        self._write_out_metrics(self.poc_lineage_as_group_SVD2_df, "poc_lineage_as_group_SVD2_df")
        self._write_out_metrics(self.poc_lineage_as_group_SVD3_df, "poc_lineage_as_group_SVD3_df")
        self._write_out_metrics(self.poc_lineage_as_group_SVD4_df, "poc_lineage_as_group_SVD4_df")
        self._write_out_metrics(self.poc_lineage_as_group_SVD5_df, "poc_lineage_as_group_SVD5_df")
        self._write_out_metrics(self.poc_lineage_as_group_higher_hierarchy_df, "poc_lineage_as_group_higher_hierarchy_df")
        self._write_out_metrics(self.por_lineage_as_group_K1_df, "por_lineage_as_group_K1_df")
        self._write_out_metrics(self.por_lineage_as_group_K2_df, "por_lineage_as_group_K2_df")
        self._write_out_metrics(self.por_lineage_as_group_K3_df, "por_lineage_as_group_K3_df")
        self._write_out_metrics(self.por_lineage_as_group_higher_hierarchy_df, "por_lineage_as_group_higher_hierarchy_df")

    def _write_out_metrics(self, df, title):
        sites_under_div_sel = len(df.loc[df["divergent"] < 0.01,:].index)
        sites_under_bal_sel = len(df.loc[df["balancing"] < 0.01,:].index)
        div_percent = sites_under_div_sel / len(df.index)
        bal_percent = sites_under_bal_sel / len(df.index)
        print(f"\t".join([title, f"{div_percent:.3f}", f"{bal_percent:.3f}"]))


    def _make_plot_from_df_list(self, df_list, title):
        # Number of rows will be the length of the df list
        fig, ax_arr = plt.subplots(ncols=1, nrows=len(df_list), figsize=(12,2))

        # Plot a line plot of the divergence and balancing selection
        if len(df_list) == 1:
            df=df_list[0]
            ax_arr.plot(df.index, -np.log10([_ if _!=0 else 0.000001 for _ in df["divergent"]]), linewidth=0.01, color="gold")
            ax_arr.plot(df.index, -np.log10([_ if _!=0 else 0.000001 for _ in df["balancing"]]), linewidth=0.01, color="blue")
            ax_arr.set_title(title)
            foo = "bar"
        else:
            for i, df in enumerate(df_list):
                # we want to plot the divergence against the index 
                ax_arr[i].plot(df.index, df["divergent"], linewidth=0.01, color="gold")
                plt.savefig("testing_plots.png", dpi=600)

PlotSelectionProfiles()
