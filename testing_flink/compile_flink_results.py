"""
UPDATE:
The read in and processing to a df was CPU limited so I wrote so multiprocessing code to 
handle the read in. For this code I did not make the complied dictionary but rather just pickled out
a dictionary for each of the scaffolds for each directory so that they could later be quickly read in and
processesed.

We have now run the flink analyses in two different structures with the lineages as either populations or as groups.
When groups, the islands are the populations.

We did this in two steps, we first ran a flink analysis estimating all of the parameters. We ran this on the contigs
that had >10000 variants (based on the flink paper saying that this was a suffcient amount to estimate parameters).
We plotted up the values that the parameters converged to and it was similar across the larger scaffolds.
This gave us confidence to use these as fixed values to then run the flink analysis on the smaller contigs.
The first estimation was done with flink.nf. The second with the fixed parameters was done with flink_fixed_parameters.nf

We did the plotting up and exploration of the first flink analysis in R locally in the script flink_testing.Rmd

It is now time to consolidate the results so that we can make a super scaffold (of all of the scaffolds combined together.)

For the lineage as population analysis there will be one superscaffold per species that shows the areas under selection and balance for the one defined group.
For the lineage as group analysis there will be 6 and 4 superscaffolds for POC and POR respectively. That is,
there will be one superscaffold per group and then one for the higher hierarchy.

The primary output from this script should be a dataframe that has the list of loci and their q value for selection and balancing.
"""

import os
import pandas as pd
import pickle
import numpy as np
import time
from multiprocessing import Pool

class CompileFlinkResults:
    def __init__(self, multiP=False):
        self.base_dir = "/home/humebc/projects/paper_4/testing_flink"
        self.out_dir = "/home/humebc/projects/paper_4/testing_flink/selection_table_results/"
        self.mp = multiP
        # Dictionaries to hold the results in 
        # Key is scaffold, value is the dictionary of coordinate to the divergent and balancing q-values
        # Population dicts
        self.poc_lineage_as_population_0_df = pd.DataFrame(columns=["divergent", "balancing"])
        self.por_lineage_as_population_0_df = pd.DataFrame(columns=["divergent", "balancing"])
        # Group dicts
        # POC
        self.poc_lineage_as_group_0_df = pd.DataFrame(columns=["divergent", "balancing"])
        self.poc_lineage_as_group_1_df = pd.DataFrame(columns=["divergent", "balancing"])
        self.poc_lineage_as_group_2_df = pd.DataFrame(columns=["divergent", "balancing"])
        self.poc_lineage_as_group_3_df = pd.DataFrame(columns=["divergent", "balancing"])
        self.poc_lineage_as_group_4_df = pd.DataFrame(columns=["divergent", "balancing"])
        self.poc_lineage_as_group_H_df = pd.DataFrame(columns=["divergent", "balancing"])
        #POR
        self.por_lineage_as_group_0_df = pd.DataFrame(columns=["divergent", "balancing"])
        self.por_lineage_as_group_1_df = pd.DataFrame(columns=["divergent", "balancing"])
        self.por_lineage_as_group_2_df = pd.DataFrame(columns=["divergent", "balancing"])
        self.por_lineage_as_group_H_df = pd.DataFrame(columns=["divergent", "balancing"])

        # start with the lineage as populations and refactor as we go
        # The results are organised per scaffold
        # We have worked out what we need to do to get the results for each scaffold in the R script flink_testing.Rmd.
        # The code here will therefore the equivalent but in python so that we can easily do this on the server (the results files are very large)
        param_list = []
        for mode in ["population", "group"]:
            for species in ["POC", "POR"]:
                for type in ["estimated", "fixed"]:
                    scaff_dir = self._get_scaff_dir(type, mode, species)
                    for scaffold in os.listdir(scaff_dir):
                        alpha_list = self._get_alpha_list(mode, species)
                        for alpha_file in alpha_list:
                            param_list.append([alpha_file, mode, species, scaff_dir, scaffold])
                            
                            
                            # self._process_alpha_file(alpha_file, mode, species, scaff_dir, scaffold)
        if self.mp: 
            # Then go through all of the scaffold directories that need processing and pickle out a dictionary
            # These dictionaries can then be collected later when running this script with multiP==False.
            with Pool(20) as p:
                p.map(self._process_alpha, param_list)
            foo = "bar"
        else:
            # Collect up all of the pickled out dictionaries
            for param_items in param_list:
                self._process_alpha(param_items)
            
            # Then its time to process each of the colleciton dictionaries into the final output
            # The output should be plain text tsv that can be read back in as a df
            # We can process these futher in a separate python script or in R.
            # We will want to sort this df by the coordinates
            foo = "bar"
            self.poc_lineage_as_population_0_df = self._sort_df_by_coord(df = self.poc_lineage_as_population_0_df)
            self.poc_lineage_as_population_0_df.to_csv(os.path.join(self.out_dir, "poc_lineage_as_population_0.tsv"), sep="\t", index=False)
            
            self.por_lineage_as_population_0_df = self._sort_df_by_coord(self.por_lineage_as_population_0_df)
            self.por_lineage_as_population_0_df.to_csv(os.path.join(self.out_dir, "por_lineage_as_population_0.tsv"), sep="\t", index=False)
            # Group dicts
            # POC
            self.poc_lineage_as_group_0_df = self._sort_df_by_coord(self.poc_lineage_as_group_0_df)
            self.poc_lineage_as_group_1_df = self._sort_df_by_coord(self.poc_lineage_as_group_1_df)
            self.poc_lineage_as_group_2_df = self._sort_df_by_coord(self.poc_lineage_as_group_2_df)
            self.poc_lineage_as_group_3_df = self._sort_df_by_coord(self.poc_lineage_as_group_3_df)
            self.poc_lineage_as_group_4_df = self._sort_df_by_coord(self.poc_lineage_as_group_4_df)
            self.poc_lineage_as_group_H_df = self._sort_df_by_coord(self.poc_lineage_as_group_H_df)
            
            self.poc_lineage_as_group_0_df.to_csv(os.path.join(self.out_dir, "poc_lineage_as_group_SVD1_df.tsv"), sep="\t", index=False)
            self.poc_lineage_as_group_1_df.to_csv(os.path.join(self.out_dir, "poc_lineage_as_group_SVD2_df.tsv"), sep="\t", index=False)
            self.poc_lineage_as_group_2_df.to_csv(os.path.join(self.out_dir, "poc_lineage_as_group_SVD3_df.tsv"), sep="\t", index=False)
            self.poc_lineage_as_group_3_df.to_csv(os.path.join(self.out_dir, "poc_lineage_as_group_SVD4_df.tsv"), sep="\t", index=False)
            self.poc_lineage_as_group_4_df.to_csv(os.path.join(self.out_dir, "poc_lineage_as_group_SVD5_df.tsv"), sep="\t", index=False)
            self.poc_lineage_as_group_H_df.to_csv(os.path.join(self.out_dir, "poc_lineage_as_group_higher_hierarchy_df.tsv"), sep="\t", index=False)
            #POR
            self.por_lineage_as_group_0_df = self._sort_df_by_coord(self.por_lineage_as_group_0_df)
            self.por_lineage_as_group_0_df.to_csv(os.path.join(self.out_dir, "por_lineage_as_group_K1_df.tsv"), sep="\t", index=False)
            self.por_lineage_as_group_1_df = self._sort_df_by_coord(self.por_lineage_as_group_1_df)
            self.por_lineage_as_group_1_df.to_csv(os.path.join(self.out_dir, "por_lineage_as_group_K2_df.tsv"), sep="\t", index=False)
            self.por_lineage_as_group_2_df = self._sort_df_by_coord(self.por_lineage_as_group_2_df)
            self.por_lineage_as_group_2_df.to_csv(os.path.join(self.out_dir, "por_lineage_as_group_K3_df.tsv"), sep="\t", index=False)
            self.por_lineage_as_group_H_df = self._sort_df_by_coord(self.por_lineage_as_group_H_df)
            self.por_lineage_as_group_H_df.to_csv(os.path.join(self.out_dir, "por_lineage_as_group_higher_hierarchy_df.tsv"), sep="\t", index=False)
            
    def _sort_df_by_coord(self, df):
        scaff_coord = [int(_.split(":")[0].split("_")[-1]) for _ in df.index]
        bp_coord = [int(_.split(":")[1]) for _ in df.index]
        df["scaff_coord"] = scaff_coord
        df["bp_coord"] = bp_coord
        # now sort first by scaff_coord then by bp_coord
        return df.sort_values(["scaff_coord", "bp_coord"], axis=0)
    
    def _process_alpha(self, param_items):
        alpha_file = param_items[0]
        mode = param_items[1]
        species = param_items[2]
        scaff_dir = param_items[3]
        scaffold = param_items[4]
        self._process_alpha_file(alpha_file, mode, species, scaff_dir, scaffold)
        
    def _get_alpha_list(self, mode, species):
        if mode == "population":
            if species == "POC":
                return ["Posterior_alphas_group_0.txt"]
            elif species == "POR":
                return ["Posterior_alphas_group_0.txt"]
            
        elif mode == "group":
            # Then we have to do an alpha for each group and for the hierarchical
            if species == "POC":
                return ["Posterior_alphas_group_0.txt", "Posterior_alphas_group_1.txt","Posterior_alphas_group_2.txt","Posterior_alphas_group_3.txt","Posterior_alphas_group_4.txt", "Posterior_A.txt"]
            elif species == "POR":
                return ["Posterior_alphas_group_0.txt", "Posterior_alphas_group_1.txt", "Posterior_alphas_group_2.txt", "Posterior_A.txt"]

    def _process_alpha_file(self, alpha_file, mode, species, scaff_dir, scaffold):
        # Look for a pickled out dictionary of the divergent/balancing selection
        df_pickle_path = os.path.join(scaff_dir, scaffold, alpha_file.replace(".txt", ".dict.p"))
        if os.path.exists(df_pickle_path):
            # selection_df = self._load_selection_df(df_pickle_path)
            if self.mp:
                # Then we are not concerned with collecting up the dictionaries
                pass
            else:
                # THen we are concerned with collecting up the dictionaries
                selection_df = self._load_selection_df(df_pickle_path=df_pickle_path)
                self._deposit_selection_df(selection_df=selection_df, species=species, mode=mode, scaffold=scaffold, alpha_file=alpha_file)

        else:
            if "fixed" in scaff_dir:
                if mode == "population":
                    flink_in_path = os.path.join(scaff_dir, scaffold, f"{species}.{scaffold}.flink.in.as.populations.txt")
                elif mode == "group":
                    flink_in_path = os.path.join(scaff_dir, scaffold, f"{species}.{scaffold}.flink.in.as.groups.txt")
            else:
                flink_in_path=os.path.join(scaff_dir, scaffold, f"{species}.{scaffold}.flink.in.txt")

            selection_df = self._make_selection_df(
                alphas_path=os.path.join(scaff_dir, scaffold, alpha_file),
                flink_in_path=flink_in_path,
                pickle_out_path=df_pickle_path
                )
        # self._deposit_selection_df(selection_df=selection_df, species=species, mode=mode, scaffold=scaffold)

    def _load_selection_df(self, df_pickle_path):
        # Then load up with dict and move on
        print(f"loading dict: {df_pickle_path}")
        with open(df_pickle_path, "rb") as f:
            return pickle.load(f)
                
    def _get_scaff_dir(self, type, mode, species):
        if type == "estimated":
            type_sub_path = "flink_results"
        else:
            type_sub_path = "flink_fixed_parameters/flink_results"
        if mode == "population":
            mode_sub_path = "lineages_as_populations"
        else:
            mode_sub_path = "lineages_as_groups"
        return os.path.join(self.base_dir, type_sub_path, mode_sub_path, species)
    
    def _deposit_selection_df(self, selection_df, species, mode, scaffold, alpha_file="Posterior_alphas_group_0.txt"):
        print(f"Depositing {mode} {species} {scaffold}")
        if mode == "population":
            if species == "POC":
                self.poc_lineage_as_population_0_df = pd.concat([self.poc_lineage_as_population_0_df, selection_df])
            else:
                self.por_lineage_as_population_0_df = pd.concat([self.por_lineage_as_population_0_df, selection_df])
        else: # group
            if species == "POC":
                if alpha_file == "Posterior_alphas_group_0.txt":
                    self.poc_lineage_as_group_0_df = pd.concat([self.poc_lineage_as_group_0_df, selection_df]) 
                elif alpha_file == "Posterior_alphas_group_1.txt":
                    self.poc_lineage_as_group_1_df = pd.concat([self.poc_lineage_as_group_1_df, selection_df])
                elif alpha_file == "Posterior_alphas_group_2.txt":
                    self.poc_lineage_as_group_2_df = pd.concat([self.poc_lineage_as_group_2_df, selection_df])
                elif alpha_file == "Posterior_alphas_group_3.txt":
                    self.poc_lineage_as_group_3_df = pd.concat([self.poc_lineage_as_group_3_df, selection_df])
                elif alpha_file == "Posterior_alphas_group_4.txt":
                    self.poc_lineage_as_group_4_df = pd.concat([self.poc_lineage_as_group_4_df, selection_df])
                elif alpha_file == "Posterior_A.txt":
                    self.poc_lineage_as_group_H_df = pd.concat([self.poc_lineage_as_group_H_df, selection_df])
            else: # POR
                if alpha_file == "Posterior_alphas_group_0.txt":
                    self.por_lineage_as_group_0_df = pd.concat([self.por_lineage_as_group_0_df, selection_df])
                elif alpha_file == "Posterior_alphas_group_1.txt":
                    self.por_lineage_as_group_1_df = pd.concat([self.por_lineage_as_group_1_df, selection_df])
                elif alpha_file == "Posterior_alphas_group_2.txt":
                    self.por_lineage_as_group_2_df = pd.concat([self.por_lineage_as_group_2_df, selection_df])
                elif alpha_file == "Posterior_A.txt":
                    self.por_lineage_as_group_H_df = pd.concat([self.por_lineage_as_group_H_df, selection_df])

    def _make_selection_df(self, alphas_path, flink_in_path, pickle_out_path):
        print(f"Reading in {alphas_path}")
        
        # The file that we are interested in is the Posterior_alphas_group_0.txt   
        start = time.process_time()
        alphas_df = pd.read_table(alphas_path, index_col=0, dtype=np.float64)
        print(f"Read complete in: {time.process_time() - start}s")
        
        div_selection = 1-((alphas_df > 0).sum(axis=0) / alphas_df.shape[0])
        bal_selection = 1-((alphas_df < 0).sum(axis=0) / alphas_df.shape[0])
        
        # The positions are given arbitrary names alpha_xxx from 0..num_pos
        # We want to put these back into coordinates
        # To do this we read in the flink input file
        with open(flink_in_path, "r") as f:
            coords = [":".join(_.rstrip().split("\t")[:2]) for _ in f][2:]
        
        selection_df = pd.DataFrame({"divergent": div_selection, "balancing": bal_selection})
        selection_df.index = coords

        with open(pickle_out_path, "wb") as f:
            pickle.dump(selection_df, f)

        return selection_df

if __name__ == '__main__':
    CompileFlinkResults(multiP=False)
                    