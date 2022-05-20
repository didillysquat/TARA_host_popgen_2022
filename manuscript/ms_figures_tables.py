"""
Figures and tables for the manuscript.

Making the table for the islands that were collected at and the samples
It is easy to lose track of which samples were used in which analysis.
So I will quickly run down how we get to the various numbers.
The /home/humebc/projects/paper_4/manuscript/input/Attrib_All_Por_18s.xls file is the file that Didier was
listing all of this SNP calls in and keeping various other information is as well.
I filtered down this information by filtering on the Assign_FG+MetaG getting rid of blank, MISS, K4 ?, K2/K3 (or all hybrids)
these are the samples in the POR_meta.txt. There are 322 samples. I then did a cross match of these samples with those found in the .vcf files he gave me
that contained 109 samples to get to the 108 samples used in the flink analysis. The one sample difference is I15S02C011POR which gave a hybrid designation
so I left it out of the flink analysis. The 109 number can be got to by filtering on Assign_MetaG_subcluster for only those non-blank samples.
  
"""


import site
import pandas as pd
import re
from collections import defaultdict

class FiguresTables:
    def __init__(self):
        self.sample_provenance_path = "/home/humebc/projects/paper_4/manuscript/input/TARA-PACIFIC_samples-provenance_20220131d.tsv"
        try:
            self.prov_df = pd.read_csv("/home/humebc/projects/paper_4/manuscript/input/sample_provenance_ms_figures_tables.csv")
            self.prov_df.set_index("sample-id_source", inplace=True, drop=True)
        except FileNotFoundError:
            self.prov_df = self._make_sample_provenance_df()
            self.prov_df.to_csv("/home/humebc/projects/paper_4/manuscript/input/sample_provenance_ms_figures_tables.csv", index=True)
        
        self.poc_table, self.por_table, self.mil_table = self.make_sample_tables()
        pass

    def _make_sample_provenance_df(self):
        df = pd.read_table(self.sample_provenance_path, sep="\t", skiprows=[0])
        df.set_index(keys='sample-id_source', drop=True, inplace=True)
        df = df.loc[:, ["sample_label", "sampling-design_label", "sampling-event_latitude_start_dd.dddddd", "sampling-event_longitude_start_ddd.dddddd"]]
        df.columns = ["sample_label", "sampling-design_label", "lat", "lon"]
        df = df.loc[df["sample_label"].str.contains("CS4L"),:]
        return df

    def make_sample_tables(self):
        """
        # Want to create a table for POR and POC where columns are sampling-design_label, sample_type (div or metaG), lineage and sublineage
        # Then we will want to have 3 sets of samples for the MS.
        # 1 - The metaG samples (samples that were metaG sampled) 109 for POR; There are 105 samples in the .vcf file but 106 in the Attrib_All_Poc_18s.xls.
        # the missing sample is I04S04C010POC it is in the excel but not in the vcf. Given that there is a genotype for it, it should be counted in the metaG
        # samples so we will say that there are 106 samples.
        
        # 2 - The additional divergent fragments 225 for POR (filter for value in Assign_FG, but no value in Assign_MetaG); 224 in POC (filter for 
        # DivFG in MetaG/DivFG)
        
        # 3 - Those samples that went into Flink analysis, 109 - 1 (I15S02C011POR hybrid) = 108; POR 105 in vcf - I09S03C010POC which was hybrid = 104.

        # Total samples for POR = 334; POC=330

        # Then there are the Milepora samples of which there were 57.
        """

        # POC
        poc_df = pd.read_csv("/home/humebc/projects/paper_4/manuscript/input/Attrib_All_Poc_18S.csv")
        poc_df.dropna(subset=["MetaG/DivFG"], inplace=True)
        
        # some of the rows have the sampling-design_label missing for some reason 
        # so work through the df rows and make the sampling_design_label
        sampling_design_label = []
        for ind in poc_df.index:
            if pd.isna(poc_df.at[ind, "sampling-design_label"]):
                sampling_design_label.append(self._extract_sampling_design_label(poc_df.at[ind, "INDV"]))
            else:
                sampling_design_label.append(poc_df.at[ind, "sampling-design_label"])
        poc_df["sampling-design_label"] = sampling_design_label
        poc_df.drop("INDV", axis=1, inplace=True)

        # POR
        por_df = pd.read_csv("/home/humebc/projects/paper_4/manuscript/input/Attrib_All_Por_18s.csv")
        por_df.dropna(subset=["Assign_FG+MetaG"], inplace=True)
        miss_inds = por_df.loc[por_df["Assign_FG+MetaG"] == "MISS", :].index
        por_df.drop(miss_inds, inplace=True)
        # Get whether the sample is MetaG DevFG or MetaG/DivFG
        metag_divfg = []
        for ind in por_df.index:
            if pd.isna(por_df.at[ind, "Assign_FG"]) and not pd.isna(por_df.at[ind, "Assign_MetaG"]):
                metag_divfg.append("MetaG")
            elif not pd.isna(por_df.at[ind, "Assign_FG"]) and pd.isna(por_df.at[ind, "Assign_MetaG"]):
                metag_divfg.append("DivFG")
            elif not pd.isna(por_df.at[ind, "Assign_FG"]) and not pd.isna(por_df.at[ind, "Assign_MetaG"]):
                metag_divfg.append("MetaG/DivFG")
        por_df["MetaG/DivFG"] = metag_divfg

        # TODO read in the Millepora table
        mil_df = pd.read_table("/home/humebc/projects/paper_4/manuscript/input/Millepora_PANAMA2021_attrib.txt")
        # Need to covert the IDV to a sampling-design_label
        mil_sampling_design_label = []
        for ind in mil_df.index:
            mil_sampling_design_label.append(self._extract_sampling_design_label(mil_df.at[ind, "INDV"]))
        mil_df["sampling-design_label"] = mil_sampling_design_label
        mil_df.drop("INDV", axis=1, inplace=True)
        # The milepora samples that were collected on Island 2 were a different species completely so we
        # will remove these from the df so that they are excluded from the paper
        mil_df = mil_df.loc[~mil_df["sampling-design_label"].str.contains("I02"), :]
        
        # TODO make the island site table here
        # Do a first pass to pick up all of the island sites
        self.island_site_holder = {}
        # We will collect the lat lons as a set. In theory we should end up with a single value in the set.
        self.lat_lon_sets = defaultdict(lambda: defaultdict(set))
        # island number is key to main dict, value is another dictionary with site number as key and defaultdict(int) as value
        # in the dict the keys are POC_MetaG POC_DivFG POC_MetaG_DivFG POR_MetaG POR_DivFG POC_MetaG_DivFG MIL
        for ind in poc_df.index:
            sampling_design_label = poc_df.at[ind, "sampling-design_label"]
            island, site = self._get_island_site(sampling_design_label)
            self._populate_island_site_holder(island, site)
            if poc_df.at[ind, "MetaG/DivFG"] == "MetaG":
                self.island_site_holder[island][site]["POC_MetaG"] += 1
            elif poc_df.at[ind, "MetaG/DivFG"] == "DivFG":
                self.island_site_holder[island][site]["POC_DivFG"] += 1
            elif poc_df.at[ind, "MetaG/DivFG"] == "MetaG/DivFG":
                self.island_site_holder[island][site]["POC_MetaG_DivFG"] += 1
            
            lat, lon = self._get_lat_lon(sampling_design_label)
            self.lat_lon_sets[island][site].add((lat, lon))

        for ind in por_df.index:
            sampling_design_label = por_df.at[ind, "sampling-design_label"]
            island, site = self._get_island_site(sampling_design_label)
            self._populate_island_site_holder(island, site)
            if por_df.at[ind, "MetaG/DivFG"] == "MetaG":
                self.island_site_holder[island][site]["POR_MetaG"] += 1
            elif por_df.at[ind, "MetaG/DivFG"] == "DivFG":
                self.island_site_holder[island][site]["POR_DivFG"] += 1
            elif por_df.at[ind, "MetaG/DivFG"] == "MetaG/DivFG":
                self.island_site_holder[island][site]["POR_MetaG_DivFG"] += 1
            
            lat, lon = self._get_lat_lon(sampling_design_label)
            self.lat_lon_sets[island][site].add((lat, lon))

        for ind in mil_df.index:
            sampling_design_label = mil_df.at[ind, "sampling-design_label"]
            island, site = self._get_island_site(sampling_design_label)
            self._populate_island_site_holder(island, site)
            self.island_site_holder[island][site]["MIL"] += 1
    
            lat, lon = self._get_lat_lon(sampling_design_label)
            self.lat_lon_sets[island][site].add((lat, lon))

        # Here we have the counts and the lat lons
        # We need to verify that each of the sets has only a single tuple in it
        for i_key in self.lat_lon_sets.keys():
            for s_key in self.lat_lon_sets[i_key].keys():
                if not len(self.lat_lon_sets[i_key][s_key]) == 1:
                    print(f"{i_key} {s_key}: {self.lat_lon_sets[i_key][s_key]}")
        # These are the island/sites that have multiple lat lon values
        # They are several hundered meters from each other so we will just take the first pair
        # 10 1: {(-14.010967, -171.843117), (-14.017767, -171.8425)}
        # 10 3: {(-13.918767, -171.541567), (-13.913517, -171.545067)}


        # Here we are ready to make the tables
        # Let's output the island table as a manually made tsv
        totals_dict = defaultdict(int)
        island_site_table = ["\t".join(["Island", "Site", "Latitude", "Longitude", "POC_MetaG", "POC_DivFG", "POC_MetaG_DivFG", "POR_MetaG", "POR_DivFG", "POR_MetaG_DivFG", "MIL"])]
        i_keys = sorted(self.island_site_holder.keys())
        for i_key in i_keys:
            s_keys = sorted(self.island_site_holder[i_key].keys())
            for s_key in s_keys:
                join_list_island_site_info = [str(_) for _ in [i_key, s_key, list(self.lat_lon_sets[i_key][s_key])[0][0], list(self.lat_lon_sets[i_key][s_key])[0][1]]]
                c_dict = self.island_site_holder[i_key][s_key]
                for k in ["POC_MetaG", "POC_DivFG", "POC_MetaG_DivFG", "POR_MetaG", "POR_DivFG", "POR_MetaG_DivFG", "MIL"]:
                    totals_dict[k] += c_dict[k]
                join_list_counts = [c_dict["POC_MetaG"], c_dict["POC_DivFG"], c_dict["POC_MetaG_DivFG"], c_dict["POR_MetaG"], c_dict["POR_DivFG"], c_dict["POR_MetaG_DivFG"], c_dict["MIL"]]
                join_list_island_site_info.extend([str(_) for _ in join_list_counts])
                island_site_table.append("\t".join(join_list_island_site_info))
        
        # Now write out the island table
        with open("/home/humebc/projects/paper_4/manuscript/tables/island_site_table.tsv", "w") as f:
            for line in island_site_table:
                f.write(f"{line}\n")
            totals = []
            for k in ["POC_MetaG", "POC_DivFG", "POC_MetaG_DivFG", "POR_MetaG", "POR_DivFG", "POR_MetaG_DivFG", "MIL"]:
                totals.append(str(totals_dict[k]))
            totals_str = "\t".join(totals)
            f.write(f"Totals\t\t\t\t{totals_str}\n")

        foo = "bar"

    def _populate_island_site_holder(self, island, site):
        if island not in self.island_site_holder.keys():
            self.island_site_holder[island] = {site:self._make_count_dict()}
        else:
            if not site in self.island_site_holder[island].keys():
                self.island_site_holder[island][site] = self._make_count_dict()

    def _make_count_dict(self):
        return {"POC_MetaG":0, "POC_DivFG":0, "POC_MetaG_DivFG":0, "POR_MetaG":0, "POR_DivFG":0, "POR_MetaG_DivFG":0, "MIL":0}

    def _get_lat_lon(self, sampling_design_label):
        # Grab the lat lon using the design_label
        lat_lon_df_of_sample = self.prov_df.loc[(self.prov_df["sampling-design_label"] == sampling_design_label) & (self.prov_df["sample_label"].str.contains("CS4L")), ["lat", "lon"]]
        
        lat_set = set(lat_lon_df_of_sample["lat"].values)
        assert(len(lat_set) == 1)
        lat = list(lat_set)[0]
        lon_set = set(lat_lon_df_of_sample["lon"].values)
        assert(len(lon_set) == 1)
        lon = list(lon_set)[0]
        return lat, lon

    def _get_island_site(self, ext_string):
        island_str = re.findall("I\d+", ext_string)
        assert(len(island_str) == 1)
        island_str = island_str[0]

        site_str = re.findall("S\d+", ext_string)
        assert(len(site_str) == 1)
        site_str = site_str[0]

        return int(island_str[1:]), int(site_str[1:])

    def _extract_sampling_design_label(self, ext_string):
        island_str = re.findall("I\d+", ext_string)
        assert(len(island_str) == 1)
        island_str = island_str[0]

        site_str = re.findall("S\d+", ext_string)
        assert(len(site_str) == 1)
        site_str = site_str[0]

        coral_str = re.findall("C\d+", ext_string)
        assert(len(coral_str) == 1)
        coral_str = coral_str[0]
        if len(coral_str) != 4:
            coral_str = coral_str[0] + "0" + coral_str[1:]

        return f"OA000-{island_str}-{site_str}-{coral_str}"

    def make_island_site_table(self):
        """
        Table should have columns island site lat lon POC_metaG POC_DF POC_total POR_metaG POR_DF POR_total MIL
        We can then manually add a totals row to the bottom of it
        """

FiguresTables()