"""
This script takes three command line inputs
1 - species: POR or POC
2 - groups_or_populations: group or population
3 - the *alleleCounts.txt.gz file output from atlas.

The aim of this script is to rearange the columns of the *alleleCounts.txt.gz
files so that they are in proper order of the groups and populations.

And to add the proper headers to the output file.

The modified output file will then be used as input to the flink analysis
"""

import sys
import gzip
import pandas as pd
import re

class ModifyAtlasOutput:
    def __init__(self):
        self.species = sys.argv[1]
        self.scaffold = sys.argv[2]
        self.mode = sys.argv[3]
        # with gzip.open(sys.argv[3], "rb") as f:
        #     self.atlas_output_lines = [_.rstrip() for _ in f]
        
        self.output = []

        if self.mode == "population":
            # Need to reorder the columns
            self.df = pd.read_table(sys.argv[4] ,compression='gzip')
            columns_names = list(self.df)[2:]
            columns_names.sort()
            new_col_order = ["chr", "pos"] + columns_names
            self.df = self.df.loc[:,new_col_order]
            if self.species == "POC":
                self.df.columns = ["chr", "pos"] + [f"SVD{_}" for _ in columns_names]
            elif self.species == "POR":
                self.df.columns = ["chr", "pos"] + [f"K{_}" for _ in columns_names]
            else:
                raise RuntimeError("unrecognised species")
            
            # Now write out with the additional columns
            with open(f"{self.species}.{self.scaffold}.flink.in.as.populations.txt", "w") as f:
                group_line = ["-", "-"] + ["Group1" for _ in list(self.df)[2:]]
                f.write("\t".join(group_line) + "\n")
                population_line = ["-", "-"] + [_ for _ in list(self.df)[2:]]
                f.write("\t".join(population_line) + "\n")
                # Then write out the data lines
                for ind in self.df.index:
                    out_line = [str(_) for _ in self.df.loc[ind,:].values]
                    f.write("\t".join(out_line) + "\n")
        elif self.mode == "group":
            # Need to reorder the columns
            self.df = pd.read_table(sys.argv[4] ,compression='gzip')
            columns_names = list(self.df)[2:]
            columns_names.sort()
            new_col_order = ["chr", "pos"] + columns_names
            self.df = self.df.loc[:,new_col_order]
            
            # Now write out with the additional columns
            with open(f"{self.species}.{self.scaffold}.flink.in.as.groups.txt", "w") as f:
                # Need to pull out the groups from the columns_names
                if self.species == "POC":
                    group_line = ["-", "-"] + [re.findall("[SVD]\d{1}", _)[0] for _ in columns_names]
                elif self.species == "POR":
                    group_line = ["-", "-"] + [re.findall("K\d{1}", _)[0] for _ in columns_names]
                else:
                    raise RuntimeError("unrecognised species")
                f.write("\t".join(group_line) + "\n")
                population_line = ["-", "-"] + [_ for _ in list(self.df)[2:]]
                f.write("\t".join(population_line) + "\n")
                # Then write out the data lines
                for ind in self.df.index:
                    out_line = [str(_) for _ in self.df.loc[ind,:].values]
                    f.write("\t".join(out_line) + "\n")
        else:
            raise RuntimeError("unknown mode")
            
ModifyAtlasOutput()