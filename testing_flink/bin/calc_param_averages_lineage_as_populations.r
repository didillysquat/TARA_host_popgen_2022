# This script will calculate the parameter averages across the non-burnin
# chains.
# It will take in command line parameters
# Pos 1 = the species i.e. POC or POR
# Pos 2 = the name of the scaffold e.g. POR.Porites_lobata_contig_26
# Pos 3 = the file containing the MCMC parameters

args = commandArgs(trailingOnly = TRUE);
species = args[1]
scaffold = args[2]
chains = read.table(args[3], header=T)

averages = list()
i = 1
for (param in names(chains)[-1]){
  average = mean(chains[,param])
  # Finally plot up an average as a redline for each of the parameters
  averages[[i]] = average
  i = i + 1
}

names(averages) = names(chains)[-1]
df = data.frame(averages)

write.table(df, file=paste0(species, ".", scaffold, ".averages.tsv"), sep="\t", row.names=F, col.names=T, quote = F)
