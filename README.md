# Tara Pacific coral host delination FLINK analysis

This repository is related to the FLINK analysis run as part of the manuscript.
The main set of FLINK analyses were run in the nextflow script: testing_flink/flink.nf
In this set, we worked with the scaffolds that were > 10 000 bp long.
We derived the population statistics from this set of runs, and then ran a
second set of flink analyses using testing_flink/flink_fixed_parameters/flink_fixed_parameters.nf where the parameters were fixed.
Details of the files used as input are detailed in the workflows themselves.

