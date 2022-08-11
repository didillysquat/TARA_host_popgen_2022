#!/bin/awk -f
{
    if (match($0, /^chr/)){
        # Then we are in the chromosome line.
        # We want to loop through each of the population names and extract the group
        # which will be the lineage.
        
        # First print the dashes before the groups
        printf "-\t-";
        # Then population by population extract the group and print
        for(i=3; i<=NF; i++){
            pops[i] = $i;
            sub(/I[0-9]+/, "", $i);
            printf "\t" $(i);  
        }
        print "";
        
        # Then do the population line
        printf "-\t-";
        for(j=3; j<=NF+1; j++){
            printf "\t" pops[j];
        }
        print "";
    }
    else{
        # we just print the line
        print $0;
    }
}