# 049_drakes_code


put folder into env variable

    nano .bahsrc
    export gbm_049_drake=/home/daniel/ubuntu/workspace/drakes_data/18-0190-049/


## start with Drake's output from Stringtie (FPKM counts)

    
    cd gbm_049_drake/expression/stringtie/ref_only
    less -S Sample6_Lane2/transcripts.gtf

   
## Use Ballgown in R for differential expression (DE) analysis (then PCA) using output from Stringtie
Perform A vs. B comparison, using all replicates, for known (reference only mode) transcripts

    mkdir -p ~workspace/drakes_data/18-0190-049/de/ballgown/ref_only
    cd $gbm_049_drake/de/ballgown/ref_only/
  

Use printf to create/print a table with ids, type (each type of sample is a type), and path to the file, as the header. Then n returns a new line.
##(note: id needs to match the file folder names created by stringtie)

Bascially, need a table that needs to look like this to feed into R:

ids type path-to-file-011_invitro_1 011 $gbm/expression/stringtie/1 011_invitro_2 011 $gbm/expression/stringtie/2 ... ...

goal is to generate a header file to load into R, for ALL samples for principal component analysis (the simplest form of multidimentional scaling), and also a file for pairwise comparisons. since we have a ton of comparisisons, might just not do this for now and only do the PCA.

file for all 049 samples for PCA: (this is how the script should look like (without the enters inbetween each line): 

printf ""ids","type","path 
"\n"Sample6_Lane2","011_slice","$gbm_049_drake/expression/stringtie/ref_only/Sample6_Lane2 "\n"2","011_slice","$gbm/expression/stringtie/ref_only/2 "\n"3","011_organoid","$gbm/expression/stringtie/ref_only/3 "\n"4","011_organoid","$gbm/expression/stringtie/ref_only/4 "\n"5","011_tissue","$gbm/expression/stringtie/ref_only/5 "\n"6","011_tissue","$gbm/expression/stringtie/ref_only/6 "\n"7","011_invitro","$gbm/expression/stringtie/ref_only/7 "\n"8","011_invitro","$gbm/expression/stringtie/ref_only/8 "\n" > GBM011_all.csv


    cd $gbm/de/ballgown/ref_only/

    printf "\"ids\",\"type\",\"path\"\n\"1\",\"011_slice\",\"$gbm/expression/stringtie/ref_only/1\"\n\"2\",\"011_slice\",\"$gbm/expression/stringtie/ref_only/2\"\n\"3\",\"011_organoid\",\"$gbm/expression/stringtie/ref_only/3\"\n\"4\",\"011_organoid\",\"$gbm/expression/stringtie/ref_only/4\"\n\"5\",\"011_tissue\",\"$gbm/expression/stringtie/ref_only/5\"\n\"6\",\"011_tissue\",\"$gbm/expression/stringtie/ref_only/6\"\n\"7\",\"011_invitro\",\"$gbm/expression/stringtie/ref_only/7\"\n\"8\",\"011_invitro\",\"$gbm/expression/stringtie/ref_only/8\"\n" > GBM011_all.csv



