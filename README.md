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
  
  
