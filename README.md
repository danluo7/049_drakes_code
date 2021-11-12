
# 049_drakes_code
below scripts are for data in /home/daniel/ubuntu/workspace/all_049/drakes_data/18-0190-049



#put folder into env variable

    nano .bahsrc
    export gbm_049_drake=/home/daniel/ubuntu/workspace/all_049/drakes_data/18-0190-049/


## start with Drake's output from Stringtie (FPKM counts)

    cd $gbm_049_drake/expression/stringtie/ref_only
    head Sample6_Lane2/transcripts.gtf

   
## Use Ballgown in R for differential expression (DE) analysis (then PCA) using output from Stringtie
Perform A vs. B comparison, using all replicates, for known (reference only mode) transcripts

    mkdir -p ~/workspace/all_049/drakes_data/18-0190-049/de/ballgown/ref_only
    cd $gbm_049_drake/de/ballgown/ref_only/
  

#Use printf to create/print a table with ids, type (each type of sample is a type), and path to the file, as the header. Then n returns a new line.
##(note: id needs to match the file folder names created by stringtie)

#Bascially, need a table that needs to look like this to feed into R:

#ids type path-to-file-011_invitro_1 011 $gbm/expression/stringtie/1 011_invitro_2 011 $gbm/expression/stringtie/2 ... ...

#goal is to generate a header file to load into R, for ALL samples for principal component analysis (the simplest form of multidimentional scaling), and also a file for pairwise comparisons. since we have a ton of comparisisons, might just not do this for now and only do the PCA.

#file for all 049 samples for PCA: 


#printf "\"ids\",\"type\",\"path
\"\n\"Sample6_Lane2\",\"049_tissue\",\"/home/daniel/ubuntu/workspace/all_049/drakes_data/18-0190-049/expression/stringtie/ref_only/Sample6_Lane2
\"\n\"Sample6_Lane3\",\"049_tissue\",\"/home/daniel/ubuntu/workspace/all_049/drakes_data/18-0190-049/expression/stringtie/ref_only/Sample6_Lane3
\"\n\"Sample17_Lane2\",\"049_slice\",\"/home/daniel/ubuntu/workspace/all_049/drakes_data/18-0190-049/expression/stringtie/ref_only/Sample17_Lane2
\"\n\"Sample17_Lane3\",\"049_slice\",\"/home/daniel/ubuntu/workspace/all_049/drakes_data/18-0190-049/expression/stringtie/ref_only/Sample17_Lane3
\"\n\"Sample18_Lane2\",\"049_neurosphere\",\"/home/daniel/ubuntu/workspace/all_049/drakes_data/18-0190-049/expression/stringtie/ref_only/Sample18_Lane2
\"\n\"Sample18_Lane3\",\"049_neurosphere\",\"/home/daniel/ubuntu/workspace/all_049/drakes_data/18-0190-049/expression/stringtie/ref_only/Sample18_Lane3
\"\n\"Sample29_Lane2\",\"049_organoid\",\"/home/daniel/ubuntu/workspace/all_049/drakes_data/18-0190-049/expression/stringtie/ref_only/Sample29_Lane3
\"\n\"Sample29_Lane3\",\"049_organoid\",\"/home/daniel/ubuntu/workspace/all_049/drakes_data/18-0190-049/expression/stringtie/ref_only/Sample29_Lane3
\"\n\"Sample26_Lane2\",\"049_invitro\",\"/home/daniel/ubuntu/workspace/all_049/drakes_data/18-0190-049/expression/stringtie/ref_only/Sample26_Lane2
\"\n\"Sample27_Lane2\",\"049_invitro\",\"/home/daniel/ubuntu/workspace/all_049/drakes_data/18-0190-049/expression/stringtie/ref_only/Sample27_Lane2
\"\n" > GBM049_all_drake.csv


	printf "\"ids\",\"type\",\"path\"\n\"
Sample6_Lane2\",\"049_tissue\",\"/home/daniel/ubuntu/workspace/all_049/drakes_data/18-0190-049/expression/stringtie/ref_only/Sample6_Lane2\"\n\"Sample6_Lane3\",\"049_tissue\",\"/home/daniel/ubuntu/workspace/all_049/drakes_data/18-0190-049/expression/stringtie/ref_only/Sample6_Lane3\"\n\"Sample17_Lane2\",\"049_slice\",\"/home/daniel/ubuntu/workspace/all_049/drakes_data/18-0190-049/expression/stringtie/ref_only/Sample17_Lane2\"\n\"Sample17_Lane3\",\"049_slice\",\"/home/daniel/ubuntu/workspace/all_049/drakes_data/18-0190-049/expression/stringtie/ref_only/Sample17_Lane3\"\n\"Sample18_Lane2\",\"049_neurosphere\",\"/home/daniel/ubuntu/workspace/all_049/drakes_data/18-0190-049/expression/stringtie/ref_only/Sample18_Lane2\"\n\"Sample18_Lane3\",\"049_neurosphere\",\"/home/daniel/ubuntu/workspace/all_049/drakes_data/18-0190-049/expression/stringtie/ref_only/Sample18_Lane3\"\n\"Sample29_Lane2\",\"049_organoid\",\"/home/daniel/ubuntu/workspace/all_049/drakes_data/18-0190-049/expression/stringtie/ref_only/Sample29_Lane3\"\n\"Sample29_Lane3\",\"049_organoid\",\"/home/daniel/ubuntu/workspace/all_049/drakes_data/18-0190-049/expression/stringtie/ref_only/Sample29_Lane3\"\n\"Sample26_Lane2\",\"049_invitro\",\"/home/daniel/ubuntu/workspace/all_049/drakes_data/18-0190-049/expression/stringtie/ref_only/Sample26_Lane2\"\n\"Sample27_Lane2\",\"049_invitro\",\"/home/daniel/ubuntu/workspace/all_049/drakes_data/18-0190-049/expression/stringtie/ref_only/Sample27_Lane2\"\n" > GBM049_all_drake.csv
	
	cat GBM049_all_drake.csv



#R script:

	R --no-restore
	library(ballgown)
	library(genefilter)
	library(dplyr)
	library(devtools)
	library(ggplot2)
	library(gplots)
	library(GenomicRanges)


	working_dir = "/home/daniel/ubuntu/workspace/all_049/drakes_data/18-0190-049/de/ballgown/ref_only"
	setwd(working_dir)
	dir()


	pheno_data = read.csv("GBM049_all_drake.csv")  


	bg = ballgown(samples=as.vector(pheno_data$path), pData=pheno_data)
	bg

	bg_table = texpr(bg, 'all')

	bg_gene_names = unique(bg_table[, 9:10])
	head(bg_gene_names)

	save(bg, file='bg.rda')
	bg


	pdf(file="GBM049_R_output.pdf")

	


#Import expression and differential expression results from the HISAT2/StringTie/Ballgown pipeline
	
	load('bg.rda')


#Load gene names for lookup later in the tutorial
	
	bg_table = texpr(bg, 'all')
	bg_gene_names = unique(bg_table[, 9:10])


#Pull the gene_expression data frame from the ballgown object

	gene_expression = as.data.frame(gexpr(bg))
	head(gene_expression)

#Get the first 3 rows of data and a selection of columns
	
	gene_expression[1:3,c(1:3,4)]


#View the column names

	colnames(gene_expression)

#View the row names by row.names(gene_expression) (will be long)

#Determine the dimensions of the dataframe.  'dim()' will return the number of rows and columns

	dim(gene_expression)



#Just for fun, check BRD4 expression across all 8 samples:

	i = row.names(gene_expression) == "NM_058243"
	gene_expression[i,]



#Load the transcript to gene index from the ballgown object. Each row of data represents a transcript. Many of these transcripts represent the same gene. Determine the numbers of transcripts and unique genes  


	transcript_gene_table = indexes(bg)$t2g
	head(transcript_gene_table)
	
	length(row.names(transcript_gene_table)) 
	length(unique(transcript_gene_table[,"g_id"])) 
	
#output of Transcript count

#< [1] 59197

#output of Unique Gene count

#< [1] 51276



#Plot #1: Plot the number of transcripts per gene. Many genes will have only 1 transcript, some genes will have several transcripts. Use the 'table()' command to count the number of times each gene symbol occurs (i.e. the # of transcripts that have each gene symbol). Then use the 'hist' command to create a histogram of these counts

	counts=table(transcript_gene_table[,"g_id"])
	c_one = length(which(counts == 1))
	c_more_than_one = length(which(counts > 1))
	c_max = max(counts)
	hist(counts, breaks=50, col="bisque4", xlab="Transcripts per gene", main="Distribution of transcript count per gene")
	legend_text = c(paste("Genes with one transcript =", c_one), paste("Genes with more than one transcript =", c_more_than_one), paste("Max transcripts for single gene = ", c_max))
	legend("topright", legend_text, lty=NULL)


#Plot #2: Plot the distribution of transcript sizes as a histogram. lengths will be those of known transcripts. Good QC step: we had a low coverage library, or other problems, we might get short 'transcripts' that are actually only pieces of real transcripts.

	full_table <- texpr(bg , 'all')
	hist(full_table$length, breaks=500, xlab="Transcript length (bp)", main="Distribution of transcript lengths", col="steelblue")

#View the summary FPKM values (minimum and maximum FPKM values) for any particular library, in this case, Sample17_Lane2
	
	min(gene_expression[,"FPKM.Sample17_Lane2"])
	max(gene_expression[,"FPKM.Sample17_Lane2"])


#Set the minimum non-zero FPKM values by one of two ways:

#coverting 0's to NA, and calculating the minimum or all non NA values
#two ways: 
#zz = fpkm_matrix[,data_columns]
#zz[zz==0] = NA
#min_nonzero = min(zz, na.rm=TRUE)
#min_nonzero


#Alternatively just set min value to 1
	
	min_nonzero=1

#Set the columns for finding FPKM and create shorter names for figures


	data_columns=c(1:10)
	short_names=c("tissue_1","tissue_2","organoid_1","organoid_2","discells_1", "discells_2","slice_1","slice_2","invitro_1","invitro_2")




#Plot #3: Plot range of values and general distribution of FPKM values for all 10 libraries

#Create boxplots using different colors by setting storing the colors of the columns in a variable called data_colors. then display on a log2 scale and add the minimum non-zero value to avoid log2(0). Note that the bold horizontal line on each boxplot is the median.


	#colors()
	data_colors=c("tomato1","tomato2","royalblue1","royalblue2","seagreen1","seagreen2","grey1","grey2","brown1","brown2")

	boxplot(log2(gene_expression[,data_columns]+min_nonzero), col=data_colors, names=short_names, las=2, ylab="log2(FPKM)", main="Distribution of FPKMs for libraries of all 10 samples")






-----------------------------------------
## plot #4 (optional): plot a pair of replicates to assess reproducibility of technical replicates. 
Tranform the data by converting to log2 scale after adding an arbitrary small value to avoid log2(0). Also add a straight line of slope 1, and intercept 0. Also calculate the correlation coefficient and display in a legend.

	x = gene_expression[,"FPKM.1"]
	y = gene_expression[,"FPKM.2"]
	plot(x=log2(x+min_nonzero), y=log2(y+min_nonzero), pch=16, col="blue", cex=0.25, xlab="FPKM (011_slice, Replicate 1)", ylab="FPKM (011_slice, Replicate 2)", main="Comparison of expression values for replicates of slice samples")
	abline(a=0,b=1)
	rs=cor(x,y)^2
	legend("topleft", paste("R squared = ", round(rs, digits=3), sep=""), lwd=1, col="black")


check organoid samples

	x = gene_expression[,"FPKM.3"]
	y = gene_expression[,"FPKM.4"]
	plot(x=log2(x+min_nonzero), y=log2(y+min_nonzero), pch=16, col="blue", cex=0.25, xlab="FPKM (011_organoid, Replicate 1)", ylab="FPKM (011_organoid, Replicate 2)", main="Comparison of expression values for replicates of organoid samples")
	abline(a=0,b=1)
	rs=cor(x,y)^2
	legend("topleft", paste("R squared = ", round(rs, digits=3), sep=""), lwd=1, col="black")

check tissue samples

	x = gene_expression[,"FPKM.5"]
	y = gene_expression[,"FPKM.6"]
	plot(x=log2(x+min_nonzero), y=log2(y+min_nonzero), pch=16, col="blue", cex=0.25, xlab="FPKM (011_tissue, Replicate 1)", ylab="FPKM (011_tissue, Replicate 2)", main="Comparison of expression values for replicates of tissue samples")
	abline(a=0,b=1)
	rs=cor(x,y)^2
	legend("topleft", paste("R squared = ", round(rs, digits=3), sep=""), lwd=1, col="black")


check in vitro samples

	x = gene_expression[,"FPKM.7"]
	y = gene_expression[,"FPKM.8"]
	plot(x=log2(x+min_nonzero), y=log2(y+min_nonzero), pch=16, col="blue", cex=0.25, xlab="FPKM (011_invitro, Replicate 1)", ylab="FPKM (011_invitro, Replicate 2)", main="Comparison of expression values for replicates of in vitro samples")
	abline(a=0,b=1)
	rs=cor(x,y)^2
	legend("topleft", paste("R squared = ", round(rs, digits=3), sep=""), lwd=1, col="black")
------------------------------------------------









## Compare the correlation distance between all replicates


#Scree plot: Determine the amount of variance coming from each principal component in a table:

	pc <- princomp(gene_expression[,data_columns],cor=TRUE,scores=TRUE)
	summary(pc)
	plot(pc,type='lines')

#Calculate the FPKM sum for all 10 libraries

	gene_expression[,"sum"]=apply(gene_expression[,data_columns], 1, sum)

#Filter out genes with a grand sum FPKM of less than 10

	i = which(gene_expression[,"sum"] > 10)


#Calculate the correlation between all pairs of data

	r=cor(gene_expression[i,data_columns], use="pairwise.complete.obs", method="pearson")
	r
	
	
## Plot MDS (plot #8)
Convert correlation to distance, and use 'multi-dimensional scaling' to plot the relative differences between libraries, by calculating 2-dimensional coordinates to plot points for each library using eigenvectors (eig=TRUE). d, k=2 means 2 dimensions
	
	d=1-r
	mds=cmdscale(d, k=2, eig=TRUE)
	par(mfrow=c(1,1))
	plot(mds$points, type="n", xlab="", ylab="", main="MDS distance plot (all non-zero genes)", xlim=c(-0.12,0.12), ylim=c(-0.12,0.12))
	points(mds$points[,1], mds$points[,2], col="grey", cex=2, pch=16)
	text(mds$points[,1], mds$points[,2], short_names, col=data_colors)


	d=1-r
	mds=cmdscale(d, k=2, eig=TRUE)
	par(mfrow=c(1,1))
	plot(mds$points, type="n", xlab="", ylab="", main="MDS distance plot (all non-zero genes)", xlim=c(-0.2,0.2), ylim=c(-0.2,0.2))
	points(mds$points[,1], mds$points[,2], col="grey", cex=2, pch=16)
	text(mds$points[,1], mds$points[,2], short_names, col=data_colors)



	dev.off()

---------------------------
#below code doesn't work anymore becuase of stattest
#Calculate the differential expression results including significance

	results_genes = stattest(bg, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
	results_genes = merge(results_genes,bg_gene_names,by.x=c("id"),by.y=c("gene_id"))



## Plot #9 - View the distribution of differential expression values as a histogram

Display only those that are significant according to Ballgown

	sig=which(results_genes$pval<0.05)
	results_genes[,"de"] = log2(results_genes[,"fc"])
	hist(results_genes[sig,"de"], breaks=50, col="seagreen", xlab="log2(Fold change) UHR vs HBR", main="Distribution of differential expression values")
	abline(v=-2, col="black", lwd=2, lty=2)
	abline(v=2, col="black", lwd=2, lty=2)
	legend("topleft", "Fold-change > 4", lwd=2, lty=2)

## Plot #10 - Display the grand expression values from UHR and HBR and mark those that are significantly differentially expressed

	gene_expression[,"brain slice"]=apply(gene_expression[,c(1:3)], 1, mean)
	gene_expression[,"in vitro"]=apply(gene_expression[,c(4:6)], 1, mean)

	x=log2(gene_expression[,"brain slice"]+min_nonzero)
	y=log2(gene_expression[,"in vitro"]+min_nonzero)
	plot(x=x, y=y, pch=16, cex=0.25, xlab="brain slice FPKM (log2)", ylab="in vitro FPKM (log2)", main="Brain Slice vs in vitro FPKMs")
	abline(a=0, b=1)
	xsig=x[sig]
	ysig=y[sig]
	points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
	legend("topleft", "Significant", col="magenta", pch=16)

Get the gene symbols for the top N (according to corrected p-value) and display them on the plot

	topn = order(abs(results_genes[sig,"fc"]), decreasing=TRUE)[1:25]
	topn = order(results_genes[sig,"qval"])[1:25]
	text(x[topn], y[topn], results_genes[topn,"gene_name"], col="black", cex=0.75, srt=45)


## Write a simple table of differentially expressed transcripts to an output file

Each should be significant with a log2 fold-change >= 2

	sigpi = which(results_genes[,"pval"]<0.05)
	sigp = results_genes[sigpi,]
	sigde = which(abs(sigp[,"de"]) >= 2)
	sig_tn_de = sigp[sigde,]

Order the output by or p-value and then break ties using fold-change

o = order(sig_tn_de[,"qval"], -abs(sig_tn_de[,"de"]), decreasing=FALSE)

	output = sig_tn_de[o,c("gene_name","id","fc","pval","qval","de")]
	write.table(output, file="SigDE_supplementary_R.txt", sep="\t", row.names=FALSE, quote=FALSE)

View selected columns of the first 25 lines of output

	output[1:25,c(1,4,5)]


You can open the file "SigDE.txt" in Excel, Calc, etc.
It should have been written to the current working directory that you set at the beginning of the R tutorial
dir()


## Plot #11 - Create a heatmap to vizualize expression differences between the eight samples

	Define custom dist and hclust functions for use with heatmaps
	mydist=function(c) {dist(c,method="euclidian")}
	myclust=function(c) {hclust(c,method="complete")}

	main_title="sig DE Transcripts"
	par(cex.main=0.8)
	sig_genes=results_genes[sig,"id"]
	sig_gene_names=results_genes[sig,"gene_name"]
	data=log2(as.matrix(gene_expression[sig_genes,data_columns])+1)
	heatmap.2(data, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="column", dendrogram="both", margins=c(6,7), Rowv=TRUE, Colv=TRUE, symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", main=main_title, cexRow=0.3, cexCol=1, labRow=sig_gene_names,labCol=short_names,col=bluered(100))

	dev.off()





  
