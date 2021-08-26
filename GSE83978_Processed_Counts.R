### Loading the series matrix files to see whats there. this code is from Geo2R. 
## It puts the data into an ExpressionSet object which can be easily subsetted to look inside
library(GEOquery)

gset <- getGEO("GSE83978", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL17021", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

#Extract phenotype data 
Phenodata <- pData(phenoData(gset))


#Create a new DF with important sample Info

colnames(Phenodata) ## Select columns from indices. I usually select a lot then cut it down after looking at it
SampleINFO <- Phenodata[, c(1,8)]

colnames(SampleINFO) <- c("title", "Status_1")

## Wasnt sure what is interesting to compare in this dataset so I made columns for all variables

SampleINFO$Donor <- ifelse(grepl("Donnor1", SampleINFO$title), "1", 
                           ifelse(grepl("Donnor2", SampleINFO$title) , "2", 
                                  ifelse(grepl("Donnor3", SampleINFO$title) , "3",
                                         "4")))



SampleINFO$Infection <- ifelse(grepl("acute", SampleINFO$Status_1), "Acute", 
                               ifelse(grepl("chronic", SampleINFO$Status_1) , "Chronic", 
                                       "Naive"))

SampleINFO$Tcf1 <- ifelse(grepl("Tcf1-", SampleINFO$title), "Negative", 
                               ifelse(grepl("Tcf1+", SampleINFO$title) , "Positive", 
                                      "NA"))

SampleINFO$GFP <- ifelse(grepl("GFP-", SampleINFO$Status_1), "Neg", 
                                      "Pos")

SampleINFO$Phenotype <- ifelse(grepl("Memory Tcf1-",SampleINFO$title),"Memory.Neg",
                         ifelse(grepl("Memory Tcf1+",SampleINFO$title),"Memory.Pos",
                                ifelse(grepl("Chronic Tcf1-",SampleINFO$title),"Chronic.Neg",
                                       ifelse(grepl("ChronicTcf1-",SampleINFO$title),"Chronic.Neg",
                                       ifelse(grepl("Chronic Tcf1+",SampleINFO$title),"Chronic.Pos",
                                              "Naive" )))))

## make factors of categorical data

CatCols<- c("Donor", "Infection", "Tcf1", "GFP", "Phenotype")
SampleINFO[, CatCols] <- lapply(SampleINFO[, CatCols], as.factor)

### the DGEList (d) has sample info stored in d$Samples. Adding my new phenodata. I dont think this is necessary though.
d_Samples<- d$samples

dSAMPLES2<- cbind(d_Samples, SampleINFO[, c(3:7)])

colnames(dSAMPLES2)

d$samples<- dSAMPLES2


## For microarrays there is usually gene annotation info in fData and the normalised expression values in exprs(gset)
## This RNAseq dataset doesnt have anything here 

## The data in the supp files I downloaded are labelled as processed counts.. not sure how they are processed.
## I tried looking at the data_processing columns for info on if the counts are normalised
## no useful information
Processing<- Phenodata[1, c(18:23)]
Processing<- t(Processing)

Processing

Phenodata[1,40]


## "Raw" data dowloaded to Data folder in my working directory, outside of R. Unzipped the folder to 14 txt.gz files
## in the folder GSE83978_RAW. These files are sepatate "processed" read counts for each sample.

## Get a vector of file names

file_list<- list.files("Data/GSE83978_RAW")


## I had to switch my WD to the "Data/GSE83978_RAW" file to read in the data as I dont know how to specify the path as well as the file_list :/  
## Have a look at one of them if you want. It might mean more to you if you know what RNAseq counts usually look like
## I got DFs with 46522 obs of 2 variables (column 1 = Gene.ID, Column 2= Read counts, headed with a number representing the sample)


Count_DF<- read.table(file_list[1], header=TRUE, sep="\t")

## Handy code I found online to automate merging all of the individual files. 
## Delete, rm(Count_DF), the above Count_DF before running so you dont get the first file twice


for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("Count_DF")){
    Count_DF <- read.table(file, header=TRUE, sep="\t")
  }

 
  
  # if the merged dataset does exist, append to it
  if (exists("Count_DF")){
    temp_dataset <-read.table(file, header=TRUE, sep="\t")
    Count_DF<-merge(Count_DF, temp_dataset)
    rm(temp_dataset)
  }
  
}

### We now have a 15 column dataframe but it has numbers as colnames. 


## Check file list to make sure sample numbers are in the same order as phenodata sample list. 
## File names have Geo asccession followed by the sample number used in the count files
file_list

## pull character vector of geo accessions
SampleIDs <- Phenodata$geo_accession

## Add "Gene.ID" to the vector so its included as the first colname
Count_Colnames <- c("Gene.ID", SampleIDs)

## Replace column names
colnames(Count_DF)<- Count_Colnames






## DE analysis with LIMMA ####
## Code partly from here https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html

BiocManager::install("edgeR")

##Load the edgeR package (which loads limma as a dependency)
library(edgeR)


## Remove summary data from Gene.ID column, keeping only ENSEMBL IDs
library(dplyr)  ## needed for filter

CountsONLY<- filter(Count_DF, grepl(pattern = 'ENSMUS', Gene.ID))


rownames(CountsONLY)<- CountsONLY$Gene.ID
counts_Matrix <- as.matrix(CountsONLY[, -1])

d0 <- DGEList(counts_Matrix)

## Calculate normalization factors

d0 <- calcNormFactors(d0)
d0

dim(d0)  ## 46517    14

##Note: calcNormFactors doesn’t normalize the data, it just calculates normalization factors for use downstream.

## Filter low-expressed genes

cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) # number of genes left = 15027   14

barplot(d$counts) ## not looking normally distributed ?? Dont know if this is needed for RNAseq

boxplot(d$counts, boxwex=0.6, notch=T, outline=FALSE, las=2 )

### Useful info on using deSeq instead of LIMMA here, seems deSeq input is non normalised count matrix so can start from here
## http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#why-un-normalized-counts



### Multidimensional scaling (MDS) plot
### I use this to make sure everything is labelled right and the factors are levelled right. 
### After cheking the sample names are coloured correctly by cross cheking the SampleInfo DF, I add in pch=16,
### This will change the samplenames to dots.

PhenotypeF<- factor(d$samples$Phenotype, levels = unique(d$samples$Phenotype))

PhenotypeF

PhenotypeF_2 <- factor (unique(PhenotypeF)) 


plotMDS(d,main= "GSE83978 MDS", xlab = "Dim1", ylab = "Dim2",col = as.integer(PhenotypeF))
legend("topright", col= as.integer(PhenotypeF_2),
       legend = c("Naive", "Memory Tcf1-","Memory Tcf1+", "Chronic Tcf1-", "Chronic Tcf1+" ) ,pch = 16, cex=0.7, inset=c(0.02))


## Fit linear model with LIMMA

### Build the design matrix for the linear modelling function.

design <- model.matrix(~0 + PhenotypeF)
design
colnames(design) <- levels(PhenotypeF)
design

## Voom transformation and calculation of variance weights
## Specify the model to be fitted. We do this before using voom since voom uses variances of the model residuals (observed - fitted)

y <- voom(d, design, plot = T)


## Aha, y is a large EList with expression data stored in y$E and all of the sample info I put in d$samples is stored to y$targets
## The expresssion data is now normalised too

boxplot(y$E, boxwex=0.6, notch=T, outline=FALSE, las=2 )

Targets<- y$targets ## Having a look to check

##lmFit fits a linear model using weighted least squares for each gene:
 
fit <- lmFit(y, design)
head(coef(fit))

colnames(fit)


contr <- makeContrasts(Memory.Neg = Memory.Neg - Naive, 
                       Memory.Pos= Memory.Pos - Naive,
                       Chronic.Neg= Chronic.Neg - Naive,
                       Chronic.Pos= Chronic.Pos - Naive, 
                       levels = design)
contr

##Estimate contrast for each gene


fit2 <- contrasts.fit(fit, contr)

colnames(fit2)

head(fit2$coefficients, 3) # After contrasts applied, coefs of fit show the contrasts as colNames

#Empirical Bayes smoothing of standard errors (shrinks standard errors that are much larger or smaller 
## than those from other genes towards the average standard error) compute statistics and table of top significant genes

fit2 <- eBayes(fit2)

colnames(fit2)

top.table <- topTable(fit2, sort.by = "F", n = Inf)
head(top.table, 20)

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05)


# Venn diagram of results
vennDiagram(dT, circle.col=palette())



# look at individual contrasts
tmp <- contrasts.fit(fit, coef = 1) # Test Memory neg vs Naive
tmp <- eBayes(tmp)
tT_1 <- topTable(tmp, sort.by = "P", n = Inf)
head(tT_1, 20)
dim(tT_1)   ### 15027   6

length(which(tT_1$adj.P.Val < 0.05)) ## 12073 genes sig DE

colnames(tT_1)
tT_1<- filter(tT_1, adj.P.Val < 0.05)  ## 12073  6
dim(tT_1) 

tT_1<- filter(tT_1, logFC <= -2  |  tT_1, logFC >= 2)

colnames(tT_1) ## LFC1 = 10484 ### alot still. LFC2 = 9517 genes


# volcano plot (log P-value vs log fold change)


volcanoplot(fit2, coef=1, main=colnames(fit2)[1], pch=20,
            highlight=length(which(dT[,1]!=0)), names=rep('+', nrow(fit2))) 





### Making annotation file. Copied from here:https://sbc.shef.ac.uk/workshops/2020-02-13-rnaseq-r/rna-seq-annotation-visualisation.nb.html#adding_annotation_to_the_deseq2_results
### Only execute if you need to install the package

install.packages("BiocManager")
BiocManager::install("org.Mm.eg.db") ## Is this is the right mouse annotation file to use? 

library(org.Mm.eg.db)
columns(org.Mm.eg.db)

## We need to filter the database by a key or set of keys in order to extract the information we want. 
## Valid names for the key can be retrieved with the keytypes function.

keytypes(org.Mm.eg.db)

### We should see ENSEMBL, which is the type of key we are going to use in this case. 
## If we are unsure what values are acceptable for the key, we can check what keys are valid with keys

keys(org.Mm.eg.db, keytype="ENSEMBL")[1:10]


anno <- AnnotationDbi::select(org.Mm.eg.db,keys= rownames(tT_1),
                              columns=c("SYMBOL","GENENAME"),
                              keytype="ENSEMBL")

## Error: select() returned 1:many mapping between keys and columns
## anno longer than top.table, dublicate values? Need to look this up
## Need to merge this with the tT file






top.tableAnn <- merge(top.table, anno, by.x= rownames, by.y= "ENSEMBL")
