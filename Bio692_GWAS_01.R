library(RColorBrewer)
library(SNPRelate)
library(devtools)
library(gplots)
library(genetics)
library(ape)
library(EMMREML)
library(compiler) 
library("scatterplot3d")
library(LDheatmap)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(tidyr)
library(ggrepel)

source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/GAPIT/emma.txt")

# set the correct working directory
setwd("~/data/...../")

#load the phenotype file
Pheno <-read.table("Whealbi_Pheno_file-1.csv", header=T) 
head(Pheno)
summary(Pheno)


# examine phenotypic data +++++++++++++++++++++++++++++

# make a histogram how phenotype values are distributed for an example accession/isolate
plot(Pheno$Bgt001, pch=20)
hist(Pheno$Bgt001, ylim=c(0,1700))


# reshape the dataframe to long format:
# ... for a few selected isolates
Pheno_long <- pivot_longer(Pheno, cols = c("Bgt001", "Bgt002", "Bgt003"), names_to = "Isolate", values_to = "Infection_score")

# plot with hard coded colors
ggplot(Pheno_long,aes(x=Infection_score, fill=Isolate)) +
  geom_histogram(color="#e9ecef", alpha=1, position = 'stack', bins = 10) +
  scale_fill_manual(values=c("blue", "orange","brown")) +
  theme_classic()

# print bars side-by-side
ggplot(Pheno_long,aes(x=Infection_score, fill=Isolate, color=Isolate)) +
  geom_histogram(color="#e9ecef", alpha=1, position = 'dodge', bins = 10) +
  scale_fill_manual(values=colors) +
  scale_color_manual(values="black") +
    theme_classic()


# now make plots with a longer list of isolates
# ... or make your own color space for longer lists
colors <- c("#0000FF", "orange","brown","#999999", "#00ff00", "#ff7700")

# ... and do the plot using your colors




# ... and the same for all isolates
# hint: extract isolate names, load into a vector and use it as variable

isolates <- colnames(Pheno)[3:dim(Pheno)[2]]





# examine genotype data +++++++++++++++++++++++++++++++++++++++++++++

#load the genotype file
hapmap <- read.table("Whealbi_Geno_File.hapmap", header = T)

# check how many missing data there is in the genotype file
NA_df <- data.frame(Accession=NA, NANr = NA)
for (i in 12:dim(hapmap)[2]){
  na_nr <- sum(is.na(hapmap[,i]), na.rm=TRUE)
  name <- colnames(hapmap)[i]
  NA_df[i,"Accession"] <- name
  NA_df[i,"NANr"] <- na_nr
}

NA_df <- na.omit(NA_df)

# plot the results with ggplot
ggplot(NA_df, aes(x=Accession, y=NANr)) +
  geom_point() +
  theme_classic() +
  theme(axis.text.x = element_blank())


# ... now it would be nicer to order accessions by amount of missing data
# order by number of NAs in each isolate/accession


# what does this tell us?



#check if the SNPs are evenly distributed
ggplot(hapmap, aes(x=pos,y=chrom)) +
  geom_point(shape=124, size=6, alpha=0.1) +
  theme_void() +
  theme(axis.text.y = element_text()) 

## add x axis labels to the plot 
# hint: refine theme by adding x-axis labels and tick etc.



#PCA of genotype data ################################################333 

# load data from vcf 
vcf.fn <- "Whealbi_Geno_File.vcf"
vcf.fn


# Transform the VCF file into gds (genomic data structure) format
showfile.gds(closeall=TRUE)

snpgdsVCF2GDS(vcf.fn, "Genotypes_ALL.gds", method="biallelic.only")
?snpgdsVCF2GDS

# open the newly created gds file
genofile <- snpgdsOpen('Genotypes_ALL.gds')
head(genofile)
genofile


# Perform PCA on the SNP data
pca <- snpgdsPCA(genofile, num.thread=4,autosome.only=FALSE )
?snpgdsPCA
head(pca)

pca_data <- data.frame(sample.id = pca$sample.id,
                       EV1 = pca$eigenvect[,1], # the first eigenvector
                       EV2 = pca$eigenvect[,2], # the second eigenvector
                       EV3 = pca$eigenvect[,3], # the third eigenvector
                       EV4 = pca$eigenvect[,4], # the fourth eigenvector
                       stringsAsFactors = FALSE)
head(pca_data)

# Plot results from PCA analysis of variant data
ggplot(pca_data,aes(x=EV1,y=EV2,label=sample.id)) +
  geom_point(size=0.2)+ theme_classic()

# add labels to some of the accessions/isolates
# since there are many dots, play with the max.overlaps parameter





### GWAS #####################################################
# now we can do GWAS for selected isolates/accessions


# Important: start out in home directory for each isolate
setwd("~/data/...")


# examine the data for individual isolates, search for those that have nice peaks

# example for Bgt001
acc <- "Bgt002"
dir.create(c(acc))
setwd(c(acc))
getwd()
head(Pheno)

myacc <- Pheno[,c(1,4)] ### select the appropriate phenotyping column
head(myacc)
myG <- read.table("../Whealbi_Geno_File.hapmap",header=F)
myG
GAPIT_JIW2 <- GAPIT(Y=myacc,G=myG) ### command to run the GWAS, this may take a while

# define name of input file to load GWAS result
infile <- paste("GAPIT.Association.GWAS_Results.MLM.",acc,".csv", sep="")
#GWAS <- read.table("GAPIT.Association.GWAS_Results.MLM.Bgt07230.csv", sep = ",", header = T)

GWAS <- read.table(infile, sep = ",", header = T)



# display unwrapped, all in one line
ggplot(GWAS, aes(x=Pos, y=-log10(P.value))) +
  geom_point(size=0.3, color="Grey40") +
  facet_wrap(~Chr, nrow = 1) +
  theme_classic()

# display wrapped
ggplot(GWAS, aes(x=Pos, y=-log10(P.value))) +
  geom_point(size=0.3, color="Grey40") +
  facet_wrap(~Chr) +
  theme_classic()



# calculate the bonferroni correction to set the sig. threshold
alpha = 0.05 # this is the normal significance threshold
test_number = dim(GWAS)[1] # this is the number of SNPs in your GWAS
bonf = alpha/test_number


# plot the GWAS results again, but this time add the bonferroni threshold:
ggplot(GWAS, aes(x=Pos, y=-log10(P.value))) +
  geom_point(size=0.1, color="Grey40") +
  geom_hline(yintercept = -log10(bonf), linetype = "dashed", linewidth = 0.3) +
  facet_wrap(~Chr) +
  theme_classic()

# make the axis more legible
ggplot(GWAS, aes(x=Pos/1000000, y=-log10(P.value))) +
  geom_point(size=0.3, color="Grey40") +
  geom_point(data=GWAS[GWAS$P.value<bonf,], size = 0.6, color="red") +
  geom_hline(yintercept = -log10(bonf), linetype = "dashed", size = 0.3) +
  scale_x_continuous(name="Position [Mb]") +
  facet_wrap(~Chr) +
  theme_classic() +
  theme(
    text = element_text(size=15)
  )


# for thos isolates/genes that give nice peaks:
## create a plot with only the chromosomes showing significant SNPs
## hint you can use a similar way of filtering the data as we did before
## when we only color significant SNPs



## Zoom in only on the area surrounding the best SNPs
## there is a little help where to put your numbers below







