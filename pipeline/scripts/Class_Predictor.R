#Make the plot for the Grant based on some of the data from Paul

#First get those libraries
library(edgeR)
library(limma)
library(gplots)
library(statmod)
library(RColorBrewer)
library(randomForest)

args = commandArgs(trailingOnly = TRUE)
fpkm_file=args[1]
rf_model=args[2]
topGene_file=args[3]
counts_file=args[4]
csv_outfile =args[5]
pdf_outfile =args[6]

#Read in FPKM for PH like and non PH like
#N.B. dont skip the first row since it includes the headers
normalisedcounts= read.csv(file=fpkm_file,sep="\t",stringsAsFactors=FALSE,row.names=1)

#Just take the sample columns, leave out the other annotational columns
ncounts = normalisedcounts[,1:289]

###############################################
##### Get Pauls Data ##########################
###############################################

#First read in the file
sample = read.csv(file=counts_file, sep="\t",stringsAsFactors=FALSE,skip=1,row.names=1)

#Extract the gene lengths
glengths = sample[,"Length",drop=FALSE]

#Just keep counts
sN = sample[,-(1:5),drop=FALSE]
#Divide every element by 2 to account for the fact that we are using double counted paired end reads (see README in PaulsDara/20150911 for why)
#sN=sN/2

#Combine samples
paul.counts = sN


#Produce DGEList
raw <- DGEList(counts=paul.counts, genes=rownames(paul.counts))
#### remove genes with very low counts

dat <- raw
dat = calcNormFactors(dat) #Normalise the raw counts
dat <- as.matrix(dat) #Transform to a matrix

#Convert to fpkm using edgeR function rpkm
gg = glengths[row.names(dat),"Length",drop=FALSE] #Select only the genelengths which are left after removing low counts
FPKM <- rpkm(dat,log=TRUE,gene.length = gg$Length)

###########################################################
##### Time to predict Pauls Data ##########################
###########################################################

#Now read in the randomForest classifier
load(rf_model)
load(topGene_file)

#One Gene name not same in Pauls set as ST.Jude data: "C10orf26" -> 
row.names(FPKM)[row.names(FPKM) == "WBP1L"] = "C10orf26"

#Check if we have all genes in FPKM that we need for the classifier
for( i in 1:42)
{
  #print(row.names(topGeneNames)[i])
  FPKM[row.names(topGeneNames)[i],1]
}
#Just take the top genes
#sort(row.names(FPKM)) #Check names alphabetically

#One Gene name not same in Pauls set as ST.Jude data: "C10orf26" -> 
row.names(FPKM)[row.names(FPKM) == "WBP1L"] = "C10orf26"

toclass = FPKM[row.names(topGeneNames),,drop=FALSE]

#Now lets predict
RF2pred <- predict(RF2,data.frame(t(toclass)))
RF2prob <- predict(RF2,data.frame(t(toclass)),type="prob")
#Now lets remove all samples which have votes below threshold
Maj = apply(RF2prob[,1:4,drop=FALSE],1,max,simplify=FALSE)
isThresh = Maj > 0.66
toclass_filter = toclass[,isThresh,drop=FALSE]
PD_ncol=ncol(toclass_filter)

#########################################################
#### Subselect St.Jude Data with genes for classifiers ##
#########################################################

nsubcounts = ncounts[row.names(topGeneNames),]
STJ_ncol = ncol(nsubcounts)

##### Combine Pauls Data and St.Judes
all.counts = cbind(toclass_filter,nsubcounts)


#Figure out the classes
#Make the groupings from Phlike and non phlike
pheno = strsplit2(colnames(all.counts),".Aligned",fixed=TRUE)[,1]
#Change underscore string to dots
pheno = gsub("___________",".",pheno)
pheno = gsub("__________",".",pheno)
pheno = gsub("_________",".",pheno)
pheno = gsub("________",".",pheno)
pheno = gsub("_______",".",pheno)

tmp = strsplit2(pheno,".",fixed=TRUE) #Split the string by .
#Make chained ifelse statement
targets = data.frame(Fusion = ifelse(grepl("Phlike", tmp[,1]) | grepl("PHALL", tmp[,1]), "Phlike", 
                                     ifelse(grepl("ERG", tmp[,1]), "ERG",
                                            ifelse(grepl("ETV", tmp[,1]),"ETV","Other")))
                     
)

#Filter the counts to remove ones with sum 1 or less
isexpr <- rowSums(all.counts) > 1
datExpr <- all.counts[isexpr,]


###############################################
#### Visualisation and Exploratory phase ######
###############################################

#Define a new colour palette with enough colours for the nine different categories
pal = brewer.pal(9,"Set1")
#Define marker set as a vector
pchvec <- c(1:5)

#First 2 dimensions
pchvec = c(rep(16,PD_ncol),rep(1,STJ_ncol))
cexvec = c(rep(1.2,PD_ncol),rep(0.8,STJ_ncol))

#Alternative sizing
#Size of dots proportional to confidence
Maj = Maj*2
cexvec = c(Maj[isThresh],rep(1,STJ_ncol))
write.table(data.frame(RF2pred, RF2prob),csv_outfile,sep=",",row.names=T,col.names=T,quote=F)

pdf(file=pdf_outfile)
if (PD_ncol > 0) {
targets$Fusion[1:PD_ncol] = RF2pred[isThresh,drop=FALSE] #Assign Pauls data as there predicted values from RF


test =plotMDS(datExpr,gene.selection="common",cex=cexvec,col=pal[targets$Fusion],pch=pchvec)



#Now annotate on sample names for Pauls Data
eTruth = c("MLM18: Pre B-ALL","MLM19: JMML","MLM20: pre B-ALL","MLM21: pre B-ALL",
           "MLM24: AML","MLM1: ETV", "MLM5: pre B-ALL", "MLM6: BCL-ABL1", 
           "MLM8: ETV6-RUNX1", "PE2: MLL-ITSN1","PE3:Acute Unclassifiable",
           "MLM10:BCR-ABL1","MLM12:MLL-TFE3", 
           "MLM15:CBFB-MYH11", "MLM16: Pax5-SPAG8", "EKL3-11: KMT2A:MLLT3",
           "EKL3-12: Failed","EKL3-13.2: Unknown","EKL3-13: Unknown","EKL3-14: P2YR8-CRLF2", 
           "EKL3-17: CDKN2A?","EKL3-1: TALL", "EKL3-21: PAX5:JAK2", "EKL3-22: Unkown","EKL3-23: Unkown", "EKL3-25: Unkown",
           "EKL3-26: Pax5:c20orf112", "EKL3-2: ETV6-RUNX1",
           "EKL3-5: Unkown" )

#text(test$x[c(1:29)],test$y[c(1:29)]+0.1,labels = eTruth,cex=0.75)
#legend('topleft',legend=unique(levels(targets$Fusion)),fill=pal)
legend('topleft',legend=c("ERG","ETV","Other","Ph-like"),fill=pal)
legend('bottomleft',legend=c("St.Judes (True Class)","MCRI (Predicted Class)"),pch=c(1,16))
}
dev.off()

