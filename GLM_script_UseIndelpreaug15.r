#!/usr/bin/env Rscript
#   PFA_VMEM="200000"  # this is used for glm reports, needed for srun
# DT[, c("z","u","v"):=NULL] #remove several columns at once

## The heart of the GLM. Uses the text files of ids generated by the naive method run
# to assign reads to categories and outputs predictions per junction to glmReports. models
# are saved into glmModels for further manual investigation.

########## FUNCTIONS ##########

require(data.table)
library(base)
set.seed(1, kind = NULL, normal.kind = NULL)

# allows for variable read length (for trimmed reads)
getOverlapForTrimmed <- function(x, juncMidpoint=150){
    if (as.numeric(x["pos"]) > juncMidpoint){
      overlap = 0
    } else if (as.numeric(x["pos"]) + as.numeric(x["readLen"]) - 1 < juncMidpoint + 1){
      overlap = 0
    } else {
      overlap = min(as.numeric(x["pos"]) + as.numeric(x["readLen"]) - juncMidpoint,
      juncMidpoint + 1 - as.numeric(x["pos"]))
    }
  
  return(overlap)
}

processScoreInput <- function(scoreFile){
  setkey(scores, id)
    return(scores)
}

addDerivedFields <- function(dt, useClass){
# correction 4/2016 of class files
if(nrow(dt) > 0){
    # calculate and add on cols for junction overlap, score adjusted for N penalties, 
    ########## now, we have info for read1 and read2 
    # calculate and add on cols for junction overlap, score adjusted for N penalties,
    dt[,`:=`(is.pos=useClass,overlap=apply(dt, 1, getOverlapForTrimmed))]  # syntax for multiple :=
    # and length-adjusted alignment score (laplace smoothing so alignment score of 0 treated different for different length reads)
    dt[, lenAdjScore:=(as.numeric(aScore) - 0.001)/as.numeric(readLen)]
    dt[,`:=`(pos=NULL, aScore=NULL, numN=NULL, readLen=NULL)]
################# repeat for read2
 ## therefore, only add length adjusted alignment score for R2 !!
    # and length-adjusted alignment score (`` smoothing" so alignment score of 0 treated different for different length reads)
    dt[, lenAdjScoreR2:=(aScoreR2 - 0.001)/readLenR2]
    dt[,`:=`(pos=NULL, aScoreR2=NULL, numNR2=NULL, readLenR2=NULL, adjScoreR2=NULL, aScore=NULL, numN=NULL, readLen=NULL, adjScore=NULL)]
}
    return(dt)
}

# the input file is just the file output by the circularRNApipeline under /ids
processClassInput <- function(classFile,my.names){

#cats = fread(classFile,  sep="\t", nrows=1000000)
cats = fread(classFile,  sep="\t")
############################################################
if ( my.names!="none"){
names(cats)=my.names
}

# syntax for changing names setnames(cats, names(cats), c("id", "R1", "R2", "class"))
  setkey(cats, id)  
  return(cats)
}

# To avoid integer underflow issue when we have too many very small or very large probabilities.
# Take inverse of posterior probability, then take log, which simplifies to sum(log(q) - /sum(log(p))
# and then reverse operations to convert answer back to a probability.
# param p: vector of p values for all reads aligning to junction
# return posterior probability that this is a circular junction based on all reads aligned to it
getPprodByJunction <- function(p ){
  out = tryCatch(
{
  x = sum(log(p))  # use sum of logs to avoid integer underflow
  return(exp(x))
},
error = function(cond){
  print(cond)
  print(p)
  return("?")
},
warning = function(cond){
  print(cond)
  print(p)
  return("-")
}
  )
return(out)
}

applyToClass <- function(dt, expr) {
  e = substitute(expr)
  dt[,eval(e),by=is.pos]
}

applyToJunction <- function(dt, expr) {
  e = substitute(expr)
  dt[,eval(e),by=junction]
}


#######################################################################
######################## BEGIN JS ADDITION ############################
####################### FIRST JS FUNCTION #############################
########################################################################
################# JS added function to FIT the GLM using arbitrary two-classes

my.glm.model<-function( linear_reads, decoy_reads,use_R2 , max.iter){
### FUNCTION TO FIT GLM TO linear READS, returns the GLM and junction predictions, 
saves = list()  # to hold all of the glms for future use
#max.iter = 2  # number of iterations updating weights and retraining glm

# set up structure to hold per-read predictions
n.neg = nrow(decoy_reads) 
n.pos = nrow(linear_reads)
n.reads = n.neg+n.pos
class.weight = min(n.pos, n.neg)

## note that this is coded as linear_reads and decoy_reads but applies to any pair class

readPredictions = rbindlist(list(linear_reads, decoy_reads))

# set initial weights uniform for class sum off all weights within any class is equal
if (n.pos >= n.neg){
  readPredictions[,cur_weight:=c(rep(n.neg/n.pos, n.pos), rep(1, n.neg))]
} else {
  readPredictions[,cur_weight:=c(rep(1, n.pos), rep(n.pos/n.neg, n.neg))]
}

# glm
for(i in 1:max.iter){
  # M step: train model based on current read assignments, down-weighting the class with more reasourcds

if (use_R2==1){
  x = glm(is.pos~overlap+lenAdjScore+qual +lenAdjScoreR2 + qualR2, data=readPredictions, family=binomial(link="logit"), weights=readPredictions[,cur_weight])

}
if (use_R2==0){
  x = glm(is.pos~overlap+lenAdjScore+qual , data=readPredictions, family=binomial(link="logit"), weights=readPredictions[,cur_weight])
}
  saves[[i]] = x

  # get CI on the output probabilities and use 95% CI
  preds = predict(x, type = "link", se.fit = TRUE)
  critval = 1.96 # ~ 95% CI
  upr = preds$fit + (critval * preds$se.fit)
  lwr = preds$fit - (critval * preds$se.fit)
  upr2 = x$family$linkinv(upr)
  lwr2 = x$family$linkinv(lwr)
  
  # use the upper 95% value for decoys and lower 95% for linear
  adj_vals = c(rep(NA, n.reads))
  adj_vals[which(readPredictions$is.pos == 1)] = lwr2[which(readPredictions$is.pos == 1)]
  adj_vals[which(readPredictions$is.pos == 0)] = upr2[which(readPredictions$is.pos == 0)]
  x$fitted.values = adj_vals  # so I don't have to modify below code
  
  # report some info about how we did on the training predictions
  totalerr = sum(abs(readPredictions[,is.pos] - round(x$fitted.values)))
  print (paste(i,"total reads:",n.reads))
  print(paste("both negative",sum(abs(readPredictions[,is.pos]+round(x$fitted.values))==0), "out of ", n.neg))
  print(paste("both positive",sum(abs(readPredictions[,is.pos]+round(x$fitted.values))==2), "out of ", n.pos))
  print(paste("classification errors", totalerr, "out of", n.reads, totalerr/n.reads ))
  print(coef(summary(x)))
  readPredictions[, cur_p:=x$fitted.values] # add this round of predictions to the running totals
  
  # calculate junction probabilities based on current read probabilities and add to junction predictions data.table

  tempDT = applyToJunction(subset(readPredictions, is.pos == 1), getPprodByJunction(cur_p))
  setnames(tempDT, "V1", paste("iter", i, sep="_")) # iter_x is the iteration of product of ps
  setkey(tempDT, junction)
  junctionPredictions = junctionPredictions[tempDT]  # join junction predictions and the new posterior probabilities
  rm(tempDT)  # clean up
  
  # E step: weight the reads according to how confident we are in their classification. Only if we are doing another loop
  if(i < max.iter){
    posScale = class.weight/applyToClass(readPredictions,sum(cur_p))[is.pos == 1,V1]
    negScale = class.weight/(n.neg - applyToClass(readPredictions,sum(cur_p))[is.pos == 0,V1])
    readPredictions[is.pos == 1,cur_weight:=cur_p*posScale]
    readPredictions[is.pos == 0,cur_weight:=((1 - cur_p)*negScale)]
  }
  setnames(readPredictions, "cur_p", paste("iter", i, sep="_")) # update names
}  

# calculate mean and variance for null distribution
## this uses a normal approximation which holds only in cases with large numbers of reads, ie the CLT only holds as the number of reads gets very large

## should be called p-predicted
read_pvals = readPredictions[,max.iter]

# rename cols to be consistent with circular glmReports, syntax below removes col. "ITER_1"
if (max.iter>1){
# cleaning up
for (myi in c(1:(max.iter-1))){
junctionPredictions[, paste("iter_",myi,sep=""):=NULL]
}
}
setnames(junctionPredictions, paste("iter_",max.iter,sep=""), "p_predicted")
list(saves, junctionPredictions) ## JS these are the outputs and done with function
}


########################################################################
###################### prediction from model ##########################
##### as a function, needs input data and model

predictNewClassP <- function(my_reads, null){ ## need not be circ_reads, just easier syntax
######### up until this point, every calculation is PER READ, now we want a function to collapse 
######### want to do hypothesis testing 
# calculate junction probabilities based on predicted read probabilities
## Use simple function-- NOTE: "p predicted" is a CI bound not the point estimate. It is still technically a consistent estimate of p predicted 
## prob of an anomaly by glm is phat/(1+phat) under 'real' 1/(1+phat) under 'decoy' junction, so the ratio of these two reduces to 1/phat. as phat -> 1, no penalty is placed on anomaly.

#merge
junctionPredictions = my_reads[, .N, by = junction] # get number of reads per junction
setnames(junctionPredictions, "N", "numReads")
setkey(junctionPredictions, junction)

my_reads[, logproduct:=sum( log (p_predicted) * (1-is.anomaly) + log( 1/(1+p_predicted) *is.anomaly)), by=junction]

## is anomaly adjusted log sum scoremm
logsum=my_reads[,sum( log ( p_predicted / (1+p_predicted*is.anomaly))), by=junction]
logsum_2=my_reads[,sum( log ( p_predicted_2 / (1+p_predicted_2*is.anomaly))), by=junction]
print ("Logsum is reported which is equal to the sum of the logs of phats-- if exponentiated, corresponds to product of ps") 

## merge these new variables to the dataframe
junctionPredictions=merge(junctionPredictions,logsum)
setnames(junctionPredictions, "V1", "logsum")

junctionPredictions=merge(junctionPredictions,logsum_2)
setnames(junctionPredictions, "V1", "logsum_2")

print (names(junctionPredictions))

########### adding quantiles of p_predicted
n.quant=2
for (qi in 1:n.quant){
my_quantiles = my_reads[,round(10*quantile(p_predicted/(1+is.anomaly* p_predicted),probs=c(0:n.quant)/n.quant,na.rm=TRUE)[qi])/10,by=junction]

# merge into junctionPredictions
print (head(my_quantiles))
setkey(my_quantiles,junction)
junctionPredictions=merge(junctionPredictions,my_quantiles)
setnames(junctionPredictions, "V1", paste("q_",qi,sep=""))
}

##################################
##  tempDT, to collapse across junctions 
# p_predicteds are the exponentiation
junctionPredictions [ ,p_predicted_2:=exp(logsum_2),by=junction]
junctionPredictions [ ,p_predicted:=exp(logsum),by=junction]

print (head(junctionPredictions[order(junction),]))

## NOTE: P VALUE IS probability of observing a posterior as extreme as it is, "getPvaluebyJunction" is a bayesian posterior
junctionPredictionsWP=assignP(junctionPredictions,null) 

rm(tempDT)  # clean up
## adding here:

unique(junctionPredictionsWP) ## returned

}
########################################################################################### ASSIGN p values through permutation
################################### 
assignP<-function(junctionPredictions,null) {
# logsum is the logged sum
# add p-value to junctionPredictions (see GB supplement with logic for this)

lognull=log(null)

use_mu = mean(lognull) # this is actually the mean of the read level predictions
use_var=var(lognull)
## for large n, 
#print ("using cdf of null distribution as "p_value" which is misnomer for convenient and replaced below ")
n.thresh.exact=15
print (n.thresh.exact)

junctionPredictions[ (numReads>n.thresh.exact) , p_value :=  pnorm((logsum - numReads*use_mu)/sqrt(numReads*use_var))]

junctionPredictions[ (numReads>n.thresh.exact) , p_value_2 :=  pnorm((logsum_2 - numReads*use_mu)/sqrt(numReads*use_var))]

## make empirical distribution of posteriors:

print ("exact calculation through sampling 10K p predicted")
my.dist=list(n.thresh.exact)
for ( tempN in 1:n.thresh.exact){ #### get distributions to convolve in next loop
n.sampled=1000 # used to compute the null distribution of posteriors
my.dist[[tempN]]=sample(lognull, n.sampled, replace=T)
}

for ( tempN in 1:n.thresh.exact){ ## use this loop to assign jncts w/ tempN
sim.readps=my.dist[[1]]
if (tempN>1){
for (tj in 2: tempN){ # loop, taking products
sim.readps=my.dist[[tj]] +  sim.readps
}
}
# convert to posterior
## fraction of time p_predicted is smaller than -- so if p_predicted is very large, the fraction of time it is smaller is big
## use the null to compute p_vals
print (head(junctionPredictions))
print (paste(tempN, "is value of readcount for exact calculation and length of sim reads is ",length(sim.readps)))

junctionPredictions [ (numReads == tempN ), p_value:= sum( exp(sim.readps)<p_predicted)/length(sim.readps),by=junction]

print ("if below table is not empty, ERROR")
print (junctionPredictions[p_value>1])


junctionPredictions [ (numReads == tempN ), p_value_2:= sum( exp(sim.readps)<p_predicted_2)/length(sim.readps), by=junction]

}
return(junctionPredictions)
}
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
######## END FUNCTIONS, BEGIN WORK #########

## command line inputs
user.input=0
no.indel.classfile=1

################# USER INPUT SHOULD BE 1 if it is used in an automated script
use.tcga.ov=0
use.uconn=0
use.normal.breast=0
use.cml=0
use.ews=0
use.simulated=0
use.ewsh5=0
use.seqc=0
## GILLIAN


########## Gillian, note the 'indels vs indel' in the normal breast names and the 'no 1' ######################
if (use.normal.breast==1){
parentdir="/scratch/PI/horence/gillian/normal_breast/circpipe/circReads/ids/"
srr="SRR1027190"
class_input = paste(parentdir,srr,"_1__output.txt",sep="")
reg_indel_class_input = paste(parentdir,srr,"_output_RegIndel.txt",sep="")

parentdir="/scratch/PI/horence/gillian/normal_breast/FarJunc_Mar9/GLM_classInput/"
FJ_indel_class_input = paste(parentdir,srr,"_output_FJIndels.txt",sep="")
fusion_class_input = paste(parentdir,srr,"_1__output_FJ.txt",sep="")
output_dir=""
sampletest=paste("normal_breast_windel_",srr,sep="")
}

if (use.tcga.ov==1){
parentdir="/srv/gsfs0/scratch/linda/TCGA_OV_Apr7/circpipe/circReads/ids/"
srr="C09BFACXX_1_AAATGC"
class_input = paste(parentdir,srr,"_1__output.txt",sep="")
reg_indel_class_input = paste(parentdir,srr,"_output_RegIndel.txt",sep="")

parentdir=" /srv/gsfs0/scratch/linda/TCGA_OV_Apr7/Machete/GLM_classInput/"

FJ_indel_class_input = paste(parentdir,srr,"_output_FJIndels.txt",sep="")
fusion_class_input = paste(parentdir,srr,"_1__output_FJ.txt",sep="")
output_dir=""
sampletest=paste("tcga_ov_",srr,sep="")
}



if (user.input==0){
args = commandArgs(trailingOnly = TRUE)
fusion_class_input=args[1]
class_input=args[2]
srr= args[3]
output_dir=args[4]
reg_indel_class_input = args[5]
FJ_indel_class_input = args[6]
print (paste("outputdir", output_dir))
## should be:
#FJ_indel_class_input = paste(parentdir,srr,"_output_FJIndels.txt",sep="")
#reg_indel_class_input = paste(parentdir,srr,"_output_RegIndel.txt",sep="")
}

max.iter=2 ## iterations for glm

if (user.input==1){
  output_dir=""}

use.indels=1
use.fusion=1

if (user.input==1){ ## THIS IF LOOP ONLY if we want to bypass command line inputs!

## GILLIAN
if(use.ewsh5==1){
parentdir="/scratch/PI/horence/alignments/EWS_FLI_bigmem/circReads/ids/"
srr="H5-AD018"
fusion_class_input = paste(parentdir,srr,"_1__output_FJ.txt",sep="")
class_input = paste(parentdir,srr,"_S9_L2L3_R1__output.txt",sep="")

FJ_indel_class_input = paste(parentdir,srr,"_output_FJIndels.txt",sep="")
reg_indel_class_input = paste(parentdir,srr,"_output_RegIndel.txt",sep="")
output_dir=""
sampletest="ews-spork"

}


if(use.simulated==1){
parentdir="/scratch/PI/horence/gillian/Engstrom/circpipe_engstrom/circReads/ids/"
srr="Human_simulated_reads2"
fusion_class_input = paste(parentdir,srr,"_1__output_FJ.txt",sep="")
class_input = paste(parentdir,srr,"_1__output.txt",sep="")

FJ_indel_class_input = paste(parentdir,srr,"_output_FJIndels.txt",sep="")
reg_indel_class_input = paste(parentdir,srr,"_output_RegIndel.txt",sep="")
output_dir=""
sampletest="engstrom1"
}
if (use.cml==1){
parentdir=""
cmli=1
srr=paste("ENCFF000HOC",cmli,sep="")
parentdir="/scratch/PI/horence/gillian/CML_test/aligned/CML/circReads/ids/"
fusion_class_input = paste(parentdir,srr,"_1__output_FJ.txt",sep="")
class_input = paste(parentdir,srr,"_1__output.txt",sep="")

FJ_indel_class_input = paste(parentdir,srr,"_output_FJIndels.txt",sep="")
reg_indel_class_input = paste(parentdir,srr,"_output_RegIndel.txt",sep="")
output_dir=""
sampletest=paste("CML",cmli,sep="")
}
if (use.ews==1){
srr="SRR1594025" #
srr="SRR1594022"
parentdir="/scratch/PI/horence/gillian/Ewing/circpipe/circReads/ids/"
fusion_class_input = paste(parentdir,srr,"_1__output_FJ.txt",sep="")
class_input = paste(parentdir,srr,"_1__output.txt",sep="")

FJ_indel_class_input = paste(parentdir,srr,"_output_FJIndels.txt",sep="")
reg_indel_class_input = paste(parentdir,srr,"_output_RegIndel.txt",sep="")
output_dir=""

output_dir=""
sampletest=paste("ews_",srr,sep="")

}

########## Gillian, note the 'indels vs indel' in the normal breast names and the 'no 1' ######################
if (use.normal.breast==1){
parentdir="/scratch/PI/horence/gillian/normal_breast/circpipe/circReads/ids/"
srr="SRR1027190"
class_input = paste(parentdir,srr,"_1__output.txt",sep="")
fusion_class_input = paste(parentdir,srr,"_1__output_FJ.txt",sep="")
#FJ_indel_class_input = paste(parentdir,srr,"_output_FJIndels.txt",sep="")
reg_indel_class_input = paste(parentdir,srr,"_output_RegIndel.txt",sep="")
output_dir=""
sampletest=paste("normal_breast_windel_",srr,sep="")
}
}

## define output class files
glm_out = paste(output_dir,srr,"_DATAOUT",sep="")
anomaly_glm_out = paste(output_dir,srr,"_AnomalyDATAOUT",sep="")
indel_glm_out = paste(output_dir,srr,"_IndelDATAOUT",sep="")
linear_juncp_out = paste(output_dir,srr,"_LINEARJUNCP_OUT",sep="")
circ_juncp_out = paste(output_dir,srr,"_CIRC_JUNCP_OUT",sep="")
fusion_juncp_out = paste(output_dir,srr,"_FUSION_JUNCP_OUT",sep="")
linearwanomaly_juncp_out = paste(output_dir,srr,"_LINEAR_W_ANOMALY_JUNCP_OUT",sep="")
fusionwanomaly_juncp_out= paste(output_dir,srr,"_FUSION_W_ANOM_JUNCPOUT", sep="")
fusionwanomaly_and_indel_juncp_out= paste(output_dir,srr,"_FUSION_W_ANOM_AND_INDEL_JUNCPOUT", sep="")
####### CURRENT GOAL IS TO UPDATE LINEAR pvalues w/ anomaly reads
linear_juncp_update_out = paste(output_dir,srr,"_LINEARJUNCP_UPDATED_OUT",sep="")
##### SINK THE SUMMARY OF THIS SCRIPT
print (paste("OUtput dir", fusionwanomaly_and_indel_juncp_out))
#sink(paste(output_dir,srr,"_glmInformation_",sep=""))

### DONE W/ LOOP


my.names="none" ## this is bc Gillians fields are not names like Lindas are
myClasses = processClassInput(class_input, my.names)

if (use.fusion==1){
print (paste("using ", fusion_class_input))
myClassesFusion = processClassInput(fusion_class_input,names(myClasses))
}
if (use.indels==1){
print (paste("using ", reg_indel_class_input," and " , FJ_indel_class_input))
myClassesRegIndel = processClassInput(reg_indel_class_input,names(myClasses))
myClassesFJIndel = processClassInput(FJ_indel_class_input,names(myClasses))
}

print(paste("class info processed", dim(myClasses)))

circ_reads = myClasses[(tolower(class) %like% 'circ'), list(id, pos, qual, aScore, numN, readLen, junction, qualR2,aScoreR2, numNR2, readLenR2),]
circ_reads = addDerivedFields(circ_reads, 1)
circ_reads [, is.anomaly:=0] ## this is not an anomaly type so WILL NOT have p value ajustment

print ("finished circ_reads")

decoy_reads = myClasses[(tolower(class) %like% 'decoy'), list(id, pos, qual, aScore, numN, readLen, junction, qualR2,aScoreR2, numNR2, readLenR2),]
decoy_reads = addDerivedFields(decoy_reads, 0)
decoy_reads [, is.anomaly:=1] ######## this IS an anomaly type 

print ("finished decoy_reads")
## was
linear_reads = myClasses[(tolower(class) %like% 'linear'), list(id, pos, qual, aScore, numN, readLen, junction, qualR2,aScoreR2, numNR2, readLenR2),]
linear_reads = addDerivedFields(linear_reads, 1)
linear_reads [, is.anomaly:=0] ## this is not an anomaly type so WILL NOT have p value ajustment

print ("finished linear_reads")

anomaly_reads = myClasses[(tolower(class) %like% 'anomaly'), list(id, pos, qual, aScore, numN, readLen, junction, qualR2,aScoreR2, numNR2, readLenR2),]
anomaly_reads [, is.anomaly:=1]
print ("finished anomaly_reads")

if (use.indels==1){
## in analogy, we first define all indels, then assign good and bad
## we will use anomaly field as a general term for 'anomaly mapping and indel'
reg_indel_reads = myClassesRegIndel[, list(id, pos, qual, aScore, numN, readLen, junction, qualR2,aScoreR2, numNR2, readLenR2),]
reg_indel_reads [, is.anomaly:=1]
print ("Finished reg indels")

FJ_indel_reads = myClassesFJIndel[, list(id, pos, qual, aScore, numN, readLen, junction, qualR2,aScoreR2, numNR2, readLenR2),]
FJ_indel_reads [, is.anomaly:=1]
FJ_indel_reads = addDerivedFields(FJ_indel_reads, 1)
print ("Finished FJ indels")
}

if (use.fusion==1){

#was na.omit
fusion_reads = (myClassesFusion[(tolower(class) %like% 'fjgood'), list(id, pos, qual, aScore, numN, readLen, junction, qualR2,aScoreR2, numNR2, readLenR2),])
fusion_reads = addDerivedFields(fusion_reads, 1)
fusion_reads [, is.anomaly:=0] ## this is not an anomaly type so WILL NOT have p value ajustment

print ("ANOMALY fusions defined as FJ bad-- of any variety")
anomaly_fusion_reads = (myClassesFusion[(tolower(class) %like% 'fjbad'), list(id, pos, qual, aScore, numN, readLen, junction, qualR2,aScoreR2, numNR2, readLenR2),]) # GILLIAN, pls comment here on what FJ bad is for the sake of documentation

anomaly_fusion_reads = addDerivedFields(anomaly_fusion_reads, 1)
anomaly_fusion_reads [, is.anomaly:=1] ## this is not an anomaly type so WILL NOT have p value ajustment
}
###############################################################################################
## CANNOT ADD DERIVED FIELDS HERE BECAUSE WE DON'T KNOW WHICH ANOMALIES ARE GOOD AND/OR BAD
##################### DERIVED FIELDS ADDED LATER ##############################################
###############################################################################################

# set up data structure to hold per-junction predictions
junctionPredictions = linear_reads[, .N, by = junction] # get number of reads per junction
setnames(junctionPredictions, "N", "numReads")
setkey(junctionPredictions, junction)

#### TRAIN EM ####
## this should be a function of any two classes; and the output will be the model

## 

n.row= dim(linear_reads)[1]
n.sample=min(n.row,10000) 

#syntax example decoy_reads[,p_predicted:=NULL]
print ("calling linear decoy model")
linearDecoyGLMoutput = my.glm.model ( linear_reads[ sample(n.row,n.sample,replace=FALSE),], decoy_reads, 1, max.iter) ## 0 does not use R2 info 

saves = linearDecoyGLMoutput[[1]]
linearJunctionPredictions =  linearDecoyGLMoutput[[2]]
#save(saves, file=glm_out)  # save models
linearDecoyGLM = saves[[max.iter]] ##### this is the glm model

## after fitting the GLM to linear vs. decoy, we want to store linear junction predictions in order to subset anomalies
######## JS ADDITION: NOTE- NOT stratifying on permutation p value, although could add this too


############################################################################
### START LINEARS
################# predict on anomaly reads -- AND TEST HOW THIS IMPACTS LINEAR PREDICTIONS
############################ linear predictions ONLY ON THE BASIS of anomalies...

preds = predict(linearDecoyGLM, newdata=linear_reads, type = "link", se.fit = TRUE) 
lwr = preds$fit - (1.96 * preds$se.fit)  # ~ lower 95% CI to be conservative 
lwr_2 = preds$fit - 2*(1.96 * preds$se.fit)  # ~ lower 99% CI to be conservative 
linear_reads[, p_predicted:= linearDecoyGLM$family$linkinv(lwr)] # add lower 95% CI prediction
linear_reads[, p_predicted_2:= linearDecoyGLM$family$linkinv(lwr_2)] # add lower 95% CI prediction
## for null
print ("Assigning null distribution for all linear reads")
null=linear_reads$p_predicted

### ASSIGN p value: 
linearJunctionPredictionsForModels = predictNewClassP(linear_reads, null)

pGoodThresh=quantile(linearJunctionPredictionsForModels$p_value,prob=.8,na.rm=TRUE)    
good.linear=linearJunctionPredictionsForModels[p_value> pGoodThresh,]
    
pBadThresh=quantile(linearJunctionPredictionsForModels$p_value,prob=.2,na.rm=TRUE)    
bad.linear=linearJunctionPredictionsForModels[p_value< pBadThresh,]    
#####################################


############## done with linears

# define two classes of anomalies: those from good vs. bad junctions
good_anomaly_reads= anomaly_reads[!is.na(match(junction, good.linear$junction)),]
bad_anomaly_reads= anomaly_reads[!is.na(match(junction, bad.linear$junction)),]
## NOW add derived fields

if (use.indels==1){
# define two classes of regular INDELS for training: those from good vs. bad junctions
good_indel_reads= reg_indel_reads[!is.na(match(junction, good.linear$junction)),]
bad_indel_reads= reg_indel_reads[!is.na(match(junction, bad.linear$junction)),]

## if either class is empty, change use.indels to zero

if (min(dim(good_indel_reads)[1],dim(bad_indel_reads)[1])<20){ 
print ("CHANGING TO INDELS=0 becuase too few reads")
use.indels=0 
}


}

##### now, re-run script training on anomalies from good vs. bad

## The "1" belowneeds to be 0 and 1 for good and bad anomalies
good_anomaly_reads = addDerivedFields(good_anomaly_reads, 1) #### 
bad_anomaly_reads = addDerivedFields(bad_anomaly_reads, 0) #### 
all_anomaly_reads=rbind(good_anomaly_reads,bad_anomaly_reads)

## in anaolgy, same for INDELS
if (use.indels==1){
good_indel_reads = addDerivedFields(good_indel_reads, 1) #### 
bad_indel_reads = addDerivedFields(bad_indel_reads, 0) #### 
all_indel_reads=rbind(good_indel_reads,bad_indel_reads)
}
######## now, RECALL GLM  FOR ANOMALY READ MAPPERS
print ("calling good anomaly bad anomaly model with .8 and .2 as thresholds")
AnomInfo = my.glm.model (good_anomaly_reads, bad_anomaly_reads,1, max.iter) ## JS these are the outputs and done with function
ANOMALYsaves=AnomInfo[[1]]
ANOMALYjunctionPredictions=AnomInfo[[2]]
AnomalyGLM = ANOMALYsaves[[max.iter]] ##### this is the glm model

if (use.indels==1){
print ("calling INDEL model with .8 and .2 as thresholds")
IndelInfo = my.glm.model (good_indel_reads, bad_indel_reads,1, max.iter) ## JS these are the outputs and done with function
INDELsaves=IndelInfo[[1]]
INDELjunctionPredictions=IndelInfo[[2]]
IndelGLM = INDELsaves[[max.iter]] ##### this is the glm model
}
if (use.indels==0){
print ("USING GLM MODEL for INDEL MODEL")
INDELsaves = linearDecoyGLMoutput[[1]]
INDELlinearJunctionPredictions =  linearDecoyGLMoutput[[2]]
IndelGLM = linearDecoyGLM
}

## save GLMS
#save(AnomalyGLM, file=anomaly_glm_out)  # save models
if (use.indels==1){
save(IndelGLM, file=indel_glm_out)  # save models
}
#save(linearDecoyGLM, file=glm_out)  # save models

preds = predict(linearDecoyGLM, newdata=decoy_reads, type = "link", se.fit = TRUE) 
lwr = preds$fit - (1.96 * preds$se.fit)  # ~ lower 95% CI to be conservative 
lwr_2 = preds$fit - 2*(1.96 * preds$se.fit)  # ~ lower 99% CI to be conservative 
decoy_reads[, p_predicted:= linearDecoyGLM$family$linkinv(lwr)] # add lower 95% CI prediction
decoy_reads[, p_predicted_2:= linearDecoyGLM$family$linkinv(lwr_2)] # add lower 95% CI prediction

preds = predict(AnomalyGLM, newdata=all_anomaly_reads, type = "link", se.fit = TRUE) 

lwr = preds$fit - (1.96 * preds$se.fit)  # ~ lower 95% CI to be conservative 
lwr_2 = preds$fit - 2*(1.96 * preds$se.fit)  # ~ lower 99% CI to be conservative 

all_anomaly_reads[, p_predicted:= AnomalyGLM$family$linkinv(lwr)] # add lower 95% CI prediction
all_anomaly_reads[, p_predicted_2:= AnomalyGLM$family$linkinv(lwr_2)] # add lower 95% CI prediction

## need to rbind anomaly junctions 
linear_and_anomaly=rbind(all_anomaly_reads[,list(qual,lenAdjScore,qualR2,lenAdjScoreR2,junction,is.pos,overlap,is.anomaly,p_predicted,p_predicted_2)], linear_reads[,list(qual,lenAdjScore,qualR2,lenAdjScoreR2,junction,is.pos,overlap,is.anomaly,p_predicted,p_predicted_2)])



linearWithAnomalyJunctionPredictions = predictNewClassP(linear_and_anomaly, null)

################# DONE WITH LINEARS
#########################################################################

#### PREDICT CIRCULAR JUNCTIONS #### SHOULD MAKE THIS MODULAR AND A FUNCTION so Farjunction and Anomalies can be used
## SIMPLE PREDICT ON CIRCLES
preds = predict(linearDecoyGLM, newdata=circ_reads, type = "link", se.fit = TRUE) 
lwr = preds$fit - (1.96 * preds$se.fit)  # ~ lower 95% CI to be conservative 
lwr_2 = preds$fit - 2*(1.96 * preds$se.fit)  # ~ lower 99% CI to be conservative 

circ_reads[, p_predicted:= linearDecoyGLM$family$linkinv(lwr)] # add lower 95% CI prediction
circ_reads[, p_predicted_2:= linearDecoyGLM$family$linkinv(lwr_2)] # add lower 95% CI prediction
circularJunctionPredictions = predictNewClassP(circ_reads, null)

#########################################################################
## start fusions ############################################################################
############################################################################
################# NOTE: Here, we are only using fusion reads not fusion anomaly reads and treating the prediction just like circle
############################################################################
############################################################################
## start prediction on good far junctions

if (use.fusion==1){
fusion_reads$overlap=as.numeric(as.vector(fusion_reads$overlap))
preds = predict(linearDecoyGLM, newdata=fusion_reads, type = "link", se.fit = TRUE) 
lwr = preds$fit - (1.96 * preds$se.fit)  # ~ lower 95% CI to be conservative 
lwr_2 = preds$fit - 2*(1.96 * preds$se.fit)  # ~ lower 99% CI to be conservative 

fusion_reads[, p_predicted:= linearDecoyGLM$family$linkinv(lwr)] # add lower 95% CI prediction
fusion_reads[, p_predicted_2:= linearDecoyGLM$family$linkinv(lwr_2)] # add lower 95% CI prediction
fusionJunctionPredictions = predictNewClassP( fusion_reads, null)


# start prediction on BAD=Anomaly mapping
print ("done with fusion normals, starting anomalies")

preds = predict(AnomalyGLM, newdata=anomaly_fusion_reads, type = "link", se.fit = TRUE) 
lwr = preds$fit - (1.96 * preds$se.fit)  # ~ lower 95% CI to be conservative
lwr_2 = preds$fit - 2*(1.96 * preds$se.fit)  # ~ lower 99% CI to be conservative 
 
anomaly_fusion_reads[, p_predicted:= AnomalyGLM$family$linkinv(lwr)] # add lower 95% CI prediction
anomaly_fusion_reads[, p_predicted_2:= AnomalyGLM$family$linkinv(lwr_2)] # add lower 95% CI prediction

linear_and_anomaly_fusions=rbind(anomaly_fusion_reads[,list(qual,lenAdjScore,qualR2,lenAdjScoreR2,junction,is.pos,overlap,is.anomaly,p_predicted,p_predicted_2)], fusion_reads[,list(qual,lenAdjScore,qualR2,lenAdjScoreR2,junction,is.pos,overlap,is.anomaly,p_predicted,p_predicted_2)])
linear_and_anomaly_fusions$overlap=as.numeric(as.vector(linear_and_anomaly_fusions$overlap))
linearWithAnomalyFusionPredictions = predictNewClassP(linear_and_anomaly_fusions, null)

# start prediction on INDELS mapping
print ("done with fusion anomalies starting indels")

if (use.indels==1){
FJ_indel_reads$overlap=as.numeric(as.vector(FJ_indel_reads$overlap))
preds = predict(IndelGLM, newdata=FJ_indel_reads, type = "link", se.fit = TRUE) 
lwr = preds$fit - (1.96 * preds$se.fit)  # ~ lower 95% CI to be conservative 
lwr_2 = preds$fit - 2*(1.96 * preds$se.fit)  # ~ lower 99% CI to be conservative 

FJ_indel_reads[, p_predicted:= IndelGLM$family$linkinv(lwr)] # add lower 95% CI prediction
FJ_indel_reads[, p_predicted_2:= IndelGLM$family$linkinv(lwr_2)] # add lower 95% CI prediction
}
}


if (use.indels==0 & no.indel.classfile==0){
## if not using indels, we assign indel reads according to the linear decoy model
preds = predict(linearDecoyGLM, newdata=FJ_indel_reads, type = "link", se.fit = TRUE) 
lwr = preds$fit - (1.96 * preds$se.fit)  # ~ lower 95% CI to be conservative 
lwr_2 = preds$fit - 2*(1.96 * preds$se.fit)  # ~ lower 99% CI to be conservative 
FJ_indel_reads[, p_predicted:= linearDecoyGLM$family$linkinv(lwr)] # add lower 95% CI prediction
FJ_indel_reads[, p_predicted_2:= linearDecoyGLM$family$linkinv(lwr_2)] # add lower 95% CI prediction

linear_and_anomaly_and_indel_fusions=rbind(FJ_indel_reads[,list(qual,lenAdjScore,qualR2,lenAdjScoreR2,junction,is.pos,overlap,is.anomaly,p_predicted,p_predicted_2)],anomaly_fusion_reads[,list(qual,lenAdjScore,qualR2,lenAdjScoreR2,junction,is.pos,overlap,is.anomaly,p_predicted,p_predicted_2)], fusion_reads[,list(qual,lenAdjScore,qualR2,lenAdjScoreR2,junction,is.pos,overlap,is.anomaly,p_predicted,p_predicted_2)])
}
if (no.indel.classfile==1){
#don't bind FJ_indel_reads bc they don't exist
linear_and_anomaly_and_indel_fusions=rbind(anomaly_fusion_reads[,list(qual,lenAdjScore,qualR2,lenAdjScoreR2,junction,is.pos,overlap,is.anomaly,p_predicted,p_predicted_2)], fusion_reads[,list(qual,lenAdjScore,qualR2,lenAdjScoreR2,junction,is.pos,overlap,is.anomaly,p_predicted,p_predicted_2)])
}

linear_and_anomaly_and_indel_fusions$overlap=as.numeric(as.vector(linear_and_anomaly_and_indel_fusions$overlap))
linearWithAnomalyAndIndelFusionPredictions = predictNewClassP(linear_and_anomaly_and_indel_fusions, null)

consolidated_fusion=merge(fusionJunctionPredictions,linearWithAnomalyFusionPredictions , all=TRUE,by="junction")
consolidated_fusion[,p_diff:=(p_predicted.x-p_predicted.y)] ## p_predicted.x should be less than p_predicted.y always so p_diff should be neegative
consolidated_fusion=data.table(unique(consolidated_fusion))

####### now add indels
#verbose adds extra output and extra info to indel table
verbose=0

print ("now adding indels and changing names to reflect accurate terminology")
consolidated_fusion_windel=merge(linearWithAnomalyFusionPredictions ,linearWithAnomalyAndIndelFusionPredictions, all=TRUE,by="junction")
consolidated_fusion_windel[,p_diff_indel:=(p_predicted.y-p_predicted.x)] ## p_predicted.x should be less than p_predicted.y always so p_diff should be neegative
consolidated_fusion_windel=data.table(unique(consolidated_fusion_windel))

setnames(consolidated_fusion_windel,"p_predicted.x", "productPhat.x")
setnames(consolidated_fusion_windel,"p_predicted.y", "productPhat.y")
setnames(consolidated_fusion_windel,"p_value.x", "junction_cdf.x")
setnames(consolidated_fusion_windel,"p_value.y", "junction_cdf.y")


setnames(consolidated_fusion_windel,"p_predicted_2.x", "productPhat_lower.x")
setnames(consolidated_fusion_windel,"p_predicted_2.y", "productPhat_lower.y")
setnames(consolidated_fusion_windel,"p_value_2.x", "junction_cdf_lower.x")
setnames(consolidated_fusion_windel,"p_value_2.y", "junction_cdf_lower.y")

setnames(consolidated_fusion_windel,"p_diff_indel", "junction_cdf_windel_diff")
if (verbose==0){
consolidated_fusion_windel[,q_1.y:=NULL]
consolidated_fusion_windel[,q_2.y:=NULL]
consolidated_fusion_windel[,q_1.x:=NULL]
consolidated_fusion_windel[,q_2.x:=NULL]
consolidated_fusion_windel[,logsum.x:=NULL]
consolidated_fusion_windel[,logsum.y:=NULL]
}

write.table(unique(consolidated_fusion_windel)[order(-junction_cdf.y),], fusionwanomaly_and_indel_juncp_out, row.names=FALSE, quote=FALSE, sep="\t")
cir.verbose=1
####################################################################################
####################################################################################
if (cir.verbose==1){
## write circle prediction

setnames(circularJunctionPredictions,"p_predicted", "productPhat.x")
setnames(circularJunctionPredictions,"p_predicted_2", "productPhat_lower.x")
setnames(circularJunctionPredictions,"p_value_2", "junction_cdf_lower.x")
setnames(circularJunctionPredictions,"p_value", "junction_cdf")

write.table(unique(linearJunctionPredictions), linear_juncp_out, row.names=FALSE, quote=FALSE, sep="\t")
write.table(unique(circularJunctionPredictions), circ_juncp_out, row.names=FALSE, quote=FALSE, sep="\t")

## with anomalies, p value should be smaller so p_diff should always be negative...
consolidated_linear=merge(linearWithAnomalyJunctionPredictions,linearJunctionPredictions,by="junction",all=TRUE) 
consolidated_linear[,p_diff:=(p_predicted.x-p_predicted.y)] ## p_predicted.x should be less than p_predicted.y always so p_diff should be neegative
consolidated_linear=data.table(unique(consolidated_linear))
write.table(unique(consolidated_linear), linearwanomaly_juncp_out, row.names=FALSE, quote=FALSE, sep="\t")

## 
my.null.quantiles=quantile(linear_reads$p_predicted,probs=c(0:10)/10,na.rm=TRUE)
## refer fusions to these quantiles; 'falsely called' vs. true will be fraction of linears (conservative estimate) ; error at this quantile can be evaluated. 

write.table(unique(fusionJunctionPredictions[order(-p_predicted),]),fusion_juncp_out, row.names=FALSE, quote=FALSE, sep="\t")
write.table(unique(consolidated_fusion)[order(-p_predicted.y),], fusionwanomaly_juncp_out, row.names=FALSE, quote=FALSE, sep="\t")

}
