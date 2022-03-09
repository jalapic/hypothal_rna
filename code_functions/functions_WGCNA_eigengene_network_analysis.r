# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/EigengeneNetwork/

#--------------------------------------------------------------------------------------------
# Set global parameters.

# This parameter controls whether my implementation of pmax should use a call to an external function
# (recommended if the requisite library is available). Otherwise an R-only implementation will be used,
# which is significantly slower (but available and stable on all R platforms).

UseCpmax = FALSE;

#--------------------------------------------------------------------------------------------

# Load the requisite libraries

WorkingDirectory = getwd();

#if (!library(MASS, logical.return=TRUE)) { # standard, no need to install
#For some reason, MASS does not seem to be installed on Titan, so we'll try to load it 
# from my own library. If that fails as well, stop.
#  if (!library(MASS, logical.return=TRUE, lib.loc="M:/Work/RLibrary")) stop()
#}   

library(MASS);
library(class)	# standard, no need to install
library(cluster)	
#library(sma)	# install it for the function plot.mat 
library(impute)# install it for imputing missing value
library(Hmisc)	# install it for the C-index calculations
library(survival)
#library(fields);

#oldwd = getwd();


#setwd(oldwd);


#####################################################################################################
# A) Assessing scale free topology and choosing the parameters of the adjacency function
#    using the scale free topology criterion (Zhang and Horvath 05)
########################################################################################################


# ===================================================
#For hard thresholding, we use the signum (step) function
if(exists("signum1") ) rm(signum1); 
signum1=function(corhelp,tau1)  {
  adjmat1= as.matrix(abs(corhelp)>=tau1)
  dimnames(adjmat1) <- dimnames(corhelp)
  diag(adjmat1) <- 0
  adjmat1}

# ===================================================
# For soft thresholding, one can use the sigmoid function 
# But we have focused on the power adjacency function in the tutorial...
if (exists("sigmoid1") ) rm(sigmoid1); sigmoid1=function(ss,mu1=0.8,alpha1=20) {
  1/(1+exp(-alpha1*(ss-mu1)))}





#This function is useful for speeding up the connectivity calculation.
#The idea is to partition the adjacency matrix into consecutive baches of a #given size.
#In principle, the larger the batch size the faster is the calculation. But #smaller batchsizes require #less memory...
# Input: gene expression data set where *rows* correspond to microarray samples #and columns correspond to genes. 
# If fewer than MinimumNoSamples contain gene expression information for a given
# gene, then its connectivity is set to missing. 
if(exists("SoftConnectivity")) rm(SoftConnectivity);
SoftConnectivity=function(datE, power=6, batchsize=1500, MinimumNoSamples=10) {
  no.genes=dim(datE)[[2]]
  no.samples=dim(datE)[[1]]
  if (no.genes<no.samples | no.genes<10 | no.samples<5 ) {stop("Error: Something seems to be wrong. Make sure that the input data frame has genes as rows and array samples as columns. Alternatively, there could be fewer than 10 genes or fewer than 5 samples. ") } else {
    sum1=function(x) sum(x,na.rm=T)
    k=rep(NA,no.genes)
    no.batches=as.integer(no.genes/ batchsize)
    if (no.batches>0) {
      for (i in 1:no.batches) {
        print(paste("batch number = ", i))
        index1=c(1:batchsize)+(i-1)* batchsize
        ad1=abs(cor(datE[,index1], datE,use="p"))^power
        ad1[is.na(ad1)]=0
        k[index1]=apply(ad1,1,sum1)
        # If fewer than MinimumNoSamples contain gene expression information for a given
        # gene, then we set its connectivity to 0.
        NoSamplesAvailable=apply(!is.na(datE[,index1]),2,sum)
        k[index1][NoSamplesAvailable< MinimumNoSamples]=NA
      } # end of for (i in 1:no.batches
    } # end of if (no.batches>0)...
    if (no.genes-no.batches*batchsize>0 ) {
      restindex=c((no.batches*batchsize+1):no.genes)
      ad1=abs(cor(datE[,restindex], datE,use="p"))^power
      ad1[is.na(ad1)]=0
      k[restindex]=apply(ad1,1,sum1)
      NoSamplesAvailable=apply(!is.na(datE[,restindex]),2,sum)
      k[restindex][NoSamplesAvailable< MinimumNoSamples]=NA
    } # end of if
  } # end of else statement
  k
} # end of function




# ===================================================
# The function PickHardThreshold can help one to estimate the cut-off value 
# when using the signum (step) function.
# The first column lists the threshold ("cut"),
# the second column lists the corresponding p-value based on the Fisher Transform 
# of the correlation. 
# The third column reports the resulting scale free topology fitting index R^2.
# The fourth column reports the slope of the fitting line, it shoud be negative for 
# biologically meaningul networks.
# The fifth column reports the fitting index for the truncated exponential model. 
# Usually we ignore it.
# The remaining columns list the mean, median and maximum resulting connectivity.
# To pick a hard threshold (cut) with the scale free topology criterion:
# aim for high scale free R^2 (column 3), high connectivity (col 6) and negative slope 
# (around -1, col 4).
# The output is a list with 2 components. The first component lists a sugggested cut-off
# while the second component contains the whole table.
# The removeFirst option removes the first point (k=0, P(k=0)) from the regression fit.
# no.breaks specifies how many intervals used to estimate the frequency p(k) i.e. the no. of points in the 
# scale free topology plot.
if (exists("PickHardThreshold")) rm(PickHardThreshold);
PickHardThreshold=function(datExpr1,RsquaredCut=0.85, cutvector=seq(0.1,0.9,by=0.05) ,removeFirst=FALSE,no.breaks=10) {
  no.genes   <- dim(datExpr1)[[2]]
  no.genes <- dim(datExpr1)[[2]]
  no.samples= dim(datExpr1)[[1]]
  colname1=c("Cut","p-value", "scale law R^2", "slope="  ,"truncated R^2","mean(k)","median(k)","max(k)")
  datout=data.frame(matrix(666,nrow=length(cutvector),ncol=length(colname1) ))
  names(datout)=colname1
  datout[,1]=cutvector
  for (i in c(1:length(cutvector) ) ){
    cut1=cutvector[i]
    datout[i,2]=2*(1-pt(sqrt(no.samples-1)*cut1/sqrt(1-cut1^2),no.samples-1))}
  if(exists("fun1")) rm(fun1)
  fun1=function(x) {
    corx=abs(cor(x,datExpr1,use="p"))
    out1=rep(NA, length(cutvector) )
    for (j in c(1:length(cutvector))) {out1[j]=sum(corx>cutvector[j])}
    out1
  } # end of fun1
  datk=t(apply(datExpr1,2,fun1))
  for (i in c(1:length(cutvector) ) ){
    nolinkshelp <- datk[,i]-1
    cut2=cut(nolinkshelp,no.breaks)
    binned.k=tapply(nolinkshelp,cut2,mean)
    freq1=as.vector(tapply(nolinkshelp,cut2,length)/length(nolinkshelp))
    # The following code corrects for missing values etc
    breaks1=seq(from=min(nolinkshelp),to=max(nolinkshelp),length=no.breaks+1)
    hist1=hist(nolinkshelp,breaks=breaks1,equidist=F,plot=FALSE,right=TRUE)
    binned.k2=hist1$mids
    binned.k=ifelse(is.na(binned.k),binned.k2,binned.k)
    binned.k=ifelse(binned.k==0,binned.k2,binned.k)
    freq1=ifelse(is.na(freq1),0,freq1)
    xx= as.vector(log10(binned.k))
    if(removeFirst) {freq1=freq1[-1]; xx=xx[-1]}
    plot(xx,log10(freq1+.000000001),xlab="log10(k)",ylab="log10(p(k))" )
    lm1= lm(as.numeric(log10(freq1+.000000001))~ xx )
    lm2=lm(as.numeric(log10(freq1+.000000001))~ xx+I(10^xx) )
    datout[i,3]=summary(lm1)$adj.r.squared 
    datout[i,4]=summary(lm1)$coefficients[2,1]  
    datout[i,5]=summary(lm2)$adj.r.squared
    datout[i,6]=mean(nolinkshelp)
    datout[i,7]= median(nolinkshelp)
    datout[i,8]= max(nolinkshelp) 
  } 
  datout=signif(datout,3) 
  print(data.frame(datout));
  # the cut-off is chosen as smallest cut with R^2>RsquaredCut 
  ind1=datout[,3]>RsquaredCut
  indcut=NA
  indcut=ifelse(sum(ind1)>0,min(c(1:length(ind1))[ind1]),indcut)
  # this estimates the cut-off value that should be used. 
  # Don't trust it. You need to consider slope and mean connectivity as well!
  cut.estimate=cutvector[indcut][[1]]
  list(cut.estimate, data.frame(datout));
} # end of function











# ===========================================================
# The function PickSoftThreshold allows one to estimate the power parameter when using
# a soft thresholding approach with the use of the power function AF(s)=s^Power
# The function PickSoftThreshold allows one to estimate the power parameter when using
# a soft thresholding approach with the use of the power function AF(s)=s^Power
# The removeFirst option removes the first point (k=1, P(k=1)) from the regression fit.
if (exists("PickSoftThreshold")) rm(PickSoftThreshold);
PickSoftThreshold=function(datExpr1,RsquaredCut=0.85, powervector=c(seq(1,10,by=1),seq(12,20,by=2)),
                           removeFirst=FALSE,no.breaks=10) {
  no.genes <- dim(datExpr1)[[2]]
  no.samples= dim(datExpr1)[[1]]
  colname1=c("Power", "scale law R^2" ,"slope", "truncated R^2","mean(k)","median(k)","max(k)")
  datout=data.frame(matrix(666,nrow=length(powervector),ncol=length(colname1) ))
  names(datout)=colname1
  datout[,1]=powervector
  if(exists("fun1")) rm(fun1)
  fun1=function(x) {
    corx=abs(cor(x,datExpr1,use="p"))
    out1=rep(NA, length(powervector) )
    for (j in c(1:length(powervector))) {out1[j]=sum(corx^powervector[j])}
    out1
  } # end of fun1
  datk=t(apply(datExpr1,2,fun1))
  for (i in c(1:length(powervector) ) ){
    nolinkshelp <- datk[,i]-1
    cut2=cut(nolinkshelp,no.breaks)
    binned.k=tapply(nolinkshelp,cut2,mean)
    freq1=as.vector(tapply(nolinkshelp,cut2,length)/length(nolinkshelp))
    # The following code corrects for missing values etc
    breaks1=seq(from=min(nolinkshelp),to=max(nolinkshelp),length=no.breaks+1)
    hist1=hist(nolinkshelp,breaks=breaks1,equidist=F,plot=FALSE,right=TRUE)
    binned.k2=hist1$mids
    binned.k=ifelse(is.na(binned.k),binned.k2,binned.k)
    binned.k=ifelse(binned.k==0,binned.k2,binned.k)
    freq1=ifelse(is.na(freq1),0,freq1)
    
    xx= as.vector(log10(binned.k))
    if(removeFirst) {freq1=freq1[-1]; xx=xx[-1]}
    plot(xx,log10(freq1+.000000001),xlab="log10(k)",ylab="log10(p(k))" )
    lm1= lm(as.numeric(log10(freq1+.000000001))~ xx )
    lm2=lm(as.numeric(log10(freq1+.000000001))~ xx+I(10^xx) )
    datout[i,2]=summary(lm1)$adj.r.squared 
    datout[i,3]=summary(lm1)$coefficients[2,1]  
    datout[i,4]=summary(lm2)$adj.r.squared
    datout[i,5]=mean(nolinkshelp)
    datout[i,6]= median(nolinkshelp)
    datout[i,7]= max(nolinkshelp) 
  } 
  datout=signif(datout,3) 
  print(data.frame(datout));
  # the cut-off is chosen as smallest cut with R^2>RsquaredCut 
  ind1=datout[,2]>RsquaredCut
  indcut=NA
  indcut=ifelse(sum(ind1)>0,min(c(1:length(ind1))[ind1]),indcut)
  # this estimates the power value that should be used. 
  # Don't trust it. You need to consider slope and mean connectivity as well!
  power.estimate=powervector[indcut][[1]]
  list(power.estimate, data.frame(datout));
}






# ===================================================
# The function ScaleFreePlot1 creates a plot for checking scale free topology
# when truncated1=T is specificed, it provides the R^2 measures for the following
# degree distributions: a) scale free topology, b) log-log R^2 and c) truncated exponential R^2

# The function ScaleFreePlot1 creates a plot for checking scale free topology
if(exists("ScaleFreePlot1")) rm(ScaleFreePlot1) ; 

ScaleFreePlot1=function(kk,no.breaks=10,AF1="" ,truncated1=FALSE, removeFirst=FALSE,cex.lab1=1){
  
  #bin data into no.breaks bins: first create the factor cut1 and code the values
  cut1=cut(kk,no.breaks)
  #now calculate the mean of each bin
  binned.k=tapply(kk,cut1,mean)
  freq1=tapply(kk,cut1,length)/length(kk)
  # The following code corrects for missing values etc
  breaks1=seq(from=min(kk, na.rm=T),to=max(kk, na.rm=T),length=no.breaks+1)
  hist1=hist(kk,breaks=breaks1,equidist=F,plot=FALSE,right=TRUE)
  binned.k2=hist1$mids
  binned.k=ifelse(is.na(binned.k),binned.k2,binned.k)
  binned.k=ifelse(binned.k==0,binned.k2,binned.k)
  freq1=ifelse(is.na(freq1),0,freq1)
  plot(log10(binned.k),log10(freq1+.000000001),xlab=paste(AF1,"log10(k)"),ylab="log10(p(k))",cex.lab=cex.lab1 )
  xx= as.vector(log10(binned.k))
  if(removeFirst) {freq1=freq1[-1]; xx=xx[-1]}
  lm1=lm(as.numeric(log10(freq1+.000000001))~ xx )
  lines(xx,predict(lm1),col=1)
  OUTPUT=data.frame(ScaleFreeRsquared=round(summary(lm1)$adj.r.squared,2),Slope=round(lm1$coefficients[[2]],2))
  if (truncated1==TRUE) { 
    lm2=lm(as.numeric(log10(freq1+.000000001))~ xx+I(10^xx) );
    OUTPUT=data.frame(ScaleFreeRsquared=round(summary(lm1)$adj.r.squared,2),Slope=round(lm1$coefficients[[2]],2),
                      TruncatedRsquared=round(summary(lm2)$adj.r.squared,2))
    print("the red line corresponds to the truncated exponential fit")
    lines(xx,predict(lm2),col=2);
    title(paste( 
      "scale free R^2=",as.character(round(summary(lm1)$adj.r.squared,2)),
      ", slope=", round(lm1$coefficients[[2]],2),
      ", trunc.R^2=",as.character(round(summary(lm2)$adj.r.squared,2))))} else { 
        title(paste("R^2=",as.character(round(summary(lm1)$adj.r.squared,2)) , 
                    " sl=", round(lm1$coefficients[[2]],2)), cex=0.4)
      }
  OUTPUT
} # end of function 









#################################################################################################################
################################################################################################################################
# B) Computing the topological overlap matrix 
#################################################################################################################
#################################################################################################################



# ===================================================
#The function TOMdist1 computes a dissimilarity 
# based on the topological overlap matrix (Ravasz et al)
# Input: an Adjacency matrix with entries in [0,1]
if(exists("TOMdist1")) rm(TOMdist1);

TOMdist1=function(adjmat1, maxADJ=FALSE) {
  diag(adjmat1)=0;
  adjmat1[is.na(adjmat1)]=0;
  maxh1=max(as.dist(adjmat1) ); minh1=min(as.dist(adjmat1) ); 
  if (maxh1>1 | minh1 < 0 ) {print(paste("ERROR: the adjacency matrix contains entries that are larger than 1 or smaller than 0!!!, max=",maxh1,", min=",minh1)) } else { 
    if (  max(c(as.dist(abs(adjmat1-t(adjmat1)))))>0   ) {print("ERROR: non-symmetric adjacency matrix!!!") } else { 
      kk=apply(adjmat1,2,sum)
      maxADJconst=1
      if (maxADJ==TRUE) maxADJconst=max(c(as.dist(adjmat1 ))) 
      Dhelp1=matrix(kk,ncol=length(kk),nrow=length(kk))
      denomTOM= pmin(as.dist(Dhelp1),as.dist(t(Dhelp1)))   +as.dist(maxADJconst-adjmat1); 
      gc();gc();
      numTOM=as.dist(adjmat1 %*% adjmat1 +adjmat1);
      #TOMmatrix=numTOM/denomTOM
      # this turns the TOM matrix into a dissimilarity 
      out1=1-as.matrix(numTOM/denomTOM) 
      diag(out1)=1
      out1
    }}
}

#---------------------------------------------------------------------------
# This is a somewhat modified TOMdist1.

SignedTOMdist = function(adjmat1, maxADJ=FALSE)
{
  diag(adjmat1)=0;
  adjmat1[is.na(adjmat1)]=0;
  collect_garbage();
  kk=apply(adjmat1,2,sum)
  collect_garbage();
  maxADJconst=1
  if (maxADJ==TRUE) maxADJconst=max(c(as.dist(adjmat1 ))) 
  collect_garbage();
  Dhelp1 = matrix(kk,ncol=length(kk),nrow=length(kk))
  collect_garbage();
  denomTOM = pmin(as.dist(Dhelp1),as.dist(t(Dhelp1))) + as.dist(maxADJconst-adjmat1); 
  rm(Dhelp1);
  collect_garbage();
  gc(); gc();
  numTOM=as.dist(adjmat1 %*% adjmat1 +adjmat1);
  collect_garbage();
  #TOMmatrix=numTOM/denomTOM
  # this turns the TOM matrix into a dissimilarity 
  out1=1-as.matrix(numTOM/denomTOM) 
  rm(numTOM); rm(denomTOM);
  collect_garbage();
  diag(out1)=1
  out1
}



# ===================================================
# This function computes a TOMk dissimilarity
# which generalizes the topological overlap matrix (Ravasz et al)
# Input: an Adjacency matrix with entries in [0,1]
# WARNING:  ONLY FOR UNWEIGHTED NETWORKS, i.e. the adjacency matrix contains binary entries...
# This function is explained in Yip and Horvath (2005)
# https://horvath.genetics.ucla.edu/html/TOM/
if(exists("TOMkdist1")) rm(TOMkdist1);
TOMkdist1 = function(adjmat1,k=1){
  maxh1=max(as.dist(adjmat1) ); minh1=min(as.dist(adjmat1) );
  if (k!=round(abs(k))) {
    stop("k must be a positive integer!!!", call.=TRUE);}
  if (maxh1>1 | minh1 < 0 ){
    print(paste("ERROR: entries of the adjacency matrix must be between inclusively 0 and 1!!!, max=",maxh1,", min=",minh1))}
  else {
    if (  max(c(as.dist(abs(adjmat1-t(adjmat1)))))>0   ) {print("ERROR: non-symmetric adjacency matrix!!!") } else { 
      
      B <- adjmat1;
      if (k>=2) {
        for (i in 2:k) {
          diag(B) <- diag(B) + 1;
          B = B %*% adjmat1;}}   # this gives the number of paths with length at most k connecting a pair
      B <- (B>0);   # this gives the k-step reachability from a node to another
      diag(B) <- 0;   # exclude each node being its own neighbor
      B <- B %*% B   # this gives the number of common k-step-neighbor that a pair of nodes share
      
      Nk <- diag(B);
      B <- B +adjmat1;   # numerator
      diag(B) <- 1;
      denomTOM=outer(Nk,Nk,FUN="pmin")+1-adjmat1;
      diag(denomTOM) <- 1;
      1 - B/denomTOM   # this turns the TOM matrix into a dissimilarity
    }}
}


# IGNORE THIS function...
# The function TOMdistROW computes the TOM distance of a gene (node)
# with that of all other genes in the network.
# WhichRow is an integer that specifies which row of the adjacency matrix
# corresponds to the gene of interest.
# Output=vector of TOM distances.
if (exists("TOMdistROW") ) rm(TOMdistROW) 
TOMdistROW=function(WhichRow=1, adjmat1, maxADJ=FALSE) {
  diag(adjmat1)=0;
  maxh1=max(as.dist(adjmat1) ); minh1=min(as.dist(adjmat1) ); 
  if (maxh1>1 | minh1 < 0 ) {print(paste("ERROR: the adjacency matrix contains entries that are larger than 1 or smaller than 0!!!, max=",maxh1,", min=",minh1)) } else { 
    kk=apply(adjmat1,2,sum)
    numTOM=adjmat1[WhichRow,] %*% adjmat1 +adjmat1[WhichRow,]; 
    numTOM[WhichRow]=1
    maxADJconst=1
    if (maxADJ==TRUE) maxADJconst=max(c(as.dist(adjmat1 ))) 
    denomTOM=pmin(kk[WhichRow],kk)+maxADJconst-adjmat1[WhichRow,]; denomTOM[WhichRow]=1
    #TOMmatrix=numTOM/denomTOM
    # this turns the TOM matrix into a dissimilarity 
    1-numTOM/denomTOM 
  }
}


#####################################################################################################
################################################################################################################################
# C) Defining gene modules using clustering procedures
#####################################################################################################
################################################################################################################################

# ===================================================
#The function modulecolor2 function assigns colors to the observations 
# in the branches of a dendrogram
# we use it to define modules....
if (exists("modulecolor2")) rm(modulecolor2);
modulecolor2=function(hier1, h1=0.9,minsize1=50) {
  # here we define modules by using a height cut-off for the branches
  labelpred= cutree(hier1,h=h1)
  sort1=-sort(-table(labelpred))
  modulename= as.numeric(names(sort1))
  modulebranch= sort1>minsize1
  no.modules=sum(modulebranch)
  # now we assume that there are fewer than a certain number of colors
  #colorcode=c("turquoise","blue","brown","yellow","green","red","black","purple","orange","pink",
  #"greenyellow","lightcyan","salmon","midnightblue","lightyellow")
  colorcode=c("turquoise","blue","brown","yellow","green","red","black","pink","magenta",
              "purple","greenyellow","tan","salmon", "midnightblue", "lightcyan","grey60",
              "lightgreen", "lightyellow", "royalblue", "darkred", "darkgreen", "darkturquoise",
              "darkgrey", "orange", "darkorange", "white" )
  
  # "grey" means not in any module;
  colorhelp=rep("grey",length(labelpred))
  if ( no.modules==0 | no.modules >length(colorcode)){ print(paste("The number of modules is problematic. Number of modules = ", as.character(no.modules)))} else { for (i in c(1:no.modules)) {colorhelp=ifelse(labelpred==modulename[i],colorcode[i],colorhelp)};
    colorhelp=factor(colorhelp,levels=c(colorcode[1:no.modules],"grey"))
  }
  factor(colorhelp, levels=unique(colorhelp[hier1$order] ))
}


#---------------------------------------------------------------------------------
#
# ModuleNumber
#
#---------------------------------------------------------------------------------
# Similar to modulecolor2 above, but returns numbers instead of colors, which is oftentimes more useful.
# 0 means unassigned.
# Return value is a simple vector, not a factor.
# Caution: the module numbers are neither sorted nor sequential, the only guarranteed fact is that grey
# probes are labeled by 0 and all probes belonging to the same module have the same number.
# If size-sorted sequential labels are required, "normalize" the result by calling NormalizeLabels
# (below).

ModuleNumber = function(HierTree, CutHeight = 0.9, MinSize = 50)
{
  Branches = cutree(HierTree, h = CutHeight);
  NOnBranches = table(Branches);
  #NOnBranches[i][[1]] somehow gives the number of elements on branch i.
  TrueBranch = NOnBranches >= MinSize;
  Branches[!TrueBranch[Branches]] = 0;
  
  #NewLabels = levels(factor(Branches));
  #for (lab in 1:length(NewLabels)) if (NewLabels[lab]!=0)
  #  Branches[Branches==NewLabels[lab]] = lab;
  
  Branches;
  
}


# The function hclustplot1 creates a barplot where the colors of the bars are sorted according to 
# a hierarchical clustering tree (hclust object)
#if (exists("hclustplot1")) rm(hclustplot1);
#hclustplot1=function(hier1,couleur,title1="Colors sorted by hierarchical clustering") 
#{
#if (length(hier1$order) != length(couleur) ) {print("ERROR: length of color vector not compatible with no. of objects in the hierarchical tree")};
#if (length(hier1$order) == length(couleur) ) {
#barplot(height=rep(1, length(couleur)), col= as.character(couleur[hier1$order]),border=F, main=title1,space=0, axes=F)}
#}

if (exists("hclustplot1")) rm(hclustplot1);
hclustplot1=function(hier1,Color1, Color2=NULL,title1="Colors sorted by hierarchical clustering") 
{
  options(stringsAsFactors=FALSE);
  if (length(hier1$order) != length(Color1) ) 
  { 
    stop("ERROR: length of color vector not compatible with no. of objects in the hierarchical tree");
  } else {
    if (is.null(Color2))
    {
      barplot(height=rep(1, length(Color1)), col= as.character(Color1[hier1$order]),
              border=F, main=title1,space=0, axes=F)
    } else if (length(Color1)==length(Color2)) {
      # height = matrix(0.5, nrow = 2, ncol = length(Color1));
      C1 = Color1[hier1$order]; C2 = Color2[hier1$order]
      step = 1/length(Color1);
      barplot(height=1, col = "white", border=F, main=title1,space=0, axes=F)
      for (i in 1:(length(Color1)))
      {
        lines(x=rep((i*step), times=2), y=c(0,0.5),  col = as.character(C1[i]));
        lines(x=rep((i*step), times=2), y=c(0.5,1),  col = as.character(C2[i]));
      } 
    }
  }
}

if (exists("hclustplotn")) rm(hclustplotn);
hclustplotn=function(hier1, Color, RowLabels=NULL, cex.RowLabels = 0.9, ...) 
{
  options(stringsAsFactors=FALSE);
  if (length(hier1$order) != dim(Color)[[1]] ) 
  { 
    stop("ERROR: length of color vector not compatible with no. of objects in the hierarchical tree");
  } else {
    No.Sets = dim(Color)[[2]];
    C = Color[hier1$order, ]; 
    step = 1/dim(Color)[[1]];
    ystep = 1/No.Sets;
    barplot(height=1, col = "white", border=F,space=0, axes=F, ...)
    for (j in 1:No.Sets)
    {
      ind = (1:(dim(C)[1]));
      xl = (ind-1) * step; xr = ind * step; 
      yb = rep(ystep*(j-1), dim(C)[1]); yt = rep(ystep*j, dim(C)[1]);
      rect(xl, yb, xr, yt, col = as.character(C[,j]), border = as.character(C[,j]));
      if (is.null(RowLabels))
      {
        text(as.character(j), pos=2, x=0, y=ystep*(j-0.5), cex=cex.RowLabels, xpd = TRUE);
      } else {
        text(RowLabels[j], pos=2, x=0, y=ystep*(j-0.5), cex=cex.RowLabels, xpd = TRUE);
      }
    }
    for (j in 1:No.Sets) lines(x=c(0,1), y=c(ystep*j,ystep*j));
  }
}


# ===================================================
# The function TOMplotn creates a TOM plot
# Inputs:  distance measure, hierarchical (hclust) object, color label=couleur
if (exists("TOMplotn")) rm(TOMplotn);
TOMplot1=function(disttom,hier1, couleur,terrainColors=FALSE) {
  no.nodes=length(couleur)
  if (no.nodes != dim(disttom)[[1]] ) {print("ERROR: number of color labels does not equal number of nodes in disttom")} else {
    labeltree=as.character(couleur)
    labelrow  = labeltree
    labelrow[hier1$order[length(labeltree):1]]=labelrow[hier1$order]
    options(expressions = 10000)
    if (terrainColors) heatmap(as.matrix(disttom),Rowv=as.dendrogram(hier1),Colv= as.dendrogram(hier1), scale="none",revC=T, ColSideColors=as.character(labeltree),RowSideColors=as.character(labelrow),
                               labRow=F, labCol=F, col = terrain.colors(1000)) else heatmap(as.matrix(disttom),Rowv=as.dendrogram(hier1),Colv= as.dendrogram(hier1), scale="none",revC=T, ColSideColors=as.character(labeltree),RowSideColors=as.character(labelrow),
                                                                                            labRow=F, labCol=F)
  }
} #end of function


# ===================================================
# The function TOMplot2 creates a TOM plot where the top and left color bars can be different
# Inputs:  distance measure, hierarchical (hclust) object, color label=couleurTop, couleurLeft
if (exists("TOMplot2")) rm(TOMplot2);
TOMplot2=function(disttom,hier1, couleurTop, couleurLeft) {
  no.nodes=length(couleurTop)
  if (no.nodes != length(couleurLeft)) {stop("ERROR: number of top color labels does not equal number of left color labels")}
  if (no.nodes != dim(disttom)[[1]] ) {stop("ERROR: number of color labels does not equal number of nodes in disttom")} else {
    labeltree = as.character(couleurTop)
    labelrow  = as.character(couleurLeft)
    labelrow[hier1$order[length(labeltree):1]]=labelrow[hier1$order]
    options(expressions = 10000)
    heatmap(as.matrix(disttom),Rowv=as.dendrogram(hier1),Colv= as.dendrogram(hier1), scale="none", revC=T, ColSideColors=as.character(labeltree),RowSideColors=as.character(labelrow),
            labRow=F, labCol=F)
  }
} #end of function



# IGNORE THIS FUNCTION...
# The function "BestHeightCut" allows one to find the best height cut-off
# for a hierarchical clustering tree when external gene information is available
# It computes a Kruskal Wallis-test p-value for each height cut-off
# based on determining whether gene significance differs across branch membership.
if(exists("BestHeightCut")) rm(BestHeightCut);
BestHeightCut=function(hier1, GeneSignif, hcut=seq(0.1,.95,by=0.01) ) {
  pvalues=rep(NA, length(hcut))
  for (i in c(1:length(hcut))) {
    colorhelp=modulecolor2(hier1,hcut[i])
    if (length(unique(colorhelp))>1 ) {pvalues[i]=kruskal.test(GeneSignif, colorhelp)$p.value}
    data.frame(hcut,pvalues)
  }}




#####################################################################################################
################################################################################################################################
# D) Summing up modules using their first principal components (first eigengene)
#####################################################################################################
################################################################################################################################

# ===================================================
#The function ModulePrinComps1 finds the first principal component (eigengene) in each 
# module defined by the colors of the input vector "couleur" (Pardon my French).
# It also reports the variances explained by the first 5 principal components.
# And it yields a measure of module conformity for each gene,
# which is highly correlated to the within module connectivity.
# The theoretical underpinnings are described in Horvath, Dong, Yip (2005)
# https://horvath.genetics.ucla.edu/html/ModuleConformity/
# This requires the R library impute
# The output is a list with 3 components: 
# 1) a data frame of module eigengenes (MEs), 
# 2) a data frame that lists the percent variance explained by the first 5 MEs of a module
# 3) a data frame that lists the module conformity for each gene. 
# The be used as alternative connectivity measure....
if(exists("ModulePrinComps1")) rm(ModulePrinComps1);
ModulePrinComps1=function(datexpr, couleur, verbose = 1, print.level = 0, Impute = FALSE,
                          GetConformity = TRUE) {
  if (is.null(datexpr))
  {  
    print("ModulePrinComps1: Error: datexpr is NULL. ");
    stop();
  }
  if (is.null(couleur))
  {  
    print("ModulePrinComps1: Error: couleur is NULL. ");
    stop()
  }
  MaxVectors = 5;
  #print(paste("datexpr dimensions:", as.character(dim(datexpr))));
  spaces = PrintSpaces(print.level);
  modlevels=levels(factor(couleur))
  PrinComps=data.frame(matrix(666,nrow=dim(datexpr)[[1]],ncol= length(modlevels))) 
  varexplained= data.frame(matrix(666,nrow= 5,ncol= length(modlevels)))
  names(PrinComps)=paste("ME",modlevels,sep="")
  for(i in c(1:length(modlevels)) )
  {
    if (verbose>0) 
      print.flush(paste(spaces, "ModulePrinComps1 : Working on ME for module ", modlevels[i], sep = ""));
    modulename    = modlevels[i]
    restrict1= as.character(couleur)== modulename
    #print(paste("length of couleur:", length(couleur), "; length of restric1:", length(restrict1)));
    # in the following, rows are genes and columns are samples     
    datModule=t(as.matrix(datexpr[, restrict1]))
    if (Impute)
    {
      saved.seed = .Random.seed;
      datModule=impute.knn(as.matrix(datModule))$data
      .Random.seed = saved.seed;
    }
    datModule=t(scale(t(datModule)));
    n = dim(datModule)[1]; p = dim(datModule)[2];
    svd1=svd(datModule, nu = min(n, p, MaxVectors), nv = min(n, p, MaxVectors));
    mtitle=paste("MEs of ", modulename," module", sep="");
    varexplained[,i]= (svd1$d[1:5])^2/sum(svd1$d^2)
    # this is the first principal component
    pc1=svd1$v[,1]
    # signh1=sign(sum(cor(pc1,  t(datModule))))
    # if (signh1 != 0)  pc1=signh1* pc1
    PrinComps[,i]= pc1
  }
  ModuleConformity= rep(666,length=dim(datexpr)[[2]])
  if (GetConformity)
  {
    for(i in 1:(dim(datexpr)[[2]])) 
      ModuleConformity[i] = abs(cor(datexpr[,i], PrinComps[,match(couleur[i], modlevels)], 
                                    use="pairwise.complete.obs"))
  } else
  {
    ModuleConformity = NULL;
  }
  
  list(PrinComps=PrinComps, varexplained=varexplained, ModuleConformity=ModuleConformity)
}



#####################################################################################################
################################################################################################################################
# E) Relating a measure of gene significance to the modules 
#####################################################################################################
################################################################################################################################

# ===================================================
# The function ModuleEnrichment1 creates a bar plot that shows whether modules are enriched with
# significant genes.
# More specifically, it reports the mean gene significance for each module.
# The gene significance can be a binary variable or a quantitative variable. 
# It also plots the 95% confidence interval of the mean (CI=mean +/- 1.96* standard error).
# It also reports a Kruskal Wallis P-value.
if( exists("ModuleEnrichment1") ) rm(ModuleEnrichment1);
ModuleEnrichment1=function(genesignif1,couleur,title1="gene significance across modules",labely="Gene Significance",boxplot=F) {
  if (length(genesignif1) != length(couleur) ) print("Error: vectors don\'t have the same lengths") else {
    if (boxplot != TRUE) {
      mean1=function(x) mean(x,na.rm=T) 
      means1=as.vector(tapply(genesignif1,couleur,mean1));
      se1= as.vector(tapply(genesignif1,couleur,stderr1))
      #par(mfrow=c(1,1))
      barplot(means1,
              names.arg=names(table(couleur) ),col= names(table(couleur) )
              ,ylab=labely)
      err.bp(as.vector(means1), as.vector(1.96*se1), two.side=T)} else {
        boxplot(split(genesignif1,couleur),notch=T,varwidth=T, col= names(table(couleur) ),ylab=labely)}
    
    title(paste(title1,", p-value=", signif(kruskal.test(genesignif1,factor(couleur))$p.value,2)))
  }
} # end of function


# IGNORE THIS...
# ===================================================
#The function fisherPvector allows one to compute Fisher exact p-values
# Thus it allows one to carry out an EASE analysis
# Output: a table of Fisher¡¯s exact p-values
# Input: annotation1 is a vector of gene annotations
# Input: couleur (French for color) denotes the module color of each gene
# Only those gene functions (levels of annotation1) that occur a certain mininum number of times
# (parameter= minNumberAnnotation) in the data set will be considered.  
if (exists("fisherPvector" ) ) rm(fisherPvector);
fisherPvector=function(couleur,annotation1,minNumberAnnotation=50) {
  levelsannotation1=levels(annotation1)
  levelscouleur=levels(factor(couleur))
  no.couleur=length(levelscouleur)
  restannotation1=table(annotation1)>minNumberAnnotation
  no.annotation=sum( restannotation1)
  datoutP=data.frame(matrix(666,nrow=no.annotation,ncol=no.couleur) )
  #datoutProp=data.frame(matrix(666,nrow=no.annotation,ncol=2*no.couleur) )
  #names(datoutProp)=paste("Prop",paste( rep(levelscouleur ,rep(2, length(levelscouleur))) ) , c("Y","N")  ,sep=".")
  datoutProp=data.frame(matrix(666,nrow=no.annotation,ncol=no.couleur) )
  names(datoutProp)=paste("Perc",levelscouleur , sep=".")
  names(datoutP)=paste("p",levelscouleur,sep=".")
  restlevelsannotation1= levelsannotation1[restannotation1]
  row.names(datoutP)= restlevelsannotation1
  for (i in c(1:no.annotation) ) {
    for (j in c(1:no.couleur) ){
      tab1=table( annotation1 !=restlevelsannotation1[i], couleur !=levelscouleur[j])
      datoutP[i,j]=signif(fisher.test(tab1)$p.value,2) 
      #datoutProp[i,2*j-1]=signif(tab1[1,1]/sum(tab1[,1] ),2)
      #datoutProp[i,2*j]= signif(tab1[1,2]/sum(tab1[,2]) ,2)
    } 
    table2=table(annotation1 !=restlevelsannotation1[i], couleur)
    datoutProp[i,]= signif(table2[1,]/apply(table2,2,sum),2)
  }
  data.frame(datoutP,datoutProp)
} # end of function fisherPvector



#####################################################################################################
################################################################################################################################
# F) Carrying out a within module analysis (computing intramodular connectivity etc) 
#####################################################################################################
################################################################################################################################

# ===================================================
#The function DegreeInOut computes for each gene 
#a) the total number of connections, 
#b) the number of connections with genes within its module, 
#c) the number of connections with genes outside its module
# When scaledToOne=TRUE, the within module connectivities are scaled to 1, i.e. the max(K.Within)=1 for each module
if (exists("DegreeInOut")) rm(DegreeInOut); DegreeInOut =function(adj1, couleur,scaledToOne=FALSE) {
  no.nodes=length(couleur)
  couleurlevels=levels(factor(couleur))
  no.levels=length(couleurlevels)
  kWithin=rep(-666,no.nodes )
  diag(adj1)=0
  for (i in c(1:no.levels) ) {
    rest1=couleur==couleurlevels[i];
    if (sum(rest1) <3 ) { kWithin[rest1]=0 } else {
      kWithin[rest1]=apply(adj1[rest1,rest1],2,sum)
      if (scaledToOne) kWithin[rest1]=kWithin[rest1]/max(kWithin[rest1])}
  }
  kTotal= apply(adj1,2,sum) 
  kOut=kTotal-kWithin
  if (scaledToOne) kOut=NA
  kDiff=kWithin-kOut
  data.frame(kTotal,kWithin,kOut,kDiff)
}


# =======================================================================
# The function WithinModuleCindex1 relates the node measures (e.g. connectivities) 
# to "external" node significance information within each  module,
# i.e. it  carries out a by module analysis.
# Output: first column reports the spearman correlation p-value between the network variable and the 
# node significance. The next columns contain the Spearman correlations between the variables.
if (exists("WithinModuleAnalysis1")) rm(WithinModuleAnalysis1);
WithinModuleAnalysis1=function(datnetwork,nodesignif, couleur) 
{
  cortesthelp=function( x ) {
    len1=dim(x)[[2]]-1
    out1=rep(666, len1);
    for (i in c(1:len1) ) {out1[i]= signif( cor.test(x[,i+1], x[,1], method="s",use="p" )$p.value ,2) }
    data.frame( variable=names(x)[-1] , NS.CorPval=out1, NS.cor=t(signif(cor (x[,1], x[,-1],use="p",method="s"),2)), 
                signif(cor(x[,-1],use="p",method="s"),2) )
  } #end of function cortesthelp
  print("IGNORE  the warnings...");
  by( data.frame(nodesignif, datnetwork), couleur, cortesthelp);
} #end of function WithinModuleAnalysis


# =======================================================================
# The function WithinModuleCindex1 relates the node measures (e.g. connectivities) 
# to "external" node significance information within each  module, 
# i.e. it  carries out a by module analysis.
# BUT it focuses on the C-index also known as area under the ROC curve
# This measure is related to Kendall's Tau statistic and Somer's D, 
# see F. Harrel (Regression Modeling Strategies). Springer. 
# It requires the following library
library(Hmisc)
# Output: the first column reports the C-index and the second, p-value 
if (exists("WithinModuleCindex1")) rm(WithinModuleCindex1);
WithinModuleCindex1=function(datnetwork,nodesignif, couleur) {
  CindexFun=function( x ) {
    len1=dim(x)[[2]]-1
    outC=rep(666, len1);
    outP=rep(666, len1);
    for (i in c(1:len1) ) {rcor1=rcorr.cens(x[,i+1], x[,1])
    outC[i]=rcor1[[1]] 
    outP[i]=1- pnorm(abs(rcor1[[2]]/rcor1[[3]]))
    }
    data.frame( variable=names(x)[-1] , C.index=outC, p.value=outP)
  } #end of function CindexFun
  #print("IGNORE  the warnings...");
  by( data.frame(nodesignif, datnetwork),couleur,CindexFun);
} #end of function WithinModuleAnalysis


# The following function allows on to plot a gene (node) significance measure versus
# connectivity.
if(exists("plotConnectivityGeneSignif1") ) rm( plotConnectivityGeneSignif1);
plotConnectivityGeneSignif1=function(degree1,genesignif1,color1="black", 
                                     title1="Gene Significance vs Connectivity" , xlab1="Connectivity", ylab1="GeneSignificance") {
  lm1=lm(genesignif1~degree1 ,na.action="na.omit")
  plot(degree1, genesignif1, col=color1,ylab=ylab1,xlab=xlab1,main=paste(title1, ", cor=",  
                                                                         signif(cor( genesignif1,degree1, method="s",use="p" )   ,2) ))
  abline(lm1)
}
