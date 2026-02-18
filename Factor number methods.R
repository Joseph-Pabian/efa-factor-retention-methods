library(foreign)
library(psych)
library(nFactors) 
libqual<-read.spss("/Users/josephpabian/Downloads/libqual-1.sav", to.data.frame=TRUE,use.value.labels=FALSE)
attach(libqual)
libqual.pa.promax.factor2<-fa(libqual[,3:14], nfactors=2, rotate="promax", SMC=TRUE, fm="pa", residuals=TRUE)
libqual.pa.promax.factor2

#SCREE PLOT#
scree(libqual[,3:14])

#OBJECTIVE SCREE PLOT MEASURES#
libqual.eigenvalues<-eigen(cor(libqual[,3:14], use="complete"))
plotnScree(nScree(x=libqual.eigenvalues$values, model="factors"))

#ZOSKI AND JURS REGRESSION#
nMreg(libqual[,3:14], model="factors")

#CHI SQUARE TEST FROM MLE 2 factors
libqual.ml.promax.factor2<-fa(libqual[,3:14], nfactors=2, rotate="promax", SMC=TRUE, fm="ml")
libqual.ml.promax.factor2
#CHI SQUARE TEST FROM MLE 3 factors
libqual.ml.promax.factor3<-fa(libqual[,3:14], nfactors=3, rotate="promax", SMC=TRUE, fm="ml")
libqual.ml.promax.factor3
#CHI SQUARE TEST FROM MLE 4 factors
libqual.ml.promax.factor4<-fa(libqual[,3:14], nfactors=4, rotate="promax", SMC=TRUE, fm="ml")
libqual.ml.promax.factor4
#PARALLEL ANALYSIS#
libqual.parallel<-fa.parallel(libqual[,3:14], fm="pa",n.iter=1000,error.bars=FALSE,SMC=TRUE)

#REVISED PARALLEL ANALYSIS#
#FIRST SUBMIT THE COMMANDS CONTAINED IN THE FILE REVISED PARALLEL ANALYSIS FUNCTIONS#
#THEN SUBMIT THE FOLLOWING COMMANDS#
libqual.cor<-cor(libqual[,3:14])
revised.pa(data=libqual.cor, N=1000)

#COMPARATIVE DATA METHOD#
#FIRST SUBMIT THE COMMANDS CONTAINED IN THE FILE COMPARISON DATA METHOD FUNCTIONS#
#THEN SUBMIT THE FOLLOWING COMMANDS#
EFA.Comp.Data(libqual[,3:14],F.Max=4,N.Pop=10000,N.Samples=500,Alpha=0.30,Graph=F,Spearman=F)

#MAP#
libqual.vss<-vss(libqual[,3:14], n=8, rotate="promax", diagonal=FALSE, fm="pa", plot=TRUE, SMC=FALSE)
libqual.vss$map

#VSS#
plot(libqual.vss)
libqual.vss

#RESIDUAL CORRELATION MATRIX#
libqual.pa.promax.factor2.resid<-factor.residuals(cor(libqual[,3:14]),libqual.pa.promax.factor2)
libqual.pa.promax.factor2.resid
