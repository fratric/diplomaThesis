#' ---
#' title: "Diplomovka"
#' author: "Peter Fatrič"
#' output:
#'    html_document:
#'      toc: true
#'      toc_depth: 5
#'      number_sections: true
#' ---
#' 
# ------ Graphical models and machine learning ------
#' 
#' # grafické modely a strojové učenie - Markovovske Siete
#' 


#' 
#' # nacitanie dat
#' 

setwd("/Users/Peter/Desktop/diplomovka/programy/programDiplomovka/lepsieDataR")


library(readxl)
df<-NULL
mldata<-NULL
mldata<-read.csv("CASP.csv")
dim(mldata)

df<-mldata
df<-df[,c(2,3,4,5,6,7,8,9,10,1)]
#df<-apply(df,2,function(x) as.numeric(as.factor(x)))
dim(df)
#df<-subset(df, df$F7 != 0)


# take a random sample of size 50 from a dataset mydata
# sample without replacement
df <- df[sample(1:nrow(df), 10000, replace=FALSE),] 

Hmisc::hist.data.frame(df$RMSD,nclass = 50)

nlen<-dim(df)[1]
clen<-dim(df)[2]

#' rozdelenie na trenovaciu a testovaciu vzortku
set.seed(20)
indT <- sample(c(TRUE,FALSE), size = nlen, replace=T, prob=c(90,10))
indE <- (1:nrow(df))[!indT]
indT <- (1:nrow(df))[indT]
nT <- length(indT)
nE <- length(indE)
model <- list() 

train<-df[indT,]
test<-df[indE,]

R<-cor(df)
library("corpcor") 
Rp <- cor2pcor(R)
library("qgraph")
qgraph(round(R,5),edge.labels=TRUE,graph="concentration") 


#'
#'  # Linear model
#'

model$linearMod <- lm(RMSD ~ ., data=train)

tmpLin<-model$linearMod$coefficients %*% t(as.matrix(unname(cbind(1,test[-clen]))))
tmpLin<-as.vector(tmpLin)
plot(test$RMSD-tmpLin,type="l",main="linear regression - residues", ylab="error", xlab="")
varReg<-var(test$RMSD-tmpLin)




#'
#' # Priprava dat na spojity nahodny vektor
#'


model$Mcdf <- lapply(train, function(x) function(y) cdfgam(y, pelgam(samlmu(x)) ))
model$Mqua <- lapply(train, function(x) function(y) quagam(y, pelgam(samlmu(x)) ))

model$Mcdf$F3 <- function(y) cdfnor(y, pelnor(samlmu(train$F3)) )
model$Mqua$F3 <- function(y) quanor(y, pelnor(samlmu(train$F3)) )

model$Mcdf$F4 <- function(y) cdfln3(y, pelln3(samlmu(train$F4)) )
model$Mqua$F4 <- function(y) qualn3(y, pelln3(samlmu(train$F4)) )

model$Mcdf$F6 <- function(y) cdfln3(y, pelln3(samlmu(train$F6)) )
model$Mqua$F6 <- function(y) qualn3(y, pelln3(samlmu(train$F6)) )



library(mixR)
library(lmom)
sfit<-mixfit(sample(train$RMSD,size = nT), ncomp = 2, family = "normal")

model$Mden$RMSD <- function(x) (sfit$pi[1]*dnorm(x, mean = sfit$mu[1], sd = sfit$sd[1])) + 
  sfit$pi[2]*dnorm(x, mean = sfit$mu[2], sd = sfit$sd[2])
#model$Mcdf$RMSD <- Vectorize(function(x) integrate(model$Mden$RMSD, lower = -Inf, upper = x, abs.tol = 0L)$value,vectorize.args = "x")
model$Mqua$RMSD <- NULL
model$Mcdf$RMSD <- function(x) (sfit$pi[1]*cdfnor(x,c(sfit$mu[1],sfit$sd[1])) + sfit$pi[2]*cdfnor(x,c(sfit$mu[2],sfit$sd[2])))

ssq<-seq(0,25,0.1)
plot(x=ssq,y=sfit$pi[1]*dnorm(ssq, mean = sfit$mu[1], sd = sfit$sd[1]),type="l",xlab="",ylab="",lty=2, main="f_RMSD as sum of two normal densities")
lines(x=ssq,y=sfit$pi[2]*dnorm(ssq, mean = sfit$mu[2], sd = sfit$sd[2]),type="l",lty=3)
lines(x=ssq,y=model$Mden$RMSD(ssq),type="l")
abline(v=6.1)

s<-train$RMSD
plot(s, copula::pobs(s), main="Empiricka a paramericka CDF")
tmpx <- seq(min(s), max(s), 0.01)
lines(x=tmpx, y=model$Mcdf$RMSD(tmpx),col="red")


#' pripravime si vysvetlujuce premenne do [0,1] pre pouzitie pri predpovedani
Predictors01 <- mapply(function(x,y) x(y), model$Mcdf[-clen], train[-clen])
Predictors01E <- mapply(function(x,y) x(y), model$Mcdf[-clen], test[-clen])
#treba davat velky pozor na poradie v ktorom su funkcie Mcdf, vzdy preverit!!! ci SVK sa nenachadza v model$Mcdf alebo dfhdp

#'
#' ## zdruzena analyza
#' 

#datT01 <- as.data.frame( apply(train, 2, rank, ties.method="random") / (nT+1) )   
datT01 <- mapply(function(x,y) x(y), model$Mcdf, train)
pairs(datT01[sample(nrow(datT01), 500),])  

#'
#' ### Elipticka trieda kopul
#' 
#library(copula)
#model$Copula$normal <- fitCopula(normalCopula(dim=ncol(datT01), dispstr="un"), as.matrix(datT01))@copula




#' 
#' ### Vine kopula
#' 
library(VineCopula)
model$Copula$vine <- RVineStructureSelect(datT01, indeptest = TRUE, type = 0, cores = 4, selectioncrit = "AIC") 

png("contours.png",width = 1200, height = 1500, pointsize = 16)
#old <- par(mfrow=c(3,3))
#plot(model$Copula$vine, type=1, edge.labels="family-tau", edge.label.cex=1, label.col="black")
contour(model$Copula$vine, cex=2.5,edge.labels="family-tau") 
#par(old)
dev.off()



integrand<-Vectorize(function(x,i){
  return(RVinePDF(c(Predictors01E[i,],model$Mcdf$RMSD(x)),model$Copula$vine)*model$Mden$RMSD(x))
},vectorize.args = "x")

xintegrand<-Vectorize(function(x,i){
  return(x*RVinePDF(c(Predictors01E[i,],model$Mcdf$RMSD(x)),model$Copula$vine)*model$Mden$RMSD(x))
},vectorize.args = "x")


hustotaCondVine <- Vectorize(function(y,i){
#  b<-integrate(function(x) integrand(x,i), lower = -Inf, upper = Inf, abs.tol = 0L)$value
#  b<-1
  return(RVinePDF(c(Predictors01E[i,],model$Mcdf$RMSD(y)),model$Copula$vine)*model$Mden$RMSD(y)) #/b
},vectorize.args = "y")


#' examples of conditional densities
seqx<-seq(0,25,0.1)

bbb<-integrate(function(x) integrand(x,3), lower = -Inf, upper = Inf, abs.tol = 0L)$value
lines(x=seqx,y=hustotaCondVine(seqx,3)/bbb,type="l")

bb<-integrate(function(x) integrand(x,16), lower = -Inf, upper = Inf, abs.tol = 0L)$value
lines(x=seqx,y=hustotaCondVine(seqx,16)/bb,type="l")

b<-integrate(function(x) integrand(x,20), lower = -Inf, upper = Inf, abs.tol = 0L)$value
plot(x=seqx,y=hustotaCondVine(seqx,20)/b,type="l", ylab="f(x_RMSD|x_1,...,x_9)", xlab="", main="Examples of conditional densities")


hist(df$RMSD,xlab="RMSD", main="histogram of RMSD", breaks = 35)

#'
#' ### clasifikacia - Bayes
#'

s<-seq(0,25,0.1)
h<-model$Mden$RMSD(seq(0,25,0.1))
plot(x=s,y=h,type="l")
lokmin<-s[splus2R::peaks(-1*h)]
lokmin
train1<-train[train$RMSD<lokmin,]
train2<-train[train$RMSD>lokmin,]

library(e1071)
trainB<-train
trainB$RMSD[trainB$RMSD < lokmin] <- 0
trainB$RMSD[trainB$RMSD >= lokmin] <- 1

testB<-test
testB$RMSD[testB$RMSD < lokmin] <- 0
testB$RMSD[testB$RMSD >= lokmin] <- 1

Naive_Bayes_Model<-naiveBayes( as.factor(RMSD) ~ ., data=trainB)

Naive_Bayes_Model

trainPred<-predict(Naive_Bayes_Model, testB)
table(trainPred, testB$RMSD)






library(lmom)
sset<-subset(trainB, RMSD == 0)[-clen]
model$Mcdf0 <- lapply(sset, function(x) function(y) cdfgam(y, pelgam(samlmu(x)) ))
model$Mqua0 <- lapply(sset, function(x) function(y) quagam(y, pelgam(samlmu(x)) ))
model$Mden0 <- lapply(sset, function(x) function(y) dgamma(y, shape = pelgam(samlmu(x))["alpha"], scale = pelgam(samlmu(x))["beta"]) )

Hmisc::hist.data.frame(sset, nclass = 20)

s<-sset$F9
plot(s, copula::pobs(s), 
     main="Empiricka a paramericka CDF age", ylab="distribucna funkcia", xlab="kvantil")
tmpx <- seq(min(s), max(s), 1)
lines(x=tmpx, y=model$Mcdf0$F9(tmpx),col="red")


sset<-subset(trainB, RMSD == 1)[-clen]
model$Mcdf1 <- lapply(sset, function(x) function(y) cdfgam(y, pelgam(samlmu(x)) ))
model$Mqua1 <- lapply(sset, function(x) function(y) quagam(y, pelgam(samlmu(x)) ))
model$Mden1 <- lapply(sset, function(x) function(y) dgamma(y, shape = pelgam(samlmu(x))["alpha"], scale = pelgam(samlmu(x))["beta"]) )

#s9<-subset(sset, F9 < 50)
s9<-sset$F9
model$Mcdf1$F9 <- function(y) cdfgno(y, pelgno(samlmu(s9)) )
model$Mqua1$F9 <- function(y) quagno(y, pelgno(samlmu(s9)) )
library(fGarch)
pelest<-pelgno(samlmu(s9))
model$Mden1$F9 <- function(y) dged(y,mean = pelest["xi"],sd = pelest["alpha"],nu = pelest["k"]) 

s<-s9
plot(s, copula::pobs(s), 
     main="Empiricka a paramericka CDF", xlab="kvantil")
tmpx <- seq(min(s), max(s), 1)
lines(x=tmpx, y=model$Mcdf1$F9(tmpx),col="red")

Hmisc::hist.data.frame(sset, nclass = 20)

sset<-subset(trainB, RMSD  == 0)[-clen]
#datT01 <- as.data.frame( apply(sset, 2, rank, ties.method="random") / (nT+1) )
datT01 <- mapply(function(x,y) x(y), model$Mcdf0, sset)
#pairs(datT01)
model$Copula$vine0 <- RVineStructureSelect(datT01, indeptest = TRUE, type = 0, selectioncrit = "AIC", cores = 4) 

sset<-subset(trainB, RMSD == 1)[-clen]
#datT01 <- as.data.frame( apply(sset, 2, rank, ties.method="random") / (nT+1) )   
datT01 <- mapply(function(x,y) x(y), model$Mcdf1, sset)
model$Copula$vine1 <- RVineStructureSelect(datT01, indeptest = TRUE, type = 0, selectioncrit = "AIC", cores = 4) 


#' prior proportions:
P1<-length(trainB$RMSD[trainB$RMSD == 0])/nT
P2<-length(trainB$RMSD[trainB$RMSD == 1])/nT
#P3<-length(df$Wine.Type[df$Wine.Type == 3])/nlen
PK<-c(P1,P2)

#kontrola
P1+P2


PKeVine <- function(i){
  predictors<-df[-clen][i,]
  predictors<-as.numeric(predictors)
  
  predictors01E0 <- mapply(function(x,y) x(y), model$Mcdf0, predictors)
  predictors01E0<-as.numeric(predictors01E0)
  
  predictors01E1 <- mapply(function(x,y) x(y), model$Mcdf1, predictors)
  predictors01E1<-as.numeric(predictors01E1)
  
  a0 <- model$Mden0$F1(predictors[1])*
    model$Mden0$F2(predictors[2])*
    model$Mden0$F3(predictors[3])*
    model$Mden0$F4(predictors[4])*
    model$Mden0$F5(predictors[5])*
    model$Mden0$F6(predictors[6])*
    model$Mden0$F7(predictors[7])*
    model$Mden0$F8(predictors[8])*
    model$Mden0$F9(predictors[9])*
    RVinePDF(predictors01E0,model$Copula$vine0)
  
  a0<-unname(a0)
  citatel<-a0*P1
  
  
  a1 <- model$Mden1$F1(predictors[1])*
    model$Mden1$F2(predictors[2])*
    model$Mden1$F3(predictors[3])*
    model$Mden1$F4(predictors[4])*
    model$Mden1$F5(predictors[5])*
    model$Mden1$F6(predictors[6])*
    model$Mden1$F7(predictors[7])*
    model$Mden1$F8(predictors[8])*
    model$Mden1$F9(predictors[9])*
    RVinePDF(predictors01E1,model$Copula$vine1)
  
  a1<-unname(a1)
  menovatel <- a0*P1 + a1*P2
  
  if(menovatel==0){
    return(NaN)
  }
  
  val1 <- citatel/menovatel #pravdepodobnost ze sme v triede 0
  val2 <- 1 - val1
  
  if(val1>val2){  #ak je pravdepodobnost ze sme v triede 0 vacia ako ze sme v triede 1
    return(0)
  }else{
    return(1)
  }
}

trainB$RMSD[2]
PKeVine(2)


#'
#'  ### clasifikacno-regresne predpovede
#'

#' Bayes classification
s<-seq(0,25,0.1)
i<-1
tmp<-vector(length = nE)


for(i in 1:nE){
  mL<-optimize(function(y) hustotaCondVine(y,i), interval = c(0,5), maximum = TRUE)$maximum
  mU<-optimize(function(y) hustotaCondVine(y,i), interval = c(5,25), maximum = TRUE)$maximum
#  minv<-optimize(function(y) hustotaCondVine(y,i), interval = c(mL$maximum,mU$maximum))$minimum
#  b<-integrate(function(x) integrand(x,i), lower = -Inf, upper = Inf, abs.tol = 0L)$value
#  p1<-integrate(function(x) integrand(x,i), lower = mL-1, upper = mL+1)$value
#  p2<-integrate(function(x) integrand(x,i), lower = mU-2, upper = mU+2)$value

  
  prd<-predict(Naive_Bayes_Model, testB[i,])
  if(prd==0){
    tmp[i]<-mL
  }else{
    tmp[i]<-mU
  }
}

plot(test$RMSD-tmp, main="Vine copula + Bayes classification",type="l", xlab="",ylab="error")
tmpVineBayes<-tmp
varVineBayes<-var(test$RMSD-tmpVineBayes)


#tmp[1:2]
#test$RMSD[1:2]

library("Metrics")
mae(test$RMSD,tmpLin)
varReg

mae(test$RMSD,tmp)
varVineBayes




library(parallel)
cl <- makeCluster(4)
clusterEvalQ(cl, library(VineCopula))
clusterExport(cl, varlist=c("model","PKeVine","testB")) 

#' Vine classification
tmp1<-vector(length = nE)
for(i in 1:nE){
  tmp1[i]<-PKeVine(indE[i])
}

stopCluster(cl)
table(tmp1,testB$RMSD)


tmp<-vector(length = nE)

cl <- makeCluster(4)
clusterEvalQ(cl, library(VineCopula))
clusterExport(cl, varlist=c("PKeVine","hustotaCondVine","tmp")) 

for(i in 1:nE){
  mL<-optimize(function(y) hustotaCondVine(y,i), interval = c(0,5), maximum = TRUE)$maximum
  mU<-optimize(function(y) hustotaCondVine(y,i), interval = c(5,25), maximum = TRUE)$maximum
  #  minv<-optimize(function(y) hustotaCondVine(y,i), interval = c(mL$maximum,mU$maximum))$minimum
  #  b<-integrate(function(x) integrand(x,i), lower = -Inf, upper = Inf, abs.tol = 0L)$value
  #  p1<-integrate(function(x) integrand(x,i), lower = mL-1, upper = mL+1)$value
  #  p2<-integrate(function(x) integrand(x,i), lower = mU-2, upper = mU+2)$value
  
  prd<-PKeVine(indE[i])
  if(prd==0){
    tmp[i]<-mL
  }else{
    tmp[i]<-mU
  }
}

stopCluster(cl)

#test$RMSD[17]
#sx<-seq(0,25,0.1)
#plot(x=sx,y=hustotaCondVine(sx,17),type = "l")
#PKeVine(indE[17])

plot(test$RMSD-VineVine, main="Vine kopula + Vine classification",type="l",xlab="",ylab="error")
VineVine<-tmp
varVine<-var(test$RMSD-VineVine)



#tmp[1:2]
#test$RMSD[1:2]

mae(test$RMSD,tmpLin)
sqrt(varReg)

mae(test$RMSD,tmpVineBayes)
sqrt(varVineBayes)

mae(test$RMSD,VineVine)
sqrt(varVine)
