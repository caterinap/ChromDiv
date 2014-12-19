##"Chromosomal diversification in tropical reef fishes relates to larval and adult traits"
##Martinez PA, Zurano JP, Amado TF, Penone C, Bidau CJ, Jacobina UP
#submitted to XX
#for information on code contact Pablo Ariel Martinez (pablo_sc82@hotmail.com) or Caterina Penone (caterina.penone@gmail.com)

##Load required packages
library("ape")
library("geiger")
library("caper")
library("mvtnorm")
library("MASS")
library("picante")

##Read tree from Near et al 2012
arbol<-read.nexus("~/multiarbol.txt", tree.names=TRUE)

##Read data table
dados<-read.table("~/Dataset.txt", header=TRUE)
head(dados)

############################Simple models & continuous variables##############################
##Example of depth range
meantrait<-"Mean_PLD"
setrait<-"ES_PLD"

##Build result table
resultados<-data.frame("Arvorenum"=rep(NA,1001),"RandomVal"=rep(NA,1001),"estim"=rep(NA,1001),
                       "rsq"=rep(NA,1001),"pval"=rep(NA,1001),"Aicc"=rep(NA,1001),"lambda"=rep(NA,1001))

conjdat<-list()

##Function to pick a random value in the interval 
funr <- function(a, b) {runif(1,a-b,a+b)}

##Choose 100 trees at random in the nexus file
Ntree<-sample(length(arbol),100,replace=F)

##Loop repeating the pgls for 100 trees ("tree replicate") and for 500 trait values ("trait replicate")
counter=1
for (i in 1:length(Ntree)) {
  for (j in 1:500){
    tryCatch({ #some models do not converge, this avoids the loop to stop
      
      #choose a random value in [mean-sd,mean+sd]
      dados$variab<-apply(dados[,c(meantrait,setrait)],1,function(x)funr(x[meantrait],x[setrait])) 
      
      #comparative data creation
      conjdat[[i]]<-comparative.data(data=dados, phy=arbol[[i]], names.col='Family', vcv=T, vcv.dim=3)
      
      #simple model
      ModeloSimple<- pgls(log(rKD)~log(variab), conjdat[[i]], lambda='ML')
      
      #extract model coefficients
      estim<-summary(ModeloSimple)$coef[2,1]
      resq<-summary(ModeloSimple)$r.squared
      pval<-summary(ModeloSimple)$coef[2,4]
      Aicc<-ModeloSimple$aicc
      lambda<-as.numeric(ModeloSimple$param.CI$lambda$opt)
      
      #write in a table
      resultados[counter,1]<- i
      resultados[counter,2]<- j    
      resultados[counter,3]<- estim
      resultados[counter,4]<- resq
      resultados[counter,5]<- pval
      resultados[counter,6]<- Aicc
      resultados[counter,7]<- lambda
      counter=counter+1
      print(paste("tree",i,"rep",j))
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
}

##Calculate statistics for each parameter
#variation due to tree choice (tree replicate)
mean_by_tree<-aggregate(.~Arvorenum, data=resultados, mean)
#variation due to continuous trait choice (trait replicate)
mean_by_randomval<-aggregate(.~RandomVal, data=resultados, mean)

statresults<-data.frame(mean_all=apply(resultados,2,mean),
                        sd_all=apply(resultados,2,sd),
                        sd_arvore=apply(mean_by_tree,2,sd),
                        sd_interval=apply(mean_by_randomval,2,sd))[-(1:2),]

#confidence interval
statresults$ICleft<-statresults$mean_all - qnorm(0.975)*statresults$sd_all/sqrt(nrow(resultados))
statresults$ICright<-statresults$mean_all + qnorm(0.975)*statresults$sd_all/sqrt(nrow(resultados))

##Save data
write.table(resultados,file="SimpleModelResults_Pld.txt")
write.csv(statresults,file="StatSimpleModel_Pld.csv")



############################Simple models & factor variables##############################
##Example of Spawning_Mode

##Build result table
resultados<-data.frame("Arvorenum"=rep(0,100),"estim"=rep(0,100),"rsq"=rep(0,100),"pval"=rep(0,100),
                       "Aicc"=rep(0,100), "lambda"=rep(NA,100))

conjdat<-list()

##Choose 100 trees at random in the nexus file
Ntree<-sample(length(arbol),100,replace=F)

##Loop repeating the pgls for 100 trees ("tree replicate")
counter=1
for (i in 1:length(Ntree)) {
  tryCatch({  #some models do not converge, this avoids the loop to stop
    
    conjdat[[i]]<-comparative.data(data=dados, phy=arbol[[i]], names.col='Family', vcv=T, vcv.dim=3)
    
    #simple model    
    ModeloSimple<- pgls(log(rKD)~Spawning_Mode, conjdat[[i]], lambda='ML') ##example of Spawning_Mode
    
    #extract model coefficients
    estim<-summary(ModeloSimple)$coef[2,1]
    resq<-summary(ModeloSimple)$r.squared
    pval<-summary(ModeloSimple)$coef[2,4]
    Aicc<-ModeloSimple$aicc
    lambda<-as.numeric(ModeloSimple$param.CI$lambda$opt)
    
    #write in a table
    resultados[counter,1]<- i
    resultados[counter,2]<- estim
    resultados[counter,3]<- resq
    resultados[counter,4]<- pval
    resultados[counter,5]<- Aicc
    resultados[counter,6]<- lambda
    counter=counter+1
    print(paste("tree",i))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

##Calculate statistics for each parameter
statresults<-data.frame(mean_all=apply(resultados,2,mean,na.rm=T),
                        sd_all=apply(resultados,2,sd,na.rm=T),
                        sd_arvore=apply(resultados,2,sd,na.rm=T))[-1,]


#confidence interval
statresults$ICleft<-statresults$mean_all - qnorm(0.975)*statresults$sd_all/sqrt(nrow(resultados))
statresults$ICright<-statresults$mean_all + qnorm(0.975)*statresults$sd_all/sqrt(nrow(resultados))

##Save data
write.table(resultados,file="SimpleModelResults_Spawn.txt")
write.csv(statresults,file="StatSimpleModel.csv")


############################Best model#########################################

##Build result table
resultados<-data.frame("Arvorenum"=rep(NA,1001),"RandomVal"=rep(NA,1001),"estim_DepthR"=rep(NA,1001),
                       "estim_Spawn"=rep(NA,1001),"estim_Dial"=rep(NA,1001),"estim_inter"=rep(NA,1001),
                       "pval_DepthR"=rep(NA,1001),"pval_spawn"=rep(NA,1001),"pval_Dial"=rep(NA,1001),
                       "pval_inter"=rep(NA,1001),"rsq"=rep(NA,1001),"pval"=rep(NA,1001), "Aicc"=rep(NA,1001), 
                       "lambda"=rep(NA,1001),"Fstat"=rep(NA,1001))

##Continuous variable of the model: depth range
meantrait<-"Depth_range"
setrait<-"ES_Depth"

conjdat<-list()

##Function to pick a random value in the interval 
funr <- function(a, b) {runif(1,a-b,a+b)}

##Choose 100 trees at random in the nexus file
Ntree<-sample(length(arbol),100,replace=F)

counter=1
for (i in 1:length(Ntree)) {
  for (j in 1:500){ 
    tryCatch({   #some models do not converge, this avoids the loop to stop
      
      conjdat[[i]]<-comparative.data(data=dados, phy=arbol[[i]], names.col='Family', vcv=T, vcv.dim=3)
      
      #choose a random value in [mean-sd,mean+sd] only for continuous variables
      dados$DR<-apply(dados[,c(meantrait,setrait)],1,function(x) funr(x[meantrait],x[setrait]))
      
      #best multiple model
      Modelo1<- pgls(log(rKD)~log(DR)*(Spawning_Mode)+Dial_Activity, conjdat[[i]], lambda='ML')
      
      #extract model coefficients
      estim_DepthR<-summary(Modelo1)$coef[2,1]
      estim_Spawn<-summary(Modelo1)$coef[3,1]
      estim_Dial<-summary(Modelo1)$coef[4,1]
      estim_inter<-summary(Modelo1)$coef[5,1]
      pval_DepthR<-summary(Modelo1)$coef[2,4]
      pval_spawn<-summary(Modelo1)$coef[3,4]
      pval_Dial<-summary(Modelo1)$coef[4,4]
      pval_inter<-summary(Modelo1)$coef[5,4]
      rsq<-summary(Modelo1)$r.squared
      pval<-pf(summary(Modelo1)$fstatistic[1],summary(Modelo1)$fstatistic[2], summary(Modelo1)$fstatistic[3],lower.tail = FALSE)
      Aicc<-Modelo1$aicc
      lambda<-as.numeric(ModeloSimple$param.CI$lambda$opt)
      Fstat<-summary(Modelo1)$fstatistic[1]
      
      #write in a table
      resultados[counter,1]<- i
      resultados[counter,2]<- j 
      resultados[counter,3]<- estim_DepthR
      resultados[counter,4]<- estim_Spawn
      resultados[counter,5]<- estim_Dial
      resultados[counter,6]<- estim_inter
      resultados[counter,7]<- pval_DepthR
      resultados[counter,8]<- pval_spawn
      resultados[counter,9]<- pval_Dial
      resultados[counter,10]<- pval_inter
      resultados[counter,11]<- rsq
      resultados[counter,12]<- pval
      resultados[counter,13]<- Aicc
      resultados[counter,14]<- lambda
      resultados[counter,15]<- Fstat
      
      counter=counter+1
      print(paste("tree",i,"rep",j))
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
}

##Calculate statistics for each parameter
#variation due to tree choice (tree replicate)
mean_by_tree<-aggregate(.~Arvorenum, data=resultados, mean,na.rm=T)
#variation due to continuous trait (trait replicate)
mean_by_randomval<-aggregate(.~RandomVal, data=resultados, mean,na.rm=T)

statresults<-data.frame(mean_all=apply(resultados,2,mean,na.rm=T),
                        sd_all=apply(resultados,2,sd,na.rm=T),
                        sd_arvore=apply(mean_by_tree,2,sd,na.rm=T),
                        sd_interval=apply(mean_by_randomval,2,sd))[-(1:2),]

#confidence interval
statresults$ICleft<-statresults$mean_all - qnorm(0.975)*statresults$sd_all/sqrt(nrow(resultados))
statresults$ICright<-statresults$mean_all + qnorm(0.975)*statresults$sd_all/sqrt(nrow(resultados))

##Save data
write.table(resultados,file="CompleteModelResults.txt")
write.csv(statresults,file="StatCompleteModel.csv")
