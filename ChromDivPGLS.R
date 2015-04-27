##"Chromosomal diversification in tropical reef fishes relates to adult traits"
##Martinez PA, Zurano JP, Amado TF, Betancur-R R, Penone C, Bidau CJ, Jacobina UP
##Submitted to XX on XX/XX/XX
##for information on code please contact Pablo Ariel Martinez (pablo_sc82@hotmail.com) or Caterina Penone (caterina.penone@gmail.com)

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
meantrait<-"Depth_range" #mean of the trait
setrait<-"ES_Depth" #standard error of the trait

##Build result table
resultados<-data.frame("RandomVal"=numeric(),"Arvorenum"=numeric(),"estim"=numeric(),
                       "rsq"=numeric(),"pval"=numeric(),"Aicc"=numeric(),"lambda"=numeric())

conjdat<-list()
variab<-list()

##Function to pick a random value in the interval 
funr <- function(a, b) {runif(1,a-b,a+b)}

##Loop repeating the pgls for 100 trees ("tree replicate") and for 500 trait values ("trait replicate")
counter=1
for (i in 1:500) {
  for (j in 1:length(arbol)){
    tryCatch({ #some models do not converge, this avoids the loop to stop
      
      #choose a random value in [mean-sd,mean+sd]
      dados$variab<-apply(dados[,c(meantrait,setrait)],1,function(x)funr(x[meantrait],x[setrait])) 
      
      #comparative data creation
      conjdat[[i]]<-comparative.data(data=dados, phy=arbol[[i]], names.col='Species_name', vcv=T, vcv.dim=3)
      
      #simple model
      ModeloSimple<- pgls(log(rKD+1)~log(variab), conjdat[[i]], lambda='ML')
      
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
write.table(resultados,file="SimpleModelResults_DepthR.txt")
write.csv(statresults,file="StatSimpleModel_DepthR.csv")



############################Simple models & factor variables##############################
  
##Build result table
resultados<-data.frame("Arvorenum"=numeric(),"estim"=numeric(),"rsq"=numeric(),
                       "pval"=numeric(),"Aicc"=numeric(),"lambda"=numeric())
                       
conjdat<-list()

##Loop repeating the pgls for 100 trees ("tree replicate")
counter=1
for (i in 1:length(arbol)) {
  tryCatch({  #some models do not converge, this avoids the loop to stop
    
    conjdat[[i]]<-comparative.data(data=dados, phy=arbol[[i]], names.col='Species_name', vcv=T, vcv.dim=3)
    
    #simple model    
    ModeloSimple<- pgls(log(rKD+1)~Spawning_Mode, conjdat[[i]], lambda='ML') #example of Spawning_Mode
    
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