
library(purrr)
library(caret)
library(scales)
library(truncnorm)

#STEP 1: SIMULATE A REAL POPULATION 
simulate.pop<-function(pop.size,mean.group,tot.areas){
  tot.pop<-rpois(1,pop.size)
  ibex.list<-c(1:tot.pop)
  ibex.list<-as.character(ibex.list)
  
  #split animals at random into the given n sectors
  ibex.areas<-split(sample(ibex.list),sample(1:tot.areas,length(ibex.list),
                                             replace=TRUE,prob=area.sizes))
  
  #split individuals into groups  of average size k
  k=0
  ibex.pop<-list()
  for(i in 1:tot.areas){
    k<-rpois(1,length(unlist(ibex.areas[i]))/mean.group)
    ibex.pop[[i]]<-split(sample(unlist(ibex.areas[i],use.names=FALSE)),
                         sample(1:k,length(unlist(ibex.areas[i])),replace=TRUE))
    names(ibex.pop[[i]])<-c(paste("g",1:length(ibex.pop[[i]]),sep=""))
  }
  
  return(list("ibex.pop"=ibex.pop,
              "tot.pop"=tot.pop))
}

# STEP 2 : PERFORM TOTAL BLOCK COUNTS 
block.count<-function(ibex.pop,tot.areas,det.p1,intra.group.det,det.group,tot.pop,detvar){
  area.count<-c()
  area.groups<-c()
  area.meansize<-c()
  census.pop.est=0
  ibex.groups=list()
  sizes=0
  
  #detectability is different between each sector
  for(t in 1:tot.areas){
    ibex.groups<-ibex.pop[[t]]
    group.sizes=c(1:length(ibex.groups))
    for(h in 1:length(ibex.groups)){
      group.sizes[h]<-sqrt(length(unlist(ibex.groups[h])))
    }
    group.sizes[length(ibex.groups)+1]<-1
    for (h in 1:length(group.sizes)){
      ifelse(group.sizes[h]>=sqrt(15),group.sizes[h]<-sqrt(15),group.sizes[h]<-group.sizes[h])
    }
    group.det<-rescale(group.sizes,c(det.group,1))
    
    #factor to correct group size 
    factor<-rnorm(1,intra.group.det,sd=0.05)
    
    #factor to correct det.p between sectors
    detp.factor<-rtruncnorm(1,0,1,det.p1,detvar)
    
    #observer perform the count detecting only some groups
    obs1.det<-c(names(ibex.groups))
    for(i in 1:length(obs1.det)){
      det.p<-detp.factor*group.det[i]
      ind<-rbernoulli(1,det.p)
      ifelse(ind==TRUE,obs1.det[i]<-obs1.det[i],obs1.det[i]<-0)
    } 
    obs1.det<-obs1.det[obs1.det!=0]
    area.count[t]<-round((sum(unlist(lapply(
      ibex.groups[obs1.det],length)))*factor),digits=0)
    area.meansize[t]<-mean(unlist(lapply(ibex.groups[obs1.det],
                                         length)))*factor
    area.groups[t]<-round(length(obs1.det),digits=0)
    
  }
  area.meansize[is.na(area.meansize)]<-0
  
  #ratio of estimated individuals with block counts 
  census.pop.est<-sum(area.count)
  census.ratio=census.pop.est/tot.pop
  
  return(list("census.ratio"=census.ratio,
              "census.abundance"=census.pop.est,
              "area.group"=area.groups,
              "area.meansize"=area.meansize,
              "area.count"=area.count))
}

#STEP 3: ESTIMATE POPULATION SIZE WITH DO METHODS 

estimate<-function(det.p1,pop.size,
                   mean.group,tot.areas,det.group,
                   intra.group.det,surv.areas,obs.effect,detvar){
  
  #simulate the real population with simulate.pop function
  ibex.pop<-simulate.pop(pop.size,mean.group,tot.areas)$ibex.pop
  
  #obtain "real" population size 
  tot.pop<-length(unlist(ibex.pop))
  
  #perform block counts 
  groundcount<-block.count(ibex.pop,tot.areas,det.p1,intra.group.det,det.group,tot.pop,detvar)
  census.ratio<-groundcount$census.ratio
  census.total<-groundcount$census.abundance
  
  #perform the full DO  
  base.sample=matrix(0,nrow=tot.areas,ncol=1)
  group.detpSR=matrix(0,nrow=tot.areas,ncol=1)
  ind.detpSR=matrix(0,nrow=tot.areas,ncol=1)
  mean.size=matrix(0,nrow=tot.areas,ncol=1)
  
  for(m in 1:tot.areas){
    ibex.groups<-ibex.pop[[m]]
    if(length(ibex.groups)>0){
      group.sizes=c(1:length(ibex.groups))
      for(h in 1:length(ibex.groups)){
        group.sizes[h]<-sqrt(length(unlist(ibex.groups[h])))}
      group.sizes[length(ibex.groups)+1]<-1
      for (h in 1:length(group.sizes)){
        ifelse(group.sizes[h]>=sqrt(15),group.sizes[h]<-sqrt(15),group.sizes[h]<-group.sizes[h])}
      group.det<-rescale(group.sizes,c(det.group,1))
      
      detp.factor<-rtruncnorm(1,0,1,det.p1,detvar)
      
      #observer 1 perform the census
      obs1.det<-c(names(ibex.groups))
      for(i in 1:length(obs1.det)){
        det.p<-detp.factor*group.det[i]
        ind<-rbernoulli(1,det.p)
        ifelse(ind==TRUE,obs1.det[i]<-obs1.det[i],obs1.det[i]<-0)
      } 
      obs1.det<-obs1.det[obs1.det!=0]
      
      #obs2 perform the census 
      obs2.det<-c(names(ibex.groups))
      for(i in 1:length(obs2.det)){
        det.p<-detp.factor*group.det[i]*obs.effect
        ind<-rbernoulli(1,det.p)
        ifelse(ind==TRUE,obs2.det[i]<-obs2.det[i],obs2.det[i]<-0)
      } 
      obs2.det<-obs2.det[obs2.det!=0]
      
    }else{
      obs1.det<-0
      obs1.det<-obs1.det[obs1.det!=0]
      obs2.det<-0
      obs2.det<-obs2.det[obs2.det!=0]
    }
    
    # DO data analysis
    tot<-union(obs1.det,obs2.det)
    factor<-rnorm(1,intra.group.det,sd=0.05)
    if(length(tot)!=0){
      mean.size[m]<-mean(unlist(lapply(ibex.groups[tot],length)))*factor
    }else
    {mean.size[m]=0}
    
    est.groups=0
    if(length(obs1.det)>1 && length(obs2.det)>0){
      both.seen<-intersect(obs1.det,obs2.det)
      only.obs1<-setdiff(obs1.det,obs2.det)
      only.obs2<-setdiff(obs2.det,obs1.det)
      est.groups<-(((length(both.seen)+length(only.obs1)+1)*
                      (length(both.seen)+length(only.obs2)+1))/ (length(both.seen)+1))-1
      base.sample[m]<-round(est.groups*mean.size[m],0)
      group.detpSR[m]<-((length(obs1.det)/est.groups)+(length(obs2.det)/est.groups))/2
      ind.detpSR[m]<-(((length(unlist(ibex.groups[obs1.det]))*factor)/base.sample[m])+
                        ((length(unlist(ibex.groups[obs2.det]))*factor)/base.sample[m]))/2
    }
  }
  
  totalDO<-sum(base.sample)/tot.pop
  sdDO.groupdetect<-sd(group.detpSR[group.detpSR!=0])
  sdDO.inddetect<-sd(ind.detpSR[ind.detpSR!=0])
  
  
  ## perform the DOAS with only some sectors (areas)
  areas<-sample(1:tot.areas,surv.areas,replace=FALSE)
  group.detpSR<-group.detpSR[areas]
  ind.detpSR<-ind.detpSR[areas]
  group.detpSR[which(is.na(group.detpSR))]<-0 
  if(length(group.detpSR)==0){
    group.detpSR<-0
  }
  ind.detpSR[which(is.na(ind.detpSR))]<-0
  if(length(ind.detpSR)==0){
    ind.detpSR<-0
  }
  estSR.groupdetp<-mean(group.detpSR)
  sdSR.groupdetp<-sd(group.detpSR)
  estSR.inddetp<-mean(ind.detpSR)
  sdSR.inddetp<-sd(ind.detpSR)
  adjSR.areagroup<-rep(0,tot.areas)
  groupadjSR.areacount<-c()
  indadjSR.areacount<-c()
  area.groups<-groundcount$area.group
  area.meansize<-groundcount$area.meansize
  area.count<-groundcount$area.count
  
  for(i in 1:tot.areas){
    adjSR.areagroup[i]<-area.groups[i]/estSR.groupdetp
    if(is.infinite(adjSR.areagroup[i])==TRUE){
      adjSR.areagroup[i]<-0
    }
    groupadjSR.areacount[i]<-adjSR.areagroup[i]*area.meansize[i]
  }
  
  indadjSR.areacount<-census.total/estSR.inddetp
  
  groupadjust.censusSR<-sum(groupadjSR.areacount)/tot.pop
  indadjust.censusSR<-sum(indadjSR.areacount)/tot.pop
  
  
  return(list("census"=census.ratio,
              "totalDO"=totalDO,
              "sdDO.groupdetect"=sdDO.groupdetect,
              "sdDO.inddetect"=sdDO.inddetect,
              "groupadj.censusSR"=groupadjust.censusSR,
              "indadj.censusSR"=indadjust.censusSR,
              "sdSR.groupdetect"=sdSR.groupdetp,
              "sdSR.inddetect"=sdSR.inddetp))
}

#STEP 4: MULTIPLE CENSUS SIMULATIONS 

area.sizes<-read.csv("Zone_sorveglianza_area.csv")
as.matrix(area.sizes)
area.sizes<-as.numeric(unlist(area.sizes[2]))

tot.areas=38
loops<-100                                  # number of simulations 
dets<-c(seq(0.3,0.9,0.1))            #detectabilities
pops=c(2500,3500,5000)           # population size
groups=c(5,10,15)                 #average group size
subareas=c(1:37)                   # number of sectors 
detg=c(0.3,0.5,0.7)               #correction factor 
intragroup=c(0.8,0.9,1)          #miscount effect
obs.effect=c(0.5,0.75,1)        #observer effect
detvar=c(0.1,0.3,0.5)            #variability of det.p between sectors

library(doParallel)
cl <- makeCluster(4,"SOCK")
registerDoParallel(cl)

totaldata<-data.frame()
result<- foreach(i.detp = 1:length(dets),.combine=rbind) %:%
  foreach(i.pop = 1:length(pops),.combine=rbind) %:%
  foreach(i.group = 1:length(groups),.combine=rbind) %:%
  foreach(i.gdet = 1:length(detg),.combine=rbind) %:%
  foreach(i.gsize = 1:length(intragroup),.combine=rbind) %:%
  foreach(i.surv = 1:length(subareas),.combine=rbind) %:%
  foreach(i.obseffect = 1:length(obs.effect),.combine=rbind) %:%
  foreach(i.detvar = 1:length(detvar),.combine=rbind) %:%
  foreach(i = 1:loops,.combine=rbind,
          .packages = c("purrr","caret",
                        "scales","truncnorm"))%dopar% {
                          detp<-dets[i.detp]
                          pop<-pops[i.pop]
                          group<-groups[i.group]
                          survs=subareas[i.surv]
                          gdet<-detg[i.gdet]
                          gsize=intragroup[i.gsize]
                          effect=obs.effect[i.obseffect]
                          detv=detvar[i.detvar]
                          b<-estimate(detp,pop,group,38,gdet,gsize,survs,effect,detv)
                          a<-data.frame("det.p"=detp,
                                        "pop.size"=pop,   
                                        "mean.group"=group,
                                        "det.groups"=gdet,
                                        "gsize"=gsize,
                                        "surv.areas"=survs,
                                        "obs.effect"=effect,
                                        "detvar"=detv,
                                        "sdDO.groupdetect"=b$sdDO.groupdetect,
                                        "sdDO.inddetect"=b$sdDO.inddetect,
                                        "totalDO"=b$totalDO,
                                        "census"=b$census,
                                        "groupadjust.censusSR"=b$groupadj.censusSR,
                                        "indadjust.censusSR"=b$indadj.censusSR,
                                        "sdSR.groupdetect"=b$sdSR.groupdetect,
                                        "sdSR.inddetect"=b$sdSR.inddetect)
                          return(a)
                        }   

totaldata<-rbind(totaldata,result)

write.csv(totaldata,"Data frames/DOAS_methods.csv")

parallel::stopCluster(cl)
