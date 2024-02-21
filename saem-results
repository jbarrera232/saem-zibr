library(ZIBR)
library(magrittr)
library(dplyr)
library(tibble)
library(nlme)
library(NBZIMM)

#### Part 1: IBD data (Lee et al., 2015) ####

### Load the raw data
PLEASE.file <- "https://raw.githubusercontent.com/chvlyl/PLEASE/master/1_Data/Raw_Data/MetaPhlAn/PLEASE/G_Remove_unclassfied_Renormalized_Merge_Rel_MetaPhlAn_Result.xls"
PLEASE.raw <- read.table(PLEASE.file,sep='\t',header=TRUE,row.names = 1,
                         check.names=FALSE,stringsAsFactors=FALSE)
taxa.raw <- t(PLEASE.raw)

### Make sure you load the data correctly
cat('samples','taxa',dim(taxa.raw),'\n')
taxa.raw[1:3,1:3]

### Load total non-human read counts
human.read.file <- 'https://raw.githubusercontent.com/chvlyl/PLEASE/master/1_Data/Raw_Data/MetaPhlAn/Human_Reads/please_combo_human_reads.xls'
human.read <- read.table(human.read.file,sep='\t',header=TRUE,
                         row.names=1,stringsAsFactors=FALSE)

### Filter low depth samples (low non human reads)
low.depth.samples <- subset(human.read,NonHumanReads<10000)
low.depth.samples[,1:5]

### Delete these samples from PLEASE data.
rownames(taxa.raw)[which(rownames(taxa.raw) %in% rownames(low.depth.samples))]

### Before deletion
dim(taxa.raw)

### After deletion
taxa.raw <- taxa.raw[-which(rownames(taxa.raw) %in% rownames(low.depth.samples)),]
dim(taxa.raw)

### Filter low abundant bacterial data
filter.index1 <- apply(taxa.raw,2,function(X){sum(X>0)>0.4*length(X)})
filter.index2 <- apply(taxa.raw,2,function(X){quantile(X,0.9)>1})
taxa.filter <- taxa.raw[,filter.index1 & filter.index2]
taxa.filter <- 100*sweep(taxa.filter, 1, rowSums(taxa.filter), FUN="/")
cat('after filter:','samples','taxa',dim(taxa.filter),'\n')
cat(colnames(taxa.filter),'\n')
head(rowSums(taxa.filter))

### 
taxa.data <- taxa.filter
dim(taxa.data)


#### Load sample information
sample.info.file <- 'https://raw.githubusercontent.com/chvlyl/PLEASE/master/1_Data/Processed_Data/Sample_Information/2015_02_13_Processed_Sample_Information.csv'
sample.info <- read.csv(sample.info.file,row.names=1)


#### create covariates, 
#### Time, antiTNF+EEN
reg.cov <-
  data.frame(Sample=rownames(taxa.data),stringsAsFactors = FALSE) %>% 
  left_join(rownames_to_column(sample.info,var = 'Sample'),by='Sample')%>%
  dplyr::filter(Treatment.Specific!='PEN') %>%
  dplyr::select(Sample,Time,Subject,Response,Treatment.Specific) %>%
  group_by(Subject) %>% summarise(count = n()) %>% dplyr::filter(count==4) %>%
  dplyr::select(Subject) %>%
  left_join(rownames_to_column(sample.info,var = 'Sample'),by='Subject') %>%
  mutate(Treat=ifelse(Treatment.Specific=='antiTNF',1,0)) %>%
  dplyr::select(Sample,Subject,Time,Response,Treat) %>%
  dplyr::mutate(Subject=paste('S',Subject,sep='')) %>%
  dplyr::mutate(Time=ifelse(Time=='1',0,ifelse(Time=='2',1,ifelse(Time=='3',4,ifelse(Time=='4',8,NA))))) %>%
  dplyr::mutate(Time.X.Treatment=Time*Treat) %>%
  as.data.frame

### take out first time point
reg.cov.t1   <-  subset(reg.cov,Time==0)
rownames(reg.cov.t1) <- reg.cov.t1$Subject
reg.cov.t234 <-  subset(reg.cov,Time!=0)
reg.cov.t234 <- data.frame(
  baseline.sample=reg.cov.t1[reg.cov.t234$Subject,'Sample'],
  baseline.subject=reg.cov.t1[reg.cov.t234$Subject,'Subject'],
  reg.cov.t234,
  stringsAsFactors = FALSE)

#### Fit ZIBR model on the real data
spe.all <- colnames(taxa.data)
p.species.list.zibr <- list()
p.species.list.saem <- list()
mod.species.list.zibr<-list()
mod.species.list.saem<-list()
for (spe in spe.all){
  #spe = 'g__Collinsella'
  ###### create covariates
  X <- data.frame(
    Baseline=taxa.data[reg.cov.t234$baseline.sample, spe]/100,
    #reg.cov.t234[,c('log.days','Delivery','Delivery.X.log.days')]
    reg.cov.t234[,c('Time','Treat')]
  )
  rownames(X) <- reg.cov.t234$Sample
  Z <- X
  subject.ind <- reg.cov.t234$Subject
  time.ind   <- reg.cov.t234$Time
  ####
  cat(spe,'\n')
  Y <- taxa.data[reg.cov.t234$Sample, spe]/100
  cat('Zeros/All',sum(Y==0),'/',length(Y),'\n')
  ####
  ## estimation with ZIBR
  if (sum(Y>0)<10 | sum(Y==0) <10 | max(Y)<0.01){
    print('skip')
    p.species.list.zibr[[spe]] <- 0
    mod.species.list.zibr[[spe]] <- 0
    p.species.list.saem[[spe]] <- 0
    mod.species.list.saem[[spe]] <- 0
    next
  }else{
    est <- zibr(logistic_cov=X,beta_cov=Z,Y=Y,
                subject_ind=subject.ind,
                time_ind=time.ind,
                quad_n=30,verbose=TRUE)
    p.species.list.zibr[[spe]] <- est$joint_p
    mod.species.list.zibr[[spe]] <- est
    
  ## estimation with SAEM
    
    X.1=Z.1=X[,-1]
    X.2=Z.2=X[,-2]
    X.3=Z.3=X[,-3]
    
    mod0<-saem_zibr(Y,X,Z,id,est$beta_v_est,est$logistic_est_table[,1],est$beta_est_table[,1],
                    232,500,5)
    mod1<-saem_zibr(Y,X.1,Z.1,id,est$beta_v_est,est$logistic_est_table[-2,1],est$beta_est_table[-2,1],
                    232,500,5)
    mod2<-saem_zibr(Y,X.2,Z.2,id,est$beta_v_est,est$logistic_est_table[-3,1],est$beta_est_table[-3,1],
                    232,500,5)
    mod3<-saem_zibr(Y,X.3,Z.3,id,est$beta_v_est,est$logistic_est_table[-4,1],est$beta_est_table[-4,1],
                    232,500,5)
    
    sls<- -2*c(mod1$loglik,mod2$loglik,mod3$loglik)+2*mod0$loglik
    p.species.list.saem[[spe]]<-pchisq(sls,2,lower.tail = F)
    mod.species.list.saem[[spe]] <- mod0
  }
  #break
}

##Fit linear transformed model for bacteroides

spe = 'g__Bacteroides'
###### create covariates
X <- data.frame(
  Baseline=taxa.data[reg.cov.t234$baseline.sample, spe]/100,
  #reg.cov.t234[,c('log.days','Delivery','Delivery.X.log.days')]
  reg.cov.t234[,c('Time','Treat')]
)
rownames(X) <- reg.cov.t234$Sample
Z <- X
subject.ind <- reg.cov.t234$Subject
time.ind   <- reg.cov.t234$Time
####
cat(spe,'\n')
Y <- taxa.data[reg.cov.t234$Sample, spe]/100
cat('Zeros/All',sum(Y==0),'/',length(Y),'\n')
####
## fit linear random effect model with arcsin square transformation on Y
tdata <- data.frame(Y.tran=asin(sqrt(Y)),X,SID=subject.ind)
lme.fit <- lme(Y.tran ~ Baseline + Time + Treat,random=~1| SID, data = tdata)
coef.mat <- summary(lme.fit)$tTable[-1,c(1,5)]
p.species.list.saem[[spe]] <- 1-coef.mat[,2]
p.species.list.zibr[[spe]] <- coef.mat[,2]
####

p.species.zibr <- t(as.data.frame(p.species.list.zibr))
p.species.zibr.adj <-
  rownames_to_column(as.data.frame(p.species.zibr),var = 'Species') %>% 
  mutate(across(-Species,~p.adjust(.x,"fdr"))) %>%
  mutate(Detection=ifelse(Treat<0.05,"Sí","No"))

p.species.saem <- t(as.data.frame(p.species.list.saem))
#p.species.saem[1,]<- 1-coef.mat[,2]
p.species.saem.adj <-
  rownames_to_column(as.data.frame(p.species.saem),var = 'Species') %>% 
  mutate(across(-Species,~p.adjust(.x,"fdr"))) %>%
  mutate(Detection=ifelse(Treat<0.05,"Sí","No"))

#### Part 2: Vaginal microbiome data (Romero et al., 2014) ####

data(Romero)

data1<-Romero$SampleData %>% 
  filter(complete.cases(Age))%>%
  mutate(ID=as.numeric(factor(Subect_ID, levels=unique(Subect_ID))),
         AGE_SC=scale(Age,center=min(Age),scale=(max(Age)-min(Age))))%>%
  mutate(Month=ifelse(pregnant == 1, 7*GA_Days/30, GA_Days/30),
         Mon_Preg=Month*pregnant)%>%
  dplyr::select(Subect_ID,ID,Month,pregnant,AGE_SC,Mon_Preg,Total.Read.Counts)

data2<-Romero$OTU%>%
  filter(complete.cases(Romero$SampleData$Age))%>%
  mutate_all(function(x) x/data1$Total.Read.Counts)%>%
  dplyr::select(where(function(x) sum(x==0)/length(x)>=0.1 & sum(x==0)/length(x)<=0.9))

taxa.zibr<-colnames(data2)

## Model 1: Pregnancy and interaction time*pregnancy
# Pregnancy effect

X.tx=Z.tx=data1[,c(3:6)]
X.tx1=Z.tx1=data1[,c(3,5:6)]
id.tx=data1$ID
res2_saem=rep(0,length(taxa.zibr))

for(i in 1:length(taxa.zibr))
{
  Y.tx=data2[,i]
  mod0=saem_zibr(Y.tx,X.tx,Z.tx,id.tx,runif(1,10,20),runif(5,-0.5,0.5),runif(5,-0.5,0.5),
                 232,500,4)
  mod1=saem_zibr(Y.tx,X.tx1,Z.tx1,id.tx,runif(1,10,20),runif(4,-0.5,0.5),runif(4,-0.5,0.5),
                 232,500,4)
  
  tst<- -2*mod1$loglik+2*mod0$loglik
  pval<- pchisq(tst,2,lower.tail = F)
  res2_saem[i]<-pval
}

## Time*pregnancy effect

X.tx1=Z.tx1=data1[,c(3:5)]
res3_saem=rep(0,length(taxa.zibr))

for(i in 1:length(taxa.zibr))
{
  Y.tx=data2[,i]
  mod0=saem_zibr(Y.tx,X.tx,Z.tx,id.tx,runif(1,10,20),runif(5,-0.5,0.5),runif(5,-0.5,0.5),
                 232,500,4)
  mod1=saem_zibr(Y.tx,X.tx1,Z.tx1,id.tx,runif(1,10,20),runif(4,-0.5,0.5),runif(4,-0.5,0.5),
                 232,500,4)
  
  tst<- -2*mod1$loglik+2*mod0$loglik
  pval<- pchisq(tst,2,lower.tail = F)
  res3_saem[i]<-pval
}

## Model 2: Only pregnancy

X.tx=Z.tx=data1[,c(3:5)]
X.tx1=Z.tx1=data1[,c(3,5)]
res4_saem=rep(0,length(taxa.zibr))

for(i in 1:length(taxa.zibr))
{
  Y.tx=data2[,i]
  mod0=saem_zibr(Y.tx,X.tx,Z.tx,id.tx,runif(1,10,20),runif(4,-0.5,0.5),runif(4,-0.5,0.5),
                 232,500,4)
  mod1=saem_zibr(Y.tx,X.tx1,Z.tx1,id.tx,runif(1,10,20),runif(3,-0.5,0.5),runif(3,-0.5,0.5),
                 232,500,4)
  
  tst<- -2*mod1$loglik+2*mod0$loglik
  pval<- pchisq(tst,2,lower.tail = F)
  res4_saem[i]<-pval
}

## Results

p.species.adj2<-data.frame(Species=taxa.zibr,
                           Pregnancy=res2_saem,
                           Preg.Time=res3_saem,
                           Pregnancy2=res4_saem) %>% 
  mutate(Detec1=ifelse(Pregnancy<0.05,T,F),
         Detec2=ifelse(Preg.Time<0.05,T,F),
         Detec3=ifelse(Pregnancy2<0.05,T,F))

colMeans(p.species.adj2[,5:7])
