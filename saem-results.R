library(ZIBR)
library(magrittr)
library(dplyr)
library(tibble)
library(nlme)
library(NBZIMM)

sessionInfo()
#R version 4.4.0 (2024-04-24 ucrt)
#Platform: x86_64-w64-mingw32/x64
#Running under: Windows 11 x64 (build 22631)

#Matrix products: default

#locale:
#[1] LC_COLLATE=Spanish_Ecuador.utf8  LC_CTYPE=Spanish_Ecuador.utf8    LC_MONETARY=Spanish_Ecuador.utf8 LC_NUMERIC=C                    
#[5] LC_TIME=Spanish_Ecuador.utf8    

#time zone: America/Santiago
#tzcode source: internal

#attached base packages:
#[1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
# [1] gtsummary_2.0.2     lubridate_1.9.3     forcats_1.0.0       stringr_1.5.1       purrr_1.0.2         readr_2.1.5         tidyr_1.3.1        
# [8] ggplot2_3.5.1       tidyverse_2.0.0     nlme_3.1-164        numDeriv_2016.8-1.1 boot_1.3-30         MASS_7.3-60.2       patchwork_1.2.0    
#[15] tibble_3.2.1        dplyr_1.1.4         magrittr_2.0.3      ZIBR_1.0.2          statmod_1.5.0       saemix_3.3          npde_3.5           
#[22] lme4_1.1-35.3       Matrix_1.7-0       

#loaded via a namespace (and not attached):
# [1] gtable_0.3.5        xfun_0.47           lattice_0.22-6      tzdb_0.4.0          vctrs_0.6.5         tools_4.4.0         generics_0.1.3     
# [8] fansi_1.0.6         pkgconfig_2.0.3     RColorBrewer_1.1-3  gt_0.11.0           readxl_1.4.3        lifecycle_1.0.4     farver_2.1.2       
#[15] compiler_4.4.0      munsell_0.5.1       sass_0.4.9          htmltools_0.5.8.1   yaml_2.3.8          pillar_1.9.0        nloptr_2.0.3       
#[22] crayon_1.5.3        mclust_6.1.1        commonmark_1.9.1    tidyselect_1.2.1    digest_0.6.35       stringi_1.8.3       labeling_0.4.3     
#[29] splines_4.4.0       fastmap_1.1.1       grid_4.4.0          colorspace_2.1-0    cli_3.6.2           cards_0.2.2         utf8_1.2.4         
#[36] withr_3.0.1         scales_1.3.0        timechange_0.3.0    rmarkdown_2.28      igraph_2.0.3        gridExtra_2.3       cellranger_1.1.0   
#[43] hms_1.1.3           evaluate_0.24.0     knitr_1.48          markdown_1.13       rlang_1.1.3         Rcpp_1.0.12         glue_1.7.0         
#[50] xaringanExtra_0.8.0 xml2_1.3.6          pkgload_1.4.0       rstudioapi_0.16.0   minqa_1.2.6         R6_2.5.1            fs_1.6.4           

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
    
   mod0<-saem_zibr(Y,X,Z,subject.ind,est$beta_v_est,est$logistic_est_table[,1],est$beta_est_table[,1],
                    86,500,5)
    mod1<-saem_zibr(Y,X.1,Z.1,subject.ind,est$beta_v_est,est$logistic_est_table[-2,1],est$beta_est_table[-2,1],
                    86,500,5)
    mod2<-saem_zibr(Y,X.2,Z.2,subject.ind,est$beta_v_est,est$logistic_est_table[-3,1],est$beta_est_table[-3,1],
                    86,500,5)
    mod3<-saem_zibr(Y,X.3,Z.3,subject.ind,est$beta_v_est,est$logistic_est_table[-4,1],est$beta_est_table[-4,1],
                    86,500,5)
    
    sls<- -2*c(mod1$loglik,mod2$loglik,mod3$loglik)+2*mod0$loglik
    p<-pchisq(sls,2,lower.tail = F)
    
    p.species.list.saem[[spe]] <- p
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
  mutate(Detection=ifelse(Treat<0.05,"Yes","No"))

p.species.saem <- t(as.data.frame(p.species.list.saem))
p.species.saem.adj <-
  rownames_to_column(as.data.frame(p.species.saem),var = 'Species') %>% 
  mutate(across(-Species,~p.adjust(.x,"fdr"))) %>%
  mutate(Detección=ifelse(Treat<0.05,"Yes","No"))

#### Part 2: Vaginal microbiome data (Romero et al., 2014) ####

data(Romero)

## Cleaning covariate data

data1<-Romero$SampleData %>% 
  filter(complete.cases(Age))%>%
  mutate(ID=as.numeric(factor(Subect_ID, levels=unique(Subect_ID))),
         AGE_SC=scale(Age,center=min(Age),scale=(max(Age)-min(Age))))%>%
  mutate(Month=ifelse(pregnant == 1, 7*GA_Days/30, GA_Days/30),
         Time=scale(Month,center=min(Month),scale=(max(Month)-min(Month))),
         Time_Preg=Time*pregnant)%>%
  dplyr::select(Subect_ID,ID,Time,pregnant,AGE_SC,Time_Preg,Total.Read.Counts,Month)

# Preparing compositional microbiome data

data2<-Romero$OTU%>%
  filter(complete.cases(Romero$SampleData$Age))%>%
  mutate_all(function(x) x/data1$Total.Read.Counts)%>%
  dplyr::select(where(function(x) sum(x==0)/length(x)>=0.1 & sum(x==0)/length(x)<=0.9))

# Discading bacteria taxa that are not present in both pregnant and non-pregnant women
                      
taxa.out<-c(31,49,50,60)
taxa.def<-setdiff(1:61,taxa.out)

## Model 1: Age, time, pregnancy
# Full model

X.tx=Z.tx=data1[,c(3:5)]
id.tx<-data1$ID

MOD1_0<-list()
 
for(tax in taxa.def)
{
  set.seed(232)
  Y.tx=data2[,tax]
  nc=ncol(X.tx)
  mod0=try(saem_zibr(Y.tx,X.tx,Z.tx,id.tx,runif(1,10,20),runif(nc+1,-0.5,0.5),runif(nc+1,-0.5,0.5),
                 232,500,5))  
  MOD1_0[[tax]]<-mod0    
}

# Model without pregnancy

X.tx=Z.tx=data1[,c(3,5)]

MOD1_1<-list()

for(tax in taxa.def)
{
  set.seed(232)
  Y.tx=data2[,tax]
  nc=ncol(X.tx)
  mod0=try(saem_zibr(Y.tx,X.tx,Z.tx,id.tx,runif(1,10,20),runif(nc+1,-0.5,0.5),runif(nc+1,-0.5,0.5),
                     232,500,5))
  MOD1_1[[tax]]<-mod0
}

## Model 2: age, time, pregnancy, interaction

# Full model

X.tx=Z.tx=data1[,c(3:6)]

MOD2_0<-list()

for(tax in taxa.def)
{
  set.seed(232)
  Y.tx=data2[,tax]
  nc=ncol(X.tx)
  mod0=try(saem_zibr(Y.tx,X.tx,Z.tx,id.tx,runif(1,10,20),runif(nc+1,-0.5,0.5),runif(nc+1,-0.5,0.5),
                     232,500,5))
  MOD2_0[[tax]]<-mod0
}

# Model without pregnancy

X.tx=Z.tx=data1[,c(3,5:6)]
                      
MOD2_1<-list()

for(tax in taxa.def)
{
  set.seed(232)
  Y.tx=data2[,tax]
  nc=ncol(X.tx)
  mod0=try(saem_zibr(Y.tx,X.tx,Z.tx,id.tx,runif(1,10,20),runif(nc+1,-0.5,0.5),runif(nc+1,-0.5,0.5),
                     232,500,5))
  MOD2_1[[tax]]<-mod0
}

# Model without interaction

X.tx=Z.tx=data1[,3:5]

MOD2_2<-list()
                      
for(tax in taxa.def)
{
  set.seed(232)
  Y.tx=data2[,tax]
  nc=ncol(X.tx)
  mod0=try(saem_zibr(Y.tx,X.tx,Z.tx,id.tx,runif(1,10,20),runif(nc+1,-0.5,0.5),runif(nc+1,-0.5,0.5),
                     232,500,5))
  MOD2_2[[tax]]<-mod0
}

## Resultas

res1<-sapply(MOD1_0,function(x) x$loglik)
res2<-sapply(MOD1_1,function(x) x$loglik)
res3<-sapply(MOD2_0,function(x) x$loglik)
res4<-sapply(MOD2_1,function(x) x$loglik)
res5<-sapply(MOD2_2,function(x) x$loglik)

p.species.romero<-data.frame(Species=taxa.zibr[taxa.def],
                             LLMOD1=res1,
                             LLMOD1_1=res2,
                             LLMOD2=res3,
                             LLMOD2_1=res4,
                             LLMOD2_2=res5) %>% 
  #mutate(across(-Species,~p.adjust(.x,"fdr"))) %>% 
  mutate(pval_Preg1=pchisq(2*(LLMOD1-LLMOD1_1),2,lower.tail=F),
         pval_Preg2=pchisq(2*(LLMOD2-LLMOD2_1),2,lower.tail=F),
         pval_Inter=pchisq(2*(LLMOD2-LLMOD2_2),2,lower.tail=F)) %>% 
  mutate(Detec_Preg1=ifelse(pval_Preg1<0.05,T,F),
         Detec_Preg2=ifelse(pval_Preg2<0.05,T,F),
         Detec_Inter=ifelse(pval_Inter<0.05,T,F))  
