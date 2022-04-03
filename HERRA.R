#####
#Maela dataset (known h^2 = ~0.8 for test)
#
#vcf="snps.vcf.gz"
#vcftools --extract-FORMAT-info GT --out vcf_matrix_Maela --maf 0.01 --max-missing 1.0 --gzvcf $vcf 
#
## replace '0/0' to 0, and '1/1' to 1
#cut -f3- vcf_matrix_Maela.GT.FORMAT| sed s/"0\/0"/"0"/g | sed s/"1\/1"/"1"/g  > vcf_matrix_Maela.txt
#
#R --slave -f HERRA_Nong.R --args vcf_matrix=vcf_matrix_Maela.txt phe_file=resistances.pheno use_cov=FALSE core=30 out=herra_allSNP_Maela.out boot=TRUE n_boot=100 prevalence=0.1 random_phe=FALSE
#####

library(glmnet)
library(doMC)
library(boot)

args <- commandArgs(trailingOnly = TRUE)
vcf_matrix <- "BP_vcf_matrix_allSNP.txt"
cov_file <- "PHENOTYPE/poppunk_1011_run2_microreact_clusters_ksnp.txt"
phe_file <- "PHENOTYPE/poppunk_binary_cluster_ksnp.txt"
out <- "herra.out"
use_cov <- TRUE
random_phe <- FALSE
core <- 8
boot <- TRUE
n_boot <- 100
prevalence <- 0.1
lambda_list <- "0.006, 0.008, 0.01, 0.0102, 0.0104"

for(i in 1:length(args)){
  temp<-strsplit(args[i],"=")
  temp<-temp[[1]]
  temp1<-temp[1]
  temp2<-temp[2]
  assign(as.character(temp1),temp2)
}
registerDoMC(cores=as.numeric(core))

cat('Select samples that available in phenotype file only!\n')
phenotype= read.table(phe_file, header=T,row.name=1, comment.char='$')#

data= read.table(vcf_matrix, header=T,sep='\t',check.names = FALSE,comment.char='')
df = data.matrix(data[,(colnames(data) %in% rownames(phenotype)) ])

#######################

d=data.matrix(phenotype)
if(random_phe == TRUE){
	library(permute)
	d=data.matrix(d[shuffle(d)])}  #random

if(use_cov == TRUE){
  covariate= read.table(cov_file, header=T,row.name=1,comment.char='$')
	#cc=data.matrix(covariate)
  cc = data.matrix(covariate[(row.names(covariate) %in% rownames(phenotype)), ])
}else{ 
	#cc=data.matrix(sample(rep(seq(2), length = nrow(phenotype)))) 
	cc=data.matrix(rep(0, nrow(phenotype)))
	rownames(cc) = rownames(d)
} #ignore covariate}  #random

log=paste(out,'.out',sep='')
cat('#### HERRA result,', as.character(Sys.time()),'####\n',file=log)
######################
fun.confounders<-function(x,y, pf, a,b,n,lam){
  model<-glmnet(x[a:b,],y[a:b],type.gaussian="naive",lambda=lam,penalty.factor=pf)
  coef<-predict(model,type="coefficients")
  n.conf = sum(pf==0)   ### num of confounders 
  nonzerocoef<-c(seq(1:(n.conf)), which(coef[-(1:(n.conf+1))]!= 0)+n.conf)
  nh = round(n/2)
  a2<-ifelse(a==1,nh+1,1)
  b2<-ifelse(b==nh,n,nh) 
  output1<-lm(y[a2:b2]~x[a2:b2,nonzerocoef])
  print(length(nonzerocoef))
  return(summary(output1)$sigma^2)
}

### main method ####
fun.herra<-function(data,idx,d,cc){  
df=data[idx,]

df=t(df)

cat('#Number of variants before screening :',ncol(df),'\n',file=log,append=TRUE)
cat('#Number of sample :',nrow(df),'\n',file=log,append=TRUE)

#### Step 1: ITRRS, iterative ridge for dimension reduction
lambda.seq = as.numeric(strsplit(lambda_list, ",")[[1]])
#lambda.seq = c(0.01) #selected 

res = matrix(0,length(lambda.seq),4) #AS
col_df = ncol(df)
row_df = nrow(df)
ITR=log2(col_df/row_df)+1
cat('itrrs numbers = ',ITR,'\n')

for (lambdai in 1:length(lambda.seq)){
  xdat.all=NULL
  #for (i in 1:22){ #ITRRS for each chromosome
    i=1
    xdat2=df #AS
    #xdat2=masnps.filter  #AS
    #xdat2=masnps.filter
    #filename = paste("...",i,".geno",sep="") #data of chromosome i
    #xdat = read.table(filename, header=T)
    #xdat2 = data.matrix(xdat)
    #rm(xdat)
    for (itrrs in 1:ITR ) {
      pf<-c(rep(0,ncol(cc)),rep(1,ncol(xdat2)))
      xdat3<-cbind(cc,xdat2)
      lam<-lambda.seq[lambdai]
      if (itrrs>3) {
        cv = cv.glmnet(xdat3,d,family="gaussian",nfold=10,alpha=0,penalty.factor=pf)
        lam<-cv$lambda.min
      }
      model<-glmnet(xdat3,d,family="gaussian",type.gaussian="naive",lambda=lam,alpha=0,penalty.factor=pf)
      coef0<-predict(model,type="coefficients")
      coef = abs(coef0[-(1:(ncol(cc)+1)),])
      inout<-ifelse(coef>median(coef),1,0)
      xdat2<-xdat2[,inout==1]
      print(c(lam,median(coef), dim(xdat2)))    
    }
    xdat.all = cbind(xdat.all,xdat2)  
  #}
  filename = paste("...",lambda.seq[lambdai],".Rdata",sep="")
  save(xdat.all,file=filename)
  #### END of Step 1
  #### Steps 2-4: lasso and OLS 
  xdat3 = cbind(cc, xdat.all)
  n = nrow(xdat3)
  nh = round(n/2)
  nfolds=10
  pf<-c(rep(0,ncol(cc)),rep(1,ncol(xdat.all)))
  cv <- cv.glmnet(xdat3[1:nh,],d[1:nh], nfolds=nfolds,penalty.factor=pf) 
  lam<-cv$lambda.min
  est.sig.e1ls<-fun.confounders(xdat3,d,pf,1,nh,n,lam)
  cv <- cv.glmnet(xdat3[(nh+1):n,],d[(nh+1):n],nfolds=nfolds,penalty.factor=pf) 
  lam<-cv$lambda.min
  est.sig.e2ls<-fun.confounders(xdat3,d,pf,nh+1,n,n,lam)
  est.sig.els<-(est.sig.e1ls+est.sig.e2ls)/2 
  vary = summary(lm(d~cc))$sigma^2
  h2.ls<-(vary - est.sig.els)/vary

  K<- as.numeric(prevalence) #disease prevalence
  p<-mean(d)
  print(K)
  z<-dnorm(qnorm(1-K))^2
  h2.liab<-h2.ls*(K*(1-K))^2/(p*(1-p))/z
  h2.est = c(h2.ls, h2.liab)   
  cat('h2',h2.est,ncol(xdat.all),nrow(xdat.all),'\n')

  res[lambdai,]=c(lambda.seq[lambdai],h2.est,ncol(xdat.all))

}
cat('#Number of variants after screening : ',ncol(xdat.all),'\n',file=log,append=TRUE)
colnames(res)=c('#lambda','h2','liability','n_variants')
res
write.table(res,sep='\t',file=log,append=TRUE,row.names=F,quote=F)
#cat(res,file=out,sep='\t')

cat('#Maximum h2 : ',max(res[,2]),'\n',file=log,append=TRUE)
return(max(res[,2]))
}

######## for boostraping only ########
if(boot==TRUE){
  result = boot(data = df, statistic=fun.herra , R = as.numeric(n_boot), cc=cc, d=d ,parallel="multicore", ncpus=as.numeric(core)) #boot method requires variants in rows, samples in columns
  print(result)
  # The bootstrap p-value can then be approximated by
  #print(sum(abs(result$t[,2]-1) > abs(result$t0[2]-1))/(1+result$R))

  confident_in= boot.ci(boot.out = result,type = c("norm", "basic"))
  cat('########### Boostraping result summary ##########\n',sep='\t',file=log,append=TRUE)  
  cat('# confident interval of h^2 at 95%\n',sep='\t',file=log,append=TRUE)
  cat('# Normal',confident_in$normal[2],confident_in$normal[3],'\n',sep='\t',file=log,append=TRUE)
  cat('# Basic',confident_in$basic[4],confident_in$basic[5],'\n',sep='\t',file=log,append=TRUE)

  pdf(paste(out,'.pdf',sep=''))
  plot(result) 
  dev.off()
  print(confident_in)


}else {
  idx=1:nrow(df)
  fun.herra(df,idx,d=d,cc=cc)
}
########

