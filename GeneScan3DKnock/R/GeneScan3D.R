#' @importFrom stats binomial dbeta gaussian glm pcauchy pchisq rbinom sd var median
utils::globalVariables(c('G_Enhancer1_surround','G_Enhancer2_surround','variants_Enhancer1_surround','variants_Enhancer2_surround','Enhancer1.pos','Enhancer2.pos','G_EnhancerAll','Z_EnhancerAll','p_EnhancerAll',"G_gene_buffer", "Z_gene_buffer", "pos_gene_buffer",'n','G_promoter','Z_promoter','G_Enhancer1','Z_Enhancer1','G_Enhancer2','Z_Enhancer2','qchisq'))

#' Data example for GeneScan3D (gene-based testing by integrating long-range chromatin interactions).
#'
#'This simulated example dataset contains outcome variable Y, covariate X, genotype and functional annotation matrices for gene and buffer region, promoter and two enhancers, positions of genetic variants in gene and buffer region.
#'
#'We generated genotypes for 2,000 individuals in a 14.5 Kb gene region, a promoter as a 0.5 Kb segment upstream of the TSS (the start point of the gene) and R = 2 enhancers with length 2 KB, which are outside the 15 Kb gene plus promoter region.
#'
#' @name Example.GeneScan3D
#' @docType data
#' @keywords data
#' @usage data("GeneScan3D.example")
#' @examples
#' data("GeneScan3D.example")
#'
#'Y=GeneScan3D.example$Y; X=GeneScan3D.example$X; n=length(Y)
#'
#'G_gene_buffer=GeneScan3D.example$G_gene_buffer
#'G_promoter=GeneScan3D.example$G_promoter
#'G_EnhancerAll=cbind(GeneScan3D.example$G_Enhancer1,GeneScan3D.example$G_Enhancer2)
#'
#'Z_gene_buffer=GeneScan3D.example$Z_gene_buffer
#'Z_promoter=GeneScan3D.example$Z_promoter
#'Z_EnhancerAll=rbind(GeneScan3D.example$Z_Enhancer1,GeneScan3D.example$Z_Enhancer2)
#'
#'pos_gene_buffer=GeneScan3D.example$pos_gene_buffer
"GeneScan3D.example"

#' The preliminary data management for GeneScan3DKnock.
#'
#' This function does the preliminary data management and fit the generalized linear model under null hypothesis for unrelated samples. The output will be used in the other GeneScan functions.
#'
#' @param Y The outcome variable, an n*1 matrix where n is the number of individuals.
#' @param X An n*d covariates matrix where d is the number of covariates.
#' @param id The subject id. This is used to match phenotype with genotype. The default is NULL, where the matched phenotype and genotype matrices are assumed.
#' @param out_type Type of outcome variable. Can be either "C" for continuous or "D" for dichotomous. The default is "C".
#' @param resampling Resampling indicator. The default is FALSE, do not conduct resampling-based moment matching when the sample size is large, especially for UK biobank-scale data.
#' @param B Number of resampling replicates. The default is 1000, only run resampling replicates when the resampling indicator is TRUE. A larger value leads to more accurate and stable p-value calculation, but requires more computing time.
#' @return It returns a list used for function GeneScan1D(), GeneScan3D() and GeneScan3D.KnockoffGeneration().
#' @examples
#' library(GeneScan3DKnock)
#'
#'# Load data example
#'data("GeneScan3D.example")
#'# Y: outcomes, n by 1 matrix for n=2000 individuals
#'# X: covariates, n by d matrix for d=1 covariate
#'Y=GeneScan3D.example$Y; X=GeneScan3D.example$X;
#'
#'# Preliminary data management
#'result.null.model=GeneScan.Null.Model(Y, X, out_type="C", resampling=FALSE)
#'
#' @export
GeneScan.Null.Model<-function(Y, X=NULL, id=NULL, out_type="C", resampling=FALSE,B=1000){
   
   Y<-as.matrix(Y);n<-nrow(Y)
   
   if(length(X)!=0){X0<-svd(as.matrix(X))$u}else{X0<-NULL}
   X0<-cbind(rep(1,n),X0)
   
   if(out_type=="C"){nullglm<-glm(Y~0+X0,family=gaussian)}
   if(out_type=="D"){nullglm<-glm(Y~0+X0,family=binomial)}
   
   if (length(id)==0){id<-1:n}
   
   mu<-nullglm$fitted.values;Y.res<-Y-mu;
   #permute the residuals for B times when sample size is small
   re.Y.res=NULL
   if(resampling==TRUE){
      index<-sapply(1:B,function(x)sample(1:length(Y)));temp.Y.res<-Y.res[as.vector(index)]
      re.Y.res<-matrix(temp.Y.res,length(Y),B)
   }
   
   #prepare invserse matrix for covariates
   if(out_type=='D'){v<-mu*(1-mu)}else{v<-rep(as.numeric(var(Y.res)),length(Y))}
   inv.X0<-solve(t(X0)%*%(v*X0))
   
   #prepare the preliminary features
   result.null.model<-list(Y=Y,id=id,n=n,mu=mu,res=Y.res,v=v,
                           X0=X0,nullglm=nullglm,out_type=out_type,
                           re.Y.res=re.Y.res,inv.X0=inv.X0)
   return(result.null.model)
}

#' Conduct GeneScan1D analysis on the gene buffer region.
#'
#' This function conducts gene-based scan test on the gene buffer region using 1D windows with sizes 1-5-10 Kb. For binary traits, we conduct SPA gene-based tests to deal with imbalance case-control issues.
#'
#' @param G The genotype matrix in the gene buffer region, which is a n*p matrix where n is the number of individuals and p is the number of genetic variants in the gene buffer region.
#' @param Z A p*q functional annotation matrix where p is the number of genetic variants in the gene buffer region and q is the number of functional annotations. If Z is NULL (do not incorporate any functional annotations), the minor allele frequency weighted dispersion and/or burden tests are applied. Specifically, Beta(MAF; 1; 25) weights are used for rare variants and weights one are used for common variants.
#' @param window.size The 1-D window sizes in base pairs to scan the gene buffer region. The recommended window sizes are c(1000,5000,10000).
#' @param pos  The positions of genetic variants in the gene buffer region, an p dimensional vector. Each position corresponds to a column in the genotype matrix and a row in the functional annotation matrix.
#' @param MAC.threshold Threshold for minor allele count. Variants below MAC.threshold are ultra-rare variants. The recommended level is 10.
#' @param MAF.threshold Threshold for minor allele frequency. Variants below MAF.threshold are rare variants. The recommended level is 0.01.
#' @param Gsub.id The subject id corresponding to the genotype matrix, an n dimensional vector. The default is NULL, where the matched phenotype and genotype matrices are assumed.
#' @param resampling Resampling indicator. The default is FALSE, do not conduct resampling-based moment matching when the sample size is large, especially for UK biobank-scale data.
#' @param result.null.model The output of function "GeneScan.Null.Model()".
#' @return \item{GeneScan1D.Cauchy.pvalue}{Cauchy combination p-values of all, common and rare variants for GeneScan1D analysis.}
#' @return \item{M}{Number of 1D scanning windows.}
#' @examples
#' library(GeneScan3DKnock)
#'
#'# Load data example
#'# Y: outcomes, n by 1 matrix for n=2000 individuals
#'# X: covariates, n by d matrix for d=1 covariate
#'# G_gene_buffer: genotype matrix of gene buffer region, n by p matrix, p=287 variants
#'# Z_gene_buffer: p by q functional annotation matrix, q=1 functional annotation
#'# pos_gene_buffer: positions of p=287 genetic variants
#'
#'data("GeneScan3D.example")
#'Y=GeneScan3D.example$Y; X=GeneScan3D.example$X;
#'G_gene_buffer=GeneScan3D.example$G_gene_buffer;Z_gene_buffer=GeneScan3D.example$Z_gene_buffer;
#'pos_gene_buffer=GeneScan3D.example$pos_gene_buffer;
#'
#'# Preliminary data management
#'set.seed(12345)
#'result.null.model=GeneScan.Null.Model(Y, X, out_type="C", resampling=FALSE)
#'
#'#Conduct GeneScan1D analysis
#'result.GeneScan1D=GeneScan1D(G=G_gene_buffer,Z=Z_gene_buffer,pos=pos_gene_buffer,
#'                             window.size=c(1000,5000,10000),
#'                             MAC.threshold=10,MAF.threshold=0.01,Gsub.id=NULL,resampling=FALSE,
#'                             result.null.model=result.null.model)
#'result.GeneScan1D$GeneScan1D.Cauchy.pvalue 
#'result.GeneScan1D$M
#'
#' @import SKAT
#' @import Matrix
#' @import WGScan
#' @import SPAtest
#' @import CompQuadForm
#' @export
GeneScan1D<-function(G=G_gene_buffer,Z=NULL,window.size=c(1000,5000,10000), pos=pos_gene_buffer,
                     MAC.threshold=10,MAF.threshold=0.01,Gsub.id=NULL,resampling=FALSE,result.null.model=result.null.model){
   
   #load preliminary features
   mu<-result.null.model$nullglm$fitted.values;
   Y.res<-result.null.model$Y-mu
   re.Y.res<-result.null.model$re.Y.res 
   X0<-result.null.model$X0
   outcome<-result.null.model$out_type
   
   impute.method='fixed'
   #match phenotype id and genotype id
   if(length(Gsub.id)==0){match.index<-match(result.null.model$id,1:nrow(G))}else{
      match.index<-match(result.null.model$id,Gsub.id)
   }
   if(mean(is.na(match.index))>0){
      msg<-sprintf("Some individuals are not matched with genotype. The rate is%f", mean(is.na(match.index)))
      warning(msg,call.=F)
   }
   #individuals ids are matched with genotype
   G=Matrix(G[match.index,])
   if(ncol(G)==0|ncol(G)==1){
      stop('Number of variants in the gene buffer region is 0 or 1')
   }
   
   #missing genotype imputation
   G[G==-9 | G==9]=NA
   N_MISS=sum(is.na(G))
   MISS.freq=apply(is.na(G),2,mean)
   if(N_MISS>0){
      msg<-sprintf("The missing genotype rate is %f. Imputation is applied.", N_MISS/nrow(G)/ncol(G))
      warning(msg,call.=F)
      G=Impute(G,impute.method)
   }
   
   #MAF filtering
   MAF<-apply(G,2,mean)/2 #MAF of nonfiltered variants
   G[,MAF>0.5 & !is.na(MAF)]<-2-G[,MAF>0.5 & !is.na(MAF)]
   MAF<-apply(G,2,mean)/2
   s<-apply(G,2,sd)
   SNP.index<-which(MAF>0 & s!=0 & !is.na(MAF)) 
   
   check.index<-which(MAF>0 & s!=0 & !is.na(MAF)  & MISS.freq<0.1)
   if(length(check.index)<=1 ){
      stop('Number of variants with missing rate <=10% in the gene plus buffer region is <=1')
   }
   
   G<-Matrix(G[,SNP.index])
   if(!is.null(Z)){Z<-Matrix(Z[SNP.index,])}
   
   pos=pos[SNP.index]
   
   #generate window matrix to specify the variants in each window
   window.matrix0_gene_buffer<-c()
   for(size in window.size){
      if (size==1){next}
      pos.tag<-seq(min(pos),max(pos),by=size*1/2)
      pos.tag<-sapply(pos.tag,function(x)pos[which.min(abs(x-pos))])
      window.matrix0_gene_buffer<-cbind(window.matrix0_gene_buffer,sapply(pos.tag,function(x)as.numeric(pos>=x & pos<x+size)))
   }
   
   window.string_gene_buffer<-apply(window.matrix0_gene_buffer,2,function(x)paste(as.character(x),collapse = ""))
   window.matrix_gene_buffer<-Matrix(window.matrix0_gene_buffer[,match(unique(window.string_gene_buffer),window.string_gene_buffer)])
   #Number of 1-D windows to scan the gene buffer region
   M_gene_buffer=dim(window.matrix_gene_buffer)[2]
   
   ##single variant score tests, using fastSPA in ScoreTest_SPA function for binary traits 
   p.single<-Get.p(G,result.null.model) 
   length(p.single) 
   #score statistics
   S=t(G)%*%Y.res
   
   GeneScan1D.Cauchy.window=matrix(NA,nrow=M_gene_buffer,ncol=3)
   for (m in 1:M_gene_buffer){
      #print(paste0('1D-window',m))
      
      #Create index for each window
      index.window<-(window.matrix_gene_buffer[,m]==1)
      G.window=G[,index.window]
      G.window=Matrix(G.window)
      if(!is.null(Z)){
         Z.window=Z[index.window,]
         Z.window=Matrix(Z.window)
      }else{
         Z.window=NULL
      }
      #if there is only 1 variant in this window, then do not conduct combined test in this window, move to the next one
      if(dim(G.window)[2]==1){
         next
      }
      
      MAF.window<-apply(G.window,2,mean)/2
      MAC.window<-apply(G.window,2,sum)   
      weight.beta_125<-dbeta(MAF.window,1,25) 
      weight.beta_1<-dbeta(MAF.window,1,1) 
      
      weight.matrix<-cbind(MAC.window<MAC.threshold,(MAF.window<MAF.threshold&MAC.window>=MAC.threshold)*weight.beta_125,(MAF.window>=MAF.threshold)*weight.beta_1) 
      #ultra-rare variants, rare and common variants
      colnames(weight.matrix)<-c('MAC<MAC.threshold','MAF<MAF.threshold&MAC>=MAC.threshold&Beta','MAF>=MAF.thresholdBeta')
      weight.matrix<-Matrix(weight.matrix)
      
      #rare variants
      if (!is.null(Z.window)){
         colnames(Z.window)<-paste0('MAF<MAF.threshold&MAC>=MAC.threshold&',1:ncol(Z.window))
         weight.matrix<-cbind(weight.matrix,(MAF.window<MAF.threshold&MAC.window>=MAC.threshold)*Z.window)
         weight.matrix<-Matrix(weight.matrix)
      } 
      
      #Single variant score test for all variants in the window, SPA p-values for binary traits
      p.single.window<-p.single[index.window]
      
      if(outcome=='D'){v<-result.null.model$v}else{v<-rep(as.numeric(var(Y.res)),nrow(G.window))}
      A<-t(G.window)%*%(v*G.window)
      B<-t(G.window)%*%(v*X0)
      C<-solve(t(X0)%*%(v*X0))
      K<-A-B%*%C%*%t(B) #covariance matrix
      
      #apply SPA gene-based tests for binary trait, deal with imbalance case-control
      if(outcome=='D'){ 
         V=diag(K)
         #adjusted variance
         v_tilde=as.vector(S^2)[index.window]/qchisq(p.single.window,df = 1, ncp = 0, lower.tail = FALSE,log.p = FALSE)
         #adjusted covariance matrix
         K_tilde=diag(sqrt(v_tilde/V))%*%K%*%diag(sqrt(v_tilde/V))
      }
      
      ##Burden test
      #for continuous traits, compute p-value of Q_Burden/Scale from chi-square 1 analytically
      #for binary traits, calculate the SPA gene-based p-value of Burden
      p.burden<-matrix(NA,1,ncol(weight.matrix))
      if(resampling==TRUE){
         for (j in 1:ncol(weight.matrix)){
            temp.window.matrix<-weight.matrix[,j]
            X<-as.matrix(G.window%*%temp.window.matrix)
            p.burden[,j]<-Get.p.base(X,result.null.model)
         }
      }else{ #do not conduct resampling-based moment matching for large sample size
         for (j in 1:ncol(weight.matrix)){
            if (sum(weight.matrix[,j]!=0)>1){ 
               #only conduct Burden test for at least 1 variants
               temp.window.matrix<-weight.matrix[,j]
               X<-as.matrix(G.window%*%temp.window.matrix)
               weights=as.vector(weight.matrix[,j])
               if(outcome=='D'){ #SPA-adjusted
                  p.burden[,j]<-pchisq(as.numeric((t(X)%*%Y.res)^2/weights%*%K_tilde%*%t(t(weights))),df=1,lower.tail=F) 
               }else{ 
                  #continuous
                  p.burden[,j]<-pchisq(as.numeric((t(X)%*%Y.res)^2/weights%*%K%*%t(t(weights))),df=1,lower.tail=F) 
               }
            }
         } 
      }
      #SKAT test
      p.dispersion<-matrix(NA,1,ncol(weight.matrix))
      score<-as.vector(S)[index.window]
      if(resampling==TRUE){
         re.score<-t(t(G.window)%*%re.Y.res) #resampling for 1000 times
         for (j in 1:ncol(weight.matrix)){
            #For extremely rare variants, do not conduct SKAT
            p.dispersion[,j]<-Get.p.SKAT(score,re.score,K,window.matrix=as.matrix(rep(1,sum(index.window))),weight=(MAC.window>=MAC.threshold)*weight.matrix[,j]) 
         }  
      }else{
         #For extremely rare variants, do not conduct SKAT, change MAC.threshold to 10, do not apply resampling-based moment matching
         weight.matrix0=(MAC.window>=MAC.threshold)*weight.matrix
         for (j in 1:ncol(weight.matrix)){
            if (sum(weight.matrix[,j]!=0)>1){ #only conduct SKAT test for at least 1 variants
               if(outcome=='D'){ 
                  #binary
                  p.dispersion[,j]<-Get.p.SKAT_noMA(score,K=K_tilde,window.matrix=as.matrix(rep(1,sum(index.window))),weight=(MAC.window>=MAC.threshold)*weight.matrix[,j])
               }else{ 
                  #continuous
                  p.dispersion[,j]<-Get.p.SKAT_noMA(score,K=K,window.matrix=as.matrix(rep(1,sum(index.window))),weight=(MAC.window>=MAC.threshold)*weight.matrix[,j]) 
               }
            }
         }
      }
      
      p.individual1<-Get.cauchy.scan(p.single.window,as.matrix((MAC.window>=MAC.threshold & MAF.window<MAF.threshold))) #rare variants
      p.individual2<-Get.cauchy.scan(p.single.window,as.matrix((MAF.window>=MAF.threshold))) #common and low frequency variants
      p.individual<-cbind(p.burden,p.dispersion,p.individual1,p.individual2);
      colnames(p.individual)<-c(paste0('burden_',colnames(weight.matrix)),paste0('dispersion_',colnames(weight.matrix)),'singleCauchy_MAF<MAF.threshold&MAC>=MAC.threshold','singleCauchy_MAF>=MAF.threshold')
      
      p.Cauchy<-as.matrix(apply(p.individual,1,Get.cauchy)) 
      #aggregated Cauchy association test
      test.common<-grep('MAF>=MAF.threshold',colnames(p.individual))
      p.Cauchy.common<-as.matrix(apply(p.individual[,test.common,drop=FALSE],1,Get.cauchy))
      p.Cauchy.rare<-as.matrix(apply(p.individual[,-test.common,drop=FALSE],1,Get.cauchy))
      GeneScan1D.Cauchy.window[m,]=c(p.Cauchy,p.Cauchy.common,p.Cauchy.rare)
   }
   GeneScan1D.Cauchy=c(Get.cauchy(GeneScan1D.Cauchy.window[,1]),Get.cauchy(GeneScan1D.Cauchy.window[,2]),Get.cauchy(GeneScan1D.Cauchy.window[,3]))
   return(list(GeneScan1D.Cauchy.pvalue=GeneScan1D.Cauchy,M=M_gene_buffer))
}



#' Conduct GeneScan3D analysis on the gene buffer region, integrating promoter and R enhancers. 
#'
#' This function conducts gene-based scan test on the gene buffer region, integrating proximal and distal regulatory elements for a gene, i.e., promoter and R enhancers. For binary traits, we conduct SPA gene-based tests to deal with imbalance case-control issues.
#'
#' @param G The genotype matrix in the gene buffer region, which is a n*p matrix where n is the number of individuals and p is the number of genetic variants in the gene buffer region.
#' @param Z A p*q functional annotation matrix, where p is the number of genetic variants in the gene buffer region and q is the number of functional annotations. If Z is NULL (do not incorporate any functional annotations), the minor allele frequency weighted dispersion and/or burden tests are applied. Specifically, Beta(MAF; 1; 25) weights are used for rare variants and weights one are used for common variants.
#' @param G.promoter The genotype matrix for promoter, which can be NULL, that is, do not integrate promoter.
#' @param Z.promoter The functional annotation matrix for promoter. Z.promoter can be NULL.
#' @param G.EnhancerAll The genotype matrix for R enhancers, by combining the genotype matrix of each enhancer by columns.
#' @param Z.EnhancerAll The functional annotation matrix for R enhancers, by combining the functional annotation matrix of each enhancer by rows. Z.EnhancerAll can be NULL.
#' @param R Number of enhancers.
#' @param p_Enhancer Number of variants in R enhancers, which is a 1*R vector.
#' @param window.size The 1-D window sizes in base pairs to scan the gene buffer region. The recommended window sizes are c(1000,5000,10000).
#' @param pos  The positions of genetic variants in the gene buffer region, an p dimensional vector. Each position corresponds to a column in the genotype matrix G and a row in the functional annotation matrix Z.
#' @param MAC.threshold Threshold for minor allele count. Variants below MAC.threshold are ultra-rare variants. The recommended level is 5.
#' @param MAF.threshold Threshold for minor allele frequency. Variants below MAF.threshold are rare variants. The recommended level is 0.01.
#' @param Gsub.id The subject id corresponding to the genotype matrix, an n dimensional vector. The default is NULL, where the matched phenotype and genotype matrices are assumed.
#' @param resampling Resampling indicator. The default is FALSE, do not conduct resampling-based moment matching when the sample size is large, especially for UK biobank-scale data.
#' @param result.null.model The output of function "GeneScan.Null.Model()".
#' @return \item{GeneScan3D.Cauchy.pvalue}{Cauchy combination p-values of all, common and rare variants for GeneScan3D analysis.}
#' @return \item{M}{Number of 1D scanning windows.}
#' @return \item{minp}{Minimum p-values of all, common and rare variants for 3D windows.}
#' @return \item{RE_minp}{The regulartory elements in the 3D windows corresponding to the minimum p-values, for all, common and rare variants. 0 represents promoter and a number from 1 to R represents promoter and r-th enhancer.}
#' @examples
#' library(GeneScan3DKnock)
#'
#'# Load data example
#'# Y: outcomes, n by 1 matrix for n=2000 individuals
#'# X: covariates, n by d matrix for d=1 covariate
#'# G_gene_buffer: genotype matrix of gene buffer region, n by p matrix, p=287 variants
#'# pos_gene_buffer: positions of p=287 genetic variants
#'# Z_gene_buffer: p by q functional annotation matrix, q=1 functional annotation
#'# G_promoter: 2000 by 6 genotype matrix of promoter
#'# Z_promoter: 6 by 1 functional annotation matrix of promoter
#'# G_EnhancerAll: 2000 by 86 genotype matrix of R=2 enhancers; 
#'# Z_EnhancerAll: 86 by 1 functional annotation matrix of R=2 enhancers
#'# p_EnhancerAll: Number of variants for R=2 enhancers.
#'
#' data("GeneScan3D.example")
#'
#'Y=GeneScan3D.example$Y; X=GeneScan3D.example$X; n=length(Y)
#'
#'G_gene_buffer=GeneScan3D.example$G_gene_buffer
#'G_promoter=GeneScan3D.example$G_promoter
#'G_EnhancerAll=cbind(GeneScan3D.example$G_Enhancer1,GeneScan3D.example$G_Enhancer2)
#'
#'Z_gene_buffer=GeneScan3D.example$Z_gene_buffer
#'Z_promoter=GeneScan3D.example$Z_promoter
#'Z_EnhancerAll=rbind(GeneScan3D.example$Z_Enhancer1,GeneScan3D.example$Z_Enhancer2)
#'
#'pos_gene_buffer=GeneScan3D.example$pos_gene_buffer
#'p_EnhancerAll=c(dim(GeneScan3D.example$G_Enhancer1)[2],dim(GeneScan3D.example$G_Enhancer2)[2])
#'
#'# Preliminary data management
#'set.seed(12345)
#'result.null.model=GeneScan.Null.Model(Y, X, out_type="C", resampling=FALSE)
#'
#'# Conduct GeneScan3D analysis
#'result.GeneScan3D=GeneScan3D(G=G_gene_buffer,Z=Z_gene_buffer,
#'                             G.promoter=G_promoter,Z.promoter=Z_promoter,
#'                             G.EnhancerAll=G_EnhancerAll,Z.EnhancerAll=Z_EnhancerAll, 
#'                             R=2,p_Enhancer=p_EnhancerAll,
#'                             pos=pos_gene_buffer,
#'                             window.size=c(1000,5000,10000),MAC.threshold=10,MAF.threshold=0.01,
#'                             result.null.model=result.null.model)
#'result.GeneScan3D$GeneScan3D.Cauchy.pvalue
#'
#' @import SKAT
#' @import Matrix
#' @import WGScan
#' @import SPAtest
#' @import CompQuadForm
#' @export
GeneScan3D<-function(G=G_gene_buffer,Z=Z_gene_buffer,G.promoter=G_promoter,Z.promoter=Z_promoter,G.EnhancerAll=G_EnhancerAll,Z.EnhancerAll=Z_EnhancerAll, R=length(p_EnhancerAll),
                     p_Enhancer=p_EnhancerAll,window.size=c(1000,5000,10000),pos=pos_gene_buffer,
                     MAC.threshold=10,MAF.threshold=0.01,Gsub.id=NULL,resampling=FALSE,result.null.model=result.null.model){
   
   #load preliminary features
   mu<-result.null.model$nullglm$fitted.values;
   Y.res<-result.null.model$Y-mu
   re.Y.res<-result.null.model$re.Y.res 
   X0<-result.null.model$X0
   outcome<-result.null.model$out_type
   
   impute.method='fixed'
   #match phenotype id and genotype id
   if(length(Gsub.id)==0){match.index<-match(result.null.model$id,1:nrow(G))}else{
      match.index<-match(result.null.model$id,Gsub.id)
   }
   if(mean(is.na(match.index))>0){
      msg<-sprintf("Some individuals are not matched with genotype. The rate is%f", mean(is.na(match.index)))
      warning(msg,call.=F)
   }
   #individuals ids are matched with genotype
   G=Matrix(G[match.index,])
   if(ncol(G)==0|ncol(G)==1){
      stop('Number of variants in the gene buffer region is 0 or 1')
   }
   
   #missing genotype imputation
   G[G==-9 | G==9]=NA
   N_MISS=sum(is.na(G))
   MISS.freq=apply(is.na(G),2,mean)
   if(N_MISS>0){
      msg<-sprintf("The missing genotype rate is %f. Imputation is applied.", N_MISS/nrow(G)/ncol(G))
      warning(msg,call.=F)
      G=Impute(G,impute.method)
   }
   
   #MAF filtering
   MAF<-apply(G,2,mean)/2 #MAF of nonfiltered variants
   G[,MAF>0.5 & !is.na(MAF)]<-2-G[,MAF>0.5 & !is.na(MAF)]
   MAF<-apply(G,2,mean)/2
   s<-apply(G,2,sd)
   SNP.index<-which(MAF>0 & s!=0 & !is.na(MAF)) 
   
   check.index<-which(MAF>0 & s!=0 & !is.na(MAF)  & MISS.freq<0.1)
   if(length(check.index)<=1 ){
      stop('Number of variants with missing rate <=10% in the gene plus buffer region is <=1')
   }
   
   G<-Matrix(G[,SNP.index])
   if(!is.null(Z)){Z<-Matrix(Z[SNP.index,])}
   pos=pos[SNP.index]
   
   #generate window matrix to specify the variants in each window
   window.matrix0_gene_buffer<-c()
   for(size in window.size){
      if (size==1){next}
      pos.tag<-seq(min(pos),max(pos),by=size*1/2)
      pos.tag<-sapply(pos.tag,function(x)pos[which.min(abs(x-pos))])
      window.matrix0_gene_buffer<-cbind(window.matrix0_gene_buffer,sapply(pos.tag,function(x)as.numeric(pos>=x & pos<x+size)))
   }
   
   window.string_gene_buffer<-apply(window.matrix0_gene_buffer,2,function(x)paste(as.character(x),collapse = ""))
   window.matrix_gene_buffer<-Matrix(window.matrix0_gene_buffer[,match(unique(window.string_gene_buffer),window.string_gene_buffer)])
   #Number of 1-D windows to scan the gene buffer region
   M_gene_buffer=dim(window.matrix_gene_buffer)[2]
   
   ##single variant score tests, using fastSPA in ScoreTest_SPA function for binary traits 
   p.single<-Get.p(G,result.null.model) 
   #score statistics
   S=t(G)%*%Y.res
   
   GeneScan1D.Cauchy.window=matrix(NA,nrow=M_gene_buffer,ncol=3)
   #Burden test: for continuous traits, compute p-value of Q_Burden/Scale from chi-square 1 analytically; for binary traits, use SPA gene- or region-based score test
   #SKAT test: for continuous traits, compute p-value use Davies and if Davies fail to converge, we use resampling moment-based adjustment (MA); for binary traits, use SPA gene- or region-based score test
   for (m in 1:M_gene_buffer){
      #print(paste0('1D-window',m))
      
      #Create index for each window
      index.window<-(window.matrix_gene_buffer[,m]==1)
      G.window=G[,index.window]
      G.window=Matrix(G.window)
      if(!is.null(Z)){
         Z.window=Z[index.window,]
         Z.window=Matrix(Z.window)
      }else{
         Z.window=NULL
      }
      #if there is only 1 variant in this window, then do not conduct combined test in this window, move to the next one
      if(dim(G.window)[2]==1){
         next
      }
      
      MAF.window<-apply(G.window,2,mean)/2
      MAC.window<-apply(G.window,2,sum)   
      weight.beta_125<-dbeta(MAF.window,1,25) 
      weight.beta_1<-dbeta(MAF.window,1,1) 
      
      weight.matrix<-cbind(MAC.window<MAC.threshold,(MAF.window<MAF.threshold&MAC.window>=MAC.threshold)*weight.beta_125,(MAF.window>=MAF.threshold)*weight.beta_1) 
      #ultra-rare variants, rare and common variants
      colnames(weight.matrix)<-c('MAC<MAC.threshold','MAF<MAF.threshold&MAC>=MAC.threshold&Beta','MAF>=MAF.thresholdBeta')
      weight.matrix<-Matrix(weight.matrix)
      
      #rare variants
      if (!is.null(Z.window)){
         colnames(Z.window)<-paste0('MAF<MAF.threshold&MAC>=MAC.threshold&',1:ncol(Z.window))
         weight.matrix<-cbind(weight.matrix,(MAF.window<MAF.threshold&MAC.window>=MAC.threshold)*Z.window)
         weight.matrix<-Matrix(weight.matrix)
      } 
      
      #Single variant score test for all variants in the window, SPA p-values for binary traits
      p.single.window<-p.single[index.window]
      
      if(outcome=='D'){v<-result.null.model$v}else{v<-rep(as.numeric(var(Y.res)),nrow(G.window))}
      A<-t(G.window)%*%(v*G.window)
      B<-t(G.window)%*%(v*X0)
      C<-solve(t(X0)%*%(v*X0))
      K<-A-B%*%C%*%t(B) #covariance matrix
      
      #apply SPA gene-based tests for binary trait, deal with imbalance case-control
      if(outcome=='D'){ 
         V=diag(K)
         #adjusted variance
         v_tilde=as.vector(S^2)[index.window]/qchisq(p.single.window,df = 1, ncp = 0, lower.tail = FALSE,log.p = FALSE)
         #adjusted covariance matrix
         K_tilde=diag(sqrt(v_tilde/V))%*%K%*%diag(sqrt(v_tilde/V))
      }
      
      ##Burden test
      #for continuous traits, compute p-value of Q_Burden/Scale from chi-square 1 analytically
      #for binary traits, calculate the SPA gene-based p-value of Burden
      p.burden<-matrix(NA,1,ncol(weight.matrix))
      if(resampling==TRUE){
         for (j in 1:ncol(weight.matrix)){
            temp.window.matrix<-weight.matrix[,j]
            X<-as.matrix(G.window%*%temp.window.matrix)
            p.burden[,j]<-Get.p.base(X,result.null.model)
         }
      }else{ #do not conduct resampling-based moment matching for large sample size
         for (j in 1:ncol(weight.matrix)){
            if (sum(weight.matrix[,j]!=0)>1){ 
               #only conduct Burden test for at least 1 variants
               temp.window.matrix<-weight.matrix[,j]
               X<-as.matrix(G.window%*%temp.window.matrix)
               weights=as.vector(weight.matrix[,j])
               if(outcome=='D'){ #SPA-adjusted
                  p.burden[,j]<-pchisq(as.numeric((t(X)%*%Y.res)^2/weights%*%K_tilde%*%t(t(weights))),df=1,lower.tail=F) 
               }else{ 
                  #continuous
                  p.burden[,j]<-pchisq(as.numeric((t(X)%*%Y.res)^2/weights%*%K%*%t(t(weights))),df=1,lower.tail=F) 
               }
            }
         } 
      }
      
      
      #SKAT test
      p.dispersion<-matrix(NA,1,ncol(weight.matrix))
      score<-as.vector(S)[index.window]
      if(resampling==TRUE){
         re.score<-t(t(G.window)%*%re.Y.res) #resampling for 1000 times
         for (j in 1:ncol(weight.matrix)){
            #For extremely rare variants, do not conduct SKAT
            p.dispersion[,j]<-Get.p.SKAT(score,re.score,K,window.matrix=as.matrix(rep(1,sum(index.window))),weight=(MAC.window>=MAC.threshold)*weight.matrix[,j]) 
         }  
      }else{
         #For extremely rare variants, do not conduct SKAT, change MAC.threshold to 10, do not apply resampling-based moment matching
         weight.matrix0=(MAC.window>=MAC.threshold)*weight.matrix
         for (j in 1:ncol(weight.matrix)){
            if (sum(weight.matrix[,j]!=0)>1){ #only conduct SKAT test for at least 1 variants
               if(outcome=='D'){ 
                  #binary
                  p.dispersion[,j]<-Get.p.SKAT_noMA(score,K=K_tilde,window.matrix=as.matrix(rep(1,sum(index.window))),weight=(MAC.window>=MAC.threshold)*weight.matrix[,j])
               }else{ 
                  #continuous
                  p.dispersion[,j]<-Get.p.SKAT_noMA(score,K=K,window.matrix=as.matrix(rep(1,sum(index.window))),weight=(MAC.window>=MAC.threshold)*weight.matrix[,j]) 
               }
            }
         }
      }
      
      p.individual1<-Get.cauchy.scan(p.single.window,as.matrix((MAC.window>=MAC.threshold & MAF.window<MAF.threshold))) #rare variants
      p.individual2<-Get.cauchy.scan(p.single.window,as.matrix((MAF.window>=MAF.threshold))) #common and low frequency variants
      p.individual<-cbind(p.burden,p.dispersion,p.individual1,p.individual2);
      colnames(p.individual)<-c(paste0('burden_',colnames(weight.matrix)),paste0('dispersion_',colnames(weight.matrix)),'singleCauchy_MAF<MAF.threshold&MAC>=MAC.threshold','singleCauchy_MAF>=MAF.threshold')
      
      p.Cauchy<-as.matrix(apply(p.individual,1,Get.cauchy)) 
      #aggregated Cauchy association test
      test.common<-grep('MAF>=MAF.threshold',colnames(p.individual))
      p.Cauchy.common<-as.matrix(apply(p.individual[,test.common,drop=FALSE],1,Get.cauchy))
      p.Cauchy.rare<-as.matrix(apply(p.individual[,-test.common,drop=FALSE],1,Get.cauchy))
      GeneScan1D.Cauchy.window[m,]=c(p.Cauchy,p.Cauchy.common,p.Cauchy.rare)
   }
   GeneScan1D.Cauchy=c(Get.cauchy(GeneScan1D.Cauchy.window[,1]),Get.cauchy(GeneScan1D.Cauchy.window[,2]),Get.cauchy(GeneScan1D.Cauchy.window[,3]))
   
   ###promoter
   if(is.null(G.promoter)){
      warning('no promoter')
      GeneScan3D.Cauchy.promoter=c()
   }else{
      
      #match phenotype id and genotype id
      G.promoter=Matrix(G.promoter[match.index,])
      
      #missing genotype imputation
      G.promoter[G.promoter==-9 | G.promoter==9]=NA
      N_MISS.promoter=sum(is.na(G.promoter))
      MISS.freq.promoter=apply(is.na(G.promoter),2,mean)
      if(N_MISS.promoter>0){
         msg<-sprintf("The missing genotype rate is %f. Imputation is applied.", N_MISS.promoter/nrow(G.promoter)/ncol(G.promoter))
         warning(msg,call.=F)
         G.promoter=Impute(G.promoter,impute.method)
      }
      
      #MAF filtering
      MAF.promoter<-apply(G.promoter,2,mean)/2 #MAF of nonfiltered variants
      G.promoter[,MAF.promoter>0.5 & !is.na(MAF.promoter)]<-2-G.promoter[,MAF.promoter>0.5 & !is.na(MAF.promoter)]
      MAF.promoter<-apply(G.promoter,2,mean)/2
      s.promoter<-apply(G.promoter,2,sd)
      SNP.index.promoter<-which(MAF.promoter>0 & s.promoter!=0 & !is.na(MAF.promoter)) 
      
      G.promoter<-Matrix(G.promoter[,SNP.index.promoter])
      if(!is.null(Z.promoter)){Z.promoter<-Matrix(Z.promoter[SNP.index.promoter,])}
      
      p_promoter=dim(G.promoter)[2] 
      
      #Obtain p-value for promoter
      if (p_promoter==0){
         warning('0 variant in promoter')
         GeneScan3D.Cauchy.promoter=c()
      }else{
         
         G.window.promoter=Matrix(G.promoter)
         if(!is.null(Z.promoter)){
            Z.window.promoter=Matrix(Z.promoter)
         }else{
            Z.window.promoter=NULL
         }
         
         MAF.window.promoter<-apply(G.window.promoter,2,mean)/2
         MAC.window.promoter<-apply(G.window.promoter,2,sum)
         
         weight.beta_125<-dbeta(MAF.window.promoter,1,25)
         weight.beta_1<-dbeta(MAF.window.promoter,1,1)
         weight.matrix<-cbind(MAC.window.promoter<MAC.threshold,(MAF.window.promoter<MAF.threshold&MAC.window.promoter>=MAC.threshold)*weight.beta_125,(MAF.window.promoter>=MAF.threshold)*weight.beta_1)
         colnames(weight.matrix)<-c('MAC<MAC.threshold','MAF<MAF.threshold&MAC>=MAC.threshold&Beta','MAF>=MAF.thresholdBeta')
         
         ##adding additional functional scores
         if (!is.null(Z.window.promoter)){
            colnames(Z.window.promoter)<-paste0('MAF<MAF.threshold&MAC>=MAC.threshold&FS',1:ncol(Z.window.promoter))
            weight.matrix<-cbind(weight.matrix,(MAF.window.promoter<MAF.threshold&MAC.window.promoter>=MAC.threshold)*Z.window.promoter)
         }
         weight.matrix<-Matrix(weight.matrix)
         
         if(outcome=='D'){v<-result.null.model$v}else{v<-rep(as.numeric(var(Y.res)),nrow(G.window.promoter))}
         A<-t(G.window.promoter)%*%(v*G.window.promoter)
         B<-t(G.window.promoter)%*%(v*X0)
         C<-solve(t(X0)%*%(v*X0))
         K<-A-B%*%C%*%t(B) #covariance matrix
         
         #SPA gene-based tests
         if(outcome=='D'){ 
            V=diag(K)
            #adjusted variance
            v_tilde=as.vector(S^2)[index.window]/qchisq(p.single.window,df = 1, ncp = 0, lower.tail = FALSE,log.p = FALSE)
            #adjusted covariance matrix
            K_tilde=diag(sqrt(v_tilde/V))%*%K%*%diag(sqrt(v_tilde/V))
         }
         
         #Burden test
         p.burden.promoter<-matrix(NA,1,ncol(weight.matrix))
         if(resampling==TRUE){
            for (j in 1:ncol(weight.matrix)){
               temp.window.matrix<-weight.matrix[,j]
               X<-as.matrix(G.window.promoter%*%temp.window.matrix)
               p.burden.promoter[,j]<-Get.p.base(X,result.null.model)
            }
         }else{ #do not conduct resampling-based moment matching for large sample size
            for (j in 1:ncol(weight.matrix)){
               if (sum(weight.matrix[,j]!=0)>1){ 
                  #only conduct Burden test for at least 1 variants
                  temp.window.matrix<-weight.matrix[,j]
                  X<-as.matrix(G.window.promoter%*%temp.window.matrix)
                  weights=as.vector(weight.matrix[,j])
                  if(outcome=='D'){ #SPA-adjusted
                     p.burden.promoter[,j]<-pchisq(as.numeric((t(X)%*%Y.res)^2/weights%*%K_tilde%*%t(t(weights))),df=1,lower.tail=F) 
                  }else{ 
                     #continuous
                     p.burden.promoter[,j]<-pchisq(as.numeric((t(X)%*%Y.res)^2/weights%*%K%*%t(t(weights))),df=1,lower.tail=F) 
                  }
               }
            } 
         }
         
         #SKAT test
         score<-as.vector(t(G.window.promoter)%*%Y.res)
         p.dispersion.promoter<-matrix(NA,1,ncol(weight.matrix))
         if(resampling==TRUE){
            re.score<-t(t(G.window.promoter)%*%re.Y.res) #resampling for 1000 times
            for (j in 1:ncol(weight.matrix)){
               #For extremely rare variants, do not conduct SKAT
               p.dispersion.promoter[,j]<-Get.p.SKAT(score,re.score,K,window.matrix=as.matrix(rep(1,dim(G.window.promoter)[2])),weight=(MAC.window.promoter>=MAC.threshold)*weight.matrix[,j]) 
            }  
         }else{
            #For extremely rare variants, do not conduct SKAT, change MAC.threshold to 10, do not apply resampling-based moment matching
            weight.matrix0=(MAC.window.promoter>=MAC.threshold)*weight.matrix
            for (j in 1:ncol(weight.matrix)){
               if (sum(weight.matrix[,j]!=0)>1){ #only conduct SKAT test for at least 1 variants
                  if(outcome=='D'){ 
                     #binary
                     p.dispersion.promoter[,j]<-Get.p.SKAT_noMA(score,K=K_tilde,window.matrix=as.matrix(rep(1,dim(G.window.promoter)[2])),weight=(MAC.window.promoter>=MAC.threshold)*weight.matrix[,j])
                  }else{ 
                     #continuous
                     p.dispersion.promoter[,j]<-Get.p.SKAT_noMA(score,K=K,window.matrix=as.matrix(rep(1,dim(G.window.promoter)[2])),weight=(MAC.window.promoter>=MAC.threshold)*weight.matrix[,j]) 
                  }
               }
            }
         }
         
         #Single variant score test for all variants in the window
         p.single.promoter<-Get.p(G.window.promoter,result.null.model)
         p.individual1.promoter<-Get.cauchy.scan(p.single.promoter,as.matrix((MAC.window.promoter>=MAC.threshold & MAF.window.promoter<MAF.threshold))) #rare variants
         p.individual2.promoter<-Get.cauchy.scan(p.single.promoter,as.matrix((MAF.window.promoter>=MAF.threshold))) #common and low frequency variants
         p.individual.promoter<-cbind(p.burden.promoter ,p.dispersion.promoter,p.individual1.promoter,p.individual2.promoter);
         colnames(p.individual.promoter)<-c(paste0('burden_',colnames(weight.matrix)),paste0('dispersion_',colnames(weight.matrix)),
                                            'singleCauchy_MAF<MAF.threshold&MAC>=MAC.threshold','singleCauchy_MAF>=MAF.threshold')
         
         #aggregated Cauchy association test
         p.Cauchy.promoter<-as.matrix(apply(p.individual.promoter,1,Get.cauchy))
         test.common<-grep('MAF>=MAF.threshold',colnames(p.individual.promoter))
         p.Cauchy.common.promoter<-as.matrix(apply(p.individual.promoter[,test.common,drop=FALSE],1,Get.cauchy))
         p.Cauchy.rare.promoter<-as.matrix(apply(p.individual.promoter[,-test.common,drop=FALSE],1,Get.cauchy))
         
         GeneScan3D.Cauchy.promoter=c(p.Cauchy.promoter,p.Cauchy.common.promoter,p.Cauchy.rare.promoter)
      }
   }
   
   ###Obtain p-values for R enhancers
   GeneScan3D.Cauchy.EnhancerAll=c()
   #Enhancer_ind=0
   if(R!=0){
      Enhancer_ind=rep(TRUE,R)
      for (r in 1:R){ #Loop for each enhancer
         #print(paste0('Enhancer',r))
         if (r==1){
            G.Enhancer=G.EnhancerAll[,1:cumsum(p_Enhancer)[r]]
         }else{
            G.Enhancer=G.EnhancerAll[,(cumsum(p_Enhancer)[r-1]+1):cumsum(p_Enhancer)[r]]
         }
         
         if(!is.null(Z.EnhancerAll)){
            if (r==1){
               Z.Enhancer=Z.EnhancerAll[1:cumsum(p_Enhancer)[r],]
            }else{
               Z.Enhancer=Z.EnhancerAll[(cumsum(p_Enhancer)[r-1]+1):cumsum(p_Enhancer)[r],]
            }
         }else{
            Z.Enhancer=NULL
         }
         
         #individuals ids are matched with genotype
         G.Enhancer=Matrix(G.Enhancer[match.index,])
         if(!is.null(Z.Enhancer)){Z.Enhancer=Matrix(Z.Enhancer)}
         
         #missing genotype imputation
         G.Enhancer[G.Enhancer==-9 | G.Enhancer==9]=NA
         N_MISS.Enhancer=sum(is.na(G.Enhancer))
         MISS.freq.Enhancer=apply(is.na(G.Enhancer),2,mean)
         if(N_MISS.Enhancer>0){
            msg<-sprintf("The missing genotype rate is %f. Imputation is applied.", N_MISS.Enhancer/nrow(G.Enhancer)/ncol(G.Enhancer))
            warning(msg,call.=F)
            G.Enhancer=Impute(G.Enhancer,impute.method)
         }
         
         #MAF filtering
         MAF.Enhancer<-apply(G.Enhancer,2,mean)/2 #MAF of nonfiltered variants
         G.Enhancer[,MAF.Enhancer>0.5 & !is.na(MAF.Enhancer)]<-2-G.Enhancer[,MAF.Enhancer>0.5 & !is.na(MAF.Enhancer)]
         MAF.Enhancer<-apply(G.Enhancer,2,mean)/2
         s.Enhancer<-apply(G.Enhancer,2,sd)
         SNP.index.Enhancer<-which(MAF.Enhancer>0 & s.Enhancer!=0 & !is.na(MAF.Enhancer)) 
         
         G.Enhancer<-Matrix(G.Enhancer[,SNP.index.Enhancer])
         if(!is.null(Z.Enhancer)){Z.Enhancer<-Matrix(Z.Enhancer[SNP.index.Enhancer,])}
         
         if(dim(G.Enhancer)[2]<1){
            Enhancer_ind[r]=FALSE
            next
         }else{
            
            G.window.Enhancer=Matrix(G.Enhancer)
            if(!is.null(Z.Enhancer)){Z.window.Enhancer=Matrix(Z.Enhancer)}else{Z.window.Enhancer=NULL}
            
            MAF.window.Enhancer<-apply(G.window.Enhancer,2,mean)/2
            MAC.window.Enhancer<-apply(G.window.Enhancer,2,sum)
            
            weight.beta_125<-dbeta(MAF.window.Enhancer,1,25)
            weight.beta_1<-dbeta(MAF.window.Enhancer,1,1)
            weight.matrix<-cbind(MAC.window.Enhancer<MAC.threshold,(MAF.window.Enhancer<MAF.threshold&MAC.window.Enhancer>=MAC.threshold)*weight.beta_125,(MAF.window.Enhancer>=MAF.threshold)*weight.beta_1)
            colnames(weight.matrix)<-c('MAC<MAC.threshold','MAF<MAF.threshold&MAC>=MAC.threshold&Beta','MAF>=MAF.thresholdBeta')
            
            ##adding additional functional scores
            if (!is.null(Z.window.Enhancer)){
               colnames(Z.window.Enhancer)<-paste0('MAF<MAF.threshold&MAC>=MAC.threshold&FS',1:ncol(Z.window.Enhancer))
               weight.matrix<-cbind(weight.matrix,(MAF.window.Enhancer<MAF.threshold&MAC.window.Enhancer>=MAC.threshold)*Z.window.Enhancer)
            }
            weight.matrix<-Matrix(weight.matrix)
            
            if(outcome=='D'){v<-result.null.model$v}else{v<-rep(as.numeric(var(Y.res)),nrow(G.window.Enhancer))}
            A<-t(G.window.Enhancer)%*%(v*G.window.Enhancer)
            B<-t(G.window.Enhancer)%*%(v*X0)
            C<-solve(t(X0)%*%(v*X0))
            K<-A-B%*%C%*%t(B) #covariance matrix
            
            #SPA gene-based tests
            if(outcome=='D'){ 
               V=diag(K)
               #adjusted variance
               v_tilde=as.vector(S^2)[index.window]/qchisq(p.single.window,df = 1, ncp = 0, lower.tail = FALSE,log.p = FALSE)
               #adjusted covariance matrix
               K_tilde=diag(sqrt(v_tilde/V))%*%K%*%diag(sqrt(v_tilde/V))
            }
            
            ##Burden test
            p.burden.Enhancer<-matrix(NA,1,ncol(weight.matrix))
            if(resampling==TRUE){
               for (j in 1:ncol(weight.matrix)){
                  temp.window.matrix<-weight.matrix[,j]
                  X<-as.matrix(G.window.Enhancer%*%temp.window.matrix)
                  p.burden.Enhancer[,j]<-Get.p.base(X,result.null.model)
               }
            }else{ #do not conduct resampling-based moment matching for large sample size
               for (j in 1:ncol(weight.matrix)){
                  if (sum(weight.matrix[,j]!=0)>1){ 
                     #only conduct Burden test for at least 1 variants
                     temp.window.matrix<-weight.matrix[,j]
                     X<-as.matrix(G.window.Enhancer%*%temp.window.matrix)
                     weights=as.vector(weight.matrix[,j])
                     if(outcome=='D'){ #SPA-adjusted
                        p.burden.Enhancer[,j]<-pchisq(as.numeric((t(X)%*%Y.res)^2/weights%*%K_tilde%*%t(t(weights))),df=1,lower.tail=F) 
                     }else{ 
                        #continuous
                        p.burden.Enhancer[,j]<-pchisq(as.numeric((t(X)%*%Y.res)^2/weights%*%K%*%t(t(weights))),df=1,lower.tail=F) 
                     }
                  }
               } 
            }
            
            #SKAT test
            score<-as.vector(t(G.window.Enhancer)%*%Y.res)
            p.dispersion.Enhancer<-matrix(NA,1,ncol(weight.matrix))
            if(resampling==TRUE){
               re.score<-t(t(G.window.Enhancer)%*%re.Y.res) #resampling for 1000 times
               for (j in 1:ncol(weight.matrix)){
                  #For extremely rare variants, do not conduct SKAT
                  p.dispersion.Enhancer[,j]<-Get.p.SKAT(score,re.score,K,window.matrix=as.matrix(rep(1,dim(G.window.Enhancer)[2])),weight=(MAC.window.Enhancer>=MAC.threshold)*weight.matrix[,j]) 
               }  
            }else{
               #For extremely rare variants, do not conduct SKAT, change MAC.threshold to 10, do not apply resampling-based moment matching
               weight.matrix0=(MAC.window.Enhancer>=MAC.threshold)*weight.matrix
               for (j in 1:ncol(weight.matrix)){
                  if (sum(weight.matrix[,j]!=0)>1){ #only conduct SKAT test for at least 1 variants
                     if(outcome=='D'){ 
                        #binary
                        p.dispersion.Enhancer[,j]<-Get.p.SKAT_noMA(score,K=K_tilde,window.matrix=as.matrix(rep(1,dim(G.window.Enhancer)[2])),weight=(MAC.window.Enhancer>=MAC.threshold)*weight.matrix[,j])
                     }else{ 
                        #continuous
                        p.dispersion.Enhancer[,j]<-Get.p.SKAT_noMA(score,K=K,window.matrix=as.matrix(rep(1,dim(G.window.Enhancer)[2])),weight=(MAC.window.Enhancer>=MAC.threshold)*weight.matrix[,j]) 
                     }
                  }
               }
            }
            
            #Single variant score test for all variants in the window
            p.single.Enhancer<-Get.p(G.window.Enhancer,result.null.model)
            
            p.individual1.Enhancer<-Get.cauchy.scan(p.single.Enhancer,as.matrix((MAC.window.Enhancer>=MAC.threshold & MAF.window.Enhancer<MAF.threshold))) #rare variants
            p.individual2.Enhancer<-Get.cauchy.scan(p.single.Enhancer,as.matrix((MAF.window.Enhancer>=MAF.threshold))) #common and low frequency variants
            p.individual.Enhancer<-cbind(p.burden.Enhancer ,p.dispersion.Enhancer,p.individual1.Enhancer,p.individual2.Enhancer);
            colnames(p.individual.Enhancer)<-c(paste0('burden_',colnames(weight.matrix)),paste0('dispersion_',colnames(weight.matrix)),'singleCauchy_MAF<MAF.threshold&MAC>=MAC.threshold','singleCauchy_MAF>=MAF.threshold')
            
            #aggregated Cauchy association test
            p.Cauchy.Enhancer<-as.matrix(apply(p.individual.Enhancer,1,Get.cauchy))
            test.common<-grep('MAF>=MAF.threshold',colnames(p.individual.Enhancer))
            p.Cauchy.common.Enhancer<-as.matrix(apply(p.individual.Enhancer[,test.common,drop=FALSE],1,Get.cauchy))
            p.Cauchy.rare.Enhancer<-as.matrix(apply(p.individual.Enhancer[,-test.common,drop=FALSE],1,Get.cauchy))
            
            GeneScan3D.Cauchy.Enhancer=c(p.Cauchy.Enhancer,p.Cauchy.common.Enhancer,p.Cauchy.rare.Enhancer)
         }
         GeneScan3D.Cauchy.EnhancerAll=rbind(GeneScan3D.Cauchy.EnhancerAll,GeneScan3D.Cauchy.Enhancer)
      } #end of the loop of R enhancers
   }
   ##Obtain 3D windows and p-values
   #M 1D windows + promoter
   GeneScan3D.window.promoter=data.frame(apply(cbind(GeneScan1D.Cauchy.window[,1],GeneScan3D.Cauchy.promoter[1]),1,Get.cauchy),
                                         apply(cbind(GeneScan1D.Cauchy.window[,2],GeneScan3D.Cauchy.promoter[2]),1,Get.cauchy),
                                         apply(cbind(GeneScan1D.Cauchy.window[,3],GeneScan3D.Cauchy.promoter[3]),1,Get.cauchy))
   colnames(GeneScan3D.window.promoter)=c('all','common','rare')
   
   #M 1D windows + promoter + Enhancer r, r=1, ..., R
   GeneScan3D.window.EnhancerAll=c()
   if(R!=0){
      for (r in 1:dim(GeneScan3D.Cauchy.EnhancerAll)[1]){
         GeneScan3D.window.enhancer=data.frame(apply(cbind(GeneScan1D.Cauchy.window[,1],GeneScan3D.Cauchy.promoter[1],GeneScan3D.Cauchy.EnhancerAll[r,1]),1,Get.cauchy),
                                               apply(cbind(GeneScan1D.Cauchy.window[,2],GeneScan3D.Cauchy.promoter[2],GeneScan3D.Cauchy.EnhancerAll[r,2]),1,Get.cauchy),
                                               apply(cbind(GeneScan1D.Cauchy.window[,3],GeneScan3D.Cauchy.promoter[3],GeneScan3D.Cauchy.EnhancerAll[r,3]),1,Get.cauchy))
         colnames(GeneScan3D.window.enhancer)=c('all','common','rare')
         GeneScan3D.window.EnhancerAll=rbind(GeneScan3D.window.EnhancerAll,GeneScan3D.window.enhancer)
      }
   }
   
   GeneScan3D.Cauchy.RE=rbind(GeneScan3D.window.promoter,GeneScan3D.window.EnhancerAll)
   
   GeneScan3D.Cauchy=c(Get.cauchy(GeneScan3D.Cauchy.RE[,1]), Get.cauchy(GeneScan3D.Cauchy.RE[,2]), Get.cauchy(GeneScan3D.Cauchy.RE[,3]))
   
   ###min-p and RE with min-p
   RE.indicator=c(rep(0,M_gene_buffer),rep((1:R)[Enhancer_ind],each=M_gene_buffer))
   
   if(is.infinite(min(GeneScan3D.Cauchy.RE[,1],na.rm=TRUE))){
      RE_minp.all=NA
   }else{
      RE_minp.all=unique(RE.indicator[which(GeneScan3D.Cauchy.RE[,1]==min(GeneScan3D.Cauchy.RE[,1],na.rm=TRUE))])
   }
   
   if(is.infinite(min(GeneScan3D.Cauchy.RE[,2],na.rm=TRUE))){
      RE_minp.common=NA
   }else{
      RE_minp.common=unique(RE.indicator[which(GeneScan3D.Cauchy.RE[,2]==min(GeneScan3D.Cauchy.RE[,2],na.rm=TRUE))])
   }
   
   if(is.infinite(min(GeneScan3D.Cauchy.RE[,3],na.rm=TRUE))){
      RE_minp.rare=NA
   }else{
      RE_minp.rare=unique(RE.indicator[which(GeneScan3D.Cauchy.RE[,3]==min(GeneScan3D.Cauchy.RE[,3],na.rm=TRUE))])
   }
   
   return(list(GeneScan3D.Cauchy.pvalue=GeneScan3D.Cauchy,M=M_gene_buffer,
               minp=c(min(GeneScan3D.Cauchy.RE[,1],na.rm=TRUE),min(GeneScan3D.Cauchy.RE[,2],na.rm=TRUE),min(GeneScan3D.Cauchy.RE[,3],na.rm=TRUE)),
               RE_minp=c(RE_minp.all,RE_minp.common,RE_minp.rare)))
}

######### Other functions #########
### p-values calculation
Get.p<-function(X,result.null.model){ 
   #single variant score test: for continuous traits, score^2/v follows chi-square 1
   #for binary traits, we use fastSPA in ScoreTest_SPA function
   X<-as.matrix(X)
   mu<-result.null.model$nullglm$fitted.values;Y.res<-result.null.model$Y-mu
   outcome<-result.null.model$out_type
   if(outcome=='D'){
      p<-ScoreTest_SPA(t(X),result.null.model$Y,result.null.model$X,method=c("fastSPA"),minmac=-Inf)$p.value
   }else{
      v<-rep(as.numeric(var(Y.res)),nrow(X))
      p<-pchisq(as.vector((t(X)%*%Y.res)^2)/(apply(X*(v*X),2,sum)-apply(t(X)%*%(v*result.null.model$X0)%*%result.null.model$inv.X0*t(t(result.null.model$X0)%*%as.matrix(v*X)),1,sum)),df=1,lower.tail=F)
   }
   return(as.matrix(p))
}
Get.p.SKAT_noMA<-function(score,K,window.matrix,weight){

   Q<-as.vector(t(score^2)%*%(weight*window.matrix)^2) #SKAT statistics
   K.temp<-weight*t(weight*K)
   
   temp<-K.temp[window.matrix[,1]!=0,window.matrix[,1]!=0]
   if(sum(temp^2)==0){p<-NA}else{
      lambda=eigen(temp,symmetric=T,only.values=T)$values #eigenvalues, mixture of chi-square
      temp.p<-SKAT_davies(Q,lambda,acc=10^(-6))$Qq
      
      if(length(temp.p)==0 || temp.p > 1 || temp.p <= 0){
         temp.p<-Get_Liu_PVal.MOD.Lambda(Q,lambda)
      }
      p<-temp.p
   }
   return(p)
}
Get.p.base<-function(X,result.null.model){
   X<-as.matrix(X)
   mu<-result.null.model$nullglm$fitted.values;Y.res<-result.null.model$Y-mu
   outcome<-result.null.model$out_type
   if(outcome=='D'){v<-mu*(1-mu)}else{v<-rep(as.numeric(var(Y.res)),nrow(X))}
   p<-pchisq((t(X)%*%Y.res)^2/(apply(X*(v*X),2,sum)-apply(t(X)%*%(v*result.null.model$X0)%*%result.null.model$inv.X0*t(t(result.null.model$X0)%*%as.matrix(v*X)),1,sum)),df=1,lower.tail=F)
   p[is.na(p)]<-NA
   return(p)
}
Get.p.moment<-function(Q,re.Q){ #Q a A*q matrix of test statistics, re.Q a B*q matrix of resampled test statistics
   re.mean<-apply(re.Q,2,mean)
   re.variance<-apply(re.Q,2,var)
   re.kurtosis<-apply((t(re.Q)-re.mean)^4,1,mean)/re.variance^2-3
   re.df<-(re.kurtosis>0)*12/re.kurtosis+(re.kurtosis<=0)*100000
   re.p<-t(pchisq((t(Q)-re.mean)*sqrt(2*re.df)/sqrt(re.variance)+re.df,re.df,lower.tail=F))
   return(re.p)
}
Get.p.SKAT<-function(score,re.score,K,window.matrix,weight){
  
   Q<-as.vector(t(score^2)%*%(weight*window.matrix)^2) #SKAT statistics
   K.temp<-weight*t(weight*K)
   
   #fast implementation by resampling based moment matching
   p0<-Get.p.moment(as.vector(t(score^2)%*%(weight*window.matrix)^2),re.score^2%*%(weight*window.matrix)^2)
   p<-p0
   for(i in which(p0<0.01 |p0>=1)){
      
      temp<-K.temp[window.matrix[,i]!=0,window.matrix[,i]!=0]
      if(sum(temp^2)==0){p[i]<-NA;next}
      
      lambda=eigen(temp,symmetric=T,only.values=T)$values
      temp.p<-SKAT_davies(Q[i],lambda,acc=10^(-6))$Qq
      
      if(length(temp.p)==0 || temp.p > 1 || temp.p <= 0){
         temp.p<-Get_Liu_PVal.MOD.Lambda(Q[i],lambda)
      }
      p[i]<-temp.p
   }
   return(as.matrix(p))
}
SKAT_davies <- function(q,lambda,h = rep(1,length(lambda)),delta = rep(0,length(lambda)),sigma=0,lim=10000,acc=0.0001) {
   r <- length(lambda)
   if (length(h) != r) warning("lambda and h should have the same length!")
   if (length(delta) != r) warning("lambda and delta should have the same length!")
   #out <- .C("qfc",lambdas=as.double(lambda),noncentral=as.double(delta),df=as.integer(h),r=as.integer(r),sigma=as.double(sigma),q=as.double(q),lim=as.integer(lim),acc=as.double(acc),trace=as.double(rep(0,7)),ifault=as.integer(0),res=as.double(0),PACKAGE="SKAT")
   out=davies(q, lambda, h = rep(1, length(lambda)), delta = rep(0,length(lambda)), sigma = 0, lim = 10000, acc = 0.0001)
   out$res <- 1 - out$res
   return(list(trace=out$trace,ifault=out$ifault,Qq=out$res))
}
Get_Liu_PVal.MOD.Lambda<-function(Q.all, lambda, log.p=FALSE){
   param<-Get_Liu_Params_Mod_Lambda(lambda)
   Q.Norm<-(Q.all - param$muQ)/param$sigmaQ
   Q.Norm1<-Q.Norm * param$sigmaX + param$muX
   p.value<- pchisq(Q.Norm1,  df = param$l,ncp=param$d, lower.tail=FALSE, log.p=log.p)
   return(p.value)
}
Get_Liu_Params_Mod_Lambda<-function(lambda){
   ## Helper function for getting the parameters for the null approximation
   
   c1<-rep(0,4)
   for(i in 1:4){
      c1[i]<-sum(lambda^i)
   }
   
   muQ<-c1[1]
   sigmaQ<-sqrt(2 *c1[2])
   s1 = c1[3] / c1[2]^(3/2)
   s2 = c1[4] / c1[2]^2
   
   beta1<-sqrt(8)*s1
   beta2<-12*s2
   type1<-0
   
   #print(c(s1^2,s2))
   if(s1^2 > s2){
      a = 1/(s1 - sqrt(s1^2 - s2))
      d = s1 *a^3 - a^2
      l = a^2 - 2*d
   } else {
      type1<-1
      l = 1/s2
      a = sqrt(l)
      d = 0
   }
   muX <-l+d
   sigmaX<-sqrt(2) *a
   
   re<-list(l=l,d=d,muQ=muQ,muX=muX,sigmaQ=sigmaQ,sigmaX=sigmaX)
   return(re)
}
##Cauchy
Get.cauchy.scan<-function(p,window.matrix){
   p[p>0.99]<-0.99
   is.small<-(p<1e-16)
   temp<-rep(0,length(p))
   temp[is.small]<-1/p[is.small]/pi
   temp[!is.small]<-as.numeric(tan((0.5-p[!is.small])*pi))
   
   cct.stat<-as.numeric(t(temp)%*%window.matrix/apply(window.matrix,2,sum))
   is.large<-cct.stat>1e+15 & !is.na(cct.stat)
   is.regular<-cct.stat<=1e+15 & !is.na(cct.stat)
   pval<-rep(NA,length(cct.stat))
   pval[is.large]<-(1/cct.stat[is.large])/pi
   pval[is.regular]<-1-pcauchy(cct.stat[is.regular])
   return(pval)
}
Get.cauchy<-function(p){
   p[p>0.99]<-0.99
   is.small<-(p<1e-16) & !is.na(p)
   is.regular<-(p>=1e-16) & !is.na(p)
   temp<-rep(NA,length(p))
   temp[is.small]<-1/p[is.small]/pi
   temp[is.regular]<-as.numeric(tan((0.5-p[is.regular])*pi))
   
   cct.stat<-mean(temp,na.rm=T)
   if(is.na(cct.stat)){return(NA)}
   if(cct.stat>1e+15){return((1/cct.stat)/pi)}else{
      return(1-pcauchy(cct.stat))
   }
}
Impute<-function(Z, impute.method){
   p<-dim(Z)[2]
   if(impute.method =="random"){
      for(i in 1:p){
         IDX<-which(is.na(Z[,i]))
         if(length(IDX) > 0){
            maf1<-mean(Z[-IDX,i])/2
            Z[IDX,i]<-rbinom(length(IDX),2,maf1)
         }
      }
   } else if(impute.method =="fixed"){
      for(i in 1:p){
         IDX<-which(is.na(Z[,i]))
         if(length(IDX) > 0){
            maf1<-mean(Z[-IDX,i])/2
            Z[IDX,i]<-2 * maf1
         }
      }
   } else if(impute.method =="bestguess") {
      for(i in 1:p){
         IDX<-which(is.na(Z[,i]))
         if(length(IDX) > 0){
            maf1<-mean(Z[-IDX,i])/2
            Z[IDX,i]<-round(2 * maf1)
         }
      }
   } else {
      stop("Error: Imputation method shoud be \"fixed\", \"random\" or \"bestguess\" ")
   }
   return(as.matrix(Z))
}
