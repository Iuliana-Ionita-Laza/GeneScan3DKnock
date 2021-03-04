#' @importFrom stats binomial dbeta gaussian glm pcauchy pchisq rbinom sd var median
utils::globalVariables(c('G_Enhancer1_surround','G_Enhancer2_surround','variants_Enhancer1_surround',
                         'variants_Enhancer2_surround','Enhancer1.pos','Enhancer2.pos','create.MK',
'KnockoffGeneration.example','GeneScan3DKnock','GeneScan3DKnock.example',
'G_EnhancerAll','Z_EnhancerAll','p_EnhancerAll',"G_gene_buffer", "Z_gene_buffer", 
"pos_gene_buffer",'n','G_promoter','Z_promoter','G_Enhancer1','Z_Enhancer1','G_Enhancer2','Z_Enhancer2'))


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


#' Data example for AR Knockoff Generation.
#'
#'This simulated example dataset contains outcome variable Y, covariate X, genotype matrices and genetic variants of surrounding regions for gene buffer and two enhancers separately, positions and functional annotations for gene buffer region, promoter and two enhancers.
#'
#'We provide genotypes of 20 Kb surrounding regions for 15 Kb gene buffer region and two 2 Kb enhancers separately. In real data analyses, the surrounding regions can increase to 200 Kb for knockoff generation.
#'
#' @name Example.KnockoffGeneration
#' @docType data
#' @keywords data
#' @usage data("KnockoffGeneration.example")
#' @examples
#' data("KnockoffGeneration.example")
#'
#'Y=KnockoffGeneration.example$Y; X=KnockoffGeneration.example$X; 
#'
#'G_gene_buffer_surround=KnockoffGeneration.example$G_gene_buffer_surround
#'variants_gene_buffer_surround=KnockoffGeneration.example$variants_gene_buffer_surround
#'G_Enhancer1_surround=KnockoffGeneration.example$G_Enhancer1_surround
#'variants_Enhancer1_surround=KnockoffGeneration.example$variants_Enhancer1_surround
#'G_Enhancer2_surround=KnockoffGeneration.example$G_Enhancer2_surround
#'variants_Enhancer2_surround=KnockoffGeneration.example$variants_Enhancer2_surround
#'
#'G_EnhancerAll_surround=cbind(G_Enhancer1_surround,G_Enhancer2_surround)
#'variants_EnhancerAll_surround=c(variants_Enhancer1_surround,variants_Enhancer2_surround)
#'p_EnhancerAll_surround=c(length(variants_Enhancer1_surround),length(variants_Enhancer2_surround))
#'
#'gene_buffer.pos=KnockoffGeneration.example$gene_buffer.pos
#'promoter.pos=KnockoffGeneration.example$promoter.pos
#'Enhancer1.pos=KnockoffGeneration.example$Enhancer1.pos
#'Enhancer2.pos=KnockoffGeneration.example$Enhancer2.pos   
#'Enhancer.pos=rbind(Enhancer1.pos,Enhancer2.pos)
#'
#'Z_gene_buffer=KnockoffGeneration.example$Z_gene_buffer
#'Z_promoter=KnockoffGeneration.example$Z_promoter
#'Z_Enhancer1=KnockoffGeneration.example$Z_Enhancer1
#'Z_Enhancer2=KnockoffGeneration.example$Z_Enhancer2
#'Z_EnhancerAll=rbind(Z_Enhancer1,Z_Enhancer2)
#'p_EnhancerAll=c(dim(Z_Enhancer1)[1],dim(Z_Enhancer2)[1])
#'
#'R=KnockoffGeneration.example$R
"KnockoffGeneration.example"


#'Data example for GeneScan3DKnock.
#'
#'This example dataset contains the original and M=5 knockoff p-values for N=100 genes. Each row presents gene id, original GeneScan3D p-value and M knockoff GeneScan3D p-values. The original and knockoff GeneScan3D p-values are generated using GeneScan3D.KnockoffGeneration() function. 
#'
#'This example dataset can be used to calculate the knockoff statistics and q-values for GeneScan3DKnock() function. 
#'
#' @name Example.GeneScan3DKnock
#' @docType data
#' @keywords data
#' @usage data("GeneScan3DKnock.example")
"GeneScan3DKnock.example"


#' The preliminary data management for GeneScan3DKnock.
#'
#' This function does the preliminary data management and fit the model under null hypothesis using all the covariates. The output will be used in the other GeneScan functions.
#'
#' @param Y The outcome variable, an n*1 matrix where n is the number of individuals.
#' @param X An n*d covariates matrix where d is the number of covariates.
#' @param id The subject id. This is used to match phenotype with genotype. The default is NULL, where the matched phenotype and genotype matrices are assumed.
#' @param out_type Type of outcome variable. Can be either "C" for continuous or "D" for dichotomous. The default is "C".
#' @param B Number of resampling replicates. The default is 1000. A larger value leads to more accurate and stable p-value calculation, but requires more computing time.
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
#'set.seed(12345)
#'result.null.model=GeneScan.Null.Model(Y, X, out_type="C", B=1000)
#'
#' @export
GeneScan.Null.Model<-function(Y, X=NULL, id=NULL, out_type="C", B=1000){
   
   Y<-as.matrix(Y);n<-nrow(Y)
   
   if(length(X)!=0){X0<-svd(as.matrix(X))$u}else{X0<-NULL}
   X0<-cbind(rep(1,n),X0)
   
   if(out_type=="C"){nullglm<-glm(Y~0+X0,family=gaussian)}
   if(out_type=="D"){nullglm<-glm(Y~0+X0,family=binomial)}
   
   if (length(id)==0){id<-1:n}
   
   mu<-nullglm$fitted.values;Y.res<-Y-mu
   #permute the residuals for B times
   index<-sapply(1:B,function(x)sample(1:length(Y)));temp.Y.res<-Y.res[as.vector(index)]
   re.Y.res<-matrix(temp.Y.res,length(Y),B)
   
   #prepare invserse matrix for covariates
   if(out_type=='D'){v<-mu*(1-mu)}else{v<-rep(as.numeric(var(Y.res)),length(Y))}
   inv.X0<-solve(t(X0)%*%(v*X0))
   
   #prepare the preliminary features
   result.null.model<-list(Y=Y,id=id,n=n,X0=X0,nullglm=nullglm,out_type=out_type,re.Y.res=re.Y.res,inv.X0=inv.X0)
   return(result.null.model)
}


#' Conduct GeneScan1D analysis on the gene buffer region.
#'
#' This function conducts gene-based scan test on the gene buffer region using 1D windows with sizes 1-5-10 Kb.
#'
#' @param G The genotype matrix in the gene buffer region, which is a n*p matrix where n is the number of individuals and p is the number of genetic variants in the gene buffer region.
#' @param Z A p*q functional annotation matrix where p is the number of genetic variants in the gene buffer region and q is the number of functional annotations. If Z is NULL (do not incorporate any functional annotations), the minor allele frequency weighted dispersion and/or burden tests are applied. Specifically, Beta(MAF; 1; 25) weights are used for rare variants and weights one are used for common variants.
#' @param window.size The 1-D window sizes in base pairs to scan the gene buffer region. The recommended window sizes are c(1000,5000,10000).
#' @param pos  The positions of genetic variants in the gene buffer region, an p dimensional vector. Each position corresponds to a column in the genotype matrix and a row in the functional annotation matrix.
#' @param MAC.threshold Threshold for minor allele count. Variants below MAC.threshold are ultra-rare variants. The recommended level is 5.
#' @param MAF.threshold Threshold for minor allele frequency. Variants below MAF.threshold are rare variants. The recommended level is 0.01.
#' @param Gsub.id The subject id corresponding to the genotype matrix, an n dimensional vector. The default is NULL, where the matched phenotype and genotype matrices are assumed.
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
#'# pos_gene_buffer: positions of p=287 genetic variants
#'# Z_gene_buffer: p by q functional annotation matrix, q=1 functional annotation
#'
#'data("GeneScan3D.example")
#'Y=GeneScan3D.example$Y; X=GeneScan3D.example$X;
#'G_gene_buffer=GeneScan3D.example$G_gene_buffer;
#'Z_gene_buffer=GeneScan3D.example$Z_gene_buffer;
#'pos_gene_buffer=GeneScan3D.example$pos_gene_buffer;
#'
#'# Preliminary data management
#'set.seed(12345)
#'result.null.model=GeneScan.Null.Model(Y, X, out_type="C", B=1000)
#'
#'#Conduct GeneScan1D analysis
#'result.GeneScan1D=GeneScan1D(G=G_gene_buffer,Z=Z_gene_buffer,pos=pos_gene_buffer,
#'                             window.size=c(1000,5000,10000),MAC.threshold=5,MAF.threshold=0.01,
#'                             result.null.model=result.null.model)
#'result.GeneScan1D$GeneScan1D.Cauchy.pvalue 
#'
#' @import SKAT
#' @import Matrix
#' @import WGScan
#' @import SPAtest
#' @import CompQuadForm
#' @export
GeneScan1D<-function(G=G_gene_buffer,Z=Z_gene_buffer,window.size=c(1000,5000,10000), pos=pos_gene_buffer,
                     MAC.threshold=5,MAF.threshold=0.01,Gsub.id=NULL,result.null.model=result.null.model){
   
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
   
   GeneScan1D.Cauchy.window=matrix(NA,nrow=M_gene_buffer,ncol=3)
   for (m in 1:M_gene_buffer){
      #Create index for each window
      index.wondow<-(window.matrix_gene_buffer[,m]==1)
      G.window=G[,index.wondow]
      G.window=Matrix(G.window)
      
      if(!is.null(Z)){
         Z.window=Z[index.wondow,]
         Z.window=Matrix(Z.window)
      }else{
         Z.window=NULL
      }
      
      #if there is only 1 variant in this window, then do not conduct combined test in this window
      #move to the next one
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
      
      #Burden test
      p.burden<-matrix(NA,1,ncol(weight.matrix))
      for (j in 1:ncol(weight.matrix)){
         temp.window.matrix<-weight.matrix[,j]
         X<-as.matrix(G.window%*%temp.window.matrix)
         p.burden[,j]<-Get.p.base(X,result.null.model)
      }
      
      #SKAT test
      p.dispersion<-matrix(NA,1,ncol(weight.matrix))
      if(outcome=='D'){v<-mu*(1-mu)}else{v<-rep(as.numeric(var(Y.res)),nrow(G.window))}
      A<-t(G.window)%*%(v*G.window)
      B<-t(G.window)%*%(v*X0)
      C<-solve(t(X0)%*%(v*X0))
      K<-A-B%*%C%*%t(B)
      score<-t(G.window)%*%Y.res;re.score<-t(t(G.window)%*%re.Y.res) #resampling for 1000 times
      for (j in 1:ncol(weight.matrix)){
         #For extremely rare variants, do not conduct SKAT
         p.dispersion[,j]<-Get.p.SKAT(score,re.score,K,window.matrix=as.matrix(rep(1,sum(index.wondow))),weight=(MAC.window>=MAC.threshold)*weight.matrix[,j],result.null.model) 
      }  
      
      #Single variant score test for all variants in the window
      p.single<-Get.p(G.window,result.null.model) 
      p.individual1<-Get.cauchy.scan(p.single,as.matrix((MAC.window>=MAC.threshold & MAF.window<MAF.threshold))) #rare variants
      p.individual2<-Get.cauchy.scan(p.single,as.matrix((MAF.window>=MAF.threshold))) #common and low frequency variants
      p.individual<-cbind(p.burden,p.dispersion,p.individual1,p.individual2);
      colnames(p.individual)<-c(paste0('burden_',colnames(weight.matrix)),paste0('dispersion_',colnames(weight.matrix)),'singleCauchy_MAF<MAF.threshold&MAC>=MAC.threshold','singleCauchy_MAF>=MAF.threshold')
      
      #aggregated Cauchy association test
      p.Cauchy<-as.matrix(apply(p.individual,1,Get.cauchy))
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
#' This function conducts gene-based scan test on the gene buffer region, integrating proximal and distal regulatory elements for a gene, i.e., promoter and R enhancers.
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
#'result.null.model=GeneScan.Null.Model(Y, X, out_type="C", B=1000)
#'
#'# Conduct GeneScan3D analysis
#'result.GeneScan3D=GeneScan3D(G=G_gene_buffer,Z=Z_gene_buffer,
#'                             G.promoter=G_promoter,Z.promoter=Z_promoter,
#'                             G.EnhancerAll=G_EnhancerAll,Z.EnhancerAll=Z_EnhancerAll, 
#'                             R=2,p_Enhancer=p_EnhancerAll,
#'                             pos=pos_gene_buffer,
#'                             window.size=c(1000,5000,10000),MAC.threshold=5,MAF.threshold=0.01,
#'                             result.null.model=result.null.model)
#'result.GeneScan3D$GeneScan3D.Cauchy.pvalue
#'
#' @import SKAT
#' @import Matrix
#' @import WGScan
#' @import SPAtest
#' @import CompQuadForm
#' @export
GeneScan3D<-function(G=G_gene_buffer,Z=Z_gene_buffer,G.promoter=G_promoter,Z.promoter=Z_promoter,G.EnhancerAll=G_EnhancerAll,Z.EnhancerAll=Z.EnhancerAll, R=length(p_EnhancerAll),
                     p_Enhancer=p_EnhancerAll,window.size=c(1000,5000,10000),pos=pos_gene_buffer,
                     MAC.threshold=5,MAF.threshold=0.01,Gsub.id=NULL,result.null.model=result.null.model){
   
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
   
   GeneScan1D.Cauchy.window=matrix(NA,nrow=M_gene_buffer,ncol=3)
   for (m in 1:M_gene_buffer){
      #Create index for each window
      index.wondow<-(window.matrix_gene_buffer[,m]==1)
      G.window=G[,index.wondow]
      G.window=Matrix(G.window)
      
      if(!is.null(Z)){
         Z.window=Z[index.wondow,]
         Z.window=Matrix(Z.window)
      }else{
         Z.window=NULL
      }
      
      #if there is only 1 variant in this window, then do not conduct combined test in this window
      #move to the next one
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
      
      #Burden test
      p.burden<-matrix(NA,1,ncol(weight.matrix))
      for (j in 1:ncol(weight.matrix)){
         temp.window.matrix<-weight.matrix[,j]
         X<-as.matrix(G.window%*%temp.window.matrix)
         p.burden[,j]<-Get.p.base(X,result.null.model)
      }
      
      #SKAT test
      p.dispersion<-matrix(NA,1,ncol(weight.matrix))
      if(outcome=='D'){v<-mu*(1-mu)}else{v<-rep(as.numeric(var(Y.res)),nrow(G.window))}
      A<-t(G.window)%*%(v*G.window)
      B<-t(G.window)%*%(v*X0)
      C<-solve(t(X0)%*%(v*X0))
      K<-A-B%*%C%*%t(B)
      score<-t(G.window)%*%Y.res;re.score<-t(t(G.window)%*%re.Y.res) #resampling for 1000 times
      for (j in 1:ncol(weight.matrix)){
         #For extremely rare variants, do not conduct SKAT
         p.dispersion[,j]<-Get.p.SKAT(score,re.score,K,window.matrix=as.matrix(rep(1,sum(index.wondow))),weight=(MAC.window>=MAC.threshold)*weight.matrix[,j],result.null.model)
      }
      
      #Single variant score test for all variants in the window
      p.single<-Get.p(G.window,result.null.model)
      p.individual1<-Get.cauchy.scan(p.single,as.matrix((MAC.window>=MAC.threshold & MAF.window<MAF.threshold))) #rare variants
      p.individual2<-Get.cauchy.scan(p.single,as.matrix((MAF.window>=MAF.threshold))) #common and low frequency variants
      p.individual<-cbind(p.burden,p.dispersion,p.individual1,p.individual2);
      colnames(p.individual)<-c(paste0('burden_',colnames(weight.matrix)),paste0('dispersion_',colnames(weight.matrix)),'singleCauchy_MAF<MAF.threshold&MAC>=MAC.threshold','singleCauchy_MAF>=MAF.threshold')
      
      #aggregated Cauchy association test
      p.Cauchy<-as.matrix(apply(p.individual,1,Get.cauchy))
      test.common<-grep('MAF>=MAF.threshold',colnames(p.individual))
      p.Cauchy.common<-as.matrix(apply(p.individual[,test.common,drop=FALSE],1,Get.cauchy))
      p.Cauchy.rare<-as.matrix(apply(p.individual[,-test.common,drop=FALSE],1,Get.cauchy))
      
      GeneScan1D.Cauchy.window[m,]=c(p.Cauchy,p.Cauchy.common,p.Cauchy.rare)
   }
   
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
         
         #Burden test
         p.burden.promoter<-matrix(NA,1,ncol(weight.matrix))
         for (j in 1:ncol(weight.matrix)){
            temp.window.matrix<-weight.matrix[,j]
            X<-as.matrix(G.window.promoter%*%temp.window.matrix)
            p.burden.promoter[,j]<-Get.p.base(X,result.null.model)
         }
         
         #SKAT test
         p.dispersion.promoter<-matrix(NA,1,ncol(weight.matrix))
         if(outcome=='D'){v<-mu*(1-mu)}else{v<-rep(as.numeric(var(Y.res)),nrow(G.window.promoter))}
         A<-t(G.window.promoter)%*%(v*G.window.promoter)
         B<-t(G.window.promoter)%*%(v*X0)
         C<-solve(t(X0)%*%(v*X0))
         K<-A-B%*%C%*%t(B)
         score<-t(G.window.promoter)%*%Y.res;re.score<-t(t(G.window.promoter)%*%re.Y.res)
         for (j in 1:ncol(weight.matrix)){
            #For extremely rare variants, do not conduct SKAT
            p.dispersion.promoter[,j]<-Get.p.SKAT(score,re.score,K,window.matrix=as.matrix(rep(1,dim(G.window.promoter)[2])),weight=(MAC.window.promoter>=MAC.threshold)*weight.matrix[,j],result.null.model)
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
         print(r)
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
            
            #Burden test
            p.burden.Enhancer<-matrix(NA,1,ncol(weight.matrix))
            for (j in 1:ncol(weight.matrix)){
               temp.window.matrix<-weight.matrix[,j]
               X<-as.matrix(G.window.Enhancer%*%temp.window.matrix)
               p.burden.Enhancer[,j]<-Get.p.base(X,result.null.model)
            }
            
            #SKAT test
            p.dispersion.Enhancer<-matrix(NA,1,ncol(weight.matrix))
            if(outcome=='D'){v<-mu*(1-mu)}else{v<-rep(as.numeric(var(Y.res)),nrow(G.window.Enhancer))}
            A<-t(G.window.Enhancer)%*%(v*G.window.Enhancer)
            B<-t(G.window.Enhancer)%*%(v*X0)
            C<-solve(t(X0)%*%(v*X0))
            K<-A-B%*%C%*%t(B)
            score<-t(G.window.Enhancer)%*%Y.res;re.score<-t(t(G.window.Enhancer)%*%re.Y.res) #resampling for 1000 times
            
            for (j in 1:ncol(weight.matrix)){
               #For extremely rare variants, do not conduct SKAT
               p.dispersion.Enhancer[,j]<-Get.p.SKAT(score,re.score,K,window.matrix=as.matrix(rep(1,dim(G.window.Enhancer)[2])),weight=(MAC.window.Enhancer>=MAC.threshold)*weight.matrix[,j],result.null.model)
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
   
   return(list(GeneScan3D.Cauchy.pvalue=GeneScan3D.Cauchy,M=M_gene_buffer,R=sum(Enhancer_ind),
               minp=c(min(GeneScan3D.Cauchy.RE[,1],na.rm=TRUE),min(GeneScan3D.Cauchy.RE[,2],na.rm=TRUE),min(GeneScan3D.Cauchy.RE[,3],na.rm=TRUE)),
               RE_minp=c(RE_minp.all,RE_minp.common,RE_minp.rare)))
}




#' GeneScan3D AR Knockoff Generation: an auto-regressive model for knockoff generation. 
#'
#' This function generates multiple knockoff genotypes for a gene and the corresponding regulatory elements based on an auto-regressive model.  Additionally, it computes p-values from the GeneScan3D test for a gene based on the original data, and each of the knockoff replicates.
#'
#' @param M Numer of multiple knockoffs.
#' @param G_gene_buffer_surround The genotype matrix of the surrounding region for gene buffer region. 
#' @param variants_gene_buffer_surround The genetic variants in the surrounding region for gene buffer region. Each position corresponds to a column in the genotype matrix G_gene_buffer_surround.
#' @param gene_buffer.pos The start and end positions of gene buffer region.
#' @param promoter.pos The start and end positions of promoter.
#' @param R Number of enhancers.
#' @param G_EnhancerAll_surround The genotype matrix of the surrounding regions for R enhancers, by combining the genotype matrix of the surrounding regions for each enhancer by columns.
#' @param variants_EnhancerAll_surround The genetic variants in the surrounding region for R enhancers. Each position corresponds to a column in the genotype matrix G_EnhancerAll_surround. 
#' @param p_EnhancerAll_surround Number of genetic variants in the surrounding region for R enhancers, which is a 1*R vector.
#' @param Enhancer.pos The start and end positions for R enhancers. One row represents one enhancer, which is a R by 2 matrix. 
#' @param p.EnhancerAll Number of genetic variants in R enhancers, which is a 1*R vector.
#' @param Z A p*q functional annotation matrix, where p is the number of genetic variants in the gene buffer region and q is the number of functional annotations. If Z is NULL (do not incorporate any functional annotations), the minor allele frequency weighted dispersion and/or burden tests are applied. Specifically, Beta(MAF; 1; 25) weights are used for rare variants and weights one are used for common variants.
#' @param Z.promoter The functional annotation matrix for promoter. Z.promoter can be NULL.
#' @param Z.EnhancerAll The functional annotation matrix for R enhancers, by combining the functional annotation matrix of each enhancer by rows. Z.EnhancerAll can be NULL.
#' @param window.size The 1-D window sizes in base pairs to scan the gene buffer region. The recommended window sizes are c(1000,5000,10000).
#' @param MAC.threshold Threshold for minor allele count. Variants below MAC.threshold are ultra-rare variants. The recommended level is 5.
#' @param MAF.threshold Threshold for minor allele frequency. Variants below MAF.threshold are rare variants. The recommended level is 0.01.
#' @param Gsub.id The subject id corresponding to the genotype matrix, an n dimensional vector. The default is NULL, where the matched phenotype and genotype matrices are assumed.
#' @param result.null.model The output of function "GeneScan.Null.Model()".
#' @return \item{GeneScan3D.Cauchy}{GeneScan3D p-values of all, common and rare variants for original genotypes.}
#' @return \item{GeneScan3D.Cauchy_knockoff}{A M by 3 GeneScan3D p-values matrix of all, common and rare variants for M knockoff genotypes.}
#' @examples
#'library(GeneScan3DKnock)
#'data(KnockoffGeneration.example)
#'Y=KnockoffGeneration.example$Y; X=KnockoffGeneration.example$X; 
#'
#'G_gene_buffer_surround=KnockoffGeneration.example$G_gene_buffer_surround
#'variants_gene_buffer_surround=KnockoffGeneration.example$variants_gene_buffer_surround
#'G_Enhancer1_surround=KnockoffGeneration.example$G_Enhancer1_surround
#'variants_Enhancer1_surround=KnockoffGeneration.example$variants_Enhancer1_surround
#'G_Enhancer2_surround=KnockoffGeneration.example$G_Enhancer2_surround
#'variants_Enhancer2_surround=KnockoffGeneration.example$variants_Enhancer2_surround
#'
#'G_EnhancerAll_surround=cbind(G_Enhancer1_surround,G_Enhancer2_surround)
#'variants_EnhancerAll_surround=c(variants_Enhancer1_surround,variants_Enhancer2_surround)
#'p_EnhancerAll_surround=c(length(variants_Enhancer1_surround),length(variants_Enhancer2_surround))
#'
#'gene_buffer.pos=KnockoffGeneration.example$gene_buffer.pos
#'promoter.pos=KnockoffGeneration.example$promoter.pos
#'Enhancer1.pos=KnockoffGeneration.example$Enhancer1.pos
#'Enhancer2.pos=KnockoffGeneration.example$Enhancer2.pos   
#'Enhancer.pos=rbind(Enhancer1.pos,Enhancer2.pos)
#'
#'Z_gene_buffer=KnockoffGeneration.example$Z_gene_buffer
#'Z_promoter=KnockoffGeneration.example$Z_promoter
#'Z_Enhancer1=KnockoffGeneration.example$Z_Enhancer1
#'Z_Enhancer2=KnockoffGeneration.example$Z_Enhancer2
#'Z_EnhancerAll=rbind(Z_Enhancer1,Z_Enhancer2)
#'p_EnhancerAll=c(dim(Z_Enhancer1)[1],dim(Z_Enhancer2)[1])
#'
#'R=KnockoffGeneration.example$R
#'
#'set.seed(12345)
#'result.null.model=GeneScan.Null.Model(Y, X, out_type="C", B=1000)
#'
#'result.GeneScan3D.KnockoffGeneration=GeneScan3D.KnockoffGeneration(
#'G_gene_buffer_surround=G_gene_buffer_surround,
#'variants_gene_buffer_surround=variants_gene_buffer_surround,
#'gene_buffer.pos=gene_buffer.pos,promoter.pos=promoter.pos,R=R, 
#'G_EnhancerAll_surround=G_EnhancerAll_surround, 
#'variants_EnhancerAll_surround=variants_EnhancerAll_surround,
#'p_EnhancerAll_surround=p_EnhancerAll_surround,
#'Enhancer.pos=Enhancer.pos,p.EnhancerAll=p_EnhancerAll,
#'Z=Z_gene_buffer,Z.promoter=Z_promoter,Z.EnhancerAll=Z_EnhancerAll, 
#'window.size=c(1000,5000,10000),
#'MAC.threshold=5,MAF.threshold=0.01,Gsub.id=NULL,result.null.model=result.null.model,M=5)
#'result.GeneScan3D.KnockoffGeneration$GeneScan3D.Cauchy 
#'result.GeneScan3D.KnockoffGeneration$GeneScan3D.Cauchy_knockoff
#' @import SKAT
#' @import Matrix
#' @import WGScan
#' @import SPAtest
#' @import CompQuadForm
#' @import abind
#' @import KnockoffScreen
#' @export
GeneScan3D.KnockoffGeneration=function(G_gene_buffer_surround=G_gene_buffer_surround,
                                       variants_gene_buffer_surround=variants_gene_buffer_surround,
                                       gene_buffer.pos=gene_buffer.pos,promoter.pos=promoter.pos,R=R,
                                       G_EnhancerAll_surround=G_EnhancerAll_surround,
                                       variants_EnhancerAll_surround=variants_EnhancerAll_surround,
                                       p_EnhancerAll_surround=p_EnhancerAll_surround,
                                       Enhancer.pos=Enhancer.pos,p.EnhancerAll=p_EnhancerAll,
                                       Z=Z_gene_buffer,Z.promoter=Z_promoter,Z.EnhancerAll=Z_EnhancerAll,
                                       window.size=c(1000,5000,10000),
                                       MAC.threshold=5,MAF.threshold=0.01,Gsub.id=NULL,result.null.model=result.null.model,M=5){
   
   
   
   mu<-result.null.model$nullglm$fitted.values;
   Y.res<-result.null.model$Y-mu
   re.Y.res<-result.null.model$re.Y.res 
   X0<-result.null.model$X0
   outcome<-result.null.model$out_type
   
   impute.method='fixed'
   ## Prelimanry checking and filtering the variants
   #match phenotype id and genotype id
   if(length(Gsub.id)==0){match.index<-match(result.null.model$id,1:nrow(G_gene_buffer_surround))}else{
      match.index<-match(result.null.model$id,Gsub.id)
   }
   if(mean(is.na(match.index))>0){
      msg<-sprintf("Some individuals are not matched with genotype. The rate is%f", mean(is.na(match.index)))
      warning(msg,call.=F)
   }
   #individuals ids are matched with genotype
   G_gene_buffer_surround=Matrix(G_gene_buffer_surround[match.index,])
   #missing genotype imputation
   G_gene_buffer_surround[G_gene_buffer_surround==-9 | G_gene_buffer_surround==9]=NA
   N_MISS=sum(is.na(G_gene_buffer_surround))
   MISS.freq=apply(is.na(G_gene_buffer_surround),2,mean)
   if(N_MISS>0){
      msg<-sprintf("The missing genotype rate is %f. Imputation is applied.", N_MISS/nrow(G_gene_buffer_surround)/ncol(G_gene_buffer_surround))
      warning(msg,call.=F)
      G_gene_buffer_surround=Impute(G_gene_buffer_surround,impute.method)
   }
   
   #MAF filtering
   MAF<-apply(G_gene_buffer_surround,2,mean)/2 #MAF of nonfiltered variants
   G_gene_buffer_surround[,MAF>0.5 & !is.na(MAF)]<-2-G_gene_buffer_surround[,MAF>0.5 & !is.na(MAF)]
   MAF<-apply(G_gene_buffer_surround,2,mean)/2
   MAC<-apply(G_gene_buffer_surround,2,sum) #minor allele count
   s<-apply(G_gene_buffer_surround,2,sd)
   SNP.index<-which(MAF>0 & s!=0 & !is.na(MAF) & MISS.freq<0.1) 
   
   check.index<-which(MAF>0 & s!=0 & !is.na(MAF)  & MISS.freq<0.1)
   if(length(check.index)<=1){
      warning('Number of variants with missing rate <=10% in the gene is <=1')
   }
   
   G_gene_buffer_surround<-Matrix(G_gene_buffer_surround[,SNP.index])
   variants_gene_buffer_surround_filter=variants_gene_buffer_surround[SNP.index]
   
   ###Generate multiple knockoffs
   G_gene_buffer_surround_knockoff<-create.MK(G_gene_buffer_surround,pos=variants_gene_buffer_surround_filter,M=M,corr_max=0.75)
   
   ##obtain knockoff genotypes for gene buffer region and promoter
   positions_gene_buffer=variants_gene_buffer_surround_filter[variants_gene_buffer_surround_filter<=gene_buffer.pos[2]&variants_gene_buffer_surround_filter>=gene_buffer.pos[1]]
   G_gene_buffer=G_gene_buffer_surround[,variants_gene_buffer_surround_filter%in%positions_gene_buffer]
   
   G_promoter=NULL
   if(!is.null(promoter.pos)){
      positions_promoter=variants_gene_buffer_surround_filter[variants_gene_buffer_surround_filter<=promoter.pos[2]&variants_gene_buffer_surround_filter>=promoter.pos[1]]
      G_promoter=G_gene_buffer_surround[,variants_gene_buffer_surround_filter%in%positions_promoter]
   }
   
   ##functional annotation
   Z_gene_buffer=NULL
   if(!is.null(Z)){
      positions_gene_buffer_nonfilter=variants_gene_buffer_surround[variants_gene_buffer_surround<=gene_buffer.pos[2]&variants_gene_buffer_surround>=gene_buffer.pos[1]]
      Z_gene_buffer=as.matrix(Z[positions_gene_buffer_nonfilter%in%positions_gene_buffer,])
   }
   Z_promoter=NULL
   if(!is.null(Z.promoter)){
      positions_promoter_nonfilter=variants_gene_buffer_surround[variants_gene_buffer_surround<=promoter.pos[2]&variants_gene_buffer_surround>=promoter.pos[1]]
      Z_promoter=as.matrix(Z.promoter[positions_promoter_nonfilter%in%positions_promoter,])
   }
   
   ## R enhancers ##
   G_EnhancerAll=c()
   p_EnhancerAll=c()
   Z_EnhancerAll=c()
   G_EnhancerAll_knockoff=c()
   
   if (R!=0){
      for (r in 1:R){
         print(r)
         ##genotype Enhancer_surround_region
         if (r==1){
            G_Enhancer_surround=G_EnhancerAll_surround[,1:cumsum(p_EnhancerAll_surround)[r]]
            positions_Enhancer_surround=variants_EnhancerAll_surround[1:cumsum(p_EnhancerAll_surround)[r]]
         }else{
            G_Enhancer_surround=G_EnhancerAll_surround[,(cumsum(p_EnhancerAll_surround)[r-1]+1):cumsum(p_EnhancerAll_surround)[r]]
            positions_Enhancer_surround=variants_EnhancerAll_surround[(cumsum(p_EnhancerAll_surround)[r-1]+1):cumsum(p_EnhancerAll_surround)[r]]
         }
         
         #individuals ids are matched with genotype
         G_Enhancer_surround=Matrix(G_Enhancer_surround[match.index,])
         #missing genotype imputation
         G_Enhancer_surround[G_Enhancer_surround==-9 | G_Enhancer_surround==9]=NA
         N_MISS=sum(is.na(G_Enhancer_surround))
         MISS.freq=apply(is.na(G_Enhancer_surround),2,mean)
         if(N_MISS>0){
            msg<-sprintf("The missing genotype rate is %f. Imputation is applied.", N_MISS/nrow(G_Enhancer_surround)/ncol(G_Enhancer_surround))
            warning(msg,call.=F)
            G_Enhancer_surround=Impute(G_Enhancer_surround,impute.method)
         }
         
         #MAF filtering
         MAF<-apply(G_Enhancer_surround,2,mean)/2 #MAF of nonfiltered variants
         G_Enhancer_surround[,MAF>0.5 & !is.na(MAF)]<-2-G_Enhancer_surround[,MAF>0.5 & !is.na(MAF)]
         MAF<-apply(G_Enhancer_surround,2,mean)/2
         MAC<-apply(G_Enhancer_surround,2,sum) #minor allele count
         s<-apply(G_Enhancer_surround,2,sd)
         SNP.index<-which(MAF>0 & s!=0 & !is.na(MAF) & MISS.freq<0.1) 
         length(SNP.index)
         
         check.index<-which(MAF>0 & s!=0 & !is.na(MAF)  & MISS.freq<0.1)
         if(length(check.index)<=1){
            warning('Number of variants with missing rate <=10% in the gene is <=1')
         }
         
         G_Enhancer_surround<-Matrix(G_Enhancer_surround[,SNP.index])
         positions_Enhancer_surround_filter=positions_Enhancer_surround[SNP.index]
         
         G_Enhancer_surround_knockoff<-create.MK(G_Enhancer_surround,pos=positions_Enhancer_surround_filter,M=M,corr_max=0.75)
         
         positions_enhancer=positions_Enhancer_surround_filter[positions_Enhancer_surround_filter<=Enhancer.pos[r,2]&positions_Enhancer_surround_filter>=Enhancer.pos[r,1]]
         G_enhancer=Matrix(G_Enhancer_surround[,positions_Enhancer_surround_filter%in%positions_enhancer])
         G_EnhancerAll=cbind(G_EnhancerAll,G_enhancer)
         
         p_Enhancer=length(positions_enhancer)
         p_EnhancerAll=c(p_EnhancerAll,p_Enhancer)
         
         G_enhancer_knockoff=array(0, dim = c(M, result.null.model$n, p_Enhancer))
         ##M knockoffs
         for(k in 1:M){
            G_Enhancer_surround_knockoff_k=G_Enhancer_surround_knockoff[k,,]
            G_enhancer_knockoff[k,,]=G_Enhancer_surround_knockoff_k[,positions_Enhancer_surround_filter%in%positions_enhancer]
         }
         G_EnhancerAll_knockoff=abind(G_EnhancerAll_knockoff,G_enhancer_knockoff)
         
         ##functional annotation
         if(!is.null(Z.EnhancerAll)){
            if (r==1){Z_Enhancer=as.matrix(Z.EnhancerAll[1:cumsum(p.EnhancerAll)[r],])}else{
               Z_Enhancer=as.matrix(Z.EnhancerAll[(cumsum(p.EnhancerAll)[r-1]+1):cumsum(p.EnhancerAll)[r],])}
            
            Z_Enhancer=as.matrix(Z_Enhancer[positions_Enhancer_surround[positions_Enhancer_surround<=Enhancer.pos[r,2]&positions_Enhancer_surround>=Enhancer.pos[r,1]]%in%positions_enhancer,])
            Z_EnhancerAll=rbind(Z_EnhancerAll,Z_Enhancer)
         }
      }
   }
   
   ####GeneScan3D: conduct gene-based test on the gene buffer region, adding one promoter and R enhancers ################
   ##original p-values
   GeneScan3D.Cauchy=GeneScan3D(G=G_gene_buffer,Z=Z_gene_buffer,G.promoter=G_promoter,Z.promoter=Z_promoter,
                                G.EnhancerAll=G_EnhancerAll,Z.EnhancerAll=Z_EnhancerAll, R=R,
                                p_Enhancer=p_EnhancerAll,window.size=c(1000,5000,10000),pos=positions_gene_buffer,
                                MAC.threshold=5,MAF.threshold=0.01,result.null.model=result.null.model,Gsub.id=row.names(G_gene_buffer))$GeneScan3D.Cauchy.pvalue
   
   #M knockoff p-values 
   GeneScan3D.Cauchy_knockoff=matrix(NA,nrow=M,ncol=3)
   for (k in 1:M){
      G_gene_buffer_surround_knockoff_k=G_gene_buffer_surround_knockoff[k,,]
      G_gene_buffer_knockoff_k=G_gene_buffer_surround_knockoff_k[,variants_gene_buffer_surround_filter%in%positions_gene_buffer]
      
      G_promoter_knockoff_k=NULL
      if(!is.null(Z.promoter)){
         G_promoter_knockoff_k=G_gene_buffer_surround_knockoff_k[,variants_gene_buffer_surround_filter%in%positions_promoter]
      }
      
      GeneScan3D.Cauchy_knockoff[k,]=GeneScan3D(G=G_gene_buffer_knockoff_k,Z=Z_gene_buffer,G.promoter=G_promoter_knockoff_k,Z.promoter=Z_promoter,
                                                G.EnhancerAll=G_EnhancerAll_knockoff[k,,],Z.EnhancerAll=Z_EnhancerAll, R=R,
                                                p_Enhancer=p_EnhancerAll,window.size=c(1000,5000,10000),pos=positions_gene_buffer,
                                                MAC.threshold=5,MAF.threshold=0.01,result.null.model=result.null.model,Gsub.id=row.names(G_gene_buffer))$GeneScan3D.Cauchy.pvalue
   }
   return(list(GeneScan3D.Cauchy=GeneScan3D.Cauchy,GeneScan3D.Cauchy_knockoff=GeneScan3D.Cauchy_knockoff))
}

#' GeneScan3DKnock: Knockoff-enhanced gene-based test for causal gene discovery (knockoff filter).
#'
#' This function performs the knockoff filter, and computes the q-value for each gene. This function takes the results from the GeneScan3D.KnockoffGeneration() function and get knockoff statistics and q-values.
#'
#' @param M Number of multiple knockoffs. 
#' @param p0 A N-dimensional vector of the original GeneScan3D p-values, calculated using GeneScan3D.KnockoffGeneration() function.
#' @param p_ko A N*M matrix of M knockoff GeneScan3D p-values, calculated using GeneScan3D.KnockoffGeneration() function.
#' @param fdr  The false discovery rate (FDR) threshold. The default is 0.1.
#' @param gene_id The genes id for N genes considered in the analysis. Usually we consider N=~20,000 protein-coding genes.
#' @return \item{W}{The knockoff statistics for each gene.}
#' @return \item{Qvalue}{The q-values for each gene.}
#' @return \item{gene_sign}{Significant genes with q-values less then the fdr threshold.}
#' @examples
#' library(GeneScan3DKnock)
#'
#'# Load data example
#'data("GeneScan3DKnock.example")
#'
#'result.GeneScan3DKnock=GeneScan3DKnock(M=5,
#'p0=GeneScan3DKnock.example$GeneScan3D.original,
#'p_ko=cbind(GeneScan3DKnock.example$GeneScan3D.ko1,
#'           GeneScan3DKnock.example$GeneScan3D.ko2,
#'           GeneScan3DKnock.example$GeneScan3D.ko3,
#'           GeneScan3DKnock.example$GeneScan3D.ko4,
#'           GeneScan3DKnock.example$GeneScan3D.ko5),fdr = 0.1,
#'           gene_id=GeneScan3DKnock.example$gene.id)
#'result.GeneScan3DKnock$W
#'result.GeneScan3DKnock$Qvalue
#'result.GeneScan3DKnock$gene_sign
#'
#' @export
GeneScan3DKnock<-function(M=5,p0=GeneScan3DKnock.example$GeneScan3D.original,
                                  p_ko=cbind(GeneScan3DKnock.example$GeneScan3D.ko1,
                                             GeneScan3DKnock.example$GeneScan3D.ko2,
                                             GeneScan3DKnock.example$GeneScan3D.ko3,
                                             GeneScan3DKnock.example$GeneScan3D.ko4,
                                             GeneScan3DKnock.example$GeneScan3D.ko5),fdr = 0.1,gene_id=GeneScan3DKnock.example$gene.id){
   
   p=cbind(p0,p_ko)
   
   #calculate knockoff statistics W, kappa, tau for given original p-value and M knockoff p-values
   T=-log10(p)
   
   W=(T[,1]-apply(T[,2:(M+1)],1,median))*(T[,1]>=apply(T[,2:(M+1)],1,max))
   kappa=apply(T,1,which.max)-1 #max T is from original data (0) or knockoff data (1 to 5)
   tau=apply(T,1,max)-apply(T,1,function(x)median(x[-which.max(x)]))
   
   Rej.Bound=10000 
   b=order(tau,decreasing=T)
   c_0=kappa[b]==0  #only calculate q-value for kappa=0
   
   #calculate ratios for top Rej.Bound tau values
   ratio<-c();temp_0<-0
   for(i in 1:length(b)){
      temp_0<-temp_0+c_0[i]
      temp_1<-i-temp_0
      temp_ratio<-(1/M+1/M*temp_1)/max(1,temp_0)
      ratio<-c(ratio,temp_ratio)
      if(i>Rej.Bound){break}
   }
   
   #calculate q value for each gene/window
   qvalue=rep(1,length(tau))
   for(i in 1:length(b)){
      qvalue[b[i]]=min(ratio[i:min(length(b),Rej.Bound)])*c_0[i]+1-c_0[i] #only calculate q-value for kappa=0, q-value for kappa!=0 is 1
      if(i>Rej.Bound){break}
   }
   
   #gene is significant if its q value less or equal than the fdr threshold
   gene_sign=as.character(gene_id[which(round(qvalue,digits=1)<=fdr)])
   
   return(list(W=W,Qvalue=qvalue,gene_sign=gene_sign))
}


### other functions
### cauculated p-values
Get.p<-function(X,result.null.model){
   X<-as.matrix(X)
   mu<-result.null.model$nullglm$fitted.values;Y.res<-result.null.model$Y-mu
   outcome<-result.null.model$out_type
   if(outcome=='D'){
      p<-ScoreTest_SPA(t(X),result.null.model$Y,result.null.model$X,method=c("fastSPA"),minmac=-Inf)$p.value
   }else{
      v<-rep(as.numeric(var(Y.res)),nrow(X))
      p<-pchisq((t(X)%*%Y.res)^2/(apply(X*(v*X),2,sum)-apply(t(X)%*%(v*result.null.model$X0)%*%result.null.model$inv.X0*t(t(result.null.model$X0)%*%as.matrix(v*X)),1,sum)),df=1,lower.tail=F)
   }
   return(as.matrix(p))
}
Get.p.moment<-function(Q,re.Q){ #Q a A*q matrix of test statistics, re.Q a B*q matrix of resampled test statistics
   re.mean<-apply(re.Q,2,mean)
   re.variance<-apply(re.Q,2,var)
   re.kurtosis<-apply((t(re.Q)-re.mean)^4,1,mean)/re.variance^2-3
   re.df<-(re.kurtosis>0)*12/re.kurtosis+(re.kurtosis<=0)*100000
   re.p<-t(pchisq((t(Q)-re.mean)*sqrt(2*re.df)/sqrt(re.variance)+re.df,re.df,lower.tail=F))
   return(re.p)
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
Get.p.SKAT<-function(score,re.score,K,window.matrix,weight,result.null.model){
   
   mu<-result.null.model$nullglm$fitted.values;Y.res<-result.null.model$Y-mu
   X0<-result.null.model$X0;outcome<-result.null.model$out_type
   
   Q<-as.vector(t(score^2)%*%(weight*window.matrix)^2)
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