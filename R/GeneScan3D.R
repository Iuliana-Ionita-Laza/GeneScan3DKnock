#' @importFrom stats binomial dbeta gaussian glm pcauchy pchisq rbinom sd var median
utils::globalVariables(c('GeneScan3DKnock.example',"G_gene_buffer", "Z_gene_buffer", "pos_gene_buffer",'n','G_promoter','Z_promoter','G_Enhancer1','Z_Enhancer1','G_Enhancer2','Z_Enhancer2'))


#' Data example for GeneScan3D (A unified framework for gene-based testing with joint analysis of coding and regulatory variation)
#'
#'The dataset contains outcome variable Y, covariate X, genotype data for gene buffer region, promoter and two enhancers, weight matrices for functional annotations and positions of genetic variants in gene buffer region as well as promoter.
#'
#' @name GeneScan3D.example
#' @docType data
#' @keywords data
#' @usage data("GeneScan3D.example")
#' @examples
#' data("GeneScan3D.example")
#'
#'Y=GeneScan3D.example$Y; X=GeneScan3D.example$X
#'
#'G_gene_buffer=GeneScan3D.example$G_gene_buffer
#'G_promoter=GeneScan3D.example$G_promoter
#'G_Enhancer1=GeneScan3D.example$G_Enhancer1
#'G_Enhancer2=GeneScan3D.example$G_Enhancer2
#'
#'Z_gene_buffer=GeneScan3D.example$Z_gene_buffer
#'Z_promoter=GeneScan3D.example$Z_promoter
#'Z_Enhancer1=GeneScan3D.example$Z_Enhancer1
#'Z_Enhancer2=GeneScan3D.example$Z_Enhancer2
#'
#'pos_gene_buffer=GeneScan3D.example$pos_gene_buffer
#'pos_promoter=GeneScan3D.example$pos_promoter
#'n=length(Y)
"GeneScan3D.example"

#'Data example for GeneScan3DKnock (Integration of knockoff statistics for causal gene identification)
#'
#'This example dataset contains the original and five knockoff p-values for N=100 genes. For each gene, there are gene id, original Cauchy3D p-value and five knockoff Cauchy3D p-values. The data can be used in GeneScan3DKnock() function to calculate the knockoff statistics and q-values for each gene. To generate knockoffs and obtain the knockoff p-values, please find the functions in R Package 'KnockoffScreen'.
#'
#' @name GeneScan3DKnock.example
#' @docType data
#' @keywords data
#' @usage data("GeneScan3DKnock.example")
"GeneScan3DKnock.example"



#' The preliminary data management for GeneScan
#'
#' This function does the preliminary data management and fit the model under null hypothesis. The output will be used in the other GeneScan functions.
#'
#' @param Y The outcome variable, an n*1 matrix where n is the total number of observations.
#' @param X An n*d covariates matrix where d is the total number of covariates.
#' @param id The subject id. This is used to match phenotype with genotype. The default is NULL, where the matched phenotype and genotype matrices are assumed.
#' @param out_type Type of outcome variable. Can be either "C" for continuous or "D" for dichotomous. The default is "C".
#' @param B Number of resampling replicates. The default is 1000. A larger value leads to more accurate and stable p-value calculation, but requires more computing time.
#' @return It returns a list used for function GeneScan1D() and GeneScan3D().
#' @examples
#' library(GeneScan3DKnock)
#'
#'# Load data example
#'# Y: outcomes, n by 1 matrix where n is the total number of observations
#'# X: covariates, n by d matrix
#'
#'data("GeneScan3D.example")
#'Y=GeneScan3D.example$Y; X=GeneScan3D.example$X;
#'
#'# Preliminary data management
#'result.prelim=GeneScan.prelim(Y, X, out_type="C", B=1000)
#'
#' @export
GeneScan.prelim<-function(Y, X=NULL, id=NULL, out_type="C", B=1000){
   ##Preliminary
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
   result.prelim<-list(Y=Y,id=id,n=n,X0=X0,nullglm=nullglm,out_type=out_type,re.Y.res=re.Y.res,inv.X0=inv.X0)
   return(result.prelim)
}

#' Conduct gene-based scan test on the gene buffer region.
#'
#' This function conduct gene-based scan test on the gene buffer region using 1D windows under different sizes, do not incorporate any regulatory elements.
#'
#' @param G The genotype matrix in the gene buffer region, which is a n*p matrix where n is the number of subjects and p is the number of genetic variants in the gene buffer region.
#' @param Z A p*q genonet matrix matrix where p is the number of genetic variables and q is the number of functional scores (weights). The default is NULL, which uses the beta(MAF; 1,25) weight.
#' @param window.size The 1-D window sizes in base pairs to scan the gene buffer region. The default is c(1000,5000,10000).
#' @param pos  The positions of genetic variants in the gene buffer region, an p dimensional vector. Each position corresponds to a column in the genotype matrix.
#' @param MAC.threshold Threshold for minor allele count. Variants below MAC.threshold are ultra-rare variants. The default is 5.
#' @param MAF.threshold Threshold for minor allele frequency. Variants below MAF.threshold are rare variants. The default is 1/sqrt(2*n).
#' @param Gsub.id The subject id corresponding to the genotype matrix, an n dimensional vector. The default is NULL, where the matched phenotype and genotype matrices are assumed.
#' @param impute.method Imputation method when there is missing genotype. Can be "random", "fixed" or "bestguess".
#' @param result.prelim The output of function "GeneScan.prelim()".
#' @return \item{GeneScan1D.Cauchy.pvalue}{Cauchy combined p-values under all, common and rare variants for GeneScan1D analysis.}
#' @return \item{M}{Number of 1D scanning windows.}
#' @examples
#' library(GeneScan3DKnock)
#'
#'# Load data example
#'# Y: outcomes, n by 1 matrix where n is the total number of observations
#'# X: covariates, n by d matrix
#'# G_gene_buffer: genotype matrix of gene buffer region, n by p matrix
#'# pos_gene_buffer: positions of genetic variants, p dimentional vector
#'# Z_gene_buffer: functional annotation matrix, p by q matrix
#'
#'data("GeneScan3D.example")
#'Y=GeneScan3D.example$Y; X=GeneScan3D.example$X;
#'G_gene_buffer=GeneScan3D.example$G_gene_buffer;
#'Z_gene_buffer=GeneScan3D.example$Z_gene_buffer;
#'pos_gene_buffer=GeneScan3D.example$pos_gene_buffer;
#'n=length(Y)
#'
#'# Preliminary data management
#'result.prelim=GeneScan.prelim(Y, X, out_type="C", B=1000)
#'
#'# Scan the gene buffer region using 1kb, 5kb and 10kb 1-D windows
#'result.GeneScan1D=GeneScan1D(G=G_gene_buffer,Z=Z_gene_buffer,window.size=c(1000,5000,10000),
#'pos=pos_gene_buffer,MAC.threshold=5,MAF.threshold=1/sqrt(2*n),
#'Gsub.id=NULL, impute.method='fixed',result.prelim=result.prelim)
#'
#' @import SKAT
#' @import Matrix
#' @import WGScan
#' @import SPAtest
#' @import CompQuadForm
#' @export
GeneScan1D<-function(G=G_gene_buffer, Z=Z_gene_buffer, window.size=c(1000,5000,10000), pos=pos_gene_buffer, MAC.threshold=5, MAF.threshold=1/sqrt(2*n), Gsub.id=NULL, impute.method='fixed', result.prelim=result.prelim){

   #load preliminary features
   mu<-result.prelim$nullglm$fitted.values;
   Y.res<-result.prelim$Y-mu
   re.Y.res<-result.prelim$re.Y.res
   X0<-result.prelim$X0
   outcome<-result.prelim$out_type

   #match phenotype id and genotype id
   if(length(Gsub.id)==0){match.index<-match(result.prelim$id,1:nrow(G))}else{
      match.index<-match(result.prelim$id,Gsub.id)
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
      Z.window=Z[index.wondow,]
      Z.window=Matrix(Z.window)

      #if there is only 1 variant in this window, then do not conduct combined test in this window
      #move to the next one
      if(dim(G.window)[2]==1){
         next
      }

      MAF.window<-apply(G.window,2,mean)/2
      MAC.window<-apply(G.window,2,sum)
      weight.beta<-dbeta(MAF.window,1,25)
      weight.matrix<-cbind(MAC.window<MAC.threshold,(MAF.window<MAF.threshold&MAC.window>=MAC.threshold)*weight.beta,(MAF.window>=MAF.threshold)*weight.beta)
      colnames(weight.matrix)<-c('MAC<MAC.threshold','MAF<MAF.threshold&MAC>=MAC.threshold&Beta','MAF>=MAF.thresholdBeta')
      weight.matrix<-Matrix(weight.matrix)

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
         p.burden[,j]<-Get.p.base(X,result.prelim)
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
         p.dispersion[,j]<-Get.p.SKAT(score,re.score,K,window.matrix=as.matrix(rep(1,sum(index.wondow))),weight=(MAC.window>=MAC.threshold)*weight.matrix[,j],result.prelim)
      }

      #Single variant score test for all variants in the window
      p.single<-Get.p(G.window,result.prelim)
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


#' Conduct gene-based scan test on the gene buffer region, adding one promoter and R enhancers.
#'
#' This function conduct gene-based scan test on the gene buffer region, incorporating the regulatory elements, i.e., one promoter and R enhancers.
#'
#' @param G The genotype matrix in the gene buffer region, which is a n*p matrix where n is the number of subjects and p is the number of genetic variants in the gene buffer region.
#' @param Z A p*q genonet matrix matrix where p is the number of genetic variables and q is the number of functional scores (weights). The default is NULL, which uses the beta(MAF; 1,25) weight.
#' @param G.promoter The genotype matrix for promoter region.
#' @param Z.promoter The genonet matrix for promoter region.
#' @param G.EnhancerAll The genotype matrix for R enhancers, combined together by columns.
#' @param Z.EnhancerAll The genonet matrix for R enhancers, combined together by rows.
#' @param p_Enhancer Number of variants in R enhancers, which is a 1*R vector.
#' @param R Number of enhancers.
#' @param window.size The 1-D window sizes in base pairs to scan the gene buffer region. The default is c(1000,5000,10000).
#' @param pos  The positions of genetic variants in the gene buffer region, an p dimensional vector. Each position corresponds to a column in the genotype matrix G.
#' @param pos_promoter The positions of genetic variants in the promoter region. Each position corresponds to a column in the genotype matrix G.promoter.
#' @param MAC.threshold Threshold for minor allele count. Variants below MAC.threshold are ultra-rare variants. The default is 5.
#' @param MAF.threshold Threshold for minor allele frequency. Variants below MAF.threshold are rare variants. The default is 1/sqrt(2*n).
#' @param Gsub.id The subject id corresponding to the genotype matrix, an n dimensional vector. The default is NULL, where the matched phenotype and genotype matrices are assumed.
#' @param impute.method Imputation method when there is missing genotype. Can be "random", "fixed" or "bestguess".
#' @param result.prelim The output of function "GeneScan.prelim()".
#' @return \item{GeneScan3D.Cauchy.pvalue}{Cauchy combined p-values under all, common and rare variants for GeneScan3D analysis.}
#' @return \item{M}{Number of 1D scanning windows.}
#' @return \item{minp}{Minimum p-values under all, common and rare variants for 3D windows.}
#' @return \item{RE_minp}{The regulartory elements in the 3D windows corresponding to the minimum p-values, under all, common and rare variants. 0 represents promoter and a number from 1 to R represents promoter plus r-th enhancer.}
#' @examples
#' library(GeneScan3DKnock)
#'
#'# Load data example
#'
#'data("GeneScan3D.example")
#'Y=GeneScan3D.example$Y; X=GeneScan3D.example$X;
#'G_gene_buffer=GeneScan3D.example$G_gene_buffer; G_promoter=GeneScan3D.example$G_promoter;
#'G_Enhancer1=GeneScan3D.example$G_Enhancer1; G_Enhancer2=GeneScan3D.example$G_Enhancer2;
#'Z_gene_buffer=GeneScan3D.example$Z_gene_buffer; Z_promoter=GeneScan3D.example$Z_promoter;
#'Z_Enhancer1=GeneScan3D.example$Z_Enhancer1; Z_Enhancer2=GeneScan3D.example$Z_Enhancer2;
#'pos_gene_buffer=GeneScan3D.example$pos_gene_buffer;
#'pos_promoter=GeneScan3D.example$pos_promoter;
#'n=length(Y)
#'
#'# Preliminary data management
#'result.prelim=GeneScan.prelim(Y, X, out_type="C", B=1000)
#'
#'# Conduct 3D gene-based scan test on the gene buffer region, adding one promoter and R enhancers
#'result.GeneScan3D=GeneScan3D(G=G_gene_buffer,Z=Z_gene_buffer,
#'G.promoter=G_promoter, Z.promoter=Z_promoter,
#'G.EnhancerAll=cbind(G_Enhancer1,G_Enhancer2),Z.EnhancerAll=rbind(Z_Enhancer1,Z_Enhancer2),
#'R=2,p_Enhancer=c(dim(G_Enhancer1)[2],dim(G_Enhancer2)[2]),window.size=c(1000,5000,10000),
#'pos=pos_gene_buffer,pos_promoter=pos_promoter,MAC.threshold=5,MAF.threshold=1/sqrt(2*n),
#'Gsub.id=NULL,impute.method='fixed',result.prelim=result.prelim)
#'
#' @import SKAT
#' @import Matrix
#' @import WGScan
#' @import SPAtest
#' @import CompQuadForm
#' @export
GeneScan3D<-function(G=G_gene_buffer,Z=Z_gene_buffer,G.promoter=G_promoter,Z.promoter=Z_promoter,G.EnhancerAll=cbind(G_Enhancer1,G_Enhancer2),Z.EnhancerAll=rbind(Z_Enhancer1,Z_Enhancer2), R=2,
                     p_Enhancer=c(dim(G_Enhancer1)[2],dim(G_Enhancer2)[2]),window.size=c(1000,5000,10000),pos=pos_gene_buffer,pos_promoter=pos_promoter,
                     MAC.threshold=5,MAF.threshold=1/sqrt(2*n),Gsub.id=NULL,impute.method='fixed',result.prelim=result.prelim){

   #load preliminary features
   mu<-result.prelim$nullglm$fitted.values;
   Y.res<-result.prelim$Y-mu
   re.Y.res<-result.prelim$re.Y.res
   X0<-result.prelim$X0
   outcome<-result.prelim$out_type

   #match phenotype id and genotype id
   if(length(Gsub.id)==0){match.index<-match(result.prelim$id,1:nrow(G))}else{
      match.index<-match(result.prelim$id,Gsub.id)
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
      pos.tag<-seq(min(pos),max(pos),by=size*1/2)
      pos.tag<-sapply(pos.tag,function(x)pos[which.min(abs(x-pos))])
      window.matrix0_gene_buffer<-cbind(window.matrix0_gene_buffer,sapply(pos.tag,function(x)as.numeric(pos>=x & pos<x+size)))
   }

   window.string_gene_buffer<-apply(window.matrix0_gene_buffer,2,function(x)paste(as.character(x),collapse = ""))
   window.matrix_gene_buffer<-Matrix(window.matrix0_gene_buffer[,match(unique(window.string_gene_buffer),window.string_gene_buffer)])
   #Number of 1-D windows for scanning the gene buffer region
   M_gene_buffer=dim(window.matrix_gene_buffer)[2]

   ###promoter
   #individuals ids are matched with genotype
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
   pos_promoter=pos_promoter[SNP.index.promoter]
   if(!is.null(Z.promoter)){Z.promoter<-Matrix(Z.promoter[SNP.index.promoter,])}

   p_promoter=dim(G.promoter)[2]

   #Set overlapped variants of promoter in window.matrix_gene_buffer matrix to 0
   window.matrix_gene_buffer[pos %in% pos_promoter,]=matrix(0,nrow=dim(as.matrix(window.matrix_gene_buffer[pos %in% pos_promoter,]))[1],
                                                            ncol=dim(as.matrix(window.matrix_gene_buffer[pos %in% pos_promoter,]))[2])

   if (p_promoter==0){
      warning('0 variant in promoter')
      GeneScan3D.Cauchy.promoter=c()
   }else{
      GeneScan3D.Cauchy.promoter=matrix(NA,nrow=M_gene_buffer,ncol=3)

      ######### GeneScan3D-promoter #########
      for (m in 1:M_gene_buffer){
         #Create index for each 3D window, using 3D window adding promoter
         index.wondow<-(window.matrix_gene_buffer[,m]==1)
         G.window.promoter=cbind(G[,index.wondow],G.promoter)
         Z.window.promoter=rbind(Matrix(Z[index.wondow,]),Z.promoter)
         G.window.promoter=Matrix(G.window.promoter)
         Z.window.promoter=Matrix(Z.window.promoter)

         #if there are less than 5 variants in this window, then do not conduct combined test in this window, move to the next one
         if(dim(G.window.promoter)[2]<5){
            next
         }

         MAF.window.promoter<-apply(G.window.promoter,2,mean)/2
         MAC.window.promoter<-apply(G.window.promoter,2,sum)

         weight.beta<-dbeta(MAF.window.promoter,1,25) #Beta(MAF; 1,25) weights
         weight.matrix<-cbind(MAC.window.promoter<MAC.threshold,(MAF.window.promoter<MAF.threshold&MAC.window.promoter>=MAC.threshold)*weight.beta,(MAF.window.promoter>=MAF.threshold)*weight.beta)
         colnames(weight.matrix)<-c('MAC<MAC.threshold','MAF<MAF.threshold&MAC>=MAC.threshold&Beta','MAF>=MAF.thresholdBeta')

         ##adding additional functional scores
         if (!is.null(Z.window.promoter)){
            colnames(Z.window.promoter)<-paste0('MAF<MAF.threshold&MAC>=MAC.threshold&',1:ncol(Z.window.promoter))
            weight.matrix<-cbind(weight.matrix,(MAF.window.promoter<MAF.threshold&MAC.window.promoter>=MAC.threshold)*Z.window.promoter)
         }
         weight.matrix<-Matrix(weight.matrix)

         #Burden test
         p.burden.promoter<-matrix(NA,1,ncol(weight.matrix))
         for (j in 1:ncol(weight.matrix)){
            temp.window.matrix<-weight.matrix[,j]
            X<-as.matrix(G.window.promoter%*%temp.window.matrix)
            p.burden.promoter[,j]<-Get.p.base(X,result.prelim)
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
            p.dispersion.promoter[,j]<-Get.p.SKAT(score,re.score,K,window.matrix=as.matrix(rep(1,dim(G.window.promoter)[2])),weight=(MAC.window.promoter>=MAC.threshold)*weight.matrix[,j],result.prelim)
         }

         #Single variant score test for all variants in the window
         p.single.promoter<-Get.p(G.window.promoter,result.prelim)
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

         GeneScan3D.Cauchy.promoter[m,]=c(p.Cauchy.promoter,p.Cauchy.common.promoter,p.Cauchy.rare.promoter)
      }
   }

   ###R enhancers
   Enhancer_ind=rep(TRUE,R)
   GeneScan3D.Cauchy.EnhancerAll=c()
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
      Z.Enhancer=Matrix(Z.Enhancer)

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

      GeneScan3D.Cauchy.Enhancer=matrix(NA,nrow=M_gene_buffer,ncol=3)

      for (m in 1:M_gene_buffer){
         #Create index for each 3D window, using 3D window adding promoter
         index.wondow<-(window.matrix_gene_buffer[,m]==1)
         #Create index for each 3D window, using 1D window adding promoter and enhancer
         if (p_promoter==0){
            G.window.Enhancer=cbind(G[,index.wondow],G.Enhancer)
            Z.window.Enhancer=rbind(Matrix(Z[index.wondow,]),Z.Enhancer)
         }else{
            G.window.Enhancer=cbind(G[,index.wondow],G.promoter,G.Enhancer)
            Z.window.Enhancer=rbind(Matrix(Z[index.wondow,]),Z.promoter,Z.Enhancer)
         }
         G.window.Enhancer=Matrix(G.window.Enhancer)
         Z.window.Enhancer=Matrix(Z.window.Enhancer)

         #if there are less than 5 variants in this window, then do not conduct combined test in this window, move to the next one
         if(dim(G.window.Enhancer)[2]<5){
            Enhancer_ind[r]=FALSE
            next
         }

         MAF.window.Enhancer<-apply(G.window.Enhancer,2,mean)/2
         MAC.window.Enhancer<-apply(G.window.Enhancer,2,sum)   #MAF=MAC/5602/2

         weight.beta<-dbeta(MAF.window.Enhancer,1,25) #Beta(MAF; 1,25) weights
         weight.matrix<-cbind(MAC.window.Enhancer<MAC.threshold,(MAF.window.Enhancer<MAF.threshold&MAC.window.Enhancer>=MAC.threshold)*weight.beta,(MAF.window.Enhancer>=MAF.threshold)*weight.beta)
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
            p.burden.Enhancer[,j]<-Get.p.base(X,result.prelim)
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
            p.dispersion.Enhancer[,j]<-Get.p.SKAT(score,re.score,K,window.matrix=as.matrix(rep(1,dim(G.window.Enhancer)[2])),weight=(MAC.window.Enhancer>=MAC.threshold)*weight.matrix[,j],result.prelim)
         }

         #Single variant score test for all variants in the window
         p.single.Enhancer<-Get.p(G.window.Enhancer,result.prelim)

         p.individual1.Enhancer<-Get.cauchy.scan(p.single.Enhancer,as.matrix((MAC.window.Enhancer>=MAC.threshold & MAF.window.Enhancer<MAF.threshold))) #rare variants
         p.individual2.Enhancer<-Get.cauchy.scan(p.single.Enhancer,as.matrix((MAF.window.Enhancer>=MAF.threshold))) #common and low frequency variants
         p.individual.Enhancer<-cbind(p.burden.Enhancer ,p.dispersion.Enhancer,p.individual1.Enhancer,p.individual2.Enhancer);
         colnames(p.individual.Enhancer)<-c(paste0('burden_',colnames(weight.matrix)),paste0('dispersion_',colnames(weight.matrix)),'singleCauchy_MAF<MAF.threshold&MAC>=MAC.threshold','singleCauchy_MAF>=MAF.threshold')

         #aggregated Cauchy association test
         p.Cauchy.Enhancer<-as.matrix(apply(p.individual.Enhancer,1,Get.cauchy))
         test.common<-grep('MAF>=MAF.threshold',colnames(p.individual.Enhancer))
         p.Cauchy.common.Enhancer<-as.matrix(apply(p.individual.Enhancer[,test.common,drop=FALSE],1,Get.cauchy))
         p.Cauchy.rare.Enhancer<-as.matrix(apply(p.individual.Enhancer[,-test.common,drop=FALSE],1,Get.cauchy))

         GeneScan3D.Cauchy.Enhancer[m,]=c(p.Cauchy.Enhancer,p.Cauchy.common.Enhancer,p.Cauchy.rare.Enhancer)
      }
      GeneScan3D.Cauchy.EnhancerAll=rbind(GeneScan3D.Cauchy.EnhancerAll,GeneScan3D.Cauchy.Enhancer)
   } #end of the loop of enhancers

   GeneScan3D.Cauchy.RE=rbind(GeneScan3D.Cauchy.promoter,GeneScan3D.Cauchy.EnhancerAll)
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



#' Integration of knockoff statistics for causal gene identification
#'
#' This function calculates the knockoff statistics and q-values after proving the original and knockoff p-values for each gene (or window). 
#'
#' @param M Number of multiple knockoffs. We use M=5 in our analysis.
#' @param p0 A N-dimensional vector of the original p-values for N genes considered in the analysis. The p-values can be obtained in GeneScan3D() function or other analysis. 
#' @param p_ko A N*M matrix of M knockoff p-values for N genes considered in the analysis. The knockoff p-values can be obtained in R package 'KnockoffScreen'.
#' @param fdr  The false discovery rate (FDR) threshold. The default is 0.1.
#' @param gene_id The genes id for N genes considered in the analysis, which can also be the windows id or an indicator vector from 1 to N. 
#' @return \item{W}{The knockoff statistics for N genes.}
#' @return \item{Qvalue}{The Q-values for N genes.}
#' @return \item{gene_sign}{Significant genes obtained in the knockoff test with Q-values less then the fdr threshold.}
#' @examples
#' library(GeneScan3DKnock)
#'
#'# Load data example
#'data("GeneScan3DKnock.example")
#'
#'result.GeneScan3DKnock=GeneScan3DKnock(M=5,p0=GeneScan3DKnock.example$Cauchy3D.all.original,
#'p_ko=cbind(GeneScan3DKnock.example$Cauchy3D.all.ko1,GeneScan3DKnock.example$Cauchy3D.all.ko2,
#' GeneScan3DKnock.example$Cauchy3D.all.ko3,GeneScan3DKnock.example$Cauchy3D.all.ko4,
#'GeneScan3DKnock.example$Cauchy3D.all.ko5),fdr = 0.1,gene_id=GeneScan3DKnock.example$gene.id)
#'
#'#Obtain knockoff statistics, q-values and the significant genes.
#'W=result.GeneScan3DKnock$W
#'Qvalue=result.GeneScan3DKnock$Qvalue
#'gene_sign=result.GeneScan3DKnock$gene_sign
#'
#' @export
GeneScan3DKnock<-function(M=5,p0=GeneScan3DKnock.example$Cauchy3D.all.original,
                          p_ko=cbind(GeneScan3DKnock.example$Cauchy3D.all.ko1,GeneScan3DKnock.example$Cauchy3D.all.ko2,
                                     GeneScan3DKnock.example$Cauchy3D.all.ko3,GeneScan3DKnock.example$Cauchy3D.all.ko4,
                                     GeneScan3DKnock.example$Cauchy3D.all.ko5),fdr = 0.1,gene_id=GeneScan3DKnock.example$gene.id){
   
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

### cauculated p-values
Get.p<-function(X,result.prelim){
   X<-as.matrix(X)
   mu<-result.prelim$nullglm$fitted.values;Y.res<-result.prelim$Y-mu
   outcome<-result.prelim$out_type
   if(outcome=='D'){
      p<-ScoreTest_SPA(t(X),result.prelim$Y,result.prelim$X,method=c("fastSPA"),minmac=-Inf)$p.value
   }else{
      v<-rep(as.numeric(var(Y.res)),nrow(X))
      p<-pchisq((t(X)%*%Y.res)^2/(apply(X*(v*X),2,sum)-apply(t(X)%*%(v*result.prelim$X0)%*%result.prelim$inv.X0*t(t(result.prelim$X0)%*%as.matrix(v*X)),1,sum)),df=1,lower.tail=F)
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
Get.p.base<-function(X,result.prelim){
   X<-as.matrix(X)
   mu<-result.prelim$nullglm$fitted.values;Y.res<-result.prelim$Y-mu
   outcome<-result.prelim$out_type
   if(outcome=='D'){v<-mu*(1-mu)}else{v<-rep(as.numeric(var(Y.res)),nrow(X))}
   p<-pchisq((t(X)%*%Y.res)^2/(apply(X*(v*X),2,sum)-apply(t(X)%*%(v*result.prelim$X0)%*%result.prelim$inv.X0*t(t(result.prelim$X0)%*%as.matrix(v*X)),1,sum)),df=1,lower.tail=F)
   p[is.na(p)]<-NA
   return(p)
}
Get.p.SKAT<-function(score,re.score,K,window.matrix,weight,result.prelim){
   
   mu<-result.prelim$nullglm$fitted.values;Y.res<-result.prelim$Y-mu
   X0<-result.prelim$X0;outcome<-result.prelim$out_type
   
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
