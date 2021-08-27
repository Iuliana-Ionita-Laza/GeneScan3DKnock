#' @importFrom stats binomial dbeta gaussian glm pcauchy pchisq rbinom sd var median
utils::globalVariables(c('G_Enhancer1_surround','G_Enhancer2_surround','variants_Enhancer1_surround',
                         'variants_Enhancer2_surround','Enhancer1.pos','Enhancer2.pos','create.MK',
'KnockoffGeneration.example','GeneScan3DKnock','GeneScan3DKnock.example',
'G_EnhancerAll','Z_EnhancerAll','p_EnhancerAll',"G_gene_buffer", "Z_gene_buffer", 
"pos_gene_buffer",'n','G_promoter','Z_promoter','G_Enhancer1','Z_Enhancer1','G_Enhancer2','Z_Enhancer2'))


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
                                       MAC.threshold=10,MAF.threshold=0.01,Gsub.id=NULL,result.null.model=result.null.model,M=5){
   
   
   
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
                                p_Enhancer=p_EnhancerAll,window.size=window.size,pos=positions_gene_buffer,
                                MAC.threshold=MAC.threshold,MAF.threshold=MAF.threshold,result.null.model=result.null.model,Gsub.id=row.names(G_gene_buffer))$GeneScan3D.Cauchy.pvalue
   
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
                                                p_Enhancer=p_EnhancerAll,window.size=window.size,pos=positions_gene_buffer,
                                                MAC.threshold=MAC.threshold,MAF.threshold=MAF.threshold,result.null.model=result.null.model,Gsub.id=row.names(G_gene_buffer))$GeneScan3D.Cauchy.pvalue
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
#' @return \item{W.threshold}{Threshold of W statistics. A gene is significant under GeneScan3DKnock if W>W.threshold or equivalently, q-value<fdr threshold. There is no significant gene if W.threshold=Inf.}
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
#'result.GeneScan3DKnock$W.threshold
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
   
   #W statistics threshold
   W.threshold=MK.threshold.byStat(kappa,tau,M=M,fdr=fdr,Rej.Bound=Rej.Bound)
   
   #gene is significant if its q value less or equal than the fdr threshold; OR W>=W.threshold
   gene_sign=as.character(gene_id[which(qvalue<=fdr)])
   return(list(W=W,W.threshold=W.threshold,Qvalue=qvalue,gene_sign=gene_sign))
}


######### Other functions #########
#knockoff filter
MK.threshold.byStat<-function (kappa,tau,M,fdr = 0.1,Rej.Bound=10000){
   b<-order(tau,decreasing=T)
   c_0<-kappa[b]==0
   ratio<-c();temp_0<-0
   for(i in 1:length(b)){
      #if(i==1){temp_0=c_0[i]}
      temp_0<-temp_0+c_0[i]
      temp_1<-i-temp_0
      temp_ratio<-(1/M+1/M*temp_1)/max(1,temp_0)
      ratio<-c(ratio,temp_ratio)
      if(i>Rej.Bound){break}
   }
   ok<-which(ratio<=fdr)
   if(length(ok)>0){
      #ok<-ok[which(ok-ok[1]:(ok[1]+length(ok)-1)<=0)]
      return(tau[b][ok[length(ok)]])
   }else{return(Inf)}
}
