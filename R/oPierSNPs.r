#' Function to prioritise genes given a list of seed SNPs together with the significance level (e.g. GWAS reported p-values)
#'
#' \code{oPierSNPs} is supposed to prioritise genes given a list of seed SNPs together with the significance level. To prioritise genes, it first defines and scores seed genes: nearby genes, eQTL genes and Hi-C genes. With seed genes and their scores, it then uses Random Walk with Restart (RWR) to calculate the affinity score of all nodes in the input graph to the seed genes. The priority score is the affinity score. It returns an object of class "pNode".
#'
#' @param data a named input vector containing the sinificance level for nodes (dbSNP). For this named vector, the element names are dbSNP ID (or in the format such as 'chr16:28525386'), the element values for the significance level (measured as p-value or fdr). Alternatively, it can be a matrix or data frame with two columns: 1st column for dbSNP, 2nd column for the significance level
#' @param include.LD additional SNPs in LD with Lead SNPs are also included. By default, it is 'NA' to disable this option. Otherwise, LD SNPs will be included based on one or more of 26 populations and 5 super populations from 1000 Genomics Project data (phase 3). The population can be one of 5 super populations ("AFR", "AMR", "EAS", "EUR", "SAS"), or one of 26 populations ("ACB", "ASW", "BEB", "CDX", "CEU", "CHB", "CHS", "CLM", "ESN", "FIN", "GBR", "GIH", "GWD", "IBS", "ITU", "JPT", "KHV", "LWK", "MSL", "MXL", "PEL", "PJL", "PUR", "STU", "TSI", "YRI"). Explanations for population code can be found at \url{http://www.1000genomes.org/faq/which-populations-are-part-your-study}
#' @param LD.customised a user-input matrix or data frame with 3 columns: 1st column for Lead SNPs, 2nd column for LD SNPs, and 3rd for LD r2 value. It is designed to allow the user analysing their pre-calculated LD info. This customisation (if provided) has the high priority over built-in LD SNPs
#' @param LD.r2 the LD r2 value. By default, it is 0.8, meaning that SNPs in LD (r2>=0.8) with input SNPs will be considered as LD SNPs. It can be any value from 0.8 to 1
#' @param significance.threshold the given significance threshold. By default, it is set to NULL, meaning there is no constraint on the significance level when transforming the significance level of SNPs into scores. If given, those SNPs below this are considered significant and thus scored positively. Instead, those above this are considered insigificant and thus receive no score
#' @param score.cap the maximum score being capped. By default, it is set to 10. If NULL, no capping is applied
#' @param distance.max the maximum distance between genes and SNPs. Only those genes no far way from this distance will be considered as seed genes. This parameter will influence the distance-component weights calculated for nearby SNPs per gene
#' @param decay.kernel a character specifying a decay kernel function. It can be one of 'slow' for slow decay, 'linear' for linear decay, and 'rapid' for rapid decay. If no distance weight is used, please select 'constant'
#' @param decay.exponent an integer specifying a decay exponent. By default, it sets to 2
#' @param GR.SNP the genomic regions of SNPs. By default, it is 'dbSNP_GWAS', that is, SNPs from dbSNP (version 146) restricted to GWAS SNPs and their LD SNPs (hg19). It can be 'dbSNP_Common', that is, Common SNPs from dbSNP (version 146) plus GWAS SNPs and their LD SNPs (hg19). Alternatively, the user can specify the customised input. To do so, first save your RData file (containing an GR object) into your local computer, and make sure the GR object content names refer to dbSNP IDs. Then, tell "GR.SNP" with your RData file name (with or without extension), plus specify your file RData path in "placeholder". Note: you can also load your customised GR object directly
#' @param GR.Gene the genomic regions of genes. By default, it is 'UCSC_knownGene', that is, UCSC known genes (together with genomic locations) based on human genome assembly hg19. It can be 'UCSC_knownCanonical', that is, UCSC known canonical genes (together with genomic locations) based on human genome assembly hg19. Alternatively, the user can specify the customised input. To do so, first save your RData file (containing an GR object) into your local computer, and make sure the GR object content names refer to Gene Symbols. Then, tell "GR.Gene" with your RData file name (with or without extension), plus specify your file RData path in "placeholder". Note: you can also load your customised GR object directly
#' @param include.TAD TAD boundary regions are also included. By default, it is 'none' to disable this option. Otherwise, inclusion of a TAD dataset to pre-filter SNP-nGene pairs (i.e. only those within a TAD region will be kept). TAD datasets can be one of "GM12878"  (lymphoblast), "IMR90" (fibroblast), "MSC" (mesenchymal stem cell) ,"TRO" (trophoblasts-like cell), "H1" (embryonic stem cell), "MES" (mesendoderm) and "NPC" (neural progenitor cell). Explanations can be found at \doi{10.1016/j.celrep.2016.10.061}
#' @param include.QTL the eQTL supported currently. By default, it is 'NA' to disable this option. Pre-built eQTL datasets are detailed in \code{\link{oDefineQTL}}
#' @param QTL.customised a user-input matrix or data frame with 4 columns: 1st column for SNPs/eQTLs, 2nd column for Genes, 3rd for eQTL mapping significance level (p-values or FDR), and 4th for contexts (required even though only one context is input). Alternatively, it can be a file containing these 4 columns. It is designed to allow the user analysing their eQTL data. This customisation (if provided) will populate built-in eQTL data
#' @param include.RGB genes linked to input SNPs are also included. By default, it is 'NA' to disable this option. Otherwise, those genes linked to SNPs will be included according to Promoter Capture HiC (PCHiC) datasets. Pre-built HiC datasets are detailed in \code{\link{oDefineRGB}}
#' @param cdf.function a character specifying a Cumulative Distribution Function (cdf). It can be one of 'exponential' based on exponential cdf, 'empirical' for empirical cdf
#' @param relative.importance a vector specifying the relative importance of nearby genes, eQTL genes and HiC genes. By default, it sets c(1/3, 1/3, 1/3)
#' @param scoring.scheme the method used to calculate seed gene scores under a set of SNPs. It can be one of "sum" for adding up, "max" for the maximum, and "sequential" for the sequential weighting. The sequential weighting is done via: \eqn{\sum_{i=1}{\frac{R_{i}}{i}}}, where \eqn{R_{i}} is the \eqn{i^{th}} rank (in a descreasing order)
#' @param network the built-in network. Currently two sources of network information are supported: the STRING database (version 10) and the Pathway Commons database (version 7). STRING is a meta-integration of undirect interactions from the functional aspect, while Pathways Commons mainly contains both undirect and direct interactions from the physical/pathway aspect. Both have scores to control the confidence of interactions. Therefore, the user can choose the different quality of the interactions. In STRING, "STRING_highest" indicates interactions with highest confidence (confidence scores>=900), "STRING_high" for interactions with high confidence (confidence scores>=700), "STRING_medium" for interactions with medium confidence (confidence scores>=400), and "STRING_low" for interactions with low confidence (confidence scores>=150). For undirect/physical interactions from Pathways Commons, "PCommonsUN_high" indicates undirect interactions with high confidence (supported with the PubMed references plus at least 2 different sources), "PCommonsUN_medium" for undirect interactions with medium confidence (supported with the PubMed references). For direct (pathway-merged) interactions from Pathways Commons, "PCommonsDN_high" indicates direct interactions with high confidence (supported with the PubMed references plus at least 2 different sources), and "PCommonsUN_medium" for direct interactions with medium confidence (supported with the PubMed references). In addition to pooled version of pathways from all data sources, the user can also choose the pathway-merged network from individual sources, that is, "PCommonsDN_Reactome" for those from Reactome, "PCommonsDN_KEGG" for those from KEGG, "PCommonsDN_HumanCyc" for those from HumanCyc, "PCommonsDN_PID" for those froom PID, "PCommonsDN_PANTHER" for those from PANTHER, "PCommonsDN_ReconX" for those from ReconX, "PCommonsDN_TRANSFAC" for those from TRANSFAC, "PCommonsDN_PhosphoSite" for those from PhosphoSite, and "PCommonsDN_CTD" for those from CTD. For direct (pathway-merged) interactions sourced from KEGG, it can be 'KEGG' for all, 'KEGG_metabolism' for pathways grouped into 'Metabolism', 'KEGG_genetic' for 'Genetic Information Processing' pathways, 'KEGG_environmental' for 'Environmental Information Processing' pathways, 'KEGG_cellular' for 'Cellular Processes' pathways, 'KEGG_organismal' for 'Organismal Systems' pathways, and 'KEGG_disease' for 'Human Diseases' pathways. 'REACTOME' for protein-protein interactions derived from Reactome pathways
#' @param STRING.only the further restriction of STRING by interaction type. If NA, no such restriction. Otherwide, it can be one or more of "neighborhood_score","fusion_score","cooccurence_score","coexpression_score","experimental_score","database_score","textmining_score". Useful options are c("experimental_score","database_score"): only experimental data (extracted from BIND, DIP, GRID, HPRD, IntAct, MINT, and PID) and curated data (extracted from Biocarta, BioCyc, GO, KEGG, and Reactome) are used
#' @param weighted logical to indicate whether edge weights should be considered. By default, it sets to false. If true, it only works for the network from the STRING database
#' @param network.customised an object of class "igraph". By default, it is NULL. It is designed to allow the user analysing their customised network data that are not listed in the above argument 'network'. This customisation (if provided) has the high priority over built-in network. If the user provides the "igraph" object with the "weight" edge attribute, RWR will assume to walk on the weighted network
#' @param seeds.inclusive logical to indicate whether non-network seed genes are included for prioritisation. If TRUE (by default), these genes will be added to the netowrk
#' @param normalise the way to normalise the adjacency matrix of the input graph. It can be 'laplacian' for laplacian normalisation, 'row' for row-wise normalisation, 'column' for column-wise normalisation, or 'none'
#' @param restart the restart probability used for Random Walk with Restart (RWR). The restart probability takes the value from 0 to 1, controlling the range from the starting nodes/seeds that the walker will explore. The higher the value, the more likely the walker is to visit the nodes centered on the starting nodes. At the extreme when the restart probability is zero, the walker moves freely to the neighbors at each step without restarting from seeds, i.e., following a random walk (RW)
#' @param normalise.affinity.matrix the way to normalise the output affinity matrix. It can be 'none' for no normalisation, 'quantile' for quantile normalisation to ensure that columns (if multiple) of the output affinity matrix have the same quantiles
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param placeholder the characters to tell the location of built-in RData files. See \code{\link{oRDS}} for details
#' @param guid a valid (5-character) Global Unique IDentifier for an OSF project. See \code{\link{oRDS}} for details
#' @return
#' an object of class "pNode", a list with following components:
#' \itemize{
#'  \item{\code{priority}: a matrix of nNode X 6 containing node priority information, where nNode is the number of nodes in the input graph, and the 6 columns are "name" (node names), "node" (1 for network genes, 0 for non-network seed genes), "seed" (1 for seeds, 0 for non-seeds), "weight" (weight values),  "priority" (the priority scores that are rescaled to the range [0,1]), "rank" (ranks of the priority scores), "description" (node description)}
#'  \item{\code{g}: an input "igraph" object}
#'  \item{\code{SNP}: a data frame of nSNP X 4 containing input SNPs and/or LD SNPs info, where nSNP is the number of input SNPs and/or LD SNPs, and the 4 columns are "SNP" (dbSNP), "Score" (the SNP score), "Pval" (the SNP p-value), "Flag" (indicative of Lead SNPs or LD SNPs)}
#'  \item{\code{Gene2SNP}: a data frame of nPair X 3 containing Gene-SNP pair info, where nPair is the number of Gene-SNP pairs, and the 3 columns are "Gene" (seed genes), "SNP" (dbSNP), "Score" (an SNP's genetic influential score on a seed gene)}
#'  \item{\code{nGenes}: if not NULL, it is a data frame containing nGene-SNP pair info}
#'  \item{\code{eGenes}: if not NULL, it is a data frame containing eGene-SNP pair info per context}
#'  \item{\code{cGenes}: if not NULL, it is a data frame containing cGene-SNP pair info per context}
#' }
#' @note The prioritisation procedure (from SNPs to target genes) consists of following steps:
#' \itemize{
#' \item{i) \code{\link{oSNPscores}} used to calculate the SNP score.}
#' \item{ii) \code{\link{oSNP2nGenes}} used to define and score the nearby genes.}
#' \item{iii) \code{\link{oSNP2eGenes}} used to define and score the eQTL genes.}
#' \item{iv) \code{\link{oSNP2cGenes}} used to define and score the HiC genes.}
#' \item{v) define seed genes as the nearby genes in ii) and the eQTL genes in iii) and the HiC genes in iv), which are then scored in an integrative manner.}
#' \item{vi) \code{\link{oPierGenes}} used to prioritise genes using an input graph and a list of seed genes and their scores from v). The priority score is the affinity score estimated by Random Walk with Restart (RWR), measured as the affinity of all nodes in the graph to the seeds.}
#' }
#' @export
#' @seealso \code{\link{oSNPscores}}, \code{\link{oSNP2nGenes}}, \code{\link{oSNP2eGenes}}, \code{\link{oSNP2cGenes}}, \code{\link{oSparseMatrix}}, \code{\link{oSM2DF}}, \code{\link{oPierGenes}}, \code{\link{oSM2DF}}, \code{\link{oDefineQTL}}, \code{\link{oDefineRGB}}
#' @include oPierSNPs.r
#' @examples
#'
#' \dontrun{
#' # a) provide the SNPs with the significance info
#' ImmunoBase <- oRDS('ImmunoBase', placeholder=placeholder)
#' gr <- ImmunoBase$AS$variants
#' AS <- as.data.frame(gr) %>% select(Variant, Pvalue)
#'
#' # b) perform priority analysis
#' pNode <- oPierSNPs(data=AS, include.QTL="eQTL_xMEdb_Monocyte", include.RGB='PCHiC_PMID27863249_Combined', network="PCommonsUN_medium", restart=0.7, placeholder=placeholder)
#'
#' # c) save to the file called 'SNPs_priority.txt'
#' pNode$priority %>% write.delim("SNPs_priority.txt", delim="\t")
#' 
#' # d) manhattan plot
#' mp <- oPierManhattan(pNode, top=20, top.label.size=1.5, y.scale="sqrt", placeholder=placeholder)
#' #pdf(file="Gene_manhattan.pdf", height=6, width=12, compress=TRUE)
#' print(mp)
#' #dev.off()
#' }

oPierSNPs <- function(data, include.LD=NA, LD.customised=NULL, LD.r2=0.8, significance.threshold=5e-8, score.cap=100, distance.max=20000, decay.kernel=c("constant","slow","linear","rapid"), decay.exponent=2, GR.SNP=c("dbSNP_Common","dbSNP_GWAS","dbSNP_Single"), GR.Gene=c("UCSC_knownGene","UCSC_knownCanonical"), include.TAD=c("none","GM12878","IMR90","MSC","TRO","H1","MES","NPC"), include.QTL=NA, QTL.customised=NULL, include.RGB=NA, cdf.function=c("empirical","exponential"), relative.importance=c(1/3,1/3,1/3), scoring.scheme=c("max","sum","sequential"), network=c("STRING_highest","STRING_high","STRING_medium","STRING_low","PCommonsUN_high","PCommonsUN_medium","PCommonsDN_high","PCommonsDN_medium","PCommonsDN_Reactome","PCommonsDN_KEGG","PCommonsDN_HumanCyc","PCommonsDN_PID","PCommonsDN_PANTHER","PCommonsDN_ReconX","PCommonsDN_TRANSFAC","PCommonsDN_PhosphoSite","PCommonsDN_CTD", "KEGG","KEGG_metabolism","KEGG_genetic","KEGG_environmental","KEGG_cellular","KEGG_organismal","KEGG_disease","REACTOME"), STRING.only=c(NA,"neighborhood_score","fusion_score","cooccurence_score","coexpression_score","experimental_score","database_score","textmining_score")[1], weighted=FALSE, network.customised=NULL, seeds.inclusive=TRUE, normalise=c("laplacian","row","column","none"), restart=0.7, normalise.affinity.matrix=c("none","quantile"), verbose=TRUE, placeholder=NULL, guid=NULL)
{

    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=TRUE)
        message("", appendLF=TRUE)
    }
    ####################################################################################
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    decay.kernel <- match.arg(decay.kernel)
    include.TAD <- match.arg(include.TAD)
    cdf.function <- match.arg(cdf.function)
    scoring.scheme <- match.arg(scoring.scheme)
    network <- match.arg(network)
    normalise <- match.arg(normalise)
    normalise.affinity.matrix <- match.arg(normalise.affinity.matrix)
    
    ####################################################################################
    
    if(verbose){
        message(sprintf("\n#######################################################"))
        message(sprintf("'oSNPscores' is being called to score SNPs (%s):", as.character(Sys.time())), appendLF=TRUE)
        message(sprintf("#######################################################"))
    }
    
	df_SNP <- oSNPscores(data=data, include.LD=include.LD, LD.customised=LD.customised, LD.r2=LD.r2, significance.threshold=significance.threshold, score.cap=score.cap, verbose=verbose, placeholder=placeholder, guid=guid)
	
	if(verbose){
        message(sprintf("#######################################################"))
        message(sprintf("'oSNPscores' has been finished (%s)!", as.character(Sys.time())), appendLF=TRUE)
        message(sprintf("#######################################################\n"))
    }
    
    ####################################################################################
    
    if(verbose){
       message(sprintf("\n#######################################################"))
        message(sprintf("'oSNP2nGenes' is being called to define nearby genes (%s):", as.character(Sys.time())), appendLF=TRUE)
        message(sprintf("#######################################################"))
    }
    
    
    if(relative.importance[1] != 0){
		df_nGenes <- oSNP2nGenes(data=df_SNP$SNP, distance.max=distance.max, decay.kernel=decay.kernel, decay.exponent=decay.exponent, GR.SNP=GR.SNP, GR.Gene=GR.Gene, include.TAD=include.TAD, verbose=verbose, placeholder=placeholder, guid=guid)
		if(include.TAD!='none'){
			TAD <- NULL
			df_nGenes <- base::subset(df_nGenes, TAD!='Excluded')
		}
	}else{
		df_nGenes <- NULL
		if(verbose){
			message(sprintf("No nearby genes are defined"), appendLF=TRUE)
		}
	}
	
	if(verbose){
        message(sprintf("#######################################################"))
        message(sprintf("'oSNP2nGenes' has been finished (%s)!", as.character(Sys.time())), appendLF=TRUE)
        message(sprintf("#######################################################\n"))
    }
    
    ####################################################################################
    
    if(verbose){
        message(sprintf("\n#######################################################"))
        message(sprintf("'oSNP2eGenes' is being called to define eQTL-containing genes (%s):", as.character(Sys.time())), appendLF=TRUE)
        message(sprintf("#######################################################"))
    }
    
    if(relative.importance[2] != 0){
		df_eGenes <- oSNP2eGenes(data=df_SNP$SNP, include.QTL=include.QTL, QTL.customised=QTL.customised, GR.Gene=GR.Gene, cdf.function=cdf.function, verbose=verbose, placeholder=placeholder, guid=guid)
	}else{
		df_eGenes <- NULL
		
		if(verbose){
			message(sprintf("No eQTL genes are defined"), appendLF=TRUE)
		}
	}
	
	if(verbose){
       message(sprintf("#######################################################"))
        message(sprintf("'oSNP2eGenes' has been finished (%s)!", as.character(Sys.time())), appendLF=TRUE)
        message(sprintf("#######################################################\n"))
    }
    ####################################################################################
    
    ####################################################################################
    
    if(verbose){
        message(sprintf("\n#######################################################"))
        message(sprintf("'oSNP2cGenes' is being called to define HiC-captured genes (%s):", as.character(Sys.time())), appendLF=TRUE)
        message(sprintf("#######################################################"))
    }
    
    if(relative.importance[3] != 0){
		df_cGenes <- oSNP2cGenes(data=df_SNP$SNP, entity="SNP", include.RGB=include.RGB, GR.SNP=GR.SNP, GR.Gene=GR.Gene, cdf.function=cdf.function, verbose=verbose, placeholder=placeholder, guid=guid)
	}else{
		df_cGenes <- NULL
		
		if(verbose){
			now <- Sys.time()
			message(sprintf("No HiC genes are defined"), appendLF=TRUE)
		}
	}
	
	if(verbose){
        message(sprintf("#######################################################"))
        message(sprintf("'oSNP2cGenes' has been finished (%s)!", as.character(Sys.time())), appendLF=TRUE)
        message(sprintf("#######################################################\n"))
    }
    ####################################################################################
    
    if(is.null(df_nGenes) & is.null(df_eGenes) & is.null(df_cGenes)){
    	G2S <- NULL
    }else{
    
		## df_SNP df_nGenes df_eGenes df_cGenes
		allGenes <- sort(base::Reduce(base::union, list(df_nGenes$Gene,df_eGenes$Gene,df_cGenes$Gene)))
		allSNPs <- sort(df_SNP$SNP)
	
		## sparse matrix of nGenes X SNPs
		G2S_n <- oSparseMatrix(df_nGenes[,c("Gene","SNP","Weight")], rows=allGenes, columns=allSNPs, verbose=FALSE)
		## sparse matrix of eGenes X SNPs
		G2S_e <- oSparseMatrix(df_eGenes[,c("Gene","SNP","Weight")], rows=allGenes, columns=allSNPs, verbose=FALSE)
		## sparse matrix of cGenes X SNPs
		G2S_c <- oSparseMatrix(df_cGenes[,c("Gene","SNP","Weight")], rows=allGenes, columns=allSNPs, verbose=FALSE)
	
		## combine both sparse matrix
		### wG2S_n
		if(is.null(G2S_n)){
			wG2S_n <- 0
		}else{
			wG2S_n <- G2S_n * relative.importance[1]
		}
		### wG2S_e
		if(is.null(G2S_e)){
			wG2S_e <- 0
		}else{
			wG2S_e <- G2S_e * relative.importance[2]
		}
		### wG2S_c
		if(is.null(G2S_c)){
			wG2S_c <- 0
		}else{
			wG2S_c <- G2S_c * relative.importance[3]
		}
		
		if(is.null(G2S_n) & is.null(G2S_e) & is.null(G2S_c)){
			G2S <- NULL
		}else{
			G2S <- wG2S_n + wG2S_e + wG2S_c
		}
    
    }
    
    #######################
    ## if NULL, return NULL
    if(is.null(G2S)){
    	return(NULL)
    }
    #######################
    
    ## consider SNP scores
    ind <- match(colnames(G2S), df_SNP$SNP)
    ########
    df_SNP <- df_SNP[ind,]
    ########
    SNP_score <- df_SNP$Score
    names(SNP_score) <- colnames(G2S)
    ## convert into matrix
    mat_SNP_score <- matrix(rep(SNP_score,each=nrow(G2S)), nrow=nrow(G2S))
    
    ## calculate genetic influence score for a gene-SNP pair
    G2S_score <- G2S * mat_SNP_score
    
    ## Gene2SNP
    Gene2SNP <- oSM2DF(data=G2S_score, verbose=FALSE)
    colnames(Gene2SNP) <- c('Gene','SNP','Score')
    Gene2SNP <- Gene2SNP[order(Gene2SNP$Gene,-Gene2SNP$Score,decreasing=FALSE),]
	
	ls_gene <- split(x=Gene2SNP$Score, f=Gene2SNP$Gene)
    ## calculate genetic influence score under a set of SNPs for each seed gene
    if(scoring.scheme=='max'){
		seeds.genes <- sapply(ls_gene, max)
		
    }else if(scoring.scheme=='sum'){
		seeds.genes <- sapply(ls_gene, sum)
		
    }else if(scoring.scheme=='sequential'){
		seeds.genes <- sapply(ls_gene, function(x){
			base::sum(x / base::rank(-x,ties.method="min"))
		})
		
    }
	
	if(verbose){
		now <- Sys.time()
		message(sprintf("%d Genes are defined as seeds and scored using '%s' scoring scheme from %d SNPs", length(seeds.genes), scoring.scheme, ncol(G2S_score)), appendLF=TRUE)
	}
    
    ######################################################################################
    
    if(verbose){
        message(sprintf("\n#######################################################"))
        message(sprintf("'oPierGenes' is being called to prioritise target genes (%s):", as.character(Sys.time())), appendLF=TRUE)
        message(sprintf("#######################################################"))
    }
    
	pNode <- suppressMessages(oPierGenes(data=seeds.genes, network=network, weighted=weighted, network.customised=network.customised, seeds.inclusive=seeds.inclusive, normalise=normalise, restart=restart, normalise.affinity.matrix=normalise.affinity.matrix, GR.Gene=GR.Gene, verbose=verbose, placeholder=placeholder, guid=guid))
	
	if(verbose){
        message(sprintf("#######################################################"))
        message(sprintf("'oPierGenes' has been finished (%s)!", as.character(Sys.time())), appendLF=TRUE)
        message(sprintf("#######################################################\n"))
    }
    
    #######################
    ## if pNode==NULL, return NULL
    if(is.null(pNode)){
    	return(NULL)
    }
    #######################
    
	if(verbose){
		now <- Sys.time()
		message(sprintf("A total of %d genes are prioritised, based on:", nrow(pNode$priority)), appendLF=TRUE)
		message(sprintf("\t%d SNPs scored positively (including %d 'Lead' and %d 'LD');", nrow(df_SNP), sum(df_SNP$Flag=='Lead'), sum(df_SNP$Flag=='LD')), appendLF=TRUE)
		if(!is.null(df_nGenes)){
			message(sprintf("\t%d nearby genes within %d(bp) genomic distance window of %d SNPs", length(unique(df_nGenes$Gene)), distance.max, length(unique(df_nGenes$SNP))), appendLF=TRUE)
		}
		if(!is.null(df_eGenes)){
			message(sprintf("\t%d eQTL genes with expression modulated by %d SNPs", length(unique(df_eGenes$Gene)), length(unique(df_eGenes$SNP))), appendLF=TRUE)
		}
		if(!is.null(df_cGenes)){
			message(sprintf("\t%d HiC genes physically interacted with %d SNP", length(unique(df_cGenes$Gene)), length(unique(df_cGenes$SNP))), appendLF=TRUE)
		}
		message(sprintf("\t%d genes defined as seeds from %d SNPs", length(seeds.genes), ncol(G2S_score)), appendLF=TRUE)
		message(sprintf("\trandomly walk the network (%d nodes and %d edges) starting from %d seed genes/nodes (with %.2f restarting prob.)", vcount(pNode$g), ecount(pNode$g), length(seeds.genes), restart), appendLF=TRUE)
	}
    
    #####
    ## SNP
    df_SNP <- df_SNP[order(df_SNP$Flag,df_SNP$Score,df_SNP$SNP,decreasing=TRUE),]
    
    ## nGenes
    if(is.null(df_nGenes)){
    	nGenes <- NULL
    }else{
		nGenes <- df_nGenes
		ind <- match(nGenes$SNP, df_SNP$SNP)
		nGenes$SNP_Flag <- df_SNP$Flag[ind]
	}
    ## eGenes
    if(is.null(df_eGenes)){
    	eGenes <- NULL
    }else{
		eGenes <- oDefineQTL(data=df_SNP$SNP, include.QTL=include.QTL, QTL.customised=QTL.customised, verbose=FALSE, placeholder=placeholder, guid=guid)
		ind <- match(eGenes$SNP, df_SNP$SNP)
		eGenes$SNP_Flag <- df_SNP$Flag[ind]
	}
    ## cGenes
    if(is.null(df_cGenes)){
    	cGenes <- NULL
    }else{
		cGenes <- oDefineRGB(data=df_SNP$SNP, entity="SNP", include.RGB=include.RGB, GR.SNP=GR.SNP, verbose=FALSE, placeholder=placeholder, guid=guid)
		#cGenes <- cGenes$df
		ind <- match(cGenes$SNP, df_SNP$SNP)
		cGenes$SNP_Flag <- df_SNP$Flag[ind]
    }
    
    #####
    ## append
    pNode[['SNP']] <- df_SNP
    pNode[['Gene2SNP']] <- Gene2SNP
    pNode[['nGenes']] <- nGenes
    pNode[['eGenes']] <- eGenes
    pNode[['cGenes']] <- cGenes
	
    ####################################################################################
    endT <- Sys.time()
    if(verbose){
        message(paste(c("\nFinish at ",as.character(endT)), collapse=""), appendLF=TRUE)
    }
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=TRUE)
    
    invisible(pNode)
}
