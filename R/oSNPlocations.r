#' Function to extract genomic locations given a list of SNPs
#'
#' \code{oSNPlocations} is supposed to extract genomic locations given a list of SNPs.
#'
#' @param data a input vector containing SNPs. SNPs should be provided as dbSNP ID (ie starting with rs). Alternatively, they can be in the format of 'chrN:xxx', where N is either 1-22 or X, xxx is genomic positional number; for example, 'chr16:28525386'
#' @param GR.SNP the genomic regions of SNPs. By default, it is 'dbSNP_GWAS', that is, SNPs from dbSNP (version 146) restricted to GWAS SNPs and their LD SNPs (hg19). It can be 'dbSNP_Common', that is, Common SNPs from dbSNP (version 146) plus GWAS SNPs and their LD SNPs (hg19). Alternatively, the user can specify the customised input. To do so, first save your RData file (containing an GR object) into your local computer, and make sure the GR object content names refer to dbSNP IDs. Then, tell "GR.SNP" with your RData file name (with or without extension), plus specify your file RData path in "RData.location". Note: you can also load your customised GR object directly
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param placeholder the characters to tell the location of built-in RData files. See \code{\link{oRDS}} for details
#' @param guid a valid (5-character) Global Unique IDentifier for an OSF project. See \code{\link{oRDS}} for details
#' @return 
#' an GR oject, with an additional metadata column called 'variant_id' storing SNP location in the format of 'chrN:xxx', where N is either 1-22 or X, xxx is genomic positional number.
#' @note none
#' @export
#' @seealso \code{\link{oSNPlocations}}
#' @include oSNPlocations.r
#' @examples
#' \dontrun{
#' # a) provide the SNPs
#' ## get lead SNPs reported in AS GWAS and their significance info (p-values)
#' gr <- ImmunoBase$AS$variant
#' data <- names(gr)
#'
#' # b) find the location
#' snp_gr <- oSNPlocations(data=data, RData.location=RData.location)
#' }

oSNPlocations <- function(data, GR.SNP=c("dbSNP_GWAS","dbSNP_Common","dbSNP_Single"), verbose=TRUE, placeholder=NULL, guid=NULL)
{
	
	## replace '_' with ':'
	data <- gsub("_", ":", data, perl=TRUE)
	## replace 'imm:' with 'chr'
	data <- gsub("imm:", "chr", data, perl=TRUE)
	
	data <- unique(data)
	
    ######################################################
    # Link to targets based on genomic distance
    ######################################################
    
  	## load positional information
	if(verbose){
		now <- Sys.time()
		message(sprintf("Load positional information for SNPs (%s) ...", as.character(now)), appendLF=TRUE)
	}
	if(is(GR.SNP,"GRanges")){
		pos_SNP <- GR.SNP
	}else{
		pos_SNP <- oRDS(GR.SNP[1], verbose=verbose, placeholder=placeholder, guid=guid)
		if(is.null(pos_SNP)){
			GR.SNP <- "dbSNP_GWAS"
			if(verbose){
				message(sprintf("Instead, %s will be used", GR.SNP), appendLF=TRUE)
			}
			pos_SNP <- oRDS(GR.SNP, verbose=verbose, placeholder=placeholder, guid=guid)
		}
	}
	
  	ind <- match(data, names(pos_SNP))
  	data_rest <- data[is.na(ind)]
  	ind <- ind[!is.na(ind)]
  	if(length(ind)){
  		gr_SNP <- pos_SNP[ind,]
  		GenomicRanges::mcols(gr_SNP) <- NULL
  	}else{
  		gr_SNP <- NULL
  	}
  	
  	#######################################################
  	### deal with data_rest having the format as: chr\w+:\d+
  	ind <- grep("^chr\\w+:\\d+", data_rest, perl=TRUE)
  	if(length(ind)>0){
		data_rest <- data_rest[ind]
		res_ls <- strsplit(data_rest, ":")
		res_df <- do.call(rbind, res_ls)
		res_gr <- GenomicRanges::GRanges(
			seqnames=res_df[,1],
			ranges = IRanges::IRanges(start=as.numeric(res_df[,2]), end=as.numeric(res_df[,2]), names=data_rest),
			strand = rep('*',nrow(res_df))
		)
	}else{
		res_gr <- NULL
	}
	### combine gr_SNP and res_gr
	if(!is.null(gr_SNP)){
		if(!is.null(res_gr)){
			gr_SNP <- c(gr_SNP, res_gr)
		}
	}else{
		if(!is.null(res_gr)){
			gr_SNP <- res_gr
		} 
	}
  	if(verbose){
		now <- Sys.time()
		message(sprintf("\tOut of %d input SNPs, %d SNPs have positional info", length(data), length(gr_SNP)), appendLF=TRUE)
  	}
	
	if(!is.null(gr_SNP)){
		tmp_df <- GenomicRanges::as.data.frame(gr_SNP, row.names=NULL)
		mcols_df <- data.frame(variant_id=paste(tmp_df[,1],':',tmp_df[,3],sep=''), stringsAsFactors=FALSE)
		GenomicRanges::mcols(gr_SNP) <- mcols_df
	}
	
    invisible(gr_SNP)
}
