#' Function to convert gene symbols to entrez geneid
#'
#' \code{oSymbol2GeneID} is supposed to convert gene symbols to entrez geneid.
#'
#' @param data an input vector containing gene symbols
#' @param org a character specifying an organism. Currently supported organisms are 'human' and 'mouse'. It can be an object 'EG'
#' @param check.symbol.identity logical to indicate whether to match the input data via Synonyms for those unmatchable by official gene symbols. By default, it sets to false
#' @param details logical to indicate whether to result in a data frame (in great details). By default, it sets to false
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to false for no display
#' @param placeholder the characters to tell the location of built-in RData files. See \code{\link{oRDS}} for details
#' @param guid a valid (5-character) Global Unique IDentifier for an OSF project. See \code{\link{oRDS}} for details
#' @return a vector containing entrez geneid with 'NA' for the unmatched if (details set to false); otherwise, a data frame is returned
#' @note If a symbol mapped many times, the one assiged as the "protein-coding" type of gene is preferred.
#' @export
#' @seealso \code{\link{oSymbol2GeneID}}
#' @include oSymbol2GeneID.r
#' @examples
#' \dontrun{
#' 
#' # a) provide the input Genes of interest (eg 100 randomly chosen human genes)
#' ## load human genes
#' org.Hs.eg <- oRDS('org.Hs.eg', placeholder=placeholder)
#' Symbol <- sample(org.Hs.eg$info$Symbol, 100)
#' Symbol
#' 
#' # b) convert into GeneID
#' GeneID <- oSymbol2GeneID(Symbol)
#' 
#' # c) convert into a data frame
#' df <- oSymbol2GeneID(Symbol, details=TRUE)
#' 
#' 
#' # advanced use
#' df <- oSymbol2GeneID(Symbol, org=org.Hs.eg, details=TRUE)
#' }

oSymbol2GeneID <- function(data, org=c("human","mouse"), check.symbol.identity=FALSE, details=FALSE, verbose=TRUE, placeholder=NULL, guid=NULL)
{
    
    if (!is.vector(data)){
        stop("The input data must be a vector.\n")
    }
    Symbol <- as.character(data)
    
    ## load Enterz Gene information
	if(is(org,"EG")){
		df_eg <- org$gene_info
		if(verbose){
			message(sprintf("Customised organism (%s)", as.character(Sys.time())), appendLF=TRUE)
		}
	}else{
		org <- org[1]
		if(org=='human'){
			#df_eg <- xRDataLoader('org.Hs.eg', RData.location=RData.location, guid=guid, verbose=verbose)$gene_info
			df_eg <- oRDS('org.Hs.eg', placeholder=placeholder, guid=guid, verbose=verbose)$info %>% as.data.frame()
		}else if(org=='mouse'){
			#df_eg <- xRDataLoader('org.Mm.eg', RData.location=RData.location, guid=guid, verbose=verbose)$gene_info
			df_eg <- oRDS('org.Mm.eg', placeholder=placeholder, guid=guid, verbose=verbose)$info %>% as.data.frame()
		}
		if(verbose){
			message(sprintf("%s organism (%s)", org, as.character(Sys.time())), appendLF=TRUE)
		}
	}
    
    ## subdived into two parts: "protein-coding" and the rest
    type_of_gene <- ''
    df_eg <- rbind(subset(df_eg,type_of_gene=='protein-coding'), subset(df_eg,type_of_gene!='protein-coding'))
    
	allGeneID <- df_eg$GeneID
	allSymbol <- as.vector(df_eg$Symbol)
	allSynonyms <- as.vector(df_eg$Synonyms)
	allDescription <- as.vector(df_eg$description)
	
    ## correct for those symbols being shown as DATE format
    if(0){
        ## for those starting with 'Mar' in a excel-input date format
        a <- Symbol
        flag <- grep("-Mar$", a, ignore.case=TRUE, perl=TRUE, value=FALSE)
        if(length(flag)>=1){
        	b <- a[flag]
        	c <- sub("-Mar$", "", b, ignore.case=TRUE, perl=TRUE)
        	d <- sub("^0", "", c, ignore.case=TRUE, perl=TRUE)
        	e <- sapply(d, function(x) paste(c("March",x), collapse=""))
        	a[flag] <- e
        	Symbol <- a
    	}

        ## for those starting with 'Sep' in a excel-input date format
        a <- Symbol
        flag <- grep("-Sep$", a, ignore.case=TRUE, perl=TRUE, value=FALSE)
        if(length(flag)>=1){
        	b <- a[flag]
            c <- sub("-Sep$", "", b, ignore.case=TRUE, perl=TRUE)
            d <- sub("^0", "", c, ignore.case=TRUE, perl=TRUE)
            e <- sapply(d, function(x) paste(c("Sept",x), collapse=""))
            a[flag] <- e
            Symbol <- a
        }
    }
    
    ## case-insensitive
    ## only keep the first matched one
    #match_flag <- match(tolower(Symbol),tolower(allSymbol))
    match_flag <- match(Symbol,allSymbol)
    
    ## match via Synonyms for those unmatchable by official gene symbols
    if(check.symbol.identity){
    	## match Synonyms (if not found via Symbol)
        na_flag <- is.na(match_flag)
        a <- Symbol[na_flag]

        ###
        tmp_flag <- is.na(match(tolower(allSymbol), tolower(Symbol)))
        tmp_Synonyms <- allSynonyms[tmp_flag]
        Orig.index <- seq(1,length(allSynonyms))
        Orig.index <- Orig.index[tmp_flag]
        ###

        b <- sapply(1:length(a), function(x){
        	tmp_pattern1 <- paste("^",a[x],"\\|", sep="")
            tmp_pattern2 <- paste("\\|",a[x],"\\|", sep="")
            tmp_pattern3 <- paste("\\|",a[x],"$", sep="")
            tmp_pattern <- paste(tmp_pattern1,"|",tmp_pattern2,"|",tmp_pattern3, sep="")
            tmp_result <- grep(tmp_pattern, tmp_Synonyms, ignore.case=TRUE, perl=TRUE, value=FALSE)
            ifelse(length(tmp_result)==1, Orig.index[tmp_result[1]], NA)
        })
        match_flag[na_flag] <- b
        
        if(verbose){
        	now <- Sys.time()
            message(sprintf("Among %d symbols of input data, there are %d mappable via official gene symbols, %d mappable via gene alias but %d left unmappable", length(Symbol), (length(Symbol)-length(a)), sum(!is.na(b)), sum(is.na(b))), appendLF=TRUE)
        }
    
    }else{
    	if(verbose){
        	now <- Sys.time()
            message(sprintf("Among %d symbols of input data, there are %d mappable via official gene symbols but %d left unmappable", length(Symbol), (sum(!is.na(match_flag))), (sum(is.na(match_flag)))), appendLF=TRUE)
        }
    }
    
	## convert into GeneID
	df_res <- df_eg[match_flag, ]
	
	if(details){
		df_res <- data.frame(Input=Symbol, df_res, stringsAsFactors=FALSE)
		return(df_res)
	}else{
		return(df_res$GeneID)
	}
}
