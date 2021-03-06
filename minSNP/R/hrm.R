#' \code{calGC} is used to calculate the number of GC content in the sequences
#' @param seq is the sequence to be calculated
#' @return will returns the sum of GC content in the sequence
#' @export
calGC <-function(seq){
    sum = 0
    for (a in seq){
        if (a=='C' || a=='c' || a == 'G' || a =='g'){
            sum = sum+1
        }
    }
    return(sum)
} 

#' \code{invalidSeq} is used to check if the sequence contains invalid character
#' @param proSeq is the sequence to be calculated
#' @param acceptedChar is the characters which are considered valid
#' @return will returns the sum of GC content in the sequence
#' @export
invalidSeq<-function(proSeq, acceptedChar = c("A", "a", "C", "c", "T", "t", "G", "g", "-")){
    invalid = FALSE
    for (i in proSeq){
        if( !is.element(i,acceptedChar) ){
            invalid = TRUE
            break
        }
    }
    return(invalid)
}

#' \code{alleles} is used to get the sequences fragments of the alleles to be analysed
#' @param loci_info is the csv file with the following header filename,allele_length,fragment_name,fragment_start,fragment_stop
#' and each rows contain relevant information with ',' seperating the field
#' @param acceptedChar is the characters which are considered valid
#' @return Will returns dataframe of different fragment with matching length, their GC content (curve number) and their fragment sequences
#' @import seqinr
#' @export
alleles <- function(loci_info, acceptedChar = c("A", "a", "C", "c", "T", "t", "G", "g", "-")) {
    info <- read.csv(loci_info, header = TRUE, sep = ",")
    filenames <- c(as.character(info$filename))
    fragment_names <- c(as.character(info$fragment_name))

    ## Initialise a list to hold the dataframes of the alleles and corresponding curve numbers
    df_list <- list()
    print('Processing loci:')
    for(i in 1:length(filenames)) {
	    x <- readFasta(filenames[i])

        ## Check the length of the allele and split up the data in to subsets
	    y <- (getLength(x) == info[i,2])
	    correct_length <- x[y == TRUE]
	    wrong_length <- x[y == FALSE]
        if(length(correct_length)==0){
            print('None of the sequences have the same length')
            next
        }
        
        ## Check for nonstandard character
        pp <- sapply(getFrag(correct_length, info[i,4], info[i,5]), invalidSeq)
        noInvalid <- correct_length[pp == FALSE]
        invalid <- correct_length[pp == TRUE]
        fragment <- getSequence(getFrag(noInvalid, info[i,4], info[i,5]), as.string=TRUE)
        fragment_STR <- sapply(fragment, `[[`, 1)
        gc_content <- sapply(getFrag(noInvalid, info[i,4], info[i,5]), calGC) # returns between 0 and 1


        ## Populate a dataframe with the locus name, allele number and curve number
        gene <- strsplit(as.character(getName(noInvalid)), '_')
        locus_name <- sapply(gene, `[[`,1)
	    allele_number <- sapply(gene, `[[`,2)

        curve <- gc_content 
        fragment_name <- fragment_names[i]
	    df_list[[i]] = data.frame(locus_name, allele_number, fragment_name, curve, fragment_STR)
        
        print(i)
        ignoredL = length(invalid) + length(wrong_length)
        print(paste("------Number of ignored alleles: ", ignoredL, "------"))
        if (length(invalid) > 0){
            print(paste("Ignored (containing invalid characters) : ", length(invalid)))
            print(getName(invalid))
        }
        if (length(wrong_length) > 0){
            print(paste("Ignored (wrong length) :", length(wrong_length)))
            print(getName(wrong_length))
        }
        }
    return(df_list)
    }

#' \code{all2curv} is a helper function that takes in a data frame consisting of ST profile and its fragment ID
#' and combine them with the fragment sequences and fragment curve number
#' @param st_loci is data frame of ST_profile and its corresponding fragment ID
#' @param alleles is the data frame generated using \code{alleles}
#' @return data frame with ST, fragments sequences, fragments curve numbers
all2curv <- function(st_loci, alleles){

    ## This code is per column
    new_ST_CURVE <- st_loci
    frag_name=lapply(alleles, function(x){levels(x[3]$fragment_name[1])})
    frag_name=unlist(frag_name)

    for (gn in colnames(new_ST_CURVE)[-1]){
        current = match(gn, frag_name)[1]
        new_ST_CURVE[,gn] <- alleles[[current]][,4][match(st_loci[,gn], alleles[[current]][,2])]
        new_col=paste(gn,'_seq', sep='')
        new_ST_CURVE[, new_col] <- NA
        new_ST_CURVE[,new_col] <- alleles[[current]][,5][match(st_loci[,gn], alleles[[current]][,2])]
    }
    return(new_ST_CURVE)
}

#' \code{combineLoci} is used to organise the fragment to its respective ST
#' @param st_profile is the location of the ST_profiles.csv that can be obtained in mlst database
#' @param alleles is the data frame generated from alleles
#' @param changeSet is the location of the csv created, see \code{changeCurve} for the format of the csv file
#' @param included is the ST to be included for analysis
#' @param excluded is the ST to be excluded from analysis
#' Note that included is used before excluded, a range can be used, e.g. to included ST 1 to 100,
#' included=c(1:100), can be similarly applied to excluded
#' @return data frame with ST, fragments sequences, fragments curve numbers
#' @export
combineLoci <- function(st_profile, alleles, changeSet=NULL, included=NULL, excluded=NULL){
    ## Load a tab-delimited .csv file of the ST profiles and save as a dataframe
    ST_profiles <- read.table(st_profile, header = TRUE, sep = "\t")
    ## Create a duplicate of 'ST_profile' for content editing
    ST_curves <- ST_profiles
    # Drop any columns that is not in the profile.
    name=lapply(alleles, function(x){levels(x[1]$locus_name[1])})
    frag_name=unlist(lapply(alleles, function(x){levels(x[3]$fragment_name)}))
    allele_name = c("ST", unlist(name))
    ST_curves <- subset(ST_curves, select=allele_name)
    col_frag_name=c('ST')

    #Start log file for analysis
    line = paste("Analysis started at: ", Sys.time())
    write(line,file="log.txt",append=TRUE)


    for(i in 2:length(colnames(ST_curves))){ 
        #The order of allele ST_profile and loci_info do not have to be the same
        for(j in frag_name){
            #Support for multiple original similar colnames separated by .
            cname = strsplit(colnames(ST_curves)[i], '\\.')[[1]][1] 
            if(grepl(cname, j)){
                col_frag_name = c(col_frag_name, j)
                frag_name=setdiff(frag_name, j)
                break
            }
        }
    }
    colnames(ST_curves) <- col_frag_name
    new_ST_curves <- all2curv(ST_curves, alleles)
    new_ST_curves <- na.omit(new_ST_curves)
    rownames(new_ST_curves) <- NULL
    if (!is.null(included)){
        new_ST_curves <- new_ST_curves[new_ST_curves$'ST' %in% included,]
    }
    if (!is.null(excluded)){
        new_ST_curves <- new_ST_curves[!new_ST_curves$'ST' %in% excluded,]
    }
    droppedST <- setdiff(ST_profiles$ST, new_ST_curves$ST)
    line = paste("ST dropped from analysis are: ", paste(droppedST, collapse=","))
    write(line,file="log.txt",append=TRUE)
    new_ST_curves <- changeCurve(changeSet, new_ST_curves)
    return(new_ST_curves)
}

#' \code{changeCurve} is used to modify the curve number if the conditions are met. To specify conditions, create a csv file with the following headings:
#' fragment_name, type, replace_with, if_match_fragment, if_match_curve,
#' where type is
#' 1. replace if match fragment
#' 2. replace if match curve number
#' for non-relevant field, use _
#' @param changeSet is the location of the csv created
#' @param c_st_all is the dataframe generated using \code{combineLoci}
#' @return the dataframe with curve number modified. All rows with NAs are also dropped. All columns with fragment sequence also dropped.
#' @export
changeCurve <- function(changeSet=NULL, c_st_all){
    if (!is.null(changeSet)){
        changeSet<-read.csv(changeSet, header = TRUE, colClasses = "character", sep = ",")
        for (change in 1:length(changeSet)){
            replaceAt <- changeSet[change, 'fragment_name']
            replaceVal <- changeSet[change, 'replace_with']
            # type 1 change - match fragment
            if(changeSet[change,]$type==1){
                mat_frag <- paste(changeSet[change, 'fragment_name'], '_seq', sep='')
                condition <- tolower(changeSet[change, 'if_match_fragment'])
                for (i in 1: nrow(c_st_all)){
                    if ((!is.na(c_st_all[i, mat_frag]))&&(c_st_all[i, mat_frag]==condition)){
                        line = paste("Changed curve value of ", c_st_all[i, 'ST'], "from ", 
                                c_st_all[i, replaceAt], " to ", replaceVal, " due to rule : ", change)
                        c_st_all[i, replaceAt] <- replaceVal
                        write(line,file="log.txt",append=TRUE)
                    }
                }
            }#type 2 change - match curve number
            else if(changeSet[change,]$type==2){
                condition <- changeSet[change, 'if_match_curve']
                for (i in 1: nrow(c_st_all)){
                    if ((!is.na(c_st_all[i, replaceAt]))&&(c_st_all[i, replaceAt]==condition)){
                        line = paste("Changed curve value of ", c_st_all[i, 'ST'], "from ", 
                                c_st_all[i, replaceAt], " to ", replaceVal, " due to rule : ", change)
                        c_st_all[i, replaceAt] <- replaceVal
                        write(line,file="log.txt",append=TRUE)
                    }
                }
            }
        }
    }
    #Dropping unnecessary columns
    m=colnames(c_st_all)[!substr(colnames(c_st_all), nchar(colnames(c_st_all))-4+1, nchar(colnames(c_st_all)))=="_seq"]
    c_st_all <- c_st_all[,m]
    #Also dropping any rows with NA
    c_st_all <- c_st_all[complete.cases(c_st_all),]
    return(c_st_all)
}

#' \code{generateID} is used to generate the MElT ID
#' @param c_st_all is the data frame from changeCurve
#' @return the same data frame with MElT ID appended to the last column
#' @export
generateID <- function(c_st_all){
    c_st_all[,'colID'] <- unlist(apply(c_st_all[,-1], 1, paste, collapse=''))
    c_st_all[, 'MElT-ID'] <- NA
    mapping = list()
    id = 1
    for (i in 1:nrow(c_st_all)){
        a = toString(c_st_all[i, 'colID'])
        if( is.null(mapping[[a]]) ){
            mapping[[a]] = id
            c_st_all[i,'MElT-ID'] = id
            id = id+1
        }else{
            c_st_all[i,'MElT-ID'] = mapping[[a]]
        }
    }
    #Drop colID
    c_st_all=c_st_all[,!names(c_st_all) %in% c('colID')]
    return(c_st_all)
}

#' \code{hrm.pattern} is used to generate the pattern for calculating simpson index
#' @param result the result from \code{generateID}
#' @return a pattern that can be input into \code{simpson.calculate}
hrm.pattern <- function(result){
type=list()
    ID = result$'MElT-ID'
	for(al in 1:length(ID)){
        name=as.character(ID[al])
        if (is.null(type[[name]])){
            type[[name]]=1}
        else{
            type[[name]]=type[[name]]+1
        }
	}
	return(type)
}

#' \code{hrm.dIndex} returns the simpson index of the result based on the MElT-ID
#' @param result is the result from \code{generateID}
#' @return will return simpson index
#' @export
hrm.dIndex <- function(result){
    pat = hrm.pattern(result)
    total = length(result$'MElT-ID')
    return(simpson.calculate(pat, total))
}

#' \code{setID} reads from csv file and out
#' @param result is the result from \code{generateID}
#' @param newLocation is the location of the csv file
#' @return will return the new dataframe with MElT-ID appended
#' @export
setID <- function (newLocation, result = NULL ) {
    n <- read.csv(newLocation, header = TRUE, sep = ",", check.names = FALSE)
    a <- colnames(result)[colnames(result) != c("ST", "MElT-ID")]
    b <- colnames(n)[colnames(n) != c("IsolateName")]
    if (length(setdiff(union(a, b), intersect(a, b))) > 0) {
        print("Problem with new data OR not result data: Check column names")
        return(n)
    }
    result$ST <- NULL
    result <- unique(result)
    newResult <- merge(x = result, y = n, by = b, all.y = TRUE)
    newResult<-newResult[c(ncol(newResult), 1:ncol(newResult)-1)]
    return(newResult)
}


#
#
#The ideal situation for us would for all the information about the alleles and STs that have been removed,
# as well as information on the manual changes that are made in the changeCurve function, output to a file.

# Would it be possible to also output the STs which have been ignored in the generation of the MelT-IDs please?

#That way we have a record we can save that documents exactly what happened for the run.