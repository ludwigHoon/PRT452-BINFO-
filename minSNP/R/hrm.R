#' \code{hrmAnalysis} is used to return the table of GC contents at different regions
#' @import seqinr
#' @param sequences a list of SeqFastadna.
#' @param regions the list of regions and their positions (start and end defined).
#' regions can be prepared in this manner:
#' region=list()
#' region[[1]]=list('s'=73,'e'=210)
#' region[[2]]=list('s'=521,'e'=617)
#' region[[3]]=list('s'=1637,'e'=1680)
#' region[[4]]=list('s'=2100,'e'=2216)
#' region[[5]]=list('s'=2316,'e'=2485)
#' region[[6]]=list('s'=2521,'e'=2644)
#' names(region)<-c('arcC78/210', 'aroE88/155', 'gmk286', 'pta294', 'tpi36', 'tpi241/243')
#' @param excluded the list of excluded STs.
#' @return Will return a matrix of GC content of STs at different regions.
#' @export
hrmAnalysis <-function(sequences, regions, excluded=NULL){
    tSeq = sequences
    for(exc in excluded){
        tSeq[[exc]]<-NULL
    }
    result = matrix(nrow=length(tSeq), ncol=length(regions))
    rownames(result) = names(tSeq)
    if(!is.null(names(regions))){
        colnames(result) = names(regions)
    }
    for(nm in names(tSeq)){
        for(re in 1:length(regions)){
            seqReg = tSeq[[nm]][regions[[re]]$s:regions[[re]]$e]
            result[nm, re] = GC(seqReg)*length(seqReg)
        }
    }
    return(result)
}

changeSets=list()
changeSets[[1]] = list(m="ATTTAATTAAAGAAATTATTTCAAAAAAAGAATTAGATGGCTTTAATATCACAATTCCTCATAAAGAGCGTATCATACCGTATTTAGATCATGTTGA", r=24.5, reg=2)
changeSets[[2]] = list(m="ATTTAATTAAAGAAATTATTTCAAAAAAAGAATTAGATGGCTTTAATATCACAATTCCCCATAAAGAACGTATCATACCGTATTTAGATTATGTTGA", r=23.5, reg=2)
changeSets[[3]] = list(m="TGCTATTTTCAAACATGGTATGACACCAATTATTTGTGTTGGTGAAACAGACGAAGAGCGTGAAAGCGGTAAAGCTAACGATGTTGTAGGTGAGCAAGTTAAGAAAGCTGTTGCAGGTTTATCTGAAGAGCAACTTAAATCAGTTGTAATTGCTTATGAACCAATCTGGG", r=67, reg=5)
#changeSets[[4]] = list(m="TGCTATTTTCAAACATGGTATGACACCAATTATTTGTGTTGGTGAAACAGACGAAGAGCGTGAAAGCGGTAAAGCTAACGATGTTGTAGGTGAGCAAGTTAAGAAAGCTGTTGCAGGTTTATCTGAAGAGCAACTTAAATCAGTTGTAATTGCTTATGAACCAATCTGGG", r=65, reg=5)
changeSets[[4]] = list(m="GCGAATGAAATGTGTGCATTTGTACGTCAAACTATTTCTGACTTATCAAGCAAAGAAGTATCAGAAGCAACTCGTATTCAATATGGTGGTAGTGTTAAACCTAACAACATTAAAGAATACATGG", r=42, reg=6)
hrmAddChanges <-function(analysis, changeSets, sequences, regions){
    for(changes in changeSets){
        for(re in 1:ncol(analysis)){
            if(re!= changes$reg){
                next
            }
            for(nm in rownames(analysis)){
                seqReg = toupper(paste(sequences[[nm]][regions[[re]]$s:regions[[re]]$e], collapse=''))
                print(seqReg==changes$m)
                if(seqReg==changes$m){ 
                    analysis[nm, re] = changes$r
                }
            }
        }
    }
    return(analysis)
}