#' getPCRsizes
#'
#' returns length of product given a GRanges span and GRangesList of transcripts
#'
#' @param pcr_span GRanges object
#' @param txlist GRangesList object
#' @param verbose logical. report intermediate output.
#'
#' @return data.frame of transcript names with detected sizes in basepairs
#'
#' @import GenomicRanges
#' @importFrom S4Vectors queryHits subjectHits elementNROWS
#' @importFrom plyr ldply
#' @export
#'
#' @author Diana Low
#'
#' @examples
#' suppressMessages(library(BSgenome.Mmusculus.UCSC.mm9))
#' bsgenome<-BSgenome.Mmusculus.UCSC.mm9
#' ## create a primer pair
#' ## for actual use, obtain primer pair from primer design (callPrimer3)
#' primer_pair <- data.frame(PRIMER_LEFT_SEQUENCE="agctcttgaaattggagctgac",
#'                           PRIMER_RIGHT_SEQUENCE="cttagaaagaacaggaaatcc",
#'                           stringsAsFactors=FALSE)
#'
#' ## confirm location
#' cp<-checkPrimer(primer_pair,bsgenome,roi)
#' cp
#'
#' ## get the PCR sizes
#' pcr_result1 <- getPCRsizes(cp,theexons)
getPCRsizes<-function(pcr_span,txlist,verbose=FALSE){
  pcr_overlap_total<-findOverlaps(pcr_span$total_span,txlist)
  pcr_hits<-txlist[subjectHits(pcr_overlap_total)]

  pcr_overlaps_left<-lapply(pcr_hits,function(x) findOverlaps(pcr_span$primer_left_span,x))
  pcr_overlaps_right<-lapply(pcr_hits,function(x) findOverlaps(pcr_span$primer_right_span,x))
  pol<-names(pcr_overlaps_left[elementNROWS(pcr_overlaps_left)!=0])
  por<-names(pcr_overlaps_right[elementNROWS(pcr_overlaps_right)!=0])

  pcr_hits<-pcr_hits[names(pcr_hits) %in% pol[pol %in% por]]

  txsizes<-ldply(pcr_hits,function(x) {
    po2<-findOverlaps(pcr_span$total_span,x)
    sum(width(overlapsRanges(ranges(pcr_span$total_span),ranges(x),po2)))
    })

  colnames(txsizes)<-c("ID","bp")
  if(verbose) message(paste(txsizes$ID,":",txsizes[,2],"bp\n"))

  return(txsizes)
}

#' splitPCRhit
#'
#' splits the PCR alignment into the two AS conditions
#'
#' @param res result from \code{\link{getPCRsizes}}
#' @param hitlist \code{\link{findCompatibleEvents}} object
#'
#' @return list of 2 data.frame objects with isoform name (ID) and
#'   length of PCR product (bp) matching Type 1 or Type 2 transcripts
#'
#' @export
#'
#' @author Diana Low
#'
#' @examples
#' ## as getPCRsizes gives you all PCR bands when the primers are used,
#' ## splitPCRhit will determine which bands are relevant to the target
#' relevant_pcr_bands<-splitPCRhit(pcr_result1,compatible_tx)
splitPCRhit<-function(res,hitlist){
  #q1<-res[res$ID %in% names(hitlist$hits[[1]]),]
  #q2<-res[res$ID %in% names(hitlist$hits[[2]]),]
  res[which(res$ID %in% names(hitlist$hits[[1]])),"type"]<-1
  res[which(res$ID %in% names(hitlist$hits[[2]])),"type"]<-2
  return(res)
}

#' getRegionDNA
#'
#' get DNA sequence give a region of interest
#'
#' @param roi \code{\link{makeROI}} object
#' @param genome BSgenome object
#' @param introns TRUE/FALSE. whether to include intronic (lowercase) DNA.
#'   By default returns only exonic (uppercase) DNA.
#'
#' @return list of \cr
#'   (1) DNA sequence
#'   (2) Junction start (for primer design)
#'
#' @export
#'
#' @author Diana Low
#'
#' @examples
#' suppressMessages(library(BSgenome.Mmusculus.UCSC.mm9))
#' bsgenome<-BSgenome.Mmusculus.UCSC.mm9
#' getRegionDNA(roi,bsgenome)
getRegionDNA<-function(roi,genome,introns=FALSE){
  #get introns
  if(roi$type %in% c("A5SS","A3SS","RI")){
    seq_pos=2
  } else {
    seq_pos=1
  }

  if(introns){
    parts<-getSeq(genome,range(roi$roi_range[[seq_pos]]))
    final_seq<-toString(parts)
  } else {
    parts<-getSeq(genome,roi$roi_range[[seq_pos]])
    final_seq<-paste(toString(parts[1]),"[]",toString(parts[length(parts)]),sep="")
  }

  jstart=unlist(gregexpr(pattern ="\\[",final_seq))
  return(list(seq=final_seq,jstart=jstart))
}

#' checkPrimer
#'
#' @param pp data.frame defining primers, or output of \code{\link{callPrimer3}}.
#'   minimal columns = PRIMER_LEFT_SEQUENCE,PRIMER_RIGHT_SEQUENCE
#' @param genome BSgenome object
#' @param roi \code{\link{makeROI}} object
#'
#' @return list of GRanges with primer locations
#'
#' @import GenomeInfoDb BSgenome.Mmusculus.UCSC.mm9
#' @export
#'
#' @author Diana Low
#'
#' @examples
#' # create a primer pair
#' roi
#' primer_pair <- data.frame(PRIMER_LEFT_SEQUENCE="agctcttgaaattggagctgac",
#'                           PRIMER_RIGHT_SEQUENCE="cttagaaagaacaggaaatcc",
#'                           stringsAsFactors=FALSE)
checkPrimer<-function(pp,genome,roi=NULL){
  roi=roi$roi
  pchr<-as.character(seqnames(roi))[1]
  strand<-as.character(strand(roi))[1]

  if(strand=="+"){
    lseq<-pp$PRIMER_LEFT_SEQUENCE
    rseq<-as.character(reverseComplement(DNAStringSet(pp$PRIMER_RIGHT_SEQUENCE)))
    aaL<-start(matchPattern(lseq,genome[[pchr]]))
    aaR<-end(matchPattern(rseq,genome[[pchr]]))
  } else{
    lseq<-as.character(reverseComplement(DNAStringSet(pp$PRIMER_LEFT_SEQUENCE)))
    rseq<-pp$PRIMER_RIGHT_SEQUENCE
    aaR<-end(matchPattern(lseq,genome[[pchr]])[1])
    aaL<-start(matchPattern(rseq,genome[[pchr]])[1])
  }
  leftspan<-matchPattern(lseq,genome[[pchr]])
  rightspan<-matchPattern(rseq,genome[[pchr]])
  return(list(total_span=GRanges(seqnames=pchr,ranges=IRanges(aaL,aaR)),
              primer_left_span=GRanges(seqnames=pchr,ranges=IRanges(start(leftspan),end(leftspan))),
              primer_right_span=GRanges(seqnames=pchr,ranges=IRanges(start(rightspan),end(rightspan)))))
}

#' callPrimer3
#'
#' call primer3 for a given set of DNAstringSet object
#'
#' modified to include SEQUENCE_TARGET as an option
#'
#' @param seq DNAstring object, one DNA string for the given amplicon
#' @param size_range default: '151-500'
#' @param Tm melting temprature parameters default:c(55,57,58)
#' @param name name of the amplicon in chr_start_end format
#' @param primer3 primer3 path
#' @param thermo.param thermodynamic parameters folder
#' @param sequence_target If one or more targets is specified then a legal
#'        primer pair must flank at least one of them.
#' @param settings text file for parameters
#'
#' @return data.frame of designed primers and parameters
#'
#' @importFrom utils read.delim read.table
#' @export
#'
#' @author Altuna Akalin's modified Arnaud Krebs' original function further modified here by Diana Low
#' @examples
#' ### NOT RUN ###
#' # primer_results<-callPrimer3(seq='')
#'
callPrimer3<-function(seq,size_range='150-500',Tm=c(57,59,62), name = "Primer1",
                       primer3="primer3-2.3.7/bin/primer3_core",
                       thermo.param="primer3-2.3.7/src/primer3_config/",
                       sequence_target=NULL,
                       settings="primer3-2.3.7/primer3web_v4_0_0_default_settings.txt") {

  #check if primer3 path/files are present
  if(sum(file.exists(primer3,thermo.param,settings))!=3){
    message("Please check your Primer3 paths!")
    return(NULL)
  }
  # make primer 3 input file
  p3.input=tempfile()
  p3.output=tempfile()
  cmmd<-paste(sprintf("SEQUENCE_ID=%s\n",name),
              sprintf("SEQUENCE_TEMPLATE=%s\n",as.character(seq)),
              "PRIMER_TASK=pick_detection_primers\n",
              "PRIMER_PICK_LEFT_PRIMER=1\n" ,
              "PRIMER_PICK_INTERNAL_OLIGO=0\n",
              "PRIMER_PICK_RIGHT_PRIMER=1\n",
              "PRIMER_EXPLAIN_FLAG=1\n",
              "PRIMER_PAIR_MAX_DIFF_TM=5\n",
              sprintf("PRIMER_MIN_TM=%s\n" ,Tm[1]),
              sprintf("PRIMER_OPT_TM=%s\n" ,Tm[2]),
              sprintf("PRIMER_MAX_TM=%s\n" ,Tm[3]),
              sprintf("PRIMER_PRODUCT_SIZE_RANGE=%s\n" ,size_range),
              sprintf("PRIMER_THERMODYNAMIC_PARAMETERS_PATH=%s\n" ,thermo.param),
              "=",
              sep='')
  if(!is.null(sequence_target)){
    cmmd<-paste(sprintf("SEQUENCE_TARGET=%s,2\n",sequence_target),cmmd)
  }
  # print(cmmd)
  write(cmmd,p3.input)
  #call primer 3 and store the output in a temporary file
  try(system(
    paste(primer3 ,p3.input, "-p3_settings_file",settings,
          ">", p3.output)
  ))

  #import and parse the output into a dataframe named designed.primers
  out <- read.delim(p3.output, sep='=', header=FALSE)

  unlink(c(p3.input,p3.output) ) # delete temp files
  returned.primers=as.numeric(as.vector(out[out[,1]=='PRIMER_PAIR_NUM_RETURNED',][,2]))
  if (length(returned.primers)==0){ warning('primers not detected for ',name,call. = FALSE);return(NA)}
  if ((returned.primers)==0){ warning('primers not detected for ',name,call. = FALSE);return(NA)}
  if (returned.primers>0){
    designed.primers=data.frame()
    for (i in seq(0,returned.primers-1,1)){
      #IMPORT SEQUENCES
      id=sprintf(  'PRIMER_LEFT_%i_SEQUENCE',i)
      PRIMER_LEFT_SEQUENCE=as.character(out[out[,1]==id,][,2])
      id=sprintf(  'PRIMER_RIGHT_%i_SEQUENCE',i)
      PRIMER_RIGHT_SEQUENCE=as.character(out[out[,1]==id,][,2])

      #IMPORT PRIMING POSITIONS
      id=sprintf(  'PRIMER_LEFT_%i',i)
      PRIMER_LEFT=as.numeric(unlist(strsplit(as.vector((out[out[,1]==id,][,2])),',')))
      #PRIMER_LEFT_LEN=as.numeric(unlist(strsplit(as.vector((out[out[,1]==id,][,2])),',')))
      id=sprintf(  'PRIMER_RIGHT_%i',i)
      PRIMER_RIGHT=as.numeric(unlist(strsplit(as.vector((out[out[,1]==id,][,2])),',')))
      #IMPORT Tm
      id=sprintf(  'PRIMER_LEFT_%i_TM',i)
      PRIMER_LEFT_TM=as.numeric(as.vector((out[out[,1]==id,][,2])),',')
      id=sprintf(  'PRIMER_RIGHT_%i_TM',i)
      PRIMER_RIGHT_TM=as.numeric(as.vector((out[out[,1]==id,][,2])),',')

      res=out[grep(i,out[,1]),]
      extra.inf=t(res)[2,,drop=FALSE]
      colnames(extra.inf)=sub( paste("_",i,sep=""),"",res[,1])
      extra.inf=extra.inf[,-c(4:9),drop=FALSE] # remove redundant columns
      extra.inf=apply(extra.inf,2,as.numeric)
      #Aggegate in a dataframe
      primer.info=data.frame(i,
                             PRIMER_LEFT_SEQUENCE,PRIMER_RIGHT_SEQUENCE,
                             PRIMER_LEFT_TM, PRIMER_RIGHT_TM,
                             PRIMER_LEFT_pos=PRIMER_LEFT[1],
                             PRIMER_LEFT_len=PRIMER_LEFT[2],
                             PRIMER_RIGHT_pos=PRIMER_RIGHT[1],
                             PRIMER_RIGHT_len=PRIMER_RIGHT[2],
                             t(data.frame(extra.inf)),stringsAsFactors = FALSE

      )
      rownames(primer.info)=NULL
      designed.primers=rbind(designed.primers, primer.info)
    }
  }
  return(designed.primers)
}
