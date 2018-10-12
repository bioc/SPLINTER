#' annotateEvents
#'
#' Gives detailed description of splicing event in terms of splicing outcome
#' post translation. Currently supports exon skipping and intron retention events.
#'
#' @param thedata list. output of extractSpliceEvents.
#' @param db TxDb object
#' @param bsgenome BSGenome object
#' @param outputdir character. relative output directory to current location.
#' @param full_output logical. writes out detailed text report and generate figures.
#' @param output_prefix character. text prefix for full_output files.
#'
#' @return list containing information on
#' (1) data.frame with splicing regions
#' (2) splice event type
#'
#'
#' @author Diana LOW
#'
#' @import GenomicFeatures googleVis ggplot2 grDevices
#' @export
annotateEvents <- function(thedata,db,bsgenome,outputdir,full_output=FALSE,output_prefix="results"){
  df<-thedata$data
  thecds<-cdsBy(db,by="tx",use.names=TRUE)
  theexons<-exonsBy(db,by="tx",use.names=TRUE)

  if(full_output){
    closeAllConnections() # to prevent writing errors
    r_file_name.ls<-c()
    p_file_name.ls<-c()
    ifelse(!dir.exists(file.path(getwd(), outputdir)), dir.create(file.path(getwd(), outputdir)), FALSE)
    ifelse(!dir.exists(file.path(getwd(), outputdir,"text")), dir.create(file.path(getwd(), outputdir,"text")), FALSE)
    ifelse(!dir.exists(file.path(getwd(), outputdir,"img")), dir.create(file.path(getwd(), outputdir,"img")), FALSE)
  }

  loops=nrow(df)
  loop_frame=list(data=cbind.data.frame(seq_id=seq.int(nrow(df[1:loops,])),
                                        df[1:loops,],
                                        event=thedata$type,
                                        VALID_TX='',VALID_CDS='',
                                        COMPATIBLE_TX='',COMPATIBLE_CDS='',
                                        ASoutcome1='',ASoutcome2='',
                                        #AALENGTH1="",AALENGTH2="",
                                        stringsAsFactors=F),type=thedata$type)

  for(i in 1:nrow(loop_frame$data)){
    if(full_output){
      r_file_name<-paste("text","/",loop_frame$data[i,]$seq_id,"_",loop_frame$data$ID[i],"_",output_prefix,"_full_report.txt",sep="")
      p_file_name<-paste("img","/",loop_frame$data[i,]$seq_id,"_",loop_frame$data$ID[i],"_",output_prefix,"_full_report.png",sep="")
      r_file_name.ls<-c(r_file_name.ls,r_file_name)
      zz <- file(paste(file.path(getwd(),outputdir),"/",r_file_name,sep=""), open="wt")
      sink(zz)
      sink(zz, type="message")
    }

    message('==================================================')
    message(paste("Checking",loop_frame$data$ID[i],loop_frame$data$geneSymbol[i]))
    message(paste(loop_frame$data[i,]$seq_id,loop_frame$data[i,]$uID,sep="_"))

    message('--> Finding valid transcripts')
    valid_tx <- findTX(id=loop_frame$data[i,]$ID,tx=theexons,db=db,valid=FALSE)

    message('--> Finding valid coding transcripts')
    valid_cds <- findTX(id=loop_frame$data[i,]$ID,tx=thecds,db=db,valid=FALSE)

    if(length(valid_tx)==0){
      message("** no valid transcripts found, skipping to next event.")
      p_file_name.ls<-c(p_file_name.ls,"")
      next
    }

    if(length(valid_cds)==0){
      message("** no valid coding transcripts found, skipping to next event.")
      p_file_name.ls<-c(p_file_name.ls,"")
      next
    }

    roi <- makeROI(loop_frame$data[i,],type=thedata$type)

    message('--> Finding compatible transcripts')
    compatible_tx_events<-findCompatibleEvents(valid_tx,roi=roi,verbose=FALSE) #use tx

    loop_frame$data[i,]$VALID_TX<-as.numeric(length(valid_tx))

    ### Check for compatible transcripts ###
    if(length(compatible_tx_events$hits)==0){
      message('.... no compatible transcripts found. skipping to next event.')
      p_file_name.ls<-c(p_file_name.ls,"")
      #loop_frame$data[i,]$ASoutcome<-"no compatible transcripts"
      next
    }
    p_file_name.ls<-c(p_file_name.ls,p_file_name)
    png(paste(file.path(getwd(),outputdir),"/",p_file_name,sep=""),width = 1000,height = 800)
    eventPlot(transcripts=valid_tx,roi_plot=roi,annoLabel=loop_frame$data[i,]$ID,rspan=2000)
    dev.off()

    message(paste("....",length(compatible_tx_events$hits),"compatible transcripts found. check if coding..."))
    loop_frame$data[i,]$COMPATIBLE_TX<-paste(length(compatible_tx_events$hits[[1]]),
                                             length(compatible_tx_events$hits[[2]]),sep=",")

    message('--> Finding compatible coding transcripts')
    compatible_cds_events<-findCompatibleEvents(valid_cds,roi=roi,verbose=FALSE) #use tx
    loop_frame$data[i,]$VALID_CDS<-as.numeric(length(valid_cds))

    ### Check for compatible cds transcripts ###
    if(length(compatible_cds_events$hits)==0){
      message("--> no compatible coding transcripts")
      loop_frame$data[i,]$ASoutcome<-"no compatible CDS"
      next
    }

    message("--> Looking for CDS transcripts with full cassette")
    if(length(compatible_cds_events$hits[[1]])!=0){
      message(paste("....",length(compatible_cds_events$hits[[1]]),"found."))
      #continue with exon skipping
      if(loop_frame$data[i,]$event=='SE') {
        newregion<-removeRegion(compatible_cds_events$hits[[1]],roi)
      } else if(loop_frame$data[i,]$event=='RI') {
        newregion<-insertRegion(compatible_cds_events$hits[[1]],roi)
      }
      ans<-eventOutcomeCompare(compatible_cds_events$hits[[1]],newregion,bsgenome,direction=loop_frame$data[i,]$IncLevelDifference>0)
      loop_frame$data[i,]$ASoutcome1<-ans$eventtypes
    } else {
      message(".... None found.")
    }

    message("**** Looking for CDS transcripts with skipped cassette")
    if(length(compatible_cds_events$hits[[2]])!=0){
      #message("**** no compatible coding transcripts of type1")
      #loop_frame$data[i,]$ASoutcome<-"no type1 CDS, possible insertion"

      #have transcripts of type 2
      #do exon insertion
      message(paste("....",length(compatible_cds_events$hits[[2]]),"found."))
      if(loop_frame$data[i,]$event=='SE') {
        newregion<-insertRegion(compatible_cds_events$hits[[2]],roi)
      } else if(loop_frame$data[i,]$event=='RI') {
        newregion<-removeRegion(compatible_cds_events$hits[[2]],roi)
      }
      ans<-eventOutcomeCompare(compatible_cds_events$hits[[2]],newregion,bsgenome,direction=loop_frame$data[i,]$IncLevelDifference>0)
      loop_frame$data[i,]$ASoutcome2<-ans$eventtypes
    } else {
      message(".... None found.")
    }



    loop_frame$data[i,]$COMPATIBLE_CDS<-paste(length(compatible_cds_events$hits[[1]]),
                                              length(compatible_cds_events$hits[[2]]),sep=",")

    if(full_output){
      sink()
      closeAllConnections()
    }
    # if(loop_frame$data[i,]$event=='SE') {
    #   newregion<-removeRegion(compatible_cds_events$hits[[1]],roi)
    # } else if(loop_frame$data[i,]$event=='RI') {
    #   newregion<-insertRegion(compatible_cds_events$hits[[1]],roi)
    # }
    #print(fpkm.vals3[fpkm.vals3$nearest_ref_id %in% names(compatible_events$hits),])

    # ans<-eventOutcomeCompare(compatible_cds_events$hits[[1]],newregion,bsgenome,direction=loop_frame$data[i,]$IncLevelDifference>0)
    #loop_frame$data[i,]$AALENGTH1<-paste(eventOutcomeTranslate(compatible_cds_events,genome=bsgenome,fullseq = FALSE)[[1]],collapse=",")
    #loop_frame$data[i,]$AALENGTH2<-paste(eventOutcomeTranslate(compatible_cds_events,genome=bsgenome,fullseq = FALSE)[[2]],collapse=",")
  }

  columns_to_take<-c("seq_id","ID","Symbol","uID",
                     "PValue","FDR","IncLevel1","IncLevel2","IncLevelDifference",
                     "VALID_TX","VALID_CDS",
                     "COMPATIBLE_TX","COMPATIBLE_CDS",
                     "ASoutcome1","ASoutcome2")
  loop_frame$data<-loop_frame$data[,columns_to_take]
  if(full_output){
    ID <- uID <- NULL # Setting the variables to NULL first
    loop_frame$data[p_file_name.ls!="",]<-transform(loop_frame$data[p_file_name.ls!="",], uID = paste('<a href = ', shQuote(p_file_name.ls[p_file_name.ls!=""]), '>', uID, '</a>',sep=''))
    loop_frame$data <- transform(loop_frame$data, ID = paste('<a href = ', shQuote(r_file_name.ls), '>', ID, '</a>',sep=''))
    x = gvisTable(loop_frame$data, options = list(allowHTML = TRUE))
    print(x,file=paste(file.path(getwd(), outputdir),"/report.html",sep=""))
  }
  sink(type="message")
  closeAllConnections()
  return(loop_frame)
}

#' extractSpliceEvents
#'
#' Extracts the location of target, upstream and downstream splice sites
#' Used for calculations and genome visualizations
#'
#' @param data character. path to file
#' @param filetype character. type of splicing output. c('mats','custom'). see Details.
#' @param splicetype character. c('SE', 'RI', 'MXE', 'A5SS', 'A3SS')
#' @param fdr numeric. false discovery rate filter range [0,1]
#' @param inclusion numeric. splicing inclusion range, takes absolute value
#' @param start0 boolean 0-base start
#'
#' @return list containing information on\cr
#'   (1) original file type\cr
#'   (2) splice event type\cr
#'   (3) data.frame with splicing regions
#'
#' @details filetype 'custom' should provide a 9-column tab-delimited text file
#'    with the following columns:
#'    ID (Ensembl gene id), Symbol (gene name), chr, strand,
#'    exonStart, exonEnd, exon2Start, exon2End, upstreamStart, upstreamEnd, downstreamStart, downstreamEnd
#'    eg. ENSG0000012345 chr1 + 3 4 5 6 1 2 7 8
#'
#' @seealso \url{http://rnaseq-mats.sourceforge.net/user_guide.htm} for MATS file definition
#' @author Diana Low
#' @export
#'
#' @examples
#' data_path<-system.file("extdata",package="SPLINTER")
#' splice_data<-extractSpliceEvents(data=paste(data_path,"/skipped_exons.txt",sep=""))
extractSpliceEvents <- function(data=NULL,filetype='mats',splicetype='SE',fdr=1,
                                inclusion=1,start0=TRUE){
  df<-read.table(data,header=TRUE,stringsAsFactors=FALSE,row.names=1)
  if(filetype=='mats'){
    if(splicetype=='MXE'){
      colnames(df)[1:12]<-c("ID","Symbol","chr","strand","exonStart","exonEnd","exon2Start","exon2End","upstreamStart","upstreamEnd","downstreamStart","downstreamEnd")
      df[,c(5,7,9,11)]<-df[,c(5,7,9,11)]+1
    } else {
      colnames(df)[1:10]<-c("ID","Symbol","chr","strand","exonStart","exonEnd","upstreamStart","upstreamEnd","downstreamStart","downstreamEnd")
      df[,c(5,7,9)]<-df[,c(5,7,9)]+1
    }
    #filter FDR
    df<-df[df$FDR<=fdr,]
    #filter range
    if(inclusion!=1) df<-df[abs(df$IncLevelDifference)>=inclusion,]
  } else if(filetype=='custom'){
    colnames(df)[1:12]<-c("ID","Symbol","chr","strand","exonStart","exonEnd","exon2Start","exon2End","upstreamStart","upstreamEnd","downstreamStart","downstreamEnd")
    if(start0){
      df[,c(5,7,9,11)]<-df[,c(5,7,9,11)]+1
    }
  } else {
    message('Invalid filetype!')
  }


  return(list(filetype=filetype,type=splicetype,data=df))
}

#' extractSpliceSites
#'
#' Extracts and formats to bed the location of target, upstream and downstream
#' splice sites
#'
#' @param df extractSpliceEvents object
#' @param target the target site to extract. See Details.
#' @param site character donor or acceptor
#' @param motif_range numeric vector of splice position to extract
#' @param start0 boolean 0-base start
#'
#' @return GRanges object
#'
#' @details target : the site to extract the sequence from. It can be either the
#'          event in question (SE, RI, MXE - first exon, MXE2 - second exon,
#'          A5SSlong, A5SSshort, A3SSlong, A3SSshort, upstream or downstream).
#'          If this function is used in conjunction with \link{shapiroDonor} or
#'          \link{shapiroAcceptor} to compute scores, then most likely it will be run
#'          twice - once for the event, and the other either up- or downstream
#'          as a comparison.
#'
#' @seealso \url{http://rnaseq-mats.sourceforge.net/user_guide.htm} for MATS file definition
#' @author Diana Low
#' @import GenomicRanges
#' @export
#'
#' @examples
#' data_path<-system.file("extdata",package="SPLINTER")
#' splice_data<-extractSpliceEvents(data=paste(data_path,"/skipped_exons.txt",sep=""))
#' splice_sites<-extractSpliceSites(splice_data,target="SE")
extractSpliceSites <- function(df,target="SE",site='donor',
                               motif_range=c(-3,6),start0=TRUE){
  df<-df$data
  if(start0) base2add=1 else base2add=0

  switch(target,
         "SE" = {
           startf<-"exonStart"
           endf<-"exonEnd"
         },
         "RI" = {
           startf<-"exonStart"
           endf<-"exonEnd"
         },
         "MXE" = {
           startf<-"exonStart"
           endf="exonEnd"
         },
         "MXE2" = {
           startf<-"exon2Start"
           endf="exon2End"
         },
         "A5SSlong" = {
           startf<-"exonStart"
           endf="exonEnd"
         },
         "A5SSshort" = {
           startf<-"exon2Start"
           endf="exon2End"
         },
         "A3SSlong" = {
           startf<-"exonStart"
           endf="exonEnd"
         },
         "A3SSshort" = {
           startf<-"exon2Start"
           endf="exon2End"
         },
         "upstream"={
           startf<-"upstreamStart"
           endf<-"upstreamEnd"
         },
         "downstream" = {
           startf<-"downstreamStart"
           endf<-"downstreamEnd"
         }
  )

  target_exons<-makeGRangesFromDataFrame(df,start.field = startf , end.field = endf, starts.in.df.are.0based=start0)
  mcols(target_exons)$names<-df$ID
  target_exons_pos<-target_exons[strand(target_exons)=="+"]
  target_exons_neg<-target_exons[strand(target_exons)=="-"]

  if(site=='donor'){
    motif_range=motif_range+1
    start(target_exons_pos)<-end(target_exons_pos)+motif_range[1]
    end(target_exons_pos)<-end(target_exons_pos)+motif_range[2]
    end(target_exons_neg)<-start(target_exons_neg)-motif_range[1]-base2add
    start(target_exons_neg)<-start(target_exons_neg)-motif_range[2]-base2add
    target_sites<-c(target_exons_pos,target_exons_neg)
  } else if(site=='acceptor'){
    motif_range=motif_range-1
    end(target_exons_pos)<-start(target_exons_pos)+motif_range[2]-base2add
    start(target_exons_pos)<-start(target_exons_pos)+motif_range[1]-base2add
    start(target_exons_neg)<-end(target_exons_neg)-motif_range[2]
    end(target_exons_neg)<-end(target_exons_neg)-motif_range[1]
    target_sites<-c(target_exons_pos,target_exons_neg)
  }

  return(target_sites)
}

#' plotting sequence logo
#'
#' Plots the sequence logo of a given set of FASTA sequences
#'
#' @param fasta_seq DNAStringSet or path to fasta-formatted file
#' @return sequence logo image
#'
#' @author Diana Low
#' @import seqLogo Biostrings
#' @export
#'
#' @examples
#' head(splice_fasta)
#' test<-Biostrings::DNAStringSet(splice_fasta$V2)
#' plot_seqlogo(test)
plot_seqlogo <-function(fasta_seq){
  if(class(fasta_seq)!="DNAStringSet")
    fasta_seq<-suppressWarnings(readDNAStringSet(fasta_seq))
  freq<-consensusMatrix(fasta_seq,as.prob=TRUE)[1:4,]
  freq<-data.frame(freq)
  seqLogo(makePWM(freq),ic.scale=FALSE) #ic.scale determines either frequency or bits
}

#' shapiroDonor
#'
#' Shapiro and Senapathy (1987) have developed a method to score the strength of a splice site
#' based on percentages of each nucleotide at each position.
#' Shapiro's score of donor site (range is from -3 [exon] to +7 [intron]) is :
#' 100 * (t - min)/ (max - min), where
#' t is the sum of percentages at positions -3 to +7,
#' min is the sum of the lowest percentages at positions -3 to +7, and
#' max is the sum of the highest percentages at positions -3 to +7.
#' @param target_fasta vector of strings or DNAStringSet of fasta to score
#' @param reference_fasta vector of strings or DNAStringSet of reference splice list
#' @return data.frame with Shapiro scores
#' @seealso \url{http://www.softberry.com/spldb/SpliceDB.html}
#' @import Biostrings
#'
#' @author Diana Low
#' @export
#'
#' @author Diana Low
#'
#' @examples
#' library(BSgenome.Mmusculus.UCSC.mm9)
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm9
#' data_path<-system.file("extdata",package="SPLINTER")
#' splice_data<-extractSpliceEvents(data=paste(data_path,"/skipped_exons.txt",sep=""))
#' splice_sites<-extractSpliceSites(splice_data)
#' donor.ss<-getSeq(bsgenome,splice_sites)
#' sdonor<-shapiroDonor(donor.m,donor.ss)
shapiroDonor<-function(reference_fasta,target_fasta){
  if(class(reference_fasta)!="matrix") {
    rpwm<-as.matrix(read.table(reference_fasta,row.names=1))/100#consensusMatrix(reference_fasta,as.prob=T)[1:4,]
  } else {
    rpwm <-reference_fasta/100
  }
  minscore<-sum(apply(rpwm,2,min))
  maxscore<-sum(apply(rpwm,2,max))

  if(class(target_fasta)!="DNAStringSet") target_fasta<-suppressWarnings(readDNAStringSet(target_fasta))

  result<-unlist(lapply(target_fasta,function(x) {
    t<-PWMscoreStartingAt(rpwm, as.character(x), starting.at=1)
    score<-(t-minscore)/(maxscore-minscore)*100
  }))
  names(result)<-as.character(target_fasta)

  return(result)
}

#' shapiroAcceptor
#'
#' Shapiro's score of acceptor site (range is from -13 [intron] to +1 [exon]) is:
#' 100 * ((t1 - l1)/(h1 - l1) + (t2 - l2)/(h2 - l2))/2, where
#' t1 is the sum of the best 8 of 10 percentages at positions -13 to -4,
#' l1 is the sum of the lowest 8 of 10 percentages at position -13 to -4,
#' h1 is the sum of the highest 8 of 10 percentages at positions -13 to -4,
#' t2 is the sum of percentages at positions -3 to +1,
#' l2 is the sum of the lowest percentages at positions -3 to +1, and
#' h2 is the sum of the highest percentages at positions -3 to +1
#'
#' @param target_fasta vector of strings or DNAStringSet of fasta to score
#' @param reference_fasta vector of strings or DNAStringSet of reference splice list
#' @return data.frame with Shapiro scores
#' @seealso \url{http://www.softberry.com/spldb/SpliceDB.html}
#' @import Biostrings
#'
#' @author Diana Low
#' @export
#'
#' @examples
#' library(BSgenome.Mmusculus.UCSC.mm9)
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm9
#' data_path<-system.file("extdata",package="SPLINTER")
#' splice_data<-extractSpliceEvents(data=paste(data_path,"/skipped_exons.txt",sep=""))
#' splice_sites<-extractSpliceSites(splice_data,site="acceptor")
#' acceptor.ss<-getSeq(bsgenome,splice_sites)
#' sacceptor<-shapiroAcceptor(acceptor.m,acceptor.ss)
shapiroAcceptor<-function(reference_fasta,target_fasta){
  #if(class(reference_fasta)!="DNAStringSet") reference_fasta<-suppressWarnings(readDNAStringSet(reference_fasta))
  if(class(target_fasta)!="DNAStringSet") target_fasta<-suppressWarnings(readDNAStringSet(target_fasta))

  #consensusMatrix(reference_fasta,as.prob=T)[1:4,]
  if(class(reference_fasta)!="matrix") {
    rpwm<-as.matrix(read.table(reference_fasta,row.names=1))/100
  } else {
    rpwm <-reference_fasta/100
  }

  #rpwm<-consensusMatrix(reference_fasta,as.prob=T)[1:4,]
  tpwm<-consensusMatrix(target_fasta,as.prob=TRUE)[1:4,]

  l1<-sum(sort(apply(rpwm[,1:10],2,min))[1:8])
  h1<-sum(sort(apply(rpwm[,1:10],2,max),decreasing=TRUE)[1:8])
  l2<-sum(apply(rpwm[,(ncol(rpwm)-3):ncol(rpwm)],2,min))
  h2<-sum(apply(rpwm[,(ncol(rpwm)-3):ncol(rpwm)],2,max))

  t1<-sum(sort(apply(tpwm[,1:10],2,max),decreasing=TRUE)[1:8])

  result<-unlist(lapply(target_fasta,function(x) {
    t2<-PWMscoreStartingAt(tpwm[,(ncol(tpwm)-3):ncol(tpwm)], substr(x,(ncol(tpwm)-3),ncol(tpwm)), starting.at=1)
    score<-100*((t1 - l1)/(h1 - l1) + (t2 - l2)/(h2 - l2))/2
  }))

  return(result)
}

#' shapiroDensity
#'
#' convenience function for plotting Shapiro score density
#'
#' @param ctrl_scores output of shapiroDonor or shapiroAcceptor
#' @param treat_scores output of shapiroDonor or shapiroAcceptor
#' @param sample samplenames
#' @return density plot of Shapiro scores
#' @author Diana Low
#'
#' @import ggplot2
#' @export
shapiroDensity<-function(ctrl_scores,treat_scores,sample=c(1,2)){
  dens<-NULL
  dat <- data.frame(dens = c(ctrl_scores,treat_scores), feature = c(rep(paste("1:",sample[1]),times=length(ctrl_scores)), rep(paste("2:",sample[2]),times=length(treat_scores))))
  plt<-ggplot(dat, aes(x = dens, fill = feature))+ geom_density(alpha = 0.9,size=1,adjust=2)+ scale_x_continuous(limits = c(50, 100))+
    xlab("Shapiro score") +  ylab("Density")+scale_fill_brewer(palette="Set1")+theme_bw(base_size=18)
  return(plt)
}

#' makeUniqueIDs
#'
#' Makes unique ID names from event location
#'
#' @param ddata extractSpliceEvents object
#'
#' @return original extractSpliceEvents list object with unique ID
#'    appended to data accessor
#'
#' @import stringr
#' @export
#'
#' @author Diana Low
#'
#' @examples
#' data_with_id<-makeUniqueIDs(splice_data)
makeUniqueIDs<-function(ddata){
  if(ddata$type=="MXE") numcol=9 else numcol=7

  uID<-apply(ddata$data,1,function(x) paste(x[3:(3+numcol)],collapse="_",sep=""))
  uID<-str_replace_all(uID, fixed(" "), "")
  ddata$data<-cbind.data.frame(ddata$data,uID,stringsAsFactors=FALSE)

  return(ddata)
}

#' addEnsemblAnnotation
#'
#' Adds annotation to \code{\link{extractSpliceEvents}} object (if not present)
#'
#' @param data \code{\link{extractSpliceEvents}} object
#' @param species character. biomaRt species passed to retrieve annotation.
#'   Common species include: 'hsapiens','mmusculus'
#'
#' @return \code{\link{extractSpliceEvents}} object with annotated genes under $geneSymbol
#'
#' @import biomaRt
#' @export
#'
#' @author Diana Low
#' @seealso \url{http://asia.ensembl.org/info/data/biomart/biomart_r_package.html#biomartexamples}
#' @examples
#' data_path<-system.file("extdata",package="SPLINTER")
#' splice_data<-extractSpliceEvents(data=paste(data_path,"/skipped_exons.txt",sep=""))
#' splice_data<-addEnsemblAnnotation(data=splice_data,species="mmusculus")
addEnsemblAnnotation<-function(data,species='hsapiens'){
  df<-data$data
  ensembl<-useMart("ensembl",dataset=paste(species,"_gene_ensembl",sep=""))
  values<-df$ID
  tt<-getBM(attributes=c('ensembl_gene_id','wikigene_name'),
            filters = 'ensembl_gene_id', values = values, mart = ensembl)
  tt <- tt[!duplicated(tt$ensembl_gene_id),]
  colnames(tt)[1]<-"ID"
  tt<-merge(df,tt,by="ID",sort=FALSE,incomparables="-",all.x=TRUE)
  tt$Symbol<-tt$wikigene_name
  tt<-tt[,-ncol(tt)]
  tt[is.na(tt$Symbol),2] <- "-"
  data$data<-tt
  return(data)
}

#' findTX
#'
#' Given an ENSEMBL id, find all transcripts that matches id
#'
#' @param id character. transcript identification (currently ENSEMBL gene names)
#' @param db TxDb object
#' @param tx GRangesList
#' @param valid logical. check if in multiples of 3 [TRUE] for CDS translation.
#' @param verbose logical. turn messages on/off.
#'
#' @return GRangesList
#'
#' @import GenomicRanges
#' @export
#'
#' @author Diana Low
#'
#' @examples
#' valid_cds <-findTX(id=splice_data$data[2,]$ID,tx=thecds,db=txdb,valid=FALSE)
findTX <- function(id,db,tx,valid=FALSE,verbose=FALSE){
  options(warn=-1)
  txid <- suppressMessages(try(select(db,keys=id,columns="TXNAME","GENEID")[["TXNAME"]],TRUE))
  if(class(txid)!="try-error") {
    tx<-tx[names(tx) %in% txid]
    if(valid){
      tx_width<-width(tx)
      tx <- tx[(sum(tx_width) %% 3==0)]
    }
  } else {
    #return empty GRanges
    tx=GRanges()
  }

  if(verbose){
    if(length(tx)>0) {
      message(paste(names(tx),"\n"))
      message(length(tx)," valid transcripts found.")
    } else {
      message("No valid transcripts found.")
    }
  }

  return(tx)
  options(warn=1)
}

#' findCompatibleExon
#'
#' Finds compatible exon in annotation with the one present in roi object
#'
#' @param tx GRangesList object of transcripts
#' @param roi \code{\link{makeROI}} object contaning event information
#' @param verbose logical. printouts and messages.
#'
#' @return list of length 3\cr
#'  (1) GRangesList hits\cr
#'  (2) Number of transcripts\cr
#'  (3) Original number of input transcripts
#'
#' @import GenomicRanges
#' @importFrom S4Vectors queryHits subjectHits elementNROWS
#' @export
#'
#' @author Diana Low
#'
#' @examples
#' compatible_exons <- findCompatibleExon(valid_cds,roi)
findCompatibleExon <-function(tx,roi,verbose=FALSE){
  fo <- findOverlaps(roi$roi,tx)
  ovreg <- tx[subjectHits(fo)]

  if(length(subjectHits(fo))>0){
    exons_of_interest<-lapply(ovreg, function(x) queryHits(findOverlaps(x,roi$roi)))
    ovreg<-ovreg[lapply(ovreg[exons_of_interest], function(x) width(x)==width(roi$roi))]
    ovreg<-ovreg[elementNROWS(ovreg) != 0]
  }

  if(verbose){
    message(length(ovreg)," match(es) from original ",length(tx)," transcripts.")
    message(paste(names(ovreg),"\n"))
  }

  return(list(hits=ovreg,ct=length(ovreg),tt=length(tx)))
}

#' findCompatibleEvents
#'
#' Which transcript contains the event?
#' Each event has 2 possibilities, as long as the transcript fulfills one,
#' it passes the test
#' Has to be exact (inner junctions)
#'
#' Seperates into event/region1 and 2 for the alternative case
#'
#' @param tx GRangesList object of transcripts
#' @param tx2 optional GRangesList object of transcripts if tx is list of cds
#' @param roi \code{\link{makeROI}} object contaning event information
#' @param sequential logical. Exons have to appear sequentially to
#'  be considered compatible
#' @param verbose logical. printouts and messages.
#'
#' @return list of length 4\cr
#'  (1) GRangesList\cr
#'  (2) Hits status [c]=coding; [nc]=non-coding\cr
#'  (3) ct - compatible transcripts\cr
#'  (4) tt - total transcripts
#'
#' @importFrom S4Vectors queryHits subjectHits
#' @export
#'
#' @author Diana Low
#'
#' @examples
#' compatible_cds <- findCompatibleEvents(valid_cds,roi=roi,verbose=TRUE)
findCompatibleEvents<-function(tx,tx2=NULL,roi,sequential=TRUE,verbose=FALSE){
  qualified_status_total=list()
  hits<-list()
  total_names<-names(tx)

  compatible_query<-findExactOverlaps(roi$roi_range,tx,sequential=sequential)
  len1<-length(compatible_query)
  len2<-0

  #additional comparison, usually when tx is cds
  if(!is.null(tx2)){
    compatible_query2<-findExactOverlaps(roi$roi_range,tx2,sequential=sequential)
    len2<-length(compatible_query2)
    status1<-"[c]"
  } else {
    status1<-""
  }

  if(len1+len2!=0){
    for(r in seq_along(roi$roi_range)){
      qualified_names=c()
      qualified_status=c()
      if(verbose) message("Checking Type ",r,".....")
      compatible_subject_names<-names(tx[subjectHits(compatible_query[queryHits(compatible_query)==r])])
      qualified_names<-c(qualified_names,compatible_subject_names)
      qualified_status<-rep(status1,length(qualified_names))

      if(!is.null(tx2)){
        compatible_subject_names2<-total_names[total_names %in% names(tx2[subjectHits(compatible_query2[queryHits(compatible_query2)==r])])]
        compatible_subject_names2<-compatible_subject_names2[!(compatible_subject_names2 %in% qualified_names)]
        qualified_names<-c(qualified_names,compatible_subject_names2)
        qualified_status<-c(qualified_status,rep("[nc]",length(compatible_subject_names2)))
      }

      qualified_status_total[[r]]<-qualified_status
      hits[[r]]<-tx[qualified_names]
      if(verbose) {
        if(length(qualified_names)==0) message("No transcripts found!")
        else message(paste(qualified_names,"\n"))
      }
    }
  } else {
    if(verbose) message("No compatible transcripts found!")
  }
  match_hits<-sum(unlist(lapply(hits,function(x) length(x))))
  if(verbose) message(match_hits," match(es) from original ",length(tx)," transcripts.")

  return(list(hits=hits,hits_status=qualified_status_total,ct=match_hits,tt=length(tx)))
}

#' findExactOverlaps
#'
#' Internal function similar to findSpliceOverlaps but only preserves internal
#' flanks
#'
#' @param query GRanges object
#' @param subject GRanges object
#' @param sequential logical. TRUE if exons are sequential.
#' @param verbose logical. report intermediate output
#'
#' @return Hits object
#'
#' @import GenomicRanges
#' @keywords internal
#'
#' @author Diana Low
findExactOverlaps<-function(query,subject,sequential=FALSE,verbose=FALSE){
  #query is roi
  #subject is tx

  delete_rows<-c()

  # finds general overlap
  fo1<-findOverlaps(query,subject,maxgap=0)

  queryhits<-query[queryHits(fo1)]
  subjecthits<-subject[subjectHits(fo1)]

  # if there's a hit, find more exact
  if(length(fo1)!=0){
    fo <- as.data.frame(fo1)
    for(i in seq_along(fo1)){
      # will have exons numbers, able to tell if sequential
      exons_of_interest<-subjectHits(findOverlaps(queryhits[[i]],subjecthits[[i]],maxgap=0))

      # GRanges of exons_of_interest
      region_of_interest<-subjecthits[[i]][exons_of_interest]

      # criteria check
      if((!sequential & matchExons(region_of_interest,queryhits[[i]])) |
         (all(diff(exons_of_interest)==1) & matchExons(region_of_interest,queryhits[[i]]))){
        if(verbose) {
          message(paste(names(subject)[fo$subjectHits[i]],"qualifies"))
        }
      } else {
        delete_rows<-c(delete_rows,i)
      }
    }

    if(length(delete_rows!=0)) fo1<-fo1[-delete_rows]
  }
  return(fo1)
}

#' matchExons
#'
#' Internal function to help match the inner coordinates of a 2/3 cassette
#' checks if reference and subject matches
#'
#' @param ref GRanges object
#' @param subject GRanges object
#'
#' @return logical. check if exons match (TRUE) or not (FALSE)
#'
#' @keywords internal
#'
#' @author Diana Low
matchExons <- function(ref,subject){
  if(as.character(strand(subject)[1])=="-") {
    ref<-rev(ref)
    subject<-rev(subject)
  }
  ref_starts<-start(ref)
  ref_ends<-end(ref)
  sub_starts<-start(subject)
  sub_ends<-end(subject)
  pass=FALSE
  if(length(subject)!=length(ref)){
    pass=FALSE
  } else if(length(subject)==3){
       if(ref_starts[2]==sub_starts[2] &
          ref_starts[3]==sub_starts[3] &
          ref_ends[1]==sub_ends[1] &
          ref_ends[2]==sub_ends[2]) pass=TRUE
  } else if(length(subject)==2) {
    if(sub_starts[2]==ref_starts[2] & ref_ends[1]==sub_ends[1]) pass=TRUE
  } else if(length(subject)==1){
    if(sub_starts[1]==ref_starts[1] & ref_ends[1]==sub_ends[1]) pass=TRUE
  }
  return(pass)
}

#' makeROI
#'
#' Creates an object to store information about the splice site (region of interest)
#' including flanking regions and alternative splice outcome
#'
#' @param df data.frame object from \code{\link{extractSpliceEvents}}
#' @param type type of splicing event c("SE","RI","MXE","A5SS","A3SS")
#'
#' @return a list containing\cr
#'   (1) type : splice type\cr
#'   (2) name : ID of transcript\cr
#'   (3) roi  : GRanges object of splice site\cr
#'   (4) flank : GRanges object of flanking exons of splice site\cr
#'   (5) roi_range : GRangesList of splice site and it's alternative outcome based on \code{type}
#'
#' @import GenomicRanges IRanges
#' @export
#'
#' @author Diana Low
#'
#' @examples
#' single_record<-splice_data$data[which(grepl("Prmt5",splice_data$data$Symbol)),]
#' roi <- makeROI(single_record,type="SE")
makeROI <-function(df,type="SE"){
  colnum<-which(colnames(df)=='chr')
  chr<-df[,colnum]

  strand<-df[,colnum+1]

  if(type=="SE" | type=="RI"){
    chr_range<-as.numeric(df[,c((colnum+2):(colnum+7))])
    roi<-GenomicRanges::GRanges(chr,IRanges::IRanges(chr_range[1],chr_range[2]),strand=strand,exon_rank=as.integer(1))
    flank<-GRanges(c(chr,chr),c(IRanges(chr_range[3],chr_range[4]),IRanges(chr_range[5],chr_range[6])),exon_rank=1:2,strand=strand)
  }
  else {
    chr_range<-as.numeric(df[,c((colnum+2):(colnum+9))])
    roi<-GRanges(c(chr,chr),c(IRanges(chr_range[1],chr_range[2]),IRanges(chr_range[3],chr_range[4])),strand=strand,exon_rank=1:2)
    if(type=="MXE"){
      flank<-GRanges(c(chr,chr),c(IRanges(chr_range[5],chr_range[6]),IRanges(chr_range[7],chr_range[8])),exon_rank=1:2,strand=strand)
    }
    else {
      flank<-GRanges(c(chr),IRanges(chr_range[5],chr_range[6]),exon_rank=as.integer(1),strand=strand)
    }
  }

  if(strand=="-") flank<-rev(flank)
  roi_ranges<-GRangesList()

  if(type!="RI"){
    for(i in 1:length(roi)){
      roi_ranges[[length(roi_ranges)+1]]<-reduce(append(flank,roi[i],after=0))
      if(strand=="-") roi_ranges[[length(roi_ranges)]]<-rev(roi_ranges[[length(roi_ranges)]])
    }
  }

  if(type=="SE" | type=="RI"){
    roi_ranges[[length(roi_ranges)+1]]<-metaremove(flank)
    if(strand=="-") roi_ranges[[length(roi_ranges)]]<-rev(roi_ranges[[length(roi_ranges)]])
  }

  if(type=="RI"){
    roi_ranges[[length(roi_ranges)+1]]<-metaremove(roi)
    if(strand=="-") roi_ranges[[length(roi_ranges)]]<-rev(roi_ranges[[length(roi_ranges)]])
  }
  return(list(type=type,name=df$ID,roi=roi,flank=flank,roi_range=roi_ranges))
}

#' extendROI
#'
#' extend the span of the current ROI by n number of up/downstream exon(s)
#' by modifying roi_range within the makeROI object while
#' retaining legacy sites by keeping $roi and $flank
#'
#' @param roi \code{\link{makeROI}} object
#' @param tx GRangesList transcript list to pull regions from
#' @param up integer. number of exons to extend upstream
#' @param down integer. number of exons to extend downstream
#' @param type integer. 1=full cassette, 2=flank only
#'
#' @return \code{\link{makeROI}} object with modified ranges
#'
#' @import GenomicAlignments
#' @export
#'
#' @examples
#' extendROI(roi,valid_tx,up=1)
extendROI<-function(roi,tx,up=0,down=0,type=1){
  # 1. find transcripts that has the current event
  if(roi$type=="SE" && type==1){
    fso<-findSpliceOverlaps(GRangesList(roi$roi),tx)
    vfso<-tx[subjectHits(fso[mcols(fso)$compatible==TRUE])]
  } else if(roi$type=="SE" && type==2){
    fso<-findSpliceOverlaps(GRangesList(roi$roi_range[[2]]),tx)
    vfso<-tx[subjectHits(fso)]
  } else if(roi$type=="RI"){
    fso<-findSpliceOverlaps(GRangesList(roi$roi_range[[1]]),tx)
    vfso<-tx[subjectHits(fso[mcols(fso)$compatible==TRUE])]
  }

  # 2. get the transcripts
  #vfso<-tx[subjectHits(fso[mcols(fso)$compatible==TRUE])]

  # 3. get the exon numbers for the compatible transcripts
  for(i in seq_along(vfso)){
    evfso<-findOverlaps(vfso[[i]],roi$flank)
    if(length(evfso)==2) {
      break
    }
  }
  # 4. check if first/last exon

  qh<-queryHits(evfso)
  message("Current exon #",qh[1],", #",qh[2])

  if(!(up==0 & down==0)){
    if(qh[1]-up<1){
      message("Extending beyond upstream limits, keeping previous range.")
      up=0
    }
    if(qh[2]+down>length(vfso[[1]])){
      message("Extending beyond downstream limits, keeping previous range. (max ",length(vfso[[1]])-qh[2],", provided ", down,")")
      down=0
    }
  }


  # 5. extend range
  message("Extending exon to #",qh[1]-up,", #",qh[2]+down)
  roi$roi_range[[2]]<-metaremove(vfso[[1]][c(queryHits(evfso)[1]-up,queryHits(evfso)[2]+down)])
  roi$roi_range[[1]]<-reduce(c(metaremove(roi$roi),roi$roi_range[[2]]))

  ss<-as.character(strand(roi$roi))
  if(ss=="-"){
    roi$roi_range[[2]]<-rev(roi$roi_range[[2]])
    roi$roi_range[[1]]<-rev(roi$roi_range[[1]])
  }
  return(roi)
}

#' removeRegion
#'
#' removes a region (exon) from a GRanges or GRangesList
#'
#' @param subject GRanges or GrangesList object
#' @param roi \code{\link{makeROI}} object containing GRanges range (to remove)
#'
#' @return GRanges object
#'
#' @import GenomicRanges
#' @importFrom S4Vectors queryHits subjectHits
#' @export
#'
#' @author Diana Low
#'
#' # Removes the exon defined in roi GRanges object from a GRanges/GRangesList
#' compatible_cds$hits[[1]]
#' region_minus_exon<-removeRegion(compatible_cds$hits[[1]],roi)
removeRegion <- function(subject,roi){
  roi<-roi$roi
  for(tx in seq_along(subject)){
    #message(names(subject)[tx])
    outcome<-subjectHits(findOverlaps(roi,subject[[tx]]))
    if(length(outcome)!=0) {
      subject[[tx]]<-subject[[tx]][-outcome] #remove particular entry
      mcols(subject[[tx]],use.names=TRUE)$exon_rank <- 1:length(subject[[tx]]) #renumber exons
    }
  }
  return(subject)
}

#' insertRegion
#'
#' inserts a region (exon or intron) into roi object
#'
#' in the case of intron retention,
#' replaces exon with intron retention range
#' reduce() the GRanges in question
#'
#' @param subject GrangesList
#' @param roi \code{\link{makeROI}} object containg region of interest (to insert).
#'   refer to makeROI().
#'
#' @return GRanges object
#'
#' @import GenomicRanges
#' @importFrom S4Vectors queryHits subjectHits DataFrame
#' @export
#'
#' @author Diana Low
#'
#' @examples
#' #Inserts the exon defined in roi GRanges object from a GRanges/GRangesList
#' region_minus_exon
#' region_with_exon<-insertRegion(region_minus_exon,roi)
insertRegion <- function(subject,roi){
  roi<-roi$roi
  s2<-subject[,"exon_rank"]
  for(tx in seq_along(s2)){
    #message(names(s2)[tx])
    current<-s2[[tx]]
    current<-append(current,roi,after=0)
    current<-reduce(current)
    if(as.character(strand(current))=="-") current<-rev(current)
    elementMetadata(current) <- DataFrame(exon_rank=seq(1,length(current)))
    s2[[tx]]<-current
  }
  return(s2)
}

#' findTermination
#'
#' Internal function to find the first stop codon that occurs in the AA sequence, returns their position
#' and the resulting truncated protein
#'
#' @param s1 character. protein sequence
#'
#' @return list containing\cr
#'   (1) stop1 : stop position\cr
#'   (2) s1 : sequence truncated to first stop
#'
#' @import Biostrings
#' @keywords internal
#'
#' @author Diana LOW
#'
findTermination <-function(s1){

  findstop1 <-vmatchPattern("*",s1)

  if(length(unlist(findstop1))!=0) {
    truncate2stop1<-subseq(s1,1,startIndex(findstop1)[[1]][1])
    stop1<-startIndex(findstop1)[[1]][1]
  } else {
    truncate2stop1<-s1
    stop1<-0
  }

  return(list(stop1=stop1,s1=truncate2stop1))
}

#' eventOutcomeCompare
#'
#' Compares two sequences and gives differences if there's a switch from 1->2
#' if seq2 is NULL, assume seq1 is a list of length 2 to compare
#'
#' @param seq1 GRangesList
#' @param seq2 GRangesList
#' @param genome BSGenome object
#' @param direction logical. Report direction of sequence change.
#' @param fullseq logical. Report full sequences.
#' @param verbose logical. turn messages on/off.
#'
#' @return list containing\cr
#'   (1) tt : PairwiseAlignmentsSingleSubject pairwise alignment\cr
#'   (2) eventtypes : string detailing primary event classification\cr
#'
#' @import Biostrings GenomicFeatures GenomicRanges
#' @export
#' @author Diana LOW
#'
#' @examples
#' suppressMessages(library(BSgenome.Mmusculus.UCSC.mm9))
#' bsgenome<-BSgenome.Mmusculus.UCSC.mm9
#' eventOutcomeCompare(seq1=compatible_cds$hits[[1]],seq2=region_minus_exon,
#'  genome=bsgenome,direction=TRUE)
eventOutcomeCompare <- function(seq1,seq2=NULL,genome,direction=TRUE,fullseq=TRUE,verbose=FALSE){
  #gapOpening=-5; gapExtension=-1
  gapOpening=10; gapExtension=0.5

  if(is.null(seq2)) seq2<-seq1
  if(direction) esign="+" else esign="-"
  tresult=''
  options( warn = -1 )
  #submat <- matrix(-1, nrow = length(AA_ALPHABET), ncol = length(AA_ALPHABET), dimnames = list(AA_ALPHABET,AA_ALPHABET))
  #diag(submat) <- 1

  prot_seq1<-Biostrings::translate(extractTranscriptSeqs(genome,seq1))
  prot_seq2<-Biostrings::translate(extractTranscriptSeqs(genome,seq2))

  for(i in seq_along(prot_seq1)){
    message()
    message(paste("###",names(prot_seq1[i]),"###"))
    main_etype=""
    sub_etype=""

    FLAG_LEN<-FALSE #if sequences have been flipped in pairwiseAlignment because of length
    FLAG_INDEL<-FALSE

    wlastexon<-width(unlist(ranges(seq1[i])))[length(granges(unlist(seq1[i])))]
    #wlastexon2<-width(unlist(ranges(seq2[i])))[length(granges(unlist(seq2[i])))]
    ft<-findTermination(prot_seq1[i])
    ft2<-findTermination(prot_seq2[i])

    if(ft$stop1==0) message("Sequence1 - no stop codon")
    if(ft2$stop1==0) message("Sequence2 - no stop codon")

    if(ft$s1==ft2$s1){
      result='(no change)'
      message("No change in protein.")
    }
    else {
      if(ft2$stop1!=0){
        distanceupstream=(ft$stop1-ft2$stop1)*3

        #if negative difference, retained intron, new protein
        if(distanceupstream<0) {
          result='(ALT)'
          message("Alternative protein.")
          }


        #if positive difference, and more than 50bp from last exon
        else if(distanceupstream-wlastexon>=50) {
          result='(NMD)'
          message("Nonsense mediated decay.")
          }


        #if positive difference, but not more than 50, consider it truncation
        else if(distanceupstream-wlastexon<50)  {
          result='(TP)'
          message("Truncated protein.")
          }

      } else {
        result='(NTC)'
        message("No termination codon.")
      }
    }

    if(result=='(NMD)'){
      #don't really care about these
      main_etype="early termination"
      elocation="-"
      message(paste(elocation," ",main_etype))
    }

    if(width(ft$s1)<width(ft2$s1)) {
      FLAG_LEN<-TRUE
      epattern<-ft$s1;esubject<-ft2$s1
      epatternstop<-ft$stop1;esubjectstop<-ft$stop2
    }
    else{
      epattern<-ft2$s1;esubject<-ft$s1
      epatternstop<-ft2$stop1;esubjectstop<-ft$stop1
    }

    tt<-pairwiseAlignment(epattern,esubject,substitutionMatrix="BLOSUM62",gapOpening=gapOpening,gapExtension=gapExtension,type="global-local") #,gapOpening=-10,gapExtension=-1,

    inserts<-insertion(indel(tt))
    deletes<-deletion(indel(tt))
    nmm<-nmismatch(tt)

    if(length(inserts[[1]])==0 & length(deletes[[1]])==0 & nmm==0 & ft$s1==ft2$s1){
      result="(NC)"
      main_etype='No change, possibly a non-transcribed isoform'
      elocation=""
      message(elocation," ",main_etype,result)
    } else {
      #insertion
      if(length(inserts[[1]])!=0){
        sub_etype="insertion"
        for(id in seq_along(length(inserts[[1]]))){
          indelstart=start(inserts[[1]][id])
          indelend=end(inserts[[1]][id])
          if(indelend!=epatternstop) elocation="middle"
          else elocation="3' end"
          message(elocation," ",sub_etype," ",
                  subseq(aligned(pattern(tt)),indelstart,indelend)," (",indelstart,"-",indelend,")")
        }

      }

      #deletion
      if(length(deletes[[1]])!=0){
        sub_etype="deletion"
        for(id in 1:length(deletes[[1]])){
          indelstart=start(deletes[[1]][id])
          indelend=end(deletes[[1]][id])
          if(indelend!=epatternstop) elocation="middle"
          else elocation="3' end"
          message(elocation," ",sub_etype," ",
                  subseq(unaligned(subject(tt)),indelstart,indelend)," (",indelstart,"-",indelend,")")
        }
      }

      #mismatches
      if(nmm!=0){
        #print(mismatchTable(tt))
        #pairwise will truncate to the shorter sequence if mismatch is 3'
        sub_etype="mismatch"
        mmt<-mismatchTable(tt)
        if(!all(diff(mmt$PatternEnd)<=2)) message("multiple mismatch sites")
        #for display purposes
        if(mismatchTable(tt)$PatternEnd[1]==width(epattern)){
          ps1<-subseq(esubject,mismatchTable(tt)$SubjectStart[1],width(esubject))#paste(as.character(mismatchTable(tt)$SubjectSubstring),collapse="")
          ps2<-subseq(epattern,mismatchTable(tt)$PatternStart[1],width(epattern))#paste(as.character(mismatchTable(tt)$PatternSubstring),collapse="")
        } else {
          ps1<-subseq(esubject,mismatchTable(tt)$SubjectStart[1],mismatchTable(tt)$SubjectEnd[nrow(mmt)])#paste(as.character(mismatchTable(tt)$SubjectSubstring),collapse="")
          ps2<-subseq(epattern,mismatchTable(tt)$PatternStart[1],mismatchTable(tt)$PatternEnd[nrow(mmt)])#paste(as.character(mismatchTable(tt)$PatternSubstring),collapse="")
        }
        mmpend<-mmt$PatternEnd[nrow(mmt)]
        mmsend<-mmt$SubjectEnd[nrow(mmt)]

        if(mmpend!=width(epattern) ) { #| (mmt$PatternEnd[nrow(mmt)]!='*' & mmt$PatternEnd[nrow(mmt)]!='*')
          elocation="middle"
        }
        else {
          elocation="3' end"
          if(mmsend!=width(esubject)) ps1<-paste(ps1,subseq(esubject,mmsend+1,width(esubject)),sep="")
        }

        if(FLAG_LEN) message(elocation," ",sub_etype," : ", ps2," to ",ps1)
        else message(elocation," ",sub_etype," : ", ps1," to ",ps2)
      }

      if(length(insertion(indel(tt))[[1]])==0 & length(deletion(indel(tt))[[1]])==0){
        elocation="3' end"
        ps1<-subseq(esubject,width(epattern)+1,width(esubject))

        #(truncated protein)/(new peptide) depending on length
        if(width(ft$s1)>width(ft2$s1)) {
          sub_etype="truncation"
        }
        else {
          sub_etype="extra peptide"
        }
        message(elocation," ",sub_etype)
        message(ps1)
      }
    }

    message("length : ",paste(width(ft$s1)-1,"AA to",width(ft2$s1)-1,"AA"))
    if(tresult!='') {
      tresult<-paste(tresult,result,sep=",")
    } else {
      tresult<-paste(result,sep="")
    }

    if(fullseq){
      message("Original\n",prot_seq1[i],"\n")
      message("Post-event\n",prot_seq2[i],"\n")
    }
  }
  return(list(alignment=tt,eventtypes=tresult))
}

#' eventOutcomeTranslate
#'
#' translates sequences, reports if NMD or NTC
#'
#' @param seq1 GRangesList
#' @param genome BSGenome object
#' @param direction logical. Report direction of sequence change.
#' @param fullseq logical. Output full AA sequence.
#' @param verbose logical. turn messages on/off.
#'
#' @return list of translated sequences
#'
#' @import Biostrings BSgenome.Mmusculus.UCSC.mm9
#'
#' @export
#' @author Diana LOW
#'
#' @examples
#' suppressMessages(library(BSgenome.Mmusculus.UCSC.mm9))
#' bsgenome<-BSgenome.Mmusculus.UCSC.mm9
#' translation_results<-eventOutcomeTranslate(compatible_cds,genome=bsgenome,
#'  direction=TRUE)
eventOutcomeTranslate<-function(seq1,genome,direction=FALSE,fullseq=TRUE,verbose=FALSE){
  totalOutcomes<-list()
  gapOpening=-5
  gapExtension=-1
  if(direction) esign="+" else esign="-"
  tresult=''
  options( warn = -1 )
  submat <- matrix(-1, nrow = length(AA_ALPHABET), ncol = length(AA_ALPHABET),
                   dimnames = list(AA_ALPHABET,AA_ALPHABET))
  diag(submat) <- 1

  for(r in seq_along(seq1$hits)){
    outcome_status=seq1$hits_status[[r]]
    outcomes<-c()
    if(verbose) message("Translating for Type ",r)
    if(length(seq1$hits[[r]])==0) {
      if(verbose) message("No transcripts available.")
      outcomes<-"*"
    } else {
      prot_seq1<-Biostrings::translate(extractTranscriptSeqs(genome,seq1$hits[[r]]))
      if(verbose) print(prot_seq1)
      for(i in 1:length(seq1$hits[[r]])){
        thetranscript<-seq1$hits[[r]][[i]]
        FLAG_LEN<-FALSE #if sequences have been flipped in pairwiseAlignment because of length
        FLAG_INDEL<-FALSE
        wlastexon<-width(thetranscript)[length(thetranscript)]
        #wlastexon2<-width(unlist(ranges(seq2[i])))[length(granges(unlist(seq2[i])))]
        ft<-findTermination(prot_seq1[i])
        if(ft$stop1==0) {
          if(verbose) message("Sequence has no stop codon")
        } else if(ft$stop1!=0){
          distanceupstream=sum(width(thetranscript))-(ft$stop1)*3

          #if negative difference, retained intron, new protein
          if(distanceupstream==0) result=paste(width(ft$s1)-1,"AA")

          #if positive difference, and more than 50bp from last exon
          else if(distanceupstream-wlastexon>=50) result='NMD'

          #if positive difference, but not more than 50, consider it truncation
          else if(distanceupstream-wlastexon<50)  result='TP'
        } else {
          result='NTC'
        }
        if(verbose) message("###",names(prot_seq1[i]),":",result,outcome_status[i],"###")
        if(fullseq){
          message("Sequence\n",prot_seq1[i])
          message()
        }
        outcomes<-c(outcomes,paste(result,outcome_status[i]))
      }
    }
  totalOutcomes[[r]]<-outcomes
  }
  return(totalOutcomes)
}

#' eventPlot
#'
#' @param transcripts GRanges object
#' @param roi_plot GRanges object region to plot
#' @param bams character vector of bam file locations
#' @param names character vector of name labels
#' @param annoLabel character. annotation label
#' @param rspan integer or NULL. number of basepairs to span from roi.
#'    if NULL, will consider whole gene of roi
#' @param pfam_dom optional GRanges object of PFAM domains from UCSC Tables.
#' @param showAll logical. TRUE = display splice junctions of entire view or
#'    FALSE = just roi.
#'
#' @return a Gviz plot of genomic region
#'
#' @import GenomicRanges graphics Gviz
#' @export
#' @author Diana Low
#'
#' @examples
#' # define BAM files
#' data_path<-system.file("extdata",package="SPLINTER")
#' mt<-paste(data_path,"/mt_chr14.bam",sep="")
#' wt<-paste(data_path,"/wt_chr14.bam",sep="")
#'
#' # plot results
#' eventPlot(transcripts=valid_tx,roi_plot=roi,bams=c(wt,mt),
#'  names=c('wt','mt'),rspan=1000)
eventPlot <-function(transcripts,roi_plot=NULL,bams=c(),names=c(),annoLabel=c('Gene A'),rspan=1000,pfam_dom=NULL,showAll=TRUE){
  #options(Gviz.ucscUrl="http://genome-asia.ucsc.edu/cgi-bin/")
  if(!is.null(roi_plot)) {
    roi1=roi_plot$roi
    chr <- as.character(unique(seqnames(roi1)))[1]
    ##subsetting the transcript list
    transcripts<-subsetByOverlaps(transcripts,roi_plot$roi_range,ignore.strand=TRUE)

    ##subsetting the pfam list
    if(!is.null(pfam_dom)){
      pfam_dom<-subsetByOverlaps(pfam_dom,transcripts)
      if(length(pfam_dom)==0) pfam_dom=NULL
    }
  }
  else {
    roi1=NULL
    chr <- as.character(unique(seqnames(transcripts)))[1]
  }

  dataTrack<-as.data.frame(unlist(transcripts),row.names=NULL)
  dataTrack<-dataTrack[,1:5]
  dataTrack<-cbind(dataTrack,feature="protein_coding",exon=1,transcript=names(unlist(transcripts)))

  ### coloring the gene models
  comp_roi1<-names(findCompatibleEvents(transcripts,roi=roi_plot,verbose=FALSE)$hits[[1]])
  comp_roi2<-names(findCompatibleEvents(transcripts,roi=roi_plot,verbose=FALSE)$hits[[2]])
  dataTrack$compatibility<-"NA"
  if(!isEmpty(comp_roi1)) dataTrack[dataTrack$transcript %in% comp_roi1,]$compatibility<-"comp_roi1"
  if(!isEmpty(comp_roi2)) dataTrack[dataTrack$transcript %in% comp_roi2,]$compatibility<-"comp_roi2"

  #single copy tracks
  gen<-unique(genome(transcripts))

  gtrack <- GenomeAxisTrack(name="",cex=1.2,col="black")
  #itrack <- IdeogramTrack(genome = gen, chromosome = chr,showId=TRUE,showBandId=FALSE)
  grtrack <- GeneRegionTrack(dataTrack, genome = gen, cex=2,col=NULL, chromosome = chr,fill="darkgoldenrod2",
                             name = "Gene Model",transcriptAnnotation="transcript",fontcolor.exon = 1,
                             fontcolor.exon = 1,fontsize=16,fontsize.group=16,
                             fontcolor.group ="black",fontcolor.item="black",
                             feature=dataTrack$compatibility, featureAnnotation="feature",alpha=0.9)
  tracklist<-list(gtrack)#list(itrack,gtrack)#
  boxsizes<-c(1)#c(1,1)#



  #scaling to particular roi
  if(!is.null(roi1)) {
    if(roi_plot$type=="RI") roi_plot$roi_range<-rev(roi_plot$roi_range)
    anno<-AnnotationTrack(roi_plot$roi_range,genome=gen,shape="box",name=annoLabel,cex.group=1,cex=0.8,
                          fontcolor.feature="black",col=NULL,
                          fontcolor.exon = 1,fontsize=16,fontsize.group=16,
                          fontcolor.group ="black",fontcolor.item="black",
                          group=rep(c("Inclusion","Skipping"),c(length(roi_plot$roi_range[[1]]),length(roi_plot$roi_range[[2]]))),
                          fill=rep(c("#C4112C", "#1155C4"),c(length(roi_plot$roi_range[[1]]),length(roi_plot$roi_range[[2]]))))
    #group(anno)<-annoLabel
    tracklist<-c(tracklist,anno)
    #tracklist<-list(itrack,gtrack,anno,grtrack,alTrack1,alTrack2)
    boxsizes<-c(boxsizes,1)
    if(is.null(rspan)) rspan=50
    #tracklist<-list(anno,grtrack,alTrack1,alTrack2)
    #boxsizes<-c(1,1,2,2)
    startrange<-start(roi_plot$roi_range[[1]])
    endrange<-start(roi_plot$roi_range[[1]])

  } else {
    rspan=NULL
    #tracklist<-list(itrack,gtrack,grtrack,alTrack1,alTrack2)
    #boxsizes<-c(1,1,2,2,2)
  }

  tracklist<-c(tracklist,grtrack)
  boxsizes<-c(boxsizes,2)

  if(!is.null(pfam_dom)){
    pfam.transcripts<-pfam_dom
    pfam.dataTrack<-as.data.frame(pfam.transcripts,row.names=NULL)
    pfam.dataTrack<-pfam.dataTrack[,1:5]
    pfam.dataTrack<-cbind(pfam.dataTrack,feature="protein_coding",exon=1,transcript=pfam.transcripts$transcript_id)
    pfam_domains <- GeneRegionTrack(pfam.dataTrack,genome=gen,chromosome=chr,cex=2,col=NULL,
                                    transcriptAnnotation="transcript",feature="gene",fill="forestgreen",
                                    name="PFAM feature", fontcolor.exon = 1,fontsize=16,fontsize.group=16,
                                    fontcolor.group ="black",fontcolor.item="black")
    # pfam_domains <- UcscTrack(genome = gen, chromosome = chr,
    #                           track = "ucscGenePfam", from = min(start(roi_plot$roi_range[[1]])), to = max(start(roi_plot$roi_range[[1]])),
    #                           trackType = "AnnotationTrack", start = "chromStart",
    #                           end = "chromEnd", id = "name", shape = "box",
    #                           fill = "#006400", name = "PFAM",showFeatureId=TRUE,featureAnnotation = "group")
    #print(pfam_domains)
    tracklist<-c(tracklist,pfam_domains)
    boxsizes<-c(boxsizes,1)
  }
  #to accomodate multibams
  if(length(bams)!=0){
    for(i in 1:length(bams)){
      alTrack <- AlignmentsTrack(bams[i], isPaired = TRUE, name=names[i],col.axis="black",cex.axis=1.1,chromosome=chr)
      tracklist<-c(tracklist,alTrack)
      boxsizes<-c(boxsizes,2)
    }
  }

  #this should auto span in the case of roi being provided
  if(!is.null(rspan)) {
    if(showAll){
      #groupAnnotation="group",
      plotTracks(tracklist,
                 from=min(startrange)-rspan,to=max(endrange)+rspan,
                 cex.title=1,col.title="black",background.title = "lightgray",lwd=2,labelPos="below",sizes=boxsizes,
                 comp_roi1="#C4112C",comp_roi2="#1155C4")
    } else {
      introns<-psetdiff(range(roi_plot$roi_range),roi_plot$roi_range)
      concatintron<-append(introns[[1]],introns[[2]])
      plotTracks(tracklist,groupAnnotation="group",type=c("coverage","sashimi"),
                 sashimiFilter=concatintron,
                 from=min(startrange)-rspan,to=max(endrange)+rspan,
                 cex.title=1,col.title="black",background.title = "lightgray",lwd=2,labelPos="below",sizes=boxsizes,
                 comp_roi1="#C4112C",comp_roi2="#1155C4")
    }

  } else {
    plotTracks(tracklist,groupAnnotation="group",type=c("coverage","sashimi"),
               cex.title=1,col.title="black",background.title = "lightgray",lwd=2,labelPos="below",sizes=boxsizes,
               comp_roi1="#C4112C",comp_roi2="#1155C4")
  }
}

#' psiPlot
#'
#' Plots percentage spliced in (PSI) values in terms of inclusion levels
#'
#' @param df data.frame containing PSI values
#' @param type character. either 'MATS' output (will read in MATS headers) or
#'    'generic' (provide 4 or 6 column data.frame)
#' @param sample_labels x-axis labels for the plot
#'
#' @return bar plot of PSI values
#'
#' @import ggplot2
#' @export
#' @author Diana Low
#'
#' @examples
#' #we give inclusion and skipped numbers as reads
#' #this will be converted into percentages
#' df<-data.frame(inclusion1=c("6,4,6"),skipped1=c("10,12,12"),inclusion2=c("15,15,15"),
#'      skipped2=c("3,3,4"),stringsAsFactors = FALSE)
#' psiPlot(df,type='generic')
#'
psiPlot <- function (df=NULL,type="MATS",sample_labels=c("Sample 1", "Sample 2")){
  value<-NULL
  condition<-NULL

  if(type=="MATS"){
    # Check file type (junctions counts?)
    if(sum(grep("IJC",colnames(df)))!=0){
      inc1=mean(as.numeric(unlist(strsplit(as.character(df$IJC_SAMPLE_1),","))))
      sk1=mean(as.numeric(unlist(strsplit(as.character(df$SJC_SAMPLE_1),","))))
      inc2=mean(as.numeric(unlist(strsplit(as.character(df$IJC_SAMPLE_2),","))))
      sk2=mean(as.numeric(unlist(strsplit(as.character(df$SJC_SAMPLE_2),","))))
    } else {
      inc1=mean(as.numeric(unlist(strsplit(as.character(df$IC_SAMPLE_1),","))))
      sk1=mean(as.numeric(unlist(strsplit(as.character(df$SC_SAMPLE_1),","))))
      inc2=mean(as.numeric(unlist(strsplit(as.character(df$IC_SAMPLE_2),","))))
      sk2=mean(as.numeric(unlist(strsplit(as.character(df$SC_SAMPLE_2),","))))
    }
    inclength=as.numeric(df$IncFormLen)
    skplength=as.numeric(df$SkipFormLen)
  } else if(type=="generic"){
    inc1=mean(as.numeric(unlist(strsplit(df[,1],","))))
    sk1=mean(as.numeric(unlist(strsplit(df[,2],","))))
    inc2=mean(as.numeric(unlist(strsplit(df[,3],","))))
    sk2=mean(as.numeric(unlist(strsplit(df[,4],","))))
    if(ncol(df)>4){
      inclength=as.numeric(df[,5])
      skplength=as.numeric(df[,6])
    } else {
      inclength=1
      skplength=1
    }
  }
  inclusion1<-(inc1/inclength)/((inc1/inclength)+(sk1/skplength))
  inclusion2<-(inc2/inclength)/((inc2/inclength)+(sk2/skplength))
  skipped1<-1-inclusion1
  skipped2<-1-inclusion2

  df.p<-data.frame(condition=c(sample_labels[1],sample_labels[1],sample_labels[2],sample_labels[2]), type=as.factor(c("Inclusion","Skipped","Inclusion","Skipped")),value=c(inclusion1,skipped1,inclusion2,skipped2))

  ggplot()+
    geom_bar(aes(y=value,x=condition,fill=type),stat="identity",data=df.p) +
    scale_fill_manual(values=c("#C4112C", "#1155C4", "#C4112C","#1155C4")) +
    theme(axis.title.x = element_text(size=14, face="bold"),
          axis.title.y = element_text(size=14, face="bold"),
          text = element_text(size=14)) +
    labs(x="Condition",y="PSI")
}

#' metaremove
#'
#' helper function to remove metadata from GRanges object
#'
#' @param x GRanges or GRangesList
#'
#' @return GRanges or GRangesList
#'
#' @import IRanges S4Vectors
#' @importFrom methods is
#' @keywords internal
metaremove <- function(x) {
  if (is(x, "GRangesList")) return(endoapply(x, metaremove))
  reduce(remvalue(x))
}

#' remvalue
#'
#' helper function to remove metadata from GRanges object
#' used within metaremove
#'
#' @param x GRanges or GRangesList
#'
#' @return GRanges or GRangesList
#'
#' @keywords internal
remvalue <- function(x) {
  values(x) <- NULL
  x
}
