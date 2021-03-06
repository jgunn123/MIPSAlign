#' Read-in Fasta
#'
#' Allows you to read-in Fasta file as dataframe
#' @param file The Fasta file location
#' @return The dataframe with "name" and "sequence" columns
#' @examples
#' fa <- ReadFasta(getwd())
#' @export

ReadFasta<-function(file) {
  fasta <- ShortRead::readFasta(file)
  DF<-data.frame(name=as.character(ShortRead::id(fasta)),sequence=as.character(ShortRead::sread(fasta)))
  return(DF)
}

#' Write Fasta File
#'
#' Saves a dataframe as fasta file
#' @param data The dataframe with "name" and "sequence" columns that you want stored
#' @param filename The name of the file output
#' @return A Fasta file in appropriate format
#' @examples
#' fa <- WriteFasta(getwd())
#' @export

WriteFasta<-function(data, filename){
  seq = data$sequence
  names(seq) = data$name
  dna = Biostrings::DNAStringSet(seq)
  Biostrings::writeXStringSet(dna, filename)
}

#' Filter Fastq Files
#'
#' Performs a variety of filtering steps to clean Fastq files
#' @param fileDir The file directory with Fastq files, default is current working directory
#' @param trimTo Trims all reads to desired length (determined by the length of the longest sequence)
#' @param leftPattern The 3' adapter to trim sequences by
#' @param rightPattern The 5' adapter to trim sequences by
#' @param minLen All sequences below this minimum length requirement are filtered out
#' @param numCores The number of cores to run on
#' @examples
#' filterSeqs(trimTo = 30)
#' @export

filterSeqs <- function(fileDir = getwd(),
                       trimTo = 75,
                       leftPattern = "",
                       rightPattern = "",
                       minLen = 20,
                       numCores = 1) {

  startTime <- Sys.time()

  gzFiles <- list.files(path = fileDir,pattern = "*.fastq.gz$")

  newDir <- paste(fileDir,"filteredSeqs",sep = "/")
  R.utils::removeDirectory(newDir,recursive = TRUE,mustExist = FALSE)
  dir.create(newDir)

  if(numCores > parallel::detectCores()) numCores <- parallel::detectCores()

  closeAllConnections()
  cl <- parallel::makeCluster(numCores,outfile = paste(newDir,"filterSeqs.out",sep = "/"))
  doSNOW::registerDoSNOW(cl)
  on.exit(parallel::stopCluster(cl))

  iterations <- length(gzFiles)
  pb <- txtProgressBar(max = iterations,style = 3)
  progress <- function(n) setTxtProgressBar(pb,n)
  opts <- list(progress = progress)

  foreach::foreach(file = gzFiles,
                   .packages = c("ShortRead","FastqCleaner"),
                   .options.snow = opts) %dopar% {

    stream <- ShortRead::FastqStreamer(paste(fileDir,file,sep = "/"))
    newName <- paste(newDir,file,sep = "/")
    newName <- paste(sub(".fastq.gz","",newName),"trunc.fastq.gz",sep = "_")

    repeat {

      fq <- yield(stream)
      if(length(fq) == 0) break
      fq <- ShortRead::narrow(fq,end = trimTo)
      fq <- ShortRead::clean(fq)
      fq <- FastqCleaner::adapter_filter(input = fq,
                                         Lpattern = leftPattern,
                                         Rpattern = rightPattern,
                                         error_rate = 0,
                                         anchored = FALSE)
      fq <- fq[ShortRead::width(fq) >= minLen]
      ShortRead::writeFastq(fq,newName,"a")
    }
  }

  endTime <- Sys.time()
  print(endTime - startTime)
}

#' Convert Bowtie Output
#'
#' Convert bowtie or bowtie2 output to counts
#' @param df Output dataframe from an aligner in MIPSAlign
#' @return Two dataframes in a list that contain the aligned and unaligned reads, respectively
#' @examples
#' newdf <- btToCounts(df)
#' @export

btToCounts <- function(df) {

  f6 <- function(X) as.data.frame(data.table::data.table(X)[,.N,keyby = X])
  test <- df %>% dplyr::select(rname) %>% f6
  perf <- df %>% dplyr::filter(NM == 0 & flag != 16) %>% dplyr::select(rname) %>% f6
  imperf <- df %>% dplyr::filter(NM != 0 & flag != 16) %>% dplyr::select(rname) %>% f6
  unaligned <- df %>% dplyr::filter(flag == 4 | flag == 16) %>% dplyr::select(seq) %>% f6
  unaligned <- unaligned[order(unaligned$N, decreasing = TRUE),]

  total <- test %>% dplyr::full_join(y = perf,by = "rname") %>% dplyr::full_join(y = imperf,by = "rname")
  colnames(total) <- c("rname","total_counts","perf_counts","imperf_counts")

  dfs <- list(total,unaligned)
  return(dfs)
}

#' Perfect Match Aligner
#'
#' Generates counts on perfect sequence matches
#' @param seqData The list of Fastq files you want to align
#' @param orfLib The reference Fasta alignment file
#' @param numCores The number of cores to run on
#' @return A list of aligned and unaligned counts dataframes, respectively
#' @examples
#' x <- matchAlign(orfLib = fastaFileName,numCores = 6)
#' @export

matchAlign <- function(seqData = getwd(),
                       orfLib = list.files(pattern = "*.fa$|.fasta$",full.names = TRUE)[[1]],
                       numCores = 1) {

  startTime <- Sys.time()

  gzFiles <- list.files(path = seqData,pattern = "*.fastq.gz$",full.names = TRUE)

  fasta <- ReadFasta(orfLib)
  name <- sub(".*/","",gzFiles)
  name <- sub(".fastq.gz","",name)
  name <- sub("_trunc","",name)

  if(numCores > parallel::detectCores()) numCores <- parallel::detectCores()

  closeAllConnections()
  cl <- parallel::makeCluster(numCores,outfile = paste(seqData,"matchAlign.out",sep = "/"))
  doSNOW::registerDoSNOW(cl)
  on.exit(parallel::stopCluster(cl))

  iterations <- length(seqData)
  pb <- txtProgressBar(max = iterations,style = 3)
  progress <- function(n) setTxtProgressBar(pb,n)
  opts <- list(progress = progress)

  aligner = foreach::foreach(file = gzFiles,
                             .combine = cbind,
                             .packages = c("ShortRead","data.table"),
                             .options.snow = opts) %dopar% {

    seqs <- c()
    stream <- open(ShortRead::FastqStreamer(file))
    on.exit(close(stream))

    repeat {
      fq <- ShortRead::yield(stream)
      if(length(fq) == 0) break
      fq <- as.character(ShortRead::sread(fq))
      seqs <- c(seqs,fq)
    }

    seqs <- c(seqs,unique(fasta$sequence))
    alignment <- data.table::data.table(seqs)[,.N,keyby = seqs]
    unaligned <- alignment[!which(as.character(unlist(alignment[,1])) %in% as.character(unlist(fasta$sequence))),]
    sumUnalign <- sum(as.numeric(unaligned$N))
    unaligned <- unaligned[order(unaligned$N,decreasing = TRUE),]
    unaligned <- as.data.frame(head(unaligned,50))

    if(nrow(unaligned) < 50) {

      numRows <- nrow(unaligned)
      unaligned[(numRows+1):50,1] <- paste(file,NA,(numRows+1):50,sep = "_")
      unaligned[(numRows+1):50,2] <- 0
    }

    unaligned$N <- as.numeric(unaligned$N)
    alignment <- alignment[which(as.character(unlist(alignment[,1])) %in% as.character(unlist(fasta$sequence))),]
    alignment <- alignment[match(as.character(unlist(fasta$sequence)),as.character(unlist(alignment[,1]))),]
    aligned <- c(sumUnalign,as.numeric(unlist(alignment[,2])) - 1)

    return(list(as.data.frame(aligned),as.data.frame(unaligned)))
  }

  counts <- as.data.frame(aligner[1,])
  QC <- as.data.frame(aligner[2,1])

  for(i in 2:(ncol(aligner))) {

    df <- as.data.frame(aligner[2,i])
    colnames(df)[1] <- colnames(QC)[1]
    QC <- dplyr::full_join(QC,df,by = colnames(QC)[1])
  }

  colnames(counts) <- name
  counts$probes <- c("unaligned",as.character(fasta$name))
  counts <- dplyr::select(counts,probes,everything())

  colnames(QC) <- c("sequence",name)
  QC[is.na(QC)] <- 0

  endTime <- Sys.time()
  print(endTime - startTime)
  print("Finished aligning...")
  return(list(counts,QC))
}

#' Bowtie Aligner
#'
#' Generates counts using the original bowtie algorithm
#' @param seqData The  directory where Fastq files are located
#' @param orfLib The reference Fasta alignment file
#' @param numCores The number of cores to run on
#' @param keepFiles TRUE or FALSE, whether you want to keep the bowtie output files
#' @param options Extra options for bowtie that feed into "alignmentParameter"
#' @return A list of aligned and unaligned counts dataframes, respectively
#' @examples
#' x <- btAlign(orfLib = fastaFileName,numCores = 6)
#' @export

btAlign <- function(seqData = getwd(),
                    orfLib = list.files(pattern = "*.fa$|*.fasta$",full.names = TRUE)[[1]],
                    aligner = c("btAlign","SpliceMap"),
                    numCores = 1,
                    keepFiles = FALSE,
                    options = NULL) {

  startTime <- Sys.time()

  tsvFileName <- paste(seqData,"nameRef.tsv",sep = "/")
  gzFiles <- list.files(path = seqData,pattern = "*.fastq.gz$",full.names = TRUE)
  name <- sub(".*/","",gzFiles)
  name <- sub(".fastq.gz","",name)
  name <- sub("_trunc","",name)
  tsvFile <- data.frame(FileName = gzFiles)
  tsvFile$SampleName <- as.character(name)
  write.table(as.data.frame(tsvFile),tsvFileName,quote = FALSE, sep = "\t", row.names = FALSE)

  if(aligner != "SpliceMap") aligner <- "Rbowtie"

  if(numCores > 3) numCores <- 3

  closeAllConnections()
  cl <- parallel::makeCluster(numCores,outfile = paste(seqData,"btAlign.out",sep = "/"))
  doSNOW::registerDoSNOW(cl)
  on.exit(parallel::stopCluster(cl))

  QuasR::qAlign(sampleFile = tsvFileName,
                genome = orfLib,
                alignmentParameter = options,
                clObj = cl)

  outBam <- list.files(path = seqData,pattern = ".bam$",full.names = TRUE)
  p <- Rsamtools::ScanBamParam(what = c("qname","flag","rname","strand","qwidth","mapq","seq"),tag = "NM")
  seqs <- as.character(ShortRead::id(ShortRead::readFasta(orfLib)))
  bView <- list()

  for(i in 1:length(outBam)) {

    checkFlag <- as.data.frame(Rsamtools::scanBam(outBam[[i]],param = Rsamtools::ScanBamParam(what = c("qname","flag","seq"))))

    if(all(checkFlag$flag == 4)) {

      bViewNew <- data.frame(matrix(nrow = nrow(checkFlag),ncol = 8))
      colnames(bViewNew) <- c("qname","flag","rname","strand","qwidth","mapq","seq","NM")
      bViewNew$qname <- checkFlag$qname
      bViewNew$flag <- checkFlag$flag
      bViewNew$seq <- checkFlag$seq
      bView <- c(bView,list(bViewNew))
    }
    else {
      bView <- c(bView,list(as.data.frame(Rsamtools::scanBam(outBam[[i]],param = p))))
    }
  }

  if(keepFiles == FALSE) {

    unlink(grep(sub(".*/","",orfLib),list.dirs(path = paste(seqData,"..",sep = "/")),value = TRUE),recursive = T)
    txtFiles <- grep("QuasR",list.files(pattern = "*.txt$",full.names= TRUE),value = TRUE)
    bamFiles <- grep(paste(name,collapse = "|"),list.files(path = seqData,pattern = "*.bam",full.names = TRUE),value = TRUE)
    faiFiles <- grep(sub(".*/","",orfLib),list.files(path = paste(seqData,"..",sep = "/"),pattern = "*.fai$",full.names = TRUE),value = TRUE)
    md5Files <- grep(sub(".*/","",orfLib),list.files(path = paste(seqData,"..",sep = "/"),pattern = "*.md5$",full.names = TRUE),value = TRUE)
    file.remove(c(txtFiles,bamFiles,faiFiles,md5Files,tsvFileName))
  }

  endTime <- Sys.time()
  print(endTime - startTime)
  print("Finished aligning...")
  return(bView)
}

#' Bowtie2 Aligner
#'
#' Generates counts using the bowtie2 algorithm
#' @param seqData The  directory where Fastq files are located
#' @param orfLib The reference Fasta alignment file
#' @param numCores The number of cores to run on
#' @param keepFiles TRUE or FALSE, whether you want to keep the bowtie output files
#' @param options Extra options that feed into bowtie2 "options"
#' @return A list of aligned and unaligned counts dataframes, respectively
#' @examples
#' x <- bt2Align(orfLib = fastaFileName,numCores = 6)
#' @export

bt2Align <- function(seqData = getwd(),
                     orfLib = list.files(pattern = "*.fa$|*.fasta$",full.names = TRUE)[[1]],
                     numCores = 1,
                     keepFiles = FALSE,
                     options = "") {

  startTime <- Sys.time()

  gzFile <- FALSE
  gzFiles <- list.files(path = seqData,pattern = "*.fastq.gz$",full.names = TRUE)

  bView <- list()

  if(numCores > parallel::detectCores()) numCores <- parallel::detectCores()

  closeAllConnections()
  cl <- parallel::makeCluster(numCores,outfile = paste(seqData,"bt2Align.out",sep = "/"))
  doSNOW::registerDoSNOW(cl)
  on.exit(parallel::stopCluster(cl))

  iterations <- length(seqData)
  pb <- txtProgressBar(max = iterations,style = 3)
  progress <- function(n) setTxtProgressBar(pb,n)
  opts <- list(progress = progress)

  bt2Out <- foreach::foreach(file = gzFiles,
                             .combine = c,
                             .packages = c("ShortRead","Rbowtie2","Rsamtools","R.utils"),
                             .options.snow = opts) %dopar% {

    if(grepl(".gz$",file)) {

      stream <- open(ShortRead::FastqStreamer(file))
      on.exit(close(stream))

      fqFile <- sub(".gz$","",file)
      file.remove(fqFile)

      repeat {
        readIn <- yield(stream)
        if(length(readIn) == 0) break
        ShortRead::writeFastq(readIn,fqFile,"a",compress = FALSE)
      }

      file <- fqFile
      gzFile <- TRUE
    }
    else if(!grepl(".fastq$",file) & !grepl(".fq$",file) & !grepl("-f",options)) {

      print("Please check your sequencing file is in the FASTQ format, or specify in options as such!")
      return()
    }

    idx <- sub(".fastq","",file)
    idx <- sub("_trunc","",idx)
    outSam <- paste(idx,"alignment.sam",sep = "_")
    outBam <- paste(idx,"alignment.bam",sep = "_")

    Rbowtie2::bowtie2_build(references = orfLib,bt2Index = idx,overwrite = TRUE)
    btResults <- Rbowtie2::bowtie2(bt2Index = idx,
                                   samOutput = outSam,
                                   seq1 = file,
                                   options = paste(options,"-t"),
                                   overwrite = TRUE)

    Rsamtools::asBam(outSam,destination = sub(".bam$","",outBam),overwrite = TRUE)
    p <- Rsamtools::ScanBamParam(what = c("qname","flag","rname","strand","qwidth","mapq","seq"),tag = "NM")

    checkFlag <- as.data.frame(Rsamtools::scanBam(outBam,param = Rsamtools::ScanBamParam(what = c("qname","flag","seq"))))

    if(all(checkFlag$flag == 4)) {

      bViewNew <- data.frame(matrix(nrow = nrow(checkFlag),ncol = 8))
      colnames(bViewNew) <-c("qname","flag","rname","strand","qwidth","mapq","seq","NM")
      bViewNew$qname <- checkFlag$qname
      bViewNew$flag <- checkFlag$flag
      bViewNew$seq <- checkFlag$seq
      bView <- bViewNew
    }
    else {
      bView <- as.data.frame(Rsamtools::scanBam(outBam,param = p))
    }

    idx <- sub(".*/","",idx)

    if (keepFiles == FALSE) {

      if (gzFile == TRUE) file.remove(file)

      bt2Files <- grep(idx,list.files(path = seqData,pattern = ".bt2",full.names = TRUE),value = TRUE)
      samFile <- grep(idx,list.files(path = seqData,pattern = ".sam$",full.names = TRUE),value = TRUE)
      bamFiles <- grep(idx,list.files(path = seqData,pattern = ".bam",full.names = TRUE),value = TRUE)
      file.remove(c(bt2Files,samFile,bamFiles))
    }
    return(list(bView))
  }

  endTime <- Sys.time()
  print(endTime - startTime)
  print("Finished aligning...")
  return(bt2Out)
}

#' MIPSAlign Trimmer and Aligner
#'
#' Combines trimming and aligning in one comprehensive function
#' @param fileDir The file directory where all the Fastq files are located
#' @param refFile The reference Fasta alignment file
#' @param filter TRUE or FALSE, if you want to use the Fastq filter
#' @param leftPattern Relevant only if filter = TRUE, the 3' adapter to trim sequences by
#' @param rightPattern Relevant only if filter = TRUE, the 5' adapter to trim sequences by
#' @param minLen Relevant only if filter = TRUE, all sequences below this minimum length threshold are filtered out
#' @param trimTo Relevant only if filter = TRUE, trims all reads to desired length (determined by the length of the longest sequence)
#' @param aligner One of "matchAlign", "btAlign", or "bt2Align"
#' @param numCores The number of cores to run on
#' @param keepFiles TRUE or FALSE, whether you want to keep the bowtie output files
#' @param options Extra options for bowtie that feed into "alignmentParameter"
#' @return A list of counts, unaligned counts, and possibly imperfect counts if a bowtie aligner was used
#' @examples
#' x <- MIPSAlign(refFile = fastaFileName,filter = TRUE,trimTo = 30,numCores = 6)
#' @export

MIPSAlign <- function(fileDir = getwd(),
                      refFile = list.files(pattern = "*.fa$|*.fasta$",full.names = TRUE)[[1]],
                      filter = FALSE,
                      leftPattern = "",
                      rightPattern = "",
                      minLen = 20,
                      trimTo = 75,
                      aligner = c("matchAlign","btAlign","bt2Align"),
                      numCores = 1,
                      keepFiles = FALSE,
                      options = "") {

  startTime <- Sys.time()

  if(filter == TRUE) {

    filterSeqs(fileDir = fileDir,
               trimTo = trimTo,
               leftPattern = leftPattern,
               rightPattern = rightPattern,
               minLen = minLen,
               numCores = numCores)

    print("Finished filtering files...")

    fileDir <- paste(fileDir,"filteredSeqs",sep = "/")
  }

  if("matchAlign" %in% aligner) {

    output <- matchAlign(seqData = fileDir,orfLib = refFile,numCores = numCores)

    counts <- output[[1]]
    unaligned_seqs <- output[[2]]

    data.table::fwrite(counts,file = paste(fileDir,"counts.csv",sep = "/"),row.names = FALSE)
    data.table::fwrite(unaligned_seqs,file = paste(fileDir,"unaligned.csv",sep = "/"),row.names = FALSE)

    endTime <- Sys.time()
    print(endTime - startTime)
    print("All done!")
    return(output)
  }

  else if(aligner == "btAlign"|aligner == "SpliceMap") {

    btOut <- btAlign(seqData = fileDir,
                     orfLib = refFile,
                     aligner = aligner,
                     numCores = numCores,
                     keepFiles = keepFiles,
                     options =  options)
  }
  else {

    btOut <- bt2Align(seqData = fileDir,
                      orfLib = refFile,
                      numCores = numCores,
                      keepFiles = keepFiles,
                      options =  options)
  }

  gzFiles <- list.files(path = fileDir,pattern = ".fastq.gz$",full.names = TRUE)

  fasta_names <- ReadFasta(refFile)

  counts <- data.frame(matrix(ncol = 0, nrow = length(unique(fasta_names$name))))
  row.names(counts) <- unique(fasta_names$name)

  imperf_counts <- counts

  unaligned_seqs <- data.frame(seq = "")

  unaligned_counts <- data.frame(matrix(ncol = 1, nrow = 1))

  for(i in 1:length(btOut)) {

    countsOut <- btToCounts(btOut[[i]])
    alignOut <- countsOut[[1]]
    unalignedSeqs <- countsOut[[2]]
    alignOut$rname <- as.character(alignOut$rname)
    name <- sub(".*/","",gzFiles[[i]])
    name <- sub(".fastq.gz","",name)
    name <- sub("_trunc","",name)
    unaligned_counts[name] <- alignOut[1,2]

    total_imperf_counts <- dplyr::select(alignOut, rname, total_counts, imperf_counts)
    completeCounts <- counts %>% dplyr::mutate(rname = row.names(counts)) %>%
      dplyr::full_join(y = total_imperf_counts ,by = "rname")
    final_unalignedSeqs <- completeCounts[!is.na(completeCounts$rname),]
    counts$probes <- final_unalignedSeqs$rname
    counts[name] <- final_unalignedSeqs$total_counts
    imperf_counts$probes <- final_unalignedSeqs$rname
    imperf_counts[name] <- final_unalignedSeqs$imperf_counts

    names(unalignedSeqs)[names(unalignedSeqs)=="N"] <- name
    unalignedSeqs <- head(unalignedSeqs,50)
    unaligned_seqs <- unaligned_seqs %>% dplyr::full_join(y = unalignedSeqs,by = "seq")
  }

  #append unaligned reads to counts file
  colnames(unaligned_counts)[1] <- "probes"
  counts <- rbind(unaligned_counts, counts)
  row.names(counts)[1] <- "unaligned"

  #Replace NA's with 0
  counts[is.na(counts)] <- 0
  counts[1,1] <- "unaligned_reads"
  imperf_counts[is.na(imperf_counts)] <- 0

  unaligned_seqs <- unaligned_seqs[-1, ]
  unaligned_seqs[is.na(unaligned_seqs)] <- 0

  data.table::fwrite(counts,file = paste(fileDir,"counts.csv",sep = "/"),row.names = FALSE)
  data.table::fwrite(unaligned_seqs,file = paste(fileDir,"unaligned.csv",sep = "/"),row.names = FALSE)
  data.table::fwrite(imperf_counts,file = paste(fileDir,"imperfCounts.csv",sep = "/"),row.names = FALSE)

  endTime <- Sys.time()
  print(endTime - startTime)
  print("All done!")
  return(list(counts,unaligned_seqs,imperf_counts))
}
