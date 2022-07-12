#' Parsing FCS files
#' @param input_events_downsampling See \code{\link{infinity_flow}}
#' @param paths Character vector of paths to store intput, intermediary results, outputs...
#' @param extra_args_read_FCS passed to flowCore:read.FCS
#' @param name_of_PE_parameter Name of the exploratory measurement
#' @description This function reads the input FCS files, downsample to a user-defined number of events if necessary, harmonize the input data names, save a concatenated expression matrix and the corresponding vector mapping events to their file of origin.
#' @param verbose Verbosity
#' @noRd
#' @importFrom flowCore read.FCS
#' @importFrom Biobase pData exprs
subsample_data <- function(
                        input_events_downsampling,
                        vector_of_file_absolute_paths,
                        paths,
                        extra_args_read_FCS,
                        name_of_PE_parameter,
                        annotation,
                        verbose=TRUE
                        ){

    ## Subsampling
    if(verbose){
        message("Parsing and subsampling input data")
        message("\tDownsampling to ",input_events_downsampling," events per input file")
    }
  new_name <- names(annotation)
  files <- vector_of_file_absolute_paths
  create_file <- function(file, new_name){
        res <- do.call(read.FCS,c(list(filename=file),extra_args_read_FCS))
        w <- sort(sample(seq_len(nrow(res)),min(input_events_downsampling,nrow(res))))
        res <- res[w,]
        write.FCS(res,file.path(paths["subset"],new_name))
    }

    for (i in 1:(length(files))){
      create_file(files[i],new_name[i])
    }


    ## convert to .Rds
    if(verbose){
        message("\tConcatenating expression matrices")
    }
    files <- list.files(paths["subset"],full.names=TRUE,recursive=FALSE,include.dirs=FALSE,pattern=".fcs")
    #Crucial to sort the file names the same way they would be sorted by the OS
    #Data names will be protected this way
    print('files before sorting')
    print(files)
    files <- gtools::mixedsort(files)
    print('files after sorting')
    print(files)
    ns <- setNames(integer(length(names(annotation))),names(annotation))
    xp <- lapply(
        files,
        function(file){
            xp <- do.call(read.FCS,c(list(filename=file),extra_args_read_FCS))
            annot <- pData(xp@parameters)
            ## ns[file]<<-nrow(xp)
            xp <- exprs(xp)
            targets <- annot$desc
            targets[is.na(targets)] <- annot$name[is.na(targets)]
            colnames(xp)[colnames(xp)!=name_of_PE_parameter] <- targets[colnames(xp)!=name_of_PE_parameter]
            colnames(xp)[colnames(xp)==name_of_PE_parameter] <- name_of_PE_parameter
            xp
        }
    )
    names(xp) <- names(annotation)
    ns <- vapply(xp, nrow, 1L)
    xp <- do.call(rbind,xp)
    ## Map which events originate from which file.
    if(verbose){
        message("\tWriting to disk")
    }
    events.code <- unlist(
        lapply(
            names(ns),
            function(x){
                rep(tail(strsplit(x,"/")[[1]],1),ns[x])
            }
        )
    )
    saveRDS(xp,file=file.path(paths["rds"],"xp.Rds"))
    saveRDS(events.code,file=file.path(paths["rds"],"pe.Rds"))
    print('done script 1')
    invisible()
}
