#' read.AsspDataObj creates an object of class dobj from a signal or parameter 
#' file readable by the ASSP Library (WAVE, SSFF, AU, ...)
#'
#' @title read.AsspDataObj from a signal/parameter file
#' @param fname filename of the signal or parameter file
#' @param begin begin time (default is in seconds) of segment to retrieve
#' @param end end time (default is in seconds) of segment to retrieve
#' @param samples (BOOL) if set to false seconds values of begin/end are sample numbers
#' @return list object containing file data
#' @author Lasse Bombien
#' @aliases getAsspDataObj
#' @useDynLib superassp, .registration = TRUE
'read.AsspDataObj' <- 'getAsspDataObj' <- function(fname, begin=0, end=0, samples=FALSE) {
  fname <- prepareFiles(fname)
  # type cast begin/end if integer
  if(inherits(begin, "integer")){
    begin = as.numeric(begin)
  } 
  if(inherits(end, "integer")){
    end = as.numeric(end)
  }
  .External("getDObj2", fname, begin=begin, end=end, samples=samples, PACKAGE="superassp")
}

#' @describeIn AsspDataObj Print a summary of the AsspDataObj (also aliased as summary.AsspDataObj).
#' @param x AsspDataObj
#' @param ... additional arguments (ignored)
#' @method print AsspDataObj
#' @export
"print.AsspDataObj" <- summary.AsspDataObj <- function(x, ...)
{
    temp <- attr(x, "filePath")
    if (is.null(temp)) {
        cat("In-memory Assp Data Object\n")
    }
    else {
        cat(paste("Assp Data Object of file ", temp, ".\n", sep=""))
    }
    cat(sprintf("Format: %s (%s)\n", AsspFileFormat(x), AsspDataFormat(x)))
    cat(paste(as.integer(n_records.AsspDataObj(x)),
              "records at", attr(x, 'sampleRate'), "Hz\n"))
    cat(sprintf("Duration: %f s\n", signal_duration.AsspDataObj(x)))
    cat(paste("Number of tracks:", length(names(x)), "\n"))
    for (track in names(x)) {
        cat('\t', track)
        cat(paste(" (", ncol(x[[track]]), " fields)\n", sep=''))
    }
    genVars <- attr(x, 'genericVars')
    if (!is.null(genVars)) {
        cat("\nGeneric variables:\n")
        for (var in names(genVars)) {
            cat(sprintf("  %s:", var))
            if (genVars[[var]]$Type %in% c("CHAR", "BYTE")) {
                cat(sprintf("\t%s\n", genVars[[var]]$Value))
            } else {
                cat(sprintf("\t%f\n", genVars[[var]]$Value))
            }
            cat(sprintf("    (%s)\n", genVars[[var]]$Type))
        }
    }
}


#' Writes an object of class AsspDataObj to a file given the meta information
#' contained in the object.
#'
#' @title write.AsspDataObj to file
#' @param dobj an object of class AsspDataObj
#' @param file file name as a character string, defaults to the
#' \code{filePath} attribute of the AsspDataObj
#' @return NULL
#' @author Lasse Bombien
#' @useDynLib superassp, .registration = TRUE
"write.AsspDataObj" <- function (dobj, file=attr(dobj, 'filePath'))
  {
    if (is.null(file))
      stop('File path not set internally. Please specify!')
    file <- path.expand(file)
    .Call("writeDObj_", dobj, file, PACKAGE="superassp")
  }

#' Internal type predicate for AsspDataObj.
#'
#' Internal use only. External callers should use `inherits(x, "AsspDataObj")`.
#' @param x Object to test.
#' @param ... Ignored.
#' @return Logical scalar.
#' @keywords internal
#' @noRd
is.AsspDataObj <- function (x, ...)
  {
    if (!inherits(x, "AsspDataObj"))
      return (FALSE)
    return (TRUE)
  }


#' Remove a track from an AsspDataObj.
#'
#' Internal class helper. Returns the object without the track named
#' `trackname`.
#'
#' @param dobj An object of class AsspDataObj
#' @param trackname the name of a track in this object
#' @return The object without the track named trackname
#' @author Lasse Bombien
#' @keywords internal
#' @noRd
delTrack <- function (dobj, trackname)
  {
    if (!is.AsspDataObj (dobj))
      stop ('First argument must be a AsspDataObj.')

    w <- which (names (dobj) == trackname)
    if (length (w) != 1)
      stop ('Invalid trackname')

    ## remove track
    dobj[[trackname]] <- NULL
    ## remove
    attr(dobj, 'trackFormats') <- attr(dobj, 'trackFormats')[-w]
    
    return (dobj)
  }

#' Add a track to an AsspDataObj.
#'
#' Internal class helper. Extends `dobj` with a new track named `trackname`.
#' If a track with that name already exists and `deleteExisting = FALSE` the
#' function errors. If TRUE the existing track is removed (see [delTrack()]).
#' `data` is coerced to a numeric matrix and must have the same number of
#' rows as existing tracks.
#'
#' @param dobj The data object to which the data is to be added
#' @param trackname The name of the new track
#' @param data a matrix with values
#' @param format format for binary writing to file (defaults to 'INT16')
#' @param deleteExisting Delete existing track with the same (default: FALSE)
#' @return the object including the new track
#' @author Lasse Bombien
#' @keywords internal
#' @noRd
addTrack <- function (dobj, trackname, data, format = 'INT16',
                      deleteExisting=FALSE) {
  if (!is.AsspDataObj(dobj))
    stop('dobj must be an AsspDataObj.')
  
  if (!is.numeric(data))
    stop('data must be a numeric matrix')
  
  if (!is.character(trackname) | length(trackname) != 1)
    stop('trackname must be an atomic string.')
  
  data <- as.matrix(data)
  
  tracks <- names(dobj)
  w <- tracks  == trackname
  if (any(w) & !deleteExisting)
    stop(paste('Track', trackname,
                'exists and will not be deleted',
                '("deleteExisting" argument)'))
  if (length(tracks) == 1 & any(w)) {
      ## this is fine: the only track will be replaced
  } else if (length(tracks) > 0) {
    if (nrow(data) != nrow(dobj[[1]]))
      stop(paste("number of rows in data must match number of rows in",
                  "existing tracks."))
  }

  dobj[[trackname]] <- data
  if (any(w))
    attr(dobj, 'trackFormats')[w] <- format
  else
    append(attr(dobj, 'trackFormats'), format)

  return(dobj)
}

#' @rdname assp_accessors
#' @exportS3Method
tracks.AsspDataObj <- function(x, ...) track_names.AsspDataObj(x, ...)

#' Get and set AsspFileFormat (internal).
#'
#' `libassp` handles a number of file formats common in speech research.
#' This helper exposes the active file format of an AsspDataObj so internal
#' code can convert between them where reasonable. Format specifiers come
#' from `AsspFileFormats` and exist as either code name or code number.
#'
#' @param x an object of class AsspDataObj
#' @return for `AsspFileFormat` the code name of the object's currently
#'   set file format
#' @author Lasse Bombien
#' @keywords internal
#' @noRd
AsspFileFormat <- function(x) {
  ## file format is in the first element (of two) in the fileInfo attribute
  xx <- x
  if (!is.AsspDataObj(xx))
    stop('Argument must be an object of class AsspDataObj')
  curFormat <- attr(xx, 'fileInfo')[1]
  ind <- match(curFormat, AsspFileFormats)
  if (is.na(ind))
    stop('Invalid file format. This AsspDataObj has been messed with!')
  return(names(AsspFileFormats)[ind])
}

#' Setter form for `AsspFileFormat` (internal).
#'
#' @param x an object of class AsspDataObj
#' @param value an integer or a string indicating the new file format
#' @return the updated object
#' @keywords internal
#' @noRd
"AsspFileFormat<-" <- function(x, value) {
  value <- value[1]
  if (!is.AsspDataObj(x))
    stop('Argument must be an object of class AsspDataObj')
  fi  <- attr(x, 'fileInfo')
  if (is.numeric(value)) {
    ind <- match(value, AsspFileFormats)
  } else if (is.character(value)) {
    ind <- match(value, names(AsspFileFormats))
  } else {
    stop ('format must be an integer or a string.')
  }
  if (is.na(ind))
    stop('format does not specify a valid file format.')
  fi[1]  <- AsspFileFormats[ind]
  attr(x, 'fileInfo')  <- as.integer(fi)
  x
}

#' Get/set the data format of an AsspDataObj (internal).
#'
#' `libassp` can store data in binary and ASCII format. This helper reports
#' (and `<-` sets) the data format. Valid values: `'ascii'` (or `1`) and
#' `'binary'` (or `2`).
#'
#' @param x an object of class AsspDataObj
#' @return a string representing the current data format
#' @author Lasse Bombien
#' @keywords internal
#' @noRd
AsspDataFormat <- function(x) {
  f <- attr(x, 'fileInfo')[2]
  if (f==1) 
    return('ascii')
  else if (f==2)
    return('binary')
  else
    stop('Invalid data format. This AsspDataObj has been messed with!')
}

#' Setter form for `AsspDataFormat` (internal).
#'
#' @param x an object of class AsspDataObj
#' @param value an integer or a string indicating the new data format
#' @return the updated object
#' @keywords internal
#' @noRd
"AsspDataFormat<-" <- function(x, value) {
  value <- value[1]
  fi <- attr(x, 'fileInfo')
  if (is.numeric(value)) {
    if (value %in% c(1,2))
      fi[2] <- value
    else
      stop('Invalid data format specified')
  } else if (is.character(value)) {
    formats <- c('ascii', 'binary')
    ind <- charmatch(tolower(value), formats)
    if (is.na(ind))
      stop('Invalid data format specified')
    fi[2] <- ind
  } else 
    stop('New value must be an integer or a string.')
  attr(x, 'fileInfo') <- as.integer(fi)
  x
}

#' @rdname assp_accessors
#' @exportS3Method
dur.AsspDataObj <- function(x, ...) signal_duration.AsspDataObj(x, ...)

#' @rdname assp_accessors
#' @exportS3Method
numRecs.AsspDataObj <- function(x, ...) n_records.AsspDataObj(x, ...)

#' @rdname assp_accessors
#' @exportS3Method
rate.AsspDataObj <- function(x, ...) sample_rate.AsspDataObj(x, ...)

#' @rdname assp_accessors
#' @exportS3Method
startTime.AsspDataObj <- function(x, ...) start_time.AsspDataObj(x, ...)

#' @rdname assp_accessors
#' @exportS3Method
tracks.JsonTrackObj    <- function(x, ...) track_names.JsonTrackObj(x, ...)

#' @rdname assp_accessors
#' @exportS3Method
dur.JsonTrackObj       <- function(x, ...) signal_duration.JsonTrackObj(x, ...)

#' @rdname assp_accessors
#' @exportS3Method
rate.JsonTrackObj      <- function(x, ...) sample_rate.JsonTrackObj(x, ...)

#' @rdname assp_accessors
#' @exportS3Method
numRecs.JsonTrackObj   <- function(x, ...) n_records.JsonTrackObj(x, ...)

#' @rdname assp_accessors
#' @exportS3Method
startTime.JsonTrackObj <- function(x, ...) start_time.JsonTrackObj(x, ...)

#' Helper function to parse unit from column name
#'
#' Extracts unit string from column names ending with "\[unit\]"
#'
#' @param col_name Character; column name to parse
#' @return Character; unit string or NA if no unit found
#' @keywords internal
.parse_unit_from_colname <- function(col_name) {
  # Pattern: match [...] at the end of the string
  pattern <- ".*\\[(.+)\\]$"

  if (grepl(pattern, col_name)) {
    # Extract unit from brackets using backreference
    unit_str <- sub(pattern, "\\1", col_name)
    return(unit_str)
  }

  return(NA_character_)
}

#' Helper function to try converting column to units
#'
#' Attempts to convert a numeric column to units. If successful, returns
#' the units object. If it fails, returns the original column and issues a warning.
#'
#' @param col Numeric vector; column data
#' @param unit_str Character; unit string to convert to
#' @param col_name Character; column name for warning messages
#' @return Numeric vector or units object
#' @keywords internal
.try_convert_to_units <- function(col, unit_str, col_name) {
  if (is.na(unit_str)) {
    return(col)
  }

  tryCatch({
    # Try to convert to units
    unit_obj <- units::as_units(unit_str)
    result <- col * unit_obj
    return(result)
  }, error = function(e) {
    warning("Could not convert column '", col_name, "' to unit '", unit_str,
            "': ", e$message, call. = FALSE)
    return(col)
  })
}

#' Old as.data.frame.AsspDataObj - REPLACED by R/assp_dataobj.R version
#' This version is kept commented out to avoid conflicts.
#' The new version in R/assp_dataobj.R supports template expansion.
# "as.data.frame.AsspDataObj" <- function(x, name.separator = "",
#                                         convert_units = TRUE, ...){
#   frame_time = seq(from = start_time.AsspDataObj(x),
#                    by = 1/sample_rate.AsspDataObj(x),
#                    length.out = n_records.AsspDataObj(x)) * 1000
#
#   all_tracks = do.call(cbind, x)
#
#   makeColumnNames <- function(a,n){
#     if(!is.null(ncol(a)) && ncol(a) > 1 ){
#       outname <- paste(n,seq(from=1,to=ncol(a),by=1),sep=name.separator)
#     }else{
#       outname <- n
#     }
#     return(outname)
#   }
#
#   colnames(all_tracks) = purrr::list_c(purrr::imap(x,makeColumnNames))
#   result_df <- as.data.frame(cbind(frame_time, all_tracks))
#
#   # Convert columns to units if requested
#   if (convert_units) {
#     for (col_name in colnames(result_df)) {
#       if (col_name == "frame_time") next  # Skip time column
#
#       unit_str <- .parse_unit_from_colname(col_name)
#       if (!is.na(unit_str)) {
#         result_df[[col_name]] <- .try_convert_to_units(
#           result_df[[col_name]],
#           unit_str,
#           col_name
#         )
#       }
#     }
#   }
#
#   return(result_df)
# }
