#' JSON Track I/O Functions
#' 
#' Read and write JSON Track Format (JSTF) files using RcppSimdJson
#' for reading and jsonlite for writing.
#'
#' @name json_track_io
NULL

#' Write JSON Track Object to File
#'
#' Writes a JsonTrackObj to a JSON file using jsonlite. Designed for
#' reliable writing with proper formatting and type handling.
#'
#' @param obj JsonTrackObj to write
#' @param file Output file path
#' @param pretty Logical, pretty-print JSON (default: FALSE for efficiency)
#' @param digits Number of decimal digits for numbers (default: 6)
#' @param auto_unbox Logical, automatically unbox single-element arrays (default: TRUE)
#'
#' @return Invisibly returns file path
#' @export
#' @examples
#' \dontrun{
#' # Create a JsonTrackObj
#' obj <- create_json_track_obj(
#'   results = list(f0_mean = 150, f0_sd = 20),
#'   function_name = "lst_example",
#'   file_path = "audio.wav",
#'   sample_rate = 16000,
#'   audio_duration = 5.0
#' )
#' 
#' # Write to file
#' write_json_track(obj, "output.json")
#' }
write_json_track <- function(obj,
                             file,
                             pretty = FALSE,
                             digits = 6,
                             auto_unbox = TRUE) {
  
  # Validate input
  stopifnot(inherits(obj, "JsonTrackObj"))
  validate_json_track(obj)
  
  # Convert to JSON using jsonlite
  json_str <- jsonlite::toJSON(
    obj,
    pretty = pretty,
    digits = digits,
    auto_unbox = auto_unbox,
    force = TRUE,
    null = "null",
    na = "null"
  )
  
  # Write to file
  writeLines(json_str, file)
  
  if (file.exists(file)) {
    message("JSON track written to: ", file)
  } else {
    stop("Failed to write file: ", file)
  }
  
  invisible(file)
}

#' Read JSON Track Object from File
#'
#' Reads a JSTF file using RcppSimdJson for high-performance parsing.
#' Falls back to jsonlite if RcppSimdJson fails.
#'
#' @param file Path to JSON track file
#' @param validate Logical, validate after reading (default: TRUE)
#'
#' @return JsonTrackObj
#' @export
#' @examples
#' \dontrun{
#' # Read JSON track
#' obj <- read_json_track("output.json")
#' 
#' # Convert to data.frame
#' df <- as.data.frame(obj)
#' }
read_json_track <- function(file, validate = TRUE) {
  
  if (!file.exists(file)) {
    stop("File not found: ", file)
  }
  
  # Try RcppSimdJson first (faster)
  obj <- tryCatch({
    read_json_track_simdjson(file)
  }, error = function(e) {
    warning("RcppSimdJson failed, falling back to jsonlite: ", e$message)
    read_json_track_jsonlite(file)
  })
  
  # Validate if requested
  if (validate) {
    validate_json_track(obj)
  }
  
  return(obj)
}

#' Read JSON using RcppSimdJson (fast)
#'
#' @param file Path to JSON file
#' @return JsonTrackObj
#' @keywords internal
read_json_track_simdjson <- function(file) {

  # Check if RcppSimdJson is available
  if (!requireNamespace("RcppSimdJson", quietly = TRUE)) {
    stop("RcppSimdJson package required but not available. Install with: install.packages('RcppSimdJson')")
  }

  # Read JSON using simdjson (no simplification - keep as nested lists)
  parsed <- RcppSimdJson::fload(file, max_simplify_lvl = "list")

  # Convert to JsonTrackObj structure
  obj <- structure(
    parsed,
    class = c("JsonTrackObj", "list")
  )

  return(obj)
}

#' Read JSON using jsonlite (fallback)
#'
#' @param file Path to JSON file
#' @return JsonTrackObj
#' @keywords internal
read_json_track_jsonlite <- function(file) {
  
  # Read JSON using jsonlite
  parsed <- jsonlite::fromJSON(
    file,
    simplifyVector = FALSE,
    simplifyDataFrame = FALSE,
    simplifyMatrix = FALSE
  )
  
  # Convert to JsonTrackObj structure
  obj <- structure(
    parsed,
    class = c("JsonTrackObj", "list")
  )
  
  return(obj)
}

#' Unified Track Reading Interface
#'
#' Reads either SSFF or JSTF files transparently based on file extension.
#'
#' @param file Path to track file (.f0, .fms, .vat, .vsj, etc.)
#' @param begin Start of region to read (seconds, or samples if
#'   \code{samples = TRUE}). Default 0 = file start. Ignored for JSTF files.
#' @param end End of region to read (seconds, or samples if
#'   \code{samples = TRUE}). Default 0 = file end. Ignored for JSTF files.
#' @param samples Logical. If \code{TRUE}, \code{begin}/\code{end} are in
#'   samples; otherwise in seconds. Ignored for JSTF files.
#' @param validate Logical, validate after reading (default: TRUE, JSTF only)
#'
#' @return AsspDataObj (for SSFF) or JsonTrackObj (for JSTF)
#' @export
#' @examples
#' \dontrun{
#' # Read SSFF pitch track
#' pitch <- read_track("audio.f0")
#'
#' # Read SSFF with time windowing
#' pitch <- read_track("audio.f0", begin = 1.0, end = 2.5)
#'
#' # Read JSTF voice quality track
#' vq <- read_track("audio.vat")
#'
#' # Both can be converted to data.frame
#' df1 <- as.data.frame(pitch)
#' df2 <- as.data.frame(vq)
#' }
read_track <- function(file, begin = 0, end = 0, samples = FALSE,
                       validate = TRUE) {

  if (!file.exists(file)) {
    stop("File not found: ", file)
  }

  # Get file extension
  ext <- tools::file_ext(file)

  # Check if it's a known JSTF extension
  jstf_extensions <- get_jstf_extensions()

  if (tolower(ext) %in% jstf_extensions) {
    # Read as JSON track
    return(read_json_track(file, validate = validate))
  } else {
    # SSFF format - use superassp's own reader
    return(read_ssff(file, begin = begin, end = end, samples = samples))
  }
}

#' Get list of JSTF extensions
#'
#' @return Character vector of JSTF file extensions
#' @keywords internal
get_jstf_extensions <- function() {
  
  # Try to read from registry
  registry_file <- system.file("extdata", "json_extensions.csv", 
                               package = "superassp")
  
  if (file.exists(registry_file)) {
    registry <- utils::read.csv(registry_file, stringsAsFactors = FALSE)
    return(registry$extension)
  }
  
  # Fallback: hardcoded list
  return(c("vat", "vsj", "dyp", "vxt", "gem", "egm", "emb", "cmp", 
           "cvq", "avq", "dsi", "vrp", "vtr", "phn"))
}

#' Get JSTF extension for function
#'
#' @param function_name Name of lst_* function
#' @return File extension (without dot)
#' @examples
#' get_jstf_extension("lst_vat")  # "vat"
#' get_jstf_extension("lst_voice_sauce")  # "vsj"
get_jstf_extension <- function(function_name) {

  # Try to read from registry (installed package)
  registry_file <- system.file("extdata", "json_extensions.csv",
                               package = "superassp")

  # Fallback for development mode (devtools::load_all)
  if (!file.exists(registry_file) || registry_file == "") {
    registry_file <- file.path("inst", "extdata", "json_extensions.csv")
  }

  if (file.exists(registry_file)) {
    registry <- utils::read.csv(registry_file, stringsAsFactors = FALSE, check.names = FALSE)
    match_row <- registry[registry[["function"]] == function_name, ]
    if (nrow(match_row) > 0) {
      return(match_row$extension[1])
    }
  }

  # Fallback: try to infer from function name
  func_name_clean <- gsub("^lst_", "", function_name)
  ext <- substr(func_name_clean, 1, 3)

  warning("Extension not found in registry for '", function_name,
          "', using inferred extension: ", ext)

  return(ext)
}
