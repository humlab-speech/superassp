% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/praat_ssff.R
\name{get_praat}
\alias{get_praat}
\title{Utility function for getting the full path of the Praat executable.}
\usage{
get_praat(praat_path = NULL)
}
\arguments{
\item{praat_path}{A character string containing the path to the executable that the function was able to find (or the executable that the function was able to verify the existance of), or NULL if no Praat executable was found and verified.}
}
\value{
A character string containing the path to the executable that the function was able to find (or the executable that the function was able to verify the existance of), or NULL if no Praat executable was found and verified.
}
\description{
This function checks the system variables and deduces where Praat is installed. On the OSX platform (Darwin) the Praat app is assumed to exist in the Applications folder,
and the actual binary inside of the application package is then used. If not OSX, then the function will search the default search paths for executables set up in the OS.
If an explicit path is given, then function will just check whether the executable is actualy present there.
}
