% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adaptconcept2_sFFLHD_R6_desfuncs.R
\name{des_func_relmax}
\alias{des_func_relmax}
\title{Relmax des func}
\usage{
des_func_relmax(mod, XX, return_se = F, N_add = 1000)
}
\arguments{
\item{mod}{mod}

\item{XX}{XX}

\item{return_se}{whether the se prediction should be returned along with
the des, all will be returned in data.frame, this will save
time if calculating the werror function since it is faster
to predict both at once instead of separately}

\item{N_add}{N_add}
}
\value{
des
}
\description{
A des func where output is scaled 0 to 1, max higher
}
