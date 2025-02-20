# <R package ACSSpack, providing ACSS, Corresponding ACSS, and GLP Algorithm>
#       Copyright (C) <2024>  <Ziqian Yang>


#' Economic data from the GLP paper 
#'
#' A list contains the five data set used in the paper of Giannone, Lenza, and Primiceri (2021).
#' Contains the following data sets: Macro1, Macro2, Micro1, Micro2, and Finance1
#'
#' @format ## `Econ_data`
#' A list contains the five lists, as the 5 data sets
#' \describe{
#'   \item{Macro1}{A list with a vector and a data frame, `Y` and `X`. `Y` is the response vector, with 659 observations in it. `X` is the data frame contains all predictors, with n=659, p=130. 
#'   It have the structure of time series data.}
#'   \item{Macro2}{A list with a vector and a data frame, `Y` and `X`. `Y` is the response vector, with 90 observations in it. `X` is the data frame contains all predictors, with n=90, p=69.
#'   It have the structure of sectional data.}
#'   \item{Micro1}{A list with a vector and a data frame, `Y` and `X`. `Y` is the response vector, with 576 observations in it. `X` is the data frame contains all predictors, with n=576, p=285.
#'   It have the structure of panel data with 48 units on 12 time points.}
#'   \item{Micro2}{A list with a vector and a data frame, `Y` and `X`. `Y` is the response vector, with 312 observations in it. `X` is the data frame contains all predictors, with n=312, p=138.
#'   It have the structure of panel data with 12 units on 26 time points.}
#'   \item{Finance1}{A list with a vector and a data frame, `Y` and `X`. `Y` is the response vector, with 68 observations in it. `X` is the data frame contains all predictors, with n=68, p=16.
#'   It have the structure of time series data.}
#' }
#' @source <https://www.econometricsociety.org/publications/econometrica/2021/09/01/economic-predictions-big-data-illusion-sparsity/supp/17842_Data_and_Programs.zip>
#' @source <https://research.stlouisfed.org/econ/mccracken/fred-databases/>
"Econ_data"