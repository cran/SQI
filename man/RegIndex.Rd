% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SQI.R
\name{RegIndex}
\alias{RegIndex}
\title{Soil Quality Index Based on Regression}
\usage{
RegIndex(DataFrame, Dep_col, OptimumValue)
}
\arguments{
\item{DataFrame}{Data set with first column as factors}

\item{Dep_col}{Dependent variable column number}

\item{OptimumValue}{Optimum value of each variable; Minimum and maximum coded as "1111" and "9999" respectively.}
}
\value{
\itemize{
  \item RegIndex: Final index
}
}
\description{
Soil Quality Index Based on Regression
}
\examples{
library("SQI")
OP<-c(7,1111,9999,9999,9999,9999,9999,9999,9999,9999,1111)
RIndex<-RegIndex(DataFrame = Data,Dep_col=7,OptimumValue = OP)
}
\references{
\itemize{
\item Bastida, F., Zsolnay, A., Hernández, T., & García, C. (2008). Past, present and future of soil quality indices: a biological perspective. Geoderma, 147(3-4), 159-171.
\item Doran, J. W., & Parkin, T. B. (1994). Defining and assessing soil quality. Defining soil quality for a sustainable environment, 35, 1-21.
\item Mukherjee, A., & Lal, R. (2014). Comparison of soil quality index using three methods. PloS one, 9(8), e105981.
}
}
