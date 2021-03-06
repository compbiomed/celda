% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/celda_G.R
\name{celdaHeatmap.celda_G}
\alias{celdaHeatmap.celda_G}
\title{Heatmap for celda_CG}
\usage{
\method{celdaHeatmap}{celda_G}(counts, celda.mod, nfeatures = 25, ...)
}
\arguments{
\item{counts}{Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.}

\item{celda.mod}{Celda object of class `celda_G`.}

\item{nfeatures}{Integer. Maximum number of features to select for each module. Default 25.}

\item{...}{Additional parameters.}
}
\value{
list A list containing the dendrograms and the heatmap grob
}
\description{
Renders an expression heatmap to visualize `celda_CG()` results. The top `nfeatures` for each module will be included in the heatmap.
}
\examples{
celdaHeatmap(celda.G.sim$counts, celda.G.mod)
}
\seealso{
`celda_G()` for clustering features and `celdaTsne()` for generating 2-dimensional coordinates
}
