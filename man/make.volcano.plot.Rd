% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make.volcano.plot.R
\name{make.volcano.plot}
\alias{make.volcano.plot}
\title{Make an enhanced volcano plot}
\usage{
make.volcano.plot(
  dat.p,
  dat.fc,
  vars,
  title,
  subtitle = NULL,
  xlim = c(-6, 8),
  ylim = c(0, 4),
  path = getwd(),
  width = 7,
  height = 8,
  pCutoff = 0.05,
  FCcutoff = 0.26,
  hline = c(0.01, 0.05),
  hlineCol = c("grey75", "grey25"),
  pointSize = 3,
  labSize = 3,
  col = c("black", "black", "blue", "red3"),
  colAlpha = 1,
  ...
)
}
\description{
Make an enhanced volcano plot
}
