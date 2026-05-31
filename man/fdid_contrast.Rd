\name{fdid_contrast}
\alias{fdid_contrast}
\title{Report a continuous-G curve contrast}
\usage{
fdid_contrast(
  x,
  g0,
  g1,
  type = c("level", "derivative"),
  inference = c("auto", "replicate", "vcov", "band", "pointwise"),
  interpolate = c("linear", "none"),
  alpha = NULL
)
}
\arguments{
\item{x}{An object returned by \code{\link{fdid}} with \code{method} equal
to \code{"kernel"} or \code{"dml_flex"}.}

\item{g0, g1}{Numeric scalar endpoints for the contrast.}

\item{type}{Curve type. \code{"level"} is the usual Notes-style contrast;
\code{"derivative"} reports a contrast between derivative values.}

\item{inference}{Inference rule. \code{"auto"} uses stored curve covariance
matrices when available, then replicate curves, then conservative
band-implied intervals, then approximate pointwise intervals.
\code{"pointwise"} ignores endpoint covariance and is labeled approximate.
\code{"band"} uses stored simultaneous or bootstrap bands when available.
\code{"replicate"} uses stored bootstrap or multiplier replicate curves.
\code{"vcov"} uses the covariance identity \eqn{V_{11} + V_{22} - 2V_{12}}
when a curve covariance matrix is stored.}

\item{interpolate}{If \code{"linear"}, linearly interpolate between stored
\code{eval_g} rows. If \code{"none"}, \code{g0} and \code{g1} must be on
the fitted grid.}

\item{alpha}{Significance level. Defaults to the fitted object's
\code{alpha}, or 0.05 for older objects.}
}
\value{
A one-row data frame with estimate, standard error when available,
confidence limits, inference label, and diagnostic note.
}
\description{
Extract a level-curve contrast, such as \eqn{\theta(g_1)-\theta(g_0)}
for \code{method = "dml_flex"} or \eqn{\mu(g_1)-\mu(g_0)} for
\code{method = "kernel"}. The helper uses stored curve values and, when
available, stored curve covariance matrices or replicate curves.
}
