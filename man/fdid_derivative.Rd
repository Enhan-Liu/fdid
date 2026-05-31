\name{fdid_derivative}
\alias{fdid_derivative}
\title{Report a continuous-G derivative value}
\usage{
fdid_derivative(
  x,
  g0,
  inference = c("auto", "replicate", "vcov", "band", "pointwise"),
  interpolate = c("linear", "none"),
  alpha = NULL
)
}
\arguments{
\item{x}{An object returned by \code{\link{fdid}} with \code{method} equal
to \code{"kernel"} or \code{"dml_flex"}.}

\item{g0}{Numeric scalar evaluation point.}

\item{inference}{Inference rule. \code{"auto"} uses stored curve covariance
matrices when available, then replicate curves, then stored band intervals,
then pointwise intervals. \code{"vcov"} uses a stored curve covariance matrix
at the requested value when available.
\code{"band"} uses stored simultaneous or bootstrap bands when available.
\code{"replicate"} uses stored bootstrap or multiplier replicate curves.}

\item{interpolate}{If \code{"linear"}, linearly interpolate between stored
\code{eval_g} rows. If \code{"none"}, \code{g0} must be on the fitted grid.}

\item{alpha}{Significance level. Defaults to the fitted object's
\code{alpha}, or 0.05 for older objects.}
}
\value{
A one-row data frame with estimate, standard error when available,
confidence limits, inference label, and diagnostic note.
}
\description{
Extract \eqn{\theta'(g_0)} for \code{method = "dml_flex"} or
\eqn{\delta(g_0)} for \code{method = "kernel"} from a fitted continuous-G
curve object.
}
