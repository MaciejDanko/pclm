# summary.pclm -------------------------------------------------------------------------

#' Summary of the fitted PCLM object
#'
#' @description
#' \emph{bold{Generic function}}
#' @param object Fitted PCLM object.
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
#' @seealso  \code{\link{plot.pclm}}
#' @keywords internal
#' @export
summary.pclm <- function(object, ...){
  message('Summary of the pclm object:')
  n1 <- diff(object$fit$CompositionMatrix$x)
  n1 <- c(n1, n1[length(n1)])
  z0 <- n1[1]/2
  n2 <- diff(object$fit$X)
  n2 <- c(n2, n2[length(n2)])
  if (!is.na(object$fit$CompositionMatrix$open.int.len)) ind <- -(length(n1):(length(n1) - 1)) else ind <- 1:length(n1)
  cat(paste(paste('PCLM total classes =', length(n2)),
            paste('Number of smoothing parameters for B-/P-splines =', object$fit$control$bs.df),
            paste('Original minimal interval length =', round(min(n1[ind]), 3)),
            paste('Original maximal interval length =', round(max(n1[ind]), 3)),
            paste('Open interval length =', round(object$fit$CompositionMatrix$open.int.len, 3)),
            paste('Fractional age/time class correction (multiple) =', object$m),
            paste('PCLM interval length =', round(min(n2), 3)),
            paste('PCLM class divider (x.div) =', object$x.div),
            paste('PCLM classes per original smallest interval length =', round(min(n1[ind]) / min(n2), 3)),
            paste('PCLM classes per original biggest interval length =', round(max(n1[ind]) / max(n2), 3)), sep='\n'))
  message('\nWarnings list:')
  W <- unlist(object$warn.list)
  print(W,  quote=F)
  cat('\n')
  invisible()
}

# plot.pclm ---------------------------------------------------------------------------

#' Diagnostic plot for PCLM object.
#' @description
#' \emph{bold{Generic function}}
#' @param object Fitted PCLM object.
#' @param type Type of PCLM plot:
#' \itemize{
#' \item{\code{"aggregated"} - Aggregated PCLM fit with interval length of \code{out.step}}.
#' See .
#' \item{\code{"nonaggregated"} - Nonaggregated (raw) PCLM fit with interval
#' of length equal to the shortest original
#' interval length divided by \code{x.div}. See \code{\link{pclm.control}}}.
#' }
#' @param xlab Optional label for the X-axis.
#' @param ylab Optional label for the Y-axis.
#' @param xlim Optional limits of the X-axis.
#' @param ylim Optional limits of the Y-axis.
#' @param legpos.x,legpos.y Position of the \code{\link{legend}}. If \code{legpos.x == NULL} then legend is not plotted.
#'
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
#' @seealso  \code{\link{summary.pclm}}
#' @keywords internal
#' @export
plot.pclm<-function(object, type = c("aggregated", "nonaggregated"), 
                    xlab, ylab, xlim, ylim, legpos.x = "topleft", legpos.y = NULL){
  if (missing(xlab)) xlab <- 'Age or time'
  
  if (length(object$exposures) == 0){
    if (missing(ylab)) ylab <- 'Counts / interval length'
    
    n1 <- diff(object$fit$CompositionMatrix$x)
    n1 <- c(n1, n1[length(n1)])
    n2 <- diff(object$fit$X)
    n2 <- c(n2, n2[length(n2)])
    n3 <- diff(object$grouped$x)
    n3 <- c(n3, n3[length(n3)])
    
    if (missing(xlim)) xlim <- range(c(object$fit$X, object$fit$CompositionMatrix$x))
    if (missing(ylim)) {
      ylim <- range(c(object$fit$Y/n2, object$fit$CompositionMatrix$y/n1, object$grouped$dx/n3))
      ylim[1] <- 0
      ylim[2] <- ylim[2]*1.2
    }
    tmp.lwd <- par('lwd'); par(lwd = 2,xaxs = 'i', yaxs = 'i')
    barplot(width = n1, space = 0, height = object$fit$CompositionMatrix$y / n1, xlab = xlab, ylab = ylab,
            col = 'gold2', border = 'white', xlim = xlim, ylim = ylim)
    par(lwd = tmp.lwd)
    
    lines(object$fit$CompositionMatrix$x, object$fit$CompositionMatrix$y/n1, type = 's')
    axis(1)
    if (tolower(type[1]) == 'nonaggregated'){
      lines(object$fit$X, object$fit$Y/n2, type = 's', col = 'red', lwd = 2)
      AType <- 'PCLM nonaggregated (raw)'
    } else if (tolower(type[1]) == 'aggregated') {
      AType <- 'PCLM aggregated'
      lines(object$grouped$x, object$grouped$dx/n3, type = 's', col = 'red', lwd = 2)
    } else stop('Unknown plotting type.')
    if (length(legpos.x) > 0) legend(x = legpos.x, y = legpos.y, legend = c('Data', AType), bty = 'n', pch = c(15, NA), lty = c(NA, 1), lwd = 2, col = c('gold2', 'red'), pt.cex = 1.8)
    box()
  } else stop('Diagnostic plots for mortality smooth not available yet.')
  invisible()
}

# head.pclm ----------------------------------------------------

#' Head function for PCLM object
#'
#' @description
#' \emph{bold{Generic function}}
#' @param object PCLM object.
#' @param n A single integer. If positive, size for the resulting object: number of rows for a life-table. If negative, all but the n last/first number of elements of x.
#' @param type which life-table  should be returned. One of \code{"aggregated"} or \code{"nonaggregated"}.
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
#' @keywords internal
#' @export
head.pclm<-function(object, n = 6L, type = c("aggregated", "nonaggregated")){
  type <- type[1]
  if (type == "aggregated")
    head(object$grouped, n = n)
  else if (type == "nonaggregated")
    head(object$raw, n = n)
  else stop('Unknown type')
}

# tail.pclm -----------------------------------------------------------

#' Tail function for PCLM object
#'
#' @description
#' \emph{bold{Generic function}}
#' @param object PCLM object.
#' @param n A single integer. If positive, size for the resulting object: number of rows for a life-table. If negative, all but the n last/first number of elements of x.
#' @param type which life-table  should be returned. One of \code{"aggregated"} or \code{"nonaggregated"}.
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
#' @keywords internal
#' @export
tail.pclm<-function(object, n=6L, type = c("aggregated", "nonaggregated")){
  type <- type[1]
  if (type == "aggregated")
    tail(object$grouped, n = n)
  else if (type == "nonaggregated")
    tail(object$raw, n = n)
  else stop('Unknown type')
}

# print.pclm ------------------------------------------------------------

#' Print function for PCLM object
#'
#' @description
#' \emph{bold{Generic function}}
#' @param object PCLM object.
#' @param type which life-table  should be returned. One of \code{"aggregated"} or \code{"nonaggregated"}.
#' @param ... other parameters passed to \code{\link{print}}.
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
#' @keywords internal
#' @export
print.pclm<-function(object, type = c("aggregated", "nonaggregated"), ...){
  type <- type[1]
  if (type == "aggregated")
    print(object$grouped, ...)
  else if (type == "nonaggregated")
    print(object$raw, ...)
  else stop('Unknown type')
}