library (devtools)
install_github("maciejdanko/PCLM")
library (PCLM)

 # Use a simple data set. Naive life-table.
 # Age: 
 x <- c(0, 1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
 # Death counts:
 dx <- c(38, 37, 17, 104, 181, 209, 452, 1190, 2436, 3164, 1852, 307, 13)
 # Survivors at the beginning of age class
 lx <- sum(dx)-c(0, cumsum(dx[-length(dx)]))
 # Interval length
 n <- diff(c(x,110))
 # Approximation of mortality per age class
 mx <- - log(1 - dx / lx) / n
 # Mid-interval vector
 xh <- x + n / 2
 # Approximated exposures
 Lx <- n * (lx - dx) + 0.5 * dx *n
 
 # *** Use PCLM
 # Ungroup dataset with out.step equal minimal interval length
 min(diff(x))
 AU10p.1a <- pclm.default(x, dx)
 print(AU10p.1a)
 plot(AU10p.1a)

 # Ungroup AU10 with out.step equal minimal interval length
 # and get good estimates of nax
 AU10p.1b <- pclm.default(x, dx, control = list(x.div = 10))
 print(AU10p.1b)
 plot(AU10p.1b)
 # This time number of internal (raw) PCLM classes was high
 # and automatically P-splines were used to prevent long computations

 # This number can be estimated before performing
 # PCLM calclualtions:
 pclm.nclasses(x, control = list(x.div = 10))
 # which is the same as in the fitted model
 length(AU10p.1b$raw$x)
 # whereas number of classes in the aggregated life-table
 # depends on out.step
 length(AU10p.1b$grouped$x) 

 # To speed-up computations we can decrease the number of P-spline knots
 AU10p.1c <- pclm.default(x, dx, control = list(x.div = 10,
                      bs.use = TRUE, bs.df.max = 100))

 # *** Diagnostic plots for fitted PCLM model
 # Aggregated PCLM fit:
 plot(AU10p.1b, type = 'aggregated')
 # Raw PCLM fit before aggregation:
 plot(AU10p.1b, type = 'nonaggregated')

 # In this PCLM fit aggregated life-table is identical
 # with nonaggregated
 plot(AU10p.1a, type = 'aggregated')
 plot(AU10p.1a, type = 'nonaggregated')

 # *** Combined summary of pash and pclm objects
 summary(AU10p.1a)
 summary(AU10p.1b)
 summary(AU10p.1c)

 # *** Smooth and aggregate data into 12-year interval
 AU10p.2 <- pclm.default(x, dx, out.step = 12)
 print(AU10p.2)
 print(AU10p.2, type = 'aggregated') # grouped PCLM life-table
 print(AU10p.2, type = 'nonaggregated') # raw PCLM life-table
 plot(AU10p.2)

 # *******************************************************************
 # Usage of PCLM methods to fit and plot mortality data
 # *******************************************************************

 AU10p.4a <- pclm.default(x, dx, control = list(x.div = 5))
 X <- AU10p.4a$grouped$x
 M <- -log(1 - AU10p.4a$grouped$dx/AU10p.4a$grouped$lx)
 plot(X, log10(M), type='l', lwd = 2,
      xlim=c(0,130), xlab='Age', ylab='log_10 mortality', col = 2)
 lines(xh, log10(mx1), type = 'p')
 tail(AU10p.4a, n = 10)
 #note that lx has standardized values

 # Improving the plot to cover more age classes
 AU10p.4b <- pclm.default(x, dx, control = list(zero.class.end = 150,
                      x.div = 4))
 X <- AU10p.4b$grouped$x
 M <- -log(1 - AU10p.4b$grouped$dx / AU10p.4b$grouped$lx)
  
 plot(X, log10(M), type='l', lwd = 2,
      xlim=c(0,130), xlab='Age', ylab='log_10 mortality', col = 2)
 lines(xh, log10(mx1), type = 'p')
 tail(AU10p.4a, n = 10)

 # The change of the order of the difference in pclm algorithm may
 # affect hte interpretation of the tail.
 # Try to check pclm.deg = 4 and 5.
 AU10p.4c <- pclm.default(x, dx, control = list(zero.class.end = 150,
                      x.div = 1, pclm.deg = 4))
 X <- AU10p.4c$grouped$x
 M <- -log(1 - AU10p.4c$grouped$dx / AU10p.4c$grouped$lx)
 plot(X, log10(M), type='l', lwd = 2,
      xlim=c(0,130), xlab='Age', ylab='log_10 mortality', col = 2)
 lines(xh, log10(mx1), type = 'p')

 # Using exposures to fit mortality, 
 # Notice that different approximation of mortality rate is used than in
 # previous cases.
 AU10p.4c <- pclm.default(x, dx, exposures = Lx, control = list(zero.class.end = 150,
                      x.div = 1, pclm.deg = 2, bs.use = FALSE))
 X <- AU10p.4c$grouped$x
 M <- AU10p.4c$grouped$mx
 plot(X, log10(M), type='l', lwd = 2,
      xlim=c(0,130), xlab='Age', ylab='log_10 mortality', col = 2)
 lines(xh, log10((dx / Lx) / n), type = 'p')

 # *******************************************************************
 # Usage of PCLM methods for more complicated dataset
 # - understanding the out.step, x.div, and interval multiple
 # *******************************************************************

 # *** Generate a dataset with varying and fractional interval lengths
 x <- c(0, 0.6, 1, 1.4, 3, 5.2, 6.4, 8.6, 11, 15,
        17.2, 19, 20.8, 23, 25, 30)
 dx <- ceiling(10000*diff(pgamma(x, shape = 3.8, rate = .4)))
 barplot(dx/diff(x), width = c(diff(x), 2)) # preview
 lx <- 10000-c(0, cumsum(dx))
 dx <- c(dx, lx[length(lx)])  

 # *** Fit PCLM with automatic out.step
 Bp1 <- pclm.default(x, dx)
 # Output interval length (out.step) is automatically set to 0.4
 # which is the minimal interval length in original data.
 min(diff(x))
 summary(Bp1) #new out.step can be also read from summary
 plot(Bp1)

 # *** Setting manually out.step
 Bp2 <- pclm.default(x, dx, out.step = 1)
 plot(Bp2, type = 'aggregated') # The fit with out.step = 1
 plot(Bp2, type = 'nonaggregated') # It is clear that
 # PCLM extended internal interval length even without changing x.div
 # It was done because of the fractional parts in x vector.
 # This is also a case for Bp1
 summary(Bp2) #PCLM interval length = 0.2
 Bp2$raw$n[1:10]

 # *** Setting manually out.step to a smaller value than
 #     the smallest original interval length
 Bp3 <- pclm.default(x, dx, out.step = 0.1)
 summary(Bp3)
 # We got a warning as out.step cannot be smaller than
 # smallest age class if x.div = 1

 # We can change x.div to make it possible
 Bp3 <- pclm.default(x, dx, out.step = 0.1, control = list(x.div = 2))
 #0.1 is two times smaller than minimal interval length
 summary(Bp3) # We were able to change the interval
 plot(Bp3)
 # NOTE: In this case x.div has not sufficient value to
 #       get good axn estimates
 Bp3$grouped$ax[1:10]

 # This can be changed by the further increase of x.div
 Bp4 <- pclm.default(x, dx, out.step = 0.1, control = list(x.div = 20))
 Bp4$grouped$ax[1:10]
 # NOTE: This time P-spline approximation was used because
 # the composition matrix was huge

 # Finally, we were able to get our assumed out.step
 Bp4$grouped$n[1:10]

 # In the fitted model the interval multiple (m) is 5.
 (m <- pclm.interval.multiple(x, control = list(x.div = 20)))
 summary(Bp4)
 # Interval multiple determines
 # the maximal interval length in raw PCLM life-table,
 (K <- 1 / m)
 # which is further divided by x.div.
 K / 20
 # Simply: 1 / (m * x.div) = 1 / (5 * 20) = 0.01
 # The interval in the raw PCLM life-table is 10 times shorter than
 # in the grouped life-table
 # interval length in aggregated PCLM life-table:
 Bp4$grouped$n[1:10]/ # divided by
 # interval length in nonaggregated PCLM life-table:
 Bp4$raw$n[1:10]
 # NOTE: The interval for the raw PCLM life-table depends
 # on original interval, m, and x.div,
 # whereas the grouped PCLM interval length is set by out.step,
 # which could be eventually increased if out.step < raw PCLM
 # interval length.

 # **** See more examples in the help for pclm.nclasses() function.