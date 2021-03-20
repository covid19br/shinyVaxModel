library(lpSolve)

#' opt_vax_rate.
#' 
#' Optimal vaccination rate function, 
#' its optimaze the problem of rollout vaccine solving a by non-linear programming solver.
#'
#' @param VAX.INITIAL.STORAGE.NUM # Number of READY vaccines on day one 'VAX.INITIAL.STORAGE.NUM'.
#' @param VAX.PRODUCTION.RATE # Number of vaccines produced per day 'VAX.PRODUCTION.RATE'.
#' @param MAX.VAC.RATE # Max number of vaccine applications per day 'MAX.VAC.RATE'.
#' @param VAX.WINDOW.DAYS # Time window between first and second vaccines 'VAX.WINDOWS.DAYS'. 
#' Its varies with the vaccine type.
#' @param SECOND.VAX.LOSS.FRAC # Fraction of people that take the first but not the second dose 'SECODNE.VAX.LOSS.FRAC'.
#' @param MAX.TIME.DAYS # time until end of vaccination schedule, in days 'MAX.TIME.DAYS'.
#' @import lpSolve
#'
#' @return
#' @export 
#'
#' @examples
opt_vax_rate <- function(VAX.INITIAL.STORAGE.NUM, VAX.PRODUCTION.RATE,
                         MAX.VAC.RATE, VAX.WINDOW.DAYS, SECOND.VAX.LOSS.FRAC,
                         MAX.TIME.DAYS) {
    times <- seq(0, MAX.TIME.DAYS)

    alpha = 1 - SECOND.VAX.LOSS.FRAC
    V.T = max(VAX.WINDOW.DAYS * (alpha * MAX.VAC.RATE - VAX.PRODUCTION.RATE), 0)
   
    # boundary-free case
    if (VAX.INITIAL.STORAGE.NUM + MAX.TIME.DAYS * (VAX.PRODUCTION.RATE - MAX.VAC.RATE) -
        (MAX.TIME.DAYS - VAX.WINDOW.DAYS) * alpha * MAX.VAC.RATE >= V.T) {
        return(rep(MAX.VAC.RATE, length(times)))
    }

    # feasibility is trivial: V(T) > 0 if VAX.RATE == 0

    # solution using LP
    N <- length(times)
    M <- matrix(0, N, N)
    M[row(M)-col(M) >= 0] = 1
    M.2 <- matrix(0, N, N)
    M.2[seq(1 + VAX.WINDOW.DAYS, N),] <- M[seq(1, N - VAX.WINDOW.DAYS),]

    M.total <- - M - alpha*M.2

    # objective
    C.T <- colSums(M.total)
    ## inequality constraints (Ax <= b)
    A <- M.total
    # diagonal matrix offset by window
    offdiag <- matrix(0, nrow = N, ncol = N)
    offdiag[row(offdiag) - col(offdiag) == VAX.WINDOW.DAYS] = 1
    A <- rbind(- M.total,
               diag(nrow = N, ncol = N) + alpha * offdiag,
               M[N,]
               )
    b <- c(VAX.INITIAL.STORAGE.NUM + VAX.PRODUCTION.RATE * times,
           rep(MAX.VAC.RATE, N),
           (VAX.INITIAL.STORAGE.NUM + VAX.PRODUCTION.RATE * (times[N] + VAX.WINDOW.DAYS))/(1+alpha)
           )
    xopt <-  lp(direction="min",
                objective.in = C.T,
                const.mat = A,
                const.dir = rep("<=", length(b)),
                const.rhs = b)
    #print(paste("Vacinas restantes no tempo final:",
    #            sum(M.total[nrow(M.total),] * xopt$solution) + b[length(times)]))
    return(xopt$solution)
}

#' plot_vac_schedule.
#' 
#' A function to plot the rollout vaccine after the optimal solution.
#' 
#' @param OPT.VAX.RATE # Optimal vaccine rate from the `opt_vax_rate` 'OPT.VAX.RATE'.
#' @param VAX.INITIAL.STORAGE.NUM # Number of READY vaccines on day one 'VAX.INITIAL.STORAGE.NUM'.
#' @param VAX.PRODUCTION.RATE # Number of vaccines produced per day 'VAX.PRODUCTION.RATE'.
#' @param MAX.VAC.RATE # Max number of vaccine applications per day 'MAX.VAC.RATE'.
#' @param VAX.WINDOW.DAYS # Time window between first and second vaccines 'VAX.WINDOWS.DAYS'. 
#' Its varies with the vaccine type.
#' @param SECOND.VAX.LOSS.FRAC # Fraction of people that take the first but not the second dose 'SECODNE.VAX.LOSS.FRAC'.
#' @param MAX.TIME.DAYS # time until end of vaccination schedule, in days 'MAX.TIME.DAYS'.
#' @import lpsolve
#'
#' @return
#' @export
#'
#' @examples
plot_vac_schedule <- function(OPT.VAX.RATE, VAX.INITIAL.STORAGE.NUM,
                              VAX.PRODUCTION.RATE, MAX.VAC.RATE,
                              VAX.WINDOW.DAYS, SECOND.VAX.LOSS.FRAC,
                              MAX.TIME.DAYS) {
    times <- seq(0, MAX.TIME.DAYS)
    alpha = 1 - SECOND.VAX.LOSS.FRAC
    V <- VAX.INITIAL.STORAGE.NUM + VAX.PRODUCTION.RATE * times -
        cumsum(OPT.VAX.RATE) -
        alpha * c(rep(0, VAX.WINDOW.DAYS),
                  cumsum(OPT.VAX.RATE[seq(1, length(times) - VAX.WINDOW.DAYS)]))

    par(mar = c(5, 4, 4, 4) + 0.3)
    plot(times, OPT.VAX.RATE,
         xlab="time (days)",
         ylab="vaccination rate (doses/day)",
         type='l')
    lines(times, alpha * c(rep(0, VAX.WINDOW.DAYS),
                           OPT.VAX.RATE[-seq(length(times)+1-VAX.WINDOW.DAYS, length(times))]),
          type='l', lty=2)
    par(new = TRUE)
    plot(times, V, type = "l", axes = FALSE, bty = "n", lty = 3, xlab = "", ylab = "", col="blue")
    axis(side=4, at = pretty(round(range(V), -5)))
    mtext("Vaccine doses stored", side=4, line=3, col="blue")
}

#' interpolate.VAX.RATE.
#' 
#' A function to interpolate vaccine rate rollout give the vaccine rate
#'
#' @param OPT.VAX.RATE # Optimal vaccine rate from the `opt_vax_rate` 'OPT.VAX.RATE'.
#'
#' @return
#' @export
#'
#' @examples
interpolate.VAX.RATE <- function(OPT.VAX.RATE) {
    VAX.RATE <- function(t){
        if (t <= 0)
            return(0)
        tu <- ceiling(t)
        x <- tu - t
        s = 0.1
        if (x > 1-s && tu > 1)
            return((x)*OPT.VAX.RATE[tu] +
                   (1-x)*(OPT.VAX.RATE[tu] + OPT.VAX.RATE[tu-1])/2)
        if (x < s)
            return((1-x/s)*OPT.VAX.RATE[tu] +
                   (x/s)*(OPT.VAX.RATE[tu] + OPT.VAX.RATE[tu+1])/2)
        return(OPT.VAX.RATE[tu])
    }
}