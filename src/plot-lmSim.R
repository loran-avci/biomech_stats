plot.lmSim <- function(obj, which=c(1L:3L), rob=FALSE, SEED=NULL, Nsim=19)
{
    ## rkst, Version 14/4/2016
    X <- model.matrix(obj)
    x.n <- nrow(X)
    yh <- fitted(obj)
    sigma <- if(rob) mad(resid(obj)) else summary(obj)$sigma
    show <- rep(FALSE, 3)
    show[which] <- TRUE
    w <- obj$weights
    if(is.null(w)) w <- rep(1, nrow(X))
    sqrtW <- sqrt(w)


    ## Tukey-Anscombe Plot
    if(show[1L]){
        res <- resid(obj)
        ylim <- extendrange(r=range(res, na.rm = TRUE), f = 0.08)
        plot(yh*sqrtW, res*sqrtW, xlab = "Fitted values",font.main = 1,
             ylab = "Residuals", main="Tukey-Anscombe", ylim = ylim, type = "n")
        if(!is.null(SEED)) set.seed(SEED)
        for(i in 1:Nsim){
            FIT <- lm.wfit(x=X, y=yh + rnorm(x.n,0,sigma/sqrtW), w=w)
            lines(lowess(fitted(FIT)*sqrtW, resid(FIT)*sqrtW,
                         f=2/3, iter=3),col="grey")#, alpha = 0.5)
        }
        abline(h = 0, lty = 3) ## , col = "gray"
        lines(lowess(yh*sqrtW, res*sqrtW, f=2/3, iter=3),
              lwd=1.5, col="red")
    }
    if(any(show[2L:3L])){
        ## Be careful, there are problems when hii==1
        hii <- lm.influence(obj, do.coef=FALSE)$hat
        if (any(isInf <- hii >= 1)) {
            warning(gettextf(paste("not plotting observations",
                                   "with leverage one:\n  %s"),
                             paste(which(isInf), collapse = ", ")),
                    call. = FALSE, domain = NA)
        }
        rdf <- obj$df.residual
        f.stdres <- function(wls) {
            r <- resid(wls)
            sr <- sqrt(w)*r/(sqrt((1 - hii) * sum(w*r^2)/rdf))
            sr[hii >= 1] <- NaN
            sr
        }
    }
    ## normal plot
    if(show[2L]){
        SIM <- matrix(NaN, ncol=Nsim, nrow=x.n)
        if(Nsim>0){
            if(!is.null(SEED)) set.seed(SEED)
            for(i in 1: Nsim){
                FIT <- lm.wfit(x=X, y=yh + rnorm(x.n,0,sigma/sqrtW), w=w)
                SIM[!isInf,i] <- sort(qqnorm(f.stdres(FIT), plot.it=FALSE)$y)
            }
        }
        RQQN <- qqnorm(f.stdres(obj), plot.it=FALSE)
        ylim <- range(c(RQQN$y, SIM), finite=TRUE)
        plot(range(RQQN$x, finite=TRUE), ylim, type="n",
             xlab = "Theoretical Quantiles", ylab = "Standardized Residuals", font.main = 1,main = "Normal Q-Q")
        if(Nsim>0)
          points(rep(sort(RQQN$x),Nsim), as.vector(SIM[!isInf,]), col="gray")#, alpha =0.5)
          points(RQQN$x, RQQN$y, lwd=2)
    }
    ## Scale-Location Plot
    if(show[3L]){
        sqrtabsR <- sqrt(abs(f.stdres(obj)))
        ylim <- c(0, max(sqrtabsR[is.finite(sqrtabsR)]))
        yl <- as.expression(substitute(sqrt(abs(YL)),
            list(YL=as.name("Standardized Residuals"))))
        plot(yh, sqrtabsR, xlab="Fitted values", ylab=yl, main="Scale-Location",font.main = 1,
             ylim=ylim, type="n")
        if(!is.null(SEED)) set.seed(SEED)
        for(i in 1:Nsim){
            FIT <- lm.wfit(x=X, y=yh + rnorm(x.n,0,sigma/sqrtW), w=w)
            h.sRes  <- sqrt(abs(f.stdres(FIT)))
            ok <- is.finite(h.sRes)
            lines(lowess(fitted(FIT)[ok], h.sRes[ok], f=2/3, iter=3),
                  col="grey")#, alpha = 0.5)
        }
        ok <- is.finite(sqrtabsR)
        lines(lowess(yh[ok], sqrtabsR[ok], f=2/3, iter=3), lwd=1.5, col="red")
    }
    invisible()
}
