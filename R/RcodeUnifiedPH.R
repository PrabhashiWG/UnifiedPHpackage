#' EM algorithm for arbitrarily censored data subject to left-truncation under the proportional hazards model
#'
#' Fits a parametric proportional hazards model (PH), proposed in Gamage et al. (2022+), to arbitrarily censored data subject to left-truncation via an EM algorithm.
#' @param d1 vector indicating whether an observation is exactly observed (1) or not (0).
#' @param d2 vector indicating whether an observation is interval-censored (1) or not (0).
#' @param d3 vector indicating whether an observation is right-censored (1) or not (0).
#' @param Li the left endpoint of the observed interval, if an observation is left-censored its corresponding entry should be 0.
#' @param Ri the right endpoint of the observed interval, if an observation is right-censored its corresponding entry should be Inf.
#' @param Ei the vector specifying the enrollment times, if no enrollment criteria is used then its corresponding entry should be 0.
#' @param Xp design matrix of predictor variables (in columns), should be specified without an intercept term.
#' @param n.int the number of interior knots to be used.
#' @param order the order of the basis functions.
#' @param g0 initial estimate of the spline coefficients; should be of length n.int+order.
#' @param b0 initial estimate of regression coefficients; should be of length dim(Xp)[2].
#' @param tol the convergence criterion of the EM algorithm, see details for further description.
#' @param t.seq an increasing sequence of points at which the cumulative baseline hazard function is evaluated.
#' @param equal logical, if TRUE knots are spaced evenly across the range of the endpoints of the observed intervals and if FALSE knots are placed at quantiles. Defaults to FALSE.
#' @details The above function fits a parametric proportional hazards model (PH), proposed in Gamage et al. (2022+), to arbitrarily censored data subject to left-truncation via an EM algorithm. For a discussion of starting values, number of interior knots, order, and further details please see Gamage et al. (2022+). The EM algorithm converges when the maximum absolute difference between consecutive parameter updates was less than the specified tol.
#' @return b: estimates of the regression coefficients.
#' @return g: estimates of the spline coefficients.
#' @return ll: the value of the maximized log-likelihood.
#' @return AIC: the Akaike information criterion.
#' @return BIC: the Bayesian information/Schwarz criterion.
#' @return bRi: I-spline basis matrix of dimension c(n.int+order, length(Ri)).
#' @return bLi: I-spline basis matrix of dimension c(n.int+order, length(Li)).
#' @return bt: I-spline basis matrix evaluated at the points t.seq.
#' @return mRi: M-spline basis matrix of dimension c(n.int+order, length(Ri)).
#' @return OPG: the variance covariance matrix of b and g
#' @import splines2
#' @import ICsurv
#' @import numDeriv
#' @references{
#' Withana Gamage, P., McMahan, C., and Wang, L. (2022+).
#' \emph{A Flexible Parametric Model for Fitting the Proportional Hazards Model}.
#' Submitted.
#' }
#' @importFrom Rdpack reprompt
#' @export

UnifiedPH.EM<-function(d1, d2, d3, Li, Ri,Ei, Xp, n.int, order, g0, b0, tol, t.seq, equal = FALSE){
 
  P <- length(b0)
  L <- length(g0)
  N <- length(d1)

  K<-n.int+order
 

  Ri[d3 == 1] <- Li[d3 == 1]
  ti <- c(Li, Ri[d3 == 0])
  #ti.new<-ti[ti>min(Ei)]
  ti.new<-ti
  if (equal == TRUE) {
    ti.new.max <- max(ti.new) + 1e-05
    ti.new.min <- min(ti.new) - 1e-05
    knots <- seq(ti.new.min, ti.new.max, length.out = (n.int + 2))
  }
  if (equal == FALSE) {
    id <- seq(0, 1, length.out = (n.int + 2))
    id <- id[-c(1, (n.int + 2))]
    ti.new.max <- max(ti.new) + 1e-05
    ti.new.min <- min(ti.new) - 1e-05
    knots <- c(ti.new.min, quantile(ti.new, id), ti.new.max)
  }
 
 
 
  ###############################################################################
  # Adding dependent packages
 
usethis::use_package("splines2")
usethis::use_package("ICsurv")
usethis::use_package("numDeriv")

  B1<-t(Ispline(x=Li,order=order,knots=knots))
  B2<-t(Ispline(x=Ri,order=order,knots=knots))
  BEi<-t(Ispline(x=Ei,order=order,knots=knots))
 
 
  bRi<-array(-1000,c(N,order+n.int))
  bLi<-array(-1000,c(N,order+n.int))
 
 # Accounting for enrollment issue

  for(i in 1:N){
    bRi[i,] <- B2[i,]-BEi[i,]
    bLi[i,] <- B1[i,]-BEi[i,]
  }


  bt <- t(Ispline(x = t.seq,order=order, knots=knots))
 
int.knots<-knots[2:(n.int+1)] # only the interior knots use in iSpline function

# getting M-splines as the derivative of I-splines
mRi <- iSpline(x=Ri, knots = int.knots, degree = (order-1), intercept = TRUE,  derivs = 1)

 
  ###############################################################################
  ###############################################################################
 
  diff0=1;
  while (diff0 > tol) {
   
    xb0 <- Xp %*% b0
    EWil=matrix(0, nrow=N, ncol=L)
    dw <- 1 - exp(-(bRi[d2==1,] %*% g0 - bLi[d2==1,] %*% g0) * exp(xb0[d2==1]))
    EWil[d2==1,] <- t(t(bRi[d2==1,] - bLi[d2==1,]) * g0) * as.vector(exp(xb0[d2==1])/dw);

    # for exactly observed
    EVil=matrix(0, nrow=N, ncol=L)
    ail<-t(t(mRi) * g0)
    dv<-as.matrix(rowSums(ail))
    EVil[d1==1,]<-as.matrix(ail[d1==1,]/dv[d1==1])

    Dil <- EVil + EWil;
    Dl <- colSums(Dil);
    Cil <-  (1-d3) * bRi + d3 * bLi;

   
    diff1=1
    b01<-b0
    while (diff1>1e-7)
    {
      xb0 <- Xp %*% b01
      Fl <- colSums(Cil*as.vector(exp(xb0)));
      deriv1<- colSums(t(Dil) %*% Xp) - colSums(t(Cil* exp(xb0)%*%(Dl /Fl)) %*%Xp)  
      deriv2 = t(as.vector(exp(xb0))*Xp) %*%( Cil%*%(as.vector(Dl/Fl^2)*(t(Cil)%*%  (as.vector(exp(xb0))*Xp) ) ))  -  t(as.vector(Cil%*%( Dl/Fl^2*( t(Cil)%*% exp(xb0))))*Xp)%*%(as.vector(exp(xb0))*Xp)
      b11= b01 - solve(deriv2)%*% deriv1
      diff1=max(abs(b11-b01))
      b01=as.vector(b11)        
    }  
    b1=as.vector(b11)
    g1= Dl/(t(Cil)%*%exp(Xp%*%b1))
    diff0=max(abs(c(b0, g0) - c(b1, g1)))
    g0=as.vector(g1)
    b0=as.vector(b1)              
    }
 
  g1=as.vector(g1)
  GRi <- bRi %*% g1
  GLi <- bLi %*% g1
  xb <- Xp %*% b1

  mGRi<-mRi %*% g1

  alllk=(mGRi*exp(xb)*exp(-GRi * exp(xb)))*d1 + (exp(-GLi * exp(xb)) - exp(-GRi * exp(xb)))*d2 + (exp(-GLi * exp(xb)))*d3
  ll<- sum(log(alllk))
  AIC <- 2 * (length(b1) + length(g1)) - 2 * ll
  BIC <- (length(b1) + length(g1)) * log(N) - 2 * ll







indc<-function(theta,d1,d2,d3,B1,B2,mR,xc){
        beta<-theta[1:P]
        g<-theta[(P+1):(P+K)]
     
        g1=as.vector(g[1:K])
    	  b1=as.vector(beta)

    GRi <- B2 %*% g1
    GLi <- B1 %*% g1
    xb <- xc %*% b1

    mGRi<-mR %*% g1

    alllk=(mGRi*exp(xb)*exp(-GRi * exp(xb)))*d1+ (exp(-GLi * exp(xb)) - exp(-GRi * exp(xb)))*d2 + (exp(-GLi * exp(xb)))*d3

    res<-sum(log(alllk))
    return(res)
      }


gr.theta<-function(theta,d1,d2,d3,B1,B2,mR,xc){
        N<-length(d1)
        grad.mat<-matrix(-99,nrow=length(theta),ncol=N)
        
        for(i in 1:N){
          grad.mat[,i]<-grad(func=indc,x=theta,d1=d1[i],d2=d2[i],d3=d3[i],B1=B1[i,],B2=B2[i,],mR=mR[i,],xc=xc[i,])
        }
        return(grad.mat)
      }
      
theta.hat<-c(b1,g1)
      gr<-gr.theta(theta.hat[1:(P+K)],d1=d1,d2=d2,d3=d3,B1=bLi,B2=bRi,mR=mRi,xc=Xp)
      Bn<-gr%*%t(gr)
      v.OPG<-solve(Bn)








  return(list("b"=b1,"g"=g1,"ll"=ll,"AIC"=AIC,"BIC"=BIC,"bRi"=bRi,"bLi"=bLi,"bt"=bt,"mRi"=mRi,"OPG"=v.OPG))
}
