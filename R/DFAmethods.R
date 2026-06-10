#' DFATools: detrended fluctuation analysis and related methods
#'
#' A toolbox of Detrended Fluctuation Analysis ('DFA') methods for long-range
#' correlations, cross-correlations and regression in nonstationary time series.
#' See \code{vignette("DFATools")} for the theory and worked examples.
#'
#' @section Function groups:
#' \describe{
#'   \item{Scaling of a single series}{\code{\link{dfa}}, \code{\link{plotdfa}}}
#'   \item{Cross-correlation}{\code{\link{rhodcca}}, \code{\link{plotrdcca}},
#'     \code{\link{plotdcca}}}
#'   \item{Partial cross-correlation}{\code{\link{rhodpcca}}}
#'   \item{Multiple cross-correlation}{\code{\link{dmc2}}}
#'   \item{Detrended fractal regression}{\code{\link{betadfa}},
#'     \code{\link{sbdfa}}, \code{\link{fracreg}}, \code{\link{effsizeDFA}}}
#'   \item{Significance tests}{\code{\link{fracreg.PStest}},
#'     \code{\link{fracreg.Ktest}}, \code{\link{fracreg.IUTest}}}
#' }
#'
#' @keywords internal
"_PACKAGE"

# Column names referenced via non-standard evaluation in plotrdcca().
utils::globalVariables(c("s", "rho"))

# ---------------------------------------------------------------------------
# Internal input-validation helpers. They stop() with a clear message on the
# common user mistakes (wrong type, NA, too few columns, series too short,
# invalid np/dpo/int/overlap/B) before any numerical work is attempted.
# ---------------------------------------------------------------------------

# Validate the shared arguments np, dpo, int and overlap.
.check_common <- function(np, dpo, int, overlap, fn) {
  if (!is.numeric(np) || length(np) != 1L || is.na(np) || np < 2 || np != round(np))
    stop(sprintf("%s(): `np` must be a single integer >= 2.", fn), call. = FALSE)
  if (!is.numeric(dpo) || length(dpo) != 1L || is.na(dpo) || dpo < 1 || dpo != round(dpo))
    stop(sprintf("%s(): `dpo` must be a single integer >= 1.", fn), call. = FALSE)
  if (!is.logical(int) || length(int) != 1L || is.na(int))
    stop(sprintf("%s(): `int` must be a single TRUE or FALSE.", fn), call. = FALSE)
  if (!is.logical(overlap) || length(overlap) != 1L || is.na(overlap))
    stop(sprintf("%s(): `overlap` must be a single TRUE or FALSE.", fn), call. = FALSE)
  invisible(TRUE)
}

# The box sizes run from a minimum of 10 up to round(n / 5); the series must be
# long enough for that maximum to exceed the minimum.
.check_length <- function(nobs, fn) {
  if (round(nobs / 5) <= 10)
    stop(sprintf(paste0("%s(): series too short (%d observations). At least ~55 ",
                        "are needed so the largest box size exceeds the smallest (10)."),
                 fn, nobs), call. = FALSE)
  invisible(TRUE)
}

# Validate a single numeric series (for dfa()).
.check_series <- function(x, fn) {
  if (!is.numeric(x) || !is.null(dim(x)))
    stop(sprintf("%s(): `data` must be a numeric vector.", fn), call. = FALSE)
  if (anyNA(x))
    stop(sprintf("%s(): `data` contains missing values (NA); remove them first.", fn),
         call. = FALSE)
  .check_length(length(x), fn)
  invisible(TRUE)
}

# Validate a numeric matrix/data frame with at least `min_cols` columns.
.check_matrix <- function(data, min_cols, fn) {
  if (is.null(dim(data)) || length(dim(data)) != 2L)
    stop(sprintf("%s(): `data` must be a matrix or data frame with at least %d columns.",
                 fn, min_cols), call. = FALSE)
  if (is.data.frame(data)) {
    ok <- vapply(data, is.numeric, logical(1))
    if (!all(ok))
      stop(sprintf("%s(): all columns of `data` must be numeric; non-numeric: %s.",
                   fn, paste(names(data)[!ok], collapse = ", ")), call. = FALSE)
  } else if (!is.numeric(data)) {
    stop(sprintf("%s(): `data` must be numeric.", fn), call. = FALSE)
  }
  if (ncol(data) < min_cols)
    stop(sprintf("%s(): `data` needs at least %d columns; got %d.",
                 fn, min_cols, ncol(data)), call. = FALSE)
  if (anyNA(data))
    stop(sprintf("%s(): `data` contains missing values (NA); remove or impute them first.", fn),
         call. = FALSE)
  .check_length(nrow(data), fn)
  invisible(TRUE)
}

# Validate the surrogate count B used by the hypothesis tests.
.check_B <- function(B, fn) {
  if (!is.numeric(B) || length(B) != 1L || is.na(B) || B < 1 || B != round(B))
    stop(sprintf("%s(): `B` must be a single integer >= 1.", fn), call. = FALSE)
  invisible(TRUE)
}

#' Multiple Fractal Regression
#'
#' Calculates the scale-dependent (DFA-based) multiple linear regression: the
#' coefficients, their variance and confidence intervals at each scale, together
#' with scale-wise collinearity diagnostics.
#'
#' @details
#' The variance of the scale-dependent coefficients follows Tilfani et al.
#' (2022, Eq. 25), \eqn{\mathrm{var}(\hat\beta_j(s)) = F^2_\varepsilon(s)
#' [F_{XX}(s)^{-1}]_{jj}}, i.e. it uses the full inverse of the detrended
#' covariance matrix of the predictors, consistent with how the coefficients
#' themselves are estimated. Under collinearity this differs from the legacy
#' "marginal" form \eqn{F^2_\varepsilon(s) / F^2_{X_j}(s)} (which under-covers);
#' the two coincide for orthogonal predictors. The default
#' \code{variance = "inv_corrected"} multiplies the inverse form by the
#' memory-correction factor \eqn{c(\widehat H) = (2\widehat H + 1)^2 / 21}
#' (Barreto et al. 2026, Eq. M4.2'), with \eqn{\widehat H} the DFA exponent of
#' the OLS residual, restoring nominal coverage under polynomial detrending;
#' \code{variance = "inv"} omits the factor and \code{variance = "marginal"}
#' reproduces the legacy form. The convention is: \eqn{\kappa(s, H) =
#' F^2(s,H)^2 / \mathrm{Var}(f^2_{X\varepsilon}(s,\nu))} is the ratio of
#' bilinears that the closed form approximates (\eqn{\kappa \in [2.3, 5]}),
#' and \eqn{c(H) = 1/\kappa(H) \approx (2H+1)^2/21} is its inverse, the
#' multiplier actually applied to the variance (\eqn{c \in [0.21, 0.43]}). The
#' applied multiplier is returned in \code{$c_factor}.
#'
#' The variance is normalised by the residual degrees of freedom
#' \eqn{T_s - k}, where \eqn{T_s = \lfloor N/s \rfloor} is the number of
#' non-overlapping boxes at scale \eqn{s} and \eqn{k} the number of predictors;
#' the same \eqn{T_s - k} is the degrees of freedom of the \code{t} quantile.
#' This differs from \eqn{T_s - k - 1} used by Tilfani et al. (2022): the
#' scale-wise intercept \eqn{\hat\beta_0(s) = \bar y - \sum_j \hat\beta_j \bar x_j}
#' is derived from the slope estimates rather than fit jointly across boxes, so
#' it does not consume an additional degree of freedom.
#' \eqn{T_s} counts disjoint boxes regardless of \code{overlap}. Scales with
#' \eqn{T_s \le k} return \code{NA} limits with a warning. The analytic interval
#' can under-cover under strong long-range dependence (estimated DFA exponent
#' above \eqn{3/4}, the Hermite-Rosenblatt threshold); a warning is then issued
#' and \code{\link{fracreg.WB}} (dependent wild bootstrap) should be preferred.
#'
#' @param data a matrix or data frame of time series; the first column is the
#'   response and the remaining columns are the predictors.
#' @param dpo detrending polynomial order.
#' @param int logical. If TRUE the integration process is applied.
#' @param np number of point scales.
#' @param overlap logical. If TRUE overlapping windows are used. Defaults to
#'   FALSE: score-based inference treats the boxes as sampling units and requires
#'   them disjoint. With \code{overlap = TRUE} only point estimates are returned
#'   (\code{variance = "none"}).
#' @param variance coefficient-variance estimator:
#'   \code{"inv_corrected"} (default) the memory-corrected inverse form
#'   \eqn{F^2_\varepsilon(s)[F_{XX}(s)^{-1}]_{jj}/(T_s-k)\cdot(2\widehat H+1)^2/21}
#'   (Barreto et al. 2026, Eq. M4.2'); \code{"inv"} the uncorrected inverse
#'   (Tilfani et al. 2022); \code{"marginal"} the legacy \eqn{1/F^2_{X_j}(s)}
#'   form (Shen 2015), which under-covers under collinearity; \code{"hc"} the
#'   heteroscedasticity-consistent sandwich; or \code{"none"} for point estimates
#'   only.
#' @param H_eps optional numeric. A pre-computed DFA exponent of the regression
#'   error to use in the memory correction; if \code{NULL} (default) it is
#'   estimated from the OLS residual.
#' @param min_boxes minimum number of non-overlapping boxes
#'   \eqn{T_s = \lfloor N/s\rfloor} required for inference at a scale (default
#'   15). Scales with \eqn{T_s < } \code{min_boxes} return \code{NA} standard
#'   errors, confidence limits and p-values, with a single warning listing them.
#'   A warning is also issued when \eqn{N < 500} (Likens et al. 2019).
#' @param abs logical. If TRUE the absolute detrended covariance is used in the
#'   cross-correlation step (more robust to outliers).
#' @return A list with, among others: scale \code{s}, detrended fluctuation
#'   \code{F}, \code{DCCA}, \code{DPCCA}, beta estimates \code{BDFA},
#'   standardized betas \code{BSDFA}, residual variance \code{UDFA}, coefficient
#'   variance \code{VDFA}, multiple correlation \code{DMC2}, \code{R2DFA},
#'   confidence limits \code{UCIB}/\code{LCIB}, \code{p.value}, critical value
#'   \code{TC}, the scale-wise diagnostics \code{VIF}, condition number
#'   \code{kappa} and adjusted \code{R2adj}, the per-series DFA exponent
#'   \code{alpha}, the residual DFA exponent \code{H_resid} and memory factor
#'   \code{c_factor} used in the correction, the chosen
#'   \code{variance_method}, the effective degrees of freedom \code{df_eff} and
#'   the box counts \code{Ts}.
#' @references
#' Barreto, I. D. C., Dore, L. H., Stosic, T. and Stosic, B. D. (2021).
#' Extending DFA-based multiple linear regression inference: application to
#' acoustic impedance models. \emph{Physica A}, 582, 126259.
#'
#' Tilfani, O., Kristoufek, L., Ferreira, P. and El Boukfaoui, M. Y. (2022).
#' Heterogeneity in economic relationships: scale dependence through the
#' multivariate fractal regression. \emph{Physica A}, 588, 126530.
#'
#' Shen, C. (2015). A new detrended semipartial cross-correlation analysis.
#' \emph{Physics Letters A}, 379(44), 2962-2969.
#'
#' Kristoufek, L. (2015). Detrended fluctuation analysis as a regression
#' framework: estimating dependence at different scales. \emph{Physical Review
#' E}, 91(2), 022802.
#' @seealso \code{\link{fracreg.diag}}, \code{\link{betadfa}},
#' \code{\link{effsizeDFA}}, \code{vignette("DFATools")}
#' @export
#' @importFrom stats pt
#' @importFrom stats qt
#' @importFrom stats cov2cor
#' @importFrom corpcor cor2pcor
#' @importFrom stats lm
#' @importFrom combinat combn
#' @importFrom combinat nCm
#' @importFrom stats quantile
#' @importFrom tseries surrogate
#' @importFrom magrittr "%>%"
#' @examples
#' set.seed(1)
#' x1 <- rnorm(400); x2 <- rnorm(400)        # stationary predictors
#' d <- data.frame(y = 0.7 * x1 - 0.5 * x2 + rnorm(400), x1 = x1, x2 = x2)
#' fit <- fracreg(d, dpo = 1, int = TRUE, np = 20, overlap = FALSE)
#' round(fit$BDFA[, 1, 10], 2)               # coefficients at the 10th scale
#' @useDynLib DFATools, .registration=TRUE

fracreg<-function(data,dpo=1,int=TRUE,np=91,overlap=FALSE,
                  variance=c("inv_corrected","inv","marginal","hc","none"),
                  H_eps=NULL,min_boxes=15,abs=FALSE){
	.check_common(np, dpo, int, overlap, "fracreg")
	.check_matrix(data, 2, "fracreg")
	variance<-match.arg(variance)
	if(!is.logical(abs)||length(abs)!=1L||is.na(abs))
		stop("fracreg(): `abs` must be a single TRUE or FALSE.", call.=FALSE)
	if(!is.null(H_eps)&&(!is.numeric(H_eps)||length(H_eps)!=1L))
		stop("fracreg(): `H_eps` must be NULL or a single number.", call.=FALSE)
	if(!is.numeric(min_boxes)||length(min_boxes)!=1L||!is.finite(min_boxes)||min_boxes<1)
		stop("fracreg(): `min_boxes` must be a single positive number.", call.=FALSE)
	if(variance!="none"&&nrow(data)<500L)
		warning("fracreg(): N = ",nrow(data)," < 500; DFA-based inference may be ",
			"unreliable (Likens et al., 2019).",call.=FALSE)
	# Score-based inference treats the boxes as sampling units and requires them
	# disjoint; overlapping windows invalidate the standard errors (paper M8).
	if(isTRUE(overlap)&&variance!="none"){
		message("fracreg(): overlap = TRUE invalidates score-based standard errors; ",
			"returning point estimates only (variance = \"none\"). ",
			"Use overlap = FALSE for valid confidence intervals.")
		variance<-"none"
	}
	data<-as.matrix(data)
	dpo0<-as.numeric(dpo)
	dpo<-as.numeric(dpo+1)
	if (int ==TRUE){int=1} else{int=0}
	nc<-ncol(data)
	pairs<-t(combn(nc,2))
	comb<-nCm(nc,2)
	sn<-matrix(,ncol=1) #Escalas
	mx<-round(nrow(data)/5)
	fn<-array(,dim=c(nc,nc,np)) #Variancia-Covarância sem tendência
	bn<-array(,dim=c((nc-1),1,np)) #Parâmetros Bdfa(s)
	uci<-array(,dim=c((nc-1),1,np)) #Limite Superior de Confiança Bdfa(s)
	lci<-array(,dim=c((nc-1),1,np)) #Limite Superior de Confiança Bdfa(s)
	un<-array(,dim=c(1,1,np)) #Variância sem tendência do resíduo
	rn<-array(,dim=c(1,1,np)) #R2(s)
	dmc2<-matrix(NA,nrow=np) #DMC2 Zebende et al (2018)
	bs<-array(,dim=c(1,(nc-1),np)) #Beta Padronizado
	vn<-array(,dim=c(1,(nc-1),np)) #Variância dos Parâmetros Bdfa(s)
	vn2<-array(,dim=c(1,(nc-1),np)) #Variância dos Parâmetros Bdfa(s)
	tn<-array(,dim=c((nc-1),1,np)) #p-valor para os Bdfas
	tnc<-array(,dim=c((nc-1),1,np)) #T crítico
	VIFn<-array(,dim=c(1,(nc-1),np)) #scale-wise VIF
	kappan<-array(,dim=c(1,1,np)) #scale-wise condition number
	r2adj<-array(,dim=c(1,1,np)) #scale-wise adjusted R2
	size=nrow(data)
#DFA
	for (i in 1:nc)
	{
seq1<-data[,i]
# ****** CONFIGURATION STRUCTURE
NPTS <- length(seq1) #size of series
NFIT	<- dpo		# order of the regression fit, plus 1
IFLAG <- as.numeric(int)		# integrate the input data if non-zero
NR 	<- as.numeric(np)		# number of box sizes
SW 	<- as.numeric(overlap)		# sliding window
MINBOX	<- 10		# minimum box size
MAXBOX	<- as.numeric(mx)		# maximum box size
cfg <- as.integer(cbind(NPTS,NFIT,IFLAG,NR,SW,MINBOX,MAXBOX))
rsi1<- numeric(NR+1)
msei1<- numeric(NR+1)
ans1<- .C("rdfa", cfg, as.numeric(seq1), as.integer(rsi1), as.numeric(msei1),PACKAGE = "DFATools")
sn<- ans1[[3]][1:NR+1]
fn[i,i,]<-ans1[[4]][1:NR+1]
}
#DCCA
#fc<-matrix(,ncol=comb)
for (i in 1:comb){
seq1<-data[,pairs[i,1]]
seq2<-data[,pairs[i,2]]
# ****** CONFIGURATION STRUCTURE
NPTS <- length(seq1) #size of series
NFIT	<- dpo		# order of the regression fit, plus 1
IFLAG <- as.numeric(int)		# integrate the input data if non-zero
NR 	<- as.numeric(np)	# number of box sizes
SW 	<- as.numeric(overlap)		# sliding window
MINBOX	<- 10		# minimum box size
MAXBOX	<- as.numeric(mx)	# maximum box size
ABSFLAG	<- as.numeric(abs)		# absolute flag for dcca
cfg <- as.integer(cbind(NPTS,NFIT,IFLAG,NR,SW,MINBOX,MAXBOX,ABSFLAG))
rs12 <- numeric(NR+1)
mse12 <- numeric(NR+1)
ans12 <- .C("rdcca", cfg, as.numeric(seq1), as.numeric(seq2), as.integer(rs12), as.numeric(mse12),PACKAGE = "DFATools")
fn[pairs[i,1],pairs[i,2],] <- ans12[[5]][1:NR+1]
fn[pairs[i,2],pairs[i,1],] <- ans12[[5]][1:NR+1]
}
#DPCCA
pn<-array(,dim=c(nc,nc,np)) #pDCCA
dp<-array(,dim=c(nc,nc,np)) #pDPCCA
for(i in 1:(np)){
	pn[,,i]<-cov2cor(fn[,,i])
	dp[,,i]<-cor2pcor(pn[,,i])
	}
###Fractal Regression
for(k in 1:np){
		bn[,1,k]<-solve(fn[(2:nc),(2:nc),k])%*%fn[(2:nc),1,k]
			}
for(k in 1:np){
u1<-data[,1]-data[,(2:nc)]%*%as.matrix(bn[,,k])
u<-scale(u1,scale=FALSE)
NPTS <- length(u) #size of series
NFIT	<- dpo		# order of the regression fit, plus 1
IFLAG <- 1		# integrate the input data if non-zero
NR 	<- 1		# number of box sizes
SW 	<- as.numeric(overlap)		# sliding window
MINBOX	<- sn[[k]]		# minimum box size
MAXBOX	<- sn[[k]]		# maximum box size
cfg <- as.integer(cbind(NPTS,NFIT,IFLAG,NR,SW,MINBOX,MAXBOX))
rsi1<- numeric(NR+1)
msei1<- numeric(NR+1)
ans1<- .C("rdfa", cfg, as.numeric(u), as.integer(rsi1), as.numeric(msei1),PACKAGE = "DFATools")
un[1,1,k]<-ans1[[4]][2]
dmc2[k]<-pn[1,(2:nc),k]%*%solve(pn[(2:nc),(2:nc),k])%*%pn[(2:nc),1,k]
bs[,,k]<-solve(pn[(2:nc),(2:nc),k])%*%pn[(2:nc),1,k]
rn[1,1,k]<-1-(un[1,1,k]/fn[1,1,k])
u<-NULL
u1<-NULL
}
for(k in 1:np){
	if(k==1L){
		# Memory-correction factor (Barreto et al. 2026, Eq. M4.2'):
		# kappa(H) = (2H+1)^2 / 21, with H the DFA-1 exponent of the OLS RESIDUAL
		# (the error's long memory drives the score CLT, not the predictors').
		if(is.null(H_eps)){
			ols_res<-as.numeric(stats::residuals(stats::lm(data[,1]~data[,-1,drop=FALSE])))
			rd<-dfa(ols_res,dpo=dpo0,int=TRUE,np=np,overlap=FALSE)
			ok<-rd$s>=max(10,min(rd$s))&rd$s<=floor(size/10)&rd$F>0
			# dfa()$F follows the Peng convention (= sqrt(F^2)), so the DFA exponent
			# is the slope of log F against log s directly.
			H_resid<-if(sum(ok)<4) NA_real_ else
				stats::coef(stats::lm(log(rd$F[ok])~log(rd$s[ok])))[[2]]
		}else H_resid<-H_eps
		H_raw<-H_resid                              # before clamp/fallback
		if(is.na(H_raw)||H_raw<0.3||H_raw>1.2){
			if(variance=="inv_corrected")
				warning("fracreg(): could not reliably estimate the residual DFA exponent; ",
					"using c = 1 (no memory correction).",call.=FALSE)
			H_resid<-NA_real_; c_factor<-1
		}else{
			H_resid<-min(max(H_raw,0.5),1-1e-6)     # clamp mild excursions into [0.5, 1)
			c_factor<-(2*H_resid+1)^2/21
		}
	}
	F2eps<-fn[1,1,k]-dmc2[k]*fn[1,1,k]   # detrended residual variance F^2_eps(s)
	# Coefficient (co)variance matrix of the predictors and its inverse (reused
	# from the beta estimate). The variance of beta_j(s) is, per Tilfani et al.
	# (2022, Eq. 25), F^2_eps(s) * [F_XX(s)^-1]_jj -- the FULL inverse, which the
	# legacy "marginal" form (1 / F^2_Xj(s)) only matches under orthogonality.
	Mk<-solve(matrix(fn[(2:nc),(2:nc),k],nc-1,nc-1))
	# Scale-wise collinearity diagnostics from the correlation matrix.
	Pk<-matrix(pn[(2:nc),(2:nc),k],nc-1,nc-1)
	VIFn[1,,k]<-diag(solve(Pk))                      # = 1/(1 - R^2_j(s))
	ev<-eigen(Pk,only.values=TRUE)$values
	kappan[1,1,k]<-max(ev)/min(ev)                   # condition number
	# Adjusted R^2 (Tilfani Eq. 26). F^2_eps/F^2_y is taken as 1 - DMC(s), the
	# model-based residual fraction, consistent with the variance estimator and
	# bounded to (-Inf, 1] (the direct residual DFA can be unstable at extreme
	# scales).
	Tsk<-floor(size/sn[[k]])                         # M3.1: T_s-based adjusted R^2
		r2adj[1,1,k]<-1-((Tsk-1)/(Tsk-(nc-1)-1))*(1-dmc2[k])
	for(i in (2:nc)){
	  dfk<-floor(size/sn[[k]])-(nc-1)
	  base_j<-switch(variance,
	    inv_corrected = Mk[(i-1),(i-1)]*c_factor,   # M4.2': inverse x kappa(H)
	    inv           = Mk[(i-1),(i-1)],                # M4.2: inverse, no correction
	    marginal      = 1/fn[i,i,k],                    # M4.1: legacy marginal
	    Mk[(i-1),(i-1)])                                # hc/none: placeholder
	  vn[,(i-1),k]<-F2eps*base_j/dfk
	  vn2[,(i-1),k]<-(un[1,1,k]%*%solve(fn[i,i,k]))
		}}
if(variance=="hc"){
	# Heteroscedasticity-consistent (sandwich) variance, per scale s:
	# Var_HC = (1/Ts) F_XX^-1 Omega F_XX^-1, Omega = (1/Ts) sum_v r_v r_v' (HC1),
	# from the per-box moment scores r_v = s_xy(s,v) - S_xx(s,v) beta_hat.
	# Non-overlapping boxes (overlap = FALSE) are recommended for the sandwich.
	pb<-.fracreg_perbox(data,dpo=dpo-1,int=as.logical(int),np=np,overlap=overlap,abs=abs)
	for(kk in 1:np){
		ps<-pb$perscale[[kk]]; Tsk<-ps$Ts
		Om<-crossprod(ps$scores)/Tsk
		if(Tsk>(nc-1)) Om<-Om*Tsk/(Tsk-(nc-1))      # HC1 finite-sample correction
		Vh<-(ps$Fxx_inv%*%Om%*%ps$Fxx_inv)/Tsk
		vn[1,,kk]<-diag(Vh)
	}
}
# Score-based inference requires both df > 0 and enough non-overlapping boxes
# (min_boxes). The two gates are combined into ok_scale[k]; suppressed scales
# return NA on SE, CI and p, and are reported in a single warning below.
Tsv_loop<-floor(size/sn); ok_scale<-(Tsv_loop>=min_boxes)&(Tsv_loop>(nc-1))
if(variance=="none") ok_scale[]<-FALSE
for(k in 1:np){
	for(i in 1:(nc-1)){
	dfk<-Tsv_loop[k]-(nc-1)
	vn[1,i,k]<-ifelse(ok_scale[k],vn[1,i,k],NA)
	uci[i,1,k]<-as.matrix(bn[i,1,k])+qt(0.975,df=ifelse(ok_scale[k],dfk,NA))*sqrt(vn[1,i,k])
	lci[i,1,k]<-as.matrix(bn[i,1,k])-qt(0.975,df=ifelse(ok_scale[k],dfk,NA))*sqrt(vn[1,i,k])
	tn [i,1,k]<-as.matrix(2*(1-pt(abs(bn[i,1,k])/sqrt(vn[1,i,k]),df=ifelse(ok_scale[k],dfk,NA))))
	tnc[i,1,k]<-as.matrix(qt(0.975,df=ifelse(ok_scale[k],dfk,NA))*sqrt(vn[1,i,k]))
	}}

	# Per-series DFA exponents (informative). The OLS-residual exponent H_resid
	# (estimated above) is what drives the memory correction and the CLT regime.
	alpha<-vapply(seq_len(nc),function(ii)
		stats::coef(stats::lm(0.5*log(fn[ii,ii,])~log(sn)))[[2]],numeric(1))
	names(alpha)<-colnames(data)
	# Gradient of advisories on the residual exponent. The closed-form factor
	# c(H) = (2H+1)^2/21 is calibrated for H in [0.5, 0.95]; below 0.5 it has
	# ~20% error and the future variance = "inv_theoretical" (table-based)
	# should be preferred; between 3/4 and 0.95 the analytic CI starts to
	# under-cover (Hermite-Rosenblatt regime); above 0.95 the series is close
	# to non-stationary; saturation near 1 is usually a sign of a missing
	# trend or omitted covariate.
	if(variance!="none"&&!is.na(H_resid)){
		Hr<-round(H_resid,2)
		if(H_resid>=0.99)
			warning("fracreg(): H_resid = ",Hr," saturated near 1 -- possible ",
				"non-stationary regime or omitted variable with long memory. ",
				"Consider (i) a higher detrending order; (ii) the dependent ",
				"wild bootstrap (fracreg.WB(weights='dependent')); ",
				"(iii) auditing the model for exogeneity.",call.=FALSE)
		else if(H_resid>=0.95)
			message("fracreg(): H_resid = ",Hr," close to the non-stationary ",
				"regime; the analytic CIs are conservative. Consider ",
				"fracreg.WB(weights='dependent') for inference.")
		else if(H_resid>0.75)
			warning("fracreg(): H_resid = ",Hr," exceeds 3/4 (Hermite-Rosenblatt ",
				"threshold); the analytic interval can under-cover under strong ",
				"long-range dependence. Prefer fracreg.WB(weights='dependent').",
				call.=FALSE)
		else if(!is.na(H_raw)&&H_raw<0.5)
			message("fracreg(): H_resid = ",round(H_raw,2)," below 0.5; the closed-",
				"form memory factor c(H) = (2H+1)^2/21 has up to ~20% error in ",
				"this range (calibrated for H in [0.5, 0.95]); H was clamped at ",
				"0.5 for the correction.")
	}
	if(variance!="none"){
		few_idx<-which(Tsv_loop<min_boxes|Tsv_loop<=nc-1)
		if(length(few_idx)){
			few_s<-sn[few_idx]
			lst<-if(length(few_s)>6) paste(c(utils::head(few_s,5),"..."),collapse=", ") else
				paste(few_s,collapse=", ")
			warning("fracreg(): scales {",lst,"} have fewer than min_boxes = ",min_boxes,
				" non-overlapping boxes (T_s = floor(N/s)); their standard errors and ",
				"intervals were set to NA. Reduce the maximum scale or increase N.",
				call.=FALSE)
		}
	}
	Tsv<-floor(size/sn)                          # non-overlapping box count per scale
	df_eff<-Tsv-(nc-1)                           # residual degrees of freedom T_s - k
	c_factor_v<-rep(c_factor,np)         # memory factor (constant across s)
fracreg<-list(sn,fn,pn,dp,bn,bs,un,vn,vn2,dmc2,rn,uci,lci,tn,tnc,VIFn,kappan,r2adj,alpha,
	H_resid,c_factor_v,variance,df_eff,Tsv)
names(fracreg)<-c("s","F","DCCA","DPCCA","BDFA","BSDFA","UDFA","VDFA","VDFA2","DMC2","R2DFA","UCIB","LCIB","p.value","TC","VIF","kappa","R2adj","alpha","H_resid","c_factor","variance_method","df_eff","Ts")

return(fracreg)
}

#' Scale-wise collinearity diagnostics for the fractal regression
#'
#' Returns scale-dependent multicollinearity diagnostics for the DFA-based
#' multiple regression: the variance inflation factors (VIF), the condition
#' number (\code{kappa}) of the predictors' scale-wise correlation matrix, and
#' the scale-wise adjusted coefficient of determination. These reveal
#' multicollinearity that depends on the time scale.
#'
#' @param data a matrix or data frame of time series; the first column is the
#'   response and the remaining columns are the predictors.
#' @param dpo detrending polynomial order.
#' @param int logical. If TRUE the integration process is applied.
#' @param np number of point scales.
#' @param overlap logical. If TRUE overlapping windows are applied.
#' @param abs logical. If TRUE the absolute detrended covariance is used.
#' @return A tibble with the scale \code{s}, one \code{VIF_*} column per
#'   predictor, the condition number \code{kappa} and the adjusted \code{R2adj},
#'   one row per scale.
#' @references
#' Tilfani, O., Kristoufek, L., Ferreira, P. and El Boukfaoui, M. Y. (2022).
#' Heterogeneity in economic relationships: scale dependence through the
#' multivariate fractal regression. \emph{Physica A}, 588, 126530.
#' @seealso \code{\link{fracreg}}, \code{vignette("DFATools")}
#' @importFrom tibble as_tibble
#' @examples
#' set.seed(1)
#' d <- data.frame(y = cumsum(rnorm(300)), x1 = cumsum(rnorm(300)),
#'                 x2 = cumsum(rnorm(300)))
#' fracreg.diag(d, np = 20)
#' @export
fracreg.diag <- function(data, dpo = 1, int = TRUE, np = 91, overlap = TRUE, abs = FALSE) {
  fit <- fracreg(data, dpo = dpo, int = int, np = np, overlap = overlap,
                 variance = "none", abs = abs)
  nc <- ncol(as.matrix(data))
  cn <- colnames(data)
  if (is.null(cn)) cn <- paste0("x", seq_len(nc))
  vif <- t(matrix(fit$VIF[1, , ], nrow = nc - 1))
  colnames(vif) <- paste0("VIF_", cn[2:nc])
  out <- cbind.data.frame(s = fit$s, vif,
                          kappa = as.vector(fit$kappa),
                          R2adj = as.vector(fit$R2adj))
  tibble::as_tibble(out)
}

# Internal: per-box detrended (co)variance scores for the fractal regression.
# Shared core for variance = "hc" and fracreg.WB(). For each scale it returns the box
# size, the number of boxes Ts, the inverse F_XX(s)^-1, the coefficient vector
# beta_hat, and the Ts x k matrix of per-box moment scores whose v-th row is
# r_v' = s_xy(s,v) - (S_xx(s,v) beta_hat)'. Uses the rdfa_box/rdcca_box C
# primitives (one detrending pass, exposed per box).
.fracreg_perbox <- function(data, dpo = 1, int = TRUE, np = 91, overlap = TRUE, abs = FALSE) {
  data <- as.matrix(data)
  n <- nrow(data); nc <- ncol(data); k <- nc - 1L
  NFIT <- as.numeric(dpo + 1); iflag <- if (isTRUE(int)) 1 else 0
  sw <- as.numeric(overlap); mx <- round(n / 5); absf <- as.numeric(abs)

  rdfaB <- function(z) {
    cfg <- as.integer(cbind(n, NFIT, iflag, np, sw, 10, mx))
    a <- .C("rdfa_box", cfg, as.numeric(z), as.integer(numeric(np + 1)),
            numeric(np + 1), numeric(np * n), PACKAGE = "DFATools")
    list(s = a[[3]][2:(np + 1)], box = matrix(a[[5]], nrow = n))
  }
  rdccaB <- function(z1, z2) {
    cfg <- as.integer(cbind(n, NFIT, iflag, np, sw, 10, mx, absf))
    a <- .C("rdcca_box", cfg, as.numeric(z1), as.numeric(z2),
            as.integer(numeric(np + 1)), numeric(np + 1), numeric(np * n),
            PACKAGE = "DFATools")
    matrix(a[[6]], nrow = n)
  }

  s <- NULL
  diagBox <- vector("list", k)
  for (j in seq_len(k)) { r <- rdfaB(data[, j + 1]); diagBox[[j]] <- r$box; s <- r$s }
  xyBox <- lapply(seq_len(k), function(j) rdccaB(data[, 1], data[, j + 1]))
  pkey <- function(i, l) paste(i, l, sep = "-")
  offBox <- list()
  if (k >= 2) for (i in 1:(k - 1)) for (l in (i + 1):k)
    offBox[[pkey(i, l)]] <- rdccaB(data[, i + 1], data[, l + 1])

  Ts <- if (overlap) n - s + 1 else floor(n / s)

  perscale <- vector("list", np)
  for (ki in seq_len(np)) {
    Tsi <- Ts[ki]; vv <- seq_len(Tsi)
    diag_v <- matrix(vapply(seq_len(k), function(j) diagBox[[j]][vv, ki], numeric(Tsi)),
                     nrow = Tsi, ncol = k)
    xy_v   <- matrix(vapply(seq_len(k), function(j) xyBox[[j]][vv, ki], numeric(Tsi)),
                     nrow = Tsi, ncol = k)
    Fxx <- matrix(0, k, k); diag(Fxx) <- colMeans(diag_v)
    if (k >= 2) for (i in 1:(k - 1)) for (l in (i + 1):k) {
      cv <- mean(offBox[[pkey(i, l)]][vv, ki]); Fxx[i, l] <- cv; Fxx[l, i] <- cv
    }
    Fxx_inv <- solve(Fxx)
    beta <- as.vector(Fxx_inv %*% colMeans(xy_v))
    Sbeta <- sweep(diag_v, 2, beta, `*`)               # diagonal part of S_xx(v) beta
    if (k >= 2) for (i in 1:(k - 1)) for (l in (i + 1):k) {
      cvv <- offBox[[pkey(i, l)]][vv, ki]
      Sbeta[, i] <- Sbeta[, i] + cvv * beta[l]
      Sbeta[, l] <- Sbeta[, l] + cvv * beta[i]
    }
    perscale[[ki]] <- list(s = s[ki], Ts = Tsi, Fxx_inv = Fxx_inv,
                           beta = beta, scores = xy_v - Sbeta)   # Ts x k
  }
  list(np = np, k = k, n = n, scales = s, perscale = perscale)
}

# Internal: Ts x B matrix of wild-bootstrap multiplier weights, E[w]=0, Var[w]=1.
.wb_weights <- function(Ts, B, type, bandwidth) {
  if (type == "rademacher") {
    matrix(sample(c(-1, 1), Ts * B, replace = TRUE), Ts, B)
  } else if (type == "mammen") {
    p <- (sqrt(5) + 1) / (2 * sqrt(5))
    matrix(ifelse(stats::runif(Ts * B) < p, -(sqrt(5) - 1) / 2, (sqrt(5) + 1) / 2),
           Ts, B)
  } else {  # "dependent" wild bootstrap (Shao 2010): Bartlett-kernel Gaussian
    l <- if (is.null(bandwidth)) max(1, Ts^(1 / 3)) else bandwidth
    d <- abs(outer(seq_len(Ts), seq_len(Ts), "-")) / l
    K <- pmax(1 - d, 0)                                   # Bartlett kernel (PSD)
    U <- chol(K + diag(1e-8, Ts))                         # K = U'U
    t(U) %*% matrix(stats::rnorm(Ts * B), Ts, B)
  }
}

#' Wild bootstrap inference for the fractal regression
#'
#' Wild-bootstrap confidence intervals and significance for the scale-dependent
#' regression coefficients. It resamples the per-box detrended-moment scores
#' (the same scores used by \code{fracreg(variance = "hc")}), so it is robust to
#' heteroscedasticity and -- with dependent weights -- to dependence between
#' boxes, and it does not recompute the DFA for each replicate.
#'
#' @details
#' For each scale \eqn{s} with \eqn{T_s} boxes, scores \eqn{r_v} and inverse
#' \eqn{F_{XX}(s)^{-1}}, each replicate is
#' \eqn{\hat\beta^*_b(s) = \hat\beta(s) + F_{XX}(s)^{-1}\,(1/T_s)\sum_v r_v w_{v,b}}.
#' The interval is the 2.5\% / 97.5\% quantiles of \eqn{\hat\beta^*_b(s)}; the
#' two-sided p-value tests \eqn{H_0\!: \beta_j(s) = 0}. The \code{"dependent"}
#' weights (Shao 2010) keep coverage when the boxes are dependent (e.g. strong
#' long memory), where the analytic t and the i.i.d. weights can fail.
#'
#' @param data a matrix or data frame; first column the response, the rest the
#'   predictors.
#' @param B number of bootstrap replicates.
#' @param weights multiplier-weight scheme: \code{"dependent"} (Shao 2010),
#'   \code{"rademacher"} or \code{"mammen"}.
#' @param bandwidth kernel bandwidth for \code{"dependent"} weights; defaults to
#'   \eqn{T_s^{1/3}} per scale.
#' @param dpo detrending polynomial order.
#' @param int logical. If TRUE the integration process is applied.
#' @param np number of point scales.
#' @param min_boxes minimum number of non-overlapping boxes \eqn{T_s =
#'   \lfloor N/s\rfloor} required for inference; at scales below this floor the
#'   bootstrap is skipped and the interval is returned as \code{NA} (default 15).
#'   Score-based inference treats the boxes as disjoint sampling units;
#'   \code{fracreg.WB()} therefore takes no \code{overlap} argument by design.
#' @param abs logical. If TRUE the absolute detrended covariance is used.
#' @return A tibble with the scale \code{s} and, per predictor, the estimate
#'   (\code{beta_*}), the lower/upper interval bounds (\code{lower_*},
#'   \code{upper_*}) and the p-value (\code{p_*}), one row per scale.
#' @references
#' Shao, X. (2010). The dependent wild bootstrap. \emph{Journal of the American
#' Statistical Association}, 105(489), 218-235.
#'
#' Mammen, E. (1993). Bootstrap and wild bootstrap for high dimensional linear
#' models. \emph{The Annals of Statistics}, 21(1), 255-285.
#'
#' Tilfani, O., Kristoufek, L., Ferreira, P. and El Boukfaoui, M. Y. (2022).
#' Heterogeneity in economic relationships: scale dependence through the
#' multivariate fractal regression. \emph{Physica A}, 588, 126530.
#' @seealso \code{\link{fracreg}}, \code{\link{fracreg.IUTest}},
#' \code{vignette("DFATools")}
#' @importFrom stats quantile rnorm runif
#' @examples
#' set.seed(1)
#' d <- data.frame(y = cumsum(rnorm(300)), x1 = cumsum(rnorm(300)),
#'                 x2 = cumsum(rnorm(300)))
#' \donttest{
#' fracreg.WB(d, B = 199, np = 15)
#' }
#' @export
fracreg.WB <- function(data, B = 999, weights = c("dependent", "rademacher", "mammen"),
                       bandwidth = NULL, dpo = 1, int = TRUE, np = 91,
                       min_boxes = 15, abs = FALSE) {
  weights <- match.arg(weights)
  # Score-based inference treats the boxes as disjoint sampling units;
  # fracreg.WB does not accept an `overlap` argument by design (paper M8).
  .check_common(np, dpo, int, overlap = FALSE, "fracreg.WB")
  .check_matrix(data, 2, "fracreg.WB")
  .check_B(B, "fracreg.WB")
  if (!is.numeric(min_boxes) || length(min_boxes) != 1L ||
      !is.finite(min_boxes) || min_boxes < 1)
    stop("fracreg.WB(): `min_boxes` must be a single positive number.", call. = FALSE)
  if (nrow(data) < 500L)
    warning("fracreg.WB(): N = ", nrow(data), " < 500; DFA-based inference may be ",
            "unreliable (Likens et al., 2019).", call. = FALSE)
  data <- as.matrix(data)
  nc <- ncol(data); k <- nc - 1L
  cn <- colnames(data); if (is.null(cn)) cn <- paste0("x", seq_len(nc))
  pred <- cn[2:nc]

  pb <- .fracreg_perbox(data, dpo = dpo, int = int, np = np, overlap = FALSE, abs = abs)
  s <- pb$scales
  beta <- lower <- upper <- pval <- matrix(NA_real_, np, k)
  for (ki in seq_len(np)) {
    ps <- pb$perscale[[ki]]; Ts <- ps$Ts
    beta[ki, ] <- ps$beta
    if (Ts < min_boxes) next                                  # too few boxes -> NA SE/CI/p
    W  <- .wb_weights(Ts, B, weights, bandwidth)              # Ts x B
    BS <- ps$beta + ps$Fxx_inv %*% (crossprod(ps$scores, W) / Ts)   # k x B
    for (j in seq_len(k)) {
      q <- stats::quantile(BS[j, ], c(0.025, 0.975), na.rm = TRUE)
      lower[ki, j] <- q[1]; upper[ki, j] <- q[2]
      pval[ki, j]  <- 2 * min(mean(BS[j, ] <= 0), mean(BS[j, ] >= 0))
    }
  }
  few_idx <- which(vapply(pb$perscale, function(z) z$Ts < min_boxes, logical(1L)))
  if (length(few_idx)) {
    lst <- if (length(few_idx) > 6) paste(c(utils::head(s[few_idx], 5), "..."), collapse = ", ")
           else paste(s[few_idx], collapse = ", ")
    warning("fracreg.WB(): scales {", lst, "} have fewer than min_boxes = ", min_boxes,
            " non-overlapping boxes; their standard errors and intervals were set ",
            "to NA. Reduce the maximum scale or increase N.", call. = FALSE)
  }
  out <- data.frame(s = s)
  for (j in seq_len(k)) {
    out[[paste0("beta_",  pred[j])]] <- beta[, j]
    out[[paste0("lower_", pred[j])]] <- lower[, j]
    out[[paste0("upper_", pred[j])]] <- upper[, j]
    out[[paste0("p_",     pred[j])]] <- pval[, j]
  }
  tibble::as_tibble(out)
}

#' Detrended Fluctuation Analysis
#'
#' Calculates the detrended fluctuation function of a single series.
#'
#' The C primitive computes the mean squared fluctuation \eqn{F^2(s)} (the
#' average of the within-box detrended residual variance over the boxes of
#' size \eqn{s}). The return uses the conventional Peng et al. (1994) form:
#' \itemize{
#'   \item \code{$F = sqrt(F^2)} -- the root mean-squared fluctuation, so that
#'         \eqn{\alpha} is the slope of \eqn{\log F(s)} vs \eqn{\log s};
#'   \item \code{$F2 = F^2} -- the squared fluctuation, the legacy quantity
#'         consumed internally by the package (\code{rhodcca}, \code{fracreg},
#'         \dots) and useful for combining DFA values across series;
#'   \item \code{$alpha} -- the estimated DFA exponent (= Hurst exponent for
#'         self-similar processes), the slope of \eqn{\log F^2(s)/2} against
#'         \eqn{\log s} over all positive scales.
#' }
#'
#' @param data a numeric vector or single-column matrix.
#' @param dpo detrending polynomial order (default 1).
#' @param int logical; if TRUE the input is integrated into the profile (the
#'   standard use for stationary inputs).
#' @param np number of scales (box sizes).
#' @param overlap logical; if TRUE overlapping windows are used.
#' @return A list with the scale vector \code{s}, the fluctuation function
#'   \code{F}, the squared fluctuation \code{F2} and the estimated DFA exponent
#'   \code{alpha}.
#' @references
#' Peng, C.-K., Buldyrev, S. V., Havlin, S., Simons, M., Stanley, H. E. and
#' Goldberger, A. L. (1994). Mosaic organization of DNA nucleotides.
#' \emph{Physical Review E}, 49(2), 1685-1689.
#'
#' Hu, K., Ivanov, P. Ch., Chen, Z., Carpena, P. and Stanley, H. E. (2001).
#' Effect of trends on detrended fluctuation analysis. \emph{Physical Review
#' E}, 64(1), 011114.
#'
#' Kantelhardt, J. W., Zschiegner, S. A., Koscielny-Bunde, E., Havlin, S.,
#' Bunde, A. and Stanley, H. E. (2002). Multifractal detrended fluctuation
#' analysis of nonstationary time series. \emph{Physica A}, 316, 87-114.
#' @seealso \code{\link{plotdfa}}, \code{vignette("DFATools")}
#' @useDynLib DFATools, .registration=TRUE
#' @examples
#' set.seed(1)
#' x <- cumsum(rnorm(300))           # random walk: alpha ~ 1.5
#' fy <- dfa(x, np = 20)
#' round(fy$alpha, 3)
#' @export

dfa<-function(data,dpo=1,int=TRUE,np=91,overlap=TRUE){
	.check_common(np, dpo, int, overlap, "dfa")
	.check_series(data, "dfa")
	data<-as.matrix(data)
	if (int ==TRUE){int=1} else{int=0}
	dpo<-as.numeric(dpo+1)
	mx<-round(nrow(data)/5)
seq1<-data[,1]
NPTS <- length(seq1)
NFIT	<- dpo
IFLAG <- as.numeric(int)
NR 	<- as.numeric(np)
SW 	<- as.numeric(overlap)
MINBOX	<- 10
MAXBOX	<- as.numeric(mx)
cfg <- as.integer(cbind(NPTS,NFIT,IFLAG,NR,SW,MINBOX,MAXBOX))
rsi1<- numeric(NR+1)
msei1<- numeric(NR+1)
ans1<- .C("rdfa", cfg, as.numeric(seq1), as.integer(rsi1), as.numeric(msei1),PACKAGE = "DFATools")
sn<- ans1[[3]][1:NR+1]
F2<- ans1[[4]][1:NR+1]                      # squared fluctuation from C
Fs<- sqrt(pmax(F2,0))                       # Peng F = sqrt(F^2)
ok<- is.finite(F2)&F2>0
alpha<- if(sum(ok)<4) NA_real_ else
	0.5*unname(stats::coef(stats::lm(log(F2[ok])~log(sn[ok])))[[2]])
list(s=sn, F=Fs, F2=F2, alpha=alpha)
}
#' rho Detrended Cross-Correlation Coefficient
#'
#' Calculates rho DCCA
#' @param data is a matrix of time series
#' @param dpo detrending polynomial order
#' @param int logical. if TRUE integration process will be applied.
#' @param np number of point scales.
#' @param overlap logical. if TRUE overlapping windows will be applied.
#' @return Scale s, Rho-DCCA
#' @references
#' Podobnik, B. and Stanley, H. E. (2008). Detrended cross-correlation analysis:
#' a new method for analyzing two nonstationary time series.
#' \emph{Physical Review Letters}, 100(8), 084102.
#'
#' Zebende, G. F. (2011). DCCA cross-correlation coefficient: quantifying level
#' of cross-correlation. \emph{Physica A}, 390(4), 614-618.
#'
#' Kwapien, J., Oswiecimka, P. and Drozdz, S. (2015). Detrended fluctuation
#' analysis made flexible to detect range of cross-correlated fluctuations.
#' \emph{Physical Review E}, 92(5), 052815.
#'
#' Cavalcanti, S. L. (2019). \emph{Aplicacoes de DFA, DCCA e DPCCA em focos
#' de calor na Amazonia Legal} (PhD thesis). Universidade Federal Rural de
#' Pernambuco.
#' @seealso \code{vignette("DFATools")}
#' @importFrom tibble as_tibble
#' @importFrom gdata upperTriangle
#' @useDynLib DFATools, .registration=TRUE
#' @examples
#' set.seed(1)
#' d <- data.frame(x = cumsum(rnorm(300)), y = cumsum(rnorm(300)))
#' rhodcca(d, np = 20)
#' @export
#'
rhodcca<-function(data,dpo=1,int=TRUE,np=91,overlap=TRUE){
	.check_common(np, dpo, int, overlap, "rhodcca")
	.check_matrix(data, 2, "rhodcca")
	data<-as.matrix(data)
	if (int ==TRUE){int=1} else{int=0}
	dpo<-as.numeric(dpo+1)
	nc<-ncol(data)
	pairs<-t(combn(nc,2))
	comb<-nCm(nc,2)
	sn<-matrix(,ncol=1) #Escalas
	mx<-round(nrow(data)/5)
	fn<-array(,dim=c(nc,nc,np)) #Variancia-Covarância sem tendência
	size=nrow(data)
#DFA
	for (i in 1:nc)
	{
seq1<-data[,i]
# ****** CONFIGURATION STRUCTURE
NPTS <- length(seq1) #size of series
NFIT	<- dpo		# order of the regression fit, plus 1
IFLAG <- as.numeric(int)		# integrate the input data if non-zero
NR 	<- as.numeric(np)		# number of box sizes
SW 	<- as.numeric(overlap)		# sliding window
MINBOX	<- 10		# minimum box size
MAXBOX	<- as.numeric(mx)		# maximum box size
cfg <- as.integer(cbind(NPTS,NFIT,IFLAG,NR,SW,MINBOX,MAXBOX))
rsi1<- numeric(NR+1)
msei1<- numeric(NR+1)
ans1<- .C("rdfa", cfg, as.numeric(seq1), as.integer(rsi1), as.numeric(msei1),PACKAGE = "DFATools")
sn<- ans1[[3]][1:NR+1]
fn[i,i,]<-ans1[[4]][1:NR+1]
}
#DCCA
#fc<-matrix(,ncol=comb)
for (i in 1:comb){
seq1<-data[,pairs[i,1]]
seq2<-data[,pairs[i,2]]
# ****** CONFIGURATION STRUCTURE
NPTS <- length(seq1) #size of series
NFIT	<- dpo		# order of the regression fit, plus 1
IFLAG <- as.numeric(int)		# integrate the input data if non-zero
NR 	<- as.numeric(np)	# number of box sizes
SW 	<- as.numeric(overlap)		# sliding window
MINBOX	<- 10		# minimum box size
MAXBOX	<- as.numeric(mx)	# maximum box size
ABSFLAG	<- 0		# absolute flag for dcca
cfg <- as.integer(cbind(NPTS,NFIT,IFLAG,NR,SW,MINBOX,MAXBOX,ABSFLAG))
rs12 <- numeric(NR+1)
mse12 <- numeric(NR+1)
ans12 <- .C("rdcca", cfg, as.numeric(seq1), as.numeric(seq2), as.integer(rs12), as.numeric(mse12),PACKAGE = "DFATools")
fn[pairs[i,1],pairs[i,2],] <- ans12[[5]][1:NR+1]
fn[pairs[i,2],pairs[i,1],] <- ans12[[5]][1:NR+1]
}
#DPCCA
pn<-array(,dim=c(nc,nc,np)) #pDCCA
for(i in 1:(np-1)){
	pn[,,i]<-cov2cor(fn[,,i])
}
ifelse(nc!=2,rhodcca<-tibble::as_tibble(cbind.data.frame(sn,t(apply(pn, 3L, gdata::upperTriangle)))),
       rhodcca<-tibble::as_tibble(cbind.data.frame(sn,apply(pn, 3L, gdata::upperTriangle))))
names(rhodcca)<-c("s",paste("DCCA",apply(t(combn(nc,2)),1,paste0,collapse = ""),sep=""))
return(rhodcca)
}

#' rho Detrended Partial Cross-Correlation Coefficient
#'
#' Calculates DPCCA
#' @param data is a matrix of time series
#' @param dpo detrending polynomial order
#' @param int logical. if TRUE integration process will be applied.
#' @param np number of point scales.
#' @param overlap logical. if TRUE overlapping windows will be applied.
#' @return Scale s, Rho DPCCA
#' @references
#' Yuan, N. et al. (2015). Detrended partial-cross-correlation analysis: a new
#' method for analyzing correlations in complex system. \emph{Scientific
#' Reports}, 5, 8143.
#' @seealso \code{\link{rhodcca}}, \code{vignette("DFATools")}
#' @useDynLib DFATools, .registration=TRUE
#' @examples
#' set.seed(1)
#' d <- data.frame(x = cumsum(rnorm(300)), y = cumsum(rnorm(300)),
#'                 z = cumsum(rnorm(300)))
#' rhodpcca(d, np = 20)
#' @export
#'
rhodpcca<-function(data,dpo=1,int=TRUE,np=91,overlap=TRUE){
	.check_common(np, dpo, int, overlap, "rhodpcca")
	.check_matrix(data, 3, "rhodpcca")
	data<-as.matrix(data)
	if (int ==TRUE){int=1} else{int=0}
	dpo<-as.numeric(dpo+1)
	nc<-ncol(data)
	pairs<-t(combn(nc,2))
	comb<-nCm(nc,2)
	sn<-matrix(,ncol=1) #Escalas
	mx<-round(nrow(data)/5)
	fn<-array(,dim=c(nc,nc,np)) #Variancia-Covarância sem tendência
	size=nrow(data)
#DFA
	for (i in 1:nc)
	{
seq1<-data[,i]
# ****** CONFIGURATION STRUCTURE
NPTS <- length(seq1) #size of series
NFIT	<- dpo		# order of the regression fit, plus 1
IFLAG <- as.numeric(int)		# integrate the input data if non-zero
NR 	<- as.numeric(np)		# number of box sizes
SW 	<- as.numeric(overlap)		# sliding window
MINBOX	<- 10		# minimum box size
MAXBOX	<- as.numeric(mx)		# maximum box size
cfg <- as.integer(cbind(NPTS,NFIT,IFLAG,NR,SW,MINBOX,MAXBOX))
rsi1<- numeric(NR+1)
msei1<- numeric(NR+1)
ans1<- .C("rdfa", cfg, as.numeric(seq1), as.integer(rsi1), as.numeric(msei1),PACKAGE = "DFATools")
sn<- ans1[[3]][1:NR+1]
fn[i,i,]<-ans1[[4]][1:NR+1]
}
#DCCA
#fc<-matrix(,ncol=comb)
for (i in 1:comb){
seq1<-data[,pairs[i,1]]
seq2<-data[,pairs[i,2]]
# ****** CONFIGURATION STRUCTURE
NPTS <- length(seq1) #size of series
NFIT	<- dpo		# order of the regression fit, plus 1
IFLAG <- as.numeric(int)		# integrate the input data if non-zero
NR 	<- as.numeric(np)	# number of box sizes
SW 	<- as.numeric(overlap)		# sliding window
MINBOX	<- 10		# minimum box size
MAXBOX	<- as.numeric(mx)	# maximum box size
ABSFLAG	<- 0		# absolute flag for dcca
cfg <- as.integer(cbind(NPTS,NFIT,IFLAG,NR,SW,MINBOX,MAXBOX,ABSFLAG))
rs12 <- numeric(NR+1)
mse12 <- numeric(NR+1)
ans12 <- .C("rdcca", cfg, as.numeric(seq1), as.numeric(seq2), as.integer(rs12), as.numeric(mse12),PACKAGE = "DFATools")
fn[pairs[i,1],pairs[i,2],] <- ans12[[5]][1:NR+1]
fn[pairs[i,2],pairs[i,1],] <- ans12[[5]][1:NR+1]
}
#DPCCA
pn<-array(,dim=c(nc,nc,np)) #pDCCA
dp<-array(,dim=c(nc,nc,np)) #pDPCCA
for(i in 1:(np-1)){
	pn[,,i]<-cov2cor(fn[,,i])
	dp[,,i]<-cor2pcor(pn[,,i])
}
  rdpcca<-tibble::as_tibble(cbind.data.frame(sn,t(apply(dp, 3L, gdata::upperTriangle))))
  names(rdpcca)<-c("s",paste("DPCCA",apply(t(combn(nc,2)),1,paste0,collapse = ""),sep=""))
  return(rdpcca)
}

#' rho Detrended Multiple-Correlation Coefficient
#'
#' Calculates DMC²
#' @param data is a matrix of time series
#' @param dpo detrending polynomial order
#' @param int logical. if TRUE integration process will be applied.
#' @param np number of point scales.
#' @param overlap logical. if TRUE overlapping windows will be applied.
#' @return Scale s, Multiple Detrended Correlation
#' @references
#' Zebende, G. F. and da Silva Filho, A. M. (2018). Detrended multiple
#' cross-correlation coefficient. \emph{Physica A}, 510, 91-97.
#'
#' Wang, F., Xu, J. and Fan, Q. (2021). Statistical test for detrended multiple
#' cross-correlation coefficient. \emph{Communications in Nonlinear Science and
#' Numerical Simulation}, 99, 105781.
#' @seealso \code{vignette("DFATools")}
#' @useDynLib DFATools, .registration=TRUE
#' @examples
#' set.seed(1)
#' d <- data.frame(y = cumsum(rnorm(300)), x1 = cumsum(rnorm(300)),
#'                 x2 = cumsum(rnorm(300)))
#' dmc2(d, np = 20)
#' @export

dmc2<-function(data,dpo=1,int=TRUE,np=91,overlap=TRUE){
	.check_common(np, dpo, int, overlap, "dmc2")
	.check_matrix(data, 2, "dmc2")
	data<-as.matrix(data)
	if (int ==TRUE){int=1} else{int=0}
	dpo<-as.numeric(dpo+1)
	nc<-ncol(data)
	pairs<-t(combn(nc,2))
	comb<-nCm(nc,2)
	sn<-matrix(,ncol=1) #Escalas
	mx<-round(nrow(data)/5)
	fn<-array(,dim=c(nc,nc,np)) #Variancia-Covarância sem tendência
	dmc2<-matrix(NA,nrow=np) #DMC2 Zebende et al (2018)
	size=nrow(data)
#DFA
	for (i in 1:nc)
	{
seq1<-data[,i]
# ****** CONFIGURATION STRUCTURE
NPTS <- length(seq1) #size of series
NFIT	<- dpo		# order of the regression fit, plus 1
IFLAG <- as.numeric(int)		# integrate the input data if non-zero
NR 	<- as.numeric(np)		# number of box sizes
SW 	<- as.numeric(overlap)		# sliding window
MINBOX	<- 10		# minimum box size
MAXBOX	<- as.numeric(mx)		# maximum box size
cfg <- as.integer(cbind(NPTS,NFIT,IFLAG,NR,SW,MINBOX,MAXBOX))
rsi1<- numeric(NR+1)
msei1<- numeric(NR+1)
ans1<- .C("rdfa", cfg, as.numeric(seq1), as.integer(rsi1), as.numeric(msei1),PACKAGE = "DFATools")
sn<- ans1[[3]][1:NR+1]
fn[i,i,]<-ans1[[4]][1:NR+1]
}
#DCCA
#fc<-matrix(,ncol=comb)
for (i in 1:comb){
seq1<-data[,pairs[i,1]]
seq2<-data[,pairs[i,2]]
# ****** CONFIGURATION STRUCTURE
NPTS <- length(seq1) #size of series
NFIT	<- dpo		# order of the regression fit, plus 1
IFLAG <- as.numeric(int)		# integrate the input data if non-zero
NR 	<- as.numeric(np)	# number of box sizes
SW 	<- as.numeric(overlap)		# sliding window
MINBOX	<- 10		# minimum box size
MAXBOX	<- as.numeric(mx)	# maximum box size
ABSFLAG	<- 0		# absolute flag for dcca
cfg <- as.integer(cbind(NPTS,NFIT,IFLAG,NR,SW,MINBOX,MAXBOX,ABSFLAG))
rs12 <- numeric(NR+1)
mse12 <- numeric(NR+1)
ans12 <- .C("rdcca", cfg, as.numeric(seq1), as.numeric(seq2), as.integer(rs12), as.numeric(mse12),PACKAGE = "DFATools")
fn[pairs[i,1],pairs[i,2],] <- ans12[[5]][1:NR+1]
fn[pairs[i,2],pairs[i,1],] <- ans12[[5]][1:NR+1]
}
#Rho DCCA
pn<-array(,dim=c(nc,nc,np)) #pDCCA
for(i in 1:(np-1)){
	pn[,,i]<-cov2cor(fn[,,i])
	}
for(k in 1:np){
dmc2[k]<-pn[1,(2:nc),k]%*%solve(pn[(2:nc),(2:nc),k])%*%pn[(2:nc),1,k]
}
DMC2<-tibble::as_tibble(cbind.data.frame(sn,dmc2))
return(DMC2)
}

#' f² scale-wise effect sizes
#'
#' Calculates f² scale-wise effect sizes
#' @param data is a matrix of time series
#' @param dpo detrending polynomial order
#' @param int logical. if TRUE integration process will be applied.
#' @param np number of point scales.
#' @param overlap logical. if TRUE overlapping windows will be applied.
#' @return Scale s, f² scale-wise effect sizes
#' @references
#' Barreto, I. D. C., Dore, L. H., Stosic, T. and Stosic, B. D. (2021).
#' Extending DFA-based multiple linear regression inference: application to
#' acoustic impedance models. \emph{Physica A}, 582, 126259.
#'
#' Cohen, J. (1988). \emph{Statistical Power Analysis for the Behavioral
#' Sciences}, 2nd ed. Lawrence Erlbaum Associates.
#' @seealso \code{\link{fracreg}}, \code{vignette("DFATools")}
#' @useDynLib DFATools, .registration=TRUE
#' @importFrom dplyr inner_join
#' @examples
#' set.seed(1)
#' d <- data.frame(y = cumsum(rnorm(300)), x1 = cumsum(rnorm(300)),
#'                 x2 = cumsum(rnorm(300)))
#' effsizeDFA(d, np = 20)
#' @export

effsizeDFA<-function(data,dpo=1,int=TRUE,np=91,overlap=TRUE){
  .check_common(np, dpo, int, overlap, "effsizeDFA")
  .check_matrix(data, 3, "effsizeDFA")
  if(ncol(data)<3){warning("the number of regressors are lower than 2\nplease provide at least 2 independent variables")}
else{
	nC<-ncol(data)
	names<-names(data)[2:nC]
	r2<-matrix(,nrow=np,ncol=nC)
	sn<-dmc2(data,dpo=dpo,int=int,np=np,overlap=overlap)$sn
	r2[,1]<-dmc2(data,dpo=dpo,int=int,np=np,overlap=overlap)$dmc2
	for(i in c(2:nC)){
	data2<-data[,-i]
	r2[,i]<-dmc2(data2,dpo=dpo,int=int,np=np,overlap=overlap)$dmc2
	data2<-NULL
	}
	r2<-r2^2
	f2<-matrix(,nrow=np,ncol=nC-1)
	for (i in c(2:nC)){
	f2[,(i-1)]<-(r2[,1]-r2[,i])/(1-r2[,1])
	}
	names1<-c("s","R2",paste0("R2_",names))
	names2<-c("s",paste0("f2_",names))
	r2<-tibble::as_tibble(cbind.data.frame(sn,r2), column_name = names1)
	f2<-tibble::as_tibble(cbind.data.frame(sn,f2), column_name = names2)
	names(r2)<-names1
	names(f2)<-names2
	rf<- r2 %>% dplyr::inner_join(f2,by = "s")
	return(rf)
	}
  }

#' Beta DFA
#'
#' Calculates Beta DFA
#' @param data is a matrix of time series
#' @param dpo detrending polynomial order
#' @param int logical. if TRUE integration process will be applied.
#' @param np number of point scales.
#' @param overlap logical. if TRUE overlapping windows will be applied.
#' @return Scale s, Beta DFA estimates
#' @keywords internal
#' @useDynLib DFATools, .registration=TRUE

bdfa<-function(data,dpo=1,int=TRUE,np=91,overlap=TRUE){
	data<-as.matrix(data)
	if (int ==TRUE){int=1} else{int=0}
	dpo<-as.numeric(dpo+1)
	nc<-ncol(data)
	pairs<-t(combn(nc,2))
	comb<-nCm(nc,2)
	sn<-matrix(,ncol=1) #Escalas
	mx<-round(nrow(data)/5)
	fn<-array(,dim=c(nc,nc,np)) #Variancia-Covarância sem tendência
	bn<-array(,dim=c((nc-1),1,np)) #Parâmetros Bdfa(s)
	size=nrow(data)
#DFA
	for (i in 1:nc)
	{
seq1<-data[,i]
# ****** CONFIGURATION STRUCTURE
NPTS <- length(seq1) #size of series
NFIT	<- dpo		# order of the regression fit, plus 1
IFLAG <- as.numeric(int)		# integrate the input data if non-zero
NR 	<- as.numeric(np)		# number of box sizes
SW 	<- as.numeric(overlap)		# sliding window
MINBOX	<- 10		# minimum box size
MAXBOX	<- as.numeric(mx)		# maximum box size
cfg <- as.integer(cbind(NPTS,NFIT,IFLAG,NR,SW,MINBOX,MAXBOX))
rsi1<- numeric(NR+1)
msei1<- numeric(NR+1)
ans1<- .C("rdfa",cfg,as.numeric(seq1),as.integer(rsi1),as.numeric(msei1),PACKAGE = "DFATools")
sn<- ans1[[3]][1:NR+1]
fn[i,i,]<-ans1[[4]][1:NR+1]
}
#DCCA
#fc<-matrix(,ncol=comb)
for (i in 1:comb){
seq1<-data[,pairs[i,1]]
seq2<-data[,pairs[i,2]]
# ****** CONFIGURATION STRUCTURE
NPTS <- length(seq1) #size of series
NFIT	<- dpo		# order of the regression fit, plus 1
IFLAG <- as.numeric(int)		# integrate the input data if non-zero
NR 	<- as.numeric(np)	# number of box sizes
SW 	<- as.numeric(overlap)		# sliding window
MINBOX	<- 10		# minimum box size
MAXBOX	<- as.numeric(mx)	# maximum box size
ABSFLAG	<- 0		# absolute flag for dcca
cfg <- as.integer(cbind(NPTS,NFIT,IFLAG,NR,SW,MINBOX,MAXBOX,ABSFLAG))
rs12 <- numeric(NR+1)
mse12 <- numeric(NR+1)
ans12 <- .C("rdcca", cfg, as.numeric(seq1), as.numeric(seq2), as.integer(rs12), as.numeric(mse12),PACKAGE = "DFATools")
fn[pairs[i,1],pairs[i,2],] <- ans12[[5]][1:NR+1]
fn[pairs[i,2],pairs[i,1],] <- ans12[[5]][1:NR+1]
}
###Fractal Regression
for(k in 1:np){
		bn[,1,k]<-solve(fn[(2:nc),(2:nc),k])%*%fn[(2:nc),1,k]
			}
bdfa<-list(sn,bn)
names(bdfa)<-c("s","BDFA")
return(bdfa)
}

#' Beta DFA
#'
#' Calculates Beta DFA
#' @param data is a matrix of time series
#' @param dpo detrending polynomial order
#' @param int logical. if TRUE integration process will be applied.
#' @param np number of point scales.
#' @param overlap logical. if TRUE overlapping windows will be applied.
#' @return Scale s, Beta DFA estimates
#' @references
#' Kristoufek, L. (2015). Detrended fluctuation analysis as a regression
#' framework: estimating dependence at different scales. \emph{Physical Review
#' E}, 91(2), 022802.
#' @seealso \code{\link{sbdfa}}, \code{\link{fracreg}},
#' \code{vignette("DFATools")}
#' @useDynLib DFATools, .registration=TRUE
#' @examples
#' set.seed(1)
#' d <- data.frame(y = cumsum(rnorm(300)), x1 = cumsum(rnorm(300)),
#'                 x2 = cumsum(rnorm(300)))
#' betadfa(d, np = 20)
#' @export

betadfa<-function(data,dpo=1,int=TRUE,np=91,overlap=TRUE){
  .check_common(np, dpo, int, overlap, "betadfa")
  .check_matrix(data, 2, "betadfa")
  nc<-ncol(data)
  cn<-c("s",names(data)[2:nc])
  data<-as.matrix(data)
  if (int ==TRUE){int=1} else{int=0}
  dpo<-as.numeric(dpo+1)
  pairs<-t(combn(nc,2))
  comb<-nCm(nc,2)
  sn<-matrix(,ncol=1) #Escalas
  mx<-round(nrow(data)/5)
  fn<-array(,dim=c(nc,nc,np)) #Variancia-Covarância sem tendência
  bn<-array(,dim=c((nc-1),1,np)) #Parâmetros Bdfa(s)
  size=nrow(data)
  #DFA
  for (i in 1:nc)
  {
    seq1<-data[,i]
    # ****** CONFIGURATION STRUCTURE
    NPTS <- length(seq1) #size of series
    NFIT	<- dpo		# order of the regression fit, plus 1
    IFLAG <- as.numeric(int)		# integrate the input data if non-zero
    NR 	<- as.numeric(np)		# number of box sizes
    SW 	<- as.numeric(overlap)		# sliding window
    MINBOX	<- 10		# minimum box size
    MAXBOX	<- as.numeric(mx)		# maximum box size
    cfg <- as.integer(cbind(NPTS,NFIT,IFLAG,NR,SW,MINBOX,MAXBOX))
    rsi1<- numeric(NR+1)
    msei1<- numeric(NR+1)
    ans1<- .C("rdfa",cfg,as.numeric(seq1),as.integer(rsi1),as.numeric(msei1),PACKAGE = "DFATools")
    sn<- ans1[[3]][1:NR+1]
    fn[i,i,]<-ans1[[4]][1:NR+1]
  }
  #DCCA
  #fc<-matrix(,ncol=comb)
  for (i in 1:comb){
    seq1<-data[,pairs[i,1]]
    seq2<-data[,pairs[i,2]]
    # ****** CONFIGURATION STRUCTURE
    NPTS <- length(seq1) #size of series
    NFIT	<- dpo		# order of the regression fit, plus 1
    IFLAG <- as.numeric(int)		# integrate the input data if non-zero
    NR 	<- as.numeric(np)	# number of box sizes
    SW 	<- as.numeric(overlap)		# sliding window
    MINBOX	<- 10		# minimum box size
    MAXBOX	<- as.numeric(mx)	# maximum box size
    ABSFLAG	<- 0		# absolute flag for dcca
    cfg <- as.integer(cbind(NPTS,NFIT,IFLAG,NR,SW,MINBOX,MAXBOX,ABSFLAG))
    rs12 <- numeric(NR+1)
    mse12 <- numeric(NR+1)
    ans12 <- .C("rdcca", cfg, as.numeric(seq1), as.numeric(seq2), as.integer(rs12), as.numeric(mse12),PACKAGE = "DFATools")
    fn[pairs[i,1],pairs[i,2],] <- ans12[[5]][1:NR+1]
    fn[pairs[i,2],pairs[i,1],] <- ans12[[5]][1:NR+1]
  }
  ###Fractal Regression
  for(k in 1:np){
    bn[,1,k]<-solve(fn[(2:nc),(2:nc),k])%*%fn[(2:nc),1,k]
  }
  bdfa<-tibble::as_tibble(cbind.data.frame(sn,t(as.matrix(bn[,,]))),column_name= cn)
  names(bdfa)<-cn
  return(bdfa)
}
#' Standardized Beta DFA
#'
#' Calculates Standardized Beta DFA
#' @param data is a matrix of time series
#' @param dpo detrending polynomial order
#' @param int logical. if TRUE integration process will be applied.
#' @param np number of point scales.
#' @param overlap logical. if TRUE overlapping windows will be applied.
#' @return Scale s, Standardized Beta DFA estimates
#' @references
#' Barreto, I. D. C., Dore, L. H., Stosic, T. and Stosic, B. D. (2021).
#' Extending DFA-based multiple linear regression inference: application to
#' acoustic impedance models. \emph{Physica A}, 582, 126259.
#'
#' Kristoufek, L. (2015). Detrended fluctuation analysis as a regression
#' framework: estimating dependence at different scales. \emph{Physical Review
#' E}, 91(2), 022802.
#' @seealso \code{\link{betadfa}}, \code{vignette("DFATools")}
#' @useDynLib DFATools, .registration=TRUE
#' @examples
#' set.seed(1)
#' d <- data.frame(y = cumsum(rnorm(300)), x1 = cumsum(rnorm(300)),
#'                 x2 = cumsum(rnorm(300)))
#' sbdfa(d, np = 20)
#' @export

sbdfa<-function(data,dpo=1,int=TRUE,np=91,overlap=TRUE){
  .check_common(np, dpo, int, overlap, "sbdfa")
  .check_matrix(data, 2, "sbdfa")
  nc<-ncol(data)
  cn<-c("s",names(data)[2:nc])
  data<-as.matrix(data)
  if (int ==TRUE){int=1} else{int=0}
  dpo<-as.numeric(dpo+1)
  pairs<-t(combn(nc,2))
  comb<-nCm(nc,2)
  sn<-matrix(,ncol=1) #Escalas
  bs<-array(,dim=c(1,(nc-1),np)) #Beta Padronizado
  mx<-round(nrow(data)/5)
  fn<-array(,dim=c(nc,nc,np)) #Variancia-Covarância sem tendência
  bn<-array(,dim=c((nc-1),1,np)) #Parâmetros Bdfa(s)
  size=nrow(data)
  #DFA
  for (i in 1:nc)
  {
    seq1<-data[,i]
    # ****** CONFIGURATION STRUCTURE
    NPTS <- length(seq1) #size of series
    NFIT	<- dpo		# order of the regression fit, plus 1
    IFLAG <- as.numeric(int)		# integrate the input data if non-zero
    NR 	<- as.numeric(np)		# number of box sizes
    SW 	<- as.numeric(overlap)		# sliding window
    MINBOX	<- 10		# minimum box size
    MAXBOX	<- as.numeric(mx)		# maximum box size
    cfg <- as.integer(cbind(NPTS,NFIT,IFLAG,NR,SW,MINBOX,MAXBOX))
    rsi1<- numeric(NR+1)
    msei1<- numeric(NR+1)
    ans1<- .C("rdfa",cfg,as.numeric(seq1),as.integer(rsi1),as.numeric(msei1),PACKAGE = "DFATools")
    sn<- ans1[[3]][1:NR+1]
    fn[i,i,]<-ans1[[4]][1:NR+1]
  }
  #DCCA
  #fc<-matrix(,ncol=comb)
  for (i in 1:comb){
    seq1<-data[,pairs[i,1]]
    seq2<-data[,pairs[i,2]]
    # ****** CONFIGURATION STRUCTURE
    NPTS <- length(seq1) #size of series
    NFIT	<- dpo		# order of the regression fit, plus 1
    IFLAG <- as.numeric(int)		# integrate the input data if non-zero
    NR 	<- as.numeric(np)	# number of box sizes
    SW 	<- as.numeric(overlap)		# sliding window
    MINBOX	<- 10		# minimum box size
    MAXBOX	<- as.numeric(mx)	# maximum box size
    ABSFLAG	<- 0		# absolute flag for dcca
    cfg <- as.integer(cbind(NPTS,NFIT,IFLAG,NR,SW,MINBOX,MAXBOX,ABSFLAG))
    rs12 <- numeric(NR+1)
    mse12 <- numeric(NR+1)
    ans12 <- .C("rdcca", cfg, as.numeric(seq1), as.numeric(seq2), as.integer(rs12), as.numeric(mse12),PACKAGE = "DFATools")
    fn[pairs[i,1],pairs[i,2],] <- ans12[[5]][1:NR+1]
    fn[pairs[i,2],pairs[i,1],] <- ans12[[5]][1:NR+1]
  }

  #DPCCA
  pn<-array(,dim=c(nc,nc,np)) #pDCCA
  for(i in 1:(np-1)){
    pn[,,i]<-cov2cor(fn[,,i])
  }
  ###Fractal Regression
  for(k in 1:np){
    bs[,,k]<-solve(pn[(2:nc),(2:nc),k])%*%pn[(2:nc),1,k]
  }
  bsdfa<-tibble::as_tibble(cbind.data.frame(sn,t(as.matrix(bs[,,]))),column_names=cn)
  names(bsdfa)<-cn
  return(bsdfa)
}

#' Kristoufek Test
#'
#' Calculates Kristoufek Test for Beta-DFA = Beta-OLS
#' @param data is a matrix of time series
#' @param B number of surrogate series
#' @param dpo detrending polynomial order
#' @param int logical. if TRUE integration process will be applied.
#' @param np number of point scales.
#' @param overlap logical. if TRUE overlapping windows will be applied.
#' @return A matrix with scale-wise Beta-DFA and critic region of Kristoufek Test.
#' @references
#' Kristoufek, L. (2015). Detrended fluctuation analysis as a regression
#' framework: estimating dependence at different scales. \emph{Physical Review
#' E}, 91(2), 022802.
#' @seealso \code{\link{fracreg.PStest}}, \code{\link{fracreg.IUTest}},
#' \code{vignette("DFATools")}
#' @useDynLib DFATools, .registration=TRUE
#' @examples
#' set.seed(1)
#' d <- data.frame(y = cumsum(rnorm(250)), x = cumsum(rnorm(250)))
#' \donttest{
#' fracreg.Ktest(d, B = 20, np = 15)
#' }
#' @export

fracreg.Ktest<-function(data,B=100,dpo=1,int=TRUE,np=91,overlap=TRUE){
.check_common(np, dpo, int, overlap, "fracreg.Ktest")
.check_matrix(data, 2, "fracreg.Ktest")
.check_B(B, "fracreg.Ktest")
x<-ncol(data)-1
if (int ==TRUE){int=1} else{int=0}
yx<-ncol(data)
pbdfa<-array(,dim=c(x,1,np,B))
data<-as.matrix(data)
res<-lm(data[,1]~data[,2:yx])$res
coef<-lm(data[,1]~data[,2:yx])$coeff
usur<-surrogate(res,ns=B,fft=TRUE,amplitude=TRUE)
ysur<-matrix(,nrow=nrow(data),ncol=B)
if(yx==2){
	for(i in (1:B)){
		ysur[,i]<-coef[1]+data[,2:yx]*as.matrix(coef)[2:yx,]+usur[,i]
	}}else{
for(i in (1:B)){
	ysur[,i]<-coef[1]+data[,2:yx]%*%as.matrix(coef)[2:yx,]+usur[,i]
	}}
for(h in (1:B)){
	data2<-cbind(ysur[,h],data[,2:yx])
	pbdfa[,,,h]<-bdfa(data2,dpo=dpo,int=int,np=np,overlap=overlap)$BDFA
	data2<-NULL
}
ucib<-array(,dim=c(x,1,np))
lcib<-array(,dim=c(x,1,np))
for(i in 1:x){
	for(j in 1:np){
		lcib[i,1,j]<-quantile(pbdfa[i,1,j,],probs = c(0.025))
		ucib[i,1,j]<-quantile(pbdfa[i,1,j,],probs = c(0.975))
}}
Beta<-bdfa(data,dpo=dpo,int=int,np=np,overlap=overlap)$BDFA
lcimat<-matrix(,nrow=np,ncol=x)
ucimat<-matrix(,nrow=np,ncol=x)
for(i in 1:x){
	lcimat[,i]<-lcib[i,,]
	ucimat[,i]<-ucib[i,,]
}
if(x==1){beta<-as.matrix(bdfa(data,dpo=dpo,int=int,np=np,overlap=overlap)$BDFA[,1,])
	colnames(beta)<-paste("bet",1:x,sep="")
	}else{
	beta<-t(bdfa(data,dpo=dpo,int=int,np=np)$BDFA[,1,])
	colnames(beta)<-paste("bet",1:x,sep="")
	}
colnames(lcimat)<-paste("klci",1:x,sep="")
colnames(ucimat)<-paste("kuci",1:x,sep="")
s<-bdfa(data,dpo=dpo,int=int,np=np,overlap=overlap)$s
cib<-tibble::as_tibble(cbind.data.frame(beta,lcimat,ucimat,s))
return(cib)
}

#' Podobnik-Shen Test
#'
#' Calculates Podobnik-Shen Test for Beta-DFA = 0
#' @param data is a matrix of time series
#' @param B number of surrogate series
#' @param dpo detrending polynomial order
#' @param int logical. if TRUE integration process will be applied.
#' @param np number of point scales.
#' @param overlap logical. if TRUE overlapping windows will be applied.
#' @return A matrix with scale-wise Beta-DFA and critic region of Podobnik-Shen Test.
#' @references
#' Podobnik, B., Jiang, Z.-Q., Zhou, W.-X. and Stanley, H. E. (2011).
#' Statistical tests for power-law cross-correlated processes. \emph{Physical
#' Review E}, 84(6), 066118.
#'
#' Shen, C. (2015). A new detrended semipartial cross-correlation analysis.
#' \emph{Physics Letters A}, 379(44), 2962-2969.
#' @seealso \code{\link{fracreg.Ktest}}, \code{\link{fracreg.IUTest}},
#' \code{vignette("DFATools")}
#' @useDynLib DFATools, .registration=TRUE
#' @examples
#' set.seed(1)
#' d <- data.frame(y = cumsum(rnorm(250)), x = cumsum(rnorm(250)))
#' \donttest{
#' fracreg.PStest(d, B = 20, np = 15)
#' }
#' @export

fracreg.PStest<-function(data,B=100,dpo=1,int=TRUE,np=91,overlap=TRUE){
  .check_common(np, dpo, int, overlap, "fracreg.PStest")
  .check_matrix(data, 2, "fracreg.PStest")
  .check_B(B, "fracreg.PStest")
  a<-bdfa(data,dpo=dpo,int=int,np=np,overlap=overlap)
	n<-nrow(data)
	m<-ncol(data)-1
	test<-array(0,dim=c(m,1,np,B))
	usur<-surrogate(data[,1],ns=B,fft=TRUE,amplitude=TRUE)
	for(i in (1:B)){
		data1<-cbind(usur[,i],data[,-1])
		test[,,,i]<-bdfa(data1,dpo=dpo,int=int,np=np,overlap=overlap)$BDFA
		data1<-NULL
	}
	q025<-matrix(0,ncol=m,nrow=np)
	q975<-matrix(0,ncol=m,nrow=np)
	for(k in 1:np){
		for (j in (1:m)){
			q025[k,j]<-quantile(test[j,,k,],probs=0.025)
			q975[k,j]<-quantile(test[j,,k,],probs=0.975)
					}}
if(m==1){beta<-as.matrix(bdfa(data,dpo=dpo,int=int,np=np,overlap=overlap)$BDFA[,1,])
	colnames(beta)<-paste("bet",1:m,sep="")
	}else{
	beta<-t(bdfa(data,dpo=dpo,int=int,np=np,overlap=overlap)$BDFA[,1,])
	colnames(beta)<-paste("bet",1:m,sep="")
	}
colnames(q025)<-paste("slci",1:m,sep="")
colnames(q975)<-paste("suci",1:m,sep="")
s<-bdfa(data,dpo=dpo,int=int,np=np,overlap=overlap)$s
cib<-tibble::as_tibble(cbind.data.frame(beta,q025,q975,s))
return(cib)
}

#' Intersection-Union Test
#'
#' Calculates Intersection-Union Test for Beta-DFA = 0 or Beta-DFA = Beta-OLS
#' @param data is a matrix of time series
#' @param B number of surrogate series
#' @param dpo detrending polynomial order
#' @param int logical. if TRUE integration process will be applied.
#' @param np number of point scales.
#' @param overlap logical. if TRUE overlapping windows will be applied.
#' @return A matrix with scale-wise Beta-DFA and critic region of Kristoufek Test and Podobnik-Shen Test.
#' @references
#' Barreto, I. D. C., Dore, L. H., Stosic, T. and Stosic, B. D. (2021).
#' Extending DFA-based multiple linear regression inference: application to
#' acoustic impedance models. \emph{Physica A}, 582, 126259.
#' @seealso \code{\link{fracreg.PStest}}, \code{\link{fracreg.Ktest}},
#' \code{vignette("DFATools")}
#' @useDynLib DFATools, .registration=TRUE
#' @examples
#' set.seed(1)
#' d <- data.frame(y = cumsum(rnorm(250)), x = cumsum(rnorm(250)))
#' \donttest{
#' fracreg.IUTest(d, B = 20, np = 15)
#' }
#' @export

fracreg.IUTest<-function(data,B=100,dpo=1,int=TRUE,np=91,overlap=TRUE){
	.check_common(np, dpo, int, overlap, "fracreg.IUTest")
	.check_matrix(data, 2, "fracreg.IUTest")
	.check_B(B, "fracreg.IUTest")
	t1<-fracreg.PStest(data,B=B,dpo=dpo,int=int,np=np,overlap=overlap)
	t2<-fracreg.Ktest(data,B=B,dpo=dpo,int=int,np=np,overlap=overlap)
	IUTest<-list(t1,t2)
	names(IUTest)<-c("PStest","Ktest")
	return(IUTest)
}

#' Plot DFA
#'
#' Plot of Detrended Fluctuation Analysis
#' @param dfa is a dfa object
#' @param seg logical. If TRUE, alpha DFA will be calculated in 2 segments
#' @param point indicate in which point segmented alpha DFA should be calculated
#' @param main plot title
#' @return a plot of Detrended Fluctuation Analysis.
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_jitter
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 annotate
#' @importFrom ggplot2 geom_smooth
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 ggtitle
#' @importFrom stats median
#' @examples
#' set.seed(1)
#' x <- cumsum(rnorm(300))
#' plotdfa(dfa(x, np = 20))
#' @export

plotdfa <- function(dfa,seg=FALSE,point=NULL,main=NULL) {
  if(seg){
    ind<-c(rep("A",point),rep("B",length(dfa$s)-point))
    s<-log10(dfa$s)
    Fl<-log10(dfa$F)
    df<-cbind.data.frame(s,Fl,ind)
    alfa1<-round(lm(df$Fl~df$s,subset = df$ind=="A")$coefficients[[2]],3)
    alfa2<-round(lm(df$Fl~df$s,subset = df$ind=="B")$coefficients[[2]],3)
    p1<-ggplot2::ggplot(df,ggplot2::aes(x=s,y=Fl,fill=ind))+
      ggplot2::geom_jitter()+
      ggplot2::xlab(expression(log[10](s)))+
      ggplot2::ylab(expression(log[10](F[X]~(s))))+
      ggplot2::scale_x_continuous(breaks=seq(min(s),max(s),by=0.25))+
      ggplot2::theme_bw()+
      ggplot2::annotate("text",label=deparse(substitute(paste(alpha[1],"=",a),list(a=alfa1))),parse=TRUE,x = mean(s[1:point])-0.2*mean(s),y=median(Fl[1:point]))+
      ggplot2::annotate("text",label=deparse(substitute(paste(alpha[2],"=",b),list(b=alfa2))),parse=TRUE,x = mean(s[point:length(s)]),y=median(Fl[point:length(s)])-0.1)+
      ggplot2::geom_smooth(method = "lm",se=FALSE)+ggplot2::theme(legend.position = "none")+
      ggplot2::ggtitle(main)
    p1
  } else {
    s<-log10(dfa$s)
    Fl<-log10(dfa$F)
    df<-cbind.data.frame(s,Fl)
    alfa<-round(lm(df$Fl~df$s)$coefficients[[2]],3)
    p1<-ggplot2::ggplot(df,ggplot2::aes(x=s,y=Fl))+
      ggplot2::geom_jitter()+
      ggplot2::xlab(expression(log[10](s)))+
      ggplot2::ylab(expression(log[10](F[X]~(s))))+
      ggplot2::scale_x_continuous(breaks=seq(min(s),max(s),by=0.25))+
      ggplot2::theme_bw()+ggplot2::annotate("text",label=deparse(substitute(paste(alpha,"=",a),list(a=alfa))),parse=TRUE,x = unname(quantile(s,0.25)),y=median(Fl))+
      ggplot2::geom_smooth(method = "lm",se=FALSE,show.legend = FALSE)+
      ggplot2::ggtitle(main)
    p1
  }
}

#' Plot DCCA
#'
#' Plot of Detrended Cross Correlation Analysis
#' @param dcca is a fracreg object
#' @param seg logical. If TRUE, alpha DCCA will be calculated in 2 segments
#' @param point indicate in which point segmented alpha DFA should be calculated
#' @param main plot title
#' @return a plot of Detrended Cross Correlation Analysis.
#' @examples
#' # 'dcca' is a list with the box scales and the detrended covariance F^2_XY(s)
#' dcca <- list(s = c(10, 20, 40, 80, 160), Fxy = c(50, 130, 320, 780, 1900))
#' plotdcca(dcca)
#' @export

plotdcca <- function(dcca,seg=FALSE,point=NULL,main=NULL) {
  if(any(dcca[[2]]<0)){
    warning("negative detrended covariance values; cannot take the logarithm")
    return(invisible(NULL))
  } else{
    if(seg){
      ind<-c(rep("A",point),rep("B",length(dcca[[1]])-point))
      s<-log10(dcca[[1]])
      Fl<-log10(dcca[[2]])/2
      df<-cbind.data.frame(s,Fl,ind)
      alfa1<-round(lm(df$Fl~df$s,subset = df$ind=="A")$coefficients[[2]],3)
      alfa2<-round(lm(df$Fl~df$s,subset = df$ind=="B")$coefficients[[2]],3)
      p1<-ggplot2::ggplot(df,ggplot2::aes(x=s,y=Fl,fill=ind))+
        ggplot2::geom_jitter()+
        ggplot2::xlab(expression(log[10](s)))+
        ggplot2::ylab(expression(log[10](F[XY]~(s))))+
        ggplot2::theme_bw()+
        ggplot2::scale_x_continuous(breaks=seq(min(s),max(s),by=0.25))+
        ggplot2::annotate("text",label=deparse(substitute(paste(lambda[1],"=",a),list(a=alfa1))),parse=TRUE,x = mean(s[1:point])-0.2*mean(s),y=median(Fl[1:point]))+
        ggplot2::annotate("text",label=deparse(substitute(paste(lambda[2],"=",b),list(b=alfa2))),parse=TRUE,x = mean(s[point:length(s)]),y=median(Fl[point:length(s)])-0.1)+
        ggplot2::geom_smooth(method = "lm",se=FALSE)+ggplot2::theme(legend.position = "none")+
        ggplot2::ggtitle(main)
      p1
    } else {
      s<-log10(dcca[[1]])
      Fl<-log10(dcca[[2]])/2
      df<-cbind.data.frame(s,Fl)
      alfa<-round(lm(df$Fl~df$s)$coefficients[[2]],3)
      p1<-ggplot2::ggplot(df,ggplot2::aes(x=s,y=Fl))+
        ggplot2::geom_jitter()+
        ggplot2::xlab(expression(log[10](s)))+
        ggplot2::ylab(expression(log[10](F[XY]~(s))))+
        ggplot2::scale_x_continuous(breaks=seq(min(s),max(s),by=0.25))+
        ggplot2::theme_bw()+
        ggplot2::annotate("text",label=deparse(substitute(paste(lambda,"=",a),list(a=alfa))),parse=TRUE,x = unname(quantile(s,0.25)),y=median(Fl))+
        ggplot2::geom_smooth(method = "lm",se=FALSE,show.legend = FALSE)+
        ggplot2::ggtitle(main)
      p1
    }
  }
}
#' Plot rho-DCCA
#'
#' Plot of Detrended Cross Correlation Coefficient
#' @param rdcca is a rhodcca object
#' @param var character. Indicate which pair in rho dcca object you want to plot.
#' @return a plot of Detrended Cross Correlation Analysis.
#' @importFrom dplyr select filter mutate rename row_number n ends_with
#' @examples
#' set.seed(1)
#' d <- data.frame(x = cumsum(rnorm(300)), y = cumsum(rnorm(300)))
#' plotrdcca(rhodcca(d, np = 20), var = "12")
#' @export

plotrdcca <- function(rdcca,var) {
  if(nchar(var)!=2){stop("`var` must be a string of exactly 2 characters")}
  else{

    df<-rdcca %>%
      select(s,ends_with(var)) %>%
      filter(row_number() <= n()-1) %>%
      rename(rho=ends_with(var)) %>%
      mutate(s=log10(s))

    p1 <- ggplot2::ggplot(df, ggplot2::aes(x = s, y = rho)) +
      ggplot2::geom_point() +
      ggplot2::xlab(expression(log[10](s))) +
      ggplot2::ylab(expression(rho ~ DCCA)) +
      ggplot2::scale_x_continuous(breaks = seq(min(df$s),max(df$s), by = 0.25)) +
      ggplot2::theme_bw()
    p1

  }

}

