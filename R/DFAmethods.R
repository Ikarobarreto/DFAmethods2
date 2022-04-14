#' Multiple Fractal Regression
#'
#' Calculates Fractal Regression.
#' @param data is a matrix of time series.
#' @param dpo Detrending polynomial order.
#' @param int logical. if TRUE integration process will be applied.
#' @param np number of point scales.
#' @param overlap logical. if TRUE overlapping windows will be applied.
#' @return Scale s, Detrended Fluctuation Function F, Rho-DCCA, Rho DPCCA,
#' Beta DFA estimates, Standardized Beta DFA estimates, DFA Residuals, DFA Variance, DFA Upper and Lower confidence interval,
#' Multiple Detrended Correlation, DFA R² , DFA p-value and DFA Calculated T statistics.
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
#' @useDynLib DFAmethods2, .registration=TRUE

fracreg<-function(data,dpo,int,np=91,overlap=T){
	data<-as.matrix(data)
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
	tn<-array(,dim=c((nc-1),1,np)) #p-valor para os Bdfas
	tnc<-array(,dim=c((nc-1),1,np)) #T crítico
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
ans1<- .C("rdfa", cfg, as.numeric(seq1), as.integer(rsi1), as.numeric(msei1),PACKAGE = "DFAmethods2")
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
ans12 <- .C("rdcca", cfg, as.numeric(seq1), as.numeric(seq2), as.integer(rs12), as.numeric(mse12),PACKAGE = "DFAmethods2")
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
ans1<- .C("rdfa", cfg, as.numeric(u), as.integer(rsi1), as.numeric(msei1),PACKAGE = "DFAmethods2")
un[1,1,k]<-ans1[[4]][2]
dmc2[k]<-pn[1,(2:nc),k]%*%solve(pn[(2:nc),(2:nc),k])%*%pn[(2:nc),1,k]
bs[,,k]<-solve(pn[(2:nc),(2:nc),k])%*%pn[(2:nc),1,k]
rn[1,1,k]<-1-(un[1,1,k]/fn[1,1,k])
u<-NULL
u1<-NULL
}
for(k in 1:np){
	for(i in (2:nc)){
	  vn[,(i-1),k]<-((fn[1,1,k]-dmc2[k]*fn[1,1,k])/(fn[i,i,k]))*(1/(sn[[k]]-nc-2))
	  #vn[,(i-1),k]<-(un[1,1,k]/fn[i,i,k])*(1/(size))
		}}
for(k in 1:np){
	for(i in 1:(nc-1)){
	uci[i,1,k]<-as.matrix(bn[i,1,k])+qt(0.975,df=(size/sn[[k]])-nc)*sqrt(vn[1,i,k])
	lci[i,1,k]<-as.matrix(bn[i,1,k])-qt(0.975,df=(size/sn[[k]])-nc)*sqrt(vn[1,i,k])
	tn [i,1,k]<-as.matrix(1-pt(abs(bn[i,1,k])/sqrt(vn[1,i,k]),df=((size/sn[[k]])-nc)))
	tnc[i,1,k]<-as.matrix(qt(0.975,df=((size/sn[[k]])-nc))*sqrt(vn[1,i,k]))
	}}

fracreg<-list(sn,fn,pn,dp,bn,bs,un,vn,dmc2,rn,uci,lci,tn,tnc)
names(fracreg)<-c("s","F","DCCA","DPCCA","BDFA","BSDFA","UDFA","VDFA","DMC2","R2DFA","UCIB","LCIB","p.value","TC")

return(fracreg)
}

#' Detrended Fluctuation Analysis
#'
#' Calculates DFA
#' @param data is a vector of time series
#' @param dpo detrending polynomial order
#' @param int logical. if TRUE integration process will be applied.
#' @param np number of point scales.
#' @param overlap logical. if TRUE overlapping windows will be applied.
#' @return Scale s, Detrended Fluctuation Function F
#' @useDynLib DFAmethods2, .registration=TRUE
#' @export

dfa<-function(data,dpo=1,int=T,np=91,overlap=T){
	data<-as.matrix(data)
	if (int ==TRUE){int=1} else{int=0}
	dpo<-as.numeric(dpo+1)
	nc<-ncol(data)
	sn<-matrix(,ncol=1) #Escalas
	mx<-round(nrow(data)/5)
	fn<-matrix(,ncol=1,nrow=np) #Variancia-Covarância sem tendência
	size=nrow(data)
#DFA
seq1<-data[,1]
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
ans1<- .C("rdfa", cfg, as.numeric(seq1), as.integer(rsi1), as.numeric(msei1),PACKAGE = "DFAmethods2")
sn<- ans1[[3]][1:NR+1]
fn<-ans1[[4]][1:NR+1]
DFA<-tibble::as_tibble(cbind.data.frame(sn,fn))
names(DFA)<-c("s","F")
return(DFA)
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
#' @importFrom tibble as_tibble
#' @importFrom gdata upperTriangle
#' @useDynLib DFAmethods2, .registration=TRUE
#' @export
#'
rhodcca<-function(data,dpo=1,int=T,np=91,overlap=T){
	data<-as.matrix(data)
	if(ncol(data)<2){warning("A matrix with at least two columns are necessary")}
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
ans1<- .C("rdfa", cfg, as.numeric(seq1), as.integer(rsi1), as.numeric(msei1),PACKAGE = "DFAmethods2")
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
ans12 <- .C("rdcca", cfg, as.numeric(seq1), as.numeric(seq2), as.integer(rs12), as.numeric(mse12),PACKAGE = "DFAmethods2")
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
#' @useDynLib DFAmethods2, .registration=TRUE
#' @export
#'
rhodpcca<-function(data,dpo=1,int=T,np=91,overlap=T){
	data<-as.matrix(data)
	if(ncol(data)<3){warning("A matrix with at least three columns are necessary")}
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
ans1<- .C("rdfa", cfg, as.numeric(seq1), as.integer(rsi1), as.numeric(msei1),PACKAGE = "DFAmethods2")
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
ans12 <- .C("rdcca", cfg, as.numeric(seq1), as.numeric(seq2), as.integer(rs12), as.numeric(mse12),PACKAGE = "DFAmethods2")
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
#' @useDynLib DFAmethods2, .registration=TRUE
#' @export

dmc2<-function(data,dpo=1,int=T,np=91,overlap=T){
	data<-as.matrix(data)
	if(ncol(data)<3){warning("A matrix with at least three columns are necessary")}
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
ans1<- .C("rdfa", cfg, as.numeric(seq1), as.integer(rsi1), as.numeric(msei1),PACKAGE = "DFAmethods2")
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
ans12 <- .C("rdcca", cfg, as.numeric(seq1), as.numeric(seq2), as.integer(rs12), as.numeric(mse12),PACKAGE = "DFAmethods2")
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
#' @useDynLib DFAmethods2, .registration=TRUE
#' @importFrom dplyr inner_join
#' @export

effsizeDFA<-function(data,dpo=1,int=T,np=91,overlap=T){
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
#' @useDynLib DFAmethods2, .registration=TRUE

bdfa<-function(data,dpo=1,int=T,np=91,overlap=T){
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
ans1<- .C("rdfa",cfg,as.numeric(seq1),as.integer(rsi1),as.numeric(msei1),PACKAGE = "DFAmethods2")
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
ans12 <- .C("rdcca", cfg, as.numeric(seq1), as.numeric(seq2), as.integer(rs12), as.numeric(mse12),PACKAGE = "DFAmethods2")
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
#' @useDynLib DFAmethods2, .registration=TRUE
#' @export

betadfa<-function(data,dpo=1,int=T,np=91,overlap=T){
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
    ans1<- .C("rdfa",cfg,as.numeric(seq1),as.integer(rsi1),as.numeric(msei1),PACKAGE = "DFAmethods2")
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
    ans12 <- .C("rdcca", cfg, as.numeric(seq1), as.numeric(seq2), as.integer(rs12), as.numeric(mse12),PACKAGE = "DFAmethods2")
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
#' @useDynLib DFAmethods2, .registration=TRUE
#' @export

sbdfa<-function(data,dpo=1,int=T,np=91,overlap=T){
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
    ans1<- .C("rdfa",cfg,as.numeric(seq1),as.integer(rsi1),as.numeric(msei1),PACKAGE = "DFAmethods2")
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
    ans12 <- .C("rdcca", cfg, as.numeric(seq1), as.numeric(seq2), as.integer(rs12), as.numeric(mse12),PACKAGE = "DFAmethods2")
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
#' @useDynLib DFAmethods2, .registration=TRUE
#' @export

fracreg.Ktest<-function(data,B=100,dpo=1,int=T,np=91,overlap=T){
x<-ncol(data)-1
if (int ==TRUE){int=1} else{int=0}
yx<-ncol(data)
pbdfa<-array(,dim=c(x,1,np,B))
data<-as.matrix(data)
res<-lm(data[,1]~data[,2:yx])$res
coef<-lm(data[,1]~data[,2:yx])$coeff
usur<-surrogate(res,ns=B,fft=T,amplitude=T)
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
if(x==1){beta<-as.matrix(dpo(data,dpo=dpo)$BDFA[,1,])
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
#' @useDynLib DFAmethods2, .registration=TRUE
#' @export

fracreg.PStest<-function(data,B=100,dpo=1,int=T,np=91,overlap=T){
  a<-bdfa(data,dpo=dpo,int=int,np=np,overlap=overlap)
	n<-nrow(data)
	m<-ncol(data)-1
	test<-array(0,dim=c(m,1,np,B))
	usur<-surrogate(data[,1],ns=B,fft=T,amplitude=T)
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
cib<-tibble::as.tibble(cbind.data.frame(beta,q025,q975,s))
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
#' @useDynLib DFAmethods2, .registration=TRUE
#' @export

fracreg.IUTest<-function(data,B=100,dpo=1,int=T,np=91,overlap=T){
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
#' @export

plotdfa <- function(dfa,seg=F,point=NULL,main) {
  if(seg==T){
    ind<-c(rep("A",point),rep("B",length(dfa[[1]])-point))
    s<-log10(dfa[[1]])
    F<-log10(dfa[[2]])/2
    df<-cbind.data.frame(s,F,ind)
    alfa1<-round(lm(df$F~df$s,subset = df$ind=="A")$coefficients[[2]],3)
    alfa2<-round(lm(df$F~df$s,subset = df$ind=="B")$coefficients[[2]],3)
    p1<-ggplot2::ggplot(df,ggplot2::aes(x=s,y=F,fill=ind))+
      ggplot2::geom_jitter()+
      ggplot2::xlab(expression(log[10](s)))+
      ggplot2::ylab(expression(log[10](F[X]~(s))))+
      ggplot2::scale_x_continuous(breaks=seq(min(s),max(s),by=0.25))+
      ggplot2::theme_bw()+
      ggplot2::annotate("text",label=substitute(paste(alpha[1],"=",a),list(a=alfa1)),x = mean(s[1:point])-0.2*mean(s),y=median(F[1:point]))+
      ggplot2::annotate("text",label=substitute(paste(alpha[2],"=",b),list(b=alfa2)),x = mean(s[point:length(s)]),y=median(F[point:length(s)])-0.1)+
      ggplot2::geom_smooth(method = "lm",se=F)+ggplot2::theme(legend.position = "none")+
      ggplot2::ggtitle(main)
    p1
  } else {
    s<-log10(dfa[[1]])
    F<-log10(dfa[[2]])/2
    df<-cbind.data.frame(s,F)
    alfa<-round(lm(df$F~df$s)$coefficients[[2]],3)
    p1<-ggplot2::ggplot(df,ggplot2::aes(x=s,y=F))+
      ggplot2::geom_jitter()+
      ggplot2::xlab(expression(log[10](s)))+
      ggplot2::ylab(expression(log[10](F[X]~(s))))+
      ggplot2::scale_x_continuous(breaks=seq(min(s),max(s),by=0.25))+
      ggplot2::theme_bw()+ggplot2::annotate("text",label=substitute(paste(alpha,"=",a),list(a=alfa)),x = quantile(s,0.25),y=median(F))+
      ggplot2::geom_smooth(method = "lm",se=F,show.legend = F)+
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
#' @export

plotdcca <- function(dcca,seg=F,point=NULL,main) {
  if(dcca[[2]]<0){
    print("Negative Detredended Covariance Function values")
  } else{
    if(seg==T){
      ind<-c(rep("A",point),rep("B",length(dcca[[1]])-point))
      s<-log10(dcca[[1]])
      F<-log10(dcca[[2]])/2
      df<-cbind.data.frame(s,F,ind)
      alfa1<-round(lm(df$F~df$s,subset = df$ind=="A")$coefficients[[2]],3)
      alfa2<-round(lm(df$F~df$s,subset = df$ind=="B")$coefficients[[2]],3)
      p1<-ggplot2::ggplot(df,ggplot2::aes(x=s,y=F,fill=ind))+
        ggplot2::geom_jitter()+
        ggplot2::xlab(expression(log[10](s)))+
        ggplot2::ylab(expression(log[10](F[XY]~(s))))+
        ggplot2::theme_bw()+
        ggplot2::scale_x_continuous(breaks=seq(min(s),max(s),by=0.25))+
        ggplot2::annotate("text",label=substitute(paste(lambda[1],"=",a),list(a=alfa1)),x = mean(s[1:point])-0.2*mean(s),y=median(F[1:point]))+
        ggplot2::annotate("text",label=substitute(paste(lambda[2],"=",b),list(b=alfa2)),x = mean(s[point:length(s)]),y=median(F[point:length(s)])-0.1)+
        ggplot2::geom_smooth(method = "lm",se=F)+ggplot2::theme(legend.position = "none")+
        ggplot2::ggtitle(main)
      p1
    } else {
      s<-log10(dcca[[1]])
      F<-log10(dcca[[2]])/2
      df<-cbind.data.frame(s,F)
      alfa<-round(lm(df$F~df$s)$coefficients[[2]],3)
      p1<-ggplot2::ggplot(df,ggplot2::aes(x=s,y=F))+
        ggplot2::geom_jitter()+
        ggplot2::xlab(expression(log[10](s)))+
        ggplot2::ylab(expression(log[10](F[XY]~(s))))+
        ggplot2::scale_x_continuous(breaks=seq(min(s),max(s),by=0.25))+
        ggplot2::theme_bw()+
        ggplot2::annotate("text",label=substitute(paste(lambda,"=",a),list(a=alfa)),x = quantile(s,0.25),y=median(F))+
        ggplot2::geom_smooth(method = "lm",se=F,show.legend = F)+
        ggplot2::ggtitle(main)
      p1
    }
  }
}
#' Plot rho-DCCA
#'
#' Plot of Detrended Cross Correlation Coefficient
#' @param rdcca is a rhodcca object
#' @param main plot title
#' @param var numerical. Indicate which pair in rho dcca object you want to plot.
#' @return a plot of Detrended Cross Correlation Analysis.
#' @export

plotrdcca <- function(rdcca,var,main) {
  lim<-length(rdcca[[1]])-1
  s<-log10(rdcca[[1]][1:lim])
  rho<-rdcca[[var+1]][1:lim]
  df<-cbind.data.frame(s,rho)
  p1<-ggplot2::ggplot(df,ggplot2::aes(x=s,y=rho))+
    ggplot2::geom_point()+ggplot2::xlab(expression(log[10](s)))+
    ggplot2::ylab(expression(rho~DCCA))+
    ggplot2::scale_x_continuous(breaks=seq(min(s),max(s),by=0.25))+
    ggplot2::theme_bw()+ggplot2::ggtitle(label=main,subtitle = names(rdcca)[[var+1]])
  p1
}

