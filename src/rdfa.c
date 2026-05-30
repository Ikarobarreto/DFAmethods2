#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/******************************************************************************************************************************************************
rdfa.h
*******************************************************************************************************************************************************/

#define SWAP(a,b) {temp = (a); (a) = (b); (b) = temp;}

typedef struct{
int NPTS;	// number points
int NFIT;	// order of the regression fit, plus 1
int IFLAG; 	// integrate the input data if non-zero
int NR;		// number of box sizes
int SW;		// sliding window
int MINBOX;	// minimum box size
int MAXBOX;	// maximum box size
int ABSFLAG;	// calculate absolute values for dcca
}CONFIG;


int iflag=1,minbox,maxbox,nfit=2,npts,nr,sw;
int *rs;	/* box size array; passed in from R as an integer vector */
double *seq;	/* input data buffer; allocated and filled by input() */
double *mse;	/* fluctuation array; allocated by setup(), filled by dfa() */
double *seq1;	// input data
double *dif1;	// difference array
double *seq2;	// input data
double *dif2;	// difference array

//dcca.c variables
double *mse1, *mse2, *msedcca;	/* fluctuation array; allocated by setup(), filled by dfa() */
int absflag;

void setup_dcca(int);
void cleanup_dcca(void);

int *rs;


/* Function prototypes. */
void run_dfa(double *seq, int npts);
void run_dcca(double *seq1, double *seq2, int npts);

int rscale(long minbox, long maxbox, double boxratio);
void dfa(double *seq, long npts, int nfit, int *rs, int nr, int sw);
void dcca(double *seq1, double *seq2, long npts, int nfit, int *rs, int nr, int sw);
void setup(void);
void cleanup(void);
void setup_dcca(int);
void cleanup_dcca(void);
double polyfit(double **x, double *y, double *dif, long ndat, int nfit);
double *vector(long nl, long nh);
int *ivector(long nl, long nh);
long *lvector(long nl, long nh);
double **matrix(long nrl, long nrh, long ncl, long nch);
void free_vector(double *v, long nl, long nh);
void free_ivector(int *v, long nl, long nh);
void free_lvector(long *v, long nl, long nh);
void free_matrix(double **m, long nrl, long nrh, long ncl, long nch);

/******************************************************************************************************************************************************
dfa.c
*******************************************************************************************************************************************************/


void run_dfa(double *seq, int npts)
{

    /* Allocate memory for dfa() and the functions it calls. */
    setup();

    /* Measure the fluctuations of the detrended input data at each box size
       using the DFA algorithm; fill mse[] with these results. */
    dfa(seq, npts, nfit, rs, nr, sw);

    /* Release allocated memory. */
    cleanup();
}



double **x;	/* matrix of abscissas and their powers, for polyfit(). */

/* Detrended fluctuation analysis
    seq:	input data array
    npts:	number of input points
    nfit:	order of detrending (2: linear, 3: quadratic, etc.)
    rs:		array of box sizes (uniformly distributed on log scale)
    nr:		number of entries in rs[] and mse[]
    sw:		mode (0: non-overlapping windows, 1: sliding window)
   This function returns the mean squared fluctuations in mse[].
*/
void dfa(double *seq, long npts, int nfit, int *rs, int nr, int sw)
{
    long i, boxsize, inc, j;
    double stat;

    for (i = 1; i <= nr; i++) {
        boxsize = rs[i];
        if (sw) { inc = 1; stat = (int)(npts - boxsize + 1) * boxsize; }
	else { inc = boxsize; stat = (int)(npts / boxsize) * boxsize; }
        for (mse[i] = 0.0, j = 0; j <= npts - boxsize; j += inc)
            mse[i] += polyfit(x, seq + j, 0, boxsize, nfit);
        mse[i] /= stat;
    }
}

/* workspace for polyfit() */
double *beta, **covar, **covar0;
int *indxc, *indxr, *ipiv;

/* This function allocates workspace for dfa() and polyfit(), and sets
   x[i][j] = i**(j-1), in preparation for polyfit(). */
void setup()
{
    long i;
    int j;

    beta = vector(1, nfit);
    covar = matrix(1, nfit, 1, nfit);
    covar0 = matrix(1, nfit, 1, nfit);
    indxc = ivector(1, nfit);
    indxr = ivector(1, nfit);
    ipiv = ivector(1, nfit);
//    mse = vector(1, nr);
    x = matrix(1, rs[nr], 1, nfit);
    for (i = 1; i <= rs[nr]; i++) {
	x[i][1] = 1.0;
	x[i][2] = i;
	for (j = 3; j <= nfit; j++)
	    x[i][j] = x[i][j-1] * i;
    }
}

/* This function frees all memory previously allocated by this program. */
void cleanup()
{
    free_matrix(x, 1, rs[nr], 1, nfit);
//    free_vector(mse, 1, nr);
    free_ivector(ipiv, 1, nfit);
    free_ivector(indxr, 1, nfit);
    free_ivector(indxc, 1, nfit);
    free_matrix(covar0, 1, nfit, 1, nfit);
    free_matrix(covar, 1, nfit, 1, nfit);
    free_vector(beta, 1, nfit);
//    free_lvector(rs, 1, rslen);	/* allocated by rscale() */
//    free(seq);			/* allocated by input() */
}

/* polyfit() is based on lfit() and gaussj() from Numerical Recipes in C
   (Press, Teukolsky, Vetterling, and Flannery; Cambridge U. Press, 1992).  It
   fits a polynomial of degree (nfit-1) to a set of boxsize points given by
   x[1...boxsize][2] and y[1...boxsize].  The return value is the sum of the
   squared errors (chisq) between the (x,y) pairs and the fitted polynomial.
*/
double polyfit(double **x, double *y, double *res, long boxsize, int nfit)
{
    int icol = 0, irow = 0, j, k;
    double big, chisq, pivinv, temp;
    long i;
    static long pboxsize = 0L;

    /* This block sets up the covariance matrix.  Provided that boxsize
       never decreases (which is true in this case), covar0 can be calculated
       incrementally from the previous value. */
    if (pboxsize != boxsize) {	/* this will be false most of the time */
	if (pboxsize > boxsize)	/* this should never happen */
	    pboxsize = 0L;
	if (pboxsize == 0L)	/* this should be true the first time only */
	    for (j = 1; j <= nfit; j++)
		for (k = 1; k <= nfit; k++)
		    covar0[j][k] = 0.0;
	for (i = pboxsize+1; i <= boxsize; i++)
	    for (j = 1; j <= nfit; j++)
		for (k = 1, temp = x[i][j]; k <= j; k++)
		    covar0[j][k] += temp * x[i][k];
	for (j = 2; j <= nfit; j++)
	    for (k = 1; k < j; k++)
		covar0[k][j] = covar0[j][k];
	pboxsize = boxsize;
    }
    for (j = 1; j <= nfit; j++) {
	beta[j] = ipiv[j] = 0;
	for (k = 1; k <= nfit; k++)
	    covar[j][k] = covar0[j][k];
    }
    for (i = 1; i <= boxsize; i++) {
	beta[1] += (temp = y[i]);
	beta[2] += temp * i;
    }
    if (nfit > 2)
	for (i = 1; i <= boxsize; i++)
	    for (j = 3, temp = y[i]; j <= nfit; j++)
		beta[j] += temp * x[i][j];
    for (i = 1; i <= nfit; i++) {
	big = 0.0;
	for (j = 1; j <= nfit; j++)
	    if (ipiv[j] != 1)
		for (k = 1; k <= nfit; k++) {
		    if (ipiv[k] == 0) {
			if ((temp = covar[j][k]) >= big ||
			    (temp = -temp) >= big) {
			    big = temp;
			    irow = j;
			    icol = k;
			}
		    }
		    else if (ipiv[k] > 1)
				return -1;
		}
	++(ipiv[icol]);
	if (irow != icol) {
	    for (j = 1; j <= nfit; j++) SWAP(covar[irow][j], covar[icol][j]);
	    SWAP(beta[irow], beta[icol]);
	}
	indxr[i] = irow;
	indxc[i] = icol;
	if (covar[icol][icol] == 0.0) 	return -2;
	pivinv = 1.0 / covar[icol][icol];
	covar[icol][icol] = 1.0;
	for (j = 1; j <= nfit; j++) covar[icol][j] *= pivinv;
	beta[icol] *= pivinv;
	for (j = 1; j <= nfit; j++)
	    if (j != icol) {
		temp = covar[j][icol];
		covar[j][icol] = 0.0;
		for (k = 1; k <= nfit; k++) covar[j][k] -= covar[icol][k]*temp;
		beta[j] -= beta[icol] * temp;
	    }
    }
    chisq = 0.0;
    if (nfit <= 2)
	for (i = 1; i <= boxsize; i++) {
	    temp = beta[1] + beta[2] * i - y[i];
	    chisq += temp * temp;
		if(res)
			res[i]=temp;
	}
    else
	for (i = 1; i <= boxsize; i++) {
	    temp = beta[1] + beta[2] * i - y[i];
	    for (j = 3; j <= nfit; j++) temp += beta[j] * x[i][j];
	    chisq += temp * temp;
		if(res)
			res[i]=temp;
	}
    return (chisq);
}


double *vector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
    double *v = (double *)malloc((size_t)((nh-nl+2) * sizeof(double)));
    if (v == NULL)
		return NULL;
    return (v-nl+1);
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
    int *v = (int *)malloc((size_t)((nh-nl+2) * sizeof(int)));
    if (v == NULL)
		return NULL;
    return (v-nl+1);
}

long *lvector(long nl, long nh)
/* allocate a long int vector with subscript range v[nl..nh] */
{
    long *v = (long *)malloc((size_t)((nh-nl+2) * sizeof(long)));
    if (v == NULL)
		return NULL;
    return (v-nl+1);
}

double **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    long i, nrow = nrh-nrl+1, ncol = nch-ncl+1;
    double **m;

    /* allocate pointers to rows */
    m = (double **) malloc((size_t)((nrow+1) * sizeof(double*)));
    if (!m)
		return NULL;
    m += 1;
    m -= nrl;

    /* allocate rows and set pointers to them */
    m[nrl] = (double *) malloc((size_t)((nrow*ncol+1) * sizeof(double)));
    if (!m[nrl])
		return NULL;

    m[nrl] += 1;
    m[nrl] -= ncl;

    for (i = nrl+1; i <= nrh; i++) m[i] = m[i-1]+ncol;

    /* return pointer to array of pointers to rows */
    return (m);
}

void free_vector(double *v, long nl, long nh)
/* free a double vector allocated with vector() */
{
    free(v+nl-1);
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
    free(v+nl-1);
}

void free_lvector(long *v, long nl, long nh)
/* free a long int vector allocated with lvector() */
{
    free(v+nl-1);
}

void free_matrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by matrix() */
{
    free(m[nrl]+ncl-1);
    free(m+nrl-1);
}


/******************************************************************************************************************************************************
dcca.c
*******************************************************************************************************************************************************/

void run_dcca(double *seq1, double *seq2, int npts)
{

    /* Allocate memory for dfa() and the functions it calls. */
    setup_dcca(npts);

    /* Measure the fluctuations of the detrended input data at each box size
       using the DFA algorithm; fill mse[] with these results. */
    dcca(seq1, seq2, npts, nfit, rs, nr, sw);

    /* Release allocated memory. */
    cleanup_dcca();
}


double **x1;	/* matrix of abscissas and their powers, for polyfit(). */
double **x2;	/* matrix of abscissas and their powers, for polyfit(). */

/* Detrended fluctuation analysis
    seq:	input data array
    npts:	number of input points
    nfit:	order of detrending (2: linear, 3: quadratic, etc.)
    rs:		array of box sizes (uniformly distributed on log scale)
    nr:		number of entries in rs[] and mse[]
    sw:		mode (0: non-overlapping windows, 1: sliding window)
   This function returns the mean squared fluctuations in mse[].
*/
void dcca(double *seq1, double *seq2, long npts, int nfit, int *rs, int nr, int sw)
{
    long i, boxsize, inc, j, k;
    double stat, temp;

    for (i = 1; i <= nr; i++) {
        boxsize = rs[i];
        if (sw) { inc = 1; stat = (int)(npts - boxsize + 1) * boxsize; }
	else { inc = boxsize; stat = (int)(npts / boxsize) * boxsize; }

		mse[i]=0.0;
        for (j = 0; j <= npts - boxsize; j += inc)
			{
            polyfit(x1, seq1 + j, dif1+j, boxsize, nfit);
            polyfit(x2, seq2 + j, dif2+j, boxsize, nfit);
			for (k = 1; k <= boxsize; k++)
				{
				temp=dif1[j+k]*dif2[j+k];
				if(absflag)
					temp=fabs(temp);
				mse[i]+=temp;
				}
			}
        mse[i] /= stat;		
    }
}

/* workspace for polyfit() */
double *beta, **covar, **covar0;
int *indxc, *indxr, *ipiv;

/* This function allocates workspace for dfa() and polyfit(), and sets
   x[i][j] = i**(j-1), in preparation for polyfit(). */
void setup_dcca(int npts)
{
    long i;
    int j;

    beta = vector(1, nfit);
    covar = matrix(1, nfit, 1, nfit);
    covar0 = matrix(1, nfit, 1, nfit);
    indxc = ivector(1, nfit);
    indxr = ivector(1, nfit);
    ipiv = ivector(1, nfit);
 	dif1 = vector(1, npts);
 	dif2 = vector(1, npts);
    x1 = matrix(1, rs[nr], 1, nfit);
    x2 = matrix(1, rs[nr], 1, nfit);
    for (i = 1; i <= rs[nr]; i++) {
	x1[i][1]=x2[i][1] = 1.0;
	x1[i][2]=x2[i][2] = i;
	for (j = 3; j <= nfit; j++)
		{
	    x1[i][j] = x1[i][j-1] * i;
	    x2[i][j] = x2[i][j-1] * i;
		}
    }
}

/* This function frees all memory previously allocated by this program. */
void cleanup_dcca()
{
    free_matrix(x1, 1, rs[nr], 1, nfit);
    free_matrix(x2, 1, rs[nr], 1, nfit);
    free_vector(dif1, 1, nr);
    free_vector(dif2, 1, nr);
    free_ivector(ipiv, 1, nfit);
    free_ivector(indxr, 1, nfit);
    free_ivector(indxc, 1, nfit);
    free_matrix(covar0, 1, nfit, 1, nfit);
    free_matrix(covar, 1, nfit, 1, nfit);
    free_vector(beta, 1, nfit);
}




/******************************************************************************************************************************************************
rdfa.c
*******************************************************************************************************************************************************/

void rdfa(CONFIG *cfg, double *r_seq, int *r_rs,double *r_mse)
{
double alpha, scale;
int i;

iflag=cfg->IFLAG;
minbox=cfg->MINBOX;
maxbox=cfg->MAXBOX;
nfit=cfg->NFIT;
npts=cfg->NPTS;
nr=cfg->NR;
sw=cfg->SW;

rs=r_rs;
mse=r_mse;
seq=malloc((npts+1)*sizeof(double));
memcpy(seq+1,r_seq,npts*sizeof(double));
if(iflag)
	{
	for(i=2;i<=npts;i++)
		seq[i]+=seq[i-1];
	}

if(!rs[1])
	{
	alpha=pow(maxbox/(double)minbox,1.0/(nr-1));
	rs[1]=minbox;
	for(i=2;i<=nr;i++)
		{
		scale=alpha;
		while(1)
			{
			rs[i]=scale*rs[i-1];
			if(rs[i]>rs[i-1])
				break;
			scale*=alpha;
			}
		}
	}

run_dfa(seq, npts);

free(seq);
}

/************************************************************************************/

void rdcca(CONFIG *cfg, double *r_seq1, double *r_seq2, int *r_rs,double *r_mse)
{
double alpha, scale;
int i;

iflag=cfg->IFLAG;
minbox=cfg->MINBOX;
maxbox=cfg->MAXBOX;
nfit=cfg->NFIT;
npts=cfg->NPTS;
nr=cfg->NR;
sw=cfg->SW;
absflag=cfg->ABSFLAG;

rs=r_rs;
mse=r_mse;
seq1=malloc((npts+1)*sizeof(double));
memcpy(seq1+1,r_seq1,npts*sizeof(double));
seq2=malloc((npts+1)*sizeof(double));
memcpy(seq2+1,r_seq2,npts*sizeof(double));
if(iflag)
	{
	for(i=2;i<=npts;i++)
		{
		seq1[i]+=seq1[i-1];
		seq2[i]+=seq2[i-1];
		}
	}

if(!rs[1])
	{
	alpha=pow(maxbox/(double)minbox,1.0/(nr-1));
	rs[1]=minbox;
	for(i=2;i<=nr;i++)
		{
		scale=alpha;
		while(1)
			{
			rs[i]=scale*rs[i-1];
			if(rs[i]>rs[i-1])
				break;
			scale*=alpha;
			}
		}
	}

run_dcca(seq1, seq2, npts);

free(seq1);
free(seq2);
}


