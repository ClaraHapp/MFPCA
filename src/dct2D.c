// compile with R CMD SHLIB dct.c -lfftw3

# ifdef HAVE_FFTW
# include <fftw3.h>
# endif

#include <R.h> 
#include <Rmath.h>

#define rho(i) ( (i) == 0 ? 1/sqrt(M_PI) : sqrt(2/M_PI) )
#define z(i) ( (i) == 0 ? 1 : 1/sqrt(2) )

void calcCoefs(int* M, int* N, double* image, double* coefs)
{	
	# ifdef HAVE_FFTW
	// make the plan
	fftw_plan plan =  fftw_plan_r2r_2d(*M, *N, image, coefs,
	                                   FFTW_REDFT10, FFTW_REDFT10,
	                                   0);
	// run dct
	fftw_execute(plan);

	double c = M_PI * M_PI /(4 * *M * *N);

	// correct for orthonormal coefficients:
	for(unsigned int i = 0; i < *M; i++)
	{
		for(unsigned int j = 0; j < *N; j++)
			coefs[i * (*N) + j] = c * rho(i) * rho(j) * coefs[i * (*N) + j];
	}
	# else
		error("dctBasis2D requires C-library fftw3 to be installed. Check http://www.fftw.org/ for more information.");
	# endif

	return;
}


void calcImage(int* M, int* N, double* coefs, double* image)
{	
	# ifdef HAVE_FFTW
	// transform coefs
	for(unsigned int i = 0; i < *M; i++)
	{
		for(unsigned int j = 0; j < *N; j++)
			coefs[i * (*N) + j] = z(i) * z(j) * coefs[i * (*N) + j];
	}

	
	// make the plan
	fftw_plan plan =  fftw_plan_r2r_2d(*M, *N, coefs, image,
	                                   FFTW_REDFT01, FFTW_REDFT01,
	                                   0);
	// run dct
	fftw_execute(plan);

	// correct for orthonormality:
	for(unsigned int i = 0; i < *M; i++)
	{
		for(unsigned int j = 0; j < *N; j++)
			image[i * (*N) + j] = image[i * (*N) + j] / M_PI;
	}
	# else
		error("dctBasis2D requires C-library fftw3 to be installed. Check http://www.fftw.org/ for more information.");
	# endif

	return;
}
