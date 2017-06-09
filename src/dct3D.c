// compile with R CMD SHLIB dct.c -lfftw3

#include "MFPCA.h"

void calcCoefs3D(int* dim, double* image, double* coefs)
{	
	# ifdef HAVE_FFTW
	// make the plan
	fftw_plan plan =  fftw_plan_r2r_3d(dim[0], dim[1], dim[2], image, coefs,
	                                   FFTW_REDFT10, FFTW_REDFT10, FFTW_REDFT10,
	                                   0);
	// run dct
	fftw_execute(plan);

	double c = M_PI * M_PI * M_PI /(8 * dim[0] * dim[1] * dim[2]);

	// correct for orthonormal coefficients:
	for(unsigned int i0 = 0; i0 < dim[0]; i0++)
	{
		for(unsigned int i1 = 0; i1 < dim[1]; i1++)
		{
			for(unsigned int i2 = 0; i2 < dim[2]; i2++)
				coefs[dim[1]*dim[2]*i0 + dim[2]*i1 + i2] = c * rho(i0) * rho(i1) * rho(i2) * coefs[dim[1]*dim[2]*i0 + dim[2]*i1 + i2];
		}
	}
	# else
		error("dctBasis3D requires C-library fftw3 to be installed. Check http://www.fftw.org/ for more information.");
	# endif

	return;
}


void calcImage3D(int* dim, double* coefs, double* image)
{	
	# ifdef HAVE_FFTW
	// transform coefs
	for(unsigned int i0 = 0; i0 < dim[0]; i0++)
	{
		for(unsigned int i1 = 0; i1 < dim[1]; i1++)
		{
			for(unsigned int i2 = 0; i2 < dim[2]; i2++)
			coefs[dim[1]*dim[2]*i0 + dim[2]*i1 + i2] = z(i0) * z(i1) * z(i2) * coefs[dim[1]*dim[2]*i0 + dim[2]*i1 + i2];
		}
	}

	
	// make the plan
	fftw_plan plan =  fftw_plan_r2r_3d(dim[0], dim[1], dim[2], coefs, image,
	                                   FFTW_REDFT01, FFTW_REDFT01, FFTW_REDFT01,
	                                   0);
	// run dct
	fftw_execute(plan);

	// correct for orthonormality:
	for(unsigned int i0 = 0; i0 < dim[0]; i0++)
	{
		for(unsigned int i1 = 0; i1 < dim[1]; i1++)
		{
			for(unsigned int i2 = 0; i2 < dim[2]; i2++)
				image[dim[1]*dim[2]*i0 + dim[2]*i1 + i2] = image[dim[1]*dim[2]*i0 + dim[2]*i1 + i2] / (M_PI * sqrt(M_PI));
		}
	}
	# else
		error("dctBasis3D requires C-library fftw3 to be installed. Check http://www.fftw.org/ for more information.");
	# endif

	return;
}
