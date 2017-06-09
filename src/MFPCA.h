#ifndef MFPCA_H  /* make sure header is loaded only once */
#define MFPCA_H

#ifdef HAVE_FFTW
#include <fftw3.h>
#endif

#include <R.h> 
#include <Rmath.h>

#define rho(i) ( (i) == 0 ? 1/sqrt(M_PI) : sqrt(2/M_PI) )
#define z(i) ( (i) == 0 ? 1 : 1/sqrt(2) )

/* C-Functions for MFPCA package  */
void calcCoefs(int* M, int* N, double* image, double* coefs);
void calcImage(int* M, int* N, double* coefs, double* image);

void calcCoefs3D(int* dim, double* image, double* coefs);
void calcImage3D(int* dim, double* coefs, double* image);

#endif /* MFPCA_H */