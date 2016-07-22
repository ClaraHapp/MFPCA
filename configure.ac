# Run this file with autoreconf configure.ac to produce a configure file

# Initialize
AC_INIT([MFPCA], 1.0) dnl package name, version

# use the same C compiler as R CMD xxx
if test -n "${R_HOME}"; then
  CC=`"${R_HOME}/bin/R" CMD config CC`
  CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
else
  echo "could not determine R_HOME"
  exit 1
fi

AC_PROG_CPP # set C preprocessor flags
AC_PROG_CC # set C compiler 

# check if fftw3 is installed
PKG_CHECK_MODULES(FFTW, fftw3,
  [ # fftw3 is installed
    FFTW_LIBS=-lfftw3
    FFTW_CFLAGS="-D HAVE_FFTW" # define makro HAVE_FFTW
  ],
  [ # fftw3 not found 
    FFTW_LIBS=
    AC_MSG_WARN([Library fftw3 not available. Package MFPCA is installed without Fourier/Cosine basis transformation. If you would like to enable this feature, please install fftw3 from http://www.fftw.org/ and reinstall MFPCA.])
  ]) 

# substitute flags for FFTW functionalities
AC_SUBST(FFTW_CFLAGS)
AC_SUBST(FFTW_LIBS)

# substitute all other flags
AC_SUBST(CFLAGS)
AC_SUBST(LIBS)

# write output to Makevars
AC_OUTPUT(src/Makevars)