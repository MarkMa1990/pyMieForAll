%module pyMie

%{
    #define SWIG_FILE_WITH_INIT
    #include "pyMie.hpp"

%}

%include "numpy.i"

%init %{
    import_array();
%}

%apply (double *INPLACE_ARRAY1, int DIM1) {(double *data_out, int Nx)}

%include "pyMie.hpp"

