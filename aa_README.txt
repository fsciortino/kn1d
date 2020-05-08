	All source code for KN1D is included in the zip file KN1D.zip.
When you unzip the files, they all will appear at the current directory
level. The manual for KN1D is included in the file kn1d_manual.pdf.

The code is written in IDL but does have a few "call_external" calls to
FORTRAN routines. These routines are basically 2D spline interpolation routines
that are used in looking up ionization rates and line radatiation emissivities
from the Johnson-Hinnov coefficients (collisional-radiative model).
It is possible to run the code without using these routines, but I will assume
here that they are desired. Only two shared object modules need to be compiled
from the FORTRAN source code that is included in this package. The following
commands do the trick for Linux (these commands are in the file make_kn1d_so.sh):

   g77 fast_b2val.f b2val.f xerror.f i1mach.f fdump.f -shared -o fast_b2val.so -w
   g77 call_b2ink.f b2ink.f xerror.f i1mach.f fdump.f -shared -o call_b2ink.so -w

If you chose to run the IDL code in a different directory than where the .so files
reside, then you will need to add the .so file  path too your LD_LIBRARY_PATH
environment variable.

Once you are ready to try the code, I suggest that you run "test_kn1d" first:

IDL> .run test_kn1d 

By look at test_kn1d.pro (and other test_kn1d*.pro files) you can get
an idea of how to set up the code to run.

I think I included all the code you need, but if the code complains about
an undefined function/procedure then let me know.

Best of luck,

Brian LaBombard  April 21, 2003

