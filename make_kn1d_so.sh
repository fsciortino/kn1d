# make_kn1d_so
#  These commands compile fortran files for kn1d into shared object modules
#
   g77 fast_b2val.f b2val.f xerror.f i1mach.f fdump.f -shared -o fast_b2val.so -w
   g77 call_b2ink.f b2ink.f xerror.f i1mach.f fdump.f -shared -o call_b2ink.so -w

