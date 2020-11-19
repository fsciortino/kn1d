# Usage:
# make           # generate KN1D libraries
# make clean     # delete previous versions

.PHONY: all kn1d clean

FC = gfortran
FLAGS = -fPIC -w -shared

# directory where shared-object library should be created
OUT_DIR = $(KN1D_DIR)/  #/home/sciortino/kn1d
#./wkn1d

# directory of Fortran code (NB: KN1D_DIR must be defined!)
KN1DF = $(KN1D_DIR)/fortran


all: kn1d

kn1d :
	@echo "Generating KN1D shared-object libraries"
	@echo "Compiler flags: " ${FLAGS}
	@echo "Building from Fortran source in " ${KN1DF}
	@echo "OUT_DIR: " ${OUT_DIR}
	$(FC) $(FLAGS) $(KN1DF)/fast_b2val.f $(KN1DF)/b2val.f $(KN1DF)/xerror.f \
			$(KN1DF)/i1mach.f $(KN1DF)/fdump.f -o $(KN1DF)/../fast_b2val.so
	$(FC) $(FLAGS) $(KN1DF)/call_b2ink.f $(KN1DF)/b2ink.f $(KN1DF)/xerror.f \
			$(KN1DF)/i1mach.f $(KN1DF)/fdump.f -o $(KN1DF)/../call_b2ink.so  #$(OUT_DIR)/call_b2ink.so	

clean :
	@echo "Eiminating KN1D shared-object library"
	rm $(KN1DF)/../*.so    ##rm $(OUT_DIR)/*.so

