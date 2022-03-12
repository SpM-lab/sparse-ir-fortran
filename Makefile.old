FC = gfortran
#FC = ifort -I/opt/intel/composer_xe_2015.5.223/ipp/include:/opt/intel/composer_xe_2015.5.223/mkl/include:/opt/intel/composer_xe_2015.5.223/tbb/include

TARGET = test_sparse_ir
OBJS = test_sparse_ir.o sparse_ir.o sparse_ir_io.o sparse_ir_preset.o
LDLIBS = -lblas -llapack
#LDLIBS = -mkl=sequential
#LDFLAGS = -L/usr/lib/aarch64-linux-gnu
LDFLAGS = -L/opt/homebrew/Cellar/openblas/0.3.20/lib
#LDFLAGS = -static-intel

.SUFFIXES: .f90

%.o: %.f90
	$(COMPILE.f) $(OUTPUT_OPTION) $< -std=f95 -Wall -fbounds-check -g -fcheck=array-temps,bounds,do,mem,pointer,recursion
	#$(COMPILE.f) $(OUTPUT_OPTION) $< -fbounds-check
	#$(COMPILE.f) $(OUTPUT_OPTION) $< -g
	#$(COMPILE.f) $(OUTPUT_OPTION) $< -O2 -g -traceback

%.mod: %.f90 %.o
	@:

$(TARGET): $(OBJS)
	$(LINK.f) $^ $(LDLIBS) $(LDFLAGS) -o $@

test_sparse_ir.o: sparse_ir.mod sparse_ir_io.mod sparse_ir_preset.mod

test: $(TARGET)
	./test_sparse_ir

clean:
	rm -f *.o *.mod
