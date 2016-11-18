F90=ifort
FCFLAGS=-O2

TARGET= 1D_RANS
OBJECT= RANS_module.o RANS_main.o RANS_setup.o RANS_poiseuille.o \
		RANS_output.o TDMA_Solver.o

all : $(TARGET)
$(TARGET) : $(OBJECT)
	$(F90) $(FCFLAGS) -o $@ $^

.SUFFIXES. : .o .f90

%.o : %.f90
	$(F90) $(FCFLAGS) -c $<

clean : 
	rm -f *.o
