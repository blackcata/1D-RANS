F90=ifort
TARGET= RANS_1D
OBJECT= RANS_module.o RANS_main.o RANS_setup.o 

all : $(TARGET)
$(TARGET) : $(OBJECT)
	$(F90) -o $@ $^

.SUFFIXES. : .o .f90

%.o : %.f90
	$(F90) -c $<

clean : 
	rm -f *.o
