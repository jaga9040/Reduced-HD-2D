FC = gfortran
LFLAGS= -O3 -I/home/application/fftw/3.3.8/include
# for debugger -g

LIBRAIRIES = -L/home/application/fftw/3.3.8/lib -lfftw3 -lm 

OBJS = RHD2D.o 

Main     :	$(OBJS)
	$(FC) -o a.out $(LFLAGS) $(OBJS) $(LIBRAIRIES)

RHD2D.o : RHD2D.f90 
	$(FC) -c $(LFLAGS) RHD2D.f90

clean :
	/bin/rm -f *.o *~ fort.* a.out

