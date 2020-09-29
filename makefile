FLAGS=-fopenmp
LAPACK=-llapack
PLOT="set term x11 size 1000,333;
PLOT+=set xrange [-10:2110]; set yrange [-0.1:1.1]; set key off;
PLOT+=p 'exact' u 1:2 w l lc 'black' lw 3, 'exact' u 1:3 w l lc 'black' lw 3 dt 2,
PLOT+='rho.out' u 1:2 w l lc 'red' lw 2, 'rho.out' u 1:7 w l lc 'red' lw 2 dt 2"

all: ranlux.o vars.o fehr.o
	gfortran -O3 -o fehr.x $^ $(LAPACK) $(FLAGS)

%.o: %.f95
	gfortran -O3 -c $^ $(FLAGS)

clean:
	@rm -fv *.x
	@rm -fv *.o
	@rm -fv *.out
	@rm -fv *.dat
	@rm -fv *.mod
	@rm -fv *.log

plot:
	@gnuplot -p -e $(PLOT)
