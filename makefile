FLAGS=-fopenmp
LAPACK=-llapack
PLOT="set term x11 size 1500,500;
PLOT+=set xrange [-10:2110]; set yrange [-0.1:1.1]; set key off;
PLOT+=p 'exact11' u 1:2 w l lc 'black' lw 2, 
PLOT+='exact22' u 1:2 w l lc 'black' lw 2 dt 2,
PLOT+='exact33' u 1:2 w l lc 'black' lw 2 dt 3,
PLOT+='pop.out' u 1:2 w l lc 'red' lw 2, 
PLOT+='pop.out' u 1:3 w l lc 'red' lw 2 dt 2,
PLOT+='pop.out' u 1:4 w l lc 'red' lw 2 dt 3"

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
