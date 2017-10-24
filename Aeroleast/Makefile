
PATHLSIPDX = /home/jiraseka/Cprograms/Sources/lsipdx/Source/
PATHLM3L   = /home/jiraseka/Cprograms/Sources/libm3l/Source/data_util/

all:
	sed 's|PATHLS|$(PATHLSIPDX)|' Elast_Solver._h > Elast_Solver.h
	sed -i 's|PATHLM|$(PATHLM3L)|' Elast_Solver.h
	gcc -O -o Elast_Solver.out Elast_Solver_cspode.c -L$(PATHLSIPDX) -llsipdx -Wl,-rpath=$(PATHLSIPDX)  -L$(PATHLM3L) -lm3l -Wl,-rpath=$(PATHLM3L) -lpthread -lm


 
