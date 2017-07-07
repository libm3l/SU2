
PATHLSIPDX = /home/jiraseka/Cprograms/Sources/lsipdx/Source
PATHLM3L   = /home/jiraseka/Cprograms/Sources/libm3l/Source/data_util/

all:
	sed 's|PATHLS|$(PATHLSIPDX)|' Elast_Solver._cpp > Elast_Solver.cpp
	sed -i 's|PATHLM|$(PATHLM3L)|' Elast_Solver.cpp
	gcc -O -o Elast_Solver.out Elast_Solver.cpp -L$(PATHLSIPDX) -llsipdx -Wl,-rpath=$(PATHLSIPDX)  -L$(PATHLM3L) -lm3l -Wl,-rpath=$(PATHLM3L) -lpthread -lm


 
