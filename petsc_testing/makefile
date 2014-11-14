
all: main

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

OBJS = OneAndTwoDHeatBTCS.o solver.o

main: main.o ${OBJS} chkopts
	-${CLINKER} -o main main.o ${OBJS} ${PETSC_MAT_LIB}
	${RM} main.o ${OBJS}


