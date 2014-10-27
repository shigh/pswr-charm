
all: OneAndTwoDHeatBTCS

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

OneAndTwoDHeatBTCS: OneAndTwoDHeatBTCS.o  chkopts
	-${CLINKER}  -o OneAndTwoDHeatBTCS OneAndTwoDHeatBTCS.o ${PETSC_MAT_LIB}
	${RM} OneAndTwoDHeatBTCS.o


