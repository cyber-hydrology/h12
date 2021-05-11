FC    = mpiifort
SRCS  = globals.f90 2D_ADP_FDM.f90 ascout.f90 flux.f90 forward.f90 inc_cell.f90 initiald.f90 mdvel.f90 mesh_gen.f90 nei_find.f90 paraout.f90 rdt.f90 regularization.f90 sdp.f90 seed_find_divide.f90 store_id_info.f90 suisin.f90 radar.f90
OBJS  = ${SRCS:.f90=.o}
#HEAD = header.h

FLAGS   = -O2 -qopenmp 
OPTION  = -O2 -qopenmp 

all : run
####################################
run : ${OBJS}
	${FC} ${OPTION} ${OBJS} -o $@

${OBJS} : ${SRCS} 
	${FC} ${OPTION} -c ${SRCS}  

####################################
clean : 
	\rm -f *.f90~ *.mod  *.o Makefile~ run
