FFLAGS = -g
FC = f77
FOBJ=nr.o
dmemimg: nmem.o dmemimg.o ${FOBJ}
	${FC} ${FFLAGS} dmemimg.o nmem.o ${FOBJ} -lg2c -lm -o dmemimg

