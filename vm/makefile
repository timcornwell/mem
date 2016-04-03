

FFLAGS = 
#

mem: memtst.o memsubs.o
	$(FC) memtst.o memsubs.o -o mem
bs.pyf: memsubs.f 
	${F2PY} -h memsubs.pyf -m memsubs memsubs.f

memsubs.so: memsubs.f memsubs.pyf
	${F2PY} ${PYFLAGS} -c -m memsubs memsubs.f memsubs.pyf -lg2c

clean:
	rm memtst.o memsubs.o memtst memsubs.so

memtst: memsubs.o memtst.o
	${FC} ${FFLAGS} memtst.o -o memtst memsubs.o -lg2c -lm

all: memsubs.so memtst