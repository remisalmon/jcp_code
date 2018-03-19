#!/bin/bash

PC1="powerspec-remi"
PC2="Newton"
BOND="sc1-m2n6.scsystem"
JFMUCCI="sci-m8n6.scsystem"

if [ $HOSTNAME = $PC1 ] || [ $HOSTNAME = $PC2 ]; then
	rm *.txt *.out code_mpi 2>/dev/null
	cp fmatlab.* code_seq/
	cd code_seq
	echo "compile code seq..."
	gcc *.h *.c -o code_seq -lm
	cd ..
	echo "compile code mpi..."
	mpicc *.h *.c -o code_mpi -lm
	echo "compilation terminee"
fi

if [ $HOSTNAME = $BOND ] || [ $HOSTNAME = $JFMUCCI ]; then
	rm *.txt *.out *.data *.tar code_mpi 2>/dev/null
	ssh remi 'cd /home/remi/Desktop/CA/code/ && cp fmatlab.* code_seq/ && rm *.txt *.out code_mpi 2>/dev/null'
	scp -q param.h remi:/home/remi/Desktop/CA/code/
	echo "compile code seq..."
	ssh remi 'cd /home/remi/Desktop/CA/code/code_seq && gcc *.h *.c -o code_seq -lm'
	echo "compile code mpi..."
	mpicc *.h *.c -o code_mpi -O3 -lscm -lm -std=c99 -lscmpi
	echo "compilation terminee"
fi
