#include "../param.h"
#include "fmatlab.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[])
{
	//system("cd /home/remi/Desktop/CA/code && octave -qf mpi_gatherdata.m 2>/dev/null");
	system("cd /home/remi/Desktop/CA/code && matlab -nosplash -r 'mpi_gatherdata;exit'");
	
	//system("cd /home/remi/Desktop/CA/code && rm data*.txt 2>/dev/null");
	
	return(0);
}
