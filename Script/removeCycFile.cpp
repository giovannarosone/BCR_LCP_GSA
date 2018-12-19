#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdio.h>
#include <assert.h>
#include <sys/stat.h>
#include <math.h>
#include <string.h> 
#include "../Parameters.h" 


int main(int argc, char *argv[]) {
	if( argc != 3 )
    {
      std::cerr << "usage: " << argv[0] << " path #seqs" << std::endl;
	  	std::cerr << "where: " << std::endl;
	  	std::cerr << "\tpath --> the path that contains cyc.*.txt files " << std::endl;
		std::cerr << "\t#seqs --> number of cyc files (it coincides with the max length of the sequence in the collection)" << std::endl;
	    exit(1);
    }

	char *filename = new char[strlen(argv[1])+30];
	dataTypeNSeq len = atoi(argv[2]);
	
	for (int j = 0 ; j < len; j++) {
		sprintf (filename, "%scyc.%u.txt", argv[1], j);
		if (remove(filename)!=0) {
			std::cerr << "Error deleting " << filename << " file" << std::endl;
			exit(1);
		}
		else
			std::cerr << "Removing: " << filename << "\n"; 
	}

	return 0;
}