/**
 ** BCR is part of:
 ** BEETL: Burrows-Wheeler Extended Tool Library
 ** Documentation in: doc/BEETL.md
 **
 ** Copyright (c) 2011-2014 Illumina, Inc. **
 ** BEETL software package is
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** Citations: 
 ** 
 ** Markus J. Bauer, Anthony J. Cox and Giovanna Rosone
 ** Lightweight BWT Construction for Very Large String Collections.
 ** Proceedings of CPM 2011, pp.219-231
 ** 
 ** Markus J. Bauer, Anthony J. Cox, Giovanna Rosone and Marinella Sciortino
 ** Lightweight LCP Construction for Next-Generation Sequencing Datasets. 
 ** Proceedings of WABI 2012, pp 326-337, 2012
 
 ** Markus J. Bauer, Anthony J. Cox, Giovanna Rosone 
 ** Lightweight algorithms for constructing and inverting the BWT of string collections. 
 ** Theoretical Computer Science 483: 134-148 (2013)
 **  
 ** Anthony J. Cox, Fabio Garofalo, Giovanna Rosone, Marinella Sciortino
 ** Lightweight LCP construction for very large collections of strings. 
 ** Journal of Discrete Algorithms (2016) 37: 17-33
 **
 ** By Giovanna Rosone
 **
 **/
 
#include "TransposeFasta.h"
#include "Tools.h"
#include "Parameters.h"
#include <assert.h>
#include <stdio.h>
#include <sys/stat.h>
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>      // std::stringstream

#include <zlib.h>

#if KSEQ_PARSER==1
	#include "external/kseq/kseq.h"
	KSEQ_INIT(gzFile, gzread)
#endif

using std::cout;
using std::cerr;
using std::endl;

TransposeFasta::TransposeFasta()
{

}


TransposeFasta::~TransposeFasta()
{

}

bool TransposeFasta::findLengthNseq( const string& input, const string& fileOutput, string BCRprefPrev){
	//#if BCR_FROMCYC==1
	for (dataTypedimAlpha z = 0 ; z < SIZE_ALPHA-1; z++)
		freq[z]=0;
	freq[SIZE_ALPHA-1]=0;
	freq[(unsigned int)(TERMINATE_CHAR)]=1;
	//#endif

	//lengthRead = CYCLENUM;
    
	char *fileLen = new char[fileOutput.length()+10];
	sprintf(fileLen, "%s.len", fileOutput.c_str());
	
	

	string bufChar;
	
	bool lenSeq = false;
   	vector <dataTypelenSeq> lengthSeqVector;

	//Find max length and number of reads
	dataTypelenSeq charsNumber=0;
	lengthRead=0;
	nSeq = 0;
	lengthTexts = 0;       //Total length of all texts without $-symbols
	
	#if KSEQ_PARSER==0
	
		std::ifstream infile(input.c_str());
		
		if (infile.is_open() == false) {
			std::cerr << "TransposeFasta::findLengthNseq: could not open file \"" << input.c_str() << "\"!"<< std::endl;	
			exit (EXIT_FAILURE);
		}
		while (getline(infile, bufChar))  {
			if ( (bufChar.length() > 0) && ((bufChar[bufChar.length()-2] == '\r') || (bufChar[bufChar.length()-2] == '\n')) )
				bufChar[bufChar.length()-2] = '\0';
			else if ( (bufChar.length() > 0) && ((bufChar[bufChar.length()-1] == '\r') || (bufChar[bufChar.length()-1] == '\n')) )
				bufChar[bufChar.length()-1] = '\0';
			
			
			//std::cerr << "Nuova lettura: seq:" << bufChar << ". Length " << bufChar.length() << " the last char is "<< (char)bufChar[bufChar.length()-1] << std::endl;
					
			if( bufChar[0] == '>' ) {     //it is the title of fasta file
				nSeq++;
				//std::cerr << "--Length Sequence N. " << (int)(nSeq) << " is " << (int)charsNumber << std::endl;
				if (charsNumber != 0) {
					//charsNumber--;
					if (lengthRead != charsNumber) {
						if (lengthRead < charsNumber)
							lengthRead = charsNumber;
						if (nSeq > 2)
							lenSeq = true;
					}
					
					//dataTypeNSeq numchar = fwrite (&charsNumber, sizeof(dataTypelenSeq), 1 , OutFileLen);
					//assert( numchar == 1); // we should always read the same number of characters
					lengthSeqVector.push_back (charsNumber);
					lengthTexts += charsNumber;
					//if ((verboseEncode==1) || (verboseDistance ==1))
						//std::cerr << "Length Sequence N. " << (int)(nSeq-1) << " is " << (int)charsNumber << " nSeq= "<< nSeq << std::endl;
				}
									
				charsNumber = 0;
			}
			else   //no title
			{

				#if BCR_FROMCYC==1
					for (dataTypelenSeq z = 0 ; z < strlen(bufChar); z++)
						freq[(unsigned int)(bufChar[z])]=1;
				#endif

				// increase the counter of chars
				charsNumber = charsNumber + bufChar.length();
				//std::cerr << "Stessa sequenza --Length Sequence N. " << (int)(nSeq) << " is " << (int)charsNumber << std::endl;
			}
				
		}   //end-while
		//Update the max length of the reads
		if (lengthRead != charsNumber) {
			if (lengthRead < charsNumber)
				lengthRead = charsNumber;
			lenSeq = true;
		}
		lengthSeqVector.push_back (charsNumber);
		lengthTexts += charsNumber;
		infile.close();
	//charsNumber--;
	
	#else
		//kseq 2019-11-08
		gzFile fp;
		kseq_t *seq;
		//int l;
		//cout << input.c_str() << "\n";
		#if  STORE_TITLE_FASTQ == 1
		    char *fnAux = new char[fileOutput.length()+10];
		    sprintf (fnAux,"%s%s",fileOutput.c_str(),".title\0");
		    FILE* OutFile = fopen(fnAux, "wb");
		    if (OutFile==NULL) {
			std::cerr << "TransposeFasta: (title) Error opening " << fnAux << std::endl;
			exit (EXIT_FAILURE);
		    }
		#endif
	
		fp = gzopen(input.c_str(), "r");
		seq = kseq_init(fp);
		while ((kseq_read(seq)) >= 0) {
			nSeq++;
			//printf("name: %s\n", seq->name.s);
			#if  STORE_TITLE_FASTQ == 1
                		fprintf(OutFile, "%s\n", seq->name.s);
            		#endif
			//if (seq->comment.l) 
				//printf("comment: %s\n", seq->comment.s);
			
			//printf("seq: %s\n", seq->seq.s);
			
			charsNumber = strlen(seq->seq.s);
			if (lengthRead != charsNumber) {
				if (lengthRead < charsNumber)
					lengthRead = charsNumber;
				if (nSeq > 2)
					lenSeq = true;
			}
			lengthSeqVector.push_back (charsNumber);
			lengthTexts += charsNumber;
					
			//if (seq->qual.l) 
			//	printf("qual: %s\n", seq->qual.s);
		} //end-while
	
		#if  STORE_TITLE_FASTQ == 1
            		fclose(OutFile);
        	#endif
	
		kseq_destroy(seq);
		gzclose(fp);
	#endif
	
	
	#if verboseEncode==1
	   	dataTypeNChar sum=0;
		for (dataTypeNSeq m=0; m < lengthSeqVector.size(); m++) {
           		sum += lengthSeqVector[m];
       }
	#endif


	if (lenSeq == true)
		cerr << "The  (new and-or old) reads have a different length." << endl;
	else
        cerr << "The (new and-or old) reads have the same length." << endl;

	static FILE *OutFileLen;                  // output file of the end positions;
	OutFileLen = fopen(fileLen, "wb");
	if (OutFileLen==NULL) {
			std::cerr << "TransposeFasta::findLengthNseq: could not open file \"" << fileLen << "\"!"<< std::endl;
			exit (EXIT_FAILURE);
	}

	//Store vector in file
	dataTypeNSeq numchar = fwrite (&lengthSeqVector[0], sizeof(dataTypelenSeq), lengthSeqVector.size(), OutFileLen);
	
//	std::cerr << "lengthTexts= " << (lengthTexts) << " lengthRead= " << (int)lengthRead << std::endl;
//	std::cerr << "lengthSeqVector.size= " << (lengthSeqVector.size()) << " nSeq= " << nSeq << std::endl;
	assert (lengthSeqVector.size() == numchar);
	assert (lengthSeqVector.size() == nSeq);
	
	lengthSeqVector.clear();
	//lengthSeqVector.shrink_to_fit();
			
	#if ( (BUILD_BCR_FROM_BCRpartials==1) && (STORE_LENGTH_IN_FILE==1)  )   //We have to add the length of the sequences from BCR partial files
		char *fnOutLenOld = new char[BCRprefPrev.length()+100];
		sprintf(fnOutLenOld, "%s.len", BCRprefPrev.c_str());
		FILE* OutFileLenOld = fopen(fnOutLenOld, "rb");
		if (OutFileLenOld==NULL) {
			std::cerr << "TransposeFasta::findLengthNseq (BUILD_BCR_FROM_BCRpartials): Error " << fnOutLenOld << "\n";
			exit (EXIT_FAILURE);
		}
		numchar =0;
		dataTypeNChar p;
		numchar = fread (&p, sizeof(dataTypelenSeq), 1 , OutFileLenOld);
		if (numchar > 0) {
			fwrite (&p, sizeof(dataTypelenSeq), numchar , OutFileLen);
			nAddedTextEGSA_transp=1;
		}
		while (numchar != 0) {
			numchar = fread (&p, sizeof(dataTypelenSeq), 1 , OutFileLenOld);
			if (numchar > 0) {
				fwrite (&p, sizeof(dataTypelenSeq), numchar , OutFileLen);
				nAddedTextEGSA_transp++;
			}
		}
		delete [] fnOutLenOld;
		fclose(OutFileLenOld);
	#endif

	delete [] fileLen;
	fclose(OutFileLen);
	
	cerr << "Number of sequences reading/writing: " << nSeq << "\n";
	cerr << "Number of characters reading/writing: " << lengthTexts << "\n";
	
	return true;
}



bool TransposeFasta::convert( const string& input, char const * fileOutput, const string& output, dataTypeNChar ram, string BCRprefPrev )
{

	int resu=true;
	resu=findLengthNseq(input, fileOutput,BCRprefPrev);
	if (resu == false) {  //Error in the reading
		std::cerr << "Error in transpose findLengthNseq! \n";
		exit (EXIT_FAILURE);
	}
	
	#if (BUILD_BCR_FROM_BCRpartials == 1)
        cerr << "In the new collection (excluding the previous BCR files): "<< endl;
    #else
        cerr << "In the new collection, we have: "<< endl;
    #endif
	
	for (dataTypedimAlpha z = 0 ; z < SIZE_ALPHA-1; z++)
		freq[z]=0;
	freq[SIZE_ALPHA-1]=0;
	freq[(unsigned int)(TERMINATE_CHAR)]=nSeq;

	//FILE* ifile;
	//ifile = fopen(input.c_str(), "rb");
	//if( ifile == NULL ) {
	//	cerr << "TrasposeFasta convert: could not open file " << input << " !" << endl;
	//}

	//if (lengthRead > BUFFER_SIZELEN) {
	//	cerr << "Warning: length Read is: " << (int)lengthRead << " and BUFFER_SIZELEN (in Tool.h) is " << BUFFER_SIZELEN << endl;
	//}
	

	vector <FILE*> outputFiles_;
	outputFiles_.resize(lengthRead);    //One for each symbol of the read.
    // create output files   (cyc files)
    for(dataTypelenSeq i=0;i<lengthRead;i++ )
    {
        std::stringstream fn;
        fn << output <<  (int)i << ".txt";
        outputFiles_[i] = fopen( fn.str().c_str(),"w" );
        if (outputFiles_[i] == NULL) {
                std::cerr << "TrasposeFasta: could not open file "  <<  fn.str().c_str() << std::endl;
								exit (EXIT_FAILURE);
				}
        fclose(outputFiles_[i]);
    }

	#if  (USE_QS==1) 
	    vector <FILE*> outputFilesQS_;
	    outputFilesQS_.resize(lengthRead);    //One for each symbol of the read.
	    // create output files   (cyc qs files)
	    for(dataTypelenSeq i=0;i<lengthRead;i++ )
	    {
		std::stringstream fnQS;
		fnQS << output << "qs." << (int)i << ".txt";
		outputFilesQS_[i] = fopen( fnQS.str().c_str(),"w" );
		if (outputFilesQS_[i] == NULL) {
			std::cerr << "TrasposeFasta: could not open file "  <<  fnQS.str().c_str() << std::endl;
					exit (EXIT_FAILURE);
			}
		fclose(outputFilesQS_[i]);
	    }
	#endif

	/*
	ram = ram * nSeq;   //ram used in BCR for vector 
	#if ( (USE_QS==1) && (FASTQ==1) )
		SIZEBUFFERcycFiles = ram / lengthRead / 2;    //I have two buffers at the same time
	#else
		SIZEBUFFERcycFiles = ram / lengthRead;
	#endif 
	
	if ( SIZEBUFFERcycFiles > nSeq )   //number of columns
		SIZEBUFFERcycFiles = nSeq;
	else if ( SIZEBUFFERcycFiles < 1048576 )
		SIZEBUFFERcycFiles = 1048576;
	//cerr << "Size of buffer in TrasposeFasta.convert: " <<  SIZEBUFFERcycFiles << " bytes" << endl;
	//char buf[BUFFER_SIZELEN];
    //for(dataTypeNChar i=0;i<BUFFER_SIZELEN;i++ )
	//	buf[i] = '\0';
	*/
	SIZEBUFFERcycFiles = 1048576;
	
	
	//TO DO: CHECK THE CASE OF charsBuffered >= SIZEBUFFERcycFiles
	//if (nSeq > SIZEBUFFERcycFiles) { //Number of columns of the matrix buffer
	//	cerr << "Warning: Number of sequences is: " << (unsigned long)nSeq << " and SIZEBUFFERcycFiles (in TransposeFasta.h) is " << SIZEBUFFERcycFiles << endl;
	//}
	std::cerr << "TrasposeFasta: init buf_ for bases of size " << (unsigned long) lengthRead << " * " << SIZEBUFFERcycFiles << std::endl;
	#if (USE_QS==1)
		std::cerr << "TrasposeFasta: init buf_ for qs of size " << (unsigned long) lengthRead << " * " << SIZEBUFFERcycFiles << std::endl;
	#endif 
	
	vector<vector<uchar> > buf_;
	buf_.resize(lengthRead);    //For each symbol/column of the read
	for (dataTypelenSeq x = 0 ; x < lengthRead; x++)         //For each symbol/column of the read
		buf_[x].resize(SIZEBUFFERcycFiles);
	for (dataTypelenSeq x = 0 ; x < lengthRead; x++)         //For each symbol/column of the read
		for (dataTypeNChar y = 0 ; y < SIZEBUFFERcycFiles; y++)         //For each buffered symbol of the read
			buf_[x][y]=TERMINATE_CHAR_LEN;

	#if  (USE_QS==1)
		vector<vector<uchar> > bufQS_;
		bufQS_.resize(lengthRead);    //For each symbol/column of the read
		for (dataTypelenSeq x = 0 ; x < lengthRead; x++)         //For each symbol/column of the read
			bufQS_[x].resize(SIZEBUFFERcycFiles);
		for (dataTypelenSeq x = 0 ; x < lengthRead; x++)         //For each symbol/column of the read
			for (dataTypeNChar y = 0 ; y < SIZEBUFFERcycFiles; y++)         //For each buffered symbol of the read
				bufQS_[x][y]=TERMINATE_CHAR_LEN;
	#endif
	
	nSeq = 0;
	dataTypeNChar num_write = 0;
	dataTypeNChar charsBuffered = 0;
	
	#if KSEQ_PARSER==0
		std::ifstream infile(input.c_str());
		string bufChar;
		dataTypelenSeq tmpLen = 0;
		dataTypelenSeq sumLenCum = 0;	
		
		// looping through the input file, add the characters to the buffer, print buffer when it's full
	//    unsigned int num_read = 0;
		
		

		while (getline(infile, bufChar))  {
			if ((bufChar[bufChar.length()-2] == '\r') || (bufChar[bufChar.length()-2] == '\n'))
				bufChar[bufChar.length()-2] = '\0';
			else if ((bufChar[bufChar.length()-1] == '\r') || (bufChar[bufChar.length()-1] == '\n'))
				bufChar[bufChar.length()-1] = '\0';		
			
			tmpLen = bufChar.length();          //tmpLen = strlen(buf)-1;
			
			if ( (charsBuffered > 0) && ( charsBuffered-1 == SIZEBUFFERcycFiles))   //it is linked to the number of sequences
			{
				// write buffers to the files, clear buffers
				for(dataTypelenSeq i=0;i<lengthRead;i++ )  {       //For each symbol/column of the read   (for each string)
					std::stringstream fn;
					fn << output <<  (int)i << ".txt";
					outputFiles_[i] = fopen( fn.str().c_str(),"a" );
					num_write = fwrite ( &buf_[i][0],sizeof(char),charsBuffered-1,outputFiles_[i] );
					assert( num_write == charsBuffered-1 );
					for(dataTypeNChar x=0; x<charsBuffered-1; x++ ) {
						//For each buffered symbol of the read (a symbol for each string)
						if (buf_[i][x] != TERMINATE_CHAR_LEN) {
							freq[(unsigned int)(buf_[i][x])]++;   //=1
							buf_[i][x] = TERMINATE_CHAR_LEN;
						}
					}
					fclose(outputFiles_[i]);
				}
				
				charsBuffered=1;
			}

			// process the input
			if( bufChar[0] != '>' )   //no title
			{
				
				dataTypelenSeq posit = 0 + sumLenCum;
				for(dataTypelenSeq i=posit; i < posit + tmpLen;i++) {
					buf_[i][charsBuffered-1] = bufChar[i-posit];
				}

				sumLenCum += tmpLen;
				//cerr << "it is not title"<< " tmpLen " << (int)tmpLen << " sumLenCum " << (unsigned int)sumLenCum << endl;
				
				/*#if verboseEncode==1
					cerr << "Partial Buf_ " << endl;
					for(dataTypelenSeq i=0;i<lengthRead;i++ ) {
						cerr << (int)i << " " ;
						for(dataTypeNChar x=0; x<charsBuffered; x++ )
							cerr << buf_[i][x] ;
						cerr << endl;
					}
				#endif */
			}
			else  {        //it is a title
				sumLenCum = 0;
				// increase the number of sequences
				nSeq++;
				
				// increase the counter of chars buffered
				charsBuffered++;
			}
			//for(dataTypelenSeq i=0;i<lengthRead;i++ )
			//	buf[i] = '\0';
			//fgets ( buf, BUFFER_SIZELEN, ifile );
			//if ((buf[strlen(buf)-2] == '\r') || (buf[strlen(buf)-2] == '\n'))
			 //   buf[strlen(buf)-2] = '\0';
			//else if ((buf[strlen(buf)-1] == '\r') || (buf[strlen(buf)-1] == '\n'))
			//    buf[strlen(buf)-1] = '\0';

	//        cerr << "buf: " << buf << "." << endl;
			
		}  //end-while
		infile.close();
		/* 	#if verboseEncode==1
		cerr << "The last Buf_ " << endl;
		cerr << "Buf_ " << endl;
		for(dataTypelenSeq i=0;i<lengthRead;i++ ) {
			for(dataTypeNChar x=0; x<charsBuffered; x++ )
				cerr << buf_[i][x] ;
			cerr << endl;
		} */
		// write the rest
		for(dataTypelenSeq i=0;i<lengthRead;i++ )
		{
	//        num_write = fwrite ( buf_[i],sizeof(uchar),charsBuffered,outputFiles_[i] );
			std::stringstream fn;
			fn << output <<  (int)i << ".txt";
			outputFiles_[i] = fopen( fn.str().c_str(),"a" );
			num_write = fwrite ( &buf_[i][0],sizeof(char),charsBuffered,outputFiles_[i] );
			assert( num_write == charsBuffered );
			for(dataTypeNChar x=0; x<charsBuffered; x++ ) {
				if (buf_[i][x] != TERMINATE_CHAR_LEN) {
					freq[(unsigned int)(buf_[i][x])]++;  //=1
					buf_[i][x] = TERMINATE_CHAR_LEN;
				}
				//cerr << "Number of characters reading/writing: " << (int) lengthTexts << "\n";
			}
			fclose( outputFiles_[i]);
		}
		
	#else    //KSEQ_PARSER==1
		//kseq 2019-11-08
		gzFile fp;
		kseq_t *seq;
		//int l;
		//cout << input.c_str() << "\n";
		fp = gzopen(input.c_str(), "r");
		seq = kseq_init(fp);
		charsBuffered=1;
		while ((kseq_read(seq)) >= 0) {
			nSeq++;
			//printf("name: %s\n", seq->name.s);
			//if (seq->comment.l) 
				//printf("comment: %s\n", seq->comment.s);
			
			//printf("seq -: %s\n", seq->seq.s);
			
			//If the buf is full, then clear the buffer and initialize
			if ( (charsBuffered > 0) && ( charsBuffered-1 == SIZEBUFFERcycFiles))   
				//it is linked to the number of sequences
			{
				// Only if buffer is full, then write buffer to the files, clear buffers
				for(dataTypelenSeq i=0;i<lengthRead;i++ )  {       //For each symbol/column of the read   (for each string)
					std::stringstream fn;
					fn << output <<  (int)i << ".txt";
					outputFiles_[i] = fopen( fn.str().c_str(),"a" );
					num_write = fwrite ( &buf_[i][0],sizeof(char),charsBuffered-1,outputFiles_[i] );
					assert( num_write == charsBuffered-1 );
					for(dataTypeNChar x=0; x<charsBuffered-1; x++ ) {
						//For each buffered symbol of the read (a symbol for each string)
						
						if (buf_[i][x] != TERMINATE_CHAR_LEN) {
							freq[(unsigned int)(buf_[i][x])]++;   //=1
							buf_[i][x] = TERMINATE_CHAR_LEN;							
						}
					}
					fclose(outputFiles_[i]);
				}
				
				#if  (USE_QS==1)
					// write buffers to the files, clear buffers
					for(dataTypelenSeq i=0;i<lengthRead;i++ )  {       //For each symbol/column of the read   (for each string)
						std::stringstream fnQS;
						fnQS << output << "qs." << (int)i << ".txt";
						outputFilesQS_[i] = fopen( fnQS.str().c_str(),"a" );
						num_write = fwrite ( &bufQS_[i][0],sizeof(char),charsBuffered-1,outputFilesQS_[i] );
						assert( num_write == charsBuffered-1 );
						for(dataTypeNChar x=0; x<charsBuffered-1; x++ ) {
							//For each buffered symbol of the read (a symbol for each string)
							if (bufQS_[i][x] != TERMINATE_CHAR_LEN) {
								bufQS_[i][x] = TERMINATE_CHAR_LEN;
							}
						}
						fclose(outputFilesQS_[i]);
					}
				#endif
				
				charsBuffered=1;
			}
			
			
			#if BCR_SET_ALN_RH==0   //align left
				dataTypelenSeq i;
				for(i=0; i < strlen(seq->seq.s);i++) {				
					buf_[i][charsBuffered-1] = seq->seq.s[i];				
				}			
				
				//if (seq->qual.l) 
				//	printf("qual: %s\n", seq->qual.s);
				#if USE_QS==1
				for(dataTypelenSeq i=0; i < strlen(seq->qual.s);i++) {				
					bufQS_[i][charsBuffered-1] = seq->qual.s[i];				
				}
				#endif
				
			#else 
				//align right
				dataTypelenSeq index=lengthRead-strlen(seq->seq.s);
				if(index>0)	buf_[index-1][charsBuffered-1] = TERMINATE_CHAR;
				for(dataTypelenSeq i=0; i < strlen(seq->seq.s); i++) {				
					buf_[index+i][charsBuffered-1] = seq->seq.s[i];	
				}
				
				#if USE_QS==1
				for(dataTypelenSeq i=0; i < strlen(seq->qual.s);i++) {				
					bufQS_[index+i][charsBuffered-1] = seq->qual.s[i];
				}
				#endif
				
			#endif
			
			/*
			for(dataTypelenSeq i=0;i<lengthRead;i++ ) {
				cerr << "Buf_ " << (int)i << endl;
				for(dataTypeNChar x=0; x<charsBuffered; x++ )
					cerr << buf_[i][x] ;
				cerr << endl;
			}
			*/
			
			// increase the counter of chars buffered
			charsBuffered++;
					
			
		} //end-while
		kseq_destroy(seq);
		gzclose(fp);
		
		#if verboseEncode==1
			cerr << "The last Buf_ " << endl;
			cerr << "Buf_ " << endl;
			for(dataTypelenSeq i=0;i<lengthRead;i++ ) {
				for(dataTypeNChar x=0; x<charsBuffered-1; x++ )
					cerr << buf_[i][x];
				cerr << endl;
			}
		#endif 
		
		
		// write the rest
		charsBuffered--;
		//cerr << "*charsBuffered " << charsBuffered << "\n";
		for(dataTypelenSeq i=0;i<lengthRead;i++ )
		{
			std::stringstream fn;
			fn << output <<  (int)i << ".txt";
			outputFiles_[i] = fopen( fn.str().c_str(),"a" );
			num_write = fwrite ( &buf_[i][0],sizeof(char),charsBuffered,outputFiles_[i] );
			assert( num_write == charsBuffered );
			
			for(dataTypeNChar x=0; x<charsBuffered; x++ ) {
				if (buf_[i][x] != TERMINATE_CHAR_LEN) {
					freq[(unsigned int)(buf_[i][x])]++;  //=1
					//cerr << "*buf_["<< (int)i<<"]["<< (int)x << "]= " << (int)buf_[i][x] << ". " << endl;
					buf_[i][x] = TERMINATE_CHAR_LEN;
				}
			}
			fclose( outputFiles_[i]);
		}
		
		#if USE_QS==1
			for(dataTypelenSeq i=0;i<lengthRead;i++ )
			{
				std::stringstream fnQS;
				fnQS << output << "qs." << (int)i << ".txt";
				outputFilesQS_[i] = fopen( fnQS.str().c_str(),"a" );
				num_write = fwrite ( &bufQS_[i][0],sizeof(char),charsBuffered,outputFilesQS_[i] );
				assert( num_write == charsBuffered );
			
				for(dataTypeNChar x=0; x<charsBuffered; x++ ) {
					if (bufQS_[i][x] != TERMINATE_CHAR_LEN) {
						//cerr << "*bufQS_["<< (int)i<<"]["<< (int)x << "]= " << (int)bufQS_[i][x] << ". " << endl;
						bufQS_[i][x] = TERMINATE_CHAR_LEN;
					}
				}
				fclose( outputFilesQS_[i]);
			}
		#endif
		
		
	#endif   //end KSEQ_PARSER

	//Free memory
	//Requests the removal of unused capacity
	for (dataTypelenSeq x = 0 ; x < lengthRead; x++)  {
		buf_[x].clear();
		buf_[x].shrink_to_fit();
		#if  (USE_QS==1)
			bufQS_[x].clear();
			bufQS_[x].shrink_to_fit();
		#endif
	}
	buf_.clear();
	buf_.shrink_to_fit();	
	outputFiles_.clear();
	outputFiles_.shrink_to_fit();
	
	#if  (USE_QS==1)
		bufQS_.clear();
		bufQS_.shrink_to_fit();
		outputFilesQS_.clear();
	#endif
			
	sizeAlpha=0;
	for (dataTypedimAlpha i = 0; i < SIZE_ALPHA-1; ++i)
		if (freq[i] > 0) {
			sizeAlpha++;
		}
	if (freq[SIZE_ALPHA-1] > 0) {		
		sizeAlpha++;
	}
	
	char *fnAux = new char[strlen(fileOutput)+10];
	sprintf (fnAux,"%s%s",fileOutput,".info\0");
	FILE* OutFile = fopen(fnAux, "wb");
	if (OutFile==NULL) {
		std::cerr << "TransposeFasta: (lengthBWT+NSequences+sizeAlpha) Error opening " << fnAux << std::endl;
		exit (EXIT_FAILURE);
	}
	
	dataTypeNChar lenTot=lengthTexts+nSeq;
	fwrite(&lenTot,sizeof(dataTypeNChar),1,OutFile);
	fwrite(&nSeq,sizeof(dataTypeNSeq),1,OutFile);
	fwrite(&sizeAlpha,sizeof(dataTypedimAlpha),1,OutFile);
	fclose(OutFile);
	delete [] fnAux;
	
	cerr << "TransposeFasta: The max length (Read) is: " << (int)lengthRead << endl;
	cerr << "TransposeFasta: Number of reads: " << nSeq << endl;
	cerr << "TransposeFasta: Total Number of chars (without end-markers): " << 	lengthTexts << endl;	
	cerr << "TransposeFasta: Size Alpha: " << 	(int)sizeAlpha << " symbols" << endl;
		
    return true;
}



bool TransposeFasta::convertFromCycFile(const string& input, char const * fileOutput, string BCRprefPrev) {

	int res=true;
	res=findLengthNseq(input, fileOutput, BCRprefPrev);
	if (res == false) {  //Error in the reading
		std::cerr << "Error in transpose findLengthNseq! \n";
		exit (EXIT_FAILURE);
	}

  	//TO DO
	//The distribution of characters is useful
	//for alpha[SIZE_ALPHA] -->Corresponding between the alphabet, the piles and tableOcc
	//and to know sizeAlpha

	lengthTexts = lengthRead * nSeq;
	
		dataTypedimAlpha sizeAlpha=0;
	for (dataTypedimAlpha i = 0; i < SIZE_ALPHA-1; ++i)
		if (freq[i] > 0) {
			sizeAlpha++;
		}
	if (freq[SIZE_ALPHA-1] > 0) {		
		sizeAlpha++;
	}
	
	char *fnAux = new char[strlen(fileOutput)+10];
	sprintf (fnAux,"%s%s",fileOutput,".info\0");
	FILE* OutFile = fopen(fnAux, "wb");
	if (OutFile==NULL) {
		std::cerr << "TransposeFasta: (lengthBWT+NSequences+sizeAlpha) Error opening " << fnAux << std::endl;
		exit (EXIT_FAILURE);
	}
	
	dataTypeNChar lenTot=lengthTexts+nSeq;
	fwrite(&lenTot,sizeof(dataTypeNChar),1,OutFile);
	fwrite(&nSeq,sizeof(dataTypeNSeq),1,OutFile);
	fwrite(&sizeAlpha,sizeof(dataTypedimAlpha),1,OutFile);
	fclose(OutFile);
	delete [] fnAux;
	
	cerr << "TransposeFasta: The max length (Read) is: " << (int)lengthRead << endl;
	cerr << "TransposeFasta: Number of reads: " << nSeq << endl;
	cerr << "TransposeFasta: Number of chars: " << 	lengthTexts << endl;	
	cerr << "TransposeFasta: Size Alpha: " << 	(int)sizeAlpha << endl;
	
	
	
	
	
    return 1;
}


// like fgets, but return number of chars
ulong TransposeFasta::readln(char* s, int n, FILE* iop) {

  register long c=0;
  register char* cs;
  register long z = 0; // count characters
  cs=s;
  while (--n > 0 && (c = getc(iop)) != EOF) {
	  //std::cerr << "n= " << n << " c= " << c << " ";
	  //printf ("%c\n", c);
	  if ((*cs++ = c) == '\n')
		break;
	z++;
  }
  if (c == EOF && z==0)
	  z = 0; // signify end of file
  *cs = '\0';
  return z;
}




bool TransposeFasta::convertQS( const string& input, const string& output, dataTypeNChar ram )
{
	FILE* ifile;
	cerr << "\nTrasposeFasta convertQS: file " << input << " !" << endl;
	
	ifile = fopen(input.c_str(), "rb");

    if( ifile == NULL ) { 
		cerr << "TrasposeFasta convertQS: could not open file " << input << " !" << endl; 
	}


	vector <FILE*> outputFiles_;
	outputFiles_.resize(lengthRead);    //One for each symbol of the read.
    // create output files   (cyc qs files)
    for(dataTypelenSeq i=0;i<lengthRead;i++ )
    {
        std::stringstream fn;
        fn << output << "qs." << (int)i << ".txt";
        outputFiles_[i] = fopen( fn.str().c_str(),"w" );
        if (outputFiles_[i] == NULL) {
                std::cerr << "TrasposeFasta: could not open file "  <<  fn.str().c_str() << std::endl;
				exit (EXIT_FAILURE);
		}
        fclose(outputFiles_[i]);
    }

	//TO DO: CHECK THE CASE OF charsBuffered >= SIZEBUFFERcycFiles
	//if (nSeq > SIZEBUFFERcycFiles) {
	//	cerr << "Warning: Number of sequences is: " << (unsigned long)nSeq << " and SIZEBUFFERcycFiles (in TransposeFasta.h) is " << SIZEBUFFERcycFiles << endl;
	//}
	
	ram = ram * nSeq;   //ram used in BCR for vector 
	
	SIZEBUFFERcycFiles = ram / lengthRead;
	
	if ( SIZEBUFFERcycFiles > nSeq )   //number of columns
		SIZEBUFFERcycFiles = nSeq;
	else if ( SIZEBUFFERcycFiles < 1048576 )
		SIZEBUFFERcycFiles = 1048576;
	
	
	std::cerr << "TrasposeFasta.convertQS: init buf_ of size " << (unsigned long) lengthRead << " * " << SIZEBUFFERcycFiles << std::endl;
	
	vector<vector<uchar> > buf_;
	buf_.resize(lengthRead);    //For each symbol/column of the read
	for (dataTypelenSeq x = 0 ; x < lengthRead; x++)         //For each symbol/column of the read
		buf_[x].resize(SIZEBUFFERcycFiles);
	for (dataTypelenSeq x = 0 ; x < lengthRead; x++)         //For each symbol/column of the read
		for (dataTypeNChar y = 0 ; y < SIZEBUFFERcycFiles; y++)         //For each buffered symbol of the read
			buf_[x][y]=TERMINATE_CHAR_LEN;

	//char buf[SIZEBUFFERcycFiles];
    //for(dataTypeNChar i=0;i<SIZEBUFFERcycFiles;i++ )
//		buf[i] = '\0';	
	//fgets ( buf, SIZEBUFFERcycFiles, ifile );			//it starts with '>', we can ignore it
	
	std::ifstream infile(input.c_str());
	string bufChar;
	
	dataTypeNChar charsBuffered = 0;
    // looping through the input file, add the characters to the buffer, print buffer when it's full
//    unsigned int num_read = 0;
    dataTypeNChar num_write = 0;
	dataTypelenSeq tmpLen = 0;
    nSeq = 0;
	dataTypelenSeq sumLenCum = 0;

    /*while( !feof(ifile) )
    {
        if( charsBuffered-1 == SIZEBUFFERcycFiles)   //it is linked to the number of sequences
        {
            // write buffers to the files, clear buffers
            for(dataTypelenSeq i=0;i<lengthRead;i++ )  {       //For each symbol/column of the read   (for each string)
                std::stringstream fn;
                fn << output << "qs." << (int)i << ".txt";
				outputFiles_[i] = fopen( fn.str().c_str(),"a" );
				num_write = fwrite ( &buf_[i][0],sizeof(char),charsBuffered-1,outputFiles_[i] );
				assert( num_write == charsBuffered-1 );
				for(dataTypeNChar x=0; x<charsBuffered-1; x++ ) {
					//For each buffered symbol of the read (a symbol for each string)
					if (buf_[i][x] != TERMINATE_CHAR_LEN) {
						//freq[(unsigned int)(buf_[i][x])]=1;
						buf_[i][x] = TERMINATE_CHAR_LEN;
					}
				}
				fclose(outputFiles_[i]);
            }
			
	        //for (dataTypelenSeq x = 0 ; x < lengthRead; x++)         //For each symbol/column of the read
		    //    for (dataTypeNChar y = 0 ; y < BUFFERSIZE; y++)         //For each buffered symbol of the read
			//        buf_[x][y]=TERMINATE_CHAR_LEN;
			
            charsBuffered=1;
        }

        // process the input
        if( buf[0] != '>' )
        {
			dataTypelenSeq posit = 0 + sumLenCum;
			for(dataTypelenSeq i=posit; i < posit + tmpLen;i++) {
				buf_[i][charsBuffered-1] = buf[i-posit];
			}

            sumLenCum += tmpLen;

			#if verboseEncode==1
				cerr << "Partial Buf_ " << endl;
				for(dataTypelenSeq i=0;i<lengthRead;i++ ) {
					cerr << (int)i << " " ;
					for(dataTypeNChar x=0; x<charsBuffered; x++ )
						cerr << buf_[i][x] ;
					cerr << endl;
				}
			#endif
        }
		else  {
			if (nSeq % SIZEBUFFERcycFiles == 0)
				cerr << "New QS sequence - " << "Sequence N. " << (int)nSeq << std::endl;

			sumLenCum = 0;
		    // increase the number of sequences
			nSeq++;
			// increase the counter of chars buffered
			charsBuffered++;
		}
		for(dataTypelenSeq i=0;i<lengthRead;i++ )
			buf[i] = '\0';
		fgets ( buf, SIZEBUFFERcycFiles, ifile );
        if ((buf[strlen(buf)-2] == '\r') || (buf[strlen(buf)-2] == '\n'))
            buf[strlen(buf)-2] = '\0';
        else if ((buf[strlen(buf)-1] == '\r') || (buf[strlen(buf)-1] == '\n'))
            buf[strlen(buf)-1] = '\0';

		tmpLen = strlen(buf);          //tmpLen = strlen(buf)-1;
    }
	*/
	
	while (getline(infile, bufChar))  {
		if ((bufChar[bufChar.length()-2] == '\r') || (bufChar[bufChar.length()-2] == '\n'))
            bufChar[bufChar.length()-2] = '\0';
        else if ((bufChar[bufChar.length()-1] == '\r') || (bufChar[bufChar.length()-1] == '\n'))
            bufChar[bufChar.length()-1] = '\0';		
		
		tmpLen = bufChar.length();          //tmpLen = strlen(buf)-1;
		
		if ( (charsBuffered > 0) && ( charsBuffered-1 == SIZEBUFFERcycFiles))   //it is linked to the number of sequences
        {
			// write buffers to the files, clear buffers
            for(dataTypelenSeq i=0;i<lengthRead;i++ )  {       //For each symbol/column of the read   (for each string)
                std::stringstream fn;
                fn << output << "qs." << (int)i << ".txt";
				outputFiles_[i] = fopen( fn.str().c_str(),"a" );
				num_write = fwrite ( &buf_[i][0],sizeof(char),charsBuffered-1,outputFiles_[i] );
				assert( num_write == charsBuffered-1 );
				for(dataTypeNChar x=0; x<charsBuffered-1; x++ ) {
					//For each buffered symbol of the read (a symbol for each string)
					if (buf_[i][x] != TERMINATE_CHAR_LEN) {
						//freq[(unsigned int)(buf_[i][x])]++;   //=1
						buf_[i][x] = TERMINATE_CHAR_LEN;
					}
				}
				fclose(outputFiles_[i]);
            }
            charsBuffered=1;
        }

        // process the input
        if( bufChar[0] != '{' )   //no title
        {
			
			dataTypelenSeq posit = 0 + sumLenCum;
			for(dataTypelenSeq i=posit; i < posit + tmpLen;i++) {
				buf_[i][charsBuffered-1] = bufChar[i-posit];
			}

            sumLenCum += tmpLen;
			//cerr << "it is not title"<< " tmpLen " << (int)tmpLen << " sumLenCum " << (unsigned int)sumLenCum << endl;
			
 			/*#if verboseEncode==1
				cerr << "Partial Buf_ " << endl;
				for(dataTypelenSeq i=0;i<lengthRead;i++ ) {
					cerr << (int)i << " " ;
					for(dataTypeNChar x=0; x<charsBuffered; x++ )
						cerr << buf_[i][x] ;
					cerr << endl;
				}
			#endif */
        }
		else  {        //it is a title
            sumLenCum = 0;
		    // increase the number of sequences
			nSeq++;
			//cerr << "it is a title"<< nSeq << endl;
			// increase the counter of chars buffered
			charsBuffered++;
		}
		//for(dataTypelenSeq i=0;i<lengthRead;i++ )
		//	buf[i] = '\0';
		//fgets ( buf, BUFFER_SIZELEN, ifile );
        //if ((buf[strlen(buf)-2] == '\r') || (buf[strlen(buf)-2] == '\n'))
         //   buf[strlen(buf)-2] = '\0';
        //else if ((buf[strlen(buf)-1] == '\r') || (buf[strlen(buf)-1] == '\n'))
        //    buf[strlen(buf)-1] = '\0';

//        cerr << "buf: " << buf << "." << endl;
		
    }  //end-while
	
	

	/* 	#if verboseEncode==1
		cerr << "The last Buf_ " << endl;
		cerr << "Buf_ " << endl;
		for(dataTypelenSeq i=0;i<lengthRead;i++ ) {
			for(dataTypeNChar x=0; x<charsBuffered; x++ )
				cerr << buf_[i][x] ;
			cerr << endl;
		}
	#endif */

    // write the rest
    for(dataTypelenSeq i=0;i<lengthRead;i++ )
    {
        std::stringstream fn;
        fn << output << "qs." << (int)i << ".txt";
        outputFiles_[i] = fopen( fn.str().c_str(),"a" );
		num_write = fwrite ( &buf_[i][0],sizeof(char),charsBuffered,outputFiles_[i] );
		assert( num_write == charsBuffered );
		for(dataTypeNChar x=0; x<charsBuffered; x++ ) {
			if (buf_[i][x] != TERMINATE_CHAR_LEN) {
				//freq[(unsigned int)(buf_[i][x])]=1;
				buf_[i][x] = TERMINATE_CHAR_LEN;
			}
		}
        fclose( outputFiles_[i]);
    }

	cerr << "Number of quality score sequences reading/writing: " << nSeq << "\n";
	//cerr << "Number of characters reading/writing: " << lengthTexts << "\n";

    return true;
}




bool TransposeFasta::convert1Sequence(char const * filename1, dataTypeNChar ram) {
  		//TO DO
	//lengthRead = CYCLENUM;
    std::cerr << "***TransposeFasta::convert1Sequence "<< std::endl;
	//The distribution of characters is useful
	//for alpha[SIZE_ALPHA] -->Corresponding between the alphabet, the piles and tableOcc

	for (dataTypedimAlpha z = 0 ; z < SIZE_ALPHA-1; z++)
		freq[z]=0;
	freq[SIZE_ALPHA-1]=0;
	freq[(unsigned int)(TERMINATE_CHAR)]=1;

	#if BCR_INPUT_IN_MEMORY==1  	// BCR reads from string
		//if (lengthRead <= 1000000000) {
			std::cerr << "(BCR_INPUT_IN_MEMORY==1) Copy File into String"<< std::endl;
			//std::ifstream inFile;
			//inFile.open(filename1);//open the input file

			//std::stringstream strStream;
			//strStream << inFile.rdbuf();//read the file
			//strInput = strStream.str();//strInput holds the content of the file

			//#if verboseEncode==1
			//	std::cout << "Input:" << strInput << "." << std::endl;
			//#endif
			//lengthRead = strInput.length();  //maximum length of reads length of (sequence without \0)

            FILE *Infile = fopen(filename1, "rb"); // b is for binary: required by DOS
            if(Infile==NULL) {
                fprintf(stderr,"Please provide the input file name (open_files)\n");
                exit(1);
            }
            fseek(Infile,0,SEEK_END);
            lengthRead=ftell(Infile);
            rewind(Infile);
            strInput = new uchar[lengthRead];      // text
            dataTypeNChar num;
            num=fread(strInput, sizeof(uchar), lengthRead, Infile);
            if(num!=lengthRead) {
                fprintf(stderr,"Error reading the input file!\n");
                exit(1);
            }
            #if verboseEncode == 1
                fprintf(stderr,"File size: %d bytes\n",lengthRead);
            #endif
            fclose(Infile);


			nSeq = 1;

			std::cout << "Length Input (lengthRead):" << lengthRead << std::endl;

			for(dataTypeNChar x=0; x<lengthRead; x++ ) {

				 freq[(unsigned int)(strInput[x])]++;

				 if ((int)strInput[x] <  (int)TERMINATE_CHAR) {
					 std::cout << "TERMINATE_CHAR = " << TERMINATE_CHAR << " = " << (int)TERMINATE_CHAR << std::endl;
					 std::cout << "is greater than strInput[ " << x << "]=" << strInput[x] << " = "<< (int)strInput[x] << std::endl;
					 exit(1);
				 }
				 else if ((int)strInput[x] ==  (int)TERMINATE_CHAR) {
					 std::cout << "TERMINATE_CHAR = " << TERMINATE_CHAR << " = " << (int)TERMINATE_CHAR << std::endl;
					 std::cout << "is equal to strInput[ " << x << "]=" << strInput[x] << " = "<< (int)strInput[x] << std::endl;
					 std::cerr << "WARNING: TERMINATE_CHAR= "<< (int)TERMINATE_CHAR << " IS NOT THE UNIQUE SYMBOL IN THE INPUT FILE"<< std::endl;
					// exit(1);
				 }
			}

			lengthTexts = lengthRead;

			// std::reverse(strInput.begin(), strInput.end());
			// #if verboseEncode==1
				// std::cout << "Input reverse:" <<  strInput << std::endl;
			// #endif
		//}
		//else {
		//	std::cerr << "TransposeFasta.cpp - BCR_INPUT_IN_MEMORY==1 -  lengthRead > 1000000000 " << std::endl;
		//	exit (EXIT_FAILURE);
		//}
	#else
		std::cerr << "TransposeFasta.cpp - change BCR_INPUT_IN_MEMORY to 1. External memory for 1 sequence is not implemented!" << std::endl;
	#endif
	 /////////////////////////////////

	 std::cerr << "End traspose\n";

	 return 1;

  /*
  }
  else {
	printf("<--> Error: Could not open %s!\n", filename1);
	return 0;
  }
  */

  std::cout << "Number of sequences reading/writing: " << nSeq << "\n";
  std::cout << "Number of characters reading/writing: " << lengthTexts << "\n";




  std::cerr << "****************End read file "<< std::endl;
  return true;
}
