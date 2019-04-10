// mex MexDataConverter.c COMPFLAGS="$COMPFLAGS -openmp"  LINKFALGS="$LINKFALGS -openmp"

#include "mex.h"
#include <stdio.h>
#include <stdlib.h>

#include "intrin.h"

// System-based h files
#include "../clinicalParameters.h"
#include "../killed_channels.h"


void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] )
{
    ///////////////////////////////////////////////////////////////////////
    // Definitions ////////////////////////////////////////////////////////
    
    // Variables
    int v;
    int Head;
	int rosszadat;
    unsigned short int *Line, *VectorR, *VectorW;
    // Rearrange detector indices
    #if NUMBER_OF_DATA_2BYTES == 36
      int sortDSR0[NUMBER_OF_DATA_2BYTES] = {29,35,17,19,1,7,33,25,21,15,11,5,31,27,23,13,9,3,32,28,24,14,10,6,34,26,22,16,12,4,30,36,18,20,2,8};
      unsigned short int DSR_KEYWORD = 21040;
    #elif NUMBER_OF_DATA_2BYTES == 72
      int sortDSR0[NUMBER_OF_DATA_2BYTES] = {63,71,59,60,72,64,67,55,51,52,56,68,39,47,43,44,48,40,35,27,31,32,28,36,7,19,23,24,20,8,11,15,3,16,4,12,61,69,57,58,70,62,65,53,49,50,54,66,37,45,41,42,46,38,33,25,29,30,26,34,5,17,21,22,18,6,9,13,1,14,2,10}; //olasz clinical ordering
	  // int sortDSR0[NUMBER_OF_DATA_2BYTES] = {53,55,17,19,65,71,29,35,69,61,33,25,37,43,1,7,47,41,11,5,57,51,21,15,59,49,23,13,67,63,31,27,68,64,32,28,45,39,9,3,46,42,10,6,60,50,24,14,58,52,22,16,70,62,34,26,66,72,30,36,48,40,12,4,38,44,2,8,54,56,18,20}; //Zoli sajátja
      unsigned short int DSR_KEYWORD = 21041;
    #endif
      
    // Right-hand values (Inputs of mex function)
    char *pathIn, *pathOut;
    FILE *InputFile, *OutputFile;
        
    // Left-hand values (Outputs of mex function)
    double *Count;
    
    
    ///////////////////////////////////////////////////////////////////////
    // Initializations ////////////////////////////////////////////////////
    
    // Right-hand values
    pathIn = mxArrayToString(prhs[0]);
    pathOut = mxArrayToString(prhs[1]);
    
    // Left-hand values
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    Count = mxGetPr(plhs[0]);
    Count[0]=0;
    
    
    ///////////////////////////////////////////////////////////////////////
    // Memory and file handling ///////////////////////////////////////////
    
    Line = (unsigned short int*) malloc (sizeof(unsigned short int)*4);
    if (Line == NULL) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "Line");
        return;
    }
    
    VectorR = (unsigned short int*) malloc (sizeof(unsigned short int)*(NUMBER_OF_HEADER_2BYTES+NUMBER_OF_DATA_2BYTES));
    if (VectorR == NULL) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "VectorR");
        return;
    }

    VectorW = (unsigned short int*) malloc (sizeof(unsigned short int)*(NUMBER_OF_HEADER_2BYTES+NUMBER_OF_DATA_2BYTES));
    if (VectorR == NULL) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "VectorW");
        return;
    }

    InputFile = fopen ( pathIn , "rb" );
    if (InputFile==NULL) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "InputFile");
        return;
    }
    
    OutputFile = fopen ( pathOut , "wb" );
    if (InputFile==NULL) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "InputFile");
        return;
    }

    
    
    ///////////////////////////////////////////////////////////////////////
    // Data Convert ///////////////////////////////////////////////////////
    
    while( !feof(InputFile) ){
        fread (Line,sizeof(unsigned short int),4,InputFile);
        
        if ( (_byteswap_ushort(Line[2]) == 17491) && (_byteswap_ushort(Line[3]) == DSR_KEYWORD) ) { // DSR0/DSR1
            
            fread(VectorR,sizeof(unsigned short int),NUMBER_OF_HEADER_2BYTES+NUMBER_OF_DATA_2BYTES,InputFile);
            
            Head = (int)(VectorR[0] / 8 - 1);
            //printf("%d-",Head);
            
            VectorW[0] = 17492; // "TD"
            for ( v=1; v<NUMBER_OF_HEADER_2BYTES; v++ )
                VectorW[v] = 0;

            for ( v=NUMBER_OF_HEADER_2BYTES; v<NUMBER_OF_HEADER_2BYTES+NUMBER_OF_DATA_2BYTES; v++ )
                VectorW[v] = _byteswap_ushort(VectorR[sortDSR0[v-4]+4-1]);
            
            // for ( v=40; v<68; v++ )
                 // VectorW[v] = 13;
            
            //-> Ez egy szures, mert az elso 18 es az utolso 18 csatorna ugyanazt a 2^x erteket adta vissza
            // tunet: fenyes pixel lesz
			rosszadat=0;

			for ( v=NUMBER_OF_HEADER_2BYTES; v<NUMBER_OF_HEADER_2BYTES+NUMBER_OF_DATA_2BYTES; v++ )
			{
                
             if (  VectorW[v]<32 ||  VectorW[v]>4080 )
			 {   // ( VectorW[NUMBER_OF_HEADER_2BYTES]!=VectorW[NUMBER_OF_HEADER_2BYTES+1] )
                 //printf("%d",v);
                 
                     if (FLAG_KILL==1)
                         {
                                if (!((Head==NODE_WITH_KILL) && (v==KILLED_CH))) 
                                {
                                        rosszadat=1;
                                 }
                         }
                     else if (FLAG_KILL==2)
                        {
                            if (!(((Head==NODE_WITH_KILL) && (v==KILLED_CH)) || (Head==NODE_WITH_KILL2) && (v==KILLED_CH2))) 
                                {
                                        rosszadat=1;
                                 }
                         }
                     else 
                        {
                            rosszadat=1;
                        }
                
             }
            }
			if (rosszadat!=1)
			{
				Count[0]+=1;
                fwrite(VectorW,sizeof(unsigned short int),NUMBER_OF_HEADER_2BYTES+NUMBER_OF_DATA_2BYTES,OutputFile);
			}
        }
    }
    
    
    // Destruct...
    free(Line);
    free(VectorR);
    free(VectorW);
    
    fclose (InputFile);
    fclose (OutputFile);
    
}