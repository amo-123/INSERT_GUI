// mex MexSPEngine_10insertUCECLin.c COMPFLAGS="$COMPFLAGS -openmp"  LINKFALGS="$LINKFALGS -openmp"


//-> Be van varrva a kodba, hogy 68*2 byte az adatmezo. Rakeresni 68-ra, 64-re es atirni, mint itt alul.
// for (m=0; m<StreamDataLength; m++){ // StreamDataLength=25e2, as it can be seen in the main function.
//        M=m*68;
//        if ( (DataStream[M+0] == 17492) && (DataStream[M+64] <256 ) ){

// Az adatmezo hosszat ki akarod vezetni csinald ugy, mint a NumPMT-vel:
//void mexFunction( int nlhs, mxArray *plhs[],
//        int nrhs, const mxArray *prhs[] )
//{
// [...]
//    int NumPMT;                 // Number of PMTs.
// [...]
// NumPMT = (int)mxGetScalar(prhs[4]);
// prhs: right-hand arguments

// 4 szalra parhuzamositva, 4XStreamDataLength hosszu darabkakat olvasok be, a tobbit eldobom. StreamDataLength kodba bevarrva:
// int StreamDataLength=25e2;        // StreamDataLength * 4 = 1e4 is the size of the data bunches to be solved in a single outer iteration step.

            
#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <omp.h>

// System-based h files
#include "../clinicalParameters.h"


// Calculation of the Bhattacharyya coefficient or its square root
// Bhattacharyya distance is a strictily monotoneously decreasing function 
// of Bhattacharyya coefficinet, therefore optimum is considered to be 
// found in the maximum.
double Distance(double *LRFpIJ, double *PMTvectP, int I, int J, int NumPMT)
{
    int p;
    double Dist;
    Dist=0.;
    
    for (p=0; p<NumPMT; p++)
        Dist+=LRFpIJ[p+I*(NumPMT+1)+J*(NumPMT+1)*1024]*PMTvectP[p];
        // LRF and PMT detetctor values have already been square rooted.
        
    return Dist;
}


///////////////////////////////////////////////////////////////////////////
// Coordinate estimation
void Run1Mcount(double *LRFpIJ, int Loop, double *PMTxy, int NumPMT, double SpectrumWindow, double *EC, double *UC, double *LinX, double *LinY, double PE, double BaseLine, double *PicOut, double *Count, double *CountEw, unsigned short *DataStream, int StreamDataLength)
{
    
    double Energy;             // Energy of the PMT vector, i.e. sum of the vector elements.
    double *PMTvectP;          // PMTvect(p).
    double PMT;

    int I, J, Iold, Jold;       // Coordinate parameters in pixel dimension.
    double Dist,Dist2;          // Dist refers to distance, in fact it is the Bhattacharyya coefficient.
    int step;                   // Step in pixes.
    int Iter;                   // Iteration loop variable.
    int m,M,p;                  // Loop variables.
    
    double ECvalue,UCvalue;     // EC and UC correction values.
    
    
    PMTvectP = (double*) malloc (sizeof(double)*NumPMT);
    if (PMTvectP == NULL) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "PMTvect");
        return;
    }
    
    
    for (m=0; m<StreamDataLength; m++){ // StreamDataLength=25e2, as it can be seen in the main function.
        M=m*(NUMBER_OF_HEADER_2BYTES+NUMBER_OF_DATA_2BYTES);
        
        if (DataStream[M+0]== 17492){ // If it is a PMT vector data packet.
        //if ( (DataStream[M+0] == 17492) && (DataStream[M+(NUMBER_OF_HEADER_2BYTES+NUMBER_OF_DATA_2BYTES-4)] <256 ) ){ //-> Ezt a 64-et nem tudom, hogy kell-e változtatni
            
            double PMTmax=0;
            
            Count[0]+=1;
            
            for (p=0; p<NumPMT; p++)
                PMTvectP[p]=(double)DataStream[M+4+p]; // From the 5th element of the data packet, the vector is read out.

            for (p=0; p<NumPMT; p++){
                PMTvectP[p]-=BaseLine;
                if (PMTvectP[p]<0)
                    PMTvectP[p]=0;
            }
            
            Energy=0.;
            for (p=0; p<NumPMT; p++)
                Energy+=PMTvectP[p]; // Calculating energy.
            
            I=0; J=0; // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            for (p=0; p<NumPMT; p++){
                J+=(int)(PMTxy[p+NumPMT*1]*PMTvectP[p]);
                I+=(int)(PMTxy[p+NumPMT*0]*PMTvectP[p]);
            }
            I=(int)(I/Energy);
            J=(int)(J/Energy); // Extremely simple centroid method for calculating some starting coordinates.
            

            for (p=0; p<NumPMT; p++)
                PMTvectP[p]=sqrt(sqrt(PMTvectP[p]/Energy)); // This does the trick! We get the sqrt twice, see documentation!
            
            ///////////////////////////////////////////////////////////////
            if ( (Energy>0.1*PE)&&(Energy<5.2*PE) ) { // 1st ENERGY WINDOW // Now it is practically meaningless.
                
                Dist=Distance(LRFpIJ, PMTvectP, I, J, NumPMT); // Distance in the starting point.
                Iold=-1;
                Jold=-1;
                Iter=0;
                // CountEw1[0]+=1; // We do not count since 1st window is meaningless.
                
                while ( (Iold!=I)&&(Jold!=J)||(Iter==1) ) // Loop until the actual and previous coordinates are equal and at least two loops has been run.
                {
                    Iold=I;
                    Jold=J;
                    Iter++;
                    if(Iter==1) { step=5; } else { step=1; } // Int the 1st loop step is 5, otherwise it is 1.
                    
                    Dist2=Distance(LRFpIJ, PMTvectP, I+step, J, NumPMT); // Whether it is increasing in the +I direction?
                    if (Dist2>Dist){
                        while (Dist2>Dist){
                            Dist=Dist2;
                            if (I<1014) {I=I+step;}
                            Dist2=Distance(LRFpIJ, PMTvectP, I+step, J, NumPMT);
                        }
                    }
                    else{
                        Dist2=Distance(LRFpIJ, PMTvectP, I-step, J, NumPMT); // If not, then we go to the -I direction.
                        while (Dist2>Dist){
                            Dist=Dist2;
                            if (I>10) {I=I-step;}
                            Dist2=Distance(LRFpIJ, PMTvectP, I-step, J, NumPMT);
                        }
                    }
                    
                    Dist2=Distance(LRFpIJ, PMTvectP, I, J+step, NumPMT); // Whether it is increasing in the +J direction?
                    if (Dist2>Dist){
                        while (Dist2>Dist){
                            Dist=Dist2;
                            if (J<1014) {J=J+step;}
                            Dist2=Distance(LRFpIJ, PMTvectP, I, J+step, NumPMT);
                        }
                    }
                    else{
                        Dist2=Distance(LRFpIJ, PMTvectP, I, J-step, NumPMT); // If not, then we go to the -J direction.
                        while (Dist2>Dist){
                            Dist=Dist2;
                            if (J>10) {J=J-step;}
                            Dist2=Distance(LRFpIJ, PMTvectP, I, J-step, NumPMT);
                        }
                    }
                }
                
                ECvalue = EC[I+J*1024];
                
                // If the coordinates are within the image and we are in the energy window.
                if ( (I>2) && (I<1022) && (J>2) && (J<1022) && (Energy<ECvalue*(1+SpectrumWindow)) && (Energy>ECvalue*(1-SpectrumWindow)) ){ // Ha egyáltalán benn maradtunk a képmátrixban.
                    
                    // Within the following part, the coordinates are transformed from integer pixel values into real pixel values.
                    int In,Jn;
                    double Idbl,Jdbl,DI,DJ;
                    
                    double Zmm,Zm0,Zmp,Z0m,Z00,Z0p,Zpm,Zp0,Zpp;
                    double a0,a1,a2,a3,a4,a5,a6,a7,a8;
                    
                    double cx,cy,N1,N2,N3,N4,ec1,ec2,ec3,ec4;
                    
                    
                    
                    
                    In=I-1;
                    
                    Zmm=0; Jn=J-1;
                    for(p=0; p<NumPMT; p++)
                        Zmm+=LRFpIJ[p+In*(NumPMT+1)+Jn*(NumPMT+1)*1024]*PMTvectP[p];
                    
                    Zm0=0; Jn=J;
                    for(p=0; p<NumPMT; p++)
                        Zm0+=LRFpIJ[p+In*(NumPMT+1)+Jn*(NumPMT+1)*1024]*PMTvectP[p];
                    
                    Zmp=0; Jn=J+1;
                    for(p=0; p<NumPMT; p++)
                        Zmp+=LRFpIJ[p+In*(NumPMT+1)+Jn*(NumPMT+1)*1024]*PMTvectP[p];
                    
                    In=I;
                    
                    Z0m=0; Jn=J-1;
                    for(p=0; p<NumPMT; p++)
                        Z0m+=LRFpIJ[p+In*(NumPMT+1)+Jn*(NumPMT+1)*1024]*PMTvectP[p];
                    
                    Z00=0; Jn=J;
                    for(p=0; p<NumPMT; p++)
                        Z00+=LRFpIJ[p+In*(NumPMT+1)+Jn*(NumPMT+1)*1024]*PMTvectP[p];
                    
                    Z0p=0; Jn=J+1;
                    for(p=0; p<NumPMT; p++)
                        Z0p+=LRFpIJ[p+In*(NumPMT+1)+Jn*(NumPMT+1)*1024]*PMTvectP[p];
                    
                    In=I+1;
                    
                    Zpm=0; Jn=J-1;
                    for(p=0; p<NumPMT; p++)
                        Zpm+=LRFpIJ[p+In*(NumPMT+1)+Jn*(NumPMT+1)*1024]*PMTvectP[p];
                    
                    Zp0=0; Jn=J;
                    for(p=0; p<NumPMT; p++)
                        Zp0+=LRFpIJ[p+In*(NumPMT+1)+Jn*(NumPMT+1)*1024]*PMTvectP[p];
                    
                    Zpp=0; Jn=J+1;
                    for(p=0; p<NumPMT; p++)
                        Zpp+=LRFpIJ[p+In*(NumPMT+1)+Jn*(NumPMT+1)*1024]*PMTvectP[p];
                    
                    
                    a0=                                   Z00;
                    a1=         -0.5*Zm0                                     +0.5*Zp0;
                    a2=                          -0.5*Z0m    +0.5*Z0p;
                    a3= 0.25*Zmm        -0.25*Zmp                    -0.25*Zpm        +0.25*Zpp;
                    a4=          0.5*Zm0                 -Z00                +0.5*Zp0;
                    a5=                           0.5*Z0m-Z00+0.5*Z0p;
                    a6=-0.25*Zmm        +0.25*Zmp+0.5*Z0m    -0.5*Z0p-0.25*Zpm        +0.25*Zpp;
                    a7=-0.25*Zmm+0.5*Zm0-0.25*Zmp                    +0.25*Zpm-0.5*Zp0+0.25*Zpp;
                    a8= 0.25*Zmm-0.5*Zm0+0.25*Zmp-0.5*Z0m+Z00-0.5*Z0p+0.25*Zpm-0.5*Zp0+0.25*Zpp;
// 0. iteration
                    DI=-a1/a4/2;
                    DJ=-a2/a5/2;
// // 1. iteration
//                     DI=-(a1+a3*DJ+a7*DJ*DJ)/(a4+a6*DJ+a8*DJ*DJ)/2;
//                     DJ=-(a2+a3*DI+a6*DI*DI)/(a5+a7*DI+a8*DI*DI)/2;
// // 2. iteration
//                     DI=-(a1+a3*DJ+a7*DJ*DJ)/(a4+a6*DJ+a8*DJ*DJ)/2;
//                     DJ=-(a2+a3*DI+a6*DI*DI)/(a5+a7*DI+a8*DI*DI)/2;
                    
                    //Idbl=(double)I;
                    //Jdbl=(double)J;
                    
                    // Somethibg was wrong with it, I could have not found it out. The below simple smoothing is inserted instead.
                    Idbl=(double)I + 2*DI; // Should be Idbl=(double)I + DI;
                    Jdbl=(double)J + 2*DJ; // Should be Jdbl=(double)J + DJ;
                    
                    
                    // Below, the linearity correction is performed.
                    cx = Idbl-floor(Idbl);
                    cy = Jdbl-floor(Jdbl);
                    
                    N1 = (1.0-cx) * (1.0-cy);
                    N2 = (cx)     * (1.0-cy);
                    N3 = (cx)     * (cy);
                    N4 = (1.0-cx) * (cy);
                    
                    ec1 = LinY[ I   + 1024*J ] / 64;
                    ec2 = LinY[ I+1 + 1024*J ] / 64;
                    ec3 = LinY[ I+1 + 1024*(J+1) ] / 64;
                    ec4 = LinY[ I   + 1024*(J+1) ] / 64;
                    
                    Idbl = ( Idbl + (ec1*N1 + ec2*N2 + ec3*N3 + ec4*N4) ); //+0.5
                    
                    ec1 = LinX[ I   + 1024*J ] / 64;
                    ec2 = LinX[ I+1 + 1024*J ] / 64;
                    ec3 = LinX[ I+1 + 1024*(J+1) ] / 64;
                    ec4 = LinX[ I   + 1024*(J+1) ] / 64;
                    
                    Jdbl = ( Jdbl + (ec1*N1 + ec2*N2 + ec3*N3 + ec4*N4) ); //+0.5
                    
                    
                    I = (int)(Idbl);
                    J = (int)(Jdbl);
                    
                    // Pixel value is increased with UC correction value.
                    UCvalue = UC[I+J*1024];
                    
                    if ( (I>0) && (I<1024) && (J>0) && (J<1024) )
                        PicOut[I+J*1024]+=UCvalue;
                    
                    // Counter CountEw is increased.
                    CountEw[0]+=1;
                    
                }
            }
        }
        
    }
    free (PMTvectP);
}
// End of coordinate estimation
///////////////////////////////////////////////////////////////////////////




void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] )
{
    ///////////////////////////////////////////////////////////////////////
    // Definitions ////////////////////////////////////////////////////////
    
    // Right-hand values
    double *LRFpIJ;             // LRF(p,I,J).
    FILE *pFile;                // File pointer.
    char *path;                 // Path to file. This is the right hand value itself.
    int Loop;                   // Number of the 136-bytes-long data packets in the file.
    double *PMTxy;              // Coordintes of the PMTs: PMT(x,y)
    int NumPMT;                 // Number of PMTs.
    double SpectrumWindow;      // Energy window. 0.1->10%
    double *EC;                 // EC table
    double *UC;                 // UC table
    double *LinX;               // X lin table
    double *LinY;               // Y lin table
    double PE;                  // Energy peak channel number.
	double BaseLine;

    // Left-hand values
    double *PicOut;             // Output image.
    double *Count;              // Counter for digested data packets.
    double *CountEw;            // Counter for count got through the energy window.
    mwSize dims[3];             // It is not left-hand value but needed for creating the matrices.

    // Left-hand clones created for parallelization purposes.
    double *PicOut01, *PicOut02, *PicOut03, *PicOut04;
    double *Count01, *Count02, *Count03, *Count04;
    double *CountEw01, *CountEw02, *CountEw03, *CountEw04;
    
    
    unsigned short *DataStream01, *DataStream02, *DataStream03, *DataStream04;
    // 136-bytes-long data packets to be read into these streams.


    int l,k,p,LoopM;                  // Loop variales.
    int StreamDataLength=25e2;        // StreamDataLength * 4 = 1e4 is the size of the data bunches to be solved in a single outer iteration step.
    
    ///////////////////////////////////////////////////////////////////////
    // Initializing input and output fields and values ////////////////////
    // Special MEX technology, see MatLab help. ///////////////////////////
    
    // Right-hand values
    LRFpIJ = mxGetPr(prhs[0]);
    path = mxArrayToString(prhs[1]);
    Loop = (int)mxGetScalar(prhs[2]);
    PMTxy = mxGetPr(prhs[3]);
    NumPMT = (int)mxGetScalar(prhs[4]);
    SpectrumWindow = (double)mxGetScalar(prhs[5]);
    EC = mxGetPr(prhs[6]);
    UC = mxGetPr(prhs[7]);
    LinX = mxGetPr(prhs[8]);
    LinY = mxGetPr(prhs[9]);
    PE=(double)mxGetScalar(prhs[10]);
	BaseLine=(double)mxGetScalar(prhs[11]);
    
    // Left-hand values and their clones
    plhs[0] = mxCreateDoubleMatrix(1024, 1024, mxREAL);
    PicOut = mxGetPr(plhs[0]);
    PicOut01 = mxGetPr(mxCreateDoubleMatrix(1024, 1024, mxREAL));
    PicOut02 = mxGetPr(mxCreateDoubleMatrix(1024, 1024, mxREAL));
    PicOut03 = mxGetPr(mxCreateDoubleMatrix(1024, 1024, mxREAL));
    PicOut04 = mxGetPr(mxCreateDoubleMatrix(1024, 1024, mxREAL));
    
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    Count = mxGetPr(plhs[1]);
    Count[0]=0.;
    Count01 = mxGetPr(mxCreateDoubleMatrix(1, 1, mxREAL));
    Count02 = mxGetPr(mxCreateDoubleMatrix(1, 1, mxREAL));
    Count03 = mxGetPr(mxCreateDoubleMatrix(1, 1, mxREAL));
    Count04 = mxGetPr(mxCreateDoubleMatrix(1, 1, mxREAL));
    Count01[0]=0.;Count02[0]=0.;Count03[0]=0.;Count04[0]=0.;
    
    plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
    CountEw = mxGetPr(plhs[2]);
    CountEw[0]=0.;
    CountEw01 = mxGetPr(mxCreateDoubleMatrix(1, 1, mxREAL));
    CountEw02 = mxGetPr(mxCreateDoubleMatrix(1, 1, mxREAL));
    CountEw03 = mxGetPr(mxCreateDoubleMatrix(1, 1, mxREAL));
    CountEw04 = mxGetPr(mxCreateDoubleMatrix(1, 1, mxREAL));
    CountEw01[0]=0.;CountEw02[0]=0.;CountEw03[0]=0.;CountEw04[0]=0.;
    

    // File and datastream
    pFile = fopen ( path , "rb" );
    if (pFile==NULL) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "path");
        return;
    }
  
    DataStream01 = (int*) malloc (sizeof(unsigned short)*(NUMBER_OF_HEADER_2BYTES+NUMBER_OF_DATA_2BYTES)*StreamDataLength);
    if (DataStream01 == NULL) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "DataStream01");
        return;
    }
    DataStream02 = (int*) malloc (sizeof(unsigned short)*(NUMBER_OF_HEADER_2BYTES+NUMBER_OF_DATA_2BYTES)*StreamDataLength);
    if (DataStream02 == NULL) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "DataStream02");
        return;
    }
    DataStream03 = (int*) malloc (sizeof(unsigned short)*(NUMBER_OF_HEADER_2BYTES+NUMBER_OF_DATA_2BYTES)*StreamDataLength);
    if (DataStream03 == NULL) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "DataStream03");
        return;
    }
    DataStream04 = (int*) malloc (sizeof(unsigned short)*(NUMBER_OF_HEADER_2BYTES+NUMBER_OF_DATA_2BYTES)*StreamDataLength);
    if (DataStream04 == NULL) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "DataStream04");
        return;
    }
    
    
    
    ///////////////////////////////////////////////////////////////////////
    // Here begins the main parallelized loop for data evaluation.
    // Parallellization for 4 threads is fixed!!!!
    
    for (LoopM=0; LoopM<Loop; LoopM+=StreamDataLength*4){
        
        
        #pragma omp parallel sections shared(LRFpIJ, pFile, Loop, PMTxy, NumPMT, SpectrumWindow, EC, UC, LinX, LinY, PE, Count, CountEw, StreamDataLength)
        {
            
            #pragma omp section
            {
                fread (DataStream01,sizeof(unsigned short),(NUMBER_OF_HEADER_2BYTES+NUMBER_OF_DATA_2BYTES)*StreamDataLength,pFile); // StreamDataLength*(NUMBER_OF_HEADER_2BYTES+NUMBER_OF_DATA_2BYTES)X2-bytes data packet has been read.
                Run1Mcount(LRFpIJ, Loop, PMTxy, NumPMT, SpectrumWindow, EC, UC, LinX, LinY, PE, BaseLine,  PicOut01, Count01, CountEw01, DataStream01, StreamDataLength);
            }
            #pragma omp section
            {
                fread (DataStream02,sizeof(unsigned short),(NUMBER_OF_HEADER_2BYTES+NUMBER_OF_DATA_2BYTES)*StreamDataLength,pFile); // StreamDataLength*(NUMBER_OF_HEADER_2BYTES+NUMBER_OF_DATA_2BYTES)X2-bytes data packet has been read.
                Run1Mcount(LRFpIJ, Loop, PMTxy, NumPMT, SpectrumWindow, EC, UC, LinX, LinY, PE, BaseLine, PicOut02, Count02, CountEw02, DataStream02, StreamDataLength);
            }
            #pragma omp section
            {
                fread (DataStream03,sizeof(unsigned short),(NUMBER_OF_HEADER_2BYTES+NUMBER_OF_DATA_2BYTES)*StreamDataLength,pFile); // StreamDataLength*(NUMBER_OF_HEADER_2BYTES+NUMBER_OF_DATA_2BYTES)X2-bytes data packet has been read.
                Run1Mcount(LRFpIJ, Loop, PMTxy, NumPMT, SpectrumWindow, EC, UC, LinX, LinY, PE, BaseLine, PicOut03, Count03, CountEw03, DataStream03, StreamDataLength);
            }
            #pragma omp section
            {
                fread (DataStream04,sizeof(unsigned short),(NUMBER_OF_HEADER_2BYTES+NUMBER_OF_DATA_2BYTES)*StreamDataLength,pFile); // StreamDataLength*(NUMBER_OF_HEADER_2BYTES+NUMBER_OF_DATA_2BYTES)X2-bytes data packet has been read.
                Run1Mcount(LRFpIJ, Loop, PMTxy, NumPMT, SpectrumWindow, EC, UC, LinX, LinY, PE, BaseLine, PicOut04, Count04, CountEw04, DataStream04, StreamDataLength);
            }
            
        }
    }
        
    fclose (pFile);

    free (DataStream01);
    free (DataStream02);
    free (DataStream03);
    free (DataStream04);
    mxFree(path);
    // Matrices created with mxCreated should be freed as well?
    
    
    // Summing the values stored in the clones.
    for(l=0;l<1024;l++)
        for(k=0;k<1024;k++)
            PicOut[l+k*1024]=PicOut01[l+k*1024] + PicOut02[l+k*1024] + PicOut03[l+k*1024] + PicOut04[l+k*1024];
    
    
    Count[0]=Count01[0] + Count02[0] + Count03[0] + Count04[0];
    CountEw[0]=CountEw01[0] + CountEw02[0] + CountEw03[0] + CountEw04[0];
    
}


