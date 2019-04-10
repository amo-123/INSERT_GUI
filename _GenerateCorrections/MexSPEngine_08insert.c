// mex MexSPEngine_08insert.c COMPFLAGS="$COMPFLAGS -openmp"  LINKFALGS="$LINKFALGS -openmp"

//-> Be van varrva a kodba, hogy 68*2 byte az adatmezo. Rakeresni 68-ra, 64-re es atirni, mint itt alul.

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
void Run1Mcount(double *LRFpIJ, int Loop, double *PMTxy, double *SpectrumMax, double SpectrumWindow, double *LUT, mwSize NumXlines, mwSize NumYlines, int NumPMT, double PE, double BaseLine, double *PicOut, double *Count, double *CountEw, double *LRFout, double *Spectrum, unsigned short *DataStream, int StreamDataLength)
{
    
    double Energy;             // Energy of the PMT vector, i.e. sum of the vector elements.
    double *PMTvectP;          // PMTvect(p).
    double PMT;

    int I, J, Iold, Jold;       // Coordinate parameters in pixel dimension.
    double Dist,Dist2;          // Dist refers to distance, in fact it is the Bhattacharyya coefficient.
    int step;                   // Step in pixes.
    int Iter;                   // Iteration loop variable.
    int l,k,m,M,p;              // Loop variables.
    int LUTi, LUTj;             // Look-up tables.
    
    
    PMTvectP = (double*) malloc (sizeof(double)*NumPMT);
    if (PMTvectP == NULL) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "PMTvect");
        return;
    }
    
    
    for (m=0; m<StreamDataLength; m++){
        M=m*(NUMBER_OF_HEADER_2BYTES+NUMBER_OF_DATA_2BYTES);
        
        if (DataStream[M+0]== 17492){ // If it is a PMT vector data packet
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
            
/*            for (p=0; p<NumPMT; p++){
                if (PMTvectP[p]>PMTmax){
                    PMTmax=PMTvectP[p];
                        J=(int)(PMTxy[p+NumPMT*1]);
                        I=(int)(PMTxy[p+NumPMT*0]);
                }
            }
*/            
            for (p=0; p<NumPMT; p++)
                PMTvectP[p]=sqrt(sqrt(PMTvectP[p]/Energy)); // This does the trick! We get the sqrt twice, see documentation!
            
            ///////////////////////////////////////////////////////////////
            if ( (Energy>0.5*PE)&&(Energy<2*PE) ) { // 1st ENERGY WINDOW //
                
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
                
                if ( (I>0)&&(I<1024)&&(J>0)&&(J<1024) ){ // If the coordinates are within the image.
                    
                    if (SpectrumWindow <=0)
                        PicOut[I+J*1024]-=1; // Pixel is cooled.
                    
                    if ( ( SpectrumWindow > 0 ) &&  ( Energy > PE*(1-SpectrumWindow) ) && ( Energy < PE*(1+SpectrumWindow) ) )
                        PicOut[I+J*1024]-=1;
                    
                    if (SpectrumWindow < 0){
                        if ( ( (int)(Energy/100)<500 )&& ( (int)(Energy/100)>=0 ) ){
                            Spectrum[ ((int)(Energy/100))*NumXlines*NumYlines ] += 1; // Count is put into the Spectrum as well.
                        }
                    }
                    
                    if (SpectrumWindow >=0){
                        LUTj=((int)LUT[I+J*1024+0]-1);           // LUT(I,J,0): coordinate "i", near to one of the crosses.
                        LUTi=((int)LUT[I+J*1024+1*1024*1024]-1); // LUT(I,J,1): coordinate "j", near to one of the crosses.
                        if ( ( LUTi >=0 ) && ( LUTj >=0 )  &&  ( ( SpectrumWindow <=0 ) || ( ( Energy > SpectrumMax[LUTi + LUTj*NumXlines]*(1-SpectrumWindow) ) && ( Energy < SpectrumMax[LUTi + LUTj*NumXlines]*(1+SpectrumWindow) ) ) ) ) {
                    // This condition ensures the following:
                    // - If SpectrumWindow >= 0 AND the LUT table value is positive, which signals that we are
                    // near to one of the crosses, then the count is registered in the LRFout and Spectrum
                    // matrices.
                    // - If SpectrumWindow == 0, then the corresponding pixel value of PicOut is increased.
                    // This make possible to check the LUT tables.
                    // - If SpectrumWindow < 0, then the only output is the image itself, the following part is skipped.
                            if (SpectrumWindow==0) 
                                PicOut[I+J*1024]+=2;
                            
                            CountEw[0]+=1;
                            
                            for (p=0; p<NumPMT; p++){
                                LRFout[ LUTi + LUTj*NumXlines + p*NumXlines*NumYlines ] += PMTvectP[p]*PMTvectP[p];
                            } // The proper element of LRF matrix is increased by the vector. (It has been sqrt-ed, therefor here it is on the 2nd power.)
                            LRFout[ LUTi + LUTj*NumXlines + NumPMT*NumXlines*NumYlines ] += 1; // Here we count the counts near to the given cross.
                            
                            if ( ( (int)(Energy/100)<500 )&& ( (int)(Energy/100)>=0 ) ){
                                Spectrum[ LUTi + LUTj*NumXlines + ((int)(Energy/100))*NumXlines*NumYlines ] += 1; // Count is put into the Spectrum as well.
                            }
                            
                        }
                    }
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
    double *SpectrumMax;        // Energy table. SpectrumMax(LineI,LineJ), where LineI and LineJ are the indices of phantom lines.
    double SpectrumWindow;      // Energy window. 0.1->10%
    double *LUT;                // LUT table of 1024X1024X2. 
                                // LUT(i,j,0/1)=0, if we are not near to any of the crosses.
                                // LUT(i,j,0)=LinesI, LUT(i,j,1)=LinesJ, when we are near to the cross of LinesI and LinesJ.
    mwSize NumXlines,NumYlines; // Number of phantom lines
    int NumPMT;                 // Number of PMTs.
    double PE=17700;            // Energy peak channel number.
	double BaseLine;

    // Left-hand values
    double *PicOut;             // Output image.
    double *Count;              // Counter for digested data packets.
    double *CountEw;            // Counter for count got through the energy window.
    double *LRFout;             // Output LRF matrix of size: LRF[NumXlines,NumYlines,NumPMT+1].
    double *Spectrum;           // Matrix for storing the spectra in all of the crosses.Size: Spectrum[NumXlines,NumYlines,500].
    mwSize dims[3];             // It is not left-hand value but needed for creating the matrices.

    // Left-hand clones created for parallelization purposes.
    double *PicOut01, *PicOut02, *PicOut03, *PicOut04;
    double *Count01, *Count02, *Count03, *Count04;
    double *CountEw01, *CountEw02, *CountEw03, *CountEw04;
    double *LRFout01, *LRFout02, *LRFout03, *LRFout04;
    double *Spectrum01, *Spectrum02, *Spectrum03, *Spectrum04;
    
    
    unsigned short *DataStream01, *DataStream02, *DataStream03, *DataStream04;
    // 136-bytes-long data packets to be read into these streams. //-> 136? - 68*2
    
    
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
    SpectrumMax = mxGetPr(prhs[4]);
    SpectrumWindow = (double)mxGetScalar(prhs[5]);
    LUT = mxGetPr(prhs[6]);
    NumXlines = (mwSize)mxGetScalar(prhs[7]);
    NumYlines = (mwSize)mxGetScalar(prhs[8]);
    NumPMT = (int)mxGetScalar(prhs[9]);
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
    
    dims[0]=NumXlines; dims[1]=NumYlines; dims[2]=NumPMT+1;
    plhs[3] = mxCreateNumericArray(3,&dims, mxDOUBLE_CLASS, mxREAL);
    LRFout = mxGetPr(plhs[3]);
    LRFout01 = mxGetPr(mxCreateNumericArray(3,&dims, mxDOUBLE_CLASS, mxREAL));
    LRFout02 = mxGetPr(mxCreateNumericArray(3,&dims, mxDOUBLE_CLASS, mxREAL));
    LRFout03 = mxGetPr(mxCreateNumericArray(3,&dims, mxDOUBLE_CLASS, mxREAL));
    LRFout04 = mxGetPr(mxCreateNumericArray(3,&dims, mxDOUBLE_CLASS, mxREAL));
    
    dims[0]=NumXlines; dims[1]=NumYlines; dims[2]=500;
    plhs[4] = mxCreateNumericArray(3,&dims, mxDOUBLE_CLASS, mxREAL);
    Spectrum = mxGetPr(plhs[4]);
    Spectrum01 = mxGetPr(mxCreateNumericArray(3,&dims, mxDOUBLE_CLASS, mxREAL));
    Spectrum02 = mxGetPr(mxCreateNumericArray(3,&dims, mxDOUBLE_CLASS, mxREAL));
    Spectrum03 = mxGetPr(mxCreateNumericArray(3,&dims, mxDOUBLE_CLASS, mxREAL));
    Spectrum04 = mxGetPr(mxCreateNumericArray(3,&dims, mxDOUBLE_CLASS, mxREAL));
    
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
        
        
        #pragma omp parallel sections shared(LRFpIJ, pFile, Loop, PMTxy, SpectrumMax, SpectrumWindow, LUT, NumXlines, NumYlines, NumPMT, PE, Count, CountEw, LRFout, Spectrum, StreamDataLength)
        {
            
            #pragma omp section
            {
                fread (DataStream01,sizeof(unsigned short),(NUMBER_OF_HEADER_2BYTES+NUMBER_OF_DATA_2BYTES)*StreamDataLength,pFile); // StreamDataLength*(NUMBER_OF_HEADER_2BYTES+NUMBER_OF_DATA_2BYTES)X2bájtos adatcsomag kiolvasva.
                Run1Mcount(LRFpIJ, Loop, PMTxy, SpectrumMax, SpectrumWindow, LUT, NumXlines, NumYlines, NumPMT, PE, BaseLine, PicOut01, Count01, CountEw01, LRFout01, Spectrum01, DataStream01, StreamDataLength);
            }
            #pragma omp section
            {
                fread (DataStream02,sizeof(unsigned short),(NUMBER_OF_HEADER_2BYTES+NUMBER_OF_DATA_2BYTES)*StreamDataLength,pFile); // StreamDataLength*(NUMBER_OF_HEADER_2BYTES+NUMBER_OF_DATA_2BYTES)X2bájtos adatcsomag kiolvasva.
                Run1Mcount(LRFpIJ, Loop, PMTxy, SpectrumMax, SpectrumWindow, LUT, NumXlines, NumYlines, NumPMT, PE, BaseLine, PicOut02, Count02, CountEw02, LRFout02, Spectrum02, DataStream02, StreamDataLength);
            }
            #pragma omp section
            {
                fread (DataStream03,sizeof(unsigned short),(NUMBER_OF_HEADER_2BYTES+NUMBER_OF_DATA_2BYTES)*StreamDataLength,pFile); // StreamDataLength*(NUMBER_OF_HEADER_2BYTES+NUMBER_OF_DATA_2BYTES)X2bájtos adatcsomag kiolvasva.
                Run1Mcount(LRFpIJ, Loop, PMTxy, SpectrumMax, SpectrumWindow, LUT, NumXlines, NumYlines, NumPMT, PE, BaseLine, PicOut03, Count03, CountEw03, LRFout03, Spectrum03, DataStream03, StreamDataLength);
            }
            #pragma omp section
            {
                fread (DataStream04,sizeof(unsigned short),(NUMBER_OF_HEADER_2BYTES+NUMBER_OF_DATA_2BYTES)*StreamDataLength,pFile); // StreamDataLength*(NUMBER_OF_HEADER_2BYTES+NUMBER_OF_DATA_2BYTES)X2bájtos adatcsomag kiolvasva.
                Run1Mcount(LRFpIJ, Loop, PMTxy, SpectrumMax, SpectrumWindow, LUT, NumXlines, NumYlines, NumPMT, PE, BaseLine, PicOut04, Count04, CountEw04, LRFout04, Spectrum04, DataStream04, StreamDataLength);
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
    
    for(l=0;l<NumXlines;l++)
        for(k=0;k<NumYlines;k++)
            for(p=0;p<=NumPMT;p++)
                LRFout[l+k*NumXlines+p*NumXlines*NumYlines]=LRFout01[l+k*NumXlines+p*NumXlines*NumYlines] + LRFout02[l+k*NumXlines+p*NumXlines*NumYlines] + LRFout03[l+k*NumXlines+p*NumXlines*NumYlines] + LRFout04[l+k*NumXlines+p*NumXlines*NumYlines];
    
    for(l=0;l<NumXlines;l++)
        for(k=0;k<NumYlines;k++)
            for(p=0;p<500;p++)
                Spectrum[l+k*NumXlines+p*NumXlines*NumYlines]=Spectrum01[l+k*NumXlines+p*NumXlines*NumYlines] + Spectrum02[l+k*NumXlines+p*NumXlines*NumYlines] + Spectrum03[l+k*NumXlines+p*NumXlines*NumYlines] + Spectrum04[l+k*NumXlines+p*NumXlines*NumYlines];
    
    
    Count[0]=Count01[0] + Count02[0] + Count03[0] + Count04[0];
    CountEw[0]=CountEw01[0] + CountEw02[0] + CountEw03[0] + CountEw04[0];
    
}


