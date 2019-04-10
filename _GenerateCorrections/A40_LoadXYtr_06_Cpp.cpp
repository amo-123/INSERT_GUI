// mex A40_LoadXYtr_06_Cpp.cpp COMPFLAGS="$COMPFLAGS -openmp"  LINKFALGS="$LINKFALGS -openmp"

#include <stdio.h>
#include <math.h>

#include "mex.h"
#include "matrix.h"

// Kisegíto függvények %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// A0_LTGfunctions.m file-ban tárolva az eredetik, meg hogy hol használjuk
// melyiket. A párhuzamosítás miatt, ebben a modulban, a kódba is be kellett
// varrni ezeket a függvényeket.

int y2i( float y, float dy )
{
	return ceil( y / dy );
}

int x2j( double x, double dx )
{
	return ceil( x / dx );
}

double i2y( int i, double dy )
{
	return dy * ( ( double ) i - 0.5 );
}

double j2x( int j, double dx)
{
	return dx * ( ( double ) j - 0.5 );
}

double Jp2x( int Jpix, double UpScale, double dx )
{
	return dx * ( (double) Jpix - 0.5) / UpScale;
}

double Ip2y( int Ipix, double UpScale, double dy )
{
	return dy * ( Ipix - 0.5 ) / UpScale;
}

int IL2I( int IL, int HalfIline, int Dy, double dy, int Iimage, double UpScale )
{
	double y;
    y = Dy * ( IL - HalfIline ) + dy * ( Iimage / 2 - 1 );
	return ceil( y * UpScale / dy );
}

int JL2J( int JL, int HalfJline, int Dx, double dx, int Jimage, double UpScale  )
{
	double x;
    x = Dx * ( JL - HalfJline ) + dx * ( Jimage / 2 - 1 );
	return ceil( x * UpScale / dx );
}

int Jline2x( int Jline, int HalfJline, int Dx, double dx, int Jimage  )
{
	return Dx * ( Jline - HalfJline ) + dx * ( Jimage / 2 - 1 );
}
 
int Iline2y( int Iline, int HalfIline, int Dy, double dy, int Iimage )
{
	return Dy *( Iline - HalfIline ) + dy * ( Iimage / 2 - 1 );
}

// Nincs szétválasztva i-ben és j-ben!!!!
void CompleteMatrix( double *A, int M, int N, int Iimage )
{
	int i, j;
	for ( i = 1; i <= ( Iimage / 2 - 1); i++ )
	{
		for ( j = ( Iimage / 2 + 1 ); j <= ( Iimage / 2 + 1 + i ); j++ )
		{
			if ( A[ Iimage / 2 - i - 1 + M * ( j - 1 ) ] == 0 )
				 A[ Iimage / 2 - i - 1 + M * ( j - 1 ) ]         = A[ Iimage / 2 - i - 1 + 1 + M * ( j - 1 ) ];
				
			if ( A[ Iimage / 2 + 1 + i - 1 + M * ( j - 1 ) ] == 0 )
				 A[ Iimage / 2 + 1 + i - 1 + M * ( j - 1 ) ]     = A[ Iimage / 2 + 1 + i - 1 - 1 + M * ( j - 1 ) ];
				
			if ( A[ j - 1 + M * ( Iimage / 2 - i - 1 ) ] == 0 )
				 A[ j - 1 + M * ( Iimage / 2 - i - 1 ) ]     = A[ j - 1 + M * ( Iimage / 2 - i - 1 + 1 ) ];
				
			if ( A[ j - 1 + M * ( Iimage / 2 + 1 + i - 1 ) ] == 0 )
				 A[ j - 1 + M * ( Iimage / 2 + 1 + i - 1 ) ] = A[ j - 1 + M * ( Iimage / 2 + 1 + i - 1 - 1 ) ];
				
		}
		
		for ( j = ( Iimage / 2 ); j >= ( Iimage / 2 - i ); j-- )
		{
			
			if ( A[ Iimage / 2 - i - 1 + M * ( j - 1 ) ] == 0 )
				 A[ Iimage / 2 - i - 1 + M * ( j - 1 ) ] = A[ Iimage / 2 - i - 1 + 1 + M * ( j - 1 ) ];
				
			if ( A[ Iimage / 2 + 1 + i - 1 + M * ( j - 1 ) ] == 0 )
				 A[ Iimage / 2 + 1 + i - 1 + M * ( j - 1 ) ] = A[ Iimage / 2 + 1 + i - 1 - 1 + M * ( j - 1 ) ];
				
			if ( A[ j - 1 + M * ( Iimage / 2 - i - 1 ) ] == 0 )
				 A[ j - 1 + M * ( Iimage / 2 - i - 1 ) ] = A[ j - 1 + M * ( Iimage / 2 - i - 1 + 1 ) ];
				
			if ( A[ j - 1 + M * ( Iimage / 2 + 1 + i - 1 ) ] == 0 )
				 A[ j - 1 + M * ( Iimage / 2 + 1 + i - 1 ) ] = A[ j - 1 + M * ( Iimage / 2 + 1 + i - 1 - 1 ) ];
		}
	}
}


/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
			     int nrhs, const mxArray *prhs[])
{
	double *ThetaX;
	int   ThetaX_n;
	int   ThetaX_m;
	
	double *ThetaY;
    int   ThetaY_n;
	int   ThetaY_m;
	
	double *M;        
	int   M_n;       
	int   M_m;       
	
	int UpScale;   
	double dy;        
	double dx;        
	double HalfIline; 
	double HalfJline; 
	double Dx;        
	double Dy;        
	int NumIlines; 
	int NumJlines; 
	int Iimage;   
	int Jimage;

	double *Out;  

	int IimageUp;
	int JimageUp;
	
	int DimJL, DimIL, DimIJ;
	
	double *Tx;
	double *Ty;   
	double *alfaX; 
	double *alfaY;
	
	double *LinTblX;
	double *LinTblY;
	
	double Xtr;
	double Ytr;
	
	//Indexes
	int IL, JL;
	
	int Ibegin, Iend, Jbegin,Jend;
	int IL3, JL3;
	int indexMatrixI, indexMatrixJ;
	double sum;
	int I, J;
	
	double x, Xwv;
	double y, Ywv;
	int i, j;
	
	double *XwvStored;
	double *YwvStored;
	double iCenter;
	double jCenter;
	
	int indexX, indexY;

	if( nrhs != 14 )
	{
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "14 inputs required.");
	}
	
	if( nlhs != 2 ) 
	{
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs", "2 outputs required.");
	}
	
	ThetaX    = mxGetPr( prhs[0] );
	ThetaX_n  = mxGetN(  prhs[0] );
	ThetaX_m  = mxGetM(  prhs[0] );
	
	ThetaY    = mxGetPr( prhs[1] );
	ThetaY_n  = mxGetN(  prhs[1] );
	ThetaY_m  = mxGetM(  prhs[1] );
	
	M         = mxGetPr( prhs[2] );
	M_n       = mxGetN(  prhs[2] );
	M_m       = mxGetM(  prhs[2] );
	
	UpScale   = (int)mxGetScalar( prhs[3] );
	dy        = (double)mxGetScalar( prhs[4] );
	dx        = (double)mxGetScalar( prhs[5] );
	HalfIline = (double)mxGetScalar( prhs[6] );
	HalfJline = (double)mxGetScalar( prhs[7] );
	Dx        = (double)mxGetScalar( prhs[8] );
	Dy        = (double)mxGetScalar( prhs[9] );
	NumIlines = (int)mxGetScalar( prhs[10] );
	NumJlines = (int)mxGetScalar( prhs[11] );
	Iimage    = (int)mxGetScalar( prhs[12] );
	Jimage    = (int)mxGetScalar( prhs[13] );
	
	
	IimageUp = Iimage * UpScale; 
	JimageUp = Jimage * UpScale;
	
	DimJL = NumJlines + 2;
	DimIL = NumIlines + 2;
	
	
	if( HalfIline == 0 )
	{
		HalfIline = (double)DimIL / 2.0 + 0.5;
	}else{
		HalfIline = HalfIline + 1;
	}
	
	if( HalfJline == 0 )
	{
		HalfJline = (double)DimJL / 2.0 + 0.5;
	}else{
		HalfJline = HalfJline + 1;
	}
	
	printf("HalfIline: %lf HalfJline: %lf\n", HalfIline, HalfJline);
	//printf( "%d %lf %d %lf\n", DimIL, HalfIline, DimJL, HalfJline );
	
	plhs[0] = mxCreateDoubleMatrix( Iimage, Jimage, mxREAL );
	plhs[1] = mxCreateDoubleMatrix( Iimage, Jimage, mxREAL );
	/* get a pointer to the real data in the output matrix */
	LinTblX = mxGetPr( plhs[0] );
	LinTblY = mxGetPr( plhs[1] );
	
	Tx    = (double*)malloc( sizeof( double ) * 16 );
	Ty    = (double*)malloc( sizeof( double ) * 16 );
	alfaX = (double*)malloc( sizeof( double ) * 16 );
	alfaY = (double*)malloc( sizeof( double ) * 16 );
	
	XwvStored  = (double*)malloc( sizeof( double ) * Iimage * Jimage );
	YwvStored  = (double*)malloc( sizeof( double ) * Iimage * Jimage );
	
	for ( i = 0; i < Iimage; i++ )
	{
		for ( j = Jbegin; j < Jimage; j++ )
		{	
			LinTblX[ i + Jimage * j ] = 0.0;
			LinTblY[ i + Jimage * j ] = 0.0;
			XwvStored[ i + Jimage * j ] = 0.0;
			YwvStored[ i + Jimage * j ] = 0.0;
		}
	}
	
	for( IL = 1; IL <= ( DimIL - 1 ); IL++  )
	{
		Ibegin = ( int )ceil( ( (double)Dy * (   (double)IL - HalfIline )         + dy * ( (double)Iimage / 2.0 - 1.0 ) ) * (double)UpScale / dy );
		Iend   = ( int )ceil( ( (double)Dy * ( ( (double)IL + 1.0 ) - HalfIline ) + dy * ( (double)Iimage / 2.0 - 1.0 ) ) * (double)UpScale / dy ) - 1;
		
		//if( IL == 0 )
		//	printf("%d %d\n", Ibegin, Iend );
		
		if ( ( Ibegin > 0 ) || ( Iend <= IimageUp ) )
		{
		
			IL3 = IL * 3 - 2;
		
			for ( JL = 1; JL <= ( DimJL - 1 ); JL++  )
			{
				Jbegin = ( int )ceil( ( (double)Dx * (double)(   JL - HalfJline )       + dx * ( (double)Jimage / 2.0 - 1.0 ) ) * (double)UpScale / dx );
				Jend   = ( int )ceil( ( (double)Dx * (double)( ( JL + 1 ) - HalfJline ) + dx * ( (double)Jimage / 2.0 - 1.0 ) ) * (double)UpScale / dx ) - 1;
				
				if ( ( Jbegin > 0 ) || ( Jend <= JimageUp ) )
				{
				
					JL3 = JL * 3 - 2;
					// A külso ciklusokat a lineáris fantom vonalak által kijelölt
					// elemek szintjén szervezzük meg.
					
					// Kiszámítjuk az adott elemben az együttható vektort.
					//Tx(1:4)   = ThetaX( IL3    , JL3:(JL3+3));
					
					Tx[ 0 ]  = ThetaX[ IL3 - 1     + ThetaX_m * ( JL3 - 1 ) ]; Tx[ 1 ]  = ThetaX[ IL3 - 1 +     ThetaX_m * ( JL3 - 1 + 1) ]; Tx[ 2 ]  = ThetaX[ IL3 - 1     + ThetaX_m * ( JL3 - 1 + 2) ]; Tx[ 3 ]  = ThetaX[ IL3 - 1     + ThetaX_m * ( JL3 - 1 + 3) ];
					//Tx(5:8)   = ThetaX( IL3 + 1, JL3:(JL3+3));
					
					Tx[ 4 ]  = ThetaX[ IL3 - 1 + 1 + ThetaX_m * ( JL3 - 1 ) ]; Tx[ 5 ]  = ThetaX[ IL3 - 1 + 1 + ThetaX_m * ( JL3 - 1 + 1) ]; Tx[ 6 ]  = ThetaX[ IL3 - 1 + 1 + ThetaX_m * ( JL3 - 1 + 2) ]; Tx[ 7 ]  = ThetaX[ IL3 - 1 + 1 + ThetaX_m * ( JL3 - 1 + 3) ];
					//Tx(9:12)  = ThetaX( IL3 + 2, JL3:(JL3+3));
					
					Tx[ 8 ]  = ThetaX[ IL3 - 1 + 2 + ThetaX_m * ( JL3 - 1 ) ]; Tx[ 9 ]  = ThetaX[ IL3 - 1 + 2 + ThetaX_m * ( JL3 - 1 + 1) ]; Tx[ 10 ] = ThetaX[ IL3 - 1 + 2 + ThetaX_m * ( JL3 - 1 + 2) ]; Tx[ 11 ] = ThetaX[ IL3 - 1 + 2 + ThetaX_m * ( JL3 - 1 + 3) ];
					
					//Tx(13:16) = ThetaX( IL3 + 3, JL3:(JL3+3));
					Tx[ 12 ] = ThetaX[ IL3 - 1 + 3 + ThetaX_m * ( JL3 - 1 ) ]; Tx[ 13 ] = ThetaX[ IL3 - 1 + 3 + ThetaX_m * ( JL3 - 1 + 1) ]; Tx[ 14 ] = ThetaX[ IL3  - 1+ 3 + ThetaX_m * ( JL3 - 1 + 2) ]; Tx[ 15 ] = ThetaX[ IL3 - 1 + 3 + ThetaX_m * ( JL3 - 1 + 3) ];
					
					//alfaX = Mi * Tx';
					for( indexMatrixI = 0; indexMatrixI < 16; indexMatrixI++ )
					{
						sum = 0.0;
						for( indexMatrixJ = 0; indexMatrixJ < 16; indexMatrixJ++ )
						{
							sum += Tx[ indexMatrixJ ] * M[ indexMatrixI + M_m * indexMatrixJ ];
						}
						alfaX[ indexMatrixI ] = sum;
						
						//if( JL == 0 )
						//	printf( "%.30lf\t", sum );
							//printf( "T:%lf Sum:%lf\t", Tx[ indexMatrixI ], sum );
					}
					
					//if( JL == 0 )	
					//	printf( "%\n" );
					
					//Ty(1:4)   = ThetaY( IL3    , JL3:(JL3+3));
					Ty[ 0 ]  = ThetaY[ IL3 - 1 + ThetaY_m * ( JL3 - 1 ) ];     Ty[ 1 ]  = ThetaY[ IL3 - 1 + ThetaY_m * ( JL3 - 1 + 1) ];     Ty[ 2 ]  = ThetaY[ IL3 - 1 + ThetaY_m * ( JL3 - 1 + 2) ];     Ty[ 3 ]  = ThetaY[ IL3 - 1 + ThetaY_m * ( JL3 - 1 + 3) ];
					
					//Ty(5:8)   = ThetaY( IL3 + 1, JL3:(JL3+3));
					Ty[ 4 ]  = ThetaY[ IL3 - 1 + 1 + ThetaY_m * ( JL3 - 1 ) ]; Ty[ 5 ]  = ThetaY[ IL3 - 1 + 1 + ThetaY_m * ( JL3 - 1 + 1) ]; Ty[ 6 ]  = ThetaY[ IL3 - 1 + 1 + ThetaY_m * ( JL3 - 1 + 2) ]; Ty[ 7 ]  = ThetaY[ IL3 - 1 + 1 + ThetaY_m * ( JL3 - 1 + 3) ];
					
					//Ty(9:12)  = ThetaY( IL3 + 2, JL3:(JL3+3));
					Ty[ 8 ]  = ThetaY[ IL3 - 1 + 2 + ThetaY_m * ( JL3 - 1 ) ]; Ty[ 9 ]  = ThetaY[ IL3 - 1 + 2 + ThetaY_m * ( JL3 - 1 + 1) ]; Ty[ 10 ] = ThetaY[ IL3 - 1 + 2 + ThetaY_m * ( JL3 - 1 + 2) ]; Ty[ 11 ] = ThetaY[ IL3 - 1 + 2 + ThetaY_m * ( JL3 - 1 + 3) ];
					
					//Ty(13:16) = ThetaY( IL3 + 3, JL3:(JL3+3));
					Ty[ 12 ] = ThetaY[ IL3 - 1 + 3 + ThetaY_m * ( JL3 - 1 ) ]; Ty[ 13 ] = ThetaY[ IL3 - 1 + 3 + ThetaY_m * ( JL3 - 1 + 1) ]; Ty[ 14 ] = ThetaY[ IL3 - 1 + 3 + ThetaY_m * ( JL3 - 1 + 2) ]; Ty[ 15 ] = ThetaY[ IL3 - 1 + 3 + ThetaY_m * ( JL3 - 1 + 3) ];
					
					//alfaY = Mi * Ty';
					for( indexMatrixI = 0; indexMatrixI < 16; indexMatrixI++ )
					{
						sum = 0;
						for( indexMatrixJ = 0; indexMatrixJ < 16; indexMatrixJ++ )
						{
							sum += Ty[ indexMatrixJ ] * M[ indexMatrixI + M_m * indexMatrixJ ];
						}
						alfaY[ indexMatrixI ] = sum;
					}
					
					// Az elemen belül végimegyünk minden pixelen és feltöltjük Xtr és
					// Ytr mezoket.
					for ( I = Ibegin; I <= Iend; I++ )  //I=IL2I(IL):(IL2I(IL+1)-1)
					{
						for ( J = Jbegin; J <= Jend; J++ )  //J=JL2J(JL):(JL2J(JL+1)-1) 
						{	
							// x és y értéke az elem sarkától számítva. Ez a feltétele
							// annak, hogy egy mátrixszal írjuk le az összes elemet.
							// x=Jp2x(J)-Jline2x(JL); !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
							// y=Ip2y(I)-Iline2y(IL); !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
							x = dx * ( (double)J - 0.5 ) / (double)UpScale - ( Dx * ( (double)JL - HalfJline ) + dx * ( (double)Jimage / 2 - 1 ) );
							y = dy * ( (double)I - 0.5 ) / (double)UpScale - ( Dy * ( (double)IL - HalfIline ) + dy * ( (double)Iimage / 2 - 1 ) );
							
							//if( JL == 0 )
							//	printf( "x:%lf y:%lf\t", x, y );
					
							Xtr = alfaX[ 0  ]                  + alfaX[ 1  ] * x                 + alfaX[ 2  ] * y                 + alfaX[ 3  ] * x * x + 
								  alfaX[ 4  ] * x * y          + alfaX[ 5  ] * y * y             + alfaX[ 6  ] * x * x * x         + alfaX[ 7  ] * x * x * y + 
								  alfaX[ 8  ] * x * y * y      + alfaX[ 9  ] * y * y * y         + alfaX[ 10 ] * x * x * x * y     + alfaX[ 11 ] * x * x * y * y + 
								  alfaX[ 12 ] * x * y * y * y  + alfaX[ 13 ] * x * x * x * y * y + alfaX[ 14 ] * x * x * y * y * y + alfaX[ 15 ] * x * x * x * y * y * y;
							
							Ytr = alfaY[ 0  ]                  + alfaY[ 1  ] * x                 + alfaY[ 2  ] * y                 + alfaY[ 3  ] * x * x + 
								  alfaY[ 4  ] * x * y          + alfaY[ 5  ] * y * y             + alfaY[ 6  ] * x * x * x         + alfaY[ 7  ] * x * x * y + 
								  alfaY[ 8  ] * x * y * y      + alfaY[ 9  ] * y * y * y         + alfaY[ 10 ] * x * x * x * y     + alfaY[ 11 ] * x * x * y * y + 
								  alfaY[ 12 ] * x * y * y * y  + alfaY[ 13 ] * x * x * x * y * y + alfaY[ 14 ] * x * x * y * y * y + alfaY[ 15 ] * x * x * x * y * y * y;
							
							Xwv = Jp2x( J, UpScale, dx ) + Xtr;
							Ywv = Ip2y( I, UpScale, dy ) + Ytr;
							i = y2i( Ywv, dy ); j = x2j( Xwv, dx );
							
							//if( JL == 0 )
							//	printf( "x:%d y:%d\t", i, j );
							
							if ( ( i > 0 ) && ( i <= Iimage ) && ( j > 0 ) && ( j <= Jimage ) )
							{
								iCenter = i2y( i, dy );
								jCenter = j2x( j, dx );
								if ( ( ( iCenter - Ywv ) * ( iCenter - Ywv ) + ( jCenter - Xwv ) * ( jCenter - Xwv ) ) < ( ( iCenter - YwvStored[ i + Jimage * j ] ) * ( iCenter - YwvStored[ i + Jimage * j ] ) + ( jCenter - XwvStored[ i + Jimage * j ] ) * ( jCenter - XwvStored[ i + Jimage * j ] ) ) )
								{
									LinTblX[ i - 1 + Jimage * ( j - 1 ) ] =- Xtr;
									LinTblY[ i - 1 + Jimage * ( j - 1 ) ] =- Ytr;
									XwvStored[ i + Jimage * j ]=Xwv;
									YwvStored[ i + Jimage * j ]=Ywv;
								}
							}
							
						}
					}
				}
			}
		}
	}
	
	for( indexX = 0; indexX < Iimage; indexX++ )
	{
		for( indexY = 0; indexY < Jimage; indexY++ )
		{
			LinTblX[ indexX + indexY * Jimage ] *= 64.0 / dx;
			LinTblY[ indexX + indexY * Jimage ] *= 64.0 / dy;
		}
	}
	
	free( Tx );
	free( Ty );
	free( alfaX );
	free( alfaY );
	free( XwvStored );
	free( YwvStored );
	
	CompleteMatrix( LinTblX, Jimage, Iimage, Iimage );
	CompleteMatrix( LinTblY, Jimage, Iimage, Iimage );
	
	return;
}
