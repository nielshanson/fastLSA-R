#include <omp.h>
#include <R.h>
#include "pnorm.h"

//int tot; // grand total of matches over all threads

// processes row pairs (i, i+1), (i, i+2), ...

// int procpairs(int i, int *m, int n) {
//   int j,k,sum=0;
//   for (j = i+1; j < n; j++) {
//     for (k = 0; k < n; k++) {
//       //find m[i][k] * m[j][k] but remember R uses col-major order
//       sum += m[n*k+i] * m[n*k+j];
//     }
//   }
//   return sum;
// }

// function to standardize each time-series row
void standardize_rows(double **mat2, int m, int n) {
	// mat - point to inputmatrix
	// m - number of rows (time-series)
	// n - number of columns (time-steps)
	
	for (int i=0; i < m; i++) { 
		  // get each rowsum
  	  double row_mean = 0;
    	for (int j = 0; j < n; j++){
    	    row_mean += mat2[i][j];
      }
		  row_mean = row_mean / n;
		  
		  // calculate standard deviation
		  float sd = 0.0 ;
		  for( int j = 0 ; j < n ; j++){
			    sd = sd + pow( mat2[i][j] - row_mean , 2 ) ;
		   }
		   sd = sqrt( sd / ( n - 1 ) ) ;
		   
		   // ask Evan about this standardization
		   for( int j = 0 ; j < n ; j++ ){
			     mat2[i][j] = ( mat2[i][j] - row_mean ) / sd  ;
		   }
	}
}

// given two row pointers to matrix m, calculate the LSA statistic
// of the two rows
double lsa (double *x , double *y , int n , int d , int *D) {
	double Pmax = 0.0 ;
	double Nmax = 0.0 ;
	
	// P[i][0] , N[i][0]
	for( int i = 0 ; i <= d ; i ++ ) {
		double Ptemp = 0.0 ;
		double Ntemp = 0.0 ;
		
		for( int k = 0 ; k + i < n ; k ++ ) {
			double temp = Ptemp + x[i+k] * y[k] ;
			if(temp > 0.0 ) {
				Ptemp = temp ;
			}
			temp = Ntemp - x[i+k] * y[k] ;
			if( temp > 0.0 ) {
				Ntemp = temp ;
			}
			if( Ptemp > Pmax ) {
				Pmax = Ptemp;
				if( Pmax > Nmax ) {
					*D = i;
				}
			}
			if( Ntemp > Nmax ) {
				Nmax = Ntemp ;
				if( Nmax > Pmax ) {
					*D = i ;
				}
			}
		}
	}
		
	// P[0][j] , N[0][j]
	for( int j = 1 ; j <= d ; j ++ ) {
		double Ptemp = 0.0;
		double Ntemp = 0.0;
		
		for( int k = 0 ; k + j < n ; k ++ ) {
			float temp = Ptemp + x[k] * y[j+k];
			if( temp > 0.0 ) {
				Ptemp = temp;
			}
			temp = Ntemp - x[k] * y[j+k];
			if( temp > 0.0 ) {
				Ntemp = temp;
			}
			if( Ptemp > Pmax ) {
				Pmax = Ptemp;
				*D = -j;
			}
			if( Ntemp > Nmax ) {
				Nmax = Ntemp;
				*D = -j;
			}
		}
	}
	
	if( Pmax > Nmax ) {
		return( Pmax / n );
	} else {
		return ( - Nmax / n );
	}
	
	return( 0.0 );
}

double lsaPBound( double x , int n , int minij , int d , double variance , PNorm * p )
{
	if ( x < 0 )
		x = -x ;
	return(  2.0 * ( n*n - ( n - d - 1 ) * ( n - d ) ) * ( 1 - p->getCumProbability( x * sqrt( n / variance ) ) )  ) ;
}


extern "C" void lsarun(double *mat, int *m, int *n, int *d,
                       double *alpha, double *minLSA, int *rez, 
                       double *lsa_result, double *pval_result, int *nthreads) {
  
  int mval = *m;
  int nval = *n;
  int dval = *d;
  double alphaval = *alpha;
  double minLSAval = *minLSA;
  int rezval = *rez;
  int nthreadsval = *nthreads;
  
  double *R = mat;
  // create index matrix for C: mat2[row][col]
  double **mat2 = new double*[mval];
  for (int i=0; i < mval; i++) {
    mat2[i] = R + (i*(nval));
  }
  
  standardize_rows(mat2, mval, nval);
  
  
  #pragma omp parallel
  {
    int i=0, j=0, D;
    int me = omp_get_thread_num();
    int nth = omp_get_num_threads();
    double mylsa = 0.0;
    double mypvaluebound = 0.0;
    PNorm *p = new PNorm( ((double)1.0)/((double)rezval) , 0.01 , 5.0 ) ;
    for (i = me; i < mval; i += nth) {
      // this thread's series
      for (j = i+1; j < mval; j++) {
        // all possible pairs
        if(i==j){
          lsa_result[(i*nval) + j] = 1.0;
        } else {
          mylsa = lsa(&mat2[i][0], &mat2[j][0], nval, dval, &D);
          mypvaluebound = lsaPBound( mylsa , nval , 1 , dval , 1.0 , p ) ;
          // lsa_result[(i*nval) + j] = mylsa;
          // pval_result[(i*nval) + j] = mypvaluebound;
          
          if( mypvaluebound < alphaval && ( mylsa > minLSAval || -mylsa > minLSAval ) ) {
            // approved, report result
            lsa_result[(i*nval) + j] = mylsa;
            pval_result[(i*nval) + j] = mypvaluebound;
          }
          
        }
      }
    }
    // #pragma omp atomic
    // tot += mysum;
    // *nthreads = nth;
  }
  
  
  //result = mat;
  
  // calculate row sums
  // for (int i = 0; i < mval; i++) {
  // 	  double row_sum = 0;
  // 	  for (int j = 0; j < nval; j++){
  // 	      row_sum += mat[(i*nval)+j];
  // 	  }
  // 	  result[i] = row_sum;
  // }
  
	//   for( int i = 0 ; i < n ; i ++ )
	//   {
	// data.setRow( i , standardize1(data.get(i)) ) ;
	//   }
	//   
	//   #pragma omp parallel
	//   {
	//     int i, mysum = 0;
	//     int me = omp_get_thread_num();
	//     int nth = omp_get_num_threads();
	//     for (i = me; i < nval; i += nth) {
	//       mysum += procpairs(i,m,nval);
	//     }
	//     #pragma omp atomic
	//     tot += mysum;
	//     *num_threads = nth;
	//   }
	//   int divisor = nval * (nval-1) / 2;
	//   *mlmean = ((float) tot) / divisor;
}