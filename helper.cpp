#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double meanC(NumericVector x) {
  int n = x.size();
  double total = 0;

  for(int i = 0; i < n; ++i) {
    total += x[i] / n;
  }
  return total;
}

// [[Rcpp::export]]
double lsa ( NumericVector x , NumericVector y , int n , int d , int D )
{
	double Pmax = 0.0 ;
	double Nmax = 0.0 ;
	
	// P[i][0] , N[i][0]
	for( int i = 0 ; i <= d ; i ++ )
	{
		double Ptemp = 0.0 ;
		double Ntemp = 0.0 ;
		
		for( int k = 0 ; k + i < n ; k ++ )
		{
			if( Ptemp + x[i+k] * y[k] > 0.0 )
			{
				Ptemp = Ptemp + x[i+k] * y[k] ;
			}
			if( Ntemp - x[i+k] * y[k] > 0.0 )
			{
				Ntemp = Ntemp - x[i+k] * y[k] ;
			}
			if( Ptemp > Pmax )
			{
				Pmax = Ptemp ;
				if( Pmax > Nmax )
				{
					D = i ;
				}
			}
			if( Ntemp > Nmax )
			{
				Nmax = Ntemp ;
				if( Nmax > Pmax )
				{
					D = i ;
				}
			}
		}
	}
	
	// P[0][j] , N[0][j]
	for( int j = 1 ; j <= d ; j ++ )
	{
		double Ptemp = 0.0 ;
		double Ntemp = 0.0 ;
		
		for( int k = 0 ; k + j < n ; k ++ )
		{
			if( Ptemp + x[k] * y[j+k] > 0.0 )
			{
				Ptemp = Ptemp + x[k] * y[j+k] ;
			}
			if( Ntemp - x[k] * y[j+k] > 0.0 )
			{
				Ntemp = Ntemp - x[k] * y[j+k] ;
			}
			if( Ptemp > Pmax )
			{
				Pmax = Ptemp ;
				D = -j ;
			}
			if( Ntemp > Nmax )
			{
				Nmax = Ntemp ;
				D = -j ;
			}
		}
	}
	
	if( Pmax > Nmax )
	{
		return( Pmax / n ) ;
	}
	else 
	{
		return ( - Nmax / n ) ;
	}
	
	return( 0.0 ) ;
}

// [[Rcpp::export]]
NumericVector standardize1 ( NumericVector v )
{
	int n = v.size();
	double mean = 0.0 ;
	for( int i = 0 ; i < n ; i ++ )
	{
		mean = mean +  v[i] ;
	}
	mean = mean / ((double) n ) ;
	double sd = 0.0 ;
	for( int i = 0 ; i < n ; i ++ )
	{
		sd = sd + pow( v[i] - mean , 2 ) ;
	}
	sd = sqrt( sd / ( n - 1 ) ) ;
	for( int i = 0 ; i < n ; i ++ )
	{
		v[i] = ( v[i] - mean ) / sd  ;
	}
	return( v );
}
