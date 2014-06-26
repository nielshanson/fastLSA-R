#include "pnorm.h"

using namespace std;
double PNorm::norm ( double x ) {
	return( sqrt( ((double)2.0) / ((double)3.14159265358979323846264338328) ) * exp( - pow( x ,(double)2.0) / ((double)2.0) ) ) ;
}

/*
PNorm::PNorm( unsigned int resolution ) {
      
	key.push_back ( 0.0 ) ;
    float x;

	n = resolution;
	for(unsigned  int i = 1 ; i < n ; i++ ) {
		x =   1.0 / ( (float)n * norm( key[i-1] )  ) + key[i-1]  ;
		key.push_back( x ) ;
	}
}
*/

PNorm::PNorm( double resolution , double step , double max )
{
	double total = 0.0 ;
	double place = 0.0 ;
	double next = step ;
	key.push_back ( place ) ;
	value.push_back ( total ) ;
	n = resolution ;
	
	while( place < max )
	{
		place = place + resolution ;
		total = total + resolution * norm( place ) ;
		
		if( place > next )
		{
			next = next + step ;
			key.push_back ( place ) ;
			value.push_back ( total ) ;
		}
	}
}
	
void PNorm::printValues() {

    cout << "Calling print values ";
	for(unsigned  int i = 1 ; i < n ; i++ ) {
//		x= getCumProbability(key[i]) + (float) 0.5/n;
		cout << " i "  << i << "   " << key[i] << "    "  << getCumProbability(key[i]) << endl;
    }

}
double PNorm::getResolution() {

   return n;
}


double PNorm::getCumProbability ( double x ) {
	  int a = 0;
    int b = key.size()-1 ;
    int i = (a + b)/2;
	
    if( x >= key[b] ) return 1.0;
	  while (a < b-1 ) {
          i = (a + b)/2;
          if( x < key[i] ) {
             b = i;
          }
          else if( x > key[i] ) {
             a = i;
          }
          else {
             break;
          }
  	}
	
	return ( (float)value[i] ) ;
}
