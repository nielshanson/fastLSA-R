#ifndef PNORM_H
#define  PNORM_H

#include <vector>
#include <math.h>
#include <iostream>

using namespace std;
class PNorm {
private:
	vector<double> key;
	vector<double> value;
	double n;
	
public:
	double norm ( double x ) ;
	PNorm( double resolution , double step , double max ) ;
	double getCumProbability ( double x ) ;
	double getResolution() ;
	void printValues() ;
	
};


#endif 
