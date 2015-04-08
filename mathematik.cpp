#include <stdlib.h>
#include <cmath>
#include "mathematik.h"
#include <limits>

using namespace std;

bool compare_abs (float & a, float & b) {

	if ( abs(abs(a/b)-1.0) <= numeric_limits<float>::epsilon() ) return true; else return false;

}

int rand_r(int max, char c) {

	float b=max/2;
	if (c==2) return  (2*int(drand48()*b)); 
	if (c==1)  if ((2*int(drand48()*b)-1) > 0) return (2*int(drand48()*b)-1); else return 1;
	return int(drand48()*max);
}

int rand_n(int max, char c, int kf) {

	int k=kf;//400
	//if (max==30) k=2000;//3000
	float r1=2.*drand48()-1., r2=2.*drand48()-1.;
	float s=r1*r1+r2*r2;
	int x1, x2;
	if (s>=1) return(rand_n(max,c,k)); else
	{
		r1=r1*sqrt(-2*log(s)/s), r2=r2*sqrt(-2*log(s)/s);
		x1=r1*k+max/4; x2=r2*k+max/4;
		if (c==2) {
		x1=x1*2; x2=x2*2;
		if ((x1<max) && (x1>=0)) return x1;
				else if ((x2<max) && (x2>=0)) return x2;
					else return(rand_n(max,c,k));
					}
		if (c==1) {
		x1=x1*2+1; x2=x2*2+1;
		if (x1<max && x1>=0) return x1;
				else if (x2<max && x2>=0) return x2;
					else return(rand_n(max,c,k));
					}
		x1=r1*k+max/2; x2=r2*k+max/2;
		if (x1<max && x1>=0) return x1;
				else if (x2<max && x2>=0) return x2;
					else return(rand_n(max,c,k));			
		
	}
	return 0;
}
