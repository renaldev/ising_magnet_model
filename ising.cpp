#include "ising.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cmath>
#include <string.h>
#include <sstream>
#include "mathematik.h"

using namespace std;

//////////////////////////////////////конструктор модели///////////////////////////////////////////

Ising::Ising(int l, double t, double b, int c, char rasp, int koef)  // формирование куба
{ 
   L = l; T = t; B = b; M = 0; // размер, температура, поле, намагниченность 

   int lat = L*L*L;
   sites = new float[lat]; // куб спинов
   lat = lat/2;
   sitesB = new float* [lat];
   sitesA = new float* [lat];
   num = 0;
   nam = 0;
      
   for(int x=0; x<L; x++)
     for(int y=0; y<L; y++)
        for(int z=0; z<L; z++)
   	    	
   	    	if (z%2) {
         		if (x%2)
         			if  (y%2) {spin(x,y,z)= - 1; sitesA[nam]=&spin(x,y,z);  nam++;}
         			else {spin(x,y,z)= 0;  sitesB[num]=&spin(x,y,z);  num++; }   //-0.8
				else
					if  (y%2)  { spin(x,y,z)= 0; sitesB[num]=&spin(x,y,z);  num++; } //-0.6
					else {spin(x,y,z)= - 1; sitesA[nam]=&spin(x,y,z);  nam++;}
				}
			else {
				if (x%2)
					if  (y%2)  {  spin(x,y,z)= 0; sitesB[num]=&spin(x,y,z);  num++;} //-0.6
					else {spin(x,y,z)= - 1; sitesA[nam]=&spin(x,y,z);  nam++;}
				else
					if  (y%2) {spin(x,y,z)= - 1; sitesA[nam]=&spin(x,y,z);  nam++;}
					else { spin(x,y,z)= 0; sitesB[num]=&spin(x,y,z);  num++; }  //-0.8
				}	
	cout << L;
	include (c,rasp, koef);
	getM();
	getE();
					
}

void Ising::include (int c, char r, int k)  //формирование состава куба
 {
 	int conc = 0, i = 0;
 	
 	while (conc < 6751) { //  случайная растановка спинов хрома
 	//i = rand_n (num, 0, 1350);
 	i = rand_r (num);
 	//cout << i << "  ";
 	if ( abs(*sitesB[i]) < 0.5) { *sitesB[i] = 1.5; conc++; } //-1.2
 	}
 	
 	for (int i=0; i<num; i++) if ( abs(*sitesB[i]) < 0.3) *sitesB[i] = 1.5; // растановка спинов железа
 	
 	conc = 0, i = 0;
 	while (conc < c) { //случайная растановка спинов магния
 	//i = rand_n (num, 0, 1700);
 	if (r=='r') i = rand_r (num); else i = rand_n (num, 0, k); //k=1700
 	if ( abs(*sitesB[i]) > 0.8) { *sitesB[i] = 0; conc++; }
 	} 
 }

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


inline float & Ising::spin(int index1, int index2, int index3) // запрос спина по координатам
{  
   if(index1<0) index1 += L;      // these operations guarantee periodic boundary conditions
   if(index1>=L) index1 -= L;
   if(index2<0) index2 += L;
   if(index2>=L) index2 -= L;
   if(index3<0) index3 += L;
   if(index3>=L) index3 -= L;

   int index = index1*L*L+ index2*L + index3;
   return sites[index];
}

inline void Ising::getM(void)  //расчет полный магнитный момент
{  
   for(int i=0; i<L; i++)
      for(int j=0; j<L; j++)
         for(int k=0; k<L; k++)  	   
         		M+=spin(i,j,k);
}


inline void Ising::getE(void)
{  int nn_count = 0;    // nearest neighbour count
	int Fnn_count = 0;
	float cr_spinF=1.5;

   for(int i=0; i<L; i++)
     for(int j=0; j<L; j++)
       for(int k=0; k<L; k++)
	{  
		nn_count += ( J*spin(i, j,k)*( spin(i+1, j, k) + spin(i, j+1,k) + spin(i, j,k+1) ) ); 
		if (compare_abs(spin(i, j, k), cr_spinF)) {
			Fnn_count += ( F*spin(i, j,k)*( spin(i-1,j-1,k)+spin(i+1,j+1,k)+spin(i-1,j-1,k-1)+spin(i+1,j+1,k-1)+spin(i-1,j-1,k+1)+spin(i+1,j+1,k+1)+spin(i-1,j+1,k)+spin(i+1,j-1,k)+spin(i-1,j+1,k-1)+spin(i+1,j-1,k-1)+spin(i-1,j+1,k+1)+spin(i+1,j-1,k+1)+spin(i-1,j,k-1)+spin(i+1,j,k-1)+spin(i-1,j,k+1)+spin(i+1,j,k+1)+spin(i,j+1,k-1)+spin(i,j-1,k-1)+spin(i,j+1,k+1)+spin(i,j-1,k+1) ) );
		}


	}


	E = -nn_count - Fnn_count/2 - u0*B*M;
}






void Ising::onestep() //шаг MCS
{   int size = L*L*L;

    for(int i=0; i<size; i++)
       flip();

}

bool Ising::flip(void)
{
    int index1 = int(drand48()*L);  int index2 = int(drand48()*L); int index3 = int(drand48()*L);
    
    float cr_spin=0.33333333; //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

    int delta_M = -2*spin(index1, index2, index3);      
    float nn_count = 0;
    if ( compare_abs(spin(index1, index2, index3), cr_spin) ) nn_count = Jcr * (spin(index1-1, index2, index3) + spin(index1+1, index2, index3)
                         				+ spin(index1, index2-1, index3) + spin(index1, index2+1, index3)
                         				+ spin(index1, index2, index3-1) + spin(index1, index2, index3+1)); 	else	{
                         				
			if  ( compare_abs(spin(index1-1, index2, index3), cr_spin) )  nn_count += Jcr * spin(index1 - 1, index2, index3);
			else nn_count += J*spin(index1 - 1, index2, index3);
			
			if  ( compare_abs(spin(index1+1, index2, index3), cr_spin) )  nn_count += Jcr * spin(index1 + 1, index2, index3);
			else nn_count += J*spin(index1 + 1, index2, index3);
			
			if  ( compare_abs(spin(index1, index2-1, index3), cr_spin) ) nn_count += Jcr * spin(index1, index2 - 1, index3);
			else nn_count += J*spin(index1, index2 - 1, index3);
			
			if  ( compare_abs(spin(index1, index2+1, index3), cr_spin) )  nn_count += Jcr * spin(index1, index2 + 1, index3);
			else nn_count += J*spin(index1, index2 + 1, index3);
			
			if  ( compare_abs(spin(index1, index2, index3-1), cr_spin) )  nn_count += Jcr * spin(index1, index2, index3 - 1);
			else nn_count += J*spin(index1, index2, index3 - 1);
			
			if  ( compare_abs(spin(index1, index2, index3+1), cr_spin) )  nn_count += Jcr * spin(index1, index2, index3 + 1);
			else nn_count += J*spin(index1, index2, index3 + 1);
			};
                   
    float cr_spinF=1.5;                       
    int delta_F = 0; if ( compare_abs(spin(index1, index2, index3), cr_spinF) )  delta_F=chrome(index1, index2, index3);
	
	delta_F = -2*spin(index1, index2, index3)*delta_F; //энергия фрустрации
	
    nn_count = -2*spin(index1, index2, index3)*nn_count; //обменная энергия

    double delta_E = -1.0*nn_count -F*delta_F -u0*B*delta_M; //изменение полной энергии

    if (delta_E<=0 || drand48()<exp(-delta_E/T))
	{  spin(index1, index2, index3) *= -1.0;     // spin(int, int) must return int & to modify its content
           M += delta_M;
           E += delta_E;
           return 1;
        }
    return 0;
}

float Ising::chrome (int i1, int i2, int i3) 
{
	float E=0;
	if 	(compare_abs(spin(i1,i2,i3), spin(i1-1,i2-1,i3)))  E=E+spin(i1-1,i2-1,i3);
	if 	(compare_abs(spin(i1,i2,i3), spin(i1+1,i2+1,i3)))  E=E+spin(i1+1,i2+1,i3);
	if 	(compare_abs(spin(i1,i2,i3), spin(i1-1,i2-1,i3-1)))  E=E+spin(i1-1,i2-1,i3-1);
	if 	(compare_abs(spin(i1,i2,i3), spin(i1+1,i2+1,i3-1)))  E=E+spin(i1+1,i2+1,i3-1);
	if 	(compare_abs(spin(i1,i2,i3), spin(i1-1,i2-1,i3+1)))  E=E+spin(i1-1,i2-1,i3+1);
	if 	(compare_abs(spin(i1,i2,i3), spin(i1+1,i2+1,i3+1)))  E=E+spin(i1+1,i2+1,i3+1);
	
	if 	(compare_abs(spin(i1,i2,i3), spin(i1-1,i2+1,i3)))  E=E+spin(i1-1,i2+1,i3);
	if 	(compare_abs(spin(i1,i2,i3), spin(i1+1,i2-1,i3)))  E=E+spin(i1+1,i2-1,i3);
	if 	(compare_abs(spin(i1,i2,i3), spin(i1-1,i2+1,i3-1)))  E=E+spin(i1-1,i2+1,i3-1);
	if 	(compare_abs(spin(i1,i2,i3), spin(i1+1,i2-1,i3-1)))  E=E+spin(i1+1,i2-1,i3-1);
	if 	(compare_abs(spin(i1,i2,i3), spin(i1-1,i2+1,i3+1)))  E=E+spin(i1-1,i2+1,i3+1);
	if 	(compare_abs(spin(i1,i2,i3), spin(i1+1,i2-1,i3+1)))  E=E+spin(i1+1,i2-1,i3+1);
	
	if 	(compare_abs(spin(i1,i2,i3), spin(i1-1,i2,i3-1)))  E=E+spin(i1-1,i2,i3-1);
	if 	(compare_abs(spin(i1,i2,i3), spin(i1+1,i2,i3-1)))  E=E+spin(i1+1,i2,i3-1);
	if 	(compare_abs(spin(i1,i2,i3), spin(i1-1,i2,i3+1)))  E=E+spin(i1-1,i2,i3+1);
	if 	(compare_abs(spin(i1,i2,i3), spin(i1+1,i2,i3+1)))  E=E+spin(i1+1,i2,i3+1);
	
	if 	(compare_abs(spin(i1,i2,i3), spin(i1,i2+1,i3-1)))  E=E+spin(i1,i2+1,i3-1);
	if 	(compare_abs(spin(i1,i2,i3), spin(i1,i2-1,i3-1)))  E=E+spin(i1,i2-1,i3-1);
	if 	(compare_abs(spin(i1,i2,i3), spin(i1,i2+1,i3+1)))  E=E+spin(i1,i2+1,i3+1);
	if 	(compare_abs(spin(i1,i2,i3), spin(i1,i2-1,i3+1)))  E=E+spin(i1,i2-1,i3+1);
	
	return E;
}

 double Ising::m_A()
 {
 		double mA=0;
 	     for (int i=0; i<nam; i++)    mA+=(*sitesA[i]);
     		return mA;
 }
 
  double Ising::m_B()
 {
 	double mB=0;
 	     for (int i=0; i<num; i++)    mB+=(*sitesB[i]);
     		return mB;
 }
 
 void Ising::snimok(double Y, int R) {

 float no=0;
 stringstream name; name <<"slise"<<Y<<"_"<<R<<".inc";
 FILE *fn = fopen (name.str().c_str(), "w+");
   for ( int i = 0; i < 30; i++ )
      for ( int j = 0; j < 30; j++ )
          for ( int k = 0; k < 30; k++ )
            {
            if ( k <= 15 || j <= -1 ) {
            if ( abs(spin(i, j, k))<0.5  ) { fprintf (fn, "sphere{<%d,%d,%d>,0.7 texture{SN}}\n",
                      i, j, k); continue;}
            
            fprintf (fn, "sphere{<%d,%d,%d>,0.7 texture{%s}}\n",
                      i, j, k, spin(i, j, k) < 0 ? "SD" : "SU");
                      }
             }
             
fclose(fn);
             
}

 
 
 
