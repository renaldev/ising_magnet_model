#ifndef _ISING_H
#define _ISING_H

#include <fstream>
#include <iostream>

const double J = -1.0;      // AB взаимодействие
const double Jcr = 0;  
const double F = - 0.2;
const double u0 = 1;     // magnetic moment

class Ising
{ 
   private:

      float * sites;         // pointer to an integer array to be allocated dynamically
      float **  sitesB;      // подрешетка A
      float **  sitesA;       // Подрешетка B
      int L;               // linear dimension of lattice
      int num, nam;
   
      float & spin(int index1, int index2, int index3); // retrieve spin of particular site
      float chrome (int i1, int i2, int i3);
      void include (int c=0, char r='r', int k=0);
      void getM(void);              // calculate overall energy
      void getE(void);
      //bool flip(void);             // attempt to flip a random spin

   public:

      Ising(int l=5, double t=10, double b=0, int c=0, char rasp='r', int koef=0);   // constructor with default arguments
      ~Ising()  {  if(sites) delete [] sites;  }    // destructor deallocating memory
      
      
	  double T, M, E, B;   // temperature, magnetization, energy, external magnetic field
	double m_A();
	double m_B();
	void snimok(double Y=0., int R=0);
      void onestep();           // one Monte Carlo step per spin
      bool flip(void);
};

#endif
