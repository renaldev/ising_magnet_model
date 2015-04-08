#include "ising.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cmath>
#include <sstream>

using namespace std;

void start(int l, int c, int k) {

		int L=l;  int size = L*L*L;
    double B= -0.25;  //поле 
    double T=5;   //максимальная температура
    double M=0; double mA=0; double mB=0;
    char rasp='r';
    int koef=k;
    
    double Msq=0; //??????????????
    
    const int conf = 20; //колличество конфигураций
    Ising experiment[conf] = Ising( L, T, B, c, rasp); //массив конфигураций
		int count=0; //число шагов
		int n=0;
		
    srand48(1234567); stringstream name; name <<"kub="<<L<<"_mg="<<c<<"B="<<B<<"frust_"<<rasp<<koef<<".dat"; ofstream output(name.str().c_str());  //!!!!!!!fgdfhgdfdhfg!!!!!!!!!!!!!!!!!!!!!!!fgdgdg!
		
		for ( double t = T; t > 0.01;t-=0.05 ) {

		  cout << "conc: " << c << "| temp: " << t << "| lat: " << L << endl; //вывод консоль!!!!!!!!!!!!!!!!!!!!!!
			double M=0; int n=0; double Msq=0; double Mna=0;
			
      for (int i=0;i<conf;i++) {	
			
				cout << flush;	experiment[i].T=t;
				
        //if (t<3.5) experiment[i].B= 0;
			
				count = 0;
				do	{  experiment[i].onestep();     // these configurations are discarded to acheive equilibrium
       			} while(++count<100);
       		
       			count =	0;
       			do 	{ experiment[i].onestep();
       			M   = M+fabs(experiment[i].M);
       			mA   = mA+fabs(experiment[i].m_A());
       			mB   = mB+fabs(experiment[i].m_B());
       			Mna   = Mna+(experiment[i].M);
       			Msq = Msq + experiment[i].M * experiment[i].M;
       			n++;      
     			} while(++count<10);
			}
			
			output<<"\n";
    		output.width(15); output<<t;
    		output.width(15); output << M/(size*n);
    		output.width(15); output<< (Msq/(size*n) - ( (M/(size*n)) * (M/(size*n)) )*size)/t;
    		output.width(15); output << Mna/(size*n);
    		output.width(15); output << 2*mA/(size*n);
    		output.width(15); output << 2*mB/(size*n);
    		//if (!(int(t*10)%10)) output<<flush;
    		//if (!(int(t*10)%10)) experiment[5].snimok(t,c);
		}
}

int main(void)
{  
 //for (int c=0; c<7000; c+=675) 	
 start(30, 3700, 0); //L,  c
 
 
		return 0;
}
