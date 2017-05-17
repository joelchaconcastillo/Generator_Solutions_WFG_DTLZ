/*
   Autor: Joel Chacón Castillo
   Date: 10-may-17

  Generates weights according to a uniform design of mixtures using the Hammersley lo-discrepancy sequence generator. 
This method were proposed by Berenguer and Coello Coello (2015).
References:
    Tan Y., Y. Jiao, H. Li, and X. Wang (2013).  "MOEA/D + uniform design:
        A new version of MOEA/D for optimization problems with many
        objectives."  Computers & Operations Research, 40, 1648-1660.
    Berenguer, J.A.M. and C.A. Coello Coello (2015).  "Evolutionary Many-
        Objective Optimization Based on Kuhn-Munkres' Algorithm."  Evolutionary
        Multi-Criterion Optimization: 8th International Conference, pp. 3-17.
 
*/
#ifndef UNIFORMGENERATOR_HPP
#define UNIFORMGENERATOR_HPP
#include <iostream>
#include <bitset>
#include <vector>
#include <cmath>
#include "WFG1.h"
#include "WFG2.h"
#include "WFG3.h"
#include "WFG4.h"
#include "WFG5.h"
#include "WFG6.h"
#include "WFG7.h"
#include "WFG8.h"
#include "WFG9.h"
#define TOL 1e-5
using namespace std;
class UniformGenerator
{
   public:
      UniformGenerator();
      vector<vector<double> > GenerateQuasiRandom(int Objectives, int N);
      vector< vector< double> > Hyperplane(int N, int dim);
      vector< vector<double> > GenerateConcaveRandomPoints(int NPoints, int M, vector<double> &scale);
      void Filter(vector<vector<double> > &N, vector<vector<double > > &R);
      vector< vector<double> > GenerateSolutionWFG1(int N, int Obj ,int k, int l);
      vector< vector<double> > GenerateSolutionWFG2(int N, int Obj ,int k, int l);
      vector< vector<double> > GenerateSolutionWFG3(int N, int Obj ,int k, int l);
      vector< vector<double> > GenerateSolutionWFG4(int N, int Obj ,int k, int l);
      vector< vector<double> > GenerateSolutionWFG5(int N, int Obj ,int k, int l);
      vector< vector<double> > GenerateSolutionWFG6(int N, int Obj ,int k, int l);
      vector< vector<double> > GenerateSolutionWFG7(int N, int Obj ,int k, int l);
      vector< vector<double> > GenerateSolutionWFG8(int N, int Obj ,int k, int l);
      vector< vector<double> > GenerateSolutionWFG9(int N, int Obj ,int k, int l);

      vector< vector<double> > GenerateSolutionDTLZ1(int N, int Obj);
      vector< vector<double> > GenerateSolutionDTLZ2_to_DTLZ4(int N, int Obj);
      vector< vector<double> > GenerateSolutionDTLZ5_to_DTLZ6(int N, int Obj);
      vector< vector<double> > GenerateSolutionDTLZ7(int N, int Obj);
      inline void setDelta(double Delta){this->Delta = Delta;}
      
   private:
      double Delta = 0.001;
      vector <int> generateFirstPrimes(int k);
      void GenerateDegenerateRandomPoints(int NPoints, int M);
      double linear(vector<double> &x, int m, int M);
      double concave(vector<double> &x, int m, int M);
      double EuclideanDistance(vector<double> &x, vector<double> &y);
      bool Dominate(vector<double> &A, vector<double> &B);
      void Non_Dominance(vector<vector<double> > &S);
};
#endif

UniformGenerator::UniformGenerator()
{

}
/**
   Create a list of firt primes with the sieve 
**/
vector <int> UniformGenerator::generateFirstPrimes(int k)
{
  vector<int> primes;
   int count = 0;
   long long _sieve_size = 1000001;
   bitset<10000010> bs;
   bs.set(); //set bits to 1
   bs[0] = bs[1] = 0;
   for(long long i = 2; i <=_sieve_size; i++)
   {
      if(bs[i])
	{
	   for(long long j = i*i; j <= _sieve_size; j +=i)bs[j]=0;
	   primes.push_back((int) i);
	   count++;
	}
	if(count > k)break;
   }
   return primes;
}
/**
  Here the hammersley method is implemented....
**/
vector<vector<double> > UniformGenerator::GenerateQuasiRandom(int Objectives, int N)
{
   vector<vector<double > > designs;
   vector<int> primes = generateFirstPrimes(Objectives);


 //     vector<double> design(Objectives, 0.0);
 //     designs.push_back(design);
   for(int i = 0 ; i < N; i++)
   {
      vector<double> design(Objectives, 0.0);
      //design[0] = (2.0*(i+1.0) - 1.0) / (2.0*N);
      design[0] = (double)( (i) % (N+1)) / (double) N;
//      design[0] = (double) (i)/N;//(double)( (i+1) % (N+1)) / (double) (N);
      for(int j = 1; j < Objectives; j++)
      {
	   double f = 1.0/primes[j-1];
	   int d = i;
	   design[j] = 0.0;
	   while(d > 0)
	   {
		design[j] += f*(d % primes[j-1]);
		d =   floor((double)d / primes[j-1]); // Here the floor is considered as part of the C++ language.....
		f = f / primes[j-1];
	   }
      }
	designs.push_back(design);
   }
//for(int j = 0; j < designs.size(); j++)
//{
//for(int i = 0; i < Objectives; i++)
//   cout << designs[j][i] << " ";
//cout << endl;
//}
 return designs;
 // The designs are transformed to weight vectors....
}
vector< vector< double> > UniformGenerator::Hyperplane(int N, int dim)
{
  vector< vector<double> > X = GenerateQuasiRandom(dim, N);
  vector< vector< double > > referenceSet;

  for(int k = 0; k < X.size(); k++)
  {
     vector<double> reference;
	
     int m =  X[k].size();
     for(int i = 1; i <= m; i++)
     {
	double fki = 1;//pow(X[k][0], 1.0/m) ;
	if(i != m)
	 fki *= (1.0 - pow(X[k][i-1], 1.0/(m-i))  );

	for(int j=1; j <= i-1; j++)
	{
	 fki *= pow(X[k][j-1],1.0/(m - j ));
	}

        reference.push_back(fki/2.0);
     }

     referenceSet.push_back(reference);
	
  }

   return referenceSet;
}
double UniformGenerator::EuclideanDistance(vector<double> &x, vector<double> &y)
{
	double Sum=0;
	for(int i = 0; i < x.size(); i++)
	Sum+= (x[i]-y[i])*(x[i]-y[i]);
	return sqrt(Sum);
}
double UniformGenerator::concave(vector<double> &x, int m, int M)
{
   double result =1.0;
   for(int i = 1; i <= M-m; i++)
	result *= sin(x[i-1] * M_PI/2.0 );
   if( m != 1)
 	result *= cos( x[M-m]*M_PI/2.0 );
   return result;

}
double UniformGenerator::linear(vector<double> &x, int m, int M)
{
   double result = 1.0;
   for(int i = 1; i <= M-m; i++)
   {
	result *= x[i-1];
   }
   if( m!=1)
   result *= 1 - x[M-m];
   return result;
}
void UniformGenerator::Filter(vector<vector<double> > &N, vector<vector<double > > &R)
{
   double Delta = this->Delta;
  
   int Size = N.size(); 
   for(int i = 0; i < Size; i++)
   {
	int flag = 0;
      for(int j = i+1 ; j < Size; j++)
//	for(int j = 0; j < N.size(); j++)
	{
///		if(i==j) continue;
		if(EuclideanDistance(N[j], N[i]) < Delta)
		{
		   N.erase(N.begin() + j);
		    Size--;
//		   flag=1;
		   //break;
		}		
	}
//	if(flag) i=0;
   }
	R = N;
}
vector< vector<double> > UniformGenerator::GenerateConcaveRandomPoints(int NPoints, int M, vector<double> &scale)
{
  
  vector< vector<double> > x = GenerateQuasiRandom(M, NPoints);
   vector<double> y(M, 0.0);//, x(M, 0.0);
   vector<vector<double> > Nurb;
   for(int i  = 0; i <  NPoints; i++)
   {
//	for(int m = 1; m <=M; m++)
//	x[m-1] = (double)rand()/RAND_MAX;
	for(int m = 1; m <=M; m++)
	{
	   //y[m-1] = concave(x, m, M);
	   y[m-1] =  concave(x[i], m, M);
	}
	Nurb.push_back(y);
   }
  vector<vector<double> >out;
  Filter(Nurb, out);

   return out;
}
void UniformGenerator::GenerateDegenerateRandomPoints(int NPoints, int M)
{
   vector<double> y(M, 0.0), x(M, 0.0);
   vector<vector<double> > Nurb;
   for(int i  = 0; i <  NPoints; i++)
   {
//	for(int m = 1; m <=M; m++)
	//x[m-1] = (double)rand()/RAND_MAX;
	x[0] = (double)rand()/RAND_MAX;

	for(int m = 1; m <=M; m++)
	{
	   y[m-1] =  linear(x, m, M)*2*m;
	}
	Nurb.push_back(y);
   }
  vector<vector<double> >out;
  Filter(Nurb, out);
  for(int i = 0; i < out.size(); i++)
  {
     for(int j = 0; j < M; j++)
	{
		cout << out[i][j] << " ";
	}
	cout << endl;
  }
}

vector< vector<double> > UniformGenerator::GenerateSolutionWFG1(int N, int Obj ,int k, int l)
{

    int a =Obj-1, b = 1;// Obj - (Obj%2);
//		k = 2 * (nObj - 1);

    vector< vector<double> > data1 = GenerateQuasiRandom(k, N); 
    vector< vector<double> > data2 = GenerateQuasiRandom(a, N); 
    vector< vector<double> > Points; 
   
    WFG1 RealObj, LightObj;

    RealObj.init(Obj, k, l); 
    LightObj.init(Obj, a, b ); //Se inicializa con el mínimo número de elementos...
     
    for(int i = 0; i < N; i++)
    {
	vector<double> fReal(Obj, 0.0), fLight(Obj,0.0);
	vector<double> Real = RealObj.WFG_1_random_soln(k, l, data1[i]);

	vector<double> Light = LightObj.WFG_1_random_soln(a, b, data2[i]);
	RealObj.evaluate(Real, fReal);
	LightObj.evaluate(Light, fLight);

	vector<double> tmp;
	for(int d = 0; d < Obj; d++)
          tmp.push_back( RealObj.getXM()-LightObj.getXM() +  fLight[d]);

	Points.push_back(tmp);
    }
	if(Delta < TOL) return Points;	
     vector<vector<double> >out;
     Filter(Points, out);
	

   return out;
}
vector< vector<double> > UniformGenerator::GenerateSolutionWFG2(int N, int Obj ,int k, int l)
{

    int a =Obj-1, b = 2;// Obj - (Obj%2);
//		k = 2 * (nObj - 1);

    vector< vector<double> > data1 = GenerateQuasiRandom(k, N); 
    vector< vector<double> > data2 = GenerateQuasiRandom(a, N); 
    vector< vector<double> > Points; 
    
    WFG2 RealObj, LightObj;

    RealObj.init(Obj, k, l); 
    LightObj.init(Obj, a, b ); //Se inicializa con el mínimo número de elementos...
     
    for(int i = 0; i < N; i++)
    {
	vector<double> fReal(Obj, 0.0), fLight(Obj,0.0);
	vector<double> Real = RealObj.WFG_2_random_soln(k, l, data1[i]);

	vector<double> Light = LightObj.WFG_2_random_soln(a, b, data2[i]);

	RealObj.evaluate(Real, fReal);
	LightObj.evaluate(Light, fLight);

	vector<double> tmp;
	for(int d = 0; d < Obj; d++)
          tmp.push_back( RealObj.getXM()-LightObj.getXM() +  fLight[d]);
	Points.push_back(tmp);
    }

     Non_Dominance(Points);


	if(Delta < TOL) return Points;	
    //limpiar 
     vector<vector<double> >out;
     Filter(Points, out);
	

   return out;
}
vector< vector<double> > UniformGenerator::GenerateSolutionWFG3(int N, int Obj ,int k, int l)
{

    int a =Obj-1, b = 2;// Obj - (Obj%2);
//		k = 2 * (nObj - 1);

    vector< vector<double> > data1 = GenerateQuasiRandom(k, N); 
    vector< vector<double> > data2 = GenerateQuasiRandom(a, N); 
    vector< vector<double> > Points; 
    
    WFG3 RealObj, LightObj;

    RealObj.init(Obj, k, l); 
    LightObj.init(Obj, a, b ); //Se inicializa con el mínimo número de elementos...
     
    for(int i = 0; i < N; i++)
    {
	vector<double> fReal(Obj, 0.0), fLight(Obj,0.0);
	vector<double> Real = RealObj.WFG_3_random_soln(k, l, data1[i]);

	vector<double> Light = LightObj.WFG_3_random_soln(a, b, data2[i]);

	RealObj.evaluate(Real, fReal);
	LightObj.evaluate(Light, fLight);

	vector<double> tmp;
	for(int d = 0; d < Obj; d++)
          tmp.push_back( RealObj.getXM()-LightObj.getXM() +  fLight[d]);
	Points.push_back(tmp);
    }
	if(Delta < TOL) return Points;	
    //limpiar 
     vector<vector<double> >out;
     Filter(Points, out);
	

   return out;
}
vector< vector<double> > UniformGenerator::GenerateSolutionWFG4(int N, int Obj ,int k, int l)
{

    vector< vector<double> > data1 = GenerateQuasiRandom(k, N); 
    vector< vector<double> > data2 = GenerateQuasiRandom(Obj, N); 
    vector< vector<double> > Points; 
    
    WFG4 RealObj, LightObj;

    RealObj.init(Obj, k, l); 
 //   LightObj.init(Obj, a, b ); //Se inicializa con el mínimo número de elementos...
     
    for(int i = 0; i < N; i++)
    {

	vector<double> fReal(Obj, 0.0);
	vector<double> Real = RealObj.WFG_4_random_soln(k, l, data1[i]);

	RealObj.evaluate(Real, fReal);

	vector<double> tmp;
	for(int m = 1; m <=Obj; m++)
	{
	   double scale = m*2;
	   tmp.push_back(  RealObj.getXM() + concave( data2[i], m, Obj)*scale);
	}
	Points.push_back(tmp);
    }
	if(Delta < TOL) return Points;	
    //limpiar 
     vector<vector<double> >out;
     Filter(Points, out);
	

   return out;
}
vector< vector<double> > UniformGenerator::GenerateSolutionWFG5(int N, int Obj ,int k, int l)
{
    vector< vector<double> > data1 = GenerateQuasiRandom(k, N); 
    vector< vector<double> > data2 = GenerateQuasiRandom(Obj, N); 
    vector< vector<double> > Points; 
    
    WFG5 RealObj, LightObj;

    RealObj.init(Obj, k, l); 
 //   LightObj.init(Obj, a, b ); //Se inicializa con el mínimo número de elementos...
     
    for(int i = 0; i < N; i++)
    {

	vector<double> fReal(Obj, 0.0);
	vector<double> Real = RealObj.WFG_5_random_soln(k, l, data1[i]);

	RealObj.evaluate(Real, fReal);

	vector<double> tmp;
	for(int m = 1; m <=Obj; m++)
	{
	   double scale = m*2;
	   tmp.push_back(  RealObj.getXM() + concave( data2[i], m, Obj)*scale);
	}
	Points.push_back(tmp);
    }
	if(Delta < TOL) return Points;	
    //limpiar 
     vector<vector<double> >out;
     Filter(Points, out);
	

   return out;
}
vector< vector<double> > UniformGenerator::GenerateSolutionWFG6(int N, int Obj ,int k, int l)
{
    vector< vector<double> > data1 = GenerateQuasiRandom(k, N); 
    vector< vector<double> > data2 = GenerateQuasiRandom(Obj, N); 
    vector< vector<double> > Points; 
    
    WFG6 RealObj, LightObj;

    RealObj.init(Obj, k, l); 
 //   LightObj.init(Obj, a, b ); //Se inicializa con el mínimo número de elementos...
     
    for(int i = 0; i < N; i++)
    {

	vector<double> fReal(Obj, 0.0);
	vector<double> Real = RealObj.WFG_6_random_soln(k, l, data1[i]);

	RealObj.evaluate(Real, fReal);

	vector<double> tmp;
	for(int m = 1; m <=Obj; m++)
	{
	   double scale = m*2;
	   tmp.push_back(  RealObj.getXM() + concave( data2[i], m, Obj)*scale);
	}
	Points.push_back(tmp);
    }
	if(Delta < TOL) return Points;	
    //limpiar 
     vector<vector<double> >out;
     Filter(Points, out);
	
   return out;
}
vector< vector<double> > UniformGenerator::GenerateSolutionWFG7(int N, int Obj ,int k, int l)
{
   vector< vector<double> > data1 = GenerateQuasiRandom(k, N); 
    vector< vector<double> > data2 = GenerateQuasiRandom(Obj, N); 
    vector< vector<double> > Points; 
    
    WFG7 RealObj, LightObj;

    RealObj.init(Obj, k, l); 
 //   LightObj.init(Obj, a, b ); //Se inicializa con el mínimo número de elementos...
     
    for(int i = 0; i < N; i++)
    {

	vector<double> fReal(Obj, 0.0);
	vector<double> Real = RealObj.WFG_7_random_soln(k, l, data1[i]);

	RealObj.evaluate(Real, fReal);

	vector<double> tmp;
	for(int m = 1; m <=Obj; m++)
	{
	   double scale = m*2;
	   tmp.push_back(  RealObj.getXM() + concave( data2[i], m, Obj)*scale);
	}
	Points.push_back(tmp);
    }
	if(Delta < TOL) return Points;	
    //limpiar 
     vector<vector<double> >out;
     Filter(Points, out);
	

   return out;
}
vector< vector<double> > UniformGenerator::GenerateSolutionWFG8(int N, int Obj ,int k, int l)
{
   vector< vector<double> > data1 = GenerateQuasiRandom(k, N); 
    vector< vector<double> > data2 = GenerateQuasiRandom(Obj, N); 
    vector< vector<double> > Points; 
    
    WFG8 RealObj, LightObj;

    RealObj.init(Obj, k, l); 
 //   LightObj.init(Obj, a, b ); //Se inicializa con el mínimo número de elementos...
     
    for(int i = 0; i < N; i++)
    {

	vector<double> fReal(Obj, 0.0);
	vector<double> Real = RealObj.WFG_8_random_soln(k, l, data1[i]);

	RealObj.evaluate(Real, fReal);

	vector<double> tmp;
	for(int m = 1; m <=Obj; m++)
	{
	   double scale = m*2;
	   tmp.push_back(  RealObj.getXM() + concave( data2[i], m, Obj)*scale);
	}
	Points.push_back(tmp);
    }
	if(Delta < TOL) return Points;	
    //limpiar 
     vector<vector<double> >out;
     Filter(Points, out);
	
   return out;
}
vector< vector<double> > UniformGenerator::GenerateSolutionWFG9(int N, int Obj ,int k, int l)
{
    vector< vector<double> > data1 = GenerateQuasiRandom(k, N); 
    vector< vector<double> > data2 = GenerateQuasiRandom(Obj, N); 
    vector< vector<double> > Points; 
    
    WFG9 RealObj, LightObj;

    RealObj.init(Obj, k, l); 
 //   LightObj.init(Obj, a, b ); //Se inicializa con el mínimo número de elementos...
     
    for(int i = 0; i < N; i++)
    {

	vector<double> fReal(Obj, 0.0);
	vector<double> Real = RealObj.WFG_9_random_soln(k, l, data1[i]);

	RealObj.evaluate(Real, fReal);

	vector<double> tmp;
	for(int m = 1; m <=Obj; m++)
	{
	   double scale = m*2;
	   tmp.push_back(  RealObj.getXM() + concave( data2[i], m, Obj)*scale);
	}
	Points.push_back(tmp);
    }

	if(Delta < TOL) return Points;	
    //limpiar 
     vector<vector<double> >out;
     Filter(Points, out);
	
   return out;
}
/**
 hyper - plane sum f_i = 0.5
**/
vector< vector<double> > UniformGenerator::GenerateSolutionDTLZ1(int N, int Obj)
{
   vector<vector<double> >out, ref = Hyperplane(N, Obj) ;
	if(Delta < TOL) return ref;	
     Filter( ref , out);
   return out;
}

vector< vector<double> > UniformGenerator::GenerateSolutionDTLZ2_to_DTLZ4(int N, int Obj)
{
  vector< vector<double> > x = GenerateQuasiRandom(Obj, N);
   vector<double> y(Obj, 0.0);//, x(M, 0.0);
   vector<vector<double> > Nurb;
   for(int i  = 0; i <  N; i++)
   {
	for(int m = 1; m <=Obj; m++)
	{
	   y[m-1] =  concave(x[i], m, Obj);
	}
	Nurb.push_back(y);
   }
  if(Delta < TOL) return Nurb;	
  vector<vector<double> >out;
  Filter(Nurb, out);
  return out;

}
vector< vector<double> > UniformGenerator::GenerateSolutionDTLZ5_to_DTLZ6(int N, int Obj)
{
  vector< vector<double> > Points; 
  vector< vector<double> > x = GenerateQuasiRandom(1, N);
  double t = M_PI / 4.0;

   for(int i = 0; i < N; i++)
    {
        vector<double> theta (Obj, 0.0);	
        vector<double> y (Obj, 1.0);	
	theta[0] = x[i][0]*(M_PI/2.0);

	for(int m = 1; m < Obj-1; m++)
	   theta[m] = t;
	for (int i = 0; i < Obj; i++)
	{
            for (int j = 0; j < Obj - (i + 1); j++)
	      y[i] *= cos(theta[j]);
	      if (i != 0){
		int aux = Obj - (i + 1);
		y[i] *= sin(theta[aux]);
	      } // if
        }
       Points.push_back(y);
    }
   return Points;
}
vector< vector<double> > UniformGenerator::GenerateSolutionDTLZ7(int N, int Obj)
{
  vector< vector<double> > Points; 
  vector< vector<double> > x = GenerateQuasiRandom(Obj, N);

   for(int i = 0; i < N; i++)
    {
	vector<double> y = x[i];
	double h = 0.0;
	for(int j = 0; j < Obj - 1; j++)
    	   h += (y[j] / 2.0) * (1.0 + sin(3.0 * M_PI * y[j]));
      	h=Obj -h;
	y[Obj-1] = 2.0*h;

	bool flag =false;
//        for(int j = 0; j < Obj; j++) if( y[j] < 0.5 ) flag=true;
//        if( y[Obj-1] > (2.0*Obj)-1.0e-2 ) continue;
//	if(flag==false)
	Points.push_back(y);
    }
   Non_Dominance(Points);

  if(Delta < TOL) return Points;	
 vector<vector<double> >out;
  Filter(Points, out);

   return out;
}

bool UniformGenerator::Dominate(vector<double> &A, vector<double> &B)
{
   for(int i = 0; i < A.size(); i++)
    {
            if( A[i] >= B[i])
            {
                return false;
            }
    }

    return true;
  int dominate1 ; // dominate1 indicates if some objective of solution1
                  // dominates the same objective in solution2. dominate2
  int dominate2 ; // is the complementary of dominate1.

  dominate1 = 0 ;
  dominate2 = 0 ;

  int flag; //stores the result of the comparison

 for (int i = 0; i < A.size(); i++) {
    double value1 = A[i];
    double value2 = B[i];
    if (value1 < value2) {
      flag = -1;
    } else if (value1 > value2) {
      flag = 1;
    } else {
      flag = 0;
    }

    if (flag == -1) {
      dominate1 = 1;
    }

    if (flag == 1) {
      dominate2 = 1;
    }
  }

  if (dominate1 == dominate2) {
    return 0; //No one dominate the other
  }
  if (dominate1 == 1) {
    return -1; // solution1 dominate
  }
  return 1;    // solution2 dominate

} 
void UniformGenerator::Non_Dominance(vector<vector<double> > &S)
{
	vector<vector<double> >P;
	
	vector<vector<double> >::iterator it_s, it_p;

	for( it_s=S.begin(); it_s != S.end() ; it_s++ )
	{
		bool dominated = false;
		it_p = P.begin();
		while(it_p != P.end() && !dominated)
		{
			if( Dominate( *it_p, *it_s ))
				dominated = true;
			else if( Dominate( *it_s ,*it_p) )
			{
				iter_swap( it_p, P.end()-1);
				P.pop_back();
			}
			else it_p++;
		}
		if(!dominated)
			P.push_back(*it_s);
	}
	S=P;
}
