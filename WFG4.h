/***********************************************************************************
 * AUTORES
 *   - Eduardo Manuel Segredo Gonzalez
 * 
 * FECHA
 *    Julio de 2008
 *
 * DESCRIPCION
 *    Individuo para la resolucion del problema de Benchmark WFG4
 *    
 *    El problema WFG4 esta definido en:
 *     "A Scalable Multi-objective Test Problem Toolkit"
 *      Simon Huband, Luigi Barone, Lyndon While and Phil Hingston
 * 
 *    Se han definido los metodos requeridos por la clase Individual: 
 *     - mutation
 *     - crossover 
 *     - evaluate
 *     - init 
 *     - clone
 *
 *    Para mas informacion acerca de la especificacion de un problema, referirse a
 *    al fichero Individual.h
 *
 * ********************************************************************************/


#ifndef __INDIVIDUALWFG4_H__
#define __INDIVIDUALWFG4_H__


#include <iostream>
#include <math.h>
using namespace std;
//#include "Individual.h"


class WFG4 {
 
  public:
 	// Funcion de inicializacion
	bool init (int nObj, int k, int l);
	
	// Evaluacion y seleccion del individuo
	void evaluate (vector<double> &X, vector<double> &Obj);
	// Rangos
	double inline getMaximum(const int i) const { return 1; }
	double inline getMinimum(const int i) const { return 0; }
	inline double getXM(){return this->Xh;}
	vector< double > WFG_4_random_soln( const int k, const int l, vector<double> &spacek );
//	unsigned int inline getOptDirection(const int i) const { return MINIMIZE; }
	inline int getDimension(){return this->nVar;}
 
 
  private:
	// Funciones auxiliares al problema
	vector<double> subvector(const vector<double> &v, const int head, const int tail);
	double correct_to_01(const double &a, const double &epsilon = 1.0e-10);
	vector<double> calculate_x (const vector<double> &t_p, const vector<short> &A);
  vector<double> calculate_f (const double &D, const vector<double> &x, const vector<double> &h,
								            	const vector<double> &S);

	// Transiciones genericas
	double s_multi (const double &y, const int A, const double &B, const double &C);
	double r_sum (const vector<double> &y, const vector<double> &w);
	
  // Superficies genericas
	double concave(const vector<double> &x, const int m);

	// Transiciones del problema WFG4
	vector<double> WFG4_t1(const vector<double> &y);
	vector<double> WFG4_t2(const vector<double> &y, const int k, const int M);
	
	// Superficie del problema WFG4
	vector<double> WFG4_shape(const vector<double> &t_p);

	// Funciones propias de los problemas WFG
	vector<short> WFG_create_A (const int M, const bool degenerate);
	vector<double> WFG_calculate_f (const vector< double >& x, const vector< double >& h);
	vector< double > WFG_normalise_z( const vector< double >& z );

	// Constantes del problema WFG4
	int k;             // Parametros relacionados con la posicion
	int l;  // Parametros relacionados con la distancia

	int nObj, nVar;
	double Xh;
};

#endif
