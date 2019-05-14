/**
   Autor: Joel Chacón Castillo
   Date: 07/may/2017
   
   This is a generator of weight vectors for MOEA/D based on the paper:
  "MOEA/D + uniform design: A new version of MOEA/D for optimization problems with many objectives"
   by
   Yan-yan Tan, Yong-chang jiao, Hong Li, Xin-kuan Wang.

**/
#include <iostream>
#include "UniformGenerator.hpp"
using namespace std;
 string verbose()
{
   string message = "";
   message += "--h --help print this summary and exit\n";
   message += "--n Number of points to generate\n";
   message += "--k Position params\n";
   message += "--l Distance params\n";
   message += "--p Problem\n";
   message += "--m Number of objectives\n";
   message += "Example:\n  ./GenratorU --n 100 --o 3 > 100.dat";
   return message;
}
int main(int argc, char * argv[])
{
   if(argc<2)
         {
	    
	    cout << "Unknown Argument.."<<endl;
	    cout << verbose();
	    exit(0);
	 }

  int Objectives =2, Vectors = 100, k =2, l=2;
   double delta= 1e-8;
   string Problem = "WFG1";
   for(int i = 1; i < argc ; i++)
    	{
		string Terminal(argv[i]);
		if( Terminal == "--h")
			cout << verbose<<endl;
		else if(Terminal == "--m" )
			Objectives = atoi(argv[++i]);
		else if(Terminal == "--n" )
			Vectors = atoi(argv[++i]);
		else if(Terminal == "--k" )
			k = atoi(argv[++i]);
		else if(Terminal == "--l" )
			l = atoi(argv[++i]);
		else if(Terminal == "--p")
			Problem = string(argv[++i]);
		else if(Terminal == "--d")
			delta = atof(argv[++i]);
		else
		{
		   cout << "Unknown Argument...";
		   exit(0);
		}
	}





  UniformGenerator Obj;
   Obj.setDelta(delta);
   vector< vector<double> >data;
  //vector< vector<double> > Wvec = Obj.Hyperplane(Vectors, Objectives);
  if( Problem == "WFG1")
     data = Obj.GenerateSolutionWFG1(Vectors, Objectives, k, l);
  else if( Problem == "WFG2")
     data = Obj.GenerateSolutionWFG2(Vectors, Objectives, k, l);
  else if( Problem == "WFG3")
     data = Obj.GenerateSolutionWFG3(Vectors, Objectives, k, l);
  else if( Problem == "WFG4")
     data = Obj.GenerateSolutionWFG4(Vectors, Objectives, k, l);
  else if( Problem == "WFG5")
     data = Obj.GenerateSolutionWFG5(Vectors, Objectives, k, l);
  else if( Problem == "WFG6")
     data = Obj.GenerateSolutionWFG6(Vectors, Objectives, k, l);
  else if( Problem == "WFG7")
     data = Obj.GenerateSolutionWFG7(Vectors, Objectives, k, l);
  else if( Problem == "WFG8")
     data = Obj.GenerateSolutionWFG8(Vectors, Objectives, k, l);
  else if( Problem == "WFG9")
     data = Obj.GenerateSolutionWFG9(Vectors, Objectives, k, l);
  else if( Problem == "DTLZ1")
     data = Obj.GenerateSolutionDTLZ1(Vectors, Objectives);
  else if( Problem == "DTLZ2" || Problem == "DTLZ3" || Problem == "DTLZ4" )
     data = Obj.GenerateSolutionDTLZ2_to_DTLZ4(Vectors, Objectives);
  else if( Problem == "DTLZ5" || Problem == "DTLZ6" )
     data = Obj.GenerateSolutionDTLZ5_to_DTLZ6(Vectors, Objectives);
  else if( Problem == "DTLZ7")
     data = Obj.GenerateSolutionDTLZ7(Vectors, Objectives);

  for(int i = 0; i < data.size(); i++)
  {
	for(int j = 0; j< Objectives; j++)
	cout << data[i][j] << " ";
        cout << endl;
  }


   return 0;
}
