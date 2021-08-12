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
   message += "Example:\n  ./GenratorU --n 100 --m 3 > 100.dat";
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

  int Objectives =2, Vectors = 100, k =2, l=2, nvar=6;
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
		else if(Terminal == "--nvar" )
			nvar = atoi(argv[++i]);
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
     data = Obj.GenerateSolutionWFG4toWFG9(Vectors, Objectives, k, l);
  else if( Problem == "WFG5")
     data = Obj.GenerateSolutionWFG4toWFG9(Vectors, Objectives, k, l);
  else if( Problem == "WFG6")
     data = Obj.GenerateSolutionWFG4toWFG9(Vectors, Objectives, k, l);
  else if( Problem == "WFG7")
     data = Obj.GenerateSolutionWFG4toWFG9(Vectors, Objectives, k, l);
  else if( Problem == "WFG8")
     data = Obj.GenerateSolutionWFG4toWFG9(Vectors, Objectives, k, l);
  else if( Problem == "WFG9")
     data = Obj.GenerateSolutionWFG4toWFG9(Vectors, Objectives, k, l);
  else if( Problem == "minusWFG1")
     data = Obj.GenerateSolutionminusWFG1(Vectors, Objectives, k, l);
  else if( Problem == "minusWFG2")
     data = Obj.GenerateSolutionminusWFG2(Vectors, Objectives, k, l);
  else if( Problem == "minusWFG3")
     data = Obj.GenerateSolutionminusWFG3(Vectors, Objectives, k, l);
  else if( Problem == "minusWFG4")
     data = Obj.GenerateSolutionminusWFG4toWFG9(Vectors, Objectives, k, l);
  else if( Problem == "minusWFG5")
     data = Obj.GenerateSolutionminusWFG4toWFG9(Vectors, Objectives, k, l);
  else if( Problem == "minusWFG6")
     data = Obj.GenerateSolutionminusWFG4toWFG9(Vectors, Objectives, k, l);
  else if( Problem == "minusWFG7")
     data = Obj.GenerateSolutionminusWFG4toWFG9(Vectors, Objectives, k, l);
  else if( Problem == "minusWFG8")
     data = Obj.GenerateSolutionminusWFG4toWFG9(Vectors, Objectives, k, l);
  else if( Problem == "minusWFG9")
     data = Obj.GenerateSolutionminusWFG4toWFG9(Vectors, Objectives, k, l);
  else if( Problem == "DTLZ1")
     data = Obj.GenerateSolutionDTLZ1(Vectors, Objectives);
  else if(Problem =="IMB2" || Problem =="IMB8" || Problem =="IMB4" || Problem == "IMB10" || Problem=="IMB6")
    {
     data = Obj.GenerateSolutionDTLZ1(Vectors, Objectives);
	for(auto &point:data)
	{
	  for(auto &var:point)
	  {
		var *=2.0;
	  }
	}
    } 
  else if( Problem == "DTLZ2" || Problem == "DTLZ3" || Problem == "DTLZ4"||Problem=="IMB5" || Problem =="IMB3"||Problem=="IMB9")
     data = Obj.GenerateSolutionDTLZ2_to_DTLZ4(Vectors, Objectives);
  else if( Problem == "DTLZ5" || Problem == "DTLZ6" )
     data = Obj.GenerateSolutionDTLZ5_to_DTLZ6(Vectors, Objectives);
  else if( Problem == "DTLZ7")
     data = Obj.GenerateSolutionDTLZ7(Vectors, Objectives);
  else if( Problem == "BT1" || Problem == "BT2" || Problem == "BT3" || Problem == "BT4" || Problem == "BT6"|| Problem == "BT7" || Problem == "BT8" || Problem == "IMB1" || Problem == "IMB7")
     data = Obj.GenerateSolutionBT1toBT4andBT6toBT8(Vectors);
  else if( Problem == "BT5")
     data = Obj.GenerateSolutionBT5(Vectors);
  else if( Problem == "BT9")
     data = Obj.GenerateSolutionBT9(Vectors);
  else if( Problem == "minusDTLZ1")
     data = Obj.GenerateSolutionminusDTLZ1(Vectors, Objectives, nvar);
  else if( Problem == "minusDTLZ2")
     data = Obj.GenerateSolutionminusDTLZ2(Vectors, Objectives, nvar);
  else if( Problem == "minusDTLZ3")
     data = Obj.GenerateSolutionminusDTLZ3(Vectors, Objectives, nvar);
  else if( Problem == "minusDTLZ4")
     data = Obj.GenerateSolutionminusDTLZ4(Vectors, Objectives, nvar);





  for(int i = 0; i < data.size(); i++)
  {
	for(int j = 0; j< Objectives; j++)
	cout << data[i][j] << " ";
        cout << endl;
  }


   return 0;
}
