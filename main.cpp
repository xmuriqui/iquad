
#include <iostream>
#include <new>

#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>


#include "iquad.hpp"
#include "IQD_qpGen.hpp"
#include "IQD_tools.hpp"


using namespace iquad;
using namespace std;



int main(int argc, char **argv)
{
    int i = 1;
    string option;
	char *outDirectory;
	char curDir[] = "./";
	IQD_IQuadProb prob;
	IQD_BranchAndBound bb;
	IQD_GeneralSolverParams sParams, mosekParams;
	
	
	IQD_welcome();
	
	
	if( argc > 3 )
		outDirectory = argv[3];
	else
		outDirectory = curDir;
	
	
	
	
    if( argc > 2 )
    {
		printf("parameters: %s %s\n", argv[1], argv[2]);
		//getchar();
		
		for(i = 0; argv[2][i] != '\0'; i++ )
		    argv[2][i] = tolower(argv[2][i]);
		    
		option = argv[2];
		
		if(option == "-biqmac")
		    IQD_biqmac(argv[1], outDirectory, true);
		else if(option == "-rbiqmac")
			IQD_biqmac(argv[1], outDirectory, false);
		else if(option == "-rmaxcut")
			IQD_biqmac(argv[1], outDirectory, false, true);
		else if(option == "-boxqp") 
		    IQD_boxqp(argv[1], outDirectory);
		else if(option == "-myqp") 
			IQD_myqp(argv[1], outDirectory);
		else if(option == "-qpprog")
			IQD_readQuadProg(argv[1], outDirectory);
		else
		    //ampl format
		    IQD_ampl(argv[1], argv, outDirectory);
    }
    else if(argc == 2)
    {
		IQD_ampl(argv[1], argv, outDirectory);
    }
    else
	{
		std::cout << "usage: " << argv[0] << ": <input problem file> [problem format]\n"
		"\twhere format can be:\n"
		"\t\t-biqmac\n"
		"\t\t-rbiqmact"
		"\t\t-rmaxcut\n"
		"\t\t-boxqp\n"
		"\t\t-myqp\n"
		"\t\t-qpprog\n"
		"\t\t-AMPL (default)\n"
		;
		
		
	}
	#if 0
    else
	{
		
		int i;
		char name[50];
		IQD_IQuadProb *prob;
		IQD_QPGenerator qpGen;
		
		//prob2 = new IQD_IQuadProb;
		
		qpGen.denv = 0.5;
		i = 1986;
		qpGen.setSeed(&i);
		qpGen.L = 20;
		qpGen.factor = 1.0;
		qpGen.lowerBoundToGenerateD = 1;
		qpGen.condd = 1;
		
		/*
		
		//prob = new IQD_IQuadProb;
		//qpGen.generateProblem(0, 2, 1, *prob);
		//qpGen.generateProblem(1, 2, 2, *prob);
		//prob->writeMIQCPModelInAMPLFile("iquad.mod");
		//prob->printProblem();
		//getchar();
		
		
		
		for(i = 20; i <= 100; i = i + 20)
		{
			for(j = 1; j <= 4; j++)
			{
				prob = new IQD_IQuadProb;
				qpGen.generateProblem(i/4, 0, i/4, *prob);
				sprintf( name, "myqp1_%03d_g1p%d_%02d_%02d_%02d.in", i, j, i/4, 0, i/4 );
				prob->writeMIQCPModelInFile( name );
				delete prob;
				
				
				prob = new IQD_IQuadProb;
				qpGen.generateProblem(0, i/4, i/4, *prob);
				sprintf( name, "myqp1_%03d_g2p%d_%02d_%02d_%02d.in", i, j, 0, i/4, i/4 );
				prob->writeMIQCPModelInFile( name );
				delete prob;
				
				
				prob = new IQD_IQuadProb;
				qpGen.generateProblem(i/6, i/6,  i/2 - 2*(i/6) , *prob);
				sprintf( name, "myqp1_%03d_g3p%d_%02d_%02d_%02d.in", i, j, i/6, i/6, i/2 - 2*(i/6) );
				prob->writeMIQCPModelInFile( name );
				delete prob;
			}
		}
		*/
		
		//generating a special problem having 50 variables
		prob = new IQD_IQuadProb;
		qpGen.generateProblem(8, 8, 9, *prob);
		sprintf(name, "myqp2_050_g3p1_8_8_9.in");
		prob->writeMIQCPModelInFile(name);
		
		delete prob;
		
		
		//generating a special problem having 50 variables
		prob = new IQD_IQuadProb;
		qpGen.generateProblem(8, 8, 9, *prob);
		sprintf(name, "myqp2_050_g3p2_8_8_9.in");
		prob->writeMIQCPModelInFile(name);
		
		delete prob;
		
		
		/*
			prob.printProblem();
			prob.writeMIQCPModelInFile("iquadProb.txt");
			getchar();
		*/
		
		
		/*
		//sParams.addIntegerParameter("CPX_PARAM_THREADS", 1 );
	    //sParams.addDoubleParameter("CPX_PARAM_TILIM", 4*60*60);

	    mosekParams.addIntegerParameter("MSK_IPAR_NUM_THREADS", 1 );
	    mosekParams.addDoubleParameter("MSK_DPAR_OPTIMIZER_MAX_TIME", 4*60*60);
	    
		bb.in_sdp_solver = iquad::IQD_SS_CSDP;
		
		bb.in_printing_frequency = 100;
		bb.in_split_way = iquad::IQD_SQ_SCHUR;
		bb.in_split_way = iquad::IQD_SQ_DIAG_SDP;
		bb.in_qcp_solver = iquad::IQD_QS_GUROBI;
		bb.in_eps_added_sd_solver_matrices = 1.0e-4;
		bb.run(*prob, &sParams, &mosekParams);
		
		printf("____________________________________________________________________________\n");
	    printf("Return code: %d obj function: %f time: %0.2f cpu time: %0.2f splitting time: %0.2f\n", bb.out_return_code, bb.out_best_obj, bb.out_clock_time, bb.out_cpu_time, bb.out_splitting_cpu_time);
	    //for(i = 0; i < n; i++)
		//printf("x[%d]: %f\n", i, bb.out_best_sol[i]);
		printf("____________________________________________________________________________\n");
		
		*/
		
		
		/*
			int i = prob.readMIQCPModelInFile("iquadProb.txt");
			printf("Problema lido: %d\n", i);
			
			prob.printProblem();
		*/
		
	}
	
	
	#endif
	
	return 0;
}








