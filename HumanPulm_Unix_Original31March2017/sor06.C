/* The sor06.C main program */

// $Id: sor06.C,v 1.10 2005/10/14 18:05:59 heine Exp $
//Last updated on October 23, 2014 by M. Umar Qureshi

#include "sor06.h"
#include "tools.h"
#include "arteries.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstring>

using namespace std;

// The vessel-network is specified and initialized, and the flow and
// pressures are to be determed. All global constants must be defined
// in the header file. That is the number of dots per cm, and the number
// of vessels.
int main(int argc, char *argv[])
{
  double tstart, tend, finaltime;
	
  double rm; double f1; double f2; double f3; 
  double fa1; double fa2; double fa3; 
  double asym; double expo; double lrrA; double lrrV; 
  bool verbosity; double numHeartBeats; int id;//double *Q0;
	
	// check command line args
  if (argc != 15) //argv[0] is the name of the program, here sor06
    {
      printf("Not enough input arguments, noargc %d and they are %s\n", argc, argv[0]);
      return 1;
    }
	
	rm = atof(argv[1]); 
	f1 = atof(argv[2]); f2 = atof(argv[3]); f3 = atof(argv[4]);
	fa1 = atof(argv[5]); fa2 = atof(argv[6]); fa3 = atof(argv[7]);
	asym = atof(argv[8]); expo = atof(argv[9]); lrrA = atof(argv[10]);
	lrrV = atof(argv[11]);
	verbosity = atoi(argv[12]); numHeartBeats = atof(argv[13]); 
	id = atoi(argv[14]);
	
	double fv1 = fa1; double fv2 = fa2; double fv3 = fa3;  

	char namepu1 [20]; char namepu2 [20]; char namepu3 [20];
	char namepu4 [20]; char namepu5 [20]; char namepu6 [20];
	char namepu7 [20]; char namepuv1 [20]; char namepuv2 [20];
	char namepuv3 [20]; char namepuv4 [20]; 
	
    sprintf(namepu1, "pu1_%d.2d", id); 
    sprintf(namepu2, "pu2_%d.2d", id);
    sprintf(namepu3, "pu3_%d.2d", id);
    sprintf(namepu4, "pu4_%d.2d", id);
    sprintf(namepu5, "pu5_%d.2d", id);
    sprintf(namepu6, "pu6_%d.2d", id);
    sprintf(namepu7, "pu7_%d.2d", id);
    sprintf(namepuv1, "puv1_%d.2d", id);
    sprintf(namepuv2, "puv2_%d.2d", id);
    sprintf(namepuv3, "puv3_%d.2d", id);
    sprintf(namepuv4, "puv4_%d.2d", id);
  
    if (verbosity)
    {
    fprintf(stdout, "Opening Files:\n");
    }
	
	FILE *fpu1 = fopen (namepu1, "w");
	if (verbosity)
	{
	if (fpu1) fprintf(stdout, "1pu OK, "); else error ("main.C","File 1pu NOT OK");
	}
	
	FILE *fpu2 = fopen (namepu2, "w");
	if (verbosity)
	{
	if (fpu2) fprintf(stdout, "2pu OK, "); else error ("main.C","File 2pu NOT OK");
	}
	
	FILE *fpu3 = fopen (namepu3, "w");
	if (verbosity)
	{
	if (fpu3) fprintf(stdout, "3pu OK, "); else error ("main.C","File 3pu NOT OK");
	}
	
	FILE *fpu4 = fopen (namepu4, "w");
	if (verbosity)
	{
	if (fpu4) fprintf(stdout, "4pu OK, "); else error ("main.C","File 4pu NOT OK");
	}
	
	FILE *fpu5 = fopen (namepu5, "w");
	if (verbosity)
	{
	if (fpu5) fprintf(stdout, "5pu OK, "); else error ("main.C","File 5pu NOT OK");
	}
	
	FILE *fpu6 = fopen (namepu6, "w");
	if (verbosity)
	{
	if (fpu6) fprintf(stdout, "6pu OK, "); else error ("main.C","File 6pu NOT OK");
	}
	
	FILE *fpu7 = fopen (namepu7, "w");
	if (verbosity)
	{
	if (fpu7) fprintf(stdout, "7pu OK, "); else error ("main.C","File 7pu NOT OK");
	}
	
    if (verbosity)
    {
    fprintf(stdout, "\n");
    }
  
	FILE *fpuv1 = fopen (namepuv1, "w");
	if (verbosity)
	{
	if (fpuv1) fprintf(stdout, "1puv OK, "); else error ("main.C","File 1pu NOT OK");
	}
	
	FILE *fpuv2 = fopen (namepuv2, "w");
	if (verbosity)
	{
	if (fpuv2) fprintf(stdout, "2puv OK, "); else error ("main.C","File 2puv NOT OK");
	}
	
	FILE *fpuv3 = fopen (namepuv3, "w");
	if (verbosity)
	{
	if (fpuv3) fprintf(stdout, "3puv OK, "); else error ("main.C","File 3puv NOT OK");
	}
	
	FILE *fpuv4 = fopen (namepuv4, "w");
	if (verbosity)
	{
	if (fpuv4) fprintf(stdout, "4puv OK, "); else error ("main.C","File 4puv NOT OK");
	}
	
	if (verbosity)
	{
  fprintf(stdout, ".\n");
	}
	

  // Workspace used by bound_match
  for(int i=0; i<8; i++) fj[i] = new double[8];

  // Workspace used by bound_bif
  for(int i=0; i<18; i++) fjac[i] = new double[18];

  clock_t c1 = clock();       // Only used when timing the program.
  nbrves    = 11;             // Total number of large arteries and veins in the network -- (temp)
  tstart    = 0.0;            // Starting time.
  finaltime = numHeartBeats*Period;      // Final end-time during a simulation.
  tend      = (numHeartBeats-1)*Period;      // Timestep before the first plot-point
                              // is reached.
	
  // The number of vessels in the network is given when the governing array of
  // vessels is declared.

  Tube   *Arteries[nbrves];                    // Array of blood vessels.

  // *---------Initialisation of the vessel network.---------*
  // init == 1 implies initial (inflow) artery.
  // init == 2 implies initial (outflow) vein.
  // init == 3 implies terminal artery to be matched to vein,
  //  (with LD == 0, RD == pointer to vein).
  // r_min == 0 implies vessel has two daughter vessels.
  // r_min != 0 implies 'terminal' vessel,
  //  (terminal arteries - one daughter, 
  //   terminal veins - no daughters).
	
// Parameters required to initiate class Tube (Length,topradius,botradius,LeftDaughter,RightDaughter,rmin, points,
                                               //init,K,f1,f2,f3,fa1,fa2,fa3,fv1,fv2,fv3,asym,expo,lrrA,lrrV);
	

  // Initialization of the Veins. Approximated dimensions: Data 1 in PV_dimensionV2.xlsx
	
  Arteries[10] = new Tube( 2.00, 0.64,0.58, 0, 0, rm, 4,2,0,f1,f2,f3,fa1,fa2,fa3,fv1,fv2,fv3, asym, expo, lrrA, lrrV);
  // LSPV, connected to LTA, approximated radii and length 2 cm
  Arteries[9]  = new Tube( 2.18, 0.97, 0.90, 0, 0, rm, 4,2,0,f1,f2,f3,fa1,fa2,fa3,fv1,fv2,fv3, asym, expo, lrrA, lrrV);
  // LIPV, connected to LIA, approximated radii and length, 2.18 cm
  Arteries[8]  = new Tube( 1.5, 0.51,  0.46, 0, 0, rm, 4,2,0,f1,f2,f3,fa1,fa2,fa3,fv1,fv2,fv3, asym, expo, lrrA, lrrV);
  // RSPV, connected to RTA, approximated radii and length, 1,5 cm,
  Arteries[7]  = new Tube( 1.24, 0.59, 0.55, 0, 0, rm, 4,2,0,f1,f2,f3,fa1,fa2,fa3,fv1,fv2,fv3, asym, expo, lrrA, lrrV);
  // RIPV, connected to RIA, approximated distal radii and length 1.24 cm,
 
  // Initialization of the Arteries.
  Arteries[6] = new Tube( 1.00, 0.58,  0.58, 0, Arteries[10], rm, 4,3,0,f1,f2,f3,fa1,fa2,fa3,fv1,fv2,fv3, asym, expo, lrrA, lrrV);
  // LTA, approximated length = 1 cm,
  Arteries[5] = new Tube( 2.25, 1.04,  0.90, 0, Arteries[9],  rm, 4,3,0,f1,f2,f3,fa1,fa2,fa3,fv1,fv2,fv3, asym, expo, lrrA, lrrV);
  // LIA, measured length = 2.18 cm
  Arteries[4] = new Tube( 1.00, 0.46,  0.46, 0, Arteries[8],  rm, 4,3,0,f1,f2,f3,fa1,fa2,fa3,fv1,fv2,fv3, asym, expo, lrrA, lrrV);
  //RTA, approximated lenghth = 1 cm,
  Arteries[3] = new Tube( 1.25, 0.57,  0.55, 0, Arteries[7],  rm, 4,3,0,f1,f2,f3,fa1,fa2,fa3,fv1,fv2,fv3, asym, expo, lrrA, lrrV);
  //RIA, measured length = 1.24 cm,
  Arteries[2] = new Tube( 2.50, 1.10,  1.08, Arteries[ 5], Arteries[ 6], 0, 4,0,0,f1,f2,f3,fa1,fa2,fa3,fv1,fv2,fv3, asym, expo, lrrA, lrrV);
  //LPA, measured length = 2.48cm,
  Arteries[1] = new Tube( 5.75, 0.93,  0.60, Arteries[ 3], Arteries[ 4], 0, 4,0,0,f1,f2,f3,fa1,fa2,fa3,fv1,fv2,fv3, asym, expo, lrrA, lrrV);
  //RPA, measured length = 5.79cm,
  Arteries[0] = new Tube( 4.50, 1.36,  1.30, Arteries[ 1], Arteries[ 2], 0, 4,1,0,f1,f2,f3,fa1,fa2,fa3,fv1,fv2,fv3, asym, expo, lrrA, lrrV);
  //MPA, measured length = 4.46cm,



  // In the next three statements the simulations are performed until
  // tstart = tend. That is this is without making any output during this
  // first  of time. If one needs output during this period, these three
  // lines should be commented out, and the entire simulation done within the
  // forthcomming while-loop.

  // Solves the equations until time equals tend.
    solver (Arteries, tstart, tend, k);
    tstart = tend;
    tend = tend + Deltat;
    
    if (verbosity)
    {
    fprintf (stdout,"plots start\n");
    }

  // The loop is continued until the final time
  // is reached. If one wants to make a plot of
  // the solution versus x, tend is set to final-
  // time in the above declaration
     while (tend <= finaltime)
  {
    for (int j=0; j<nbrves; j++)
    {
      int ArtjN = Arteries[j]->N;
      for (int i=0; i<ArtjN; i++)
      {
        Arteries[j]->Qprv[i+1] = Arteries[j]->Qnew[i+1];
        Arteries[j]->Aprv[i+1] = Arteries[j]->Anew[i+1];
      }
    }

    // Solves the equations until time equals tend.
    solver (Arteries, tstart, tend, k);
	
    if (verbosity)
	{
    fprintf (stdout,".");
	}

    // A 2D plot of P(x_fixed,t) is made. The resulting 2D graph is due to
    // the while loop, since for each call of printPxt only one point is set.
      
	Arteries[ 0] -> printPxt (fpu1, tend, 0);
	Arteries[ 1] -> printPxt (fpu2, tend, 0);
	Arteries[ 2] -> printPxt (fpu3, tend, 0);
	Arteries[ 3] -> printPxt (fpu4, tend, 0);
    Arteries[ 4] -> printPxt (fpu5, tend, 0);
	Arteries[ 5] -> printPxt (fpu6, tend, 0);
	Arteries[ 6] -> printPxt (fpu7, tend, 0);
	Arteries[ 7] -> printPxt (fpuv1, tend, 0);
	Arteries[ 8] -> printPxt (fpuv2, tend, 0);
    Arteries[ 9] -> printPxt (fpuv3, tend, 0);
	Arteries[10] -> printPxt (fpuv4, tend, 0);
    
	  
	  // The time within each print is increased.
    tstart = tend;
    tend   = tend + Deltat; // The current ending time is increased by Deltat.
  }
  
    if (verbosity)
    {
    fprintf(stdout,"\n");
    }

  // The following statements is used when timing the simulation.
  fprintf(stdout,"nbrves = %d, Lax, ", nbrves);
  clock_t c2 = clock(); // FIXME clock() may wrap after about 72 min.
  int tsec = (int) ((double) (c2-c1)/CLOCKS_PER_SEC + 0.5);
  fprintf(stdout,"cpu-time %d:%02d\n", tsec / 60, tsec % 60);
  
    if (verbosity)
    {
    fprintf(stdout,"\n");
    }

  // In order to termate the program correctly the vessel network and hence
  // all the vessels and their workspace are deleted.
  for (int i=0; i<nbrves; i++) delete Arteries[i];

  // Matrices and arrays are deleted
  for (int i=0; i<18; i++) delete[] fjac[i];
  for (int i=0; i<8;  i++) delete[] fj[i];
  
    if (verbosity)
    {
    fprintf(stdout, "Closing Files: \n");
  	if (fclose (fpu1)  != EOF) fprintf(stdout,"1pu OK, ");
    else error("main.C","Close 1pu NOT OK");
	if (fclose (fpu2)  != EOF) fprintf(stdout,"2pu OK, ");
    else error("main.C","Close 2pu NOT OK");
	if (fclose (fpu3)  != EOF) fprintf(stdout,"3pu OK, ");
    else error("main.C","Close 3pu NOT OK");
	if (fclose (fpu4)  != EOF) fprintf(stdout,"4pu OK, ");
    else error("main.C","Close 4pu NOT OK");
	if (fclose (fpu5)  != EOF) fprintf(stdout,"5pu OK, ");
    else error("main.C","Close 5pu NOT OK");
	if (fclose (fpu6)  != EOF) fprintf(stdout,"6pu OK, ");
    else error("main.C","Close 6pu NOT OK");
	if (fclose (fpu7)  != EOF) fprintf(stdout,"7pu OK, ");
    else error("main.C","Close 7pu NOT OK");
  
  fprintf(stdout, "\n");
	
	if (fclose (fpuv1)  != EOF) fprintf(stdout,"1puv OK, ");
    else error("main.C","Close 1puv NOT OK");
	if (fclose (fpuv2)  != EOF) fprintf(stdout,"2puv OK, ");
    else error("main.C","Close 2puv NOT OK");
	if (fclose (fpuv3)  != EOF) fprintf(stdout,"3puv OK, ");
    else error("main.C","Close 3puv NOT OK");
	if (fclose (fpuv4)  != EOF) fprintf(stdout,"4puv OK, ");
    else error("main.C","Close 4puv NOT OK");
	  
  fprintf(stdout, ".\n");
	}

 return 0;
}
