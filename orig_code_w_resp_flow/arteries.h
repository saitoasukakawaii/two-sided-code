/***************************************************************************/
/*                                                                         */
/*  Program: arteries.h                                                    */
/*  Version: 2.0                                                           */
/*  By: Mette Olufsen, Math-Tech                                           */
/*  Date: 14. Jan. 1997                                                    */
/*                                                                         */
/*  This header file includes the definition of the global constants, and  */
/*  the definition of the class representing one blood vessel.             */
/*                                                                         */
/***************************************************************************/

// $Id: arteries.h,v 1.9 2005/07/07 22:19:51 heine Exp $
// Last updated on October 23, 2014 by M. Umar Qureshi

#ifndef _ARTERIES_H
#define _ARTERIES_H

#include <cstdio>
#include <cmath>

// Global parameters imported from main.h

extern double   conv, rho, mu, mu_pl, nu, Lr, Lr2, Lr3, g, q, Fr2,
                Re, p0, pmean, tmst, Period, Fcst, CO, COm,
				Deltat;		//tau; 
                
extern double  *fjac[18], *fj[8];

// The class structure.
class Tube {
public:
  double L;                    // The length of the vessel
  double rtop, rbot;           // The top and bottom radii of the vessel
  Tube *LD, *RD;               // The left and right daughter-vessels. If
                               // these are set the Peripheral resistance
                               // should be set to zero.
  double rm;                             
  double pts;                  // The number of grid points per cm
  int init;
  double K_loss;
  double ff1, ff2, ff3, fa1, fa2, fa3, fv1, fv2, fv3;
  double asym, expo;
  double lrrA, lrrV;

  int N;                       // The number of grid points along the vessel
  double h;                    // The interval length of delta x
  double RLrb;                 // The peripheral resistance of the vessel

  double Ah05, Qh05;
  
  double *Qnew, *Qold, *Qh,    // The arrays needed to store data during
         *Anew, *Aold, *Ah,    // the numerical solution of the system.
         *R1, *R2, *R1h, *R2h,
         *S1, *S2, *S1h, *S2h,
         *Qprv, *Aprv,
         *pL, *y11, *y12, *y21, *y22, *QL, *Pout, *Z,
     *Q0,
     //*Ps,
	 *r0, *r0h,
	 *dr0dx, *dr0dxh,
	 *A0, *A0h, *wom,
	 *fr, *frh,
	 *dfrdr0, *dfrdr0h,
	 *p1, *p1h,
	 *dp1dr0, *dp1dr0h;

  Tube (double Length,
        double topradius, double botradius,
        Tube *LeftDaughter, Tube *RightDaughter,
        double rmin, double points, int init, double K,
        double f1, double f2, double f3,
        double fa1, double fa2, double fa3,
        double fv1, double fv2, double fv3,
        double asym, double expo,
        double lrrA, double lrrV);
                                                         // Constructor.
  ~Tube ();                                              // Destructor.


  // Prints P(x_fixed,t), A(x_fixed,t), F(x_fixed,t), or Q(x_fixed,t) for all
  // t's along the tube.
    
  void printQ0  (FILE *fd);
  //void printPs  (FILE *fd);
  void printPt  (FILE *fd, double t, int i);
  void printQt  (FILE *fd, double t, int i);
  void printAt  (FILE *fd, double t, int i);
  void printPUt (FILE *fd, double t, int i);
  void printFt  (FILE *fd, double t, int i);


  // Prints P(x,t), A(x,t), F(x,t), or Q(x,t) for all x's and t's
  // along the tube. The argument offset makes sure that the vessel
  // is located with the right offset from the inlet.
  void printPxt (FILE *fd, double t, int offset);
  void printAxt (FILE *fd, double t, int offset);
  void printFxt (FILE *fd, double t, int offset);
  void printQxt (FILE *fd, double t, int offset);

  // Prints P as a function of Q and A, respectively,  at a given location x.
  void printPQ (FILE *fd, int i);
  void printPA (FILE *fd, int i);

  // Prints dP/dx(x_fixed,t), dQ/dx(x_fixed,t), dA/dt(x_fixed,t),
  // dQ/dt(x_fixed,t), Fric(x_fixed, t), TotConRes(x_fixed,t),
  // TotMomRes(x_fixed,t) for any fixed x and for all t.
  // I.e. all terms including the sums of the momentum and
  // continuity equations.
  void printdAdt      (FILE *fd, double t, int i, double Aprev, double tmst);
  void printdQdx      (FILE *fd, double t, int i);
  void printTotConRes (FILE *fd, double t, int i, double Qprev, double tmst);

  void printdQdt      (FILE *fd, double t, int i, double Qprev, double tmst);
  void printddxQ2divA (FILE *fd, double t, int i);
  void printdPdx      (FILE *fd, double t, int i);
  void printFric      (FILE *fd, double t, int i);
  void printTotMomRes (FILE *fd, double t, int i, double Qprev, double tmst);


  // Defines P(x,A(x,t)).
  double P (int i, double A);

  // Defines dPdA(x,A(x,t)).
  double dPdA (int i, double A);

  // Defines dPdx1(x,A(x,t)).
  double dPdx1 (int i, double A);

  // Defines B(x,A(x,t)).
  double B (int i, double A);

  // Defines Bh(x,A(x,t)).
  double Bh (int i, double A);

  // Defines dBdx1(x,A(x,t)).
  double dBdx1 (int i, double A);

  // Defines dBdx1h(x,A(x,t)).
  double dBdx1h (int i, double A);

  // Defines dBdAh (x,A(x,t)).
  double dBdAh (int i, double A);

  // Defines d2BdAdxh (x, A(x,t));
  double d2BdAdxh (int i, double A);


  // Tests that the CFL-condition is valid throughout the tube.
  double CFL ();


  // Finds the flux acc. to sys. eq.
  double Rvec (int k, int i, int j, double Q, double A);

  // Finds the rhs. of system eq.
  double Svec (int k, int i, int j, double Q, double A);


  // Steps through interior points.
  void step (double k);


  // Updates left bndry. This should only be done for the inlet tube.
  void bound_left (double t, double k, double Period);

  // Updates the matching boundary.
  void bound_match (int qLnb, double t, double k, double theta, double gamma);


  // Updates right bndry. This should only be done for terminal vessels.
  double c  (int i, double A); // The wave speed through aorta.
  // double Hp (int i, double Q, double A);
  // void poschar (double theta, double &qR, double &aR, double &cR, double &HpR);
  // void bound_right (int qLnb, double k, double theta, double t);


  // Updates bifurcation conditions. Uses daughter vessels, and should
  // only be called when such exist.
  void bound_bif (double theta, double gamma);

  // Updates outflow.
  void bound_right (double t, double k, double Period);

  // In order to ensure a more efficient execution of the program the following
  // functions is made as in-line functions.

// A function returning the Friction of the system. The definition of this
// function is given according to the derivation in the mathematical model.
// The constant cst, determines the amount of damping in the system.
inline double F (double Q, double A)
{
  double tmp1 = Fcst*M_PI*Q;
  double tmp2 = A*Re;
  double tmp3 = -tmp1/tmp2;
  return(tmp3);
  // return(-Fcst*M_PI*Q/(A*Re));
}

inline double dFdQ (double A)
{
  return(-Fcst*M_PI/(A*Re));
}

inline double dFdA (double Q, double A)
{
  return(Fcst*M_PI*Q/(sq(A)*Re));
}

private:
  // The private function Q0 may only be accessed from the left boundary
  // function. It ensures a certain and given CO (defined in main.h).
  double Q0_init (double t, double k, double Period);
  //double Ps_init (double t, double k, double Period);
  double PL_init (double t, double k, double Period);
};

void solver (Tube *Arteries[], double tstart, double tend, double k);

#endif
