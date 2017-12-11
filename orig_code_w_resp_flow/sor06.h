/***************************************************************************/
/*                                                                         */
/*  Program: sor06.h                                                       */
/*  Version: 2.0                                                           */
/*  By: Mette Olufsen, Math-Tech                                           */
/*  Date: 14. Jan. 1997                                                    */
/*                                                                         */
/*  This header file defines the global parameters                         */
/*                                                                         */
/***************************************************************************/

// $Id: sor06.h,v 1.7 2005/07/07 22:19:51 heine Exp $
// Last updated on October 23, 2014 by M. Umar Qureshi

#ifndef _SOR06_H
#define _SOR06_H

#include <cmath>
#include "tools.h"

int    nbrves, N_aorta;               // Number of vessels in the tree.

int    tmstps = 32768,                 // The number of timesteps per period.
	   plts   = 1024;                 // Number of plots per period.

const char*  CO_filename = "gvPA_3_16384.dat";       // Input flow file at the heart.
//const char*  CO_filename = "Flow0.dat";         // Input flow file at the heart.
//const char*  CO_filename = "p0_in4096.dat";     //Changeed to "p0_in4096.dat" if using pressure as an input condition.
const char*  PL_filename = "zero_3.dat";            // Outflow (Static pressure) file.

double  conv   = 1333.220,             // Conversion from mmHg to SI-units.
        rho    = 1.055,                // Density of blood [g/cm^3].
        mu     = 0.049,                // Viscosity of blood [g/cm/s].
        mu_pl  = mu,                   // Viscosity of blood [g/cm/s].
        nu     = mu/rho,

        Tper   = 0.70*8,                 // The period of one heart beat [s].
       
 //     rm     = 0.005,                // The minimum radius for the structured tree [cm].

        
 // Stiffness parameters for the large vessels.
 //      f1     = 0, //1.99925e07,
 //      f2     = 1,// -25.5267,
 //      f3     = 260000,//8.65e05   //which is 0.3*k3 (k3 = 8.65e05 in Olufsen's et al (1999))


 // Stiffness parameters for the small arteries.
 //      fa1    =  0.0, // 1.99925e07,
 //      fa2    =  1.0, // -25.5267,
 //      fa3    =  50000,// 465251,


 // Stiffness parameters for the small veins.
 //      fv1    =  0.0, // 1.99925e07,
 //      fv2    =  1.0, // -25.5267,
 //      fv3    =  50000,//165000,//50000.0, // 465251,


 // To simulate the PAH and CTEPH and HLD in Qureshi et al (2014)
 // Used f3 = 0.48*k3 for HLD plots or 340000 if starting from xi=2.4 (and for PVH)
 // For PAH use 0.45*k3 or 305000 for up to 50% increase in fa3,fv3 and 336000 for up to 75% increase in fa3 and fv3


 //    asym   =  0.41,              // The asymmetry ratio of the structured tree.
 //    expo   =	 2.76,              // Exponent in radius relation. (Decrease up to 2.3 to simulate PH associated with HLD)


 //    lrrA   =  15.75,       // Length to radius ratio in small arteries.
 //    lrrV   =  14.54,       // Length to radius ratio in small veins.


 //    Fcst = 10.0,                 // Determines the damping coeff.
       Fcst   = 50, //17.7778,      // Determines the damping coeff for the friction.
       Lr     = 1.0,                // characteristic radius of the vessels in the tree [cm].
       Lr2    = sq(Lr),             // The squared radius [cm2].
       Lr3    = cu(Lr),             // The radius to the third power [cm^3].
       g      = 981.0,              // The gravitational force [cm/s^2].
       q      = 10.0*Lr2,           // The characteristic flow [cm^3/s].
       Fr2    = sq(q)/g/pow(Lr,5),  // The squared Froudes number.
       Re     = q*rho/mu/Lr,        // Reynolds number.
       Period = Tper*q/Lr3,         // The dimension-less period.
       //tau    = 0.08*q/Lr3,       // End of pulse, dimension-less.
       k      = Period/tmstps,      // Length of a timestep.
       Deltat = Period/plts,        // Interval between each point plottet.
       p0     = 0.0;//-10/rho/g/Lr*conv;    // Ensures a certain diastolic pressure. (Unstressed pressure)

double       *fjac[18],   // Work space used by bound_bif.
             *fj[8];      // Work space used by bound_match.

#endif
