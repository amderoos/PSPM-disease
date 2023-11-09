/***
  NAME
    PSPM-Disease.h
  PURPOSE
    Module computes the internal equilibrium of the the basic size-structured 
    consumer-resource model used in:

    A.M. de Roos, J.A.J. Metz & L. Persson, 2013. 
    Ontogenetic symmetry and asymmetry in energetics. J. Math. Biol. 66 (4): 889-914.

    with disease added.

    Last modification: AMdR - June 24, 2023
 ***/

#include "globals.h"


/*
 *===========================================================================
 *     DEFINITION OF PROBLEM DIMENSIONS AND NUMERICAL SETTINGS
 *===========================================================================
 */
// Dimension settings: Required
#define EQUATIONS_DIM             3
#define EXTRAOUTPUT_DIM           12
#define PARAMETER_NR              17

// Numerical settings: Optional (default values adopted otherwise)
#define DYTOL                     1.0E-7                                            // Variable tolerance
#define RHSTOL                    1.0E-8                                            // Function tolerance
#define ALLOWNEGATIVE             0                                                 // Negative solution components allowed?

#define JACOBIAN_STEP             1.0E-4

/*
 *===========================================================================
 *     DEFINITION OF ALIASES
 *===========================================================================
 */
#define RHO                       parameter[ 0]
#define RMAX                      parameter[ 1]

#define SB                        parameter[ 2]
#define SM                        parameter[ 3]

#define M                         parameter[ 4]

#define Q                         parameter[ 5]
#define H                         parameter[ 6]

#define SIGMA                     parameter[ 7]

#define TS                        parameter[ 8]
#define TI                        parameter[ 9]
#define TIJ                       parameter[10]
#define TIA                       parameter[11]

#define MUS                       parameter[12]
#define MUI                       parameter[13]
#define MUIJ                      parameter[14]
#define MUIA                      parameter[15]

#define BETA                      parameter[16]


/*
 *===========================================================================
 *     DEFINITION OF NAMES AND DEFAULT VALUES OF THE PARAMETERS
 *===========================================================================
 */
// At least two parameters should be specified in this array
char  *parameternames[PARAMETER_NR] =
    { "RHO", "Rmax", "Sb", "Sm", "M", "Q", "H", "SIGMA",
      "TS", "TI", "TIJ", "TIA",
      "MUS", "MUI", "MUIJ", "MUIA", "BETA"};

// These are the default parameters values
double  parameter[PARAMETER_NR] =
    { 0.1, 100.0, 0.1, 1.0, 1.0, 1.0, 3.0, 0.5,
      0.1, 0.0, 0.0, 0.0, 
      0.015, 0.0, 0.0, 0.0, 5.0E-4};

/*
 *===========================================================================
 *     DEFINITION OF THE SYSTEM OF EQUATIONS TO SOLVE
 *===========================================================================
 */

// Routine specifying the system of equalities from which to solve for
// R, B and IT at equilibrium
static double                     nu_J, nu_A, nu_IJ, nu_IA, nu_IJplus, nu_IAplus;
static double                     mu_J, mu_A, mu_IJ, mu_IA;
static double                     Juv_bio_S, Adu_bio_S, Juv_bio_I, Adu_bio_I;
static double                     Juv_nr_S, Adu_nr_S, Juv_nr_I, Adu_nr_I;
static double                     PCBirthRate_S, PCBirthRate_I, TotMatRate_S, TotMatRate_I;

#define R                         argument[0]                                       // Resource
#define B                         argument[1]                                       // Total birth rate in numbers
#define ITOT                      argument[2]                                       // Transmission

#define SMOOTHINGWIDTH            1.0E-3

// smooth maximum function

double smax(double a, double b)
  {
    double h1 = max(SMOOTHINGWIDTH   - fabs(a - b), 0.0);
    double h2 = max(1.0 - fabs(a - b) / SMOOTHINGWIDTH, 0.0);
    return (max(a, b) + 0.25 * h1 * h2);
  }

int    Equations(double *argument, double *result)

{
  double                          FR;
  double                          x, powzx, z;
  double                          TRANS, common = 0.0;

  TRANS = BETA * ITOT;
  // if ((CurveType == EQ) && (ITOT < 1.0E-4)) TRANS = -DYTOL;

    
  // Define the global variables
  z          = SB / SM;
  FR         = M * R / (H + R);

  nu_J       = SIGMA * (2 - Q) * FR - TS;
  nu_A       = SIGMA *      Q  * FR - TS;

  mu_J       = MUS;
  mu_A       = MUS;

  x          = (mu_J + TRANS) / nu_J;
  powzx      = pow(z, x);

  Juv_nr_S   = B * (1 - powzx) / (mu_J + TRANS);
  Juv_bio_S  = (B * SB / nu_J) * (1 - pow(z, x - 1)) / (x - 1);
  Adu_nr_S   = B * powzx / (mu_A + TRANS);
  Adu_bio_S  = SM * Adu_nr_S;

  TotMatRate_S  = B * powzx;
  PCBirthRate_S = nu_A * ((powzx / z) / (mu_A + TRANS));

  if (BETA > DYTOL)
    {
      nu_IJ     = nu_J - (TI + TIJ);
      nu_IA     = nu_A - (TI + TIA);

#if (0)
      nu_IJplus = max(nu_IJ, 0);
      nu_IAplus = max(nu_IA, 0);
#else
      nu_IJplus = smax(nu_IJ, 0);
      nu_IAplus = smax(nu_IA, 0);
#endif
      mu_IJ     = MUS + MUI + MUIJ;
      mu_IA     = MUS + MUI + MUIA;

      mu_IJ    -= nu_IJ - nu_IJplus;
      mu_IA    -= nu_IA - nu_IAplus;

      Juv_nr_I  = Juv_nr_S;
      Juv_bio_I = Juv_bio_S;

      if (nu_IJplus > 0)
        {
          double y     = mu_IJ / nu_IJplus;
          double powzy = pow(z, y);

          common       = nu_IJplus * (1 - powzy / powzx) * (mu_A + TRANS);
          common      /= nu_J * (mu_IJ - nu_IJplus * x);

          Juv_nr_I    += nu_IJplus * B *      (powzy - powzx) / (nu_J * (mu_IJ - nu_IJplus * x));
          Juv_bio_I   -= nu_IJplus * B * SB * (1 - powzy / z) / (nu_J * (mu_IJ - nu_IJplus    ));
       }

      Juv_nr_I      *= TRANS / mu_IJ;
      Juv_bio_I     *= TRANS / (mu_IJ - nu_IJplus * x);
      Adu_nr_I       = TRANS * (1 + common) * Adu_nr_S  / mu_IA;
      Adu_bio_I      = TRANS * (1 + common) * Adu_bio_S / mu_IA;

      PCBirthRate_I  = nu_IAplus * TRANS * (1 + common) / mu_IA;
      PCBirthRate_I *= ((powzx / z) / (mu_A + TRANS));
      TotMatRate_I   = TRANS * common * Adu_nr_S;
    }
  else
    {
      Juv_bio_I     = 0.0;
      Adu_bio_I     = 0.0;
      Juv_nr_I      = 0.0;
      Adu_nr_I      = 0.0;
      PCBirthRate_I = 0.0;
    }

  //==============================================================================
  // Compute the final values of the fixed point equation F(y)=0,

  if (result)
    {
      result[0]   = RHO * (RMAX - R);
      result[0]  -= M * FR * ((2 - Q) * (Juv_bio_S + Juv_bio_I) + Q * (Adu_bio_S + Adu_bio_I));
      result[1]   = (PCBirthRate_S + PCBirthRate_I - 1.0);
      if (ITOT > 1.0)
        result[2] = ((Juv_nr_I + Adu_nr_I) / ITOT - 1.0);
      else
        result[2] = ((Juv_nr_I + Adu_nr_I) - ITOT);
      // result[2]  = (log10(max(Juv_nr_I + Adu_nr_I,1.0E-5)) -  log10(max(ITOT,1.0E-5)));
    }

  return SUCCES;
}


/*===========================================================================*/

// Define all variables to be written to the output file (column-organized ASCII file)

int    DefineExtraOutput(double *argument, double *ExtraOutput)

{
  // Invoke the routine that sets the right-hand side for setting the output variables
  if (Equations(argument, NULL) == FAILURE) return FAILURE;

  ExtraOutput[ 0] = Juv_nr_S;
  ExtraOutput[ 1] = Adu_nr_S;
  ExtraOutput[ 2] = Juv_nr_I;
  ExtraOutput[ 3] = Adu_nr_I;
  ExtraOutput[ 4] = Juv_bio_S;
  ExtraOutput[ 5] = Adu_bio_S;
  ExtraOutput[ 6] = Juv_bio_I;
  ExtraOutput[ 7] = Adu_bio_I;
  ExtraOutput[ 8] = PCBirthRate_S;
  ExtraOutput[ 9] = PCBirthRate_I;
  ExtraOutput[10] = TotMatRate_S;
  ExtraOutput[11] = TotMatRate_I;

  return SUCCES;
}


/*==============================================================================*/
