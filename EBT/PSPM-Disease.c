/***
  NAME
    PSPM-Disease.c
  DESCRIPTION
    EBT implementation of the basic size-structured consumer-resource model
    used in:

    A.M. de Roos, J.A.J. Metz & L. Persson, 2013. 
    Ontogenetic symmetry and asymmetry in energetics. J. Math. Biol. 66 (4): 889-914.

    with disease added.

    Implementation without events and without boundary cohorts for the infected population

    OUTPUT VARIABLES OF THE STRUCTURED MODEL
       0 ( 2): R
       1 ( 3): Uninfected Juvenile biomass
       2 ( 4): Infected Juvenile biomass
       3 ( 5): Total Juvenile biomass

       4 ( 6): Uninfected Adult biomass
       5 ( 7): Infected Adult Biomass
       6 ( 8): Total Adult biomass

       7 ( 9): Total uninfected biomass
       8 (10): Total infected biomass
       9 (11): Total biomass

      10 (12): Uninfected Juvenile number
      11 (13): Infected Juvenile number
      12 (14): Total Juvenile number

      13 (15): Uninfected Adult number
      14 (16): Infected Adult number
      15 (17): Total Adult number

      16 (18): Total uninfected number
      17 (19): Total infected number
      18 (20): Total number

      19 (21): ingest_R
      20 (22): nu_J
      21 (23): nu_IJ
      22 (24): nu_A
      23 (25): nu_IA
      24 (26): mort_J
      25 (27): mort_IJ
      26 (28): mort_A
      27 (29): mort_IA

      28 (30): Uninfected adult fecundity (numbers)
      29 (31): Infected adult fecundity (numbers)

      30 (32): Uninfected birthrate (biomass)
      31 (33): Infected birthrate (biomass)

      32 (34): Number of cohorts in susceptible population
      33 (35): Number of cohorts in infected population

      34 (36): Current cohort limit

   Last modification: AMdR - July 11, 2023
***/
#include "escbox.h"

#ifndef LOGISTIC
#define LOGISTIC                  0                                                 // 0: Semi-chemostat; 1: Logistic
#endif


/*
 *====================================================================================================================================
 *
 *    LABELLING ENVIRONMENT AND I-STATE VARIABLES
 *
 *====================================================================================================================================
 */
#define SUSCEPT                   0
#define INFECTED                  1

#define time                      env[0]
#define R                         env[1]

#define age                       (i_state(0))
#define lnmass                    (i_state(1))

#define IDmaturity                (i_const(0))

#define TINY                      1.0E-6
#define EQUALINMASS               1.0E-4

#define ismature(p, n)            ((popIDcard[(p)][(n)][IDmaturity] > 0.5) || (pop[(p)][(n)][lnmass] >= logsm))
#define isequalmass(p, n, m)      (fabs(pop[(p)][(n)][lnmass] - pop[(p)][(m)][lnmass]) < EQUALINMASS * fabs(pop[(p)][(n)][lnmass] + pop[(p)][(m)][lnmass])) 

/*
 *====================================================================================================================================
 *
 *    DEFINING AND LABELLING CONSTANTS AND PARAMETERS
 *
 *====================================================================================================================================
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

#define LAKEVOLUME                1.0E3


/*
 *====================================================================================================================================
 *
 *    OTHER DEFINITIONS
 *
 *====================================================================================================================================
 */

static double                     Adu_S_bio, Adu_I_bio, Juv_S_bio, Juv_I_bio;
static double                     Adu_S_nr,  Adu_I_nr,  Juv_S_nr,  Juv_I_nr;
static double                     ingest_R;
static double                     TRANS;
static double                     nu_J, nu_IJ, nu_A, nu_IA, nuplus_J, nuplus_IJ, nuplus_A, nuplus_IA;
static double                     mort_J, mort_A, mort_IJ, mort_IA;
static double                     S_birthrate, I_birthrate;
static double                     logsb, logsm;
static double                     IntroInfectedNewborns = 0, IntroInfectedTime = 200;
static int                        NewBirths;
#if (SUSCEPTIBLEONLY == 0)
static int                        NewInfections;
#endif
static int                        InGradient = 0;

#if (BIFURCATION == 1)
#define MIN_S_COHORTS             200
#define MAX_S_COHORTS             800
#define MIN_I_COHORTS             300
#define MAX_I_COHORTS             1500

#define COHORTLIMITS              13
static double                     CohLimits[COHORTLIMITS] = {0.001, 0.002, 0.0025, 0.005,
                                                             0.01,  0.02,  0.025,  0.05,
                                                             0.1,   0.2,   0.25,   0.5,
                                                             1.0};
extern double                     AllAve[OUTPUT_VAR_NR], AllMin[OUTPUT_VAR_NR], AllMax[OUTPUT_VAR_NR];
#endif

/*
 *====================================================================================================================================
 *
 * USER INITIALIZATION ROUTINE ALLOWS OPERATIONS ON INITIAL POPULATIONS
 *
 *====================================================================================================================================
 */

void UserInit(int argc, char **argv, double *env, population *pop)

{
  int                             ii;

  switch (argc)
    {
      case 4:
        IntroInfectedTime     = atof(argv[3]);
      case 3:
        IntroInfectedNewborns = atof(argv[2]);
      default:
        break;
    }

  logsb = log(SB);
  logsm = log(SM) - logsb;

  for (ii = 0; ii < cohort_no[SUSCEPT]; ii++)
    {
      if (pop[SUSCEPT][ii][lnmass] > logsm * (1.0 - TINY))
        popIDcard[SUSCEPT][ii][IDmaturity] = 1.0;
      else
        popIDcard[SUSCEPT][ii][IDmaturity] = 0.0;
    }

#if (SUSCEPTIBLEONLY == 0)
  for (ii = 0; ii < cohort_no[INFECTED]; ii++)
    {
      if (pop[INFECTED][ii][lnmass] > logsm * (1.0 - TINY))
        popIDcard[INFECTED][ii][IDmaturity] = 1.0;
      else
        popIDcard[INFECTED][ii][IDmaturity] = 0.0;
    }
#endif


#if (LOGISTIC > 0)
  ReportNote("Logistic resource dynamics");
#else
  ReportNote("Semi-chemostat resource dynamics");
#endif

  ReportNote("Disease transmission proportional to numerical abundance");
  ReportNote("\n    Ratio of offspring-adult body size: z = %4.2f\n", SB / SM);

#if (BIFURCATION == 1)
  ReportNote("Cohort limit adjusted to target %d - %d cohorts in the suceptible and %d - %d in the infected population", 
             MIN_S_COHORTS, MAX_S_COHORTS, MIN_I_COHORTS, MAX_I_COHORTS);
#endif                                                                              // BIFURCATION

  return;
}


/*
 *====================================================================================================================================
 *
 *  SPECIFICATION OF THE NUMBER AND VALUES OF BOUNDARY POINTS
 *
 *====================================================================================================================================
 */

void SetBpointNo(double *env, population *pop, int *bpoint_no)

{
  bpoint_no[SUSCEPT]  = 0;

  NewBirths     = AddCohorts(pop, SUSCEPT, 1);                                      // Add new offspring cohort

  pop[SUSCEPT][NewBirths][number]             = 0.0;
  pop[SUSCEPT][NewBirths][age]                = 0.0;
  pop[SUSCEPT][NewBirths][lnmass]             = 0.0;
  popIDcard[SUSCEPT][NewBirths][IDmaturity]   = 0;

#if (SUSCEPTIBLEONLY == 0)
  int                             ii, jj;

  bpoint_no[INFECTED] = 0;

  NewInfections = AddCohorts(pop, INFECTED, cohort_no[SUSCEPT]);                    // Add new cohorts to be infected
  for (ii = NewInfections, jj = 0; ii < cohort_no[INFECTED]; ii++, jj++)
    {
      pop[INFECTED][ii][number]            = 0.0;
      pop[INFECTED][ii][age]               = pop[SUSCEPT][jj][age];
      pop[INFECTED][ii][lnmass]            = pop[SUSCEPT][jj][lnmass];
      popIDcard[INFECTED][ii][IDmaturity]  = popIDcard[SUSCEPT][jj][IDmaturity];
    }

  if ((!iszero(IntroInfectedNewborns)) && isequal(IntroInfectedTime, time))
    {
      NewInfections = AddCohorts(pop, INFECTED, 1);
      pop[INFECTED][NewInfections][number]            = IntroInfectedNewborns;
      pop[INFECTED][NewInfections][age]               = 0.0;
      pop[INFECTED][NewInfections][lnmass]            = logsb;
      popIDcard[INFECTED][NewInfections][IDmaturity]  = 0.0;
    }
#endif

  return;
}


/*===============================================================================================================================*/

void SetBpoints(double *env, population *pop, population *bpoints)

{
  return;
}


/*===============================================================================================================================*/

void SetStructVars(double *env, population *pop)

{
  int                             ii;
  double                          FR;

  Juv_S_nr  = 0.0;
  Adu_S_nr  = 0.0;
  Juv_S_bio = 0.0;
  Adu_S_bio = 0.0;
  Juv_I_nr  = 0.0;
  Adu_I_nr  = 0.0;
  Juv_I_bio = 0.0;
  Adu_I_bio = 0.0;
  TRANS     = 0.0;
  ingest_R  = 0.0;

  for (ii = 0; ii < cohort_no[SUSCEPT]; ii++)
    {
      if (ismature(SUSCEPT, ii))
        {
          Adu_S_nr  += pop[SUSCEPT][ii][number];
          Adu_S_bio += SM * pop[SUSCEPT][ii][number];
        }
      else
        {
          Juv_S_nr  += pop[SUSCEPT][ii][number];
          Juv_S_bio += SB * exp(pop[SUSCEPT][ii][lnmass]) * pop[SUSCEPT][ii][number];
        }
    }

#if (SUSCEPTIBLEONLY == 0)
  for (ii = 0; ii < cohort_no[INFECTED]; ii++)
    {
      if (ismature(INFECTED, ii))
        {
          Adu_I_nr  += pop[INFECTED][ii][number];
          Adu_I_bio += SM * pop[INFECTED][ii][number];
        }
      else
        {
          Juv_I_nr  += pop[INFECTED][ii][number];
          Juv_I_bio += SB * exp(pop[INFECTED][ii][lnmass]) * pop[INFECTED][ii][number];
        }
    }
#endif

#if (SUSCEPTIBLEONLY == 1)
  TRANS      = 0.0;
#else
#if (BIFURCATION == 1)
  TRANS      = BETA * max(Juv_I_nr + Adu_I_nr, 10.0) / LAKEVOLUME;
#else
  TRANS      = BETA * (Juv_I_nr + Adu_I_nr) / LAKEVOLUME;
#endif
#endif

  FR         = M * R / (H + R);

  nu_J       = SIGMA * (2 - Q) * FR - TS;
  nu_A       = SIGMA *      Q  * FR - TS;
  ingest_R   = ((2 - Q) * Juv_S_bio + Q * Adu_S_bio) * FR;

  nu_IJ      = SIGMA * (2 - Q) * FR - (TS + TI + TIJ);
  nu_IA      = SIGMA *      Q  * FR - (TS + TI + TIA);
  ingest_R  += ((2 - Q) * Juv_I_bio + Q * Adu_I_bio) * FR;

  nuplus_J   = max(nu_J, 0.0);
  nuplus_A   = max(nu_A, 0.0);
  nuplus_IJ  = max(nu_IJ, 0.0);
  nuplus_IA  = max(nu_IA, 0.0);

  mort_J     = MUS;
  mort_A     = MUS;

  mort_IJ    = MUS + MUI + MUIJ;
  mort_IA    = MUS + MUI + MUIA;

  mort_J    -= min(nu_J, 0.0);
  mort_A    -= min(nu_A, 0.0);

  mort_IJ   -= min(nu_IJ, 0.0);
  mort_IA   -= min(nu_IA, 0.0);

  S_birthrate = nuplus_A  * Adu_S_bio;
  I_birthrate = nuplus_IA * Adu_I_bio;

  return;
}


/*
 *====================================================================================================================================
 *
 *      SPECIFICATION OF DERIVATIVES
 *
 *====================================================================================================================================
 */

void Gradient(double *env, population *pop, population *ofs, 
              double *envgrad, population *popgrad, population *ofsgrad, population *bpoints)

{
  int                             ii;
  double                          repro = 0.0;

  SetStructVars(env, pop);

#if (BIFURCATION == 1)
  double                          tmpout[OUTPUT_VAR_NR];

  if ((rk_level == 1) && (fmod(env[0], BifPeriod) > 0.6*BifPeriod))
    {
      InGradient = 1;
      DefineOutput(env, pop, tmpout);
      InGradient = 0;
      for (ii = 0; ii < OUTPUT_VAR_NR; ii++)
        {
          if (tmpout[ii] > AllMax[ii]) AllMax[ii] = tmpout[ii];
          if (tmpout[ii] < AllMin[ii]) AllMin[ii] = tmpout[ii];
        }
    }
#endif

  for (ii = 0; ii < cohort_no[SUSCEPT]; ii++)
    {
      popgrad[SUSCEPT][ii][number] = -1.0;
      popgrad[SUSCEPT][ii][age]    =  1.0;
      popgrad[SUSCEPT][ii][lnmass] =  0.0;

      // We assume that individuals never shrink
      // Starvation is dealt with as additional mortality

      if (ismature(SUSCEPT, ii))
        {
          popgrad[SUSCEPT][ii][number]  = -(mort_A + TRANS) * pop[SUSCEPT][ii][number];
          popgrad[SUSCEPT][ii][lnmass]  = 0.0;
        }
      else
        {
          popgrad[SUSCEPT][ii][number]  = -(mort_J + TRANS) * pop[SUSCEPT][ii][number];
          popgrad[SUSCEPT][ii][lnmass]  = nuplus_J;
        }
    }

  // The derivatives for the susceptible cohort of newborn individuals
  repro = (S_birthrate + I_birthrate) / SB;
  popgrad[SUSCEPT][NewBirths][number] += repro;

  if (pop[SUSCEPT][NewBirths][number] > 1.0)
    {
      popgrad[SUSCEPT][NewBirths][age]    -= repro * pop[SUSCEPT][NewBirths][age]    / pop[SUSCEPT][NewBirths][number];
      popgrad[SUSCEPT][NewBirths][lnmass] -= repro * pop[SUSCEPT][NewBirths][lnmass] / pop[SUSCEPT][NewBirths][number];
    }

#if (SUSCEPTIBLEONLY == 0)
  int                             jj;

  for (ii = 0; ii < NewInfections; ii++)
    {
      popgrad[INFECTED][ii][number] = -1.0;
      popgrad[INFECTED][ii][age]    =  1.0;
      popgrad[INFECTED][ii][lnmass] =  0.0;

      if (ismature(INFECTED, ii))
        {
          popgrad[INFECTED][ii][number]  = -mort_IA * pop[INFECTED][ii][number];
          popgrad[INFECTED][ii][lnmass]  = 0.0;
        }
      else
        {
          popgrad[INFECTED][ii][number]  = -mort_IJ * pop[INFECTED][ii][number];
          popgrad[INFECTED][ii][lnmass]  = nuplus_IJ;
        }
    }

  for (ii = NewInfections, jj = 0; ii < cohort_no[INFECTED]; ii++, jj++)
    {
      popgrad[INFECTED][ii][number] = -1.0;
      popgrad[INFECTED][ii][age]    =  1.0;
      popgrad[INFECTED][ii][lnmass] =  0.0;
      
      repro = TRANS * pop[SUSCEPT][jj][number];

      if (ismature(INFECTED, ii))
        {
          popgrad[INFECTED][ii][number]  = -mort_IA * pop[INFECTED][ii][number] + repro;
          popgrad[INFECTED][ii][lnmass]  = 0.0;
        }
      else
        {
          popgrad[INFECTED][ii][number]  = -mort_IJ * pop[INFECTED][ii][number] + repro;
          popgrad[INFECTED][ii][lnmass]  = nuplus_IJ;
        }

      if (pop[INFECTED][ii][number] > 1.0)
        {
          popgrad[INFECTED][ii][lnmass] += repro * (pop[SUSCEPT][jj][lnmass] - pop[INFECTED][ii][lnmass]) / pop[INFECTED][ii][number];
        }
    }
#endif

  envgrad[0] = 1.0;
#if (LOGISTIC > 0)
  envgrad[1] = RHO*R*(1.0 - R/RMAX) - ingest_R/LAKEVOLUME;
#else
  envgrad[1] = RHO * (RMAX - R) - ingest_R / LAKEVOLUME;
#endif

  return;
}


/*
 *====================================================================================================================================
 *
 *    SPECIFICATION OF BETWEEN COHORT CYCLE DYNAMICS
 *
 *====================================================================================================================================
 */

void InstantDynamics(double *env, population *pop, population *ofs)

{
  int                             ii;
  int                             first_adult = -1;
  double                          total_number, total_age;

  // Lump all susceptible adult cohorts
  first_adult = -1; total_number = total_age = 0.0;
  for (ii = 0; ii < cohort_no[SUSCEPT]; ii++)
    {
      if ((pop[SUSCEPT][ii][number] < 0) || iszero(pop[SUSCEPT][ii][number])) continue;

      if (ismature(SUSCEPT, ii))
        {
          if (first_adult == -1) first_adult = ii;

          total_number += pop[SUSCEPT][ii][number];
          total_age    += pop[SUSCEPT][ii][age] * pop[SUSCEPT][ii][number];
          pop[SUSCEPT][ii][number] = 0.0;

          pop[SUSCEPT][first_adult][number]           = total_number;
          pop[SUSCEPT][first_adult][lnmass]           = logsm;
          pop[SUSCEPT][first_adult][age]              = total_age / total_number;
          popIDcard[SUSCEPT][first_adult][IDmaturity] = 1;
        }
    }

#if (SUSCEPTIBLEONLY == 0)
  int                             jj;
  int                             merge_cohorts = 0;
  double                          total_lnmass;

  // Lump all similar infected cohorts
  first_adult = -1;
  for (ii = 0; ii < cohort_no[INFECTED]; ii++)
    {
      if ((pop[INFECTED][ii][number] < 0) || iszero(pop[INFECTED][ii][number])) continue;

      if ((first_adult == -1) && ismature(INFECTED, ii)) first_adult = ii;

      total_number  = pop[INFECTED][ii][number];
      total_age     = pop[INFECTED][ii][number] * pop[INFECTED][ii][age];
      total_lnmass  = pop[INFECTED][ii][number] * pop[INFECTED][ii][lnmass];
      merge_cohorts = 0;

      for (jj = ii + 1; jj < cohort_no[INFECTED]; jj++)
        {
          if ((pop[INFECTED][jj][number] < 0) || iszero(pop[INFECTED][jj][number])) continue;

          if ((first_adult == ii) && (ismature(INFECTED, jj) || (isequalmass(INFECTED, ii, jj))))
            {
              total_number  += pop[INFECTED][jj][number];
              total_age     += pop[INFECTED][jj][number] * pop[INFECTED][jj][age];
              merge_cohorts  = 1;

              pop[INFECTED][jj][number] = 0.0;
            }
          else if (isequalmass(INFECTED, ii, jj))
            {
              total_number  += pop[INFECTED][jj][number];
              total_age     += pop[INFECTED][jj][number] * pop[INFECTED][jj][age];
              total_lnmass  += pop[INFECTED][jj][number] * pop[INFECTED][jj][lnmass];
              merge_cohorts  = 1;

              pop[INFECTED][jj][number] = 0.0;
            }
        }
      if (merge_cohorts)
        {
          pop[INFECTED][ii][number]            = total_number;
          pop[INFECTED][ii][age]               = total_age    / total_number;
          pop[INFECTED][ii][lnmass]            = total_lnmass / total_number;
        }
      if (first_adult == ii)
        {
          pop[INFECTED][ii][lnmass]            = logsm;
          popIDcard[INFECTED][ii][IDmaturity]  = 1;
        }
    }
#endif

  return;
}


/*
 *====================================================================================================================================
 *
 *      SPECIFICATION OF OUTPUT VARIABLES
 *
 *====================================================================================================================================
 */

void DefineOutput(double *env, population *pop, double *output)

{
  int                             outnr = 0;

  if (!InGradient)
    {
      LabelState(0, "T = %.2f", env[0]);
      SetStructVars(env, pop);
    }

#ifdef RPACKAGE
  if (iszero(fmod(env[0], 1.0))) fprintf(stderr, ".");
#endif

  output[outnr++] = R;                                                              //  2:
  output[outnr++] = Juv_S_bio / LAKEVOLUME;                                         //  3: Uninfected juvenile biomass
  output[outnr++] = Juv_I_bio / LAKEVOLUME;                                         //  4: Infected juvenile biomass
  output[outnr++] = (Juv_I_bio + Juv_S_bio) / LAKEVOLUME;                           //  5: Total juvenile biomass

  output[outnr++] = Adu_S_bio / LAKEVOLUME;                                         //  6: Uninfected adult biomass
  output[outnr++] = Adu_I_bio / LAKEVOLUME;                                         //  7: Infected adult biomass
  output[outnr++] = (Adu_S_bio + Adu_I_bio) / LAKEVOLUME;                           //  8: Total adult biomass

  output[outnr++] = (Juv_S_bio + Adu_S_bio) / LAKEVOLUME;                           //  9: Total uninfected biomass
  output[outnr++] = (Juv_I_bio + Adu_I_bio) / LAKEVOLUME;                           // 10: Total infected biomass
  output[outnr++] = (Juv_S_bio + Juv_I_bio + Adu_S_bio + Adu_I_bio) / LAKEVOLUME;   // 11: Total consumer biomass

  output[outnr++] = Juv_S_nr / LAKEVOLUME;                                          // 12: Uninfected juvenile number
  output[outnr++] = Juv_I_nr / LAKEVOLUME;                                          // 13: Infected juvenile number
  output[outnr++] = (Juv_S_nr + Juv_I_nr) / LAKEVOLUME;                             // 14: Total juvenile number

  output[outnr++] = Adu_S_nr / LAKEVOLUME;                                          // 15: Uninfected adult number
  output[outnr++] = Adu_I_nr / LAKEVOLUME;                                          // 16: Infected adult number
  output[outnr++] = (Adu_S_nr + Adu_I_nr) / LAKEVOLUME;                             // 17: Total adult number

  output[outnr++] = (Juv_S_nr + Adu_S_nr) / LAKEVOLUME;                             // 18: Total uninfected number
  output[outnr++] = (Juv_I_nr + Adu_I_nr) / LAKEVOLUME;                             // 19: Total infected number
  output[outnr++] = (Juv_S_nr + Juv_I_nr + Adu_S_nr + Adu_I_nr) / LAKEVOLUME;       // 20: Total consumer number

  output[outnr++] = ingest_R / LAKEVOLUME;                                          // 21:
  output[outnr++] = nu_J;                                                           // 22:
  output[outnr++] = nu_IJ;                                                          // 23:
  output[outnr++] = nu_A;                                                           // 24:
  output[outnr++] = nu_IA;                                                          // 25:
  output[outnr++] = mort_J;                                                         // 26:
  output[outnr++] = mort_IJ;                                                        // 27:
  output[outnr++] = mort_A;                                                         // 28:
  output[outnr++] = mort_IA;                                                        // 29:

  output[outnr++] = nuplus_A  * SM / SB;                                            // 30: Uninfected Fecundity (in number of offspring!!)
  output[outnr++] = nuplus_IA * SM / SB;                                            // 31: Infected Fecundity (in number of offspring!!)

  output[outnr++] = S_birthrate / LAKEVOLUME;                                       // 32: Uninfected biomass reproduction rate
  output[outnr++] = I_birthrate / LAKEVOLUME;                                       // 33: Infected biomass reproduction rate

  output[outnr++] = cohort_no[SUSCEPT];                                             // 34: Number of cohort in susceptible population
#if (SUSCEPTIBLEONLY == 1)
  output[outnr++] = 0;
#else
  output[outnr++] = cohort_no[INFECTED];                                            // 35: Number of cohort in infected population
#endif
  output[outnr++] = cohort_limit;                                                   // 36: Current cohort limit

#if (BIFURCATION == 1)
  if (!InGradient) 
    {
      if ((fmod(env[0], BifPeriod) < 0.1 * cohort_limit) || (fmod(env[0], BifPeriod) > (BifPeriod - 0.1 * cohort_limit)))
        {
          double                  Ave_CohNr_S = AllAve[32];
#if (SUSCEPTIBLEONLY == 0)
          double                  Ave_CohNr_I = AllAve[33];
#endif
          int                     clindex;

          // Find the current index of the cohort limit in CohLimits
          for (clindex = 0; clindex < COHORTLIMITS; clindex++)
            {
              if (fabs(cohort_limit - CohLimits[clindex]) < 0.01 * cohort_limit) break;
            }
          if (clindex >= COHORTLIMITS) clindex = COHORTLIMITS - 1;
          if (clindex <= 0)            clindex = 0;
#if (SUSCEPTIBLEONLY == 1)
          if (Ave_CohNr_S > MAX_S_COHORTS) clindex++;
          if (Ave_CohNr_S < MIN_S_COHORTS) clindex--;
#else
          if ((Adu_I_nr + Juv_I_nr) < 10.0)
            {
              if (Ave_CohNr_S > MAX_S_COHORTS) clindex++;
              if (Ave_CohNr_S < MIN_S_COHORTS) clindex--;              
            }
          else 
            {
              if ((Ave_CohNr_S > MAX_S_COHORTS) || (Ave_CohNr_I > MAX_I_COHORTS)) clindex++;
              if ((Ave_CohNr_S < MIN_S_COHORTS) || (Ave_CohNr_I < MIN_I_COHORTS)) clindex--;
            }
#endif
          if (clindex >= COHORTLIMITS) clindex = COHORTLIMITS - 1;
          if (clindex <= 0)            clindex = 0;

          cohort_limit = CohLimits[clindex];
        }
    }
#endif

  return;
}


/*===============================================================================================================================*/
