// AMdR - May 19, 2023

#define SUSCEPTIBLEONLY           0                                                 // With (0) or without (1) infections

#if (SUSCEPTIBLEONLY == 1)
#define POPULATION_NR             1
#else
#define POPULATION_NR             2
#endif
#define I_STATE_DIM               2
#define I_CONST_DIM               1
#define ENVIRON_DIM               2

#define OUTPUT_VAR_NR             35
#define PARAMETER_NR              17

#define EVENT_NR                  0
#define TIME_METHOD               RKTSI5
#define DYNAMIC_COHORTS           0                                                 // Exclusively dynamic (0/1)

#if (BIFURCATION == 1)
#define ADJUST_COH_LIMIT          0
#define TARGET_COHORTS            500
#endif                                                                              // BIFURCATION

/*===============================================================================================================================*/
