#include "common_structures_strat.h"

/* Define the physical constants used in the new_tephra.c code
 * units are SI except WIND_INTERVAL (in km)
 * MAX_LINE, COL_STEPS and PART_STEPS dimensionless
*/
#define LOG_FILE "node_"
/* #define _PRINT 1 */
#define pi 3.141592654
#define DEG2RAD 0.017453293
/*#define DEBUG 0 */
#define MAX_LINE 200

/*mixed diffusion model*/
/*eddy diff for small particles in m2/s (400 cm2/s) */
extern double EDDY_CONST;

/* diffusion coeff for large particles (m2/s) */
extern double DIFFUSION_COEFFICIENT; 

/*threshold for change in diffusion (seconds fall time) */
extern double FALL_TIME_THRESHOLD; 


/* density of air */
/* air density in kg/m3 */
#define AIR_DENSITY 1.293
   
/* dynamic viscosity of air */
#define AIR_VISCOSITY 0.000018325
#define GRAVITY 9.81

/*density model for the pyroclasts */
/* These are now defined in config file */
extern double LITHIC_DENSITY;
extern double PUMICE_DENSITY;

#define LITHIC_DIAMETER_THRESHOLD 7.0
#define PUMICE_DIAMETER_THRESHOLD -1.0

/* #define BETA_LIMIT 100.0, now same as COL_STEPS */

enum{DISTR_1, DISTR_2};

/* These are now defined in config file. */
extern int PLUME_MODEL;
extern int PART_STEPS;
extern int COL_STEPS;
extern double PLUME_RATIO; /* Hb/Ht[area_of_release] of the laterally spreading cloud */
extern int WIND_DAYS;
extern int WIND_COLUMNS;

/* wind data read in at 0.5km intervals - this can be changed to suit available data
 * make sure input file intervals agree */


