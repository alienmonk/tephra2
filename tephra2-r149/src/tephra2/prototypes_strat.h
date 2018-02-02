#include "parameters_strat.h"

/* define the following functions */

int get_eruptions(void);
int get_points(FILE *in_points);
int get_wind(FILE *in_wind);
int init_globals(char *in_config);
void set_global_values(FILE *log_file);
void set_eruption_values(ERUPTION *erupt, WIND *wind);

double phi2m( double xx );
double particle_density (double phi_slice);
double part_diff_time( double col_ht );
double (*pdf)(double, int, double, double);
double plume_pdf0(double x, int slice, double none, double ht_section_width);
double plume_pdf1(double x, int slice, double beta, double none);
double part_fall_time(double col_ht, double layer, double ashdiam, double part_density);
double pdf_grainsize(double part_mean_size, double part_size_slice, double part_step_width);

//double strat_average(WIND *level, double col_ht, double xspace, double yspace, double fall_time, double sigma, double elev);
double strat_average( double average_wind_direction, 
                      double average_wind_speed,             
			                double xspace, double yspace, 
			                double total_fall_time,
			                double sigma);
void tephra_calc(ERUPTION *erupt, POINT *pt, WIND *level, STATS *stats);
double part_fall_time_vent2ground (double vent_height, double layer,  double ashdiam,double part_density);
