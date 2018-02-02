#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <gc.h>
#include "prototypes_strat.h"

static double threshold;
static double AIR_VISCOSITY_x_18;
static double LITHIC_DENSITY_minus_PUMICE_DENSITY;
static double PUMICE_DIAMETER_THRESHOLD_minus_LITHIC_DIAMETER_THRESHOLD;
static double ONE_THIRD;
static double AIR_VISCOSITY_x_225;
static double GRAV_SQRD_x_4;
static double SQRT_TWO_PI; /* add new line, 12-2-2010 */
static double BETA_x_SQRT_TWO_PI;
static double TWO_BETA_SQRD;
static double PDF_GRAINSIZE_DEMON1;
static double TWO_x_PART_SIGMA_SIZE;
static double EDDY_CONST_x_8_div_5;
static TABLE **T;
static FILE *log_file;

/*
Code: tephra_calc.c
By: C.B. and L.J. Connor & T. Hincks & C. Bonadonna
Copyright (C) 2003  C.B. Connor, L.J. Connor, C. Bonadonna, and T. Hincks
See: http://www.cas.usf.edu/~cconnor/parallel/tephra/tephra.html

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

*/
/* ---------------------------------------------------------------------
 * FUNCTION:tephra_calc.c 
 *
 * Purpose: This function calculates and returns the expected accumulation
 * of volcanic ash (kg/m2) at a specific geogrphic location (x,y)
 * due to an eruption with specific input parameters. 
 * These points may be random or on  a UTM grid (m)
 * 
 * This implementation accounts for variation in wind velocity with height.
 * The model is discretized w.r.t. height and particle size. 

 * This function is called for each point (x,y,) If more than one eruption is
 * involved, for example in a probabilistic analysis, the function is called for
 * each set of eruption parameters. 
 
 * INPUTS:
 * ERUPTION *erupt: pointer to array of eruption parameters
 * POINT *pt: pointer to an array of location specific parameters, 
 * WIND *level: pointer to a day of wind data :
              height asl in m; 
              wind speed in ms-1
              wind direction in degrees N
 *	     
 * OUTPUTs:
 *   the value of the mass accumulated at the input location (northing, easting) in kg/m2 
 *
 *   a distribution of particle sizes
 *   the exact number of binss (i.e. sizes) and phi size used per bin is an integer and 
 *   is determined by (erupt->max_phi - erupt->min_phi)
 *   each bin accumulates phi sizes up to its integer size
 *   ex. bin[0] holds grainsizes [min_phi to min_phi+1)  
 ***************************************************************************/

void tephra_calc(ERUPTION *erupt, POINT *pt, WIND *level, STATS *stats) { /* tephra_calc starts ... */
 
   /**********************************************************************************
   * WIND structure:
   * level->windspeed: windspeed in m/s
   * level->wind_dir: wind direction in +/- degrees from north
   * level->wind_height: meters above sea level
   *
   * See common_structures_strat.h for structure definitions.
   **********************************************************************************/
  
  
	int i = 0,j = 0, bin = -1; 
	double new_xspace, new_yspace, ht_above_vent, cos_wind = 0.0, sin_wind = 0.0, windspeed = 0.0;
	double sigma, demon2, demon3, ash_fall, layer, fall_time_adj = 0.0, total_fall_time=0.0;
	double average_windspeed_x, average_windspeed_y, average_wind_direction, average_windspeed =0.0;
	double wind_sum_x = 0.0, wind_sum_y = 0.0;
	static double min=10e6, max=0.0;
  
#ifdef _PRINT
	fprintf(log_file, "IN tephra_calc ...");
#endif

  /* Initialize mass to zero */
	pt->mass = 0.0;
	wind_sum_x = 0.0;
	wind_sum_y = 0.0;
	
  /* Transform the volcano location coordinate to 0,0 
   */
	new_xspace = pt->northing - erupt->volcano_northing;
	new_yspace = pt->easting - erupt->volcano_easting;

  /* do the double integration over grainsize and column height 
   */
#ifdef _PRINT
	fprintf(log_file, "\nBeginning integration loops ... \n");
#endif
  
  /* Interpolate to fine the windspeed and direction below the height of the vent.
   * Find one average wind speed and wind direction between vent and grid elevation point. 
   * The first values in the wind array give the wind speed and direction at the
   * vent height.
   */
	layer = erupt->vent_height - pt->elevation;
  
	windspeed = (level[0].windspeed * pt->elevation) / erupt->vent_height;
	cos_wind = cos(level[0].wind_dir) * windspeed;
	sin_wind = sin(level[0].wind_dir) * windspeed;

	   
	for (i = 0; i < PART_STEPS; i++) { /* PART_STEPS_LOOP */
   
		fall_time_adj = 0.0;
		
		/* Accumulate the particle sizes into bins of whole numbered phi sizes */
		if (!(i % 10)) { 
			bin++;
#ifdef _PRINT
			fprintf(log_file, "PART_STEP=%d phi[%d] = %g\n", i, bin, pt->phi[bin]);
#endif
		}
  
		/* Modify the total fall time of each particle size (i) 
		   by the time that it takes particle size (i) to descend from vent height to the grid cell. 
		   This is only necessary when the elevation of the grid cell < elevation of the vent.
		*/           
		if (layer > 0) {
			fall_time_adj = 
			part_fall_time(erupt->vent_height, layer, T[i][0].ashdiam, T[i][0].part_density);
     
#ifdef DEBUG 	
			fprintf(log_file, "%d %g %g\n",  i, layer, fall_time_adj) ;
#endif	 
		}

		for (j = 0; j < COL_STEPS; j++) { /* COL_STEPS_LOOP */
     
  		total_fall_time = T[i][j].total_fall_time + fall_time_adj;
    	// fprintf(stderr, "%g %g %g ", T[i][j].total_fall_time, total_fall_time, fall_time_adj);
      
    	/* Sum the adjustments (windspeed and wind_direction) 
    	 * for each particle size  falling from each level.
    	 */
		/* removed 2 lines, 12-2-2010
   		wind_sum_x = cos_wind * fall_time_adj * windspeed;
	  	wind_sum_y = sin_wind * fall_time_adj * windspeed;
	  	*/
	  	/* change 2 lines, 12-2-2010 */
	  	wind_sum_x = cos_wind * fall_time_adj;
	  	wind_sum_y = sin_wind * fall_time_adj;
	    
	    
	    /* Now add the summed adjustments to the already summed
	     * windspeeds and directions 
	     and
    	 Account for the wind:
    	 Find the average windspeed in the x and y directions 
    	 over the total fall time.
			*/
			average_windspeed_x = 
			(T[i][j].wind_sum_x + wind_sum_x)/total_fall_time;
			
			average_windspeed_y = 
			(T[i][j].wind_sum_y + wind_sum_y)/total_fall_time;
    	    	
			/* If zero, make windspeed a very small value (cannot divide by zero in next step) */
			if (!average_windspeed_x) average_windspeed_x = .001;
			if (!average_windspeed_y) average_windspeed_y = .001;
      
			/* Find the average wind direction (direction of the velocity vector) */
			if (average_windspeed_x < 0) {
				average_wind_direction = 
				atan(average_windspeed_y/average_windspeed_x ) + pi;
			} else 
    	  average_wind_direction = 
    	  atan(average_windspeed_y/average_windspeed_x);
    	
			/* Find the average windspeed ( magnitude of the velocity vector) */
			average_windspeed = 
			sqrt(average_windspeed_x*average_windspeed_x + average_windspeed_y*average_windspeed_y);
				
			if (total_fall_time > max) max = total_fall_time;
			if (total_fall_time < min) min = total_fall_time;
				
			/* calculate the value of sigma (dispersion) based on total_fall_time  
			 * to acct for the change in the shape of the column with ht - increasing radius 
			 */
			ht_above_vent = T[i][j].particle_ht - erupt->vent_height;
			
			/* falltime for fine particles */
			if (total_fall_time >= FALL_TIME_THRESHOLD) {
				sigma = 
				EDDY_CONST_x_8_div_5 * pow((total_fall_time + T[i][j].plume_diffusion_fine_particle), 2.5);
				//fprintf(stderr,"f");
			} else { /* falltime for coarse particles */
				sigma = 
				4.0 * DIFFUSION_COEFFICIENT * (total_fall_time + T[i][j].plume_diffusion_coarse_particle);
				//fprintf(stderr, "c");
			}

		  demon2 =  pi * sigma;
      
    	/* Modify fall time by the variation of wind velocity with height */
			demon3 = 
			strat_average( average_wind_direction, 
      	             average_windspeed,             
				             new_xspace, new_yspace, 
				             total_fall_time,
				             sigma); 
/*
			if (!demon2 || isnan(demon2) || isinf(demon2) || isnan(demon3) || isinf(demon3)) {
 				fprintf(stderr, 
      	"[%d][%d] layer= %.1f totalfalltime=%g [falltimeadj=%g] demon1=%g demon2=%g demon3=%g sigma=%g\n",
      	i,j, layer,total_fall_time, fall_time_adj, T[i][j].demon1, demon2, demon3, sigma);
      	exit(-1);
			}
 */     
			ash_fall = (T[i][j].demon1 / demon2) * demon3;
			pt->mass += ash_fall;
			pt->phi[bin] += ash_fall;
			//fprintf(stderr, "\n");
		}    
	}
#ifdef _PRINT
  fprintf(log_file, "PART_STEP=%d phi[%d] = %g\n", i, bin, pt->phi[bin]);
  fprintf(log_file, "OUT\n");
#endif
  stats->min_falltime = min;
  stats->max_falltime = max;
}


/* ----------------- New Function Starts Here -------------------- */
/* function phi2m converts the ash diameter from 
   units of phi to m
*/

 double phi2m(double xx) {
   double cms;
   //printf("in phi2m");
   cms = 0.001 * pow(2, -xx);
   return cms;
 }
/* ----------------- New Function Starts Here -------------------- */
/* function particle_density calculates varying particle density based on their grain size diamete
   using a linear correlation between pumice_threshold (PHI) and lithic_threshold (PHI)
*/

 double particle_density (double phi_slice) {

  double mean_density = 0;

  if (phi_slice >= LITHIC_DIAMETER_THRESHOLD) mean_density = LITHIC_DENSITY;
  else if (phi_slice <= PUMICE_DIAMETER_THRESHOLD) mean_density = PUMICE_DENSITY;
  else if (phi_slice < LITHIC_DIAMETER_THRESHOLD && phi_slice > PUMICE_DIAMETER_THRESHOLD) 
    mean_density = 
    
      LITHIC_DENSITY - 
      LITHIC_DENSITY_minus_PUMICE_DENSITY * 
      (phi_slice - LITHIC_DIAMETER_THRESHOLD) / PUMICE_DIAMETER_THRESHOLD_minus_LITHIC_DIAMETER_THRESHOLD;

   return mean_density;
 }


/* ----------------- New Function Starts Here -------------------- */
/* function part_diff_time determines the particle diffusion time const
   in the atmosphere, where 
   col_ht = height of the particle (m) in the column w.r.t. vent. 
   THis is added to the particle fall time to account for the width of the column,
   which changes as a function of height

   c = eddy diffusivity in the atmosphere
   given in si units, e.g., 0.04 m^2 s^-1
   
   returns the particle diffusion time
*/


/* ----------------- New Function Starts Here -------------------- */
/* function part_fall_time determines the time of particle fall within each falling step
   falling steps are here:
   set = particle rising steps = ht_step_width
   
   returns the particle fall time within each falling step

This function follows the approach outlined in Bonadonna et al. (1998) 
Briefly, particle fall time is calculated based on terminal velocities in
layers that are 1000 m thick. The terminal velocity is a function of the
particle Reynolds number, which varies with grainsize, air properties.

The thickness of the first layer (closest to the ground) is equal to the vent 
height. The vent_height is in meters above sea level. The area the ash falls on
is considered to be at sea level. This leads to some assumptions (!) near the
volcano...
*/


double part_fall_time(double particle_ht, double layer, double ashdiam, double part_density) {
  
	double rho, hz, temp0, temp1;
   	double vtl, vti, vtt;
   	double reynolds_number;
   	double particle_term_vel;
   	double particle_fall_time;
 
  	particle_fall_time = 0.0;
  	hz = particle_ht;  /* height of the particle above sea level */
    
  	/*rho is the density of air (kg/m^3) at the elevation of the current particle*/
  	temp0 = -hz / 8200.0;
  	rho = AIR_DENSITY * exp(temp0);
  
	/*
   	(friction due to the air) :
    	vtl is terminal velocity (m/s) in laminar regime RE<6 
    	vti is terminal velocity (m/s) in intermediate regime 6<RE<500
    	vtt is terminal velocity (m/s) in turbulent regime RE>500
  	*/
  	vtl = part_density * GRAVITY * ashdiam * ashdiam / AIR_VISCOSITY_x_18; /* 18.0 * AIR_VISCOSITY */
  
  	/*
    	vti = ashdiam * 
    	pow(((4.0*GRAVITY*GRAVITY*erupt->part_mean_density *erupt->part_mean_density )/		(225.0*AIR_VISCOSITY*rho)),(1.0/3.0));
    	vtt=sqrt(3.1*erupt->part_mean_density *GRAVITY*ashdiam/rho);
  	*/
  
  	/*
    	RE is calculated using vtl (RE is Reynolds Number)
  	*/
  	reynolds_number = ashdiam * rho * vtl / AIR_VISCOSITY;
  	particle_term_vel = vtl;
  	temp0 = ashdiam * rho;


  	/*
    	c...if laminar RE>6 (intermediate regime), RE is calculated again considering vti
  	*/
  
	if (reynolds_number >= 6.0) {

    		/*4.0 * GRAVITY * GRAVITY * part_density * part_density / AIR_VISCOSITY * 225.0 * rho */
    		temp1 = GRAV_SQRD_x_4 * part_density * part_density / AIR_VISCOSITY_x_225 * rho; 
    		vti = ashdiam * pow(temp1, ONE_THIRD); /* ONE_THIRD = 1.0/3.0 */    
    		reynolds_number = temp0 * vti / AIR_VISCOSITY;
    		particle_term_vel = vti;
    		/*
    		c...if intermediate RE>500 (turbulent regime), RE is calculated again considering vtt 
  		*/ 
  		if (reynolds_number >= 500.0) {
    			vtt = sqrt( 3.1 * part_density * GRAVITY * ashdiam / rho);
    			reynolds_number =  temp0 * vtt / AIR_VISCOSITY; 
    			particle_term_vel = vtt;
  		}  
  	}
/* Calculate the time it takes this particle to fall through this distance=layer */
  particle_fall_time = layer / particle_term_vel;
  
/* particle fall time is in sec   */
  
  //printf("i= %d, layer = %f, hz = %f, particle_term_vel = %f, diam=%f, reynolds = %f\n", i,layer_thickness, hz, particle_term_vel, a//shdiam, reynolds_number);
  
  return particle_fall_time;
}



/* ----------------- New Function Starts Here -------------------- */
/* this function calculates the expected fraction of particles
   in a given grainsize class (part_size_slice) assuming a normal 
   distribution in phi units about the mean, dmean,
   with standard deviation sigma.
   The probability that 
*/

double pdf_grainsize(double part_mean_size, double part_size_slice, double part_step_width) {
  
  double func_rho, temp;
  double demon3, demon2;
  
  /* PDF_GRAINSIZE_DEMON1 = 1.0 / 2.506628 * erupt->part_sigma_size */
  demon3   = part_size_slice - part_mean_size;
  temp = -demon3 * demon3 / TWO_x_PART_SIGMA_SIZE; /* 2.0 * erupt->part_sigma_size * erupt->part_sigma_size */
  demon2   = exp(temp);
  func_rho = PDF_GRAINSIZE_DEMON1 * demon2 * part_step_width; 
  
  if (func_rho < 0.0) 
    fprintf(log_file, "error in ash size distribution - method pdf_grainsize");

  return func_rho;
 
} 
/* ----------------- New Function Starts Here -------------------- */
/* Function strat_average accounts for the variation in wind velocity 
   with height by using the average velocity value
   
		exp[ -5{ (x'-ut)^2 + y'^2} / {8*pi*C(t+td)} ]
		
   over the path of the particle as it falls from 
   its column release height to the ground.
    
   u = wind velocity (m/s), varies with height
   t = particle fall time 
   td = particle diffusion time
   
   The Suzuki equation has been formulated s.t. the wind is in the x direction
   We therefore need to transform the coordinates (xspace and yspace) to xprime
   and yprime, with xprime increasing in the downwind direction:
   x' = x cos a + y sin a
   y' = y cos a - x sin a

   
*/
double strat_average( double average_wind_direction, 
                      double average_windspeed,             
			                double xspace, double yspace, 
			                double total_fall_time,
			                double sigma) {
		
		double temp0, temp1, xprime, yprime, demon1, demon3;
			                  
    temp0 = cos(average_wind_direction);
    temp1 = sin(average_wind_direction);
    
    xprime = xspace * temp0 + yspace * temp1;
    yprime = yspace * temp0 - xspace * temp1;
    
    temp0 = xprime - average_windspeed * total_fall_time;
    demon1 = temp0 * temp0 + yprime * yprime;
	/* where sigma is calculated for the total fall time */
    demon3 = exp(-demon1/sigma);
    return demon3;
			          
}

/* 
   inputs:
   x: height of a particle within the plume, relative to vent height
   slice: integration step (index)
   ht_section_width: the width of an integration step
   none: not used
   
   output: the probability that a given grainsize will be released from a given height
*/

double plume_pdf0(double x, int slice, double none0, double none1) {

  double probability;
  static int num_slices_left = 0;
  static double plume_slice = 0.0;
  
  /* if (!slice) fprintf(stderr, "ENTER plume_pdf0 ....\n"); */
 
  probability = 0.0;
  if (x > threshold) {
     
    if (!num_slices_left) {
      num_slices_left = COL_STEPS - slice;
      plume_slice = 1.0 / (double)num_slices_left;
      /* fprintf(stderr, "slices left = %d\n ", num_slices_left); */
    }

    probability = plume_slice;
  }
  
  /* fprintf(stderr, "x=%g threshold=%g plume_slice=%g prob=%g\n", x, threshold, plume_slice, probability); 
  if (probability < 0.0) This only gets printed if an error occurs. 
    fprintf(stderr, "col_ht=%f prob=%f\n", x, probability); */

  return probability; 
}

/* 
   inputs:
   x: height of a particle within the plume, relative to vent height
   slice: integration step
   beta: the column beta parameter
   none: not used
   
   output: the probability that a given grainsize will be released from a given height
*/

double plume_pdf1(double x, int slice, double plume, double total) {

  double col_ht, probability, beta_limit;
  double temp1, temp0, demon1, demon2;

  beta_limit = plume;
  col_ht = beta_limit - beta_limit * x / total;
  if (col_ht <= 0.0) col_ht = 1e-9; 
  temp1 = log(col_ht);
  temp1 *= temp1;
  temp0 = -temp1/TWO_BETA_SQRD; /* 2.0 * beta * beta */
  demon1 = exp(temp0);
  demon2 = col_ht * BETA_x_SQRT_TWO_PI; /* beta * sqrt(2.0 * PI) */
  

  probability = demon1 / demon2;
  if (isnan(probability)) {
    fprintf(stderr, "ht=%g  demon1=%g demon2=%g temp0=%g\n", 
	    col_ht, demon1, demon2, temp0); 
    exit(-1);
  }
  
  if (probability < 0.0) {/* This only gets printed if an error occurs. */
    fprintf(stderr, "col_ht=%f demon1=%f demon2=%f prob=%f\n", col_ht, demon1, demon2, probability); 
		exit(-1);
	}
  return probability;
}

void set_global_values(FILE *log) {

  log_file = log;
#ifdef _PRINT  
  fprintf(log_file, "IN set_global_values ...");
#endif
  /* Set values for global static variables */
  
  AIR_VISCOSITY_x_18 = 18.0 * AIR_VISCOSITY;
  LITHIC_DENSITY_minus_PUMICE_DENSITY = LITHIC_DENSITY - PUMICE_DENSITY;
  PUMICE_DIAMETER_THRESHOLD_minus_LITHIC_DIAMETER_THRESHOLD = PUMICE_DIAMETER_THRESHOLD - LITHIC_DIAMETER_THRESHOLD;
  ONE_THIRD = 1.0 / 3.0;
  AIR_VISCOSITY_x_225 = AIR_VISCOSITY * 225.0;
  GRAV_SQRD_x_4 = 4.0 * GRAVITY * GRAVITY;
  EDDY_CONST_x_8_div_5 = 8.0 * EDDY_CONST / 5.0;
  T = NULL;
#ifdef _PRINT  
  fprintf(log_file, "OUT");
#endif
}

void set_eruption_values(ERUPTION *erupt, WIND *wind) { /* set_eruption_values */

  /* The following parameters are the properties of a eruption
   * each eruption must have all of these parameters defined:
   *
   * erupt->total_ash_mass is the total amount of ash erupted by
   * the volcano over the course of the entire eruption or calculation period
   * erupt->max_part_size is the maximum particle diameter considered
   * in the calculation. This is input in phi units (so it will likely be
   * a negative number like -5 and appear to be less than min_part_size)
   * erupt->min_part_size is the minimum particle diameter condsidered in the
   * calculation. This input is in phi units.
   *
   * Note: erupt->max/min_part_size are used to set the limits of integration
   * on the calculation. Particles outside this range are not considered at all.
   *
   * erupt->part_mean_size is the mean particle diameter erupted in phi units
   * erupt->part_sigma_size is the standard deviation in particle diameter in phi units
   * erupt-> vent_height is the elevation of the vent m.a.s.l. in meters
   * erupt->max_column_height is the eruption column height m.a.s.l. 
   * (not used) erupt->column_beta is the shape factor governing 
   * the particle size distribution in the eruption column. 
   */

  int i, j;

  double y, prob;
  double x, total_P_col, total_P_part, cum_prob_part, cum_prob_col, total_P;
  double particle_ht, cum_fall_time, wind_x, wind_y, ht_above_vent, temp;
  double col_prob, part_prob;

  double ht_section_width;
  double part_section_width;
  double ht_step_width;
  double part_step_width;
  double pmin=10e6, pmax=0.0;
  
#ifdef _PRINT
  fprintf(log_file, "IN set_eruption_values ... ");
#endif

  PART_STEPS = (erupt->max_phi - erupt->min_phi) * 10;
#ifdef _PRINT
  fprintf(log_file, "PART_STEPS=%d\n", PART_STEPS);
#endif

  /* threshold = (PLUME_RATIO * (erupt->max_plume_height - erupt->vent_height)); replaced 01-23-2011 */
  threshold = erupt->vent_height + (PLUME_RATIO * (erupt->max_plume_height - erupt->vent_height)); /*new line 01-23-2011 */
  SQRT_TWO_PI = sqrt(2.0 * pi); /* new line, 12-2-2010 */
  BETA_x_SQRT_TWO_PI = erupt->column_beta * SQRT_TWO_PI;
  TWO_BETA_SQRD = 2.0 * erupt->column_beta * erupt->column_beta;
  /* PDF_GRAINSIZE_DEMON1 = 1.0 / 2.506628 * erupt->sigma_phi; line removed, 12-2-2010 */
  PDF_GRAINSIZE_DEMON1 = 1.0 / (2.506628 * erupt->sigma_phi); /* changed 12-02-2010 */
  TWO_x_PART_SIGMA_SIZE = 2.0 * erupt->sigma_phi * erupt->sigma_phi;
  
  /*define the limits of integration */ 
  ht_section_width = erupt->max_plume_height - erupt->vent_height; 
  part_section_width = erupt->max_phi - erupt->min_phi;
  ht_step_width = ht_section_width / (double)COL_STEPS; 
  part_step_width = part_section_width / (double)PART_STEPS;

  /* steps for nomalization of probabilities */
  cum_prob_col = 0.0;
  x = erupt->vent_height;
  for (i=0; i < COL_STEPS; i++) {
    x += ht_step_width;
     
    prob = (*pdf)(x, (double)i, ht_section_width, erupt->max_plume_height);
    cum_prob_col += prob;
     
    //#ifdef _PRINT
    //fprintf(stderr, " slice_ht=%g, prob=%g, cum_prob=%g\n", x, prob, cum_prob_col); 
    //#endif
  }    
  total_P_col = cum_prob_col; 
  //fprintf( stderr, "total_P_col=%g\n ", total_P_col);

  cum_prob_part = 0.0;
  y = (erupt)->min_phi;
  for (i=0; i < PART_STEPS; i++) {
    
    prob = pdf_grainsize(erupt->mean_phi, y, part_step_width);
    cum_prob_part += prob;
    //fprintf(stderr, " grain_size=%.2f, prob=%g, cum_prob=%g\n", y, prob, cum_prob_part);
    
#ifdef _PRINT
    fprintf(log_file, " grain_size=%g, prob=%g, cum_prob=%g\n", y, prob, cum_prob_part);
#endif
    y += part_step_width;
  }
  total_P_part = cum_prob_part;
  //fprintf( stderr, "total_P_part=%g \n", total_P_part);

  /* Normalization constant */
  total_P = (total_P_col * total_P_part);


  /* End of normalization steps */  

  /* Dynamically allocated table for storing integration data.
     Used in the double integration steps below for each point considered.
  */
    if (T == NULL) {

      T = (TABLE **)GC_MALLOC((size_t)PART_STEPS * sizeof(TABLE *));
      if (T == NULL) {
        fprintf(log_file, 
        "Cannot malloc memory for Integration Table:[%s]\n", strerror(errno));
        exit(1);
      }
      for (i=0; i<PART_STEPS; i++) {
        T[i] = (TABLE *)GC_MALLOC((size_t)COL_STEPS * sizeof(TABLE));
        if (T[i] == NULL) {
          fprintf(log_file, 
          "Cannot malloc memory for Integration Table[%d]:[%s]\n", i, strerror(errno));
          exit(1);
        }
      }
    } else {
      T = (TABLE **)GC_REALLOC(T, (size_t)PART_STEPS * sizeof(TABLE *));
      if (T == NULL) {
        fprintf(log_file, 
        "Cannot malloc memory for Integration Table:[%s]\n", strerror(errno));
        exit(1);
      }
      for (i=0; i<PART_STEPS; i++) {
        T[i] = (TABLE *)GC_REALLOC(T[i], (size_t)COL_STEPS * sizeof(TABLE));
        if (T[i] == NULL) {
          fprintf(log_file, 
          "Cannot malloc memory for Integration Table[%d]:[%s]\n", i, strerror(errno));
          exit(1);
        }
      }
    }
/* fprintf(log_file,
	      "\nPart_Ht\tAsh_Diam\tPart-Den\tFalltime\tDFP\tDCP\tTFalltime\tWsumX\tWsumY\tdemon1\n");
*/	      
    /* Start with the maximum particle size */
    y = (erupt)->min_phi;
    for (i = 0; i < PART_STEPS; i++) { /* PART_STEPS_LOOP */
      /* y += part_step_width;  changed to end of loop 12-2-2010 */
      T[i][0].part_density  =  particle_density(y);    
      T[i][0].ashdiam = phi2m(y);
      
      /* the expected fraction of particles of this size based on given mean and std deviation */
      part_prob = pdf_grainsize(erupt->mean_phi, y, part_step_width);
      cum_fall_time = 0.0;
      wind_x = 0.0;
      wind_y = 0.0;
      
      /* Start at the height of the vent */
      particle_ht = erupt->vent_height;
      for (j = 0; j < COL_STEPS; j++) { /* COL_STEPS_LOOP */

	      /* define the small slice dz */
	      particle_ht += ht_step_width;
	      
	      /* Calculate the time it takes a particle to fall from its release point
	         in the column to the next column release point.
	       */    
	      T[i][j].fall_time = 
	      part_fall_time(particle_ht, ht_step_width, T[i][0].ashdiam, T[i][0].part_density); 
	      
    	 /* Particle diffusion time (seconds) */
    	 ht_above_vent = particle_ht - erupt->vent_height;
	     
	     temp = 0.2 * ht_above_vent * ht_above_vent; 
       T[i][j].plume_diffusion_fine_particle = 
       pow(temp, 0.4); /* 0.4 = 2.0/5.0 */
       
	     T[i][j].plume_diffusion_coarse_particle = 
	     0.0032 * (ht_above_vent *  ht_above_vent) / DIFFUSION_COEFFICIENT; 
	     
	      /* Sum the windspeed and wind_direction for each particle size 
	       * falling from each level. In the wind array, the first wind level
	       * gives wind speed and direction at the vent height. 
	       * Start with the next wind level, 
	       * so that we are using the wind speed and direction 
	       * starting from one step above the vent. 
	       */
     
	     wind_x += 
	     T[i][j].fall_time * wind[j+1].windspeed * cos(wind[j+1].wind_dir);
	     
	     wind_y += 
	     T[i][j].fall_time * wind[j+1].windspeed * sin(wind[j+1].wind_dir);
	      
	     T[i][j].wind_sum_x = wind_x;
       T[i][j].wind_sum_y = wind_y;
	      
	      /* Accumulate the time it takes each particle size to descend
	         from its release point down
	         to its final resting place.This part of the code just 
	         calculates the fall_time from the release point to the 
	         height of the vent.
	         The time it takes a particle to fall from the vent height 
	         to a grid cell will be calculated later. 
	       */
	      cum_fall_time += T[i][j].fall_time;
	      T[i][j].total_fall_time = cum_fall_time;
	      if (T[i][j].total_fall_time > pmax) pmax = T[i][j].total_fall_time;
			  if (T[i][j].total_fall_time < pmin) pmin = T[i][j].total_fall_time;
	      /* the probability that a given grainsize will be released from a given height */
        col_prob = 
        (*pdf)(particle_ht, j, erupt->column_beta, erupt->max_plume_height);
      
        /* Normalization is now done here */
        T[i][j].demon1 =  
        (erupt->total_ash_mass * col_prob  * part_prob)/total_P;
        
	      T[i][j].particle_ht = particle_ht;
	/*      	
	      fprintf(log_file,
	      "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
	      T[i][j].particle_ht, 
	      T[i][j].ashdiam, 
	      T[i][j].part_density, 
	      T[i][j].fall_time,
	      T[i][j].plume_diffusion_fine_particle,
	      T[i][j].plume_diffusion_coarse_particle,
	      T[i][j].total_fall_time,
	      T[i][j].wind_sum_x,
	      T[i][j].wind_sum_y,
	      T[i][j].demon1);
*/	      
      } /* END COL_STEPS_LOOP */ 
      
 /*     fprintf(log_file, "\n"); */
    	y += part_step_width; /* moved from beg of loop 12-02-2010 */
         
    } /* END PART_STEPS_LOOP */
/*  	fprintf(log_file, "OUT\n"); */
//  fprintf(stderr, "MIN particle fall time = %.1f\n", pmin);
//  fprintf(stderr, "MAX particle fall time = %.1f\n", pmax);
}
