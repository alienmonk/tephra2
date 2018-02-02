#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <gc.h>
#include "prototypes_strat.h"

/*
Code: new_tephra.c
By: C.B. & L.J. Connor, T. Hincks, and C. Bonadonna
Copyright (C) 2003  C.B. Connor, L.J. Connor, C. Bonadonna, T. Hincks
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

/* The following Global Variables are assigned some default values */

double DIFFUSION_COEFFICIENT = 200.0;
double FALL_TIME_THRESHOLD = 180.0;
double EDDY_CONST = .04;
double LITHIC_DENSITY = 2350.0;
double PUMICE_DENSITY = 1000.0;
int PART_STEPS = 100;
int COL_STEPS = 100;
int PLUME_MODEL = DISTR_1;
double PLUME_RATIO = 0.1;
double WIND_INTERVAL;
int WIND_DAYS = 1;
int WIND_COLUMNS = 3;
double PLUME_HEIGHT = 5.0;
double ERUPTION_MASS = 1e10;
double MAX_GRAINSIZE = -7.0;
double MIN_GRAINSIZE = 7.0;
double MEDIAN_GRAINSIZE = 1.5;
double STD_GRAINSIZE = 2.0;
double VENT_EASTING = 0.0;
double VENT_NORTHING = 0.0;
double VENT_ELEVATION = 0.0;

/*define the following data structures for the code
  full descriptions are found in common_structures.h */
static ERUPTION *erupt;
static WIND **W;
static POINT *pt;

/*define the following variables passed among functions in this file */
static int num_pts = 0; /*total number of points used in the analysis */
static int num_eruptions = 0; /*total number of eruptions used in the analysis */
/* static int num_wind_data = 0; number of wind data points in the analysis */

static int local_n;

//FILE *in_eruptions;
FILE *in_points;
FILE *in_wind;
FILE *log_file;


int cmp(double *x, double *y) {
  if (*x < *y) return -1;
  if (*x > *y) return 1;
  return 0;
}

void exit_now(int e) {
  
  //(void) fclose(in_eruptions);
  (void) fclose(in_points);
  (void) fclose(in_wind);
#ifdef _PRINT
  (void) fclose(log_file);
#endif
#ifdef DEBUG
  (void) fclose(log_file);
#endif
  exit(e);
}

int main(int argc, char *argv[]) { /* MAIN CODE STARTS HERE */
//FILE *Output;
//Output=fopen("Tephra2out.dat","w");
  
  int i, j, bin, phi_bins;
  double val;
  STATS stats;
  char log_name[25];


  /* Check for correct number of comand line arguments */
  if (argc != 4) {
      fprintf(stderr, 
	      "Missing comand line arguments,\nUSAGE: <program name> <config file> <points file> <wind file>\n\n");
    exit(1);
  }
  

  /* Each node opens a file for logging */
  sprintf(log_name, "%s", LOG_FILE);
  
//  fprintf(stderr, "%s\n", log_name);
  
  log_file  = fopen(log_name, "w+");
  if (log_file == NULL) {
    fprintf(stderr, "Cannot open LOG file=[%s]:[%s]. Exiting.\n", 
	    log_name, strerror(errno));
    exit(1);
  }

  
  /* Initialize the global variables (see top of file) with inputs from the configuration file. */
  if ( init_globals(argv[1]) ) {
    exit(1);
  }
  
  /*make sure the eruptions file exists 
  in_eruptions = fopen(argv[2], "r");
  if (in_eruptions == NULL) {

    fprintf(stderr, "Cannot open eruptions  file=[%s]:[%s]. Exiting.\n", 
	    argv[2], strerror(errno));
    exit_now(1);
  }
  */
  
#ifdef _PRINT
  fflush(log_file); 
#endif

  /*make sure the points file exists*/
  in_points= fopen(argv[2], "r");
  if (in_points == NULL) {
    fprintf(stderr, "Cannot open points  file=[%s]:[%s]. Exiting.\n", 
	    argv[2], strerror(errno));
    exit_now(1);
  }
  
  /* Input the data points from a file using the
     get_points function. 
  */
  
  if (get_points(in_points) ) {
    exit_now(1);
  }

#ifdef _PRINT
  fflush(log_file); 
#endif
  
  /*make sure the wind file exists*/
  in_wind= fopen(argv[3], "r");
  if (in_wind == NULL) {
    fprintf(stderr, "Cannot open wind file=[%s]:[%s]. Exiting.\n", 
	    argv[3], strerror(errno));
    exit_now(1);
  }
  
  if (get_wind(in_wind) ) {
    exit_now(1);
  }
  
  /* 
     
  Note: "local_n" is the number of points from the
  input file. Therefore local_n is
  declared a static variable and is assigned a value in the
  get_points function 
  
  Note: tephra_calc is the main computational part of the
  code.
  
  Input the data points from wind file using the get_wind function. 
  
 
  */
  
#ifdef _PRINT
  fflush(log_file); 
#endif
  
  set_global_values(log_file);
  
#ifdef _PRINT
  fflush(log_file); 
#endif
 
/* each node gets all of the eruption data */ 
  if (get_eruptions() ) {
    exit_now(1);
  }

/* Calculating an accumulation map */
    
  for (j=0;  j < num_eruptions; j++) { /* For each eruptive senario */
    set_eruption_values(erupt+j, W[j]);
    //fprintf(stderr, "[%d]PARTICLE STEPS=%d ", j, PART_STEPS);
      for (i = 0;i < local_n; i++) {  /* For each location */
	      /* Note: W[j]  means that if there are multiple eruptions, 
	      there should be multiple WIND_DAYS in the wind.in file, 
	      1 WIND_DAY for each eruption line */
	      tephra_calc(erupt+j, pt+i, W[j], &stats); 
	      (pt+i)->cum_mass += (pt+i)->mass;
      }
    } j--;
 //   fprintf(stderr,"Min Particle Fall Time = %.2g(s)\nMax Particle Fall Time = %.2g(s)\n\n",stats.min_falltime, stats.max_falltime);
    
	/* Here I am assuming that each eruptive senario has the same min and max particle size range
	   So, the grainsize distribution contains the same number of bins. The amount in each bin accumulates
	   between eruptive senarios.
	*/
	phi_bins = (int)((erupt+j)->max_phi - (erupt+j)->min_phi);
	fprintf(stdout, "#Easting Northing Elev. Mass/Area  "); 
	for (bin = 0; bin < phi_bins; bin++) 
	  fprintf(stdout,"[%d->%d) ", 
		(int)((erupt+j)->min_phi)+bin, 
		(int)((erupt+j)->min_phi)+bin+1);
	
	fprintf(stdout,"\n");
	
//	fprintf(stderr, "\nPART_STEPS=%d phi_bins=%d\n", PART_STEPS, phi_bins);
	for (i=0; i < num_pts; i++) {
	  fprintf(stdout, "%.0f %.0f %.0f %.2g  ", 
		(pt+i)->easting, 
		(pt+i)->northing,
		(pt+i)->elevation, 
		(pt+i)->cum_mass);
	  for (bin=0; bin < phi_bins; bin++) {
	  	val = ((pt+i)->phi[bin]/(pt+i)->cum_mass) * 100.0;
	    fprintf(stdout, "%.2g ", val);
	  }
	  fprintf(stdout, "\n");
	}



  fprintf(log_file, "Finished.\n");
 // (void) fclose(Output);
  exit_now(0);
  return 1;
}

/**************************************************************
FUNCTION:  get_eruptions
DESCRIPTION:  This function reads eruption data into the
ERUPTION array. Each node stores 
all of the eruption  parameters which are varied and then 
used in calculating the mass loading value at each point. 
INPUTS: (IN) FILE *in  (file handle from which to read)
OUTPUTS: int -1=error, 0=no error
***************************************************************/
int get_eruptions(void) {
  
  //char line[MAX_LINE];
  int i;
  
#ifdef _PRINT
  fprintf(log_file,"ENTER[get_eruptions]\n");
#endif
  
  /* while (fgets(line, MAX_LINE, in) != NULL) {
    if (line[0] == '#' || line[0] == '\n') continue;
    num_eruptions++;
  } */
  num_eruptions = 1;
  erupt = (ERUPTION *)GC_MALLOC((size_t)(num_eruptions +1) * sizeof(ERUPTION));
  
  if (erupt == NULL) {
    fprintf(stderr, "[%d-of-%d]\tCannot malloc memory for eruptions:[%s]\n",
            0,0, strerror(errno));
    return -1;
  } 
  
#ifdef _PRINT 
  fprintf(log_file,"\t%d eruptions in file.\n", num_eruptions);
#endif  
  
  //rewind(in);

  i=0;
 /* while (fgets(line, MAX_LINE, in) != NULL) {
    if (line[0] == '#' || line[0] == '\n') continue;
    else {
      while (ret = sscanf(line,
			  "%lf %lf %lf %lf %lf %lf %lf %lf %lf",*/
			  (erupt+i)->volcano_easting = VENT_EASTING;
			  (erupt+i)->volcano_northing = VENT_NORTHING;
			  (erupt+i)->total_ash_mass = ERUPTION_MASS;
			  (erupt+i)->min_phi = MAX_GRAINSIZE;
			  (erupt+i)->max_phi= MIN_GRAINSIZE;
			  (erupt+i)->mean_phi = MEDIAN_GRAINSIZE;
			  (erupt+i)->sigma_phi= STD_GRAINSIZE;
			  (erupt+i)->vent_height = VENT_ELEVATION;	
			  (erupt+i)->max_plume_height = PLUME_HEIGHT;
//), ret != 9) { 
// &(erupt+i)->column_beta,
/*	
	if (ret == EOF && errno == EINTR) continue;
	fprintf(stderr, "[line=%d,ret=%d] Did not read in 9 parameters:[%s]\n", i+1,ret, strerror(errno));
	return -1;
      }
      i++;
    }
  }
 */ 
#ifdef _PRINT
  fprintf(log_file, "EXIT[get_eruptions].\n");	  
#endif
  return 0;
}


/**************************************************************
FUNCTION:  get_points
DESCRIPTION:  This function reads eruption data into the
ERUPTION array. Each node stores 
all of the eruption  parameters which are varied and then 
used in calculating the mass loading value at each point. 
INPUTS: (IN) FILE *in  (file handle from which to read)
OUTPUTS: int -1=error, 0=no error
***************************************************************/
int get_points(FILE *in) {
  
  char line[MAX_LINE];
  int i, j, ret, my_start=0, pts_read=0;
  
#ifdef _PRINT
  fprintf(log_file,"ENTER[get_points]\n");
#endif
  while (fgets(line, MAX_LINE, in) != NULL)  {
    if (line[0] == '#' || line[0] == '\n') continue;
    num_pts++;
  }
  rewind(in);
  

  local_n = num_pts;
  

#ifdef _PRINT  
    fprintf(log_file, "Total locations: %d.\n", num_pts);
#endif    
  
  pt = (POINT *)GC_MALLOC((size_t)local_n * sizeof(POINT));
  if (pt == NULL) {
    fprintf(stderr, "Cannot malloc memory for my points:[%s]\n", strerror(errno));
    return -1;
  } 
  
#ifdef _PRINT  
  fprintf(log_file, "\tReading in %d locations, starting at line %d.\n",
	  local_n, my_start);
#endif    
   
  i=0;
  while (i < num_pts) {
    fgets(line, MAX_LINE, in);
    if (line[0] == '#' || line[0] == '\n') continue;
    else {
      while (ret = sscanf(line,
			  "%lf %lf %lf",
			  &(pt+pts_read)->easting,
			  &(pt+pts_read)->northing,
			  &(pt+pts_read)->elevation),
	     ret != 3) { 
	
	if (ret == EOF && errno == EINTR) continue;
	fprintf(stderr, "[function:get_points, line=%d,ret=%d] Did not read in 3 coordinates\n", i+1,ret);
	return -1;
      }
      /* Initialize some values in the point structure */
      (pt+pts_read)->cum_mass = 0.0;
      (pt+pts_read)->mass_pt = NULL;
      for (j=0; j<20; j++)
	(pt+pts_read)->phi[j] = 0.0;
      if (i >= my_start) {
	pts_read++;
	if (pts_read == local_n) break;
      } 
    }
    i++;
  }
  
  
#ifdef _PRINT  
  fprintf(log_file, "EXIT[get_points].\n");
#endif		  
  return 0;
}

/**************************************************************
FUNCTION:  get_wind
DESCRIPTION:  This function reads wind data into the
WIND array. Each node stores all of the wind data. 

INPUTS: (IN) FILE *in  (file handle from which to read)
OUTPUTS: int -1=error, 0=no error
***************************************************************/
int get_wind(FILE *in) {
  
  int i=0, j=0, ret;
  char line[MAX_LINE];
  double wind_height, wind_dir, windspeed, dir0, ht0, sp0;
  double level;
  
#ifdef _PRINT
  fprintf(log_file,"ENTER[get_wind].\n");
#endif
  
  WIND_INTERVAL = (PLUME_HEIGHT - VENT_ELEVATION)/COL_STEPS;
  
  W = (WIND**)GC_MALLOC(WIND_DAYS * sizeof(WIND *));
  if (W == NULL) {
    fprintf(stderr, "Cannot malloc memory for wind columns:[%s]\n", strerror(errno));
    return -1;
  } else {
    for (i=0; i < WIND_DAYS; i++) {
      W[i] = (WIND *)GC_MALLOC((COL_STEPS+1) * sizeof(WIND));
      if (W[i] == NULL) {
	fprintf(stderr, "Cannot malloc memory for wind rows %d:[%s]\n", i, strerror(errno));
	return -1;
      }
    }
  }
  /* Assume one wind day */
  i=0;
  /* start at the vent */ 
    level = VENT_ELEVATION;
    
  /* Do for each column step */
  /* j = 0 is for the interval between the vent and the ground.
   * Here we set the wind speed and direction to be at the level of the vent;
   * The values used in the calculations change for each location and are 
   * set in the tephra_calc routine when the point elevation is known. 
   * The last interval ends at the top of the column. 
   */
  for (j=0; j <= COL_STEPS; j++) { 
    W[i][j].wind_height = 0.0;
    ht0 = 0.0;
    dir0 = 0.0;
    sp0 = 0.0;
    
    /* Find wind elevation just greater than current level */
    /* Start scanning the wind file for the best match.
     * Each new level starts scanning the file from the beginning.
     */
    while (NULL != fgets(line, MAX_LINE, in)) {
	    if (line[0] == '#' || strlen(line) < WIND_COLUMNS) continue;
	    else {
	      while (ret = sscanf(line,
			      "%lf %lf %lf",
			      &wind_height,
			      &windspeed,
			      &wind_dir), ret != 3) { 
	    
	        if (ret == EOF && errno == EINTR) continue;
	        
	        fprintf(stderr, 
	        "[line=%d,ret=%d] Did not read in 3 parameters:[%s]\n", 
	        i+1,ret, strerror(errno));
	        
	        return -1;
	      }
	    }
	    
	    /* This is the case where we find the first height that is equal to
	     * or greater that the level that we are assigning.
	     */
	    if (wind_height >= level) {
	      if(wind_height == level) {
	        W[i][j].wind_dir = wind_dir;
	        W[i][j].windspeed = windspeed;
	        
	      } else { /* interpolate */
	        W[i][j].wind_dir = 
	        ((wind_dir - dir0) * (level - ht0) / (wind_height - ht0)) + dir0;
	      
	        W[i][j].windspeed = 
	        ((windspeed - sp0) * (level - ht0) / (wind_height - ht0)) + sp0;
	      }
	      W[i][j].wind_height = level;
	      fprintf(log_file, 
	      "%f %f %f\n", 
	      W[i][j].wind_height, W[i][j].windspeed, W[i][j].wind_dir);
	      W[i][j].wind_dir *= DEG2RAD; /* change to radians */
	      break; /* ready to rescan the file for a match for the next level */
	    }
	    /* This is the case where the scanned height is less than the level
	     * we are assigning.
	     */
	    else {
	      /* Maintain the scanned values for possible interpolation 
	       * at the next level.
	       */
	      ht0 = wind_height;
	      dir0 = wind_dir;
	      sp0 = windspeed;
	    }
	  }
	  /* If we finish scanning the file and all heights are below the level we are
	   * currently assigning, then just use the direction and speed
	   * at the upper-most height.
	   */
	  if (!W[i][j].wind_height) {
	    W[i][j].wind_height = level;
	    W[i][j].windspeed = sp0;
	    W[i][j].wind_dir = dir0;
	  }
	  /* Go to the next column height */
	  rewind(in); 
	  level += WIND_INTERVAL; 
  } 

  
#ifdef _PRINT
  fprintf(log_file, "\tRead %d wind days with %d wind levels per day.\n", i, j);
  fprintf(log_file, "EXIT[get_wind].\n");
#endif	  	  
  return 0;
}

/**************************************************************
FUNCTION:  get_config_data
DESCRIPTION:  This function reads the configuration file,
and sets some global variables.

INPUTS: (IN) FILE *in  (file handle from which to read)
OUTPUTS: int -1=error, 0=no error
***************************************************************/

int init_globals(char *config_file) {

  FILE *in_config;
  char buf[1][30], **ptr1;
  char line[MAX_LINE];
  char space[4] = "\n\t ";
  char *token;
  
#ifdef _PRINT  
  fprintf(log_file, "ENTER[init_globals].\n");
#endif
 
  in_config = fopen(config_file, "r");
  if (in_config == NULL) {
    fprintf(stderr, 
	    "Cannot open configuration file=[%s]:[%s]. Exiting.\n", config_file, strerror(errno)); 
    return 1;
  }
  
  ptr1 = (char **)&buf[0];
  while (fgets(line, MAX_LINE, in_config) != NULL) {
    /* fprintf(stderr, "%s\n", line); */
    if (line[0] == '#' || line[0] == '\n' || line[0] == ' ') {
      /* fprintf(stderr, "%s\n", line); */
      continue;
    }
    token = strtok_r(line, space, ptr1);
    if (!strncmp(token, "DIFFUSION_COEFFICIENT", strlen("DIFFUSION_COEFFICIENT"))) {
      token = strtok_r(NULL,space,ptr1);
      DIFFUSION_COEFFICIENT = strtod(token, NULL);
      /* DIFFUSION_COEFFICIENT can never be 0 as it is used in divisions */
      if (!DIFFUSION_COEFFICIENT) DIFFUSION_COEFFICIENT = 1.0;
//      fprintf(stderr, "DIFFUSION_COEFFICIENT=%.1f\n", DIFFUSION_COEFFICIENT);

    }else if (!strncmp(token, "EDDY_CONST", strlen("EDDY_CONST"))) {
      token = strtok_r(NULL,space,ptr1);
      EDDY_CONST = strtod(token, NULL);
//      fprintf(stderr, "EDDY_CONST=%g\n", EDDY_CONST);
    }
    else if (!strncmp(token, "FALL_TIME_THRESHOLD", strlen("FALL_TIME_THRESHOLD"))) {
      token = strtok_r(NULL,space,ptr1);
      FALL_TIME_THRESHOLD = strtod(token, NULL);
//      fprintf(stderr, "FALL_TIME_THRESHOLD=%.1f\n", FALL_TIME_THRESHOLD);
    }
    else if (!strncmp(token, "LITHIC_DENSITY", strlen("LITHIC_DENSITY"))) {
      token = strtok_r(NULL,space,ptr1);
      LITHIC_DENSITY = strtod(token, NULL);
//      fprintf(stderr, "LITHIC_DENSITY=%.1f\n", LITHIC_DENSITY);
    }
    else if (!strncmp(token, "PUMICE_DENSITY", strlen("PUMICE_DENSITY"))) {
      token = strtok_r(NULL,space,ptr1);
      PUMICE_DENSITY = strtod(token, NULL);
//      fprintf(stderr, "PUMICE_DENSITY=%.1f\n", PUMICE_DENSITY);
    }
    else if (!strncmp(token, "PART_STEPS", strlen("PART_STEPS"))) {
      token = strtok_r(NULL, space, ptr1);
      PART_STEPS = (int)atoi(token);
//      fprintf(stderr, "PART_STEPS = %d\n", PART_STEPS);
    }
    else if (!strncmp(token, "COL_STEPS", strlen("COL_STEPS"))) {
      token = strtok_r(NULL, space, ptr1);
      COL_STEPS = (int)atoi(token);
//      fprintf(stderr, "COL_STEPS = %d\n", COL_STEPS);
    }
    else if (!strncmp(token, "PLUME_MODEL", strlen("PLUME_MODEL"))) {
      token = strtok_r(NULL,space,ptr1);
      PLUME_MODEL = (int)atoi(token);
      
      if (PLUME_MODEL == 0) {
	      pdf = plume_pdf0;
//	      fprintf(stderr, "PLUME_MODEL=[%d]%s\n", PLUME_MODEL, "Uniform Distribution with threshold");
      }
      else if (PLUME_MODEL == 1) {
	      pdf = plume_pdf1;
//	      fprintf(stderr, "PLUME_MODEL=[%d]%s\n", PLUME_MODEL, "log-normal Distribution using beta");
      }
    }
    else if (!strncmp(token, "PLUME_RATIO", strlen("PLUME_RATIO"))) {
      token = strtok_r(NULL, space, ptr1);
      PLUME_RATIO = strtod(token, NULL);
//      if (!PLUME_MODEL) fprintf(stderr, "PLUME_RATIO = %.2f\n", PLUME_RATIO);
    }
    else if (!strncmp(token, "WIND_DAYS", strlen("WIND_DAYS"))) {
      token = strtok_r(NULL, space, ptr1);
      WIND_DAYS = (int)atoi(token);
//      fprintf(stderr, "WIND_DAYS = %d\n", WIND_DAYS);
    }
    else if (!strncmp(token, "WIND_COLUMNS", strlen("WIND_COLUMNS"))) {
      token = strtok_r(NULL, space, ptr1);
      WIND_COLUMNS = (int)atoi(token);
//      fprintf(stderr, "WIND_COLUMNS = %d\n", WIND_COLUMNS);
    }
    else if (!strncmp(token, "PLUME_HEIGHT", strlen("PLUME_HEIGHT"))) {
      token = strtok_r(NULL, space, ptr1);
      PLUME_HEIGHT = strtod(token, NULL);
//      fprintf(stderr, "PLUME_HEIGHT = %.1f\n", PLUME_HEIGHT);
    }
    else if (!strncmp(token, "ERUPTION_MASS", strlen("ERUPTION_MASS"))) {
      token = strtok_r(NULL, space, ptr1);
      ERUPTION_MASS = strtod(token, NULL);
//      fprintf(stderr, "ERUPTION_MASS = %g\n", ERUPTION_MASS);
    }
    else if (!strncmp(token, "MAX_GRAINSIZE", strlen("MAX_GRAINSIZE"))) {
      token = strtok_r(NULL, space, ptr1);
      MAX_GRAINSIZE = strtod(token, NULL);
//      fprintf(stderr, "MAX_GRAINSIZE = %.0f\n", MAX_GRAINSIZE);
    }
    else if (!strncmp(token, "MIN_GRAINSIZE", strlen("MIN_GRAINSIZE"))) {
      token = strtok_r(NULL, space, ptr1);
      MIN_GRAINSIZE = strtod(token, NULL);
//      fprintf(stderr, "MIN_GRAINSIZE = %.0f\n", MIN_GRAINSIZE);
    }
    else if (!strncmp(token, "MEDIAN_GRAINSIZE", strlen("MEDIAN_GRAINSIZE"))) {
      token = strtok_r(NULL, space, ptr1);
      MEDIAN_GRAINSIZE = strtod(token, NULL);
//      fprintf(stderr, "MEDIAN_GRAINSIZE = %.2f\n", MEDIAN_GRAINSIZE);
    }
    else if (!strncmp(token, "STD_GRAINSIZE", strlen("STD_GRAINSIZE"))) {
      token = strtok_r(NULL, space, ptr1);
      STD_GRAINSIZE = strtod(token, NULL);
//      fprintf(stderr, "STD_GRAINSIZE = %.2f\n", STD_GRAINSIZE);
    }
    else if (!strncmp(token, "VENT_EASTING", strlen("VENT_EASTING"))) {
      token = strtok_r(NULL, space, ptr1);
      VENT_EASTING = strtod(token, NULL);
//      fprintf(stderr, "VENT_EASTING = %.1f\n", VENT_EASTING);
    }
    else if (!strncmp(token, "VENT_NORTHING", strlen("VENT_NORTHING"))) {
      token = strtok_r(NULL, space, ptr1);
      VENT_NORTHING = strtod(token, NULL);
//      fprintf(stderr, "VENT_NORTHING = %.1f\n", VENT_NORTHING);
    }
    else if (!strncmp(token, "VENT_ELEVATION", strlen("VENT_ELEVATION"))) {
      token = strtok_r(NULL, space, ptr1);
      VENT_ELEVATION = strtod(token, NULL);
//      fprintf(stderr, "VENT_ELEVATION = %.1f\n", VENT_ELEVATION);
    }
    else continue;
  }
  (void) fclose(in_config);

#ifdef _PRINT
  fprintf(log_file, "EXIT[init_globals].\n");
#endif


  return 0;
}
