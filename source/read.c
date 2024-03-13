#include "global.h"
#include "prototypes.h"

/*..............................................................................*/

/*
Reads the run parameters from the external file "options".  See INSTRUCTIONS
for a list of keywords.
*/

void read_options(long *npart, struct vector *box, long *nsweeps,
   long *report, int *periodic,
   long *equilibrate, char *configurations_filename,
   long *sourcetype, char *potential_filename, long *num_test_particles, double *kt,
   int *grid)
{
   int i;
   char command[20];
   char option[20];
   char error[200];
   char fname[80]; //name of file where configurations are read from
   long tmp_grid; // Temporary variable for loading grid value

   FILE *infile = NULL;
   FILE *datafile = NULL; //LAMMPS file


   /* Prototypes for keywords library */
   int  read_line(FILE *);
   int  get_string(char [], int);
   void upper_case(char []);
   int  get_int(long int *);
   int  get_double(double *);


   /*--- 0. Defaults ---*/
   box->x = box->y -1.0;              /* Box dimensions */
   *npart = 10;                       /* Number of rods */
   *nsweeps = 100;                    /* Number of accumulation sweeps */
   *periodic = 1;                     /* Wrapping rather than spanning percolation criterion */
   *report = 100;                     /* Number of sweeps between statistics reports */
   seed = -1;                         /* Random number seed */
   *equilibrate = 0;                  /* Default equilibration sweeps is 0 */
   *sourcetype = 0;                   /* 0 = LAMMPS, 1 = MC */
   *num_test_particles = 0;           /* Number of test particles to insert per frame */
   *kt = 0;                           /* kT */
   tmp_grid = 0;                     /* By default, use random points */

   /*--- 1. Read in values ---*/

   infile = fopen("options", "r");
   if (infile == NULL) die ("Could not open \"options\" file");
   printf ("Reading run options from the \"options\" file\n\n");

   while ( read_line(infile) ) {
      get_string(command, sizeof(command));
      upper_case(command);

      if (*command == '#') {
         continue;

      } else if (strcmp(command, "FILENAME") == 0) {
         if(!get_string(fname, sizeof(fname))) die ("Could not read filename after FILENAME");
         strcpy(configurations_filename, fname);
      } else if (strcmp(command, "POTENTIAL_FILE") == 0) {
         if(!get_string(fname, sizeof(fname))) die ("Could not read filename after POTENTIAL_FILENAME}");
         strcpy(potential_filename, fname);
      } else if (strcmp(command, "DEFINITION") == 0) {
         get_string(option, sizeof(option));
         upper_case(option);
         if (strcmp(option, "SPANNING") == 0) {
            *periodic = 0;
         } else if (strcmp(option, "WRAPPING") == 0) {
            *periodic = 1;
         } else {
            sprintf (error, "Unrecognised option after DEFINITION keyword: %s", option);
            die (error);
         }

      } else if (strcmp(command, "EQUILIBRATE") == 0){
         if (!get_int(equilibrate)) die ("Could not read equilibration sweeps after EQUILIBRATE");
      } else if (strcmp(command, "SEED") == 0) {
         if (!get_int(&seed)) die ("Could not read random number seed after SEED");
         seed = -abs(seed);   
      } else if (strcmp(command, "STATISTICS") == 0) {
         if (!get_int(report)) die ("Could not read number sweeps between statistics reports after STATISTICS");
      } else if (strcmp(command, "SWEEPS") == 0) {
         if (!get_int(nsweeps)) die ("Could not read number of sweeps after SWEEPS");
      } else if (strcmp(command, "SOURCETYPE") == 0) {
         if (!get_int(sourcetype)) die("Could not read source file type (LAMMPS or MC).");
      } else if (strcmp(command, "NUM_TEST_PARTICLES") == 0) {
         if (!get_int(num_test_particles)) die("Could not read number of test particles after NUM_TEST_PARTICLES");
      } else if (strcmp(command, "KT") == 0) {
         if (!get_double(kt)) die("Could not read value of kt after KT");
      } else if (strcmp(command, "GRID") == 0) {
         if (!get_int(&tmp_grid)) die ("Could not read grid after GRID");
         if (tmp_grid == 0) {*grid = 0;}
         else {*grid = 1;}

      } else {
         sprintf (error, "Unrecognised keyword: %s", command);
         die (error);
      }
   }

   fclose(infile);

   datafile = fopen(configurations_filename, "r");
   if (datafile == NULL) die ("Could not open configurations file");
   
   for (i = 0; i < 9; i++) {
      read_line(datafile);

      if (i==3) {
         if(!get_int(npart)) die ("Could not read number of particles from LAMMPS configuration file");
      }

      if (i==5) {
         double tmp_xa, tmp_xb;
         if (!get_double(&tmp_xa)) die ("Could not read first box size in x direction from LAMMPS configuration file");
         if (!get_double(&tmp_xb)) die ("Could not read second box size in x direction from LAMMPS configuration file");
         box->x = tmp_xb - tmp_xa;
      }

      if (i==6) {
         double tmp_ya, tmp_yb;
         if (!get_double(&tmp_ya)) die ("Could not read first box size in y direction from LAMMPS configuration file");
         if (!get_double(&tmp_yb)) die ("Could not read second box size in y direction from LAMMPS configuration file");
         box->y = tmp_yb - tmp_ya;
      }
   }

   fclose(datafile);
   
   

   /*--- 2. Validity checks ---*/

   if (*npart < 1) {
      die ("The number of discs must be at least 1.");
   }

   if (*equilibrate < 0) {
      die ("Equilibration time must be greater than or equal to 0.");
   }

   if (*num_test_particles < 1) {
      die ("Need at least 1 test particle per frame");
   }

   if (seed == 0) {
      die ("The random seed must be a negative integer (not zero).");
   }

   if (*kt <= 0.0) {
      die ("kt must be greater than 0.0");
   }

   if (*report > *nsweeps) *report=*nsweeps;



   /*--- 3. Summarize results on standard output ---*/
   printf (" Configurations file:                      %s\n", configurations_filename);
   printf (" Simulation cell dimensions:               %.8lf, %.8lf\n",
           box->x, box->y);
   printf (" Total number of particles:                %ld\n", *npart);
   printf (" Total number of test particles per frame: %ld\n", *num_test_particles);
   //sweeps
   printf (" Total sweeps in data:                     %ld\n", *nsweeps);
   printf (" Equilibration sweeps:                     %ld\n", *equilibrate);
   printf (" Total sweeps for statistics:              %ld\n", *nsweeps - *equilibrate);
   printf (" Sweeps between statistics reports:        %ld\n", *report);
   
   printf (" Random number seed:                       %ld\n", seed);
   printf (" kt:                                       %lf\n", *kt);
   printf (" Insertion positions:                      %s\n", *grid ? "Grid" : "Random");
   printf ("\n");

}

/*..............................................................................*/
