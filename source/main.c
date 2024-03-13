#include "global.h"
#include "prototypes.h"

/*..............................................................................*/

int main()
{
   char configurations_filename[80];          /*Name of file where LAMMPS configs are */
   char potential_filename[80];/* Name of file for reading in the potential*/
   double dr;                  /* step in r for binning RDF */
   double **potential;         /* External potential array */
   double max_potential_distance; /* largest separation of particles for pairwise potential */
   double kt;
   int periodic;               /* 1=wrapping criterion, 0=spanning */
   int grid;
   long sourcetype;             /* 0 = LAMMPS, 1 = MC sim */
   long i;
   long npart;                 /* Number of particles */
   long nsweeps;               /* Number of production sweeps */
   long report;                /* Number of sweeps between statistics reports */
   long equilibrate;           /* Number of sweeps to wait before sampling statistics */
   long potential_length;            /* Number of points in the potential file */
   long num_test_particles;    /* Number of test particles per frame */
   struct disc *particle;  /* Configuration of entire system */
   struct vector box;          /* Simulation cell dimensions */
   FILE *potential_file;
   //struct mystat st;

   printf ("\nTest Particle Insertion - 2D Disks");
   printf ("\n--------------------------\n\n");
   

   /* Get user parameters */
   read_options(&npart, &box, &nsweeps, &report,
                &periodic, 
                &equilibrate, configurations_filename,
                &sourcetype, potential_filename, &num_test_particles, &kt,
                &grid);
   
   /* Set aside memory for the configuration */
   particle = (struct disc *)malloc(npart * sizeof(struct disc));
   for (i=0; i<npart; i++) { particle[i].idx=i; }
   potential = NULL;

   // Set up array for potential
   potential_file = fopen(potential_filename, "r");
   if (!potential_file) die(strcat("Could not open potential with filename: ", potential_filename));

   potential_length = get_file_length(potential_file);

   potential = (double **) malloc(potential_length * sizeof(double *));
   for (i=0; i<potential_length; i++){
      potential[i] = (double *) malloc(2* sizeof(double));
      potential[i][0] = potential[i][1] = 0.0;
   }

   // Load Potential
   load_potential(potential_length, potential, potential_file, &max_potential_distance, &dr);

   fclose(potential_file);

   test_interpolated_potential(potential, potential_length, dr);

   printf ("Starting simulation\n\n");
   fflush (stdout);

   simulate(npart, &box, nsweeps, periodic, particle, equilibrate, configurations_filename, sourcetype, max_potential_distance, dr, potential, potential_length, kt, num_test_particles, grid);
   

   // Free memory

   printf ("\nDone\n\n");


   return 0;
}

/*..............................................................................*/

/*
Print error message and exit.
*/

void die(char string[])
{
   fprintf (stderr, "\nERROR: %s\n\n", string);
   exit (1);
}

/*................................................................................*/
