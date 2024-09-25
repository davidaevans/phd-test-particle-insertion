#include "global.h"
#include "prototypes.h"

/*..............................................................................*/

void simulate(long npart, struct vector *box, long nsweeps,
   int periodic, struct disc *particle, long equilibrate, char *configurations_filename, 
   long sourcetype, double max_potential_distance, double dr, double **potential, long potential_length,
   double kt, long num_test_particles, int grid)
{
   double maxsep;      /* Maximum separation for binning of RDF */
   double shorter;
   double tmp_energy; /* Temporary vairabel to store the value of sum of energies*/
   double tmp_exp;    /* Temporary variable to store the value of exp(-E/kt) */
   double r;
   double gridspacingx, gridspacingy;
   double boltzmann_factor_mean;
   double boltzmann_factor_sum;
   long ngridx, ngridy;
   long i, j;
   long bin;
   long ncellx, ncelly;      /* Number of cell-list cells in each direction */
   long ncells;        /* Total number of cell-list cells */
   long **neighbour = NULL;   /* List of neighbouring cells for each cell */
   long sweep;         /* Current sweep number */
   long numbins;        /* Number of bins in RDF histogram */
   long *rdfhist;      /* Array for RDF histogram */
   long *cell;         /* For looping over neighbouring cell */
   long test_cell;     /* Cell that the test inserted particle would be in */
   long tmp_neighbours;  /* Counter for number of neighbours */
   long tmp_cells; /* Counter for number of cells iterated through */
   struct disc *neighbouring_particle; /* Particle that neighbours the test inserted particle */
   struct disc **cfirst = NULL;    /* Array of pointers to first particle in cell list */
   struct vector test_particle; /* x/y position of the particle */
   struct vector r_cm;  /* Vector for separation of particles centres of mass */
   struct mystat *rdf;         /* Array containing mystats for rdf bins */
   
   struct mystat nullstat = {0.0, 0.0, 0, 0.0, 0.0};
   FILE *configurations_file = NULL;
   
   /*=== Initialise counters etc. ===*/

   boltzmann_factor_mean = boltzmann_factor_sum = 0.0;
   numbins = 0;
   maxsep = 0;
   rdfhist = NULL;

   shorter = box->y > box->x ? box->x : box->y;
   maxsep = shorter/2;
   numbins = maxsep/dr + 1;
   rdf = (struct mystat *)malloc(numbins * sizeof(struct mystat));
   //set all elements to zero as no guarantee they will all be accessed
   for (i=0;i<numbins;i++) rdf[i] = nullstat;
   


   /*=== Initialise cell list ===*/
   ncellx = (long)(box->x / (max_potential_distance)); 
   ncelly = (long)(box->y / (max_potential_distance));
   if (ncellx > 3 || ncelly > 3) {
      ncells = ncellx * ncelly;
      cfirst = (struct disc **)malloc(sizeof(struct disc *) * ncells);
      neighbour = (long **)malloc(sizeof(long *) * ncells);
      for (i=0; i<ncells; i++) {
         neighbour[i] = (long *)malloc(sizeof(long) * 10);
      }
      /* Work out neighbouring cells for each cell by pointer */
      /* Interior of box */
      for (i=0; i<ncellx; i++) {
         for (j=0; j<ncelly; j++) {
            neighbour[j*ncellx+i][0] = ((j-1) + (j==0?ncelly:0))*ncellx + ((i-1) + (i==0?ncellx:0));
            neighbour[j*ncellx+i][1] = ((j-1) + (j==0?ncelly:0))*ncellx + i;
            neighbour[j*ncellx+i][2] = ((j-1) + (j==0?ncelly:0))*ncellx + ((i+1) - (i==ncellx-1?ncellx:0));
            neighbour[j*ncellx+i][3] = j*ncellx + ((i-1) + (i==0?ncellx:0));
            neighbour[j*ncellx+i][4] = j*ncellx + i;
            neighbour[j*ncellx+i][5] = j*ncellx + ((i+1) - (i==ncellx-1?ncellx:0));
            neighbour[j*ncellx+i][6] = ((j+1) - (j==ncelly-1?ncelly:0))*ncellx + ((i-1) + (i==0?ncellx:0));
            neighbour[j*ncellx+i][7] = ((j+1) - (j==ncelly-1?ncelly:0))*ncellx + i;
            neighbour[j*ncellx+i][8] = ((j+1) - (j==ncelly-1?ncelly:0))*ncellx + ((i+1) - (i==ncellx-1?ncellx:0));
            neighbour[j*ncellx+i][9] = -1;  /* end token */
         }
      }
      if (!periodic) {
         printf("OVERWRITING BOUNDARIES WITHOUT PBC\n");
         /* Overwrite periodic results along the boundaries */
         /* Edges */
         for (i=1; i<ncellx-1; i++) {
            /* top */
            neighbour[i][0] = i-1;
            neighbour[i][1] = i;
            neighbour[i][2] = i+1;
            neighbour[i][3] = ncellx + (i-1);
            neighbour[i][4] = ncellx + i;
            neighbour[i][5] = ncellx + (i+1);
            neighbour[i][6] = -1;
            /* bottom */
            neighbour[(ncelly-1)*ncellx + i][0] = (ncelly-2)*ncellx + (i-1);
            neighbour[(ncelly-1)*ncellx + i][1] = (ncelly-2)*ncellx + i;
            neighbour[(ncelly-1)*ncellx + i][2] = (ncelly-2)*ncellx + (i+1);
            neighbour[(ncelly-1)*ncellx + i][3] = (ncelly-1)*ncellx + (i-1);
            neighbour[(ncelly-1)*ncellx + i][4] = (ncelly-1)*ncellx + i;
            neighbour[(ncelly-1)*ncellx + i][5] = (ncelly-1)*ncellx + (i+1);
            neighbour[(ncelly-1)*ncellx + i][6] = -1;
         }
         for (j=1; j<ncelly-1; j++) {
            /* left */
            neighbour[j*ncellx][0] = (j-1)*ncellx;
            neighbour[j*ncellx][1] = (j-1)*ncellx + 1;
            neighbour[j*ncellx][2] = j*ncellx;
            neighbour[j*ncellx][3] = j*ncellx + 1;
            neighbour[j*ncellx][4] = (j+1)*ncellx;
            neighbour[j*ncellx][5] = (j+1)*ncellx + 1;
            neighbour[j*ncellx][6] = -1;
            /* right */
            neighbour[(j+1)*ncellx-1][0] = j*ncellx-2;
            neighbour[(j+1)*ncellx-1][1] = j*ncellx-1;
            neighbour[(j+1)*ncellx-1][2] = (j+1)*ncellx-2;
            neighbour[(j+1)*ncellx-1][3] = (j+1)*ncellx-1;
            neighbour[(j+1)*ncellx-1][4] = (j+2)*ncellx-2;
            neighbour[(j+1)*ncellx-1][5] = (j+2)*ncellx-1;
            neighbour[(j+1)*ncellx-1][6] = -1;
         }
         /* Corners */
         /* Top left */
         neighbour[0][0] = 0;
         neighbour[0][1] = 1;
         neighbour[0][2] = ncellx;
         neighbour[0][3] = ncellx+1;
         neighbour[0][4] = -1;
         /* Top right */
         neighbour[ncellx-1][0] = ncellx-2;
         neighbour[ncellx-1][1] = ncellx-1;
         neighbour[ncellx-1][2] = 2*ncellx-2;
         neighbour[ncellx-1][3] = 2*ncellx-1;
         neighbour[ncellx-1][4] = -1;
         /* Bottom left */
         neighbour[ncellx*(ncelly-1)][0] = ncellx*(ncelly-2);
         neighbour[ncellx*(ncelly-1)][1] = ncellx*(ncelly-2) + 1;
         neighbour[ncellx*(ncelly-1)][2] = ncellx*(ncelly-1);
         neighbour[ncellx*(ncelly-1)][3] = ncellx*(ncelly-1) + 1;
         neighbour[ncellx*(ncelly-1)][4] = -1;
         /* Bottom right */
         neighbour[ncellx*ncelly-1][0] = ncellx*(ncelly-1) - 2;
         neighbour[ncellx*ncelly-1][1] = ncellx*(ncelly-1) - 1;
         neighbour[ncellx*ncelly-1][2] = ncellx*ncelly - 2;
         neighbour[ncellx*ncelly-1][3] = ncellx*ncelly - 1;
         neighbour[ncellx*ncelly-1][4] = - 1;
      }
      printf ("Cell list grid: %ld x %ld\n", ncellx, ncelly);
      printf ("Cell size x:    %lf\n", box->x / ncellx);
      printf ("Cell size y:    %lf\n\n", box->y / ncelly);
   } else {
      printf ("Too few cells for cell lists\n");
      ncells = 0;
   }

   // printf("\nCheck Cell Lists\n");
   // test_cell_lists(neighbour, ncells);
   // printf("\n\n");

   // Open configurations file

   configurations_file = fopen(configurations_filename, "r");
   if (!configurations_file) die("Could not open configurations file");

   // Set up grid spacing, if using
   if (grid) {
      ngridx = ngridy = ceil(sqrt(num_test_particles));
      gridspacingx = (double) box->x/ngridx;
      gridspacingy = (double) box->y/ngridy;

      printf("Using grid for insertion positions\n");
      printf("%ld x %ld grid points\nSpacing x: %lf, y: %lf\n\n", ngridx, ngridy, gridspacingx, gridspacingy);

      // for (i=0;i<num_test_particles;i++){
      //    test_particle = get_grid_point(ngridx, ngridy, gridspacingx, gridspacingy, i);
      //    printf("%lf %lf\n", test_particle.x, test_particle.y);
      // }
   } else {
      ngridx = ngridy = 0;
      gridspacingx = gridspacingy = 0.0;
   }

   //loop through configurations
   for (sweep = 1; sweep <= nsweeps; sweep++){
      printf("Frame: %ld\n", sweep);
      //skip equilibrate steps
      if (sweep <= equilibrate){
         printf("Equilibration frame - ignoring\n");
         fflush(stdout);
         load(configurations_file, particle, npart, ncellx, ncelly, cfirst, *box, sourcetype);
         continue;
      }


      //reset pointers for cell lists and histogram data
      for (i=0; i<npart; i++) { particle[i].next = NULL; }
      for (i=0; i<ncells; i++) { cfirst[i] = NULL; }

      //printf("Sweep: %ld\n", sweep);      
      load(configurations_file, particle, npart, ncellx, ncelly, cfirst, *box, sourcetype);

      if (ncells == 0){
         calculate_gr_no_cell_lists(num_test_particles, box, potential, potential_length, dr, kt, rdf, npart, particle, maxsep, grid, ngridx, ngridy, gridspacingx, gridspacingy, &boltzmann_factor_sum);
      } else {
         calculate_gr_cell_lists(num_test_particles, box, ncellx, ncelly, neighbour, cfirst, potential, potential_length, dr, kt, rdf, npart, particle, maxsep, grid, ngridx, ngridy, gridspacingx, gridspacingy, &boltzmann_factor_sum);
      }
      

   } //end of configurations loop
 

   //close LAMMPs file so it can be reset for next shell thickness
   fclose(configurations_file);

   //free up arrays
   free(cfirst);
   for (i=0; i<ncells; i++) {
      free(neighbour[i]);
   }
   free(neighbour);

   boltzmann_factor_mean = (double) boltzmann_factor_sum / (nsweeps * num_test_particles);

   normalise_and_write_rdf(rdf, numbins, npart, nsweeps, equilibrate, dr, *box, boltzmann_factor_mean);
   
}


/*..............................................................................*/

/*
Accumulate a value into the statistics and update the mean and rms values.
*/
void accumulate(struct mystat *q, double x)
{
   (*q).sum += x;
   (*q).sum2 += x*x;
   (*q).samples++;
   (*q).mean = (*q).sum / (*q).samples;
   (*q).rms = sqrt(fabs((*q).sum2 / (*q).samples -
                        (*q).sum * (*q).sum / (*q).samples / (*q).samples));
}

/*................................................................................*/

/*
Accumulate a batch of values into the statistics and update the mean and rms values.
*/
void baccumulate(struct mystat *q, double x, double x2, long n)
{
   if (LONG_MAX - (*q).samples < n) {
      fprintf (stderr, "ERROR: Overflow of sample number in baccumulate\n");
      exit (1);
   }
   (*q).sum += x;
   (*q).sum2 += x2;
   (*q).samples += n;
   if ((*q).samples > 0) {
      (*q).mean = (*q).sum / (*q).samples;
      (*q).rms = sqrt(fabs((*q).sum2 / (*q).samples -
                        (*q).sum * (*q).sum / (*q).samples / (*q).samples));
   }
}

/*................................................................................*/


void test_cell_lists(long **neighbours, long ncells) {
   long i, j;

   for (i=0;i<ncells;i++){

      printf("Cell: %ld\n", i);
      for (j=0;j<9; j++) {
         printf(" %ld ", neighbours[i][j]);
      }
      printf("\n");

   }
}

void calculate_gr_cell_lists(long num_test_particles, struct vector *box, long ncellx, long ncelly, long **neighbour, struct disc **cfirst,
                             double **potential, long potential_length, double dr, double kt, struct mystat *rdf, long npart, struct disc *particle,
                             double maxsep, int grid, long ngridx, long ngridy, double gridspacingx, double gridspacingy, double *boltzmann_factor_sum) {
   double tmp_energy, tmp_exp;
   double r;
   long i, j,k;
   long *cell, test_cell;
   long tmp_neighbours, tmp_cells;
   long bin;
   struct vector test_particle;
   struct vector r_cm;
   struct disc *neighbouring_particle;

   for (i = 0; i < num_test_particles; i++) {
      // Get test particle location
      if (grid) {
         test_particle = get_grid_point(ngridx, ngridy, gridspacingx, gridspacingy, i);
      } else {
         test_particle = get_test_particle_location(box);
      }
      // printf("\n\ntp-x: %lf y: %lf\n", test_particle.x, test_particle.y);
      // fflush(stdout);
      test_cell = getcell(test_particle, ncellx, ncelly, *box);
      // printf("test cell: %ld\n", test_cell); 

      // Loop through "neighbouring" particles to calculate new energies around test particle
      tmp_energy = 0.0;
      tmp_neighbours = 0;
      tmp_cells = 0;

      cell = &neighbour[test_cell][0];
      while (*cell >= 0) {
         // printf("    cell: %ld\n", *cell);
         // fflush(stdout);
         tmp_cells++;         
         neighbouring_particle = cfirst[*cell];
         
         while (neighbouring_particle) {
            // Don't need to test if neighbouring particle = test particle as haven't added the test particle to the cell
            r_cm = image(test_particle, neighbouring_particle->pos, *box);
            r = sqrt(DOT(r_cm, r_cm));
            tmp_energy += calculate_energy(potential, potential_length, r, dr);
            tmp_neighbours++;      

            neighbouring_particle = neighbouring_particle->next;
         }
         cell++;
         // printf("    cell++: %ld\n", *cell);
         // fflush(stdout);
      }


      // Calculate temporary value of exponent around this test particle
      tmp_exp = exp(-tmp_energy/kt);
      *boltzmann_factor_sum += tmp_exp;
      //printf("%le %le %le\n", tmp_energy, kt, tmp_exp);

      // Loop over all other particles up to cutoff of half length of box
      for (j = 0; j < npart; j++) {
         r_cm = image(test_particle, particle[j].pos, *box);
         r = sqrt(DOT(r_cm, r_cm));

         if (r > maxsep) {continue;}

         // Calculate rdf bin from separation
         bin = (long) floor(r/dr);
         // printf("r: %lf bin: %ld\n", r, bin);

         // Add energy to statistics
         accumulate(&(rdf[bin]), tmp_exp);

      }

      // Go to next test particle
   }

}

void calculate_gr_no_cell_lists(long num_test_particles, struct vector *box, double **potential, long potential_length, double dr, double kt, struct mystat *rdf, long npart, struct disc *particle,
                                double maxsep, int grid, long ngridx, long ngridy, double gridspacingx, double gridspacingy, double *boltzmann_factor_sum) {
   double tmp_energy, tmp_exp;
   double r;
   long i, j;
   long tmp_neighbours;
   long bin;
   struct vector test_particle;
   struct vector r_cm;

   for (i = 0; i < num_test_particles; i++) {
      // Get test particle location
      if (grid) {
         test_particle = get_grid_point(ngridx, ngridy, gridspacingx, gridspacingy, i);
      } else {
         test_particle = get_test_particle_location(box);
      }
      // printf("\n\ntp-x: %lf y: %lf\n", test_particle.x, test_particle.y);
      // fflush(stdout);

      // Loop through "neighbouring" particles to calculate new energies around test particle
      tmp_energy = 0.0;
      tmp_neighbours = 0;

      for (j = 0; j < npart; j++) {
         
         // Don't need to test if neighbouring particle = test particle as haven't added the test particle to the cell
         r_cm = image(test_particle, particle[j].pos, *box);
         r = sqrt(DOT(r_cm, r_cm));
         tmp_energy += calculate_energy(potential, potential_length, r, dr);
         tmp_neighbours++;      

      }


      // Calculate temporary value of exponent around this test particle
      tmp_exp = exp(-tmp_energy/kt);
      *boltzmann_factor_sum += tmp_exp;
      //printf("%le %le %le\n", tmp_energy, kt, tmp_exp);

      // Loop over all other particles up to cutoff of half length of box
      for (j = 0; j < npart; j++) {
         r_cm = image(test_particle, particle[j].pos, *box);
         r = sqrt(DOT(r_cm, r_cm));

         if (r > maxsep) {continue;}

         // Calculate rdf bin from separation
         bin = (long) floor(r/dr);
         // printf("r: %lf bin: %ld\n", r, bin);

         // Add energy to statistics
         accumulate(&(rdf[bin]), tmp_exp);

      }

      // Go to next test particle
   }

}