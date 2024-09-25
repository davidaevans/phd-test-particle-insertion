
void accumulate(struct mystat *q, double x);

double anint(double arg);

void baccumulate(struct mystat *q, double x, double x2, long n);

void die(char string[]);

void draw(FILE *outfile, struct vector box, long npart,
          struct disc *particle);

long getcell(struct vector pos, long ncellx, long ncelly, struct vector box);

struct vector image(struct vector r1, struct vector r2, struct vector box);

double imagesep(struct vector r1, struct vector r2, struct vector box);


double ran2(long *idum);

void read_options(long *npart, struct vector *box, long *nsweeps,
   long *report,
   int *periodic,
   long *equilibrate, char *configurations_filename,
   long *sourcetype, char *potential_filename, long *num_test_particles, double *kt, int *grid);

void simulate(long npart, struct vector *box, long nsweeps,
   int periodic, struct disc *particle, long equilibrate, char *configurations_filename, 
   long sourcetype, double max_potential_distance, double dr, double **potential, long potential_length, double kt,
   long num_test_particles, int grid);



void rdf(struct disc *particle,
         long *rdfhist,
         long npart,
         long nsweeps,
         long equilibrate,
         long numbins,
         double dr,
         double maxsep,
         struct vector box);

void normalise_and_write_rdf(struct mystat *rdf,
                             long numbins,
                             long npart,
                             long nsweeps,
                             long equilibrate,
                             double dr,
                             struct vector box,
                             double boltzmann_factor_mean);

void load(FILE *configurations_file, struct disc *particle, long npart,
         long ncellx, long ncelly, struct disc **cfirst, struct vector box,
         long sourcetype);

void load_LAMMPS(FILE *configurations_file, struct disc *particle, long npart,
         long ncellx, long ncelly, struct disc **cfirst, struct vector box);

void load_MC(FILE *configurations_file, struct disc *particle, long npart,
         long ncellx, long ncelly, struct disc **cfirst, struct vector box);

void load_potential(long potential_length, double **potential, FILE *potential_filename, double *max_potential_distance, double *dr);

struct vector get_test_particle_location(struct vector *box);

double interpolate(double r1, double r2, double v1, double v2, double r_target);

double calculate_energy(double **potential, long potential_length, double r, double dr);

long get_file_length(FILE *fp);

void test_cell_lists(long **neighbours, long ncells);

void test_interpolated_potential(double **potential, long potential_length, double dr);

void calculate_gr_cell_lists(long num_test_particles, struct vector *box, long ncellx, long ncelly, long **neighbour, struct disc **cfirst,
                             double **potential, long potential_length, double dr, double kt, struct mystat *rdf, long npart, struct disc *particle,
                             double maxsep, int grid, long ngridx, long ngridy, double gridspacingx, double gridspacingy, double *boltzmann_factor_mean);

void calculate_gr_no_cell_lists(long num_test_particles, struct vector *box, double **potential, long potential_length, 
                                double dr, double kt, struct mystat *rdf, long npart, struct disc *particle, double maxsep,
                                int grid, long ngridx, long ngridy, double gridspacingx, double gridspacingy, double *boltzmann_factor_mean);

struct vector get_grid_point(long ngridx, long ngridy, double gridspacingx, double gridspacingy, long index);