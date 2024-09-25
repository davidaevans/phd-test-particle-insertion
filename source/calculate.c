#include "global.h"
#include "prototypes.h"


// Calculates the number of lines in a file and then resets the pointer back to the start
long get_file_length(FILE *fp){
    long lines = 0;
    int ch;
    while(!feof(fp)) {
        ch = fgetc(fp);
        if(ch == '\n')
        {
            lines++;
        }
    }
    rewind(fp);
    return lines;
}

/*
Returns a random location in the box as the test particle location
*/ 
struct vector get_test_particle_location(struct vector *box) {
    struct vector testp;

    testp.x = ran2(&seed) * box->y;
    testp.y = ran2(&seed) * box->x;

    return testp;
}

/*
Returns next grid point for test particle location
*/
struct vector get_grid_point(long ngridx, long ngridy, double gridspacingx, double gridspacingy, long index) {
    struct vector testp;
    long x, y;

    x = index % ngridx;
    y = floor(index / ngridx);

    if (y > ngridy) die("Calculated y grid index greater than number of grid points / get_grid_point");

    testp.x = x * gridspacingx;
    testp.y = y * gridspacingy;

    return testp;
}

/* 
Return value of potential via linear interpolation 
*/

double interpolate(double r1, double r2, double v1, double v2, double r_target)
{
   double grad, val;

   //printf("x1: %lf x2: %lf y1: %lf y2: %lf\n", r1, r2, v1, v2);

   grad = (v2-v1)/(r2-r1);

   val = v1 + (r_target-r1) * grad;
   //printf("target: %lf pot: %lf\n",r_target,val);
   return val;
}


/*
Returns energy of particle separation based of linear interpolation of potential
*/

double calculate_energy(double **potential, long potential_length, double r, double dr) {
    //Assumes potential values are left edge of bin
    long i;
    double r1, r2, v1, v2;
    double tmp_energy;
    
    // v1 = v2 = r1 = r2 = 0.0;

    if (r < potential[0][0]) {return 10000000000.0;} //die("Particle separation less than first potential spacing.");
    if (r >= (potential[potential_length-1][0] + dr)) return 0.0;

    i = floor((r - potential[0][0])/dr);
    // Deal with if between last bin and cutoff
    // printf("potential length: %lf, bin: %ld, r: %lf ", potential[potential_length-1][0] + dr, i, r);
    if (i == potential_length - 1) {i--;}
    if (i > potential_length) {die("Bin error in interpolation");}

    // Account for potential being in the middle of the bin but stored as left edges
    if (i > 0 && r < (potential[i][0] + dr / 2.0 )){
        i--;
    }
    
    r1 = potential[i][0] + dr / 2.0;
    // r1 = potential[i][0];
    v1 = potential[i][1];

    r2 = potential[i+1][0] + dr / 2.0;
    // r2 = potential[i+1][0];
    v2 = potential[i+1][1];

    // printf("left edge: %lf, right edge: %lf\n", r1, r2);

    tmp_energy = interpolate(r1, r2, v1, v2, r);
    // printf("  tmp_energy: %.15lf", tmp_energy);
    return tmp_energy;
}

void test_interpolated_potential(double **potential, long potential_length, double dr) {

    long i;
    double tmp_e, tmp_r;
    FILE *outfile = NULL;

    outfile = fopen("interpolated-potential.dat", "w");
    if (!outfile) {die ("Could not open interpolated-potential.dat for writing");}

    for (i = 0; i < 50000; i++) {
        tmp_r = (double) i / 1000.0;
        tmp_e = calculate_energy(potential, potential_length, tmp_r, dr);

        fprintf(outfile, "%.10lf %.15lf\n", tmp_r, tmp_e);
    }   
    fclose(outfile);
}