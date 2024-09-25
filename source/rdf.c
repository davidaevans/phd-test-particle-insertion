#include "global.h"
#include "prototypes.h"


// FOLLOWING ALGORITH 7, pg.86 from Frenkel and Smit 2nd ed.
void rdf(struct disc *particle,
         long *rdfhist,
         long npart,
         long nsweeps,
         long equilibrate,
         long numbins,
         double dr,
         double maxsep,
         struct vector box) {

    
    long i,j;
    double sep, sep2;             /* separation of particles */ 
    int bin;

    for (i=0;i<npart-1;i++) {
        for (j=i+1; j<npart;j++){
            sep2 = imagesep(particle[i].pos, particle[j].pos, box);
            sep = sqrt(sep2);

            if(sep < maxsep) {
                bin = sep/dr;
                rdfhist[bin] += 2;
            }

        }
        //printf("i: %ld\n", i);
    }

}

//normalisation of the rdf called at the end
void normalise_and_write_rdf(struct mystat *rdf,
                             long numbins,
                             long npart,
                             long nsweeps,
                             long equilibrate,
                             double dr,
                             struct vector box,
                             double boltzmann_factor_mean) {

    long i;
    long sweeps;
    long nparteff;
    long n_samples;
    double *normalised_rdf;
    double tmp_denom;
    double denom;
    FILE *outfile;

    normalised_rdf = (double *)malloc(numbins * sizeof(double));
    

    for (i=0;i<numbins;i++) {
        normalised_rdf[i] = (double) rdf[i].mean / boltzmann_factor_mean;
    }

    //write the resulting normalised histogram to the rdf.dat file 
    // Use the centre of the bin as r value for printing g(r)
    outfile = fopen("rdf.dat", "w");
    if (outfile == NULL) die ("Could not open rdf file for writing");
    for (i = 0; i < numbins; i++) {
        fprintf(outfile, "%.10lf %.15lf %ld\n",
                (i+0.5)*dr,
                normalised_rdf[i],
                rdf[i].samples);

        fflush(outfile);
    }
    fclose(outfile);
}