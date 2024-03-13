#include "global.h"
#include "prototypes.h"

long seed; 

int  read_line(FILE *);
int  get_string(char [], int);
void upper_case(char []);
int  get_int(long int *);
int  get_double(double *);

void load(FILE *configurations_file, struct disc *particle, long npart,
         long ncellx, long ncelly, struct disc **cfirst, struct vector box, long sourcetype){
    if (sourcetype == 0){
        load_LAMMPS(configurations_file, particle, npart, ncellx, ncelly, cfirst, box);
    } else if (sourcetype == 1) {
        load_MC(configurations_file, particle, npart, ncellx, ncelly, cfirst, box);
    } else {
        die ("Could not determine source file type on load.");
    }
}

void load_LAMMPS(FILE *configurations_file, struct disc *particle, long npart,
         long ncellx, long ncelly, struct disc **cfirst, struct vector box){
    
    long i;
    //int j;
    double xpos, xpos_scaled;
    double ypos, ypos_scaled;
    //double xmin, xmax, ymin, ymax;
    long cell;
    long tmp_id;
    long type;
    int headerln;

    //xmin = xmax = ymin = ymax = 0;

    for (headerln = 0; headerln < 9; headerln++) {
        read_line(configurations_file);
    }

    for(i = 0; i<npart;i++) {
        
        if(read_line(configurations_file)){
            //id type x y....
            if(!get_int(&tmp_id)) die ("Could not read particle ID");
            if(!get_int(&type)) die ("Could not read particle type");
            if(!get_double(&xpos)) die ("Could not read particle x position");
            if(!get_double(&ypos)) die ("Could not read particle y position");

            xpos_scaled = xpos/box.x - floor(xpos/box.x);
            ypos_scaled = ypos/box.y - floor(ypos/box.y);

            //printf("%lf %lf\n", xpos_scaled, ypos_scaled);

            //input file is in scaled coords 0 to 1, convert to -0.5 to 0.5
            particle[i].pos.x = xpos_scaled * box.x;
            particle[i].pos.y = ypos_scaled * box.y;

            //for testing with scaled coords in LAMMPs file
            //particle[i].pos.x = xpos-0.5;
            //particle[i].pos.y = ypos-0.5;
            
            //type of structure the particle belongs to
            particle[i].structure = type;

            //if particle positions are still out of bounds, throw error
            if (particle[i].pos.x > box.x || particle[i].pos.x < 0.0){
		        printf("Particle %ld: x position out of bounds\n", i+1);
		        printf("Particle: %ld, x: %lf, y: %lf\n", i+1, particle[i].pos.x, particle[i].pos.y);
		        die("Particle x position outside of scaled coordinates bounds");
            }
	        if (particle[i].pos.y > box.y || particle[i].pos.y < 0.0) {
		        printf("Particle %ld: y position out of bounds\n", i+1); 
		        printf("Particle: %ld, x: %lf, y: %lf\n", i+1, particle[i].pos.x, particle[i].pos.y);
		        die("Particle y position outside of scaled coordinates bounds");
            }

            //read in neighbouring particle ids
            /*
            for (j = 0; j < 10; j++) {
                long tmp_n; 
                if (!get_int(&tmp_n)) die ("Could not read in neighbouring particle ID");
                particle[i].neighbours[j] = tmp_n - 1;
            }
            */
            //get cells

            if (ncellx > 3 && ncelly > 3) {
                cell = getcell(particle[i].pos, ncellx, ncelly, box);
                particle[i].cell = cell;
                particle[i].next = cfirst[cell];
                cfirst[cell] = &particle[i];
            }

        //debugging tool for getting max and min positions in each direction to test out of bounds
        /* 
	    if (particle[i].pos.x > xmax) {xmax = particle[i].pos.x;}
	    if (particle[i].pos.x < xmin) {xmin = particle[i].pos.x;}
	    if (particle[i].pos.y > ymax) {ymax = particle[i].pos.y;}
	    if (particle[i].pos.y < ymin) {ymin = particle[i].pos.y;}
        */

        } else {
            die("Configuration could not be read in. Load.c: ln107");
        }
        
    //printf("i: %ld, x: %lf, y: %lf\n", i, particle[i].pos.x, particle[i].pos.y);    	
    }
    //printf("xmin: %lf, xmax: %lf\n", xmin, xmax);
    //printf("ymin: %lf, ymax: %lf\n", ymin, ymax);

    //printf("config read\n");
    

    
}


void load_MC(FILE *configurations_file, struct disc *particle, long npart,
         long ncellx, long ncelly, struct disc **cfirst, struct vector box){
    
    long i;
    //int j;
    double xpos;
    //double xpos_scaled;
    double ypos;
    //double ypos_scaled;
    //double xmin, xmax, ymin, ymax;
    long cell;
    long tmp_id;
    //long type;
    int headerln;

    //xmin = xmax = ymin = ymax = 0;

    for (headerln = 0; headerln < 9; headerln++) {
        read_line(configurations_file);
    }

    for(i = 0; i<npart;i++) {
        //printf("%ld\n", i);
        if(read_line(configurations_file)){
            //id type x y....
            if(!get_int(&tmp_id)) die ("Could not read particle ID");
            if(!get_double(&xpos)) die ("Could not read particle x position");
            if(!get_double(&ypos)) die ("Could not read particle y position");
            
            //xpos_scaled = xpos/box.x - floor(xpos/box.x);
            //ypos_scaled = ypos/box.y - floor(ypos/box.y);

            //printf("%lf %lf\n", xpos_scaled, ypos_scaled);

            //Some LAMMPS coorinates are out of scaled bounds so reset them to be at boundary
            /*
            if (xpos_scaled >= 1.0) {xpos_scaled = 0.999999;}
            if (xpos_scaled < 0.0) {xpos_scaled = 0.0;}
            if (ypos_scaled >= 1.0) {ypos_scaled = 0.999999;}
            if (ypos_scaled < 0.0) {ypos_scaled = 0.0;}
            */
           
            //input file is in scaled coords 0 to 1, convert to -0.5 to 0.5
            //particle[i].pos.x = xpos_scaled-0.5;
            //particle[i].pos.y = ypos_scaled-0.5;

            //for testing with scaled coords in LAMMPs file
            particle[i].pos.x = xpos;
            particle[i].pos.y = ypos;
            
            //type of structure the particle belongs to
            particle[i].structure = -1;

            //if particle positions are still out of bounds, throw error
            if (particle[i].pos.x > 0.5 || particle[i].pos.x < -0.5){
		        printf("Particle %ld: x position out of bounds\n", i+1);
		        printf("Particle: %ld, x: %lf, y: %lf\n", i+1, particle[i].pos.x, particle[i].pos.y);
		        die("Particle x position outside of scaled coordinates bounds");
            }
	        if (particle[i].pos.y > 0.5 || particle[i].pos.y < -0.5) {
		        printf("Particle %ld: y position out of bounds\n", i+1); 
		        printf("Particle: %ld, x: %lf, y: %lf\n", i+1, particle[i].pos.x, particle[i].pos.y);
		        die("Particle y position outside of scaled coordinates bounds");
            }

            //read in neighbouring particle ids
            /*
            for (j = 0; j < 10; j++) {
                long tmp_n; 
                if (!get_int(&tmp_n)) die ("Could not read in neighbouring particle ID");
                particle[i].neighbours[j] = tmp_n - 1;
            }
            */
            //get cells
            if (ncellx > 3 && ncelly > 3) {
                cell = getcell(particle[i].pos, ncellx, ncelly, box);
                particle[i].cell = cell;
                particle[i].next = cfirst[cell];
                cfirst[cell] = &particle[i];
            }


        //debugging tool for getting max and min positions in each direction to test out of bounds
        /* 
	    if (particle[i].pos.x > xmax) {xmax = particle[i].pos.x;}
	    if (particle[i].pos.x < xmin) {xmin = particle[i].pos.x;}
	    if (particle[i].pos.y > ymax) {ymax = particle[i].pos.y;}
	    if (particle[i].pos.y < ymin) {ymin = particle[i].pos.y;}
        */

        } else {
            die("Configuration could not be read in. Load.c: ln96");
        }
        
    //printf("i: %ld, x: %lf, y: %lf\n", i, particle[i].pos.x, particle[i].pos.y);    	
    }
    //printf("xmin: %lf, xmax: %lf\n", xmin, xmax);
    //printf("ymin: %lf, ymax: %lf\n", ymin, ymax);

    //printf("config read\n");
    

    
}


void load_potential(long potential_length, double **potential, FILE *potential_file, double *max_potential_distance, double *dr) {
    long i;

    // Load potential from file
    for (i = 0; i < potential_length; i++) {
        if (read_line(potential_file)) {
            if (!get_double(&(potential[i][0]))) die ("Could not read value of potential when loading potential.");
            if (!get_double(&(potential[i][1]))) die ("Could not read value of position when loading potential.");
        } else {
            die("Error reading line in potential file");
        }
    }

    *dr = (double) (potential[1][0] - potential[0][0]);
    *max_potential_distance = potential[(potential_length)-1][0];

    printf("Spacing in r: %lf\n", *dr);
    printf("Maximum separation: %lf\n", *max_potential_distance);

    // Print potential as test
    for (i = 0; i < potential_length; i++) {

        printf("%.10lf %.15lf\n", potential[i][0], potential[i][1]);

    }
}