#include "global.h"
#include "prototypes.h"

/*..............................................................................*/

/*
Dumps a configuration to the supplied file handle.  The director is normalised
to the length of the rod.
*/

void draw(FILE *outfile, struct vector box, long npart,
          struct disc *particle)
{
   long i;

   for (i=0; i<npart; i++) {
      fprintf (outfile,
         "%15.8le %15.8le 0.0   0.0 0.0 0.0\n",
         box.x * (particle[i].pos.x - anint(particle[i].pos.x)),
         box.y * (particle[i].pos.y - anint(particle[i].pos.y))
         );
   }
}

/*..............................................................................*/


/*
Writes the configuration as above, but in a different format for debugging whether the
percolating cluster detection is working correctly

Format:

sweep# xpos ypos radius percolating

*/

void write_config(FILE *testfile, struct vector box, long npart,
                  struct disc *particle, long sweep, int percolating)
{
   long i;

   for (i=0; i<npart; i++) {
      fprintf (testfile,
         "%ld %15.8le %15.8le %lf %d %d\n",
         sweep, 
         box.x * (particle[i].pos.x - anint(particle[i].pos.x)),
         box.y * (particle[i].pos.y - anint(particle[i].pos.y)),
         particle[i].diameter,
         percolating,
         particle[i].structure
      );
   }
}

/*..............................................................................*/
