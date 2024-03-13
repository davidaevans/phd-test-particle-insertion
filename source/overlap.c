#include "global.h"
#include "prototypes.h"

/*..............................................................................*/

/*
Returns the vector pointing from the centre of mass of particle 2 to the
centre of mass of the closest image of particle 1 in 2D.

If the standard separation in the simulation box in either direction is 
greater than 0.5, the image must be closer - that's what the anint 
subtraction takes care of

*/

struct vector image(struct vector r1, struct vector r2, struct vector box)
{
   struct vector r12;

   r12.x = (r1.x - r2.x)/box.x;
   r12.y = (r1.y - r2.y)/box.y;

   r12.x = box.x * (r12.x - anint(r12.x));
   r12.y = box.y * (r12.y - anint(r12.y));

   return r12;
}

/*..............................................................................*/

/*
Returns the scalar separation from the centre of mass of particle 2 to the
centre of mass of the closest image of particle 1 in 2D.

If the standard separation in the simulation box in either direction is 
greater than 0.5, the image must be closer - that's what the anint 
subtraction takes care of

*/

double imagesep(struct vector r1, struct vector r2, struct vector box)
{
   struct vector r12;

   r12.x = r1.x - r2.x;
   r12.y = r1.y - r2.y;

   r12.x = box.x * (r12.x - anint(r12.x));
   r12.y = box.y * (r12.y - anint(r12.y));
   //printf("x: %lf, y: %lf\n", r12.x, r12.y);

   return DOT(r12,r12);
}

/*..............................................................................*/
 
/*
Returns the nearest integer to its argument as a double precision number. e.g.
anint(-0.49) = 0.0 and anint(-0.51) = -1.0. Equivalent to the Fortran intrinsic
ANINT.
*/
 
double anint(double arg)
{
   if (arg < 0) {
      return (double)( (long)(arg-0.5) );
   } else {
      return (double)( (long)(arg+0.5) );
   }
}
                                
/*..............................................................................*/
