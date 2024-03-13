#include "global.h"
#include "prototypes.h"

/*..............................................................................*/

/*
Returns the area of a 2D disc of given diameter and shell
thickness "shell" in units of D^2
*/

double disc_area(double diameter)
{
   return SQ(diameter/2) * M_PI;
}

/*..............................................................................*/

/*
Returns the cell number for a position in reduced coordinated
(on the unit square -0.5<=x<0.5, -0.5<=y<0.5).
*/

long getcell(struct vector pos, long ncellx, long ncelly, struct vector box)
{
   return (long)( (pos.y/box.x) * ncelly) * ncellx + (long)( (pos.x/box.y) * ncellx);
}

/*..............................................................................*/
