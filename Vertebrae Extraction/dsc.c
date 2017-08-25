/*********************************************************************/
/* dsc - computes the Dice coefficient score.				         */
/*********************************************************************/

#include "VisXV4.h"          /* VisionX structure include file       */
#include "Vutil.h"           /* VisionX utility header files         */

VXparam_t par[] =            /* command line structure               */
{
    {    "if=",    0,   " binary mask file (i.e. solution) "},
    {    "ig=",    0,   " user output file "},
    {     0,       0,    0}
};

#define  MVAL   par[0].val
#define  IVAL   par[1].val

V3fstruct (soln);
V3fstruct (user);

int
main(argc, argv)
int argc;
char *argv[];
{
    VXparse(&argc, &argv, par); /* parse the command line    */
    V3fread(&soln, MVAL);
    V3fread(&user, IVAL);

    int         x,y,z;           /* index counters            */
    float       X = 0;
    float       Y = 0;
    float	U = 0;	

    for (z = soln.zlo; z<= soln.zhi; z++) {
        for (y = soln.ylo; y <= soln.yhi; y++) {
            for (x = soln.xlo; x <= soln.xhi; x++) {
		if (soln.u[z][y][x] == 255)
			X += 1;
		if (user.u[z][y][x] == 255)
			Y += 1;
                if (soln.u[z][y][x] == 255 && user.u[z][y][x] == 255)
                    	U += 1;
            }
        }
    }

	printf("%f\n", 2*U/(X+Y));
    
    exit(0);
}
