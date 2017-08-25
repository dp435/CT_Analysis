/****************************************************************/
/* csv-conv                                                     */
/*            Converts binary .vx image into .csv file.         */
/****************************************************************/

#include "VisXV4.h"          /* VisionX structure include file       */
#include "Vutil.h"           /* VisionX utility header files         */
#include <stdio.h>
#include <stdlib.h>

VXparam_t par[] =            /* command line structure               */
{
{    "if=",    0,   " input file  v3dmean: compute local mean"},
{    "of=",    0,   " output file "},
{     0,       0,    0}
};
#define  IVAL   par[0].val
#define  OVAL   par[1].val

int
main(argc, argv)
int argc;
char *argv[];
{
  FILE * fp;
  V3fstruct (im);
  int        x,y,z;               /* index counters                 */

  VXparse(&argc, &argv, par);     /* parse the command line         */

  V3fread( &im, IVAL);            /* read 3D image                  */
  if ( im.type != VX_PBYTE || im.chan != 1) {     /* check  format  */
     fprintf (stderr, "image not byte type or single channel\n");
     exit (1);
  }   

  fp = fopen(OVAL,"w");

  /* for all pixels */
  for (z = im.zlo; z <= im.zhi; z++) {
    for (y = im.ylo; y <= im.yhi; y++) {
      for (x = im.xlo; x <= im.xhi; x++) {
        if (im.u[z][y][x] == 255)
          fprintf(fp, "%d,%d,%d\n",x,y,z);
      }
    }
  }
  
  fclose(fp);
  exit(0);
}