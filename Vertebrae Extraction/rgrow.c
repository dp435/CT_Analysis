/****************************************************************/
/* VisX4 program rgrow                                          */
/* Region grows bone region from initialized seedpoints.        */
/*                                                              */
/****************************************************************/
#include "VisXV4.h"          /* VisX structure include file     */
#include "Vutil.h"           /* VisX utility header files       */

extern char *VisXhist;

char * pname = "v3tpl";

VisXfile_t *VXin,*VXin2,     /* input file structure            */
           *VXout;           /* output file structure           */
VisXelem_t *VXlist,**VXpt;   /* VisX data structure             */
VisXelem_t *mlist;           /* VisX data structure             */

VisX3dim_t  sim;             /* source image structure          */
VisX3dim_t  rim;             /* result image structure          */
VisX3dim_t  bim;             /* bone boundary image structure   */
VisXiinfo_t imginfo;
float xres,yres,zres,ri,rs;
int thresh, edgeThresh, intensityDiff, k;
int i,j,k,xs,ys,zs;

void VX3frameset(VisX3dim_t *is, VisX3dim_t *ir);
void region_grow(int x, int y, int z, int zSeed);
void remove_discs(int lowerBound, int upperBound);
void disconnect_ribs(int lowerBound, int upperBound);

VXparam_t par[] = {
    {"if=",    0,   "input file     v3tpl: Apply mask and threshold"},
    {"of=",    0,   "output file"},
    {"th=",    0,   "threshold value (default: 1100)"},
    {"xs=",    0,   "X-coordinate of seed."},
    {"ys=",    0,   "Y-coordinate of seed."},
    {"zs=",    0,   "Z-coordinate of seed."},
    {"-v",     0,   "verbose flag"},
    { 0,       0,   0},
};

/* Command line parameters are accessed in code by vars below */
#define  IVAL    par[0].val
#define  OVAL    par[1].val
#define  TVAL    par[2].val
#define  XVAL    par[3].val
#define  YVAL    par[4].val
#define  ZVAL    par[5].val
#define  VFLAG   par[6].val

int
main(argc, argv)
int argc;
char *argv[];
{
  VisXelem_t *vptr = NULL, *mptr = NULL;
  
  VXparse(&argc, &argv, par);    /* parse the command line      */

  VXin  = VXopen(IVAL, 0);       /* open input file             */
  VXout = VXopen(OVAL, 1);       /* open the output file        */

  VXlist = VXread(VXin);         /* read input file             */

  if (TVAL) {
    thresh = atoi(TVAL);
  } else {
    fprintf(stderr, "Threshold not specified; defaulting to 1100.\n");
    thresh = 1100;
  }

  if (XVAL) {
    xs = atoi(XVAL);
  } else {
    fprintf(stderr,"X-seed point missing, exiting,\n");
    exit(1);
  }

  if (YVAL) {
    ys = atoi(YVAL);
  } else {
    fprintf(stderr,"Y-seed point missing, exiting,\n");
    exit(1);
  }

  if (ZVAL) {
    zs = atoi(ZVAL);
  } else {
    fprintf(stderr,"Z-seed point missing, exiting,\n");
    exit(1);
  }

  fprintf(stderr,"SEED @ X:%d,Y:%d,Z:%d\n",xs,ys,zs);

  if(VXNIL == (vptr = VXfind(VXlist, VX_PSHORT))){
    fprintf(stderr, "%s: no acceptable input image found, exiting.\n",pname);
    exit(1);
  }

  /* Initialize input image structure */
  VXset3dim(&sim, vptr, VXin);
  if(sim.chan != 1){
    fprintf(stderr, "%s: Multi-channel images are not supported.\n",pname);
    exit(1);
  }

  /* Create boundary image structure */
  VXmake3dim(&bim, VX_PBYTE, sim.bbx, sim.chan);

  /* Create result image structure */
  VXmake3dim(&rim, VX_PBYTE, sim.bbx, sim.chan);


  edgeThresh = 1200;
  intensityDiff = 100;

  region_grow(xs,ys,zs,zs);
  int upperBound = MIN(zs + 20,sim.zhi-2);
  int lowerBound = MAX(zs - 20,sim.zlo+2);
  remove_discs(lowerBound, upperBound);
  disconnect_ribs(lowerBound, upperBound);


  VX3frameset(&sim,&rim);
  
  VXwrite(VXout, rim.list);     /* write data                   */
  VXclose(VXin);                /* close files                  */
  VXclose(VXout);

  exit(0);
}

void region_grow(int x, int y, int z, int zSeed) {
  /* limit growth in Z direction. */
  int upperBound = MIN(zs + 20,sim.zhi-2);
  int lowerBound = MAX(zs - 20,sim.zlo+2);
  if (z > upperBound || z < lowerBound)
    return;

  rim.u[z][y][x] = 255;

  if (z+1 <= sim.zhi && rim.u[z+1][y][x] == 0 && sim.s[z+1][y][x] > thresh) {
    if (sim.s[z][y][x] > edgeThresh) {
        if (sim.s[z+1][y][x] > edgeThresh) {
            region_grow(x, y, z+1, zSeed);
        }
    } else
        region_grow(x, y, z+1, zSeed);
  }

  if (z-1 >= sim.zlo && rim.u[z-1][y][x] == 0 && sim.s[z-1][y][x] > thresh) {
    if (sim.s[z][y][x] > edgeThresh) {
        if (sim.s[z-1][y][x] > edgeThresh) {
            region_grow(x, y, z-1, zSeed);
        }
    } else
        region_grow(x, y, z-1, zSeed);
  }

  if (y+1 <= sim.yhi && rim.u[z][y+1][x] == 0 && sim.s[z][y+1][x] > thresh) {
    if (sim.s[z][y][x] > edgeThresh) {
        if (sim.s[z][y+1][x] > edgeThresh) {
            region_grow(x, y+1, z, zSeed);
        }
    } else
        region_grow(x, y+1, z, zSeed);
  }

  if (y-1 >= sim.ylo && rim.u[z][y-1][x] == 0 && sim.s[z][y-1][x] > thresh) {
    if (sim.s[z][y][x] > edgeThresh) {
        if (sim.s[z][y-1][x] > edgeThresh) {
            region_grow(x, y-1, z, zSeed);
        }
    } else
      region_grow(x, y-1, z, zSeed);
  }

  if (x+1 <= sim.xhi && rim.u[z][y][x+1] == 0 && sim.s[z][y][x+1] > thresh) {
    if (sim.s[z][y][x] > edgeThresh) {
        if (sim.s[z][y][x+1] > edgeThresh) {
            region_grow(x+1, y, z, zSeed);
        }
    } else
    region_grow(x+1, y, z, zSeed);
  }

  if (x-1 >= sim.xlo && rim.u[z][y][x-1] == 0 && sim.s[z][y][x-1] > thresh) {
    if (sim.s[z][y][x] > edgeThresh) {
        if (sim.s[z][y][x-1] > edgeThresh) {
            region_grow(x-1, y, z, zSeed);
        }
    } else
    region_grow(x-1, y, z, zSeed);
  }
}

void remove_discs(int lowerBound, int upperBound) {
  for (k=lowerBound+3; k<=upperBound-3; k++) {
    for (j=rim.ylo; j<=rim.yhi; j++) {
      for (i=rim.xlo; i<=rim.xhi; i++) {
        if (rim.u[k][j][i] == 255) {
          /* intensity: hi,(lo),hi*/
          if (sim.s[k-1][j][i] >= edgeThresh && sim.s[k][j][i] < edgeThresh
            && sim.s[k+1][j][i] >= edgeThresh)
            rim.u[k][j][i] = 0;

          /* intensity: hi,(lo),lo,hi */
          for (zs=0; zs<=1; zs++) {
            if (sim.s[k-1-zs][j][i] >= edgeThresh && sim.s[k-zs][j][i] < edgeThresh
            && sim.s[k+1-zs][j][i] < edgeThresh && sim.s[k+2-zs][j][i] >= edgeThresh)
              rim.u[k][j][i] = 0;
          }

          /* intensity: hi,(lo),lo,lo,hi */
          for (zs=0; zs<=2; zs++) {
            if (sim.s[k-1-zs][j][i] >= edgeThresh && sim.s[k-zs][j][i] < edgeThresh
            && sim.s[k+1-zs][j][i] < edgeThresh && sim.s[k+2-zs][j][i] < edgeThresh
            && sim.s[k+3-zs][j][i] >= edgeThresh)
              rim.u[k][j][i] = 0;
          }

          /* intensity: hi,(lo),lo,lo,lo,hi */
          for (zs=0; zs<=3; zs++) {
            if (sim.s[k-1-zs][j][i] >= edgeThresh && sim.s[k-zs][j][i] < edgeThresh
            && sim.s[k+1-zs][j][i] < edgeThresh && sim.s[k+2-zs][j][i] < edgeThresh
            && sim.s[k+3-zs][j][i] < edgeThresh && sim.s[k+4-zs][j][i] >= edgeThresh)
            rim.u[k][j][i] = 0;
          }

          /* intensity: hi,(lo),lo,lo,lo,lo,hi */
          for (zs=0; zs<=4; zs++) {
            if (sim.s[k-1-zs][j][i] >= edgeThresh && sim.s[k-zs][j][i] < edgeThresh
            && sim.s[k+1-zs][j][i] < edgeThresh && sim.s[k+2-zs][j][i] < edgeThresh
            && sim.s[k+3-zs][j][i] < edgeThresh && sim.s[k+4-zs][j][i] < edgeThresh
            && sim.s[k+5-zs][j][i] >= edgeThresh)
            rim.u[k][j][i] = 0;
          }

        }
      }
    }
  }
}

void disconnect_ribs(int lowerBound, int upperBound) {
    /* calculate centroid of each frame. */
  int x,y,z;
  for (z = lowerBound; z <= upperBound; z++) {
    double x_acc = 0;
    double y_acc = 0;
    double area = 0;
    for (y = rim.ylo; y <= rim.yhi; y++) {
      for (x = rim.xlo; x <= rim.xhi; x++) {
        if (rim.u[z][y][x] == 255) {
          x_acc += x;
          y_acc += y;
          area++;
        }
      }
    }
    if (x_acc != 0 || y_acc != 0) {
      int x_centroid = (int) (x_acc/area);
      int y_centroid = (int) (y_acc/area);

      for (j=rim.ylo; j<=rim.yhi; j++) {
        for (i=rim.xlo; i<=rim.xhi; i++) {
          if (i > x_centroid + 115 || i < x_centroid - 115)
            rim.u[z][j][i] = 0;
          if (j > y_centroid + 100 || j < y_centroid - 125)
            rim.u[z][j][i] = 0;
        }
      }
    }
  }
}


void VX3frameset(VisX3dim_t *is, VisX3dim_t *ir) {
/* VX3frameset: inserts frame markers in ir in the same location
   as in is

   This function assumes that both 3D images have the same number
   of images. If not, it will print an error and quit.
*/
    VisXelem_t *src, *dest;     
    VisXelem_t *fptr,*sptr;    /* ptrs into src list */
    VisXelem_t *dptr;          /* ptrs into dst list */

    /* ensure at the beginning of each list */
    src=VXfirst(is->list);
    dest=VXfirst(ir->list);

    if ( VXNIL == (fptr = VXfind(src, VX_FRAME)) ) {
        /* No frames, don't do anything */
        return;
    }
    /* start searching after the first frame marker */
    sptr = src;
    dptr = dest;
    while(VXNIL != (sptr = VXfind(sptr, VX_BBX)) ) {
        if (VXNIL == (dptr = VXfind(dptr, VX_BBX)) ) {
            fprintf(stderr, "Error: Image count not equal!\n");
            exit(1);
        }
        fptr = VXbfind(sptr, VX_FRAME);
        /* Go back one element from BBX to add frame marker */
        dptr = VXaddelem(dptr->prev,VX_FRAME,"",fptr->size);

        fptr = VXfind(fptr, VX_EFRAME);
        dptr = dptr->next->next;
        dptr = VXaddelem(dptr, VX_EFRAME, "", fptr->size);
        sptr = sptr->next;
    }
    if (VXNIL != VXfind(dptr, VX_BBX)) {
        fprintf(stderr, "Error: Image count not equal at end!\n");
        exit(1);
    }
}

