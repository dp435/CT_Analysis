/****************************************************************/
/* VisX4 program maxcut                                         */
/*                                */
/* Finds max-cut plane that passes between each vertebral   */
/* column, and uses that to find the center of each vertebra.   */
/****************************************************************/
#include "VisXV4.h"          /* VisX structure include file     */
#include "Vutil.h"           /* VisX utility header files       */
#include "math.h"

extern char *VisXhist;

char * pname = "maxcut";

VisXfile_t *VXin,*VXin2, *VXin3;      /* input file structure            */
FILE * fp; VisXfile_t *VXout;         /* output file structure           */
VisXelem_t *VXlist,**VXpt;            /* VisX data structure             */
VisXelem_t *mlist, *splist;           /* VisX data structure             */

VisX3dim_t  sim;                      /* source image structure          */
VisX3dim_t  rim;                      /* result image structure          */
VisX3dim_t  mim;                      /* mask image structure            */
VisX3dim_t  spim;                     /* spcanal mask image structure    */
VisXiinfo_t imginfo;
float xres,yres,zres,ri,rs;

void VX3frameset(VisX3dim_t *is, VisX3dim_t *ir);
void quicksort(int cand[], int first, int last);
void rotateX(float arr[3], float angle);
int findPlaneCoords(float normal[3], int point[3], int x, int y);
void limitSearchSpace(int searchSpace[], int lowerBound, int upperBound);

VXparam_t par[] = {
  {"if=",    0,   "input CT image"},
  {"ig=",    0,   "vertebrae mask file"},
  {"sp=",    0,   "spinal canal mask file"},
  {"of=",    0,   "output file"},
  {"-v",     0,   "verbose flag"},
  { 0,       0,   0},
};

/* Command line parameters are accessed in code by vars below */
#define  IVAL    par[0].val
#define  MVAL    par[1].val
#define  SPVAL   par[2].val
#define  OVAL    par[3].val
#define  VFLAG   par[4].val

int
main(argc, argv)
int argc;
char *argv[];
{
  int i,j,k,x,y,z;
  int xmin,xmax,ymin,ymax,zmin,zmax;

  int sp_arr[400][3];

  int cand[13];                 /* array holding Z-coords of top/bottom
                                   surfaces of vertebral body. */
  int candIdx = 0;

  /* control param for how close candidates can be. */
  int restrictionWindow = 15;

  /* control params for plane width & length*/
  int width = 75;
  int length = 75;

  VisXelem_t *vptr = NULL, *mptr = NULL, *sptr = NULL;

  VXparse(&argc, &argv, par);    /* parse the command line      */

  VXin  = VXopen(IVAL, 0);       /* open input file             */
  VXout = VXopen(OVAL, 1);       /* open the output file        */

  VXlist = VXread(VXin);         /* read input file             */
  VXgetresinfo( &imginfo);
  VXgetrescale( &imginfo);
  xres = imginfo.xres;
  yres = imginfo.yres;
  zres = imginfo.zres;
  ri = imginfo.ri;
  rs = imginfo.rs;
  if( VFLAG) {
    fprintf(stderr, "img res = %f x %f x %f ri %f rs %f\n",
      imginfo.xres, imginfo.yres, imginfo.zres, imginfo.ri, imginfo.rs);
  }

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

  if ( SPVAL ) {
    if ( VFLAG ) {
      fprintf(stderr, "%s: Mask file specified, applying mask...\n",pname);
    }
    VXin3 = VXopen(SPVAL, 0);       /* open mask file              */

    /* Read spinal canal mask file */
    splist = VXread(VXin3);
    if (VXNIL == (sptr = VXfind(splist,VX_PBYTE)) ) {
      fprintf(stderr, "%s: Invalid format for mask file.\n",pname);
      exit(1);
    }

    VXset3dim(&spim, sptr, VXin3);

    /* Check if image and mask have same bounding box, warn if not 
       Note that this is not a problem for this program, so we don't
       do anything.                                                */
    if ( (sim.xlo != spim.xlo) || (sim.xhi != spim.xhi) ||
         (sim.ylo != spim.ylo) || (sim.yhi != spim.yhi) ||
         (sim.zlo != spim.zlo) || (sim.zhi != spim.zhi) ) {
      fprintf(stderr, "%s: bounding boxes do not match.\n",pname);
    }

    /* Load spinal canal coords into array. */
    int sp_idx = 0;
    for (k = spim.zlo; k <= spim.zhi; k++) {
      for (j = spim.ylo; j <= spim.yhi; j++) {
        for (i = spim.xlo; i <= spim.xhi; i++) {
          if (spim.u[k][j][i] == 255) {
            sp_arr[k][0] = i;
            sp_arr[k][1] = j;
            sp_arr[k][2] = k;
          }
        }
      }
    }
  }

  /* Apply mask if mask image is specified*/
  if ( MVAL ) {
    if ( VFLAG ) {
      fprintf(stderr, "%s: Mask file specified, applying mask...\n",pname);
    }
    VXin2 = VXopen(MVAL, 0);       /* open mask file              */

    /* Read mask file */
    mlist = VXread(VXin2);
    if (VXNIL == (mptr = VXfind(mlist,VX_PBYTE)) ) {
      fprintf(stderr, "%s: Invalid format for mask file.\n",pname);
      exit(1);
    }

    VXset3dim(&mim, mptr, VXin2);

    /* Check if image and mask have same bounding box, warn if not 
       Note that this is not a problem for this program, so we don't
       do anything.                                                */
    if ( (sim.xlo != mim.xlo) || (sim.xhi != mim.xhi) ||
         (sim.ylo != mim.ylo) || (sim.yhi != mim.yhi) ||
         (sim.zlo != mim.zlo) || (sim.zhi != mim.zhi) ) {
      fprintf(stderr, "%s: bounding boxes do not match.\n",pname);
    }

    /* Determine what regions overlap */ 
    if (sim.xlo > mim.xlo) {
        xmin = sim.xlo;
    } else {
        xmin = mim.xlo;
    }
    if (sim.xhi < mim.xhi) {
        xmax = sim.xhi;
    } else {
        xmax = mim.xhi;
    }
    if (sim.ylo > mim.ylo) {
        ymin = sim.ylo;
    } else {
        ymin = mim.ylo;
    }
    if (sim.yhi < mim.yhi) {
        ymax = sim.yhi;
    } else {
        ymax = mim.yhi;
    }
    if (sim.zlo > mim.zlo) {
        zmin = sim.zlo;
    } else {
        zmin = mim.zlo;
    }
    if (sim.zhi < mim.zhi) {
        zmax = sim.zhi;
    } else {
        zmax = mim.zhi;
    }
    /* Apply mask to source image */
    for (k = zmin; k <= zmax; k++) {
      for (j = ymin; j <= ymax; j++) {
        for (i = xmin; i <= xmax; i++) {
          if (mim.u[k][j][i] == 0) {
            /* set pixels to lowest value for short */
            sim.s[k][j][i] = -32767;
          }
        }
      }
    }
  }

  /* Create result image structure */
  VXmake3dim(&rim, VX_PBYTE, sim.bbx, sim.chan);

  /* Find first frame of interest. */
  int isFound = 0;
  int frameStart = 0;

  for (k = zmin; k <= zmax; k++) {
    x = sp_arr[k][0];
    y = sp_arr[k][1];
    z = sp_arr[k][2];
    for (j=0; j<length; j++) {
      for (i=-width/2; i<= width/2; i++) {
        if (mim.u[z][y+j][x+i] == 255 && isFound == 0) {
          frameStart = k;
          isFound = 1;
        }
        if (isFound == 1)
          break;
      }
      if (isFound == 1)
        break;
    }
    if (isFound == 1)
      break;
  }

  int searchSpace[spim.zhi];
  /* Initialize potential search space to all true.*/
  for (i=0; i<spim.zhi; i++)
    searchSpace[i] = 1;

  /* Add first candidate. */
  cand[0] = sp_arr[frameStart][2];
  limitSearchSpace(searchSpace, 0, frameStart+restrictionWindow);

  int px,py,pz;                         /* plane coords.            */
  float angle;                          /* angle of tile of plane.  */
  signed long int acc;                  /* accumulator for intensity of voxels
                                           that the plane intersects. */    
  signed long int maxcut[spim.zhi];     /* max-cut array */

  /* initialize maxcut array with max signed long. */
  for (k=0; k<spim.zhi; k++)
    maxcut[k] = -2147483647;

  for (k=frameStart+1; k<spim.zhi; k++) {
    x = sp_arr[k][0];
    y = sp_arr[k][1];
    z = sp_arr[k][2];

    for (angle=-M_PI/8.0; angle<=M_PI/8.0; angle+=M_PI/64.0) {
      acc = 0;
      float normal[3] = {0.0,0.0,1.0};
      rotateX(normal, angle);

      for (j=0; j<length; j++) {
        px = sp_arr[k][0];
        py = sp_arr[k][1]+j;
        pz = findPlaneCoords(normal, sp_arr[k], px, py);

        /* Out-of-bounds: just skip iteration. */
        if (pz < frameStart+1 || pz >= sim.zhi) {
          acc = -2147483647;
          break;
        }

        for (i=-width/2; i<= width/2; i++) {
          /* Skip non-vertebrae voxels. */
          if (mim.u[pz][py][px+i] != 0) {
            acc += 1;
          }
        }
      }
      maxcut[k] = MAX(maxcut[k], acc);
    }
    if (VFLAG)
      fprintf(stderr,"FRAME:%d, MAXCUT:%d\n",k,maxcut[k]);
  }

  /* Search for local maximas. */
  for (i=1; i<=13; i++) {
    acc = -2147483647;
    isFound = 0;
    for (k=0; k<spim.zhi; k++) {
      if (searchSpace[k] == 1 && maxcut[k] > acc) {
        acc = maxcut[k];
        candIdx = k;
        isFound = 1;
      }
    }
    if (isFound == 1) {
      cand[i] = candIdx;
      if (((float) candIdx) > ((float)spim.zhi)*0.5) {
        limitSearchSpace(searchSpace, MAX(0,candIdx-20), 
          MIN(spim.zhi,candIdx+18));
      }
      else
        limitSearchSpace(searchSpace, MAX(0,candIdx-restrictionWindow),
          MIN(spim.zhi,candIdx+restrictionWindow));
    } else {
      cand[i] = -2147483647;
    }
  }

  quicksort(cand,0,12);
  int numMissed = 0;            /* number of local maximas missed. */
  for (k=0; k<13; k++) {
    if (cand[k] < 0)
      numMissed++;
  }

  int missedDist, potentialCand, maxCand;
  /* Extrapolate unidentified maximas. */
  for (i=0; i<numMissed; i++) {
    missedDist = -2147483647;
    candIdx = 0;
    /* Search upper portions of vertebrae. */
    for (k=0; k<6; k++) {
      if (cand[k] < 0)
        continue;
      if (cand[k+1] - cand[k] > missedDist) {
        missedDist = cand[k+1] - cand[k];
        candIdx=k;
      }
    }

    /* Include searching for lower portions of vertebrae. */
    for (k=0; k<12; k++) {
      if (cand[k] < 0)
        continue;
      if (cand[k+1] - cand[k] > 35) {
        missedDist = cand[k+1] - cand[k];
        candIdx=k;
      }
    }

    /* Calculate midpoint of widest gap, and check nearby pixels of the mdpt. */
    potentialCand = (cand[candIdx]+cand[candIdx+1])/2;
    maxCand = potentialCand;
    for (k=potentialCand-5; k>= potentialCand+5; k++) {
      if (maxcut[k] < maxcut[maxCand])
        maxCand = k;
    }
    cand[0] = maxCand;
    quicksort(cand,0,12);
  }

  numMissed = 0;            /* number of local maximas missed. */
  for (k=0; k<3; k++) {
    if (cand[k+1] - cand[k] > 28)
      numMissed++;
  }

  for (i=0; i<numMissed; i++) {
    candIdx = -1;
    /* Search upper portions of vertebrae. */
    for (k=0; k<3; k++) {
      if (cand[k+1] - cand[k] > 28) {
        candIdx=k;
      }
    }
    potentialCand = (cand[candIdx]+cand[candIdx+1])/2;
    maxCand = potentialCand;
    for (k=potentialCand-5; k>= potentialCand+5; k++) {
      if (maxcut[k] < maxcut[maxCand])
        maxCand = k;
    }
    cand[12] = maxCand;
    quicksort(cand,0,12);
  }


  if (VFLAG)
    for (k=0; k<13; k++)
      fprintf(stderr,"CANDIDATE ID:%d, Z:%d\n",k,cand[k]);

  fp = fopen(OVAL,"w");

  for (k = 0; k < 13-1; k++) {
    candIdx = (cand[k+1]+cand[k])/2;
    x = sp_arr[candIdx][0];
    y = sp_arr[candIdx][1];
    z = sp_arr[candIdx][2];
    fprintf(fp, "%d,%d,%d\n",x,y,z);
  }
  fclose(fp);

  VX3frameset(&sim,&rim);

//  VXwrite(VXout, rim.list);     /* write data                   */
  VXclose(VXout);

  VXclose(VXin);                /* close files                  */

  exit(0);
}

/* Method to sort 1D-array in increasing order. */
void quicksort(int cand[], int first, int last) {
  int pivot,j,temp,i;

  if(first<last){
    pivot=first;
    i=first;
    j=last;

    while(i<j){
      while(cand[i]<=cand[pivot] && i<last)
        i++;
      while(cand[j]>cand[pivot])
        j--;
      if(i<j) {
        temp=cand[i];
        cand[i]=cand[j];
        cand[j]=temp;
      }
    }

    temp=cand[pivot];
    cand[pivot]=cand[j];
    cand[j]=temp;

    quicksort(cand,first,j-1);
    quicksort(cand,j+1,last);
    }
}

/* Rotates a 3D vector by <angle> degrees (in radians) about the X-axis. */
void rotateX(float arr[3], float angle) {
  float y = arr[1];
  float z = arr[2];
  arr[1] = cos(angle) * y - sin(angle) * z;
  arr[2] = sin(angle) * y + cos(angle) * z;
}

/* Finds plane with normal that intersects point, and then find z-coords at (x,y) on the plane. */
int findPlaneCoords(float normal[3], int point[3], int px, int py) {
  float intercept = (normal[0] * ((float) point[0])
    + normal[1] * ((float) point[1]) + normal[2] * ((float) point[2]));
  return (int) ((intercept - normal[0] * ((float) px) - normal[1] * ((float) py))/normal[2]);
}

/* Restrict search space for potential candidates. */
void limitSearchSpace(int searchSpace[], int lowerBound, int upperBound) {
  int idx;
  for (idx = lowerBound; idx <= upperBound; idx++)
    searchSpace[idx] = 0;
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
        /* Set pointer to the image (bbx, then image) */
        /* TODO: Implement variable number of element skipping */
        dptr = dptr->next->next;
        dptr = VXaddelem(dptr, VX_EFRAME, "", fptr->size);
        sptr = sptr->next;
    }
    if (VXNIL != VXfind(dptr, VX_BBX)) {
        fprintf(stderr, "Error: Image count not equal at end!\n");
        exit(1);
    }
}
