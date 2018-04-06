/******************************************************************************  
   Ntools header file

   14-11-00: Major overhaul from utrecht.h related to new INOVA installations
             First Ntools version
             
   27-06-01: Inclusion of "on the fly" shaped pulse generation routines.
   
   20-09-01: First version for V500, V600 and V800: specific Ntools_sys.h
             file whith define V500, V600, or V800 Needed
             
   29-01-02: readPbox_shape routine; adaptation of the SHAPE structure and 
             S[TXYZ]PULSE routines.
             
   02-07-02: Inclusion of event code (from event.h) for (partially) overlapping 
             pulses.
             
   13-09-02: Fixed error in Pbox_shape routine: call to Pbox had args -p and -l
             interchanged.

   05-11-02: Improved error reporting in event routines. 
   
   14-07-04: delay-pulse-delay function extended with blank and unblank event
   
   12-10-04: Bug fix: replaced POWER_DELAY with PRRF_DELAY in FINE_POWER
             
   12-10-04: Changed implementation S(TXYZ)PULSE:
             - remainder of rof1 time carried over into predelay of shaped pulse
             - if rof2>WFG_STOP_DELAY, rof2-WFG_STOP_DELAY is carried into post-
               delay of pulse
               
   27-01-05: Bug fix: replaced one POWER_DELAY with PWRF_DELAY in POWER, 
             DECOUPLE_POWER and PRESAT_POWER routines

   27-01-05: Moved spectrometer-specific things (like probe defs) to a separate
             file (e.g. Ntools_V800.h) which should be linked to Ntools_sys.h
             
   08-03-05: Added GRADZ definition.
   
   xx-11-06: Cpd structure and cpdon, cpdoff functions.

   15-10-07: Initialisation of ampfactor in shaped pulse routines.

   17-10-07: Adjustments of hard_shape, sinc_shape, hsinc_shape and hsinctr_shape: 
	     now all use new_RFshape for initialisation.
	     Implementation of TripleResPulses; on-the-fly initialization of all
	     possible triple-res pulses.
	     
   28-10-07: Adjusted inplementation of TripleResPulses.
   	     shaka6_shape routine

*******************************************************************************/

#include <stdio.h>
#include <math.h>
#include <group.h>


/*
#include <dirent.h> 
#include <sys/stat.h>: inclusion gives compile warnings involving definition of SEC 
*/
#define S_IRWXU         00700   /* read, write, execute: owner */



/***********************************/
/********* handies to have *********/
/***********************************/

#define MIN(a,b)        ((a<b) ? a : b)
#define MAX(a,b)        ((a>b) ? a : b)
#define TRUE            1
#define FALSE           0
#define VERBOSE         1
#define NOVERBOSE       0
#define ERROR           1
#define NOERROR         0

#ifndef PI
#define PI              3.1416
#endif

#define ERROR_FILE      stdout

#define FINEPOWER_MAXVAL 4095.0

char *strdup();


/***********************************/
/*********** Probe & power defs ****/
/***********************************/
#define PROBE_MAX_RF_CHANNEL 4
#define PROBE_MAX_GRAD_CHANNEL 3
#define PROBE_GRADX 0
#define PROBE_GRADY 1
#define PROBE_GRADZ 2
typedef struct ProbeDef {

         char   *name;                                  /* probe name */
         int    id;                                     /* id of probe */
         double max_power[PROBE_MAX_RF_CHANNEL];        /* max power on rf-channels */
         double max_decpower[PROBE_MAX_RF_CHANNEL];     /* max decouple power on rf-channels */
         double max_satpower[PROBE_MAX_RF_CHANNEL];     /* max presat power on rf-channels */
         double grad_factor[PROBE_MAX_GRAD_CHANNEL];    /* 1% DAC value of gradient channels */

        } ProbeDef;

static ProbeDef *probe=NULL;                            /* pointer to current probe */

/*****************************************************/
/*********** some specific defs **********************/
/*****************************************************/
#include "Ntools_sys.h"


#define MIN_DELAY       0.2e-6          /* Minimal executable duration of delay */
#define LOCK_ON         lk_sample()     /* lock statements */
#define LOCK_OFF        lk_hold()

/***************************************/
/*********** functionalities ***********/
/***************************************/

#define us( a ) ((double) a * 1.0e-6)
#define ms( a ) ((double) a * 1.0e-3)

#define NI1 ((int) (d2*getval("sw1")))  /* number of increment in first indirect */
#define NI2 ((int) (d3*getval("sw2")))  /* number of increment in second indirect */
#define NI3 ((int) (d4*getval("sw3")))  /* number of increment in third indirect */

#define SQRT15 3.8729833
#define SQRT3  1.7320508
double pw90CaCo = SQRT15/(4.0*(174.0-57.0)*C13FREQ);    /* 90 degree on Ca or CO with zero excitation on CO or Ca */
double pw180CaCo = SQRT3/(2.0*(174.0-57.0)*C13FREQ);    /* 180 degree on Ca or CO with zero excitation on CO or Ca */

double pw90CCo = SQRT15/(4.0*(174.0-35.0)*C13FREQ);     /* 90 degree on C or CO with zero excitation on CO or C */
double pw180CCo = SQRT3/(2.0*(174.0-35.0)*C13FREQ);     /* 180 degree on C or CO with zero excitation on CO or C */


#define GRADZ(shape, gt, glx, gly, glz)  \
        {setProbe(); \
         zgradpulse( (glz)*probe->grad_factor[PROBE_GRADZ], gt )\
        }


/***************************************/
/***** delay routine  ******************/
/***************************************/
void udelay(time, file, line, statement)
double time; char *file; int line; char *statement;
/* delay with checks for minimal delay times and reports value
   fid #, file and line # in case of invalid delay time
*/
{if (time > MIN_DELAY)
        G_Delay(DELAY_TIME, time, 0);
 else
        fprintf(ERROR_FILE,"Warning: %s line %d\nFID %d:  %s = %.7f\n",
                file, (int) line, (int) ix, statement, time);
}

/* redefine delay */
#undef delay
#define delay(time) udelay((time), __FILE__, __LINE__, #time)


/***************************************/
/**** checking routines ****************/
/***************************************/
#define DECOUPLE_CHECK( period, duration ) \
{if ((dm[period] == 'y' || dm2[period] == 'y' || dm2[period] == 'y' || \
      dm[period] == 's' || dm2[period] == 's' || dm2[period] == 's' \
     )  \
      && at > duration ) \
        {fprintf(ERROR_FILE,"Reduce acquisition (%f) to <= %f \n", at, duration); \
         exit( 1 ); \
        } \
}


#define CHECK( s1, s2 ) check_string( s1, #s1, s2 )
int check_string( string1, string1_name, string2 )
char *string1, *string2, *string1_name;
/* character by character check of string1 against string2
   string2 can contain stars to indicate arbitrary values
   string1 is expanded with last character if shorter than string2
   Ex. check_string( dmm, "dmm", "ccc*");
*/
{char *s1=string1, *s2=string2;

 while (*s2 != '\0')
        {if (*s2 == '*' || *s1 == *s2)
                {s2++;
                 if (*(s1+1) != '\0') {s1++;}
                }
         else
                {fprintf(ERROR_FILE,"Error: invalid %s (%s, required %s)\n", string1_name, string1, string2);
                 exit(1);
                }
        }
 return 0;
}
/*****************************************************/
/*********** on the fly shaped-pulse stuff ***********/
/*****************************************************/

#define MAXLEN 256
#define SKIP_SPACES(line) {while (*line != '\0' && (*line == ' ' || *line == '\t')) line++;}

typedef struct Shape {
        int    type;            /* shape type flag SHAPE_RF or SHAPE_DEC */
        int    nsteps;          /* number of elements in shape */
        double pw;              /* duration of the shape */
        double ampfactor;       /* amplitude scaling factor */
        double pwr, rf, B1max;  /* rf power, fine power level, and B1max */
        double dres, dmf, dcyc; /* Additional Pbox parameters read from header */
        double offset, mod;     /* offset (Hz) and modulation offset factor (0.0 - 1.0) of the shape */
        double *amp, *phase;    /* pointers to amplitude and phase arrays */
        char   *name;           /* VNMR shapelib filename */
        char   *waveType;       /* Shape waveform type; e.g. g3, sinc etc. */
	int    device;		/* Channel/device on which the pulse will be given */
        } SHAPE;

#define SHAPE_RF                0
#define SHAPE_DEC               1
static char *shape_type[2] = {"RF","DEC"};

#define SHAPE_RESOLUTION        0.2e-6
#define SHAPE_MAXAMPLITUDE      1023.0
#define SHAPE_PWR_NOT_DEFINED   -99.99


int list_shape( shape , fp )
SHAPE *shape; FILE *fp;
{if (shape != NULL) {
/* NB: line cannot be to long because pulsetool will crash */
        fprintf(fp,
                "--- %s (%s) ---\nn: %4d, pw: %5.1f us, off.:%7.1f Hz, mod.: %4.2f, fac.: %4.2f\n",
                shape->name, 
                shape->waveType,
                shape->nsteps, 
                shape->pw*1e6, 
                shape->offset,
                shape->mod,
                shape->ampfactor
                );
 	 return NOERROR;
 } else {
 	 fprintf(ERROR_FILE,"Error list_shape: undefined shape (NULL)\n");
	 return ERROR;
 }
}


SHAPE *new_shape()
{SHAPE *tmp;

 if ((tmp = (SHAPE *) malloc(sizeof(SHAPE))) == NULL)
        {return NULL;
        }
        
 tmp->type     = SHAPE_RF;
 tmp->nsteps   = 0;
 tmp->pw       = 0.0;
 tmp->pwr      = SHAPE_PWR_NOT_DEFINED;
 tmp->rf       = FINEPOWER_MAXVAL;
 tmp->B1max    = 0.0;
 tmp->ampfactor= 1.0;
 tmp->amp      = NULL;
 tmp->phase    = NULL;
 tmp->offset   = 0.0;
 tmp->mod      = 0.0;
 tmp->name     = "not defined";
 tmp->waveType = "not defined";
 tmp->device   = -1; /* denotes undefined */

 return tmp;
}

SHAPE *new_RFshape( char *waveType, int nsteps )
/* return a new SHAPE instance with allocated arrays amp, phase for nsteps 
   Fill the amplitude array with SHAPE_MAXAMPLITUDE, phase array with zero's
*/
{SHAPE *shape;
 int i;

 if ((shape = new_shape()) == NULL)
        {fprintf(ERROR_FILE, "ERROR - new_RFshape: no memory init (%s)\n", waveType);
         exit( 1 );
        }
        
 shape->type     = SHAPE_RF;
 shape->nsteps   = nsteps;
 shape->waveType = strdup( waveType );

 if ((shape->amp = (double *) malloc( sizeof(double) * shape->nsteps)) == NULL)
        {fprintf(ERROR_FILE, "ERROR - new_RFshape: no memory amplitude array (%s,%d)\n", waveType, shape->nsteps);
         exit( 1 );
        }
 if ((shape->phase = (double *) malloc( sizeof(double) * shape->nsteps)) == NULL)
        {fprintf(ERROR_FILE, "ERROR - new_RFshape: no memory phase array (%s,%d)\n", waveType, shape->nsteps);
         exit( 1 );
        }
 for (i=0; i != shape->nsteps; i++) {
 	 shape->amp[i]   = SHAPE_MAXAMPLITUDE;
	 shape->phase[i] = 0.0;
 	}
 shape->ampfactor = 1.0;

 return shape;
}

SHAPE *write_shape( shape, shapename )
SHAPE *shape; char *shapename;
{
 char   outname[MAXLEN], name[MAXLEN], *home;
 FILE   *fout;
 char   *getenv();
 int    i;

 if ((home = getenv("HOME")) == NULL)
        {fprintf(ERROR_FILE,"ERROR - shift_shape: retrieving home\n");
         return NULL;
        }
        
 sprintf(name, "%s.%s", shapename, shape_type[shape->type]);
 sprintf(outname,"%s/vnmrsys/shapelib/%s", home, name);
 if ((fout = fopen(outname,"w")) == NULL)
        {fprintf(ERROR_FILE,"ERROR - write_shape: opening %s\n", outname);
         return NULL;
        }
 shape->name = strdup(name);

 fprintf(fout,"# NTools generated %s (%s)\n", shape->name, shape->waveType );
 fprintf(fout,"# pw = %7.1f\n", shape->pw );
 fprintf(fout,"# pwr = %6.2f\n",shape->pwr );
 fprintf(fout,"# rf = %7.1f \n",shape->rf );
 fprintf(fout,"# B1max = %7.2f\n",shape->B1max );
 fprintf(fout,"# ampfactor = %4.2f\n", shape->ampfactor );

 for (i=0; i!= shape->nsteps; i++)
        {fprintf( fout, "%10.3f %10d.0    1.0 \n", shape->phase[i], (int) (shape->amp[i]+0.5));
        }

 fclose(fout);

 return shape;

}

SHAPE *convolute_shape( shape, offset, mod )
SHAPE *shape; double offset, mod;
{
 double phase, phase_step, roundf;
 int i;

/* phase shifter in 4096 steps (12 bit resolution) */
 roundf = 360.0/ 4096.0;
 
 phase_step = 360.0 * offset * shape->pw / ((double) (shape->nsteps-1));

 phase = (-mod)*phase_step*((double) (shape->nsteps-1));
 while (phase < 0.0) phase += 360.0;
 while (phase > 360.0) phase -= 360.0;

 for (i=0; i!=shape->nsteps; i++)
        { /* now round phase and store */
         shape->phase[i] = ((int)((phase/roundf)+0.5))*roundf;
         phase += phase_step;
         while (phase < 0.0) phase += 360.0;
         while (phase > 360.0) phase -= 360.0;
        }
        
 shape->offset = offset;
 shape->mod = mod;
        
 return shape;

}

double calcrf_shape( shape, pwref, comp )
/* calculate rf power for fine attenuator
   pwref: reference pw; has to have the same flipangle as pw of the shape
   comp: amplifier compression factor
*/
SHAPE *shape; double pwref, comp;
{
 shape->rf = (comp*FINEPOWER_MAXVAL*pwref) / shape->pw * 1.0/shape->ampfactor;
 shape->rf = (int) (shape->rf + 0.5);
 return shape->rf;
}               


SHAPE *hard_shape(pw, offset, mod, pwref, comp, name, verbose)
double pw, offset, mod, pwref, comp; char *name; int verbose;
/*      pw:     pulse width
        offset: offset (in Hz) of the pulse
        mod:    zero point of phase modulation pattern; begin (mod=0.0) to end (mod=1.0)
        pwref:  reference pw for calculation of the fineattenuator value; has to have the same flipangle as pw of the shape
        comp:   amplifier compression factor
        name:   VNMR shapelib name
        verbose: verbose flag
*/

{SHAPE *shape; int nsteps;

 nsteps    = (int) (pw/SHAPE_RESOLUTION + 0.5);
 shape     = new_RFshape("hard", nsteps);   
 shape->pw = SHAPE_RESOLUTION*nsteps;

 shape->ampfactor = 1.0;
 calcrf_shape( shape, pwref, comp );
 
 if (offset != 0.0)
        { convolute_shape( shape, offset, mod );
        }
 
 write_shape( shape, name );

 if (verbose)
        list_shape( shape, stdout );

 return shape;
}       


SHAPE *comp270_90_shape(double pw90, double pwr, char *name, int verbose)
/*      pw90:    90 pulse width
        pwr:     power corresponding to pw90
        name:    VNMR shapelib name
        verbose: verbose flag
*/
{SHAPE *shape; int i;

 shape = new_RFshape("comp270_90",  5 );
        
 shape->pw     = pw90*shape->nsteps;
 shape->pwr    = pwr;
 
/* compensating evolution of duration pw90 */
 shape->amp[0]   = 0.0;
 shape->phase[0] = 0.0;

/* 270 x-pulse */
 for (i=1; i!=4; i++)
        {shape->amp[i]   = 1.0;
	 shape->phase[i] = 0.0;
        }
/* 90 -x pulse */
 shape->amp[4]   = 1.0;
 shape->phase[4] = 180.0;
 
 write_shape( shape, name );

 if (verbose)
        list_shape( shape, stdout );

 return shape;
}       

SHAPE *sinc_shape(pw, offset, mod, pwref, comp, name, verbose)
double pw, offset, mod, pwref, comp; char *name; int verbose;
/*      pw:      pulse width
        offset:  offset (in Hz) of the pulse
        mod:     zero point of phase modulation pattern; begin (mod=0.0) to end (mod=1.0)
        pwref:   reference pw for calculation of the fineattenuator value; has to have the same flipangle as pw of the shape
        comp:    amplifier compression factor
        name:    VNMR shapelib name
        verbose: verbose flag
*/
{SHAPE *shape; int i, nsteps;
 double stepsize, x, amp;

 nsteps    = (int) (pw/SHAPE_RESOLUTION + 0.5);
 shape     = new_RFshape("sinc", nsteps);   
 shape->pw = SHAPE_RESOLUTION*nsteps;

/* sin(x)/x from -pi to pi (centerlobe) */
 stepsize = PI*2.0/((double) (shape->nsteps-1));
 x=-PI;
 shape->ampfactor = 0.0;
 for (i=0; i!=shape->nsteps; i++)
        {if (x > -0.01  && x < 0.01)
                {amp =1.0;
                }
         else
                {amp = sin(x)/x;
                }
         shape->ampfactor += amp;
         shape->amp[i] = amp * SHAPE_MAXAMPLITUDE;
         shape->phase[i] = 0.0;
         x += stepsize;
        }
 shape->ampfactor /= shape->nsteps;
 calcrf_shape( shape, pwref, comp );

 if (offset != 0.0)
        { convolute_shape( shape, offset, mod );
        }
 
 write_shape( shape, name );

 if (verbose)
        list_shape( shape, stdout );

 return shape;
}       


SHAPE *hsinc_shape(pw, offset, mod, pwref, comp, name, verbose)
double pw, offset, mod, pwref, comp; char *name; int verbose;
/*      pw:      pulse width
        offset:  offset (in Hz) of the pulse
        mod:     zero point of phase modulation pattern; begin (mod=0.0) to end (mod=1.0)
        pwref:   reference pw for calculation of the fineattenuator value; has to have the same flipangle as pw of the shape
        comp:    amplifier compression factor
        name:    VNMR shapelib name
        verbose: verbose flag
*/
{SHAPE *shape; int i, nsteps;
 double stepsize, x, amp;

 nsteps    = (int) (pw/SHAPE_RESOLUTION + 0.5);
 shape     = new_RFshape("hsinc", nsteps);   
 shape->pw = SHAPE_RESOLUTION*nsteps;

/* half sin(x)/x from -pi to 0 (half centerlobe) */
 stepsize = PI/((double) (shape->nsteps-1));
 x=-PI;
 shape->ampfactor = 0.0;
 for (i=0; i!=shape->nsteps; i++)
        {if (x > -0.01  && x < 0.01)
                {amp =1.0;
                }
         else
                {amp = sin(x)/x;
                }
         shape->ampfactor += amp;
         shape->amp[i] = amp * SHAPE_MAXAMPLITUDE;
         shape->phase[i] = 0.0;
         x += stepsize;
        }
 shape->ampfactor /= shape->nsteps;
 calcrf_shape( shape, pwref, comp );

 if (offset != 0.0)
        { convolute_shape( shape, offset, mod );
        }
 
 write_shape( shape, name );

 if (verbose)
        list_shape( shape, stdout );

 return shape;
}       


SHAPE *hsinctr_shape(pw, offset, mod, pwref, comp, name, verbose)
double pw, offset, mod, pwref, comp; char *name; int verbose;
/*      pw:      pulse width
        offset:  offset (in Hz) of the pulse
        mod:     zero point of phase modulation pattern; begin (mod=0.0) to end (mod=1.0)
        pwref:   reference pw for calculation of the fineattenuator value; has to have the same flipangle as pw of the shape
        comp:    amplifier compression factor
        name:    VNMR shapelib name
        verbose: verbose flag
*/
{SHAPE *shape; int i, nsteps;
 double stepsize, x, amp;

 nsteps    = (int) (pw/SHAPE_RESOLUTION + 0.5);
 shape     = new_RFshape("hsinctr", nsteps);   
 shape->pw = SHAPE_RESOLUTION*nsteps;

/* half sin(x)/x from -pi to 0 (half centerlobe) */
 stepsize = PI/((double) (shape->nsteps-1));
 x=0;
 shape->ampfactor = 0.0;
 for (i=0; i!=shape->nsteps; i++)
        {if (x > -0.01  && x < 0.01)
                {amp =1.0;
                }
         else
                {amp = sin(x)/x;
                }
         shape->ampfactor += amp;
         shape->amp[i] = amp * SHAPE_MAXAMPLITUDE;
         shape->phase[i] = 0.0;
         x += stepsize;
        }
 shape->ampfactor /= shape->nsteps;
 calcrf_shape( shape, pwref, comp );

 if (offset != 0.0)
        { convolute_shape( shape, offset, mod );
        }
 
 write_shape( shape, name );

 if (verbose)
        list_shape( shape, stdout );

 return shape;
}       

SHAPE *readPbox_shape(shapename, verbose)          
/* retrieve Pbox parameters from .RF file header 
   adapted from EK's getRsh routine in the Pbox_psh.h file
*/
char  *shapename; int verbose;
{
  SHAPE    *rshape;
  FILE     *inpf;
  int      j;
  char     ch, str[MAXLEN];
  char     outname[MAXLEN], name[MAXLEN], *home;
  char    *getenv();

  if ((rshape = new_shape()) == NULL)
        {fprintf(ERROR_FILE,"ERROR - readPbox_shape: allocating shape memory\n");
         return NULL;
        }

  if ((home = getenv("HOME")) == NULL)
        {fprintf(ERROR_FILE,"ERROR - readPbox_shape: retrieving home\n");
         return NULL;
        }
        
  sprintf(name, "%s.%s", shapename, shape_type[rshape->type]);
  sprintf(outname,"%s/vnmrsys/shapelib/%s", home, name);

  if ((inpf = fopen(outname,"r")) == NULL)
        {fprintf(ERROR_FILE,"ERROR - readPbox_shape: opening %s\n", outname);
         exit(1);
        }
  rshape->name = strdup(name);

  
  j = fscanf(inpf, "%c %s %lf %lf %lf %lf %lf %lf %lf\n", 
                   &ch, str, 
                   &(rshape->pw), 
                   &(rshape->pwr),  
                   &(rshape->rf), 
                   &(rshape->dres), 
                   &(rshape->dmf), 
                   &(rshape->dcyc), 
                   &(rshape->B1max)                
            );
            
  fclose(inpf);

  if ((j > 4) && (ch == '#') && (str[0] == 'P') && (str[1] == 'b') && 
      (str[2] == 'o') && (str[3] == 'x') && (rshape->pwr < 63.1))
  {
    rshape->pw /= 1000000.0;
  }
  else
        {fprintf(ERROR_FILE,"ERROR - readPbox_shape: retrieving parameters\n");
         exit(1);
        }
  rshape->waveType = strdup("Pbox");
        
  if (verbose)
        list_shape( rshape, stdout );


  return rshape; 
}

SHAPE *Pbox_shape(waveType, bandwidth, offset, pwref, refpwr, options, name, verbose)          
/* use Pbox to make a .RF shape.
   read name.RF file header into the structure.
   return structure pointer upon successfull completion.
   
        wavetype:       Pbox wave type; eg. g3, sinc, ...etc
        bandwidth:      Bandwidth of the pulse          
        offset:         offset (in Hz) of the pulse
        pwref:          reference pw90 for calculation of the power settings
        refpwr:         reference power 
        options:        optional Pbox flags/commands/parameters
        name:           VNMR shapelib name (no .RF extension)
        verbose:        verbose flag
      
*/
char  *waveType, *name, *options; double bandwidth, offset, pwref, refpwr; int verbose;
{ char    string[MAXLEN];
  char    *getenv();
  SHAPE   *shape;
 
/* generate the Pbox command */
  sprintf(string, "Pbox %s -w \"%s %.3f %.3f\" -l %.2f -p %.2f -stepsize %0.2f %s",
                   name,
                   waveType,
                   bandwidth,
                   offset,
                   pwref*1.0e6,
                   refpwr,
                   SHAPE_RESOLUTION*1e6,
                   options
         );
  if (verbose>1)
        {fprintf(ERROR_FILE,">>%s<<\n",string);
         fflush( ERROR_FILE );
         system( string );
        }
  else
        {sprintf(string, "%s > /dev/null", string);
         system(string);
        }
  
  if ((shape = readPbox_shape(name, 0)) == NULL)
        {fprintf(ERROR_FILE,"ERROR - Pbox_shape: retrieving header information\n");
         return NULL;
        }
  shape->offset = offset;
  shape->waveType = strdup(waveType);

  if (verbose)
        list_shape( shape, stdout );
  
  return shape;
}

SHAPE *shaka6_shape( pw90, pw90pwr, name, verbose )
double pw90, pw90pwr; char *name; int verbose;
/*      pw90:    basic pw90 to derive shaka6 from
	pw90pwr: power for pw90
        name:   VNMR shapelib name
        verbose: verbose flag
*/

{SHAPE *shape; int nsteps;

/* shaka 6 definition, flip angles in degrees; Chem. Phys. Lett. 120, 201 (1985) */
 double shaka6[6] = { 158.0, 171.2, 342.8, 145.5, 81.2, 85.3}; 
 double sum;
 int    count, n, i, j;

 nsteps    = (int) (pw90*(984.0/90.0)/SHAPE_RESOLUTION + 0.5);
 shape     = new_RFshape("shaka6", nsteps);   
 shape->pw = SHAPE_RESOLUTION*nsteps;
 
 /* set phase according to shaka6 pattern */
 sum = 0.0; count = 0;
 for (i=0; i<6; i++) {
 	/* prevent cummulative rounding errors */
 	n = (int) ( (sum+shaka6[i])/90.0*pw90 / SHAPE_RESOLUTION + 0.5) - count; 
	
	for (j=count; j<count+n; j++) {
		shape->phase[j] = (i%2)*180.0;
	}
	count += n;
	sum += shaka6[i];
 }

 write_shape( shape, name );

 if (verbose)
        list_shape( shape, stdout );

 return shape;
}       

SHAPE *comp244_shape( pw90, pw90pwr, name, verbose )
double pw90, pw90pwr; char *name; int verbose;
/*      90(x),2.44*90(y),90(x) composite 180

	pw90:    basic pw90 to derive shaka6 from
	pw90pwr: power for pw90
        name:   VNMR shapelib name
        verbose: verbose flag
*/

{SHAPE *shape; int nsteps;

/* comp244 definition, flip angles in degrees; phase in degrees */
 double comp244[3] = { 90.0, 2.44*90.0, 90.0 }; 
 double phase[3] =   {  0.0,      90.0,  0.0 }; 
 double sum;
 int    count, n, i, j;

 nsteps    = (int) (pw90*4.44/SHAPE_RESOLUTION + 0.5);
 shape     = new_RFshape("comp244", nsteps);   
 shape->pw = SHAPE_RESOLUTION*nsteps;
 
 /* set phase according to shaka6 pattern */
 sum = 0.0; count = 0;
 for (i=0; i<3; i++) {
 	/* prevent cummulative rounding errors */
 	n = (int) ( (sum+comp244[i])/90.0*pw90 / SHAPE_RESOLUTION + 0.5) - count; 
	
	for (j=count; j<count+n; j++) {
		shape->phase[j] = phase[i];
	}
	count += n;
	sum += comp244[i];
 }

 write_shape( shape, name );

 if (verbose)
        list_shape( shape, stdout );

 return shape;
}       

SHAPE *watergate180_shape( pw1, pw2, pw90pwr, wgdelay, name, verbose )
double pw1, pw2, pw90pwr, wgdelay; char *name; int verbose;
/*      pw1, pw2: basic pw90 to derive watergate from
	pw90pwr:  power for pw1/pw2
	wgdelay:  watergate delay
        name:     VNMR shapelib name
        verbose:  verbose flag
*/

{SHAPE *shape; int nsteps;

/* watergate definition, */
 double wgdef[11] = { pw1/13.0*3.0,  wgdelay, pw1/13.0*9.0, wgdelay, pw1/13.0*19.0, wgdelay,
                      pw2/13.0*19.0, wgdelay, pw2/13.0*9.0, wgdelay, pw2/13.0*3.0
		    }; 
 double sum;
 int    count, n, i, j;

 nsteps    = (int) ((pw1/13.0*31.0+wgdelay*3+pw2/13.0*31.0+wgdelay*2)/SHAPE_RESOLUTION + 0.5);
 shape     = new_RFshape("watergate180", nsteps);   
 shape->pw = SHAPE_RESOLUTION*nsteps;
 
 /* set amplitudes according to wgdef pattern */
 sum = 0.0; count = 0;
 
 /* first '90', phase 0 */
 for (i=0; i<5; i++) {
 	/* prevent cummulative rounding errors */
 	n = (int) ( (sum+wgdef[i]) / SHAPE_RESOLUTION + 0.5) - count; 
	
	for (j=count; j<count+n; j++) {
		shape->amp[j] = (i%2) ? 0.0 : SHAPE_MAXAMPLITUDE;
		shape->phase[j] = 0.0;
	}
	count += n;
	sum += wgdef[i];
 }
 /* second '90', phase 180 */
 for (i=5; i<11; i++) {
 	/* prevent cummulative rounding errors */
 	n = (int) ( (sum+wgdef[i]) / SHAPE_RESOLUTION + 0.5) - count; 
	
	for (j=count; j<count+n; j++) {
		shape->amp[j] = (i%2) ? 0.0 : SHAPE_MAXAMPLITUDE;
		shape->phase[j] = 180.0;
	}
	count += n;
	sum += wgdef[i];
 }

 write_shape( shape, name );

 if (verbose)
        list_shape( shape, stdout );

 return shape;
}       

/**********************************************************/
/*********** triple res stuff *****************************/
/**********************************************************/

/* Pulse file names are defined as:
  
  pulseFileName :== dirname / prefix + shapeName
  
  shapeName :== shapePart + resonancePart

  dirname:       directory name for storage of all shapes; gets created in the vnmrsys/shapelib directory
  prefix: 	 user-defined prefix to the name; e.g. "CACO" to indicate Carrier on Ca, off-resonance Co.
  shapePart:     describes kind of pulse; e.g. hard, g3, sincCoff, zeroCoff.
  resonancePart: describes the frequency, e.g H90 for proton 90 on resonance, C180 for carbon 180 on resonance,
                 Coff90, carbon 90 frequency-shifted to CoffCar.

*/
typedef struct TripleResPulses {

         char *dirname;                         /* directory name for storage of all shapes */
         char *prefix;                          /* prefix for all shapes: dirname/prefix+shapeName */

/* TODO: make nomenclature consistent */
        /* proton */
         double Hfrq, pwHref, HpwrRef, compH;   /* proton freq, ref pw, ref power, compH */
         SHAPE  *hardH90, *hardH180;            /* hard proton 90, 180 pulses on resonance */
         SHAPE  *sincH2OH90;                    /* sinc pulse for H2O */
         
        /* carbon */
         double Cfrq, pwCref, refpwrC, compC;   /* carbon freq, ref pw, ref power, compC */
         double Ccar, CoffCar;                  /* carbon carrier (ppm), Coff in ppm */
         SHAPE  *hardC90,     *hardC180;        /* hard carbon 90, 180 pulses on resonance */
	 SHAPE  *shaka6C180,  *comp244C180;	/* shaka6 composite 180, 90(x),2.44*90(y),90(x) comp 180 */
         SHAPE  *hardCoff90,  *hardCoff180;     /* hard carbon 90, 180 pulses off resonance */
         SHAPE  *zeroCoffC90, *zeroCoffC180;    /* hard pulse, zero rotation Coff, 90, 180 on resonance; NB only if |Coff-Car| > MINCARDIFF */
         SHAPE  *zeroCCoff90, *zeroCCoff180;    /* hard pulse, zero rotation carrier, 90, 180 off resonance; NB only if |Coff-Car| > MINCARDIFF */
         SHAPE  *sincCoffC90, *sincCoffC180;    /* sinc zero at Coff, 90, 180 pulses on resonance; NB only if |Coff-Car| > MINCARDIFF  */
         SHAPE  *sincCCoff90, *sincCCoff180;    /* sinc zero at carrier, 90, 180 pulses off resonance; NB only if |Coff-Car| > MINCARDIFF  */
         SHAPE  *g3C180;                        /* g3 180 pulse on-resonance */
         SHAPE  *g3Coff180;                     /* g3 180 pulse at Coff */

        /* nitrogen */
         double Nfrq, pwNref, NpwrRef, compN;   /* nitrogen freq, ref pw, ref power, compN */
         SHAPE  *hardN90, *hardN180;            /* hard nitrogen 90, 180 pulses on resonance */
 
        } TRIPLERESPULSES;


TRIPLERESPULSES *newTripleResPulses() 
/* initialize a new TRIPERESPULSES instance */
{TRIPLERESPULSES *tmp;

 if ((tmp = (TRIPLERESPULSES *) malloc(sizeof(TRIPLERESPULSES))) == NULL)
        {return NULL;
        }
 tmp->dirname        = "";
 tmp->prefix         = "";
        
 tmp->Hfrq           = 0.0;
 tmp->pwHref         = 0.0;
 tmp->HpwrRef        = 0.0;
 tmp->compH          = 0.0;
 tmp->hardH90        = NULL;
 tmp->hardH180       = NULL;
 tmp->sincH2OH90     = NULL;

 tmp->Cfrq           = 0.0;
 tmp->pwCref         = 0.0;
 tmp->refpwrC        = 0.0;
 tmp->compC          = 0.0;
 tmp->Ccar           = 0.0;
 tmp->CoffCar        = 0.0;
 
 tmp->hardC90        = NULL;
 tmp->hardC180       = NULL;
 tmp->shaka6C180     = NULL;
 tmp->comp244C180    = NULL;
 tmp->hardCoff90     = NULL;
 tmp->hardCoff180    = NULL;
 
 tmp->zeroCoffC90    = NULL;
 tmp->zeroCoffC180   = NULL;
 tmp->zeroCCoff90    = NULL;
 tmp->zeroCCoff180   = NULL;
 
 tmp->sincCoffC90    = NULL;
 tmp->sincCoffC180   = NULL;
 tmp->sincCCoff90    = NULL;
 tmp->sincCCoff180   = NULL;
 
 tmp->g3C180         = NULL;
 tmp->g3Coff180      = NULL;

 tmp->Nfrq           = 0.0;
 tmp->pwNref         = 0.0;
 tmp->NpwrRef        = 0.0;
 tmp->compN          = 0.0;
 tmp->hardN90        = NULL;
 tmp->hardN180       = NULL;

 return tmp;
}

TRIPLERESPULSES *TripleResPulses(
/* initialize the TRIPERESPULSES instance
   create directory dirname
   store parameters
   initialize shapes 
   
   9 April:: remember to call with dirname AND prefix strings
*/
        dirname, prefix, 
        Hdevice, pwHref, HpwrRef, compH,
        Cdevice, pwCref, refpwrC, compC, Ccar, CoffCar,
        Ndevice, pwNref, NpwrRef, compN,
        verbose
)
char *dirname, *prefix;
int  Hdevice; double pwHref, HpwrRef, compH;
int  Cdevice; double pwCref, refpwrC, compC, Ccar, CoffCar;
int  Ndevice; double pwNref, NpwrRef, compN;
int  verbose;
 
{
 TRIPLERESPULSES *tmp;

 char    name[MAXLEN], *home;
 char   *getenv();

 double Hfrq = sfrq,
 	Cfrq = dfrq,
 	Nfrq = dfrq2; /* TODO; make these device dependent */
 double offresFreq = (CoffCar-Ccar) * Cfrq;

/* make a root directory for the shapes, based on dirname */
 if ((home = getenv("HOME")) == NULL)
        {fprintf(ERROR_FILE,"ERROR - TripleResPulses: retrieving home\n");
         return NULL;
        }        
 sprintf(name,"%s/vnmrsys/shapelib/%s", home, dirname);
 mkdir(name, S_IRWXU );

 if ((tmp = newTripleResPulses()) == NULL) {return NULL;}
 
 tmp->prefix  = strdup(prefix);
 tmp->dirname = strdup(dirname);

/* generating shapes */

/***** proton *****/
 tmp->Hfrq    = Hfrq;
 tmp->pwHref  = pwHref;
 tmp->HpwrRef = HpwrRef;
 tmp->compH   = compH;

/* proton hard pulses */
 sprintf(name, "%s/%s%s", dirname, prefix, "hardH90");
 if ((tmp->hardH90 = hard_shape( pwHref, 0.0, 0.0, pwHref, 1.0, name, verbose )) == NULL)  
 	{fprintf(ERROR_FILE,"Error generating pulse %s\n", name); return NULL;}
 tmp->hardH90->pwr    = HpwrRef;
 tmp->hardH90->device = Hdevice;
 sprintf(name, "%s/%s%s", dirname, prefix, "hardH180");
 if ((tmp->hardH180 = hard_shape( pwHref*2.0, 0.0, 0.0, pwHref*2.0, 1.0, name, verbose )) == NULL)  
 	{fprintf(ERROR_FILE,"Error generating pulse %s\n", name); return NULL;}
 tmp->hardH180->pwr    = HpwrRef;
 tmp->hardH180->device = Hdevice;
 /* Proton sinc for H2O */
 sprintf(name, "%s/%s%s", dirname, prefix, "sincH2OH90");
 if ((tmp->sincH2OH90 = sinc_shape( us(1700)*500.0/Hfrq, 0.0, 0.0, pwHref, compH, name, verbose )) == NULL) {return NULL;}
 tmp->sincH2OH90->pwr    = HpwrRef;
 tmp->sincH2OH90->device = Hdevice;

/***** carbon *****/
 tmp->Cfrq    = Cfrq;
 tmp->pwCref  = pwCref;
 tmp->refpwrC = refpwrC;
 tmp->compC   = compC;
 tmp->Ccar    = Ccar;
 tmp->CoffCar = CoffCar;
 
 /* Carbon hard pulses */
 sprintf(name, "%s/%s%s", dirname, prefix, "hardC90");
 if ((tmp->hardC90 = hard_shape( pwCref, 0.0, 0.0, pwCref, 1.0, name, verbose )) == NULL) 
 	{fprintf(ERROR_FILE,"Error generating pulse %s\n", name); return NULL;}
 tmp->hardC90->pwr = refpwrC;
 tmp->hardC90->device = Cdevice;
 sprintf(name, "%s/%s%s", dirname, prefix, "hardC180");
 if ((tmp->hardC180 = hard_shape( pwCref*2.0, 0.0, 0.0, pwCref*2.0, 1.0, name, verbose )) == NULL)  
 	{fprintf(ERROR_FILE,"Error generating pulse %s\n", name); return NULL;}
 tmp->hardC180->pwr = refpwrC;
 tmp->hardC180->device = Cdevice;
 sprintf(name, "%s/%s%s", dirname, prefix, "shaka6C180");
 if ((tmp->shaka6C180 = shaka6_shape( pwCref, refpwrC, name, verbose )) == NULL)  
 	{fprintf(ERROR_FILE,"Error generating pulse %s\n", name); return NULL;}
 tmp->shaka6C180->pwr = refpwrC;
 tmp->shaka6C180->device = Cdevice;
  sprintf(name, "%s/%s%s", dirname, prefix, "comp244C180");
 if ((tmp->comp244C180 = comp244_shape( pwCref, refpwrC, name, verbose )) == NULL)  
 	{fprintf(ERROR_FILE,"Error generating pulse %s\n", name); return NULL;}
 tmp->comp244C180->pwr = refpwrC;
 tmp->comp244C180->device = Cdevice;
sprintf(name, "%s/%s%s", dirname, prefix, "hardCoff90");
 if ((tmp->hardCoff90 = hard_shape( pwCref, offresFreq, 0.0, pwCref, 1.0, name, verbose )) == NULL)  
 	{fprintf(ERROR_FILE,"Error generating pulse %s\n", name); return NULL;}
 tmp->hardCoff90->pwr = refpwrC;
 tmp->hardCoff90->device = Cdevice;
 sprintf(name, "%s/%s%s", dirname, prefix, "hardCoff180");
 if ((tmp->hardCoff180 = hard_shape( pwCref*2.0, offresFreq, 0.0, pwCref*2.0, 1.0, name, verbose )) == NULL)  
 	{fprintf(ERROR_FILE,"Error generating pulse %s\n", name); return NULL;}
 tmp->hardCoff180->pwr = refpwrC;
 tmp->hardCoff180->device = Cdevice;
 
 /* carbon hard with zero rotation at an offset Co, Ca, or Cali */
#define MINCARDIFF 20.0
 if (fabs(CoffCar-Ccar) > MINCARDIFF) {
        sprintf(name, "%s/%s%s", dirname, prefix, "zeroCoffC180");
        if ((tmp->zeroCoffC180 = hard_shape( (double) fabs(2.0*sqrt(3.0)/(4.0*offresFreq)), 0.0, 0.0, pwCref*2.0, compC, name, verbose )) == NULL)  
 		{fprintf(ERROR_FILE,"Error generating pulse %s\n", name); return NULL;}
        tmp->zeroCoffC180->pwr = refpwrC;
        tmp->zeroCoffC180->device = Cdevice;
        sprintf(name, "%s/%s%s", dirname, prefix, "zeroCoffC90");
        if ((tmp->zeroCoffC90 = hard_shape( (double) fabs(sqrt(15.0)/(4.0*offresFreq)), 0.0, 0.0, pwCref, compC, name, verbose )) == NULL)  
 		{fprintf(ERROR_FILE,"Error generating pulse %s\n", name); return NULL;}
        tmp->zeroCoffC90->pwr = refpwrC;
        tmp->zeroCoffC90->device = Cdevice;
        sprintf(name, "%s/%s%s", dirname, prefix, "zeroCCoff180");
        if ((tmp->zeroCCoff180 = hard_shape( (double) fabs(2.0*sqrt(3.0)/(4.0*offresFreq)), offresFreq, 0.0, pwCref*2.0, compC, name, verbose )) == NULL)  
 		{fprintf(ERROR_FILE,"Error generating pulse %s\n", name); return NULL;}
        tmp->zeroCCoff180->pwr = refpwrC;
        tmp->zeroCCoff180->device = Cdevice;
        sprintf(name, "%s/%s%s", dirname, prefix, "zeroCCoff90");
        if ((tmp->zeroCCoff90 = hard_shape( (double) fabs(sqrt(15.0)/(4.0*offresFreq)), offresFreq, 0.0, pwCref, compC, name, verbose )) == NULL)  
 		{fprintf(ERROR_FILE,"Error generating pulse %s\n", name); return NULL;}
        tmp->zeroCCoff90->pwr = refpwrC;
        tmp->zeroCCoff90->device = Cdevice;
} else {
        fprintf(ERROR_FILE, "Warning: skipping zero-exicitation pulses; carriers too close (CoffCar-Ccar = %.1f)\n", CoffCar-Ccar);
}

/* sinc's with zero's at specific frequencies */
 if (fabs(CoffCar-Ccar) > MINCARDIFF) {
        sprintf(name, "%s/%s%s", dirname, prefix, "sincCoffC180");
        if ((tmp->sincCoffC180 = sinc_shape( (double) fabs( us(100)*14300.0/offresFreq ), 0.0, 0.0, pwCref*2.0, compC, name, verbose )) == NULL)  
 		{fprintf(ERROR_FILE,"Error generating pulse %s\n", name); return NULL;}
        tmp->sincCoffC180->pwr = refpwrC;
        tmp->sincCoffC180->device = Cdevice;
        sprintf(name, "%s/%s%s", dirname, prefix, "sincCoffC90");
        if ((tmp->sincCoffC90 = sinc_shape( (double) fabs( us(100)*15800.0/offresFreq), 0.0, 0.0, pwCref, compC, name, verbose )) == NULL)  
 		{fprintf(ERROR_FILE,"Error generating pulse %s\n", name); return NULL;}
        tmp->sincCoffC90->pwr = refpwrC;
        tmp->sincCoffC90->device = Cdevice;
        sprintf(name, "%s/%s%s", dirname, prefix, "sincCCoff180");
        if ((tmp->sincCCoff180 = sinc_shape( (double) fabs( us(100)*14300.0/offresFreq ), offresFreq, 0.0, pwCref*2.0, compC, name, verbose )) == NULL)  
 		{fprintf(ERROR_FILE,"Error generating pulse %s\n", name); return NULL;}
        tmp->sincCCoff180->pwr = refpwrC;
        tmp->sincCCoff180->device = Cdevice;
        sprintf(name, "%s/%s%s", dirname, prefix, "sincCCoff90");
        if ((tmp->sincCCoff90 = sinc_shape( (double) fabs( us(100)*15800.0/offresFreq), offresFreq, 0.0, pwCref, compC, name, verbose )) == NULL)  
 		{fprintf(ERROR_FILE,"Error generating pulse %s\n", name); return NULL;}
        tmp->sincCCoff90->pwr = refpwrC;
        tmp->sincCCoff90->device = Cdevice;
} else {
        fprintf(ERROR_FILE, "Warning: skipping sinc pulses; carriers too close (CoffCar-Ccar = %.1f)\n", CoffCar-Ccar);
}
  
 /* carbon g3's */
 sprintf(name, "%s/%s%s", dirname, prefix, "g3C180");
 if ((tmp->g3C180 = Pbox_shape("g3", 10000.0*sfrq/800.0, 0.0, pwCref, refpwrC,"",name, verbose )) == NULL)  
 	{fprintf(ERROR_FILE,"Error generating pulse %s\n", name); return NULL;}
 tmp->g3C180->device = Cdevice;
 sprintf(name, "%s/%s%s", dirname, prefix, "g3Coff180");
 if ((tmp->g3Coff180 = Pbox_shape("g3", 10000.0*sfrq/800.0, offresFreq, pwCref, refpwrC,"",name, verbose )) == NULL)  
 	{fprintf(ERROR_FILE,"Error generating pulse %s\n", name); return NULL;}
 tmp->g3Coff180->device = Cdevice;

/***** nitrogen *****/
 tmp->Nfrq    = Nfrq;
 tmp->pwNref  = pwNref;
 tmp->NpwrRef = NpwrRef;
 tmp->compN   = compN;
 
 sprintf(name, "%s/%s%s", dirname, prefix, "hardN90");
 if ((tmp->hardN90 = hard_shape( pwNref, 0.0, 0.0, pwNref, 1.0, name, verbose )) == NULL)  {return NULL;}
 tmp->hardN90->pwr = NpwrRef;
 tmp->hardN90->device = Ndevice;
 sprintf(name, "%s/%s%s", dirname, prefix, "hardN180");
 if ((tmp->hardN180 = hard_shape( 2.0*pwNref, 0.0, 0.0, pwNref*2.0, 1.0, name, verbose )) == NULL)  {return NULL;}
 tmp->hardN180->pwr = NpwrRef;
 tmp->hardN180->device = Ndevice;
 
 return tmp;
}


/**********************************************************/
/*********** Probe & RF-channel dependent stuff ***********/
/**********************************************************/

#define TDEV            10
#define XDEV            11
#define YDEV            12
#define ZDEV            13
#define MINDEVICE       10
#define MAXDEVICE       13


int setProbe()
/* set current probe */
{double Id;
 int iId;
 
 if (probe == NULL) {
/* Patch to disable probing probeId */
#ifndef DISABLE_PROBE_CHECK
        /* get probeId from global parameter tree */
        if (getparm("probeId","real",GLOBAL, (char *) &Id, 1))
                {fprintf(ERROR_FILE, "ERROR: get_probeId, parameter undefined\n");
                 exit(1);
                }
#else
        Id  = 0.0;
#endif
        /* check the Id number for validity */
        iId = (int) Id;
        if (iId < 0 || iId >= MAX_PROBES)
                {fprintf(ERROR_FILE, "ERROR: get_probeId, invalid probeId (%d)\n", iId);
                 exit(1);
                }
        
         probe = probes  + iId;
         fprintf(ERROR_FILE,"Current probe: %s (id %d)\n", probe->name, probe->id);
        }
 return NOERROR;
}

#define POWER(dev, val, t) \
{setProbe(); setpower( dev, val, probe->max_power[dev-MINDEVICE] ); delay(t - POWER_DELAY - PWRF_DELAY);}

#define PRESAT_POWER(dev, val, t) \
{setProbe(); setpower( dev, val, probe->max_satpower[dev-MINDEVICE] ); delay(t - POWER_DELAY - PWRF_DELAY);}

#define DECOUPLE_POWER(dev, val, t) \
{setProbe(); setpower( dev, val, probe->max_decpower[dev-MINDEVICE] ); delay(t - POWER_DELAY - PWRF_DELAY);}

#define LN2 0.693147

void setpower(dev, value, maxval)
int dev; double value, maxval;
/* routine can also set non-integer dB values using fine attenuator */
{int    ival;
 double val, remain, pwrl;
 double exp();

 ival = (int) value;
 val  = (double) ival; 
 remain = value - val;

/* 
printf(">>%8.1f %8d %8.1f %8.1f\n", value, ival, val, remain);
*/

 if (remain == 0.0)
        /* set fine power to full value since it is an 'integer' number of dB's */
        {pwrl = 4095.0;
        }
 else
        /* ok, there is a remainder, set the 'integer' dB value one unit higher */
        {val += 1.0;
         remain = 1.0 - remain;         
        /* convert to 0, 4095 range of linear attenuator, dividing by 2 equals 6 dB */
        /* NB pow() function gives compile errors */
        
/*      
         remain = exp( (double) (remain / 6.0 * LN2)) - 1.0;
         pwrl = 4095.0 - remain * 2048.0;
*/
         pwrl = 4095.0/(exp( (double) (remain / 6.0 * LN2)));
         
        /* round value */
         pwrl = (double) ((int) (pwrl + 0.5));
        }

/*
printf(">>%8.1f %8.1f\n", val, pwrl);
*/

 switch (dev) { 
        case TDEV: {if (val > maxval) 
                        {fprintf(ERROR_FILE, "ERROR: setpower, too much transmitter RF power (%.1f, Max. %.1f)\n", 
                         (float) value, (float) maxval); 
                         exit(1); 
                        } 
                    obspower( val ); 
                    obspwrf( pwrl );
                    break; 
                   } 
        case XDEV: {if (val > maxval) 
                        {fprintf(ERROR_FILE, "ERROR: setpower, too much X-channel (Dec) RF power (%.1f, Max. %.1f)\n", 
                         (float) value, (float) maxval); 
                         exit(1); 
                        } 
                    decpower( val );
                    decpwrf( pwrl ); 
                    break; 
                   } 
        case YDEV: {if (val > maxval) 
                        {fprintf(ERROR_FILE, "ERROR: setpower, too much Y-channel (Dec2) RF power (%.1f, Max. %.1f)\n", 
                         (float) value, (float) maxval); 
                         exit(1); 
                        } 
                    dec2power( val );
                    dec2pwrf( pwrl ); 
                    break; 
                   } 
        case ZDEV: {if (val > maxval) 
                        {fprintf(ERROR_FILE, "ERROR: setpower, too much Z-channel (Dec3) RF power (%.1f, Max. %.1f)\n", 
                         (float) value, (float) maxval); 
                         exit(1); 
                        } 
                    dec3power( val );
                    dec3pwrf( pwrl ); 
                    break; 
                   } 
        default: {fprintf(ERROR_FILE, "ERROR: setpower, invalid device\n"); exit(1);} 
        }  
  }
  
int NTsetGradient( axis, level )
/* set gradient using -100-100% scale */
int axis;
double level;
{
 setProbe();
 switch (axis) {
         case PROBE_GRADX: { rgradient('x', level*probe->grad_factor[PROBE_GRADX]);}
         case PROBE_GRADY: { rgradient('y', level*probe->grad_factor[PROBE_GRADY]);}
         case PROBE_GRADZ: { rgradient('z', level*probe->grad_factor[PROBE_GRADZ]);} 
        }
}

#define FINE_POWER(dev, val, t) \
 {switch (dev) { \
        case TDEV: {obspwrf( val ); break;} \
        case XDEV: {decpwrf( val ); break;} \
        case YDEV: {dec2pwrf( val ); break;} \
        case ZDEV: {dec3pwrf( val ); break;} \
        default: {fprintf(ERROR_FILE, "ERROR: FINE_POWER invalid device\n"); exit(1);} \
        } \
   delay(t - PWRF_DELAY); \
  }
 
#define OFFSET(dev, val, t) \
 {switch (dev) { \
        case TDEV: {obsoffset( val ); break;} \
        case XDEV: {decoffset( val ); break;} \
        case YDEV: {dec2offset( val ); break;} \
        case ZDEV: {dec3offset( val ); break;} \
        default: {fprintf(ERROR_FILE, "ERROR: OFFSET, invalid device\n"); exit(1);} \
        } \
   delay(t - OFFSET_DELAY); \
  }
 

#define STEPSIZE(dev, val) \
 {switch (dev) { \
        case TDEV: {obsstepsize( val ); break;} \
        case XDEV: {decstepsize( val ); break;} \
        case YDEV: {dec2stepsize( val ); break;} \
        case ZDEV: {dec3stepsize( val ); break;} \
        default: {fprintf(ERROR_FILE, "ERROR: STEPSIZE, invalid device\n"); exit(1);} \
        } \
 }

#define SAPHASE(dev, val, t) \
 {switch (dev) { \
        case TDEV: {xmtrphase( val ); break;} \
        case XDEV: {dcplrphase( val ); break;} \
        case YDEV: {dcplr2phase( val ); break;} \
        case ZDEV: {dcplr3phase( val ); break;} \
        default: {fprintf(ERROR_FILE, "ERROR: SAPHASE, invalid device\n"); exit(1);} \
        } \
   delay(t - SAPS_DELAY); \
 }
 
#define TRON(dev) \
 {switch (dev) { \
        case TDEV: {obsunblank(); xmtron(); break;} \
        case XDEV: {decunblank(); decon(); break;} \
        case YDEV: {dec2unblank(); dec2on(); break;} \
        case ZDEV: {dec3unblank(); dec3on(); break;} \
        default: {fprintf(ERROR_FILE, "ERROR: TRON, invalid device\n"); exit(1);} \
        } \
 }
 
#define TRDON(dev,rof1) \
 {switch (dev) { \
        case TDEV: {obsunblank(); delay(rof1); xmtron(); break;} \
        case XDEV: {decunblank(); delay(rof1); decon(); break;} \
        case YDEV: {dec2unblank(); delay(rof1); dec2on(); break;} \
        case ZDEV: {dec3unblank(); delay(rof1); dec3on(); break;} \
        default: {fprintf(ERROR_FILE, "ERROR: TRDON, invalid device\n"); exit(1);} \
        } \
 }
 
#define TROFF(dev) \
 {switch (dev) { \
        case TDEV: {xmtroff(); obsblank(); break;} \
        case XDEV: {decoff();  decblank(); break;} \
        case YDEV: {dec2off(); dec2blank(); break;} \
        case ZDEV: {dec3off(); dec3blank(); break;} \
        default: {fprintf(ERROR_FILE, "ERROR: TROFF, invalid device\n"); exit(1);} \
        } \
 }
 
#define TRDOFF(dev,rof2) \
 {switch (dev) { \
        case TDEV: {xmtroff(); delay(rof2); obsblank(); break;} \
        case XDEV: {decoff(); delay(rof2);  decblank(); break;} \
        case YDEV: {dec2off(); delay(rof2); dec2blank(); break;} \
        case ZDEV: {dec3off(); delay(rof2); dec3blank(); break;} \
        default: {fprintf(ERROR_FILE, "ERROR: TRDOFF, invalid device\n"); exit(1);} \
        } \
 }

#define TRMTRON(dev) \
 {switch (dev) { \
        case TDEV: { xmtron(); break;} \
        case XDEV: { decon(); break;} \
        case YDEV: { dec2on(); break;} \
        case ZDEV: { dec3on(); break;} \
        default: {fprintf(ERROR_FILE, "ERROR: TRMTRON, invalid device\n"); exit(1);} \
        } \
 }
 
#define TRMTROFF(dev) \
 {switch (dev) { \
        case TDEV: {xmtroff(); break;} \
        case XDEV: {decoff();  break;} \
        case YDEV: {dec2off(); break;} \
        case ZDEV: {dec3off(); break;} \
        default: {fprintf(ERROR_FILE, "ERROR: TRMTROFF, invalid device\n"); exit(1);} \
        } \
 }

#define BLANK(dev) \
 {switch (dev) { \
        case TDEV: {obsblank(); break;} \
        case XDEV: {decblank(); break;} \
        case YDEV: {dec2blank(); break;} \
        case ZDEV: {dec3blank(); break;} \
        default: {fprintf(ERROR_FILE, "ERROR: BLANK, invalid device\n"); exit(1);} \
        } \
 }
 
#define UNBLANK(dev) \
 {switch (dev) { \
        case TDEV: {obsunblank(); break;} \
        case XDEV: {decunblank(); break;} \
        case YDEV: {dec2unblank(); break;} \
        case ZDEV: {dec3unblank(); break;} \
        default: {fprintf(ERROR_FILE, "ERROR: UNBLANK, invalid device\n"); exit(1);} \
        } \
 }
 
#define PHASE(dev, ph) \
 {switch (dev) { \
        case TDEV: {txphase(ph); break;} \
        case XDEV: {decphase(ph); break;} \
        case YDEV: {dec2phase(ph); break;} \
        case ZDEV: {dec3phase(ph); break;} \
        default: {fprintf(ERROR_FILE, "ERROR: PHASE, invalid device\n"); exit(1);} \
        } \
 }
 
#define PRGON(dev, pattern, pw, dres, dly) \
 {switch (dev) { \
        case TDEV: {obsprgon( pattern, pw, dres ); if (dly>0) delay(dly-PRG_START_DELAY); break;} \
        case XDEV: {decprgon( pattern, pw, dres ); if (dly>0) delay(dly-PRG_START_DELAY); break;} \
        case YDEV: {dec2prgon( pattern, pw, dres ); if (dly>0) delay(dly-PRG_START_DELAY); break;} \
        case ZDEV: {dec3prgon( pattern, pw, dres ); if (dly>0) delay(dly-PRG_START_DELAY); break;} \
        default: {fprintf(ERROR_FILE, "ERROR: PHASE, invalid device\n"); exit(1);} \
        } \
 }
 
#define PRGOFF(dev, dly) \
 {switch (dev) { \
        case TDEV: {obsprgoff();  if (dly>0) delay(dly-PRG_STOP_DELAY); break;} \
        case XDEV: {decprgoff();  if (dly>0) delay(dly-PRG_STOP_DELAY); break;} \
        case YDEV: {dec2prgoff(); if (dly>0) delay(dly-PRG_STOP_DELAY); break;} \
        case ZDEV: {dec3prgoff(); if (dly>0) delay(dly-PRG_STOP_DELAY); break;} \
        default: {fprintf(ERROR_FILE, "ERROR: PHASE, invalid device\n"); exit(1);} \
        } \
 }
 
typedef struct Cpd {
        int     device;         /* RF device/channel */
	char   *pattern;	/* cpd decoupling waveform pattern */
        double  pw;             /* 90 degree pulse */
        double  pwr;  		/* rf power for decoupling */
        double  dres;  		/* resolution of the cpd pattern */
	int     ison;		/* defines if cpd is on */
        } CPD;

CPD *Cpd( int device, char *pattern, double dres, double pw, double pwr )
/* Initialize a CPD instance, rturn NULL pointer on error */
{CPD *tmp;
 if ((tmp = (CPD *) malloc(sizeof(CPD))) == NULL)
        {return NULL;
        }
 tmp->device  = device;
 tmp->pattern = pattern;
 tmp->pw      = pw;
 tmp->dres    = dres;
 tmp->pwr     = pwr;
 tmp->ison    = FALSE;
 return tmp; 
}

/* 
int cpdon(int dev, char *pat, double p90, double res, codeint ph, double dly) 
 {switch (dev) { 
        case TDEV: {txphase(ph);   obsunblank();  delay(dly-PRG_START_DELAY); xmtron(); obsprgon( pat, p90, res);  break;} 
        case XDEV: {decphase(ph);  decunblank();  delay(dly-PRG_START_DELAY); decon();  decprgon( pat, p90, res);  break;} 
        case YDEV: {dec2phase(ph); dec2unblank(); delay(dly-PRG_START_DELAY); dec2on(); dec2prgon( pat, p90, res); break;} 
        case ZDEV: {dec3phase(ph); dec3unblank(); delay(dly-PRG_START_DELAY); dec3on(); dec3prgon( pat, p90, res); break;} 
        default: {fprintf(ERROR_FILE, "ERROR: cpdon, invalid device\n"); exit(1);} 
        } 
 }
*/
int cpdon( CPD *cpd, codeint ph, double dly ) 
{PHASE(cpd->device, ph);
 UNBLANK( cpd->device );
 DECOUPLE_POWER( cpd->device, cpd->pwr, dly-PRG_START_DELAY );
 TRMTRON( cpd->device );
 PRGON( cpd->device, cpd->pattern, cpd->pw, cpd->dres, 0.0 );  /* 0.0 results in No delay, time is compensated for above */
 cpd->ison = TRUE;
 return 0;
}

int cpdoff( CPD *cpd, double dly )
{TRMTROFF( cpd->device );
 PRGOFF( cpd->device, dly-us(2) );
 BLANK( cpd->device );
 cpd->ison = FALSE;
 return 0;
}

#define CPDon(cpd,ph,dly) \
{PHASE(cpd->device, ph);\
 UNBLANK( cpd->device );\
 DECOUPLE_POWER( cpd->device, cpd->pwr, us(2) );\
 TRMTRON( cpd->device );\
 PRGON( cpd->device, cpd->pattern, cpd->pw, cpd->dres, dly-us(2) );  \
}\


#define CPDoff(cpd,dly) TRCoff(cpd->device,dly)


#define TRCoff(dev, dly) \
 {switch (dev) { \
        case TDEV: {xmtroff(); obsprgoff();  delay(dly-PRG_STOP_DELAY); obsblank();  break;} \
        case XDEV: {decoff();  decprgoff();  delay(dly-PRG_STOP_DELAY); decblank();  break;} \
        case YDEV: {dec2off(); dec2prgoff(); delay(dly-PRG_STOP_DELAY); dec2blank(); break;} \
        case ZDEV: {dec3off(); dec3prgoff(); delay(dly-PRG_STOP_DELAY); dec3blank(); break;} \
        default: {fprintf(ERROR_FILE, "ERROR: TRCoff, invalid device\n"); exit(1);} \
        } \
 }
 

  
/* ============ */

#define STPULSE( shape, ph, rof1, rof2 ) \
{double roftmp=rof1; \
 txphase( ph ); \
 if (shape->pwr > SHAPE_PWR_NOT_DEFINED) { \
        obspower( shape->pwr );  \
        roftmp -= POWER_DELAY; \
        } \
 obspwrf( shape->rf ); \
 shaped_pulse( shape->name, shape->pw, ph, roftmp-PWRF_DELAY-WFG_START_DELAY, MAX(0.0,rof2-WFG_STOP_DELAY) ); \
}

#define SXPULSE( shape, ph, rof1, rof2 ) \
{double roftmp=rof1; \
 decphase( ph ); \
 if (shape->pwr > SHAPE_PWR_NOT_DEFINED) { \
        decpower( shape->pwr );  \
        roftmp -= POWER_DELAY; \
        } \
 decpwrf( shape->rf ); \
 decshaped_pulse( shape->name, shape->pw, ph, roftmp-PWRF_DELAY-WFG_START_DELAY, MAX(0.0,rof2-WFG_STOP_DELAY) ); \
}
 
#define SYPULSE( shape, ph, rof1, rof2 ) \
{double roftmp=rof1; \
 dec2phase( ph ); \
 if (shape->pwr > SHAPE_PWR_NOT_DEFINED) { \
        dec2power( shape->pwr );  \
        roftmp -= POWER_DELAY; \
        } \
 dec2pwrf( shape->rf ); \
 dec2shaped_pulse( shape->name, shape->pw, ph, roftmp-PWRF_DELAY-WFG_START_DELAY, MAX(0.0,rof2-WFG_STOP_DELAY) ); \
}

#define SZPULSE( shape, ph, rof1, rof2 ) \
{double roftmp=rof1; \
 dec3phase( ph ); \
 if (shape->pwr > SHAPE_PWR_NOT_DEFINED) { \
        dec3power( shape->pwr );  \
        roftmp -= POWER_DELAY; \
        } \
 dec3pwrf( shape->rf ); \
 dec3shaped_pulse( shape->name, shape->pw, ph, roftmp-PWRF_DELAY-WFG_START_DELAY, MAX(0.0,rof2-WFG_STOP_DELAY) ); \
}

double s3pulse( SHAPE *shapeT, SHAPE *shapeX, SHAPE *shapeY, codeint phT, codeint phX, codeint phY, 
                double rof1, double rof2 
	       )
/* Execute a shapedpulse on arbitrary channels 
   return the MAX of the pulsewidth
*/
{double  roftmp=rof1, rof2tmp=rof2; 
 double  pw1,    pw2,    pw3; 
 char   *name1, *name2, *name3; 
 
 if (shapeT!=NULL) { 
 	pw1 = shapeT->pw; name1 = shapeT->name; 
 	txphase( phT ); 
 	if (shapeT->pwr > SHAPE_PWR_NOT_DEFINED) { 
        	obspower( shapeT->pwr );  roftmp -= POWER_DELAY; 
        	} 
 	obspwrf( shapeT->rf ); roftmp -= PWRF_DELAY; 
 } else {pw1 = 0.0; name1 = ""; } 
 
 if (shapeX!=NULL) { 
 	pw2 = shapeX->pw; name2 = shapeX->name; 
 	decphase( phX ); 
 	if (shapeX->pwr > SHAPE_PWR_NOT_DEFINED) { 
        	decpower( shapeX->pwr );  roftmp -= POWER_DELAY; 
        	} 
 	decpwrf( shapeX->rf ); roftmp -= PWRF_DELAY; 
 } else {pw2 = 0.0; name2 = ""; } 

 if (shapeY!=NULL) { 
 	pw3 = shapeY->pw; name3 = shapeY->name; 
 	dec2phase( phY ); 
 	if (shapeY->pwr > SHAPE_PWR_NOT_DEFINED) { 
        	dec2power( shapeY->pwr );  roftmp -= POWER_DELAY; 
        	} 
 	dec2pwrf( shapeY->rf ); roftmp -= PWRF_DELAY; 
 } else {pw3 = 0.0; name3 = ""; } 
 
 roftmp -= WFG3_START_DELAY;
 rof2tmp -= WFG3_STOP_DELAY;
 
 if (roftmp < 0.0) {
 	fprintf(ERROR_FILE,"WARNING s3pulse: FID %d, roftmp = %.7f\n", 
		           (int) ix, roftmp
	       );
 }
 if (rof2tmp < 0.0) {
 	fprintf(ERROR_FILE,"WARNING s3pulse: FID %d, rof2tmp = %.7f\n", 
		           (int) ix, rof2tmp
	       );
 }
 
 sim3shaped_pulse( name1, name2, name3, 
		   pw1,   pw2,   pw3, 
		   phT,   phX,   phY, 
		   MAX(roftmp, 0.0), MAX(0.0,rof2tmp) 
		  ); 
 return MAX(MAX(pw1,pw2),pw3);
}
  



/*************************************************************************/
/* Event code                                                            */
/*************************************************************************/

#define QPHASE(dev, ph) \
 {switch (dev) { \
        case TDEV: {txphase( ph ); break;} \
        case XDEV: {decphase( ph ); break;} \
        case YDEV: {dec2phase( ph ); break;} \
        case ZDEV: {dec3phase( ph ); break;} \
        default: {fprintf(ERROR_FILE, "ERROR: SET_PHASE, invalid device\n"); exit(1);} \
        } \
 }


/*************************************************************************/


#define DELAY_EVENT     0
#define BLANK_EVENT     1
#define UNBLANK_EVENT   2
#define XMTRON_EVENT    3
#define XMTROFF_EVENT   4
#define SET_PHASE_EVENT 5
#define NUM_EVENT       6

char *eventList[NUM_EVENT] = {
        "delay",
        "blank",
        "unblank",
        "xmtr on",
        "xmtr off",
        "set phase"
};

typedef struct Event {
         int    chan;
         int    type;
         double duration;
         codeint phase;
        } EVENT;
        
#define MIN_DURATION 25e-9      /* minimum duration */
        
                
/*************************************************************************/

int print_event( EVENT *e )
{
 fprintf( ERROR_FILE, "Event chan %2d, %-8s (%1d), duration %11.8f\n", 
                       e->chan - MINDEVICE, eventList[e->type], e->type, e->duration );
 return NOERROR;
}

int set_event( EVENT *e, int chan, int type, double duration )
{
 e->chan = chan;
 e->type = type;
 e->duration = duration;
 return NOERROR ;
}

int execute_event( EVENT *e )
{
 switch (e->type)
        {case DELAY_EVENT:
         break;

         case BLANK_EVENT: BLANK( e->chan );
         break;

         case UNBLANK_EVENT: UNBLANK( e->chan );
         break;

         case XMTRON_EVENT: TRMTRON( e->chan );
         break;

         case XMTROFF_EVENT:TRMTROFF( e->chan );
         break;

         case SET_PHASE_EVENT: QPHASE( e->chan, e->phase );
         break;

         default: return ERROR;
         break;
        }
 return NOERROR;
}


void sim_dpd( int chan1, double d1a, double p1, codeint ph1, double d1b,
              int chan2, double d2a, double p2, codeint ph2, double d2b
            )
/* delay-pulse-delay on two channels simultaneously */

#define N_EVENT 8

{EVENT channel1[N_EVENT], channel2[N_EVENT];
 int ec1, ec2;

 if ( fabs( (d1a+p1+d1b) - (d2a+p2+d2b)) > MIN_DURATION )
        {fprintf( ERROR_FILE, "ERROR - SIM_DPD: unequal total durations, channel %d: %11.8f  channel %d: %11.8f \n",
                  chan1,  (d1a+p1+d1b), chan2, (d2a+p2+d2b)
                );
         exit( ERROR );
        }
        

/* store the delay and pulse durations */
 
 set_event( &(channel1[0]), chan1, SET_PHASE_EVENT, 0.0 ); channel1[0].phase = ph1;
 set_event( &(channel1[1]), chan1, DELAY_EVENT,     d1a );
 set_event( &(channel1[2]), chan1, UNBLANK_EVENT,   0.0 );
 set_event( &(channel1[3]), chan1, XMTRON_EVENT,    0.0 );
 set_event( &(channel1[4]), chan1, DELAY_EVENT,     p1 );
 set_event( &(channel1[5]), chan1, XMTROFF_EVENT,   0.0 );
 set_event( &(channel1[6]), chan1, BLANK_EVENT,     0.0 );
 set_event( &(channel1[7]), chan1, DELAY_EVENT,     d1b );


 set_event( &(channel2[0]), chan2, SET_PHASE_EVENT, 0.0 ); channel2[0].phase = ph2;
 set_event( &(channel2[1]), chan2, DELAY_EVENT,     d2a );
 set_event( &(channel2[2]), chan2, UNBLANK_EVENT,   0.0 );
 set_event( &(channel2[3]), chan2, XMTRON_EVENT,    0.0 );
 set_event( &(channel2[4]), chan2, DELAY_EVENT,     p2 );
 set_event( &(channel2[5]), chan2, XMTROFF_EVENT,   0.0 );
 set_event( &(channel2[6]), chan2, BLANK_EVENT,     0.0 );
 set_event( &(channel2[7]), chan2, DELAY_EVENT,     d2b );


/* execute this by continuously executing the shortes event */
 ec1 = 0; ec2 = 0;
 while (ec1 < N_EVENT && ec2 < N_EVENT )
        {
        
/*
print_event( &(channel1[ec1]) );        
print_event( &(channel2[ec2]) );        
*/      
        
         execute_event( &(channel1[ec1]) );
         execute_event( &(channel2[ec2]) );

         if (fabs(channel1[ec1].duration - channel2[ec2].duration) <= MIN_DURATION)
                {if (channel1[ec1].duration > 0) delay( channel1[ec1].duration );
                 ec1++;
                 ec2++;
                }
        
         else if (channel2[ec2].duration - channel1[ec1].duration > MIN_DURATION)
                {if (channel1[ec1].duration > 0) delay( channel1[ec1].duration );
                 channel2[ec2].type = DELAY_EVENT; /* true action has been taken, we only need to complete the duration */
                 channel2[ec2].duration -= channel1[ec1].duration;
                 ec1++;
                }
         else if (channel1[ec1].duration - channel2[ec2].duration > MIN_DURATION)
                {if (channel2[ec2].duration > 0) delay( channel2[ec2].duration );
                 channel1[ec1].type = DELAY_EVENT; /* true action has been taken, we only need to complete the duration */
                 channel1[ec1].duration -= channel2[ec2].duration;
                 ec2++;
                }
                
         else
                {fprintf( ERROR_FILE, "ERROR - SIM_DPD: this should not occur\n");
                 exit( ERROR );
                }
        }

 if (ec1 < N_EVENT )
        {fprintf( ERROR_FILE, "ERROR - SIM_DPD: Fid %d: unfinished event handeling (stack 1: %d)\n", (int) ix, ec1);
         print_event( &(channel1[ec1]) );
         exit( ERROR );
        }

 if (ec2 < N_EVENT )
        {fprintf( ERROR_FILE, "ERROR - SIM_DPD: Fid %d: unfinished event handeling (stack 2: %d)\n", (int) ix, ec2);
         print_event( &(channel2[ec2]) );
         exit( ERROR );
        }

}


/*************************************************************************/


