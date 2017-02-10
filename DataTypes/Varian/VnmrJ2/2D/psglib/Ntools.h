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
             
   07-03-06: Implemented DISABLE_PROBE_CHECK directive. Defining this parameter
             with '#define' before inclusion of Ntools.h, disable the probe check
             and values from 'probes[0]' (the generic one) are taken.
	     
   05-11-07: Adapted for redhad LINUX using ANSI statement and moving routine comments
             into the body (after the "{"), because DPS does not like it.

*******************************************************************************/

#include <stdio.h>
#include <math.h>
#include <group.h>

#define ANSI

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

#ifndef ANSI 
char *strdup();
#endif

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
{
/* character by character check of string1 against string2
   string2 can contain stars to indicate arbitrary values
   string1 is expanded with last character if shorter than string2
   Ex. check_string( dmm, "dmm", "ccc*");
*/

 char *s1=string1, *s2=string2;

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
        char   *waveType;       /* Shape wavefor type; e.g. g3, sinc etc. */
        } SHAPE;

#define SHAPE_RF                0
#define SHAPE_DEC               1
static char *shape_type[2] = {"RF","DEC"};

#define SHAPE_RESOLUTION        0.2e-6
#define SHAPE_MAXAMPLITUDE      1023.0
#define SHAPE_PWR_NOT_DEFINED   -99.99

#ifdef ANSI
int list_shape( SHAPE *shape, FILE *fp )
#else
int list_shape( shape , fp )
SHAPE *shape; FILE *fp;
#endif
{if (shape != NULL)
/* NB: line cannot be to long because pulsetool will crash */
        fprintf(fp,
                "%s (%s): n %d, pw %.1f us, off. %.1f Hz, mod. %.2f, fac. %.2f\n",
                shape->name, 
                shape->waveType,
                shape->nsteps, 
                shape->pw*1e6, 
                shape->offset,
                shape->mod,
                shape->ampfactor
                );
 return NOERROR;
}


SHAPE *new_shape()
{SHAPE *tmp;

 if ((tmp = (SHAPE *) malloc(sizeof(SHAPE))) == NULL)
        {return NULL;
        }
        
 tmp->type = SHAPE_RF;
 tmp->nsteps = 0;
 tmp->pw = 0.0;
 tmp->pwr = SHAPE_PWR_NOT_DEFINED;
 tmp->rf  = FINEPOWER_MAXVAL;
 tmp->B1max = 0.0;
 tmp->amp = NULL;
 tmp->phase = NULL;
 tmp->offset =0.0;
 tmp->mod = 0.0;
 tmp->name = NULL;
 tmp->waveType = NULL;

 return tmp;
}


#ifdef ANSI
SHAPE *write_shape( SHAPE *shape, char *shapename )
#else
SHAPE *write_shape( shape, shapename )
SHAPE *shape; char *shapename;
#endif
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

 fprintf(fout,"# NTools ");
 list_shape( shape, fout );

 for (i=0; i!= shape->nsteps; i++)
        {fprintf( fout, "%10.3f %10d.0    1.0 \n", shape->phase[i], (int) (shape->amp[i]+0.5));
        }

 fclose(fout);

 return shape;

}

#ifdef ANSI
SHAPE *convolute_shape( SHAPE *shape, double offset, double mod)
#else
SHAPE *convolute_shape( shape, offset, mod )
SHAPE *shape; double offset, mod;
#endif
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

#ifdef ANSI
double calcrf_shape( SHAPE *shape, double pwref, double comp )
#else
double calcrf_shape( shape, pwref, comp )
SHAPE *shape; double pwref, comp;
#endif
{
/* calculate rf power for fine attenuator
   pwref: reference pw; has to have the same flipangle as pw of the shape
   comp: amplifier compression factor
*/

 shape->rf = (comp*FINEPOWER_MAXVAL*pwref) / shape->pw * 1.0/shape->ampfactor;
 shape->rf = (int) (shape->rf + 0.5);
 return shape->rf;
}               


SHAPE *hard_shape(pw, offset, mod, pwref, comp, name, verbose)
double pw, offset, mod, pwref, comp; char *name; int verbose;
{
/*      pw:     pulse width
        offset: offset (in Hz) of the pulse
        mod:    zero point of phase modulation pattern; begin (mod=0.0) to end (mod=1.0)
        pwref:  reference pw for calculation of the fineattenuator value; has to have the same flipangle as pw of the shape
        comp:   amplifier compression factor
        name:   VNMR shapelib name
        verbose: verbose flag
*/

 SHAPE *shape; int i;

 if ((shape = new_shape() ) == NULL)
        {fprintf(ERROR_FILE, "ERROR - hard_shape: no memory\n");
         exit( 1 );
        }
        
 shape->nsteps = (int) (pw/SHAPE_RESOLUTION + 0.5);
 shape->pw = SHAPE_RESOLUTION*shape->nsteps;
 
 if ((shape->amp = (double *) malloc( sizeof(double) * shape->nsteps)) == NULL)
        {fprintf(ERROR_FILE, "ERROR - hard_shape: no memory\n");
         exit( 1 );
        }
 if ((shape->phase = (double *) malloc( sizeof(double) * shape->nsteps)) == NULL)
        {fprintf(ERROR_FILE, "ERROR - hard_shape: no memory\n");
         exit( 1 );
        }
        
 for (i=0; i!=shape->nsteps; i++)
        {shape->amp[i] = 1.0 * SHAPE_MAXAMPLITUDE;
         shape->phase[i] = 0.0;
        }
 shape->ampfactor = 1.0;
 calcrf_shape( shape, pwref, comp );
 shape->waveType = strdup("hard");
 
 if (offset != 0.0)
        { convolute_shape( shape, offset, mod );
        }
 
 write_shape( shape, name );

 if (verbose)
        list_shape( shape, stdout );

 return shape;
}       


SHAPE *sinc_shape(pw, offset, mod, pwref, comp, name, verbose)
double pw, offset, mod, pwref, comp; char *name; int verbose;
{
/*      pw:      pulse width
        offset:  offset (in Hz) of the pulse
        mod:     zero point of phase modulation pattern; begin (mod=0.0) to end (mod=1.0)
        pwref:   reference pw for calculation of the fineattenuator value; has to have the same flipangle as pw of the shape
        comp:    amplifier compression factor
        name:    VNMR shapelib name
        verbose: verbose flag
*/
 SHAPE *shape; int i;
 double stepsize, x, amp;

 if ((shape = new_shape() ) == NULL)
        {fprintf(ERROR_FILE, "ERROR - sinc_shape (new): no memory\n");
         exit( 1 );
        }
        
 shape->nsteps = (int) (pw/SHAPE_RESOLUTION + 0.5);
 shape->pw = SHAPE_RESOLUTION*shape->nsteps;
 
 if ((shape->amp = (double *) malloc( sizeof(double) * shape->nsteps)) == NULL)
        {fprintf(ERROR_FILE, "ERROR - sinc_shape (amp): no memory\n");
         exit( 1 );
        }
 if ((shape->phase = (double *) malloc( sizeof(double) * shape->nsteps)) == NULL)
        {fprintf(ERROR_FILE, "ERROR - sinc_shape (phase): no memory\n");
         exit( 1 );
        }

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
 shape->waveType = strdup("sinc");

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
{
/*      pw:      pulse width
        offset:  offset (in Hz) of the pulse
        mod:     zero point of phase modulation pattern; begin (mod=0.0) to end (mod=1.0)
        pwref:   reference pw for calculation of the fineattenuator value; has to have the same flipangle as pw of the shape
        comp:    amplifier compression factor
        name:    VNMR shapelib name
        verbose: verbose flag
*/
 SHAPE *shape; int i;
 double stepsize, x, amp;

 if ((shape = new_shape() ) == NULL)
        {fprintf(ERROR_FILE, "ERROR - sinc_shape (new): no memory\n");
         exit( 1 );
        }
        
 shape->nsteps = (int) (pw/SHAPE_RESOLUTION + 0.5);
 shape->pw = SHAPE_RESOLUTION*shape->nsteps;
 
 if ((shape->amp = (double *) malloc( sizeof(double) * shape->nsteps)) == NULL)
        {fprintf(ERROR_FILE, "ERROR - sinc_shape (amp): no memory\n");
         exit( 1 );
        }
 if ((shape->phase = (double *) malloc( sizeof(double) * shape->nsteps)) == NULL)
        {fprintf(ERROR_FILE, "ERROR - sinc_shape (phase): no memory\n");
         exit( 1 );
        }

/* half sin(x)/x from -pi to 0 (centerlobe) */
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
 shape->waveType = strdup("hsinc");

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
{
/*      pw:      pulse width
        offset:  offset (in Hz) of the pulse
        mod:     zero point of phase modulation pattern; begin (mod=0.0) to end (mod=1.0)
        pwref:   reference pw for calculation of the fineattenuator value; has to have the same flipangle as pw of the shape
        comp:    amplifier compression factor
        name:    VNMR shapelib name
        verbose: verbose flag
*/
 SHAPE *shape; int i;
 double stepsize, x, amp;

 if ((shape = new_shape() ) == NULL)
        {fprintf(ERROR_FILE, "ERROR - sinc_shape (new): no memory\n");
         exit( 1 );
        }
        
 shape->nsteps = (int) (pw/SHAPE_RESOLUTION + 0.5);
 shape->pw = SHAPE_RESOLUTION*shape->nsteps;
 
 if ((shape->amp = (double *) malloc( sizeof(double) * shape->nsteps)) == NULL)
        {fprintf(ERROR_FILE, "ERROR - sinc_shape (amp): no memory\n");
         exit( 1 );
        }
 if ((shape->phase = (double *) malloc( sizeof(double) * shape->nsteps)) == NULL)
        {fprintf(ERROR_FILE, "ERROR - sinc_shape (phase): no memory\n");
         exit( 1 );
        }

/* half sin(x)/x from -pi to 0 (centerlobe) */
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
 shape->waveType = strdup("hsinctr");

 if (offset != 0.0)
        { convolute_shape( shape, offset, mod );
        }
 
 write_shape( shape, name );

 if (verbose)
        list_shape( shape, stdout );

 return shape;
}       

SHAPE *readPbox_shape(shapename, verbose)          
char  *shapename; int verbose;
{
/* retrieve Pbox parameters from .RF file header 
   adapted from EK's getRsh routine in the Pbox_psh.h file
*/
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

SHAPE *Pbox_shape(waveType, bandwidth, offset, pwref, refpwr, myresol, options, name, verbose)          
char  *waveType, *name, *options; double bandwidth, offset, pwref, refpwr,
myresol; int verbose;
{ 
/* use Pbox to make a .RF shape.
   read name.RF file header into the structure.
   return structure pointer upon successfull completion.
   
        wavetype:       Pbox wave type; eg. g3, sinc, ...etc
        bandwidth:      Bandwidth of the pulse          
        offset:         offset (in Hz) of the pulse
        pwref:          reference pw90 for calculation of the power settings
        refpwr:         reference power 
        myresol:        time resolution (in microseconds) 
        options:        optional Pbox flags/commands/parameters
        name:           VNMR shapelib name (no .RF extension)
        verbose:        verbose flag
      
*/

  char    string[MAXLEN];
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
                   myresol,
                   options
         );
  if (verbose)
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


SHAPE *myPbox_shape(waveType, bandwidth, offset, mod, pwref, refpwr, myresol, options, name, verbose)          
char  *waveType, *name, *options; double bandwidth, offset,mod, pwref, refpwr,
myresol; int verbose;
{ 
/* use Pbox to make a .RF shape.
   read name.RF file header into the structure.
   return structure pointer upon successfull completion.
   
        wavetype:       Pbox wave type; eg. g3, sinc, ...etc
        bandwidth:      Bandwidth of the pulse          
        offset:         offset (in Hz) of the pulse
	mod		0.0 for exc, 0.5 for refoc. 1.0 for de-exc
        pwref:          reference pw90 for calculation of the power settings
        refpwr:         reference power 
        options:        optional Pbox flags/commands/parameters
        name:           VNMR shapelib name (no .RF extension)
        verbose:        verbose flag
      
*/

  char    string[MAXLEN];
  char    *getenv();
  SHAPE   *shape;
 
/* generate the Pbox command */
  sprintf(string, "Pbox %s -w \"%s %.3f %.3f %.f \" -l %.2f -p %.2f -stepsize %0.2f %s",
                   name,
                   waveType,
                   bandwidth,
                   offset,
		   mod,
                   pwref*1.0e6,
                   refpwr,
                   myresol,
                   options
         );
  if (verbose)
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


SHAPE *myPbox_shape2(waveType1, bandwidth1, offset1, mod1,waveType2, bandwidth2,
offset2, mod2, pwref, refpwr, myresol, options, name, verbose)          
char  *waveType1, *waveType2, *name, *options; double bandwidth1, bandwidth2,
offset1, offset2, mod1, mod2, pwref, refpwr, myresol; int verbose;
{ 
/* use Pbox to make a .RF shape.
   read name.RF file header into the structure.
   return structure pointer upon successfull completion.
   
        wavetype1:       Pbox wave type 1; eg. g3, sinc, ...etc
        bandwidth1:      Bandwidth 1 of the pulse          
        offset1:         offset 1 (in Hz) of the pulse
	mod1		0.0 for exc, 0.5 for refoc. 1.0 for de-exc
        wavetype2:       Pbox 2 wave type; eg. g3, sinc, ...etc
        bandwidth2:      Bandwidth 2 of the pulse          
        offset2:         offset 2(in Hz) of the pulse
	mod2		0.0 for exc, 0.5 for refoc. 1.0 for de-exc
       pwref:          reference pw90 for calculation of the power settings
        refpwr:         reference power 
        options:        optional Pbox flags/commands/parameters
        name:           VNMR shapelib name (no .RF extension)
        verbose:        verbose flag
      
*/

  char    string[MAXLEN];
  char    *getenv();
  SHAPE   *shape;
 
/* generate the Pbox command */
  sprintf(string, "Pbox %s -w \"%s %.3f %.3f %.f \" \"%s %.3f %.3f %.f \" -l %.2f -p %.2f -stepsize %0.2f %s",
                   name,
                   waveType1,
                   bandwidth1,
                   offset1,
		   mod1,
                   waveType2,
                   bandwidth2,
                   offset2,
		   mod2,
                   pwref*1.0e6,
                   refpwr,
                   myresol,
                   options
         );
  if (verbose)
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
  shape->offset = offset1;
  shape->waveType = strdup(waveType1);

  if (verbose)
        list_shape( shape, stdout );
  
  return shape;
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
{
/* set current probe */

 double Id;
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
{
/* routine can also set non-integer dB values using fine attenuator */
 
 int    ival;
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
{
/* delay-pulse-delay on two channels simultaneously */

#define N_EVENT 8

 EVENT channel1[N_EVENT], channel2[N_EVENT];
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


