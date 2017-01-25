/*********** V500 ***********/
#define V500

#define MAX_PROBES 4
static ProbeDef probes[MAX_PROBES] = {
	{"Generic", 			/* name of probe */
	  0,				/* probe id */
	 {60.0, 63.0, 63.0, 63.0},	/* max powers of Transmitter, Dec, Dec2, Dec3 */
	 {50.0, 48.0, 48.0, 48.0},	/* max decouple powers of Transmitter, Dec, Dec2, Dec3 */
	 {10.0, 10.0, 10.0, 10.0},	/* max presat powers of Transmitter, Dec, Dec2, Dec3 */
         {327.67, 327.67, 327.67}       /* 1% DAC-level of x, y, z gradient */
	},
	{"HCN-PFG", 
	  1,
	 {60.0, 63.0, 63.0, 63.0},
	 {50.0, 48.0, 48.0, 48.0},
	 {10.0, 10.0, 10.0, 10.0},	/* max presat powers of Transmitter, Dec, Dec2, Dec3 */
         {0.0, 0.0, 327.67}       	/* 1% DAC-level of x, y, z gradient */
	},
	{"HCP-PFG", 
	  2,
	 {60.0, 63.0, 63.0, 63.0},
	 {50.0, 48.0, 48.0, 48.0},
	 {10.0, 10.0, 10.0, 10.0},	/* max presat powers of Transmitter, Dec, Dec2, Dec3 */
         {0.0, 0.0, 327.67}             /* 1% DAC-level of x, y, z gradient */
	},
	{"HCN-Cold", 
	  3,
	 {54.0, 59.0, 59.0, 59.0},
	 {44.0, 45.0, 42.0, 42.0},
	 {5.0, 10.0, 10.0, 10.0},	/* max presat powers of Transmitter, Dec, Dec2, Dec3 */
         {0.0, 0.0, 327.67}             /* 1% DAC-level of x, y, z gradient */
	}
};

/* 26-03-2003: remove shaped gradient; still not working
#define GTIME( time ) 			(time + WFG_START_DELAY + WFG_STOP_DELAY)
#define GRAD(shape, gt, glx, gly, glz) 	\
        {setProbe(); shapedgradient( shape, gt, glz*probe->grad_factor[PROBE_GRADZ], 'z', 1, 1 );}
*/
#define GTIME( time ) 			(time + 2.0*GRADIENT_DELAY)
#define GRAD(shape, gt, glx, gly, glz)  \
        {setProbe(); zgradpulse( (glz)*probe->grad_factor[PROBE_GRADZ], gt);}

#define C13FREQ 125.71

/*********** end V500 ***********/
