/*********** Generic probe def file ***********/

#define MAX_PROBES 1
static ProbeDef probes[MAX_PROBES] = {
        {"Generic",                     /* name of probe */
          0,                            /* probe id */
         {63.0, 63.0, 63.0, 63.0},      /* max powers of Transmitter, Dec, Dec2, Dec3 */
         {50.0, 48.0, 48.0, 48.0},      /* max decouple powers of Transmitter, Dec, Dec2, Dec3 */
         {10.0, 10.0, 10.0, 10.0},      /* max presat powers of Transmitter, Dec, Dec2, Dec3 */
         {327.67, 327.67, 327.67}       /* 1% DAC-level of x, y, z gradient */
        }
};
#define DISABLE_PROBE_CHECK

/****** adjust for your spectrometer ******/

/* simple rectangular Z-gradient */
#define GTIME( time )                   (time + 2.0*GRADIENT_DELAY)
#define GRAD(shape, gt, glx, gly, glz)  {setProbe(); zgradpulse( (glz)*probe->grad_factor[PROBE_GRADZ], gt);}

/* triple axis gradient 
#define GTIME( time ) 			(time + 6.0*GRADIENT_DELAY)
#define GRAD(shape, gt, glx, gly, glz)  \
	{setProbe(); \
	 rgradient('z', (glz)*probe->grad_factor[PROBE_GRADZ]); \
	 rgradient('x', (glx)*probe->grad_factor[PROBE_GRADX]); \
	 rgradient('y', (gly)*probe->grad_factor[PROBE_GRADY]); \
	 delay( gt ); \
	 rgradient('z', 0.0); rgradient('x', 0.0); rgradient('y', 0.0); \
	}
*/

#define C13FREQ 125.71

/*********** end generic probe def file ***********/
