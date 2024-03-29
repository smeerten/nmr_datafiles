/* Precompiled from vmb_2D_hydr2ar_dev_mt on Wed Apr 22 15:47:32 2015
 */
/******************************************************************************* 

SEQUENCE:	vmb_2D_hydr2ar_dev_mt
AIM:		Acquisition of 2D CT-COSY spectrum between hydrides (from pH2, 
                 evolution period) and aromatic/aliphatic protons of the 
		substrates in asymmetric Iridium-H2-(sub)-(cosub) complexes.
		
		The experiment consists of the following steps:
		1) Bubbling pH2 in the NMR tube: when p-H2 binds to the metal,
		   two non-chemically equivalent hydrides (typically resonating 
		   between -15 ppm and -30 ppm) are obtained
		2) Due to the chemical shift difference of the hydrides and the
		   uncoherent binding of p-H2 to the metal  center, the singlet
		   state of the hydrides is converted to longitudinal spin order
		    in a few ms.
		3) The initial block consists of selective 90-180-90 pulses.
		   The selectivity of the pulses is crucial; particularly, the
		   first 90 deg pulse should excite only one hydride or, at most
		   high field hydrides only. The 180 must refocus/invert ALL
		   hydrides and invert the aromatics or aliphatics of the 
		   substrate. The selectivity of the last 90 deg pulse is less
		   important; in the current implementation, for simmetry
		   reasons the pulse is the same as the starting 90 deg pulse.
		   The duration of this block is tuned to (n+1)/2JHH. If several
		   complexes are present in the sample and more than one hydride 
		   have been excited, an average value of the JHH coupling constants
		   should be used. At the end of this period hydrides magnetization 
		   is converted to longitudinal spin order between hydride and 
		   a proton of the substrate. Alternatively, a square low
		   power 45 deg pulse covering the whole hydride region can be 
		   used for excitation. In this case the starting block consists
		   of (square)45 - (shaped)180 - (square)90 pulses. At the end
		   of such block enhanced longitudinal spin order for ALL 
		   hydrides and relative substrates is obtained.
		   In addition longitudinal spin order  DQ + ZQ of hydrides is 
		   present, but it should not contribute to the measured signal.
		   During the refocusing/dephasing period dtauhh chemical shift
		   evolution of the hydrides takes place in Constant Time
		   fashion. The experiment can be acquired with gradient coherence
		   selection (echo-anti-echo)(RECOMMENDED).  
		 4)Purging gradient
		 5)Selective excitation of the aromatic/aliphatic protons.
		   In this block (selective 90-180-90) of duration
		   dtauhh2 antiphase between aromatics and hydrides is
		   refocused. The  pulses are selective on the whole
		   aromatics region. At the end of the period dtauhh2 antiphase
		   coherence is partially converted to in-phase transverse 
		   magnetization of aromatics, that can be detected. To destroy 
		   possible contribution to the acquired signal originating 
		   from residual antiphase, an optional purging block can be 
		   applied, which consists of an off-resonance 90 deg pulse on 
		   hydride protons, flanked by a pair of gradients applied in 
		   opposite direction (RECOMMENDED).

		The experiment can be run in conventional Varian fashion (set 
		by the flag fpseudo3Dmode='n' i.e. setting ni = desired number 
		of complex points, phase=1,2 for quadrature detection, and nt 
		defining the  number of scans per increment).
		When this option is selected, bubbling of p-H2 occurs at the
		beginning of each individual scan. Fluctuations in the field
		homogeneity typically generates undesired t1-noise ridge in the
		2D spectrum, only partially alleviated by gradient coherence 
		selection. 
		In alternative, the experiment can be acquired in pseuso 3D
		fashion (flag fpseudo3Dmode='y'). In this case, bubbling occurs
		only at the beginning of a 2D experiment.When bubbling stops,
		the signals of all increments is rapidly acquired, only one scan
		per increment. When the last increment has been acquired, pH2 is
		again bubbled in the NMR tube, and the quadrature components for
		each increment are acquired. This is repeated a number of times
		(corresponding to the number of transient in a regular 2D) to 
		increase the signal-to-noise ratio and for phase cycling.In order
		to acquire a 2D experiment in this fashion these parameters should
		be set as folloing:
		ni = 0
		phase = 1
		ni2 = desired number of scans
		phase2 = 1,2
		d2 = 0,1/sw1,2/sw1,3/sw1,4/sw1,5/sw1 ...
		nbubble = 2 (or 4),0,0,0,0,0,0,0, ...
		array = 'phase2,(d2,nbubble)'
		The 2D datasets can not be processed and displayed in vnmrj:
		first it must be rearranged to recover the standard Varian
		format. A macro in matlab or Octave is available.
		Note that even if gradient coherence selection is set, some
		undesired coherence pathways are possible, due to exchange
		of pH2 and substrate. Such pathways are largely reduced 
		(practically cancelled) by phase cycling, even if the scans
		that are to be combined are acquired several minutes from each
		other, under different bubbling. Longest phase cycling consists
		of 16 elements.
		
		Important: when fpseudo3Dmode='y', the exchange of  pH2
		and substrate should be optimized by changing the sample 
		temperature: if exchange is too fast, pH2 will be consumed in a
		few seconds and only a few increments in the 2D will provide
		some signals. On the other hand, when exchange is too slow,
		binding of new pH2 to the complex might take too long, 
		determining a strong reduction of the signal enhancement.

		Finally: these settings should be used for the determination of
		J couplings:
		
		FOR ARRAY OF 1D SPECTRA 

		fpseudo3Dmode = 'y' 		(IMPORTANT)
		f2D = 'n'
		fPN = 'y'			(IMPORTANT)
		fpurge_evol = 'y'		(IMPORATNT)
		nt = 1
		ni2 = desired number of scans (multiple of 2)
		ni = 1
		phase = phase2 = 1
		d2 = 0	
		dtauhh2 = fixed (typically between 120 and 300 ms)
		dtauhh arrayed, typically in steps of 50 ms, starting from 
			the longest duration down to the shortst (typ. 50 ms)
		nbubble arrayed (same number of elements as dtauhh, only the
			furst element is non-zero)
		array = '(dtauhh,nbubble)'	


		FOR ARRAY OF 2D SPECTRA 

		fpseudo3Dmode = 'y' 		(IMPORTANT)
		f2D = 'Y'
		fPN = 'y'			(IMPORTANT)
		fpurge_evol = 'y'		(IMPORATNT)
		nt = 1
		ni2 = desired number of scans (multiple of 2)
		ni = 1
		phase = 1
		phase2 = 1,2
		dtauhh2 = fixed (typically between 120 and 300 ms)
		dtauhh arrayed, typically in steps of 50 ms, starting from 
			the longest duration down to the shortst (typ. 50 ms)
		d2 = 0,1/sw1,2/sw1,3/sw1,4/sw1,...	
		nbubble arrayed (same number of elements as d2, only the
			first elent is non-zero)
		array = 'dtauhh,phase2,(d2,nbubble)'	
		

REFERENCE:    	

		
IMPLEMENTED:  	Marco Tessari
HISTORY:        12-12-2014

ADJUSTABLE PARAMETERS:
		tpwr		power level for optional square 45 deg pulse
		pw		optional low power square 45 deg pulse
		bw_hyr		bandwidth entire hydride region.
		bw_hyr_sel	bandwidth hydride region of interest
		tof_hydr	center of acquisition (hydride) region
		bw_aro		bandwidth aromatic region
		tof_aro		center of evolution region (i.e. of spectral 
				region ovserved in the indirect dimension).
		dtauhh		delay to generate antiphase hydride-substrate
				protons
		dtauhh2		delay to refocus antiphase hydride-substrate
				protons
		fshaped		flag for shaped excitation of hydrides in the first block
		fshaped2	flag for shaped excitation of hydrides @ end of the first block
				
		fpseudo3Dmode	flag for selection of bubbling scheme	
		fPN		flag for selection of gradient coherence
				selection
		fpurge_evol	flag for removal of residual antiphase to
				substrate protons before acquisition

PROCESSING:	

*******************************************************************************/







#include <standard.h>
#include "Ntools.h"

static SHAPE *eburp_hydr_zx,
	     *eburp_hydr_xz,
	     *eburp_hydr_zx_sel,
	     *eburp_hydr_xz_sel,
	     *eburp_hydr_off,
	     *reburp_hydr,
	     *eburp_aro_zx,
	     *eburp_aro_xz,
	     *eburp_aro_off,
	     *reburp_aro,
	     *reburp_INEPT_hydr,
	     *reburp_INEPT_aro,
	     *iburp_hydr_off,
	     *iburp_hydr_off_dum,
	     *eburp_purge;

pulsesequence()


{
 double back_pres_delay;	 /* system flush  */
 double bubble_open_delay;	 /* system flush  */
 double bw_aro;	 /* bandwidht for entire arom/aliph region */
 double bw_aro_sel;	 /* bandwidht for arom/aliph of interest */
 double bw_hydr;	 /* bandwidht entire hydride region */
 double bw_hydr_sel;	 /* bandwidht for hydrides of interest */
 double bw_purge;	 /* bandwidht for Boltzmann purge */
 double dcheck;	 /* preparation period  */
 double dtau;	 /* dephase delay */
 double dtauhh;	 /* delay to create antiphase (INEPT) */
 double dtauhh2;	 /* delay to refocus antiphase (revINEPT) */
 char   fCOSY[MAXSTR];	 /*flag for purge Boltzmann aromatics*/
 char   fPN[MAXSTR];	 /*flag for P/N acquisition  */
 char   fpseudo3Dmode[MAXSTR];	 /*flag for 2D spectrum  in 3D mode*/
 char   fpurge_evol[MAXSTR];	 /*flag for purge pulse on evol. region before acq  */
 char   fshaped[MAXSTR];	 /*flag for shaped 1st 90 deg on hydride */
 char   fshaped2[MAXSTR];	 /*flag for shaped 2nd 90 deg on hydrides  */
 char   fstartpurge[MAXSTR];	 /*flag for purge Boltzmann aromatics*/
 double grdelay;	 /* grad delay */
 char   gshape[MAXSTR];	 /* general gradient shape */
 double gt0;	 /* z-gradient pulse 0 */
 double gt1;	 /* z-gradient pulse 1 */
 double gt10;	 /* z-gradient pulse 9 */
 double gt2;	 /* z-gradient pulse 2 */
 double gtH;	 /* H-gradient pulse duration */
 double gx0;	 /* x-gradient level 0 */
 double gx1;	 /* x-gradient level 1 */
 double gx10;	 /* x-gradient level 10 */
 double gx2;	 /* x-gradient level 2 */
 double gxH;	 /* H x-gradient level 11 */
 double gy0;	 /* y-gradient level 0 */
 double gy1;	 /* y-gradient level 1 */
 double gy10;	 /* y-gradient level 10 */
 double gy2;	 /* y-gradient level 2 */
 double gyH;	 /* H y-gradient level 11 */
 double gz0;	 /* z-gradient level 0 */
 double gz1;	 /* z-gradient level 1 */
 double gz10;	 /* z-gradient level 10 */
 double gz2;	 /* z-gradient level 2 */
 double gzH;	 /* H z-gradient level 11 */
 double linesh_rec_delay;	 /* system flush  */
 int    nbubble;	 /* number of p-H2 bubbling loops */
 int    nloops;	 /* number of DIPSI cycles */
 double pw2;	 /* 90 deg T pulse */
 double tof_aro;	 /* offset evolution region */
 double tof_hydr;	 /* offset evolution region */
 double tof_purge;	 /* offset purge */
 double vent_open_delay;	 /* preparation period  */
 int    verbose;	 /* verbose flag (0/1) */

 back_pres_delay	= getval("back_pres_delay"); 
 bubble_open_delay	= getval("bubble_open_delay"); 
 bw_aro	= getval("bw_aro"); 
 bw_aro_sel	= getval("bw_aro_sel"); 
 bw_hydr	= getval("bw_hydr"); 
 bw_hydr_sel	= getval("bw_hydr_sel"); 
 bw_purge	= getval("bw_purge"); 
 d1	= getval("d1"); 
 d2	= getval("d2"); 
 d3	= getval("d3"); 
 dcheck	= getval("dcheck"); 
 dpwr	= getval("dpwr"); 
 dpwr2	= getval("dpwr2"); 
 dtau	= getval("dtau"); 
 dtauhh	= getval("dtauhh"); 
 dtauhh2	= getval("dtauhh2"); 
 (void) getstr("fCOSY",fCOSY); 
 (void) getstr("fPN",fPN); 
 (void) getstr("fpseudo3Dmode",fpseudo3Dmode); 
 (void) getstr("fpurge_evol",fpurge_evol); 
 (void) getstr("fshaped",fshaped); 
 (void) getstr("fshaped2",fshaped2); 
 (void) getstr("fstartpurge",fstartpurge); 
 grdelay	= getval("grdelay"); 
 (void) getstr("gshape",gshape); 
 gt0	= getval("gt0"); 
 gt1	= getval("gt1"); 
 gt10	= getval("gt10"); 
 gt2	= getval("gt2"); 
 gtH	= getval("gtH"); 
 gx0	= getval("gx0"); 
 gx1	= getval("gx1"); 
 gx10	= getval("gx10"); 
 gx2	= getval("gx2"); 
 gxH	= getval("gxH"); 
 gy0	= getval("gy0"); 
 gy1	= getval("gy1"); 
 gy10	= getval("gy10"); 
 gy2	= getval("gy2"); 
 gyH	= getval("gyH"); 
 gz0	= getval("gz0"); 
 gz1	= getval("gz1"); 
 gz10	= getval("gz10"); 
 gz2	= getval("gz2"); 
 gzH	= getval("gzH"); 
 linesh_rec_delay	= getval("linesh_rec_delay"); 
 nbubble	= (int) (getval("nbubble") + 0.5); 
 nloops	= (int) (getval("nloops") + 0.5); 
 phase2	= (int) (getval("phase2") + 0.5); 
 pw	= getval("pw"); 
 pw2	= getval("pw2"); 
 tof_aro	= getval("tof_aro"); 
 tof_hydr	= getval("tof_hydr"); 
 tof_purge	= getval("tof_purge"); 
 tpwr	= getval("tpwr"); 
 tpwr2	= getval("tpwr2"); 
 vent_open_delay	= getval("vent_open_delay"); 
 verbose	= (int) (getval("verbose") + 0.5); 

 (void) loadtable("vmb_2D_hydr2ar_dev_mt.tab"); 

{int icosel;
double myfact = 1.0;
double ref_pw90=104.25e-6;
double ref_pwr=36;

int mytest = (d3*getval("sw2")) ;

int mytest2 = (mytest / 2) % 2;

int mytest3 = (mytest / 4) % 2;

int mytest4 = (mytest / 8) % 2;


/* States or PN -TPPI indirect-1 */
 
if ( *fpseudo3Dmode == 'y' )
   {

    if ( *fPN == 'n')
       {
        if (phase2 == 2)
           {
	    tsadd(t1, 1, 4);
	    myfact = 0.0;
           }
       }
    else
       {
        myfact = 1.0;
        if (phase2 == 2)
            {icosel = 1.0;
            }
        else
            {icosel = -1.0;
	    }
       }

     if (( (int)(d2*getval("sw1")+1e-6) - (int) ((int)(d2*getval("sw1")+1e-6) / 2) * 2))
          {
          tsadd(t1, 2, 4);
          tsadd(t60, 2, 4); 
         }
      
     if ( mytest % 2 )
          {
           tsadd(t1, 2, 4); 
           tsadd(t60, 2, 4); 
          }
 
     if ( mytest2 )
          {
           tsadd(t2, 2, 4); 
           tsadd(t60, 2, 4); 
          }
 
     if ( mytest3 )
           {
           tsadd(t3, 2, 4); 
           tsadd(t60, 2, 4); 
          }

     if ( mytest4 )
          {
           tsadd(t5, 1, 4); 
          }
    }
else 
    {
     if ( *fPN == 'n')
        {
	  myfact = 0.0;
/* States-TPPI indirect-1 */
         if (phase1 == 2)
            {
	      tsadd(t1, 1, 4);
            }
        }
    else
       {
        myfact = 1.0;
/* PN-TPPI indirect-1 */
        if (phase1 == 1)
            {icosel = -1.0;
            }
        else
            {icosel = 1.0;
            }
       }
        if (NI1 % 2)
           { 
   	    tsadd(t60, 2, 4);
            tsadd(t1, 2, 4);
           }
     } 
 
setreceiver( t60 );


if ( ix == 1) 
    {	
	 eburp_hydr_zx_sel = myPbox_shape("eburp1",
	  			   bw_hydr_sel,
				   0.0,
				   0.0,
				   ref_pw90,
				   ref_pwr,
				   10.0,
				   "",
                       	   	   "eburp_hydr_zx_sel", 
                                   verbose 
                           	  );
				  
	 eburp_hydr_xz_sel = myPbox_shape("eburp1",
	  			   bw_hydr_sel,
				   0.0,
				   1.0,
				   ref_pw90,
				   ref_pwr,
				   10.0,
				   "",
                       	   	   "eburp_hydr_xz_sel", 
                                   verbose 
                           	  );

	 eburp_hydr_zx = myPbox_shape("eburp1",
	  			   bw_hydr,
				   0.0,
				   0.0,
				   ref_pw90,
				   ref_pwr,
				   10.0,
				   "",
                       	   	   "eburp_hydr_zx", 
                                   verbose 
                           	  );
				  
	 eburp_hydr_xz = myPbox_shape("eburp1",
	  			   bw_hydr,
				   0.0,
				   1.0,
				   ref_pw90,
				   ref_pwr,
				   10.0,
				   "",
                       	   	   "eburp_hydr_xz", 
                                   verbose 
                           	  );
	 eburp_hydr_off = myPbox_shape("eburp1",
	  			   bw_hydr,
				   (tof_hydr-tof_aro),
				   0.0,
				   ref_pw90,
				   ref_pwr,
				   10.0,
				   "",
                       	   	   "eburp_hydr_off", 
                                   verbose 
                           	  );
	 reburp_hydr = myPbox_shape("reburp",
	  			   bw_hydr,
				   0.0,
				   0.5,
				   ref_pw90,
				   ref_pwr,
				   10.0,
				   "",
                       	   	   "reburp_hydr", 
                                   verbose 
                           	  );

	 eburp_aro_zx = myPbox_shape("eburp1",
	  			   bw_aro,
				   0.0,
				   0.0,
				   ref_pw90,
				   ref_pwr,
				   10.0,
				   "",
                       	   	   "eburp_aro_zx", 
                                   verbose 
                           	  );

	 eburp_aro_xz = myPbox_shape("eburp1",
	  			   bw_aro,
				   0.0,
				   1.0,
				   ref_pw90,
				   ref_pwr,
				   10.0,
				   "",
                       	   	   "eburp_aro_xz", 
                                   verbose 
                           	  );

	 eburp_aro_off = myPbox_shape("eburp1",
	  			   bw_aro,
				   (-tof_hydr+tof_aro),
				   0.0,
				   ref_pw90,
				   ref_pwr,
				   10.0,
				   "",
                       	   	   "eburp_aro_off", 
                                   verbose 
                           	  );

	 reburp_INEPT_hydr = myPbox_shape2("reburp",
	  			   bw_hydr,
				   0.0,
				   0.5,
				   "reburp",
	  			   bw_aro,
				   (-tof_hydr+tof_aro),
				   0.0,
				   ref_pw90,
				   ref_pwr,
				   10.0,
				   "",
                       	   	   "reburp_INEPT_hydr", 
                                   verbose 
                           	  );
				  
	 reburp_INEPT_aro = myPbox_shape2("reburp",
	  			   bw_aro,
				   0.0,
				   0.5,
				   "reburp",
	  			   bw_hydr,
				   (-tof_aro+tof_hydr),
				   0.0,
				   ref_pw90,
				   ref_pwr,
				   10.0,
				   "",
                       	   	   "reburp_INEPT_aro", 
                                   verbose 
                           	  );

	 reburp_aro = myPbox_shape("reburp",
	  			   bw_aro,
				   0.0,
				   0.5,
				   ref_pw90,
				   ref_pwr,
				   10.0,
				   "",
                       	   	   "reburp_aro", 
                                   verbose 
                           	  );

	 iburp_hydr_off = myPbox_shape("iburp2",
	  			   bw_hydr,
				   (-tof_aro+tof_hydr),
				   0.0,
				   ref_pw90,
				   ref_pwr,
				   10.0,
				   "",
                       	   	   "iburp_hydr_off", 
                                   verbose 
                           	  );
	  
	 iburp_hydr_off_dum = myPbox_shape("iburp2",
	  			   bw_hydr,
				   (-tof_aro-tof_hydr),
				   0.0,
				   ref_pw90,
				   ref_pwr,
				   10.0,
				   "",
                       	   	   "iburp_hydr_off_dum", 
                                   verbose 
                           	  );
				  
	 eburp_purge = myPbox_shape("eburp1",
	  			   bw_purge,
				   tof_purge,
				   0.0,
				   ref_pw90,
				   ref_pwr,
				   10.0,
				   "",
                       	   	   "eburp_purge", 
                                   verbose 
                           	  );
       



  }
/* initialise hardloop counter */
   initval((double) nbubble, v10);

/* initialise hardloop counter */
 initval((double) nloops, v1);
if (ix == 1)
{
printf("> Mixing time=%f\n",(115.104 *pw2)*nloops ); 
}


status( A );

LOCK_ON; 
 DECOUPLE_POWER( XDEV, dpwr, us(10) );
 DECOUPLE_POWER( YDEV, dpwr2, us(10) );
 
 recoff();
 delay(d1);
 POWER( TDEV,tpwr, us(10) );
 
       if ( nbubble > 0 )
         { 
          starthardloop( v10 );
             sp1on();
	       delay( vent_open_delay ); 
             sp1off();
LOCK_OFF;

             sp2on();
	       delay( bubble_open_delay);
             sp2off();

             sp3on();
	       delay( back_pres_delay);
             sp3off();

               delay(linesh_rec_delay);
          endhardloop( );
	  delay(dcheck);
         }
LOCK_OFF;
	 
    if ( *fstartpurge == 'y')
        {
          OFFSET( TDEV, tof_purge, us(5) );
           STPULSE(eburp_purge, zero, rof1, 0.0);
           GRAD( gshape, gt10, gx10, gy10, gz10);
           delay(grdelay);
        }

 status( B );
 
    OFFSET( TDEV, tof_hydr, us(5) );
    POWER( TDEV, tpwr, us(10.0));	 

/* Hydrides excitation, refocusing of antiphase between hydrides, CT t1
evolution e COSY transfer to aromatics */

	 

if ( *fshaped == 'y' )
    {
	STPULSE( eburp_hydr_zx_sel, t1, rof1, 0.0 );
   }
else
    {
         rgpulse(pw,t1,rof1,0.0);
    }	  

           delay(dtauhh*0.5 + d2*0.5 - (gt0 + grdelay) );
           GRAD( gshape, gt0, gx0, gy0, gz0);
	   delay(grdelay - rof1);
	   
        STPULSE(reburp_INEPT_hydr, one, rof1, 0.0);

           GRAD( gshape, gt0, gx0, gy0, gz0 );
           delay(dtauhh*0.5 - d2*0.5 - gt0 - us(10) - (gtH+grdelay) );
           GRAD( gshape, gtH, gxH*icosel*myfact , gyH*icosel*myfact, gzH*icosel*myfact );
           delay(grdelay - rof1);
if ( *fshaped2 == 'y' )
    {
           delay(us(10.0));
	STPULSE( eburp_hydr_xz_sel, t2, rof1, 0.0 );
    }
else
    {
            POWER( TDEV, tpwr,us(10) );  
         rgpulse(pw*2.0,t2,rof1,0.0);
    }

/* PURGING GRADIENT*/
 
            OFFSET( TDEV, tof_aro, us(5) );
            GRAD( gshape, gt1, gx1, gy1, gz1 );
            delay( grdelay);

/* ON AROMATICS:  */

       STPULSE(eburp_aro_zx, t3, rof1, 0.0);
       
 if ( *fCOSY == 'y' ) 
     {
	  delay(dtauhh2*0.5 - (gt2 + grdelay));
          GRAD( gshape, gt2, gx2, gy2, gz2 );
          delay(grdelay - rof1 );
	 
       STPULSE(reburp_aro, t4, rof1, 0.0);      
       
         GRAD( gshape, gt2, gx2, gy2, gz2 );
         delay( grdelay  );
	 delay(dtauhh2*0.5  - (gt2 + grdelay)  );

        STPULSE(eburp_aro_xz, t5, rof1, 0.0);
 
      if ( *fPN == 'y' )
         {
          GRAD( gshape, gtH, gxH*myfact, gyH*myfact, gzH*myfact );
	  delay(grdelay);
         }
       }

 if ( *fCOSY == 'n' ) 
     {
 POWER( TDEV, tpwr2, us(5) );

      if ( *fPN == 'y' )
         {
          GRAD( gshape, gtH, gxH*myfact, gyH*myfact, gzH*myfact );
	  delay(grdelay);
         }

       if (nloops >= 1.0)
        {
         starthardloop(v1);

            rgpulse(pw2*2.000, zero, rof1, 0.0);
              
            rgpulse(pw2*1.556, zero, rof1, 0.0);
            rgpulse(pw2*3.556, two, rof1, 0.0);
              

            rgpulse(pw2*1.000, two, rof1, 0.0);
            rgpulse(pw2*3.000, zero, rof1, 0.0);
              
              
            rgpulse(pw2*0.222, zero, rof1, 0.0);
            rgpulse(pw2*2.222, two, rof1, 0.0);
              
              
            rgpulse(pw2*0.944, two, rof1, 0.0);
            rgpulse(pw2*0.333, zero, rof1, 0.0);
            rgpulse(pw2*1.389, two, rof1, 0.0);
              

            rgpulse(pw2*1.333, two, rof1, 0.0);
            rgpulse(pw2*3.333, zero, rof1, 0.0);
              
              
            rgpulse(pw2*0.833, zero, rof1, 0.0);
            rgpulse(pw2*2.833, two, rof1, 0.0);
              
              
            rgpulse(pw2*0.111, two, rof1, 0.0);
            rgpulse(pw2*2.111, zero, rof1, 0.0);
              
              
            rgpulse(pw2*2.000, zero, rof1, 0.0);
              
            
            rgpulse(pw2*2.000, two, rof1, 0.0);
              

            rgpulse(pw2*1.556, two, rof1, 0.0);
            rgpulse(pw2*3.556, zero, rof1, 0.0);
              

            rgpulse(pw2*1.000, zero, rof1, 0.0);
            rgpulse(pw2*3.000, two, rof1, 0.0);
              
              
            rgpulse(pw2*0.222, two, rof1, 0.0);
            rgpulse(pw2*2.222, zero, rof1, 0.0);
              
              
            rgpulse(pw2*0.944, zero, rof1, 0.0);
            rgpulse(pw2*0.333, two, rof1, 0.0);
            rgpulse(pw2*1.389, zero, rof1, 0.0);
              

            rgpulse(pw2*1.333, zero, rof1, 0.0);
            rgpulse(pw2*3.333, two, rof1, 0.0);
              
              
            rgpulse(pw2*0.833, two, rof1, 0.0);
            rgpulse(pw2*2.833, zero, rof1, 0.0);
              
              
            rgpulse(pw2*0.111, zero, rof1, 0.0);
            rgpulse(pw2*2.111, two, rof1, 0.0);
              
              
            rgpulse(pw2*2.000, two, rof1, 0.0);
              
            
            rgpulse(pw2*2.000, two, rof1, 0.0);
              

            rgpulse(pw2*1.556, two, rof1, 0.0);
            rgpulse(pw2*3.556, zero, rof1, 0.0);
              

            rgpulse(pw2*1.000, zero, rof1, 0.0);
            rgpulse(pw2*3.000, two, rof1, 0.0);
              
              
            rgpulse(pw2*0.222, two, rof1, 0.0);
            rgpulse(pw2*2.222, zero, rof1, 0.0);
              
              
            rgpulse(pw2*0.944, zero, rof1, 0.0);
            rgpulse(pw2*0.333, two, rof1, 0.0);
            rgpulse(pw2*1.389, zero, rof1, 0.0);
              

            rgpulse(pw2*1.333, zero, rof1, 0.0);
            rgpulse(pw2*3.333, two, rof1, 0.0);
              
              
            rgpulse(pw2*0.833, two, rof1, 0.0);
            rgpulse(pw2*2.833, zero, rof1, 0.0);
              
              
            rgpulse(pw2*0.111, zero, rof1, 0.0);
            rgpulse(pw2*2.111, two, rof1, 0.0);
              
              
            rgpulse(pw2*2.000, two, rof1, 0.0);
              
            
            rgpulse(pw2*2.000, zero, rof1, 0.0);
              

            rgpulse(pw2*1.556, zero, rof1, 0.0);
            rgpulse(pw2*3.556, two, rof1, 0.0);
              

            rgpulse(pw2*1.000, two, rof1, 0.0);
            rgpulse(pw2*3.000, zero, rof1, 0.0);
              
              
            rgpulse(pw2*0.222, zero, rof1, 0.0);
            rgpulse(pw2*2.222, two, rof1, 0.0);
              
              
            rgpulse(pw2*0.944, two, rof1, 0.0);
            rgpulse(pw2*0.333, zero, rof1, 0.0);
            rgpulse(pw2*1.389, two, rof1, 0.0);
              

            rgpulse(pw2*1.333, two, rof1, 0.0);
            rgpulse(pw2*3.333, zero, rof1, 0.0);
              
              
            rgpulse(pw2*0.833, zero, rof1, 0.0);
            rgpulse(pw2*2.833, two, rof1, 0.0);
              
              
            rgpulse(pw2*0.111, two, rof1, 0.0);
            rgpulse(pw2*2.111, zero, rof1, 0.0);
              
              
            rgpulse(pw2*2.000, zero, rof1, 0.0);
              
         endhardloop();
        }

	  delay(dtau*0.5 - (gt2 + grdelay));
          GRAD( gshape, gt2, gx2, gy2, gz2 );
          delay(grdelay - rof1 );
	 
       STPULSE(reburp_aro, t4, rof1, 0.0);      
       
          GRAD( gshape, gt2, gx2, gy2, gz2 );
          delay( grdelay  );
  	  delay(dtau*0.5  - (gt2 + grdelay)  );

       }
LOCK_ON; 

       
       recon();
 status(C);

}
}
