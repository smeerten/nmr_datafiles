/* Precompiled from zgtest_mt8 on Fri Aug 24 13:23:19 2012
 */
/******************************************************************************* 

SEQUENCE:	JRselinv_mt
AIM:		1D exchange experiment.
		Scheme: selective inversion by jump-return - exchange delay -
		90 degree pulse - acquisition 
REFERENCE:    	M. Piotto, V. Saudek & V. Sklenar, J. Biomol. NMR 2, 661 - 666 (1992)

		V. Sklenar, M. Piotto, R. Leppik $ V. Saudek, J. Magn. Reson. Series A 102, 241 -245 (1993)

IMPLEMENTED:  	Marco Tessari
HISTORY:      	16-03-2012

ADJUSTABLE PARAMETERS:
		satflg		presat on ('ynn') or off ('nnn')
		satfrq		presat frequency
		satpwr		presat power
		pw		1H 90 at tpwr
		pw2		1H 90 at tpwr , typically pw-0.1
		tpwr
		fshape		shaped water-flipback pulse
		pw3		water flip-back pulse at tpwr3
		wgdelay		watergate delay
		dm,dmm		decoupling on X-channel
		dm2,dmm2 	decoupling on Y-channel
		dm2,dmm2 	decoupling on Z-channel

PROCESSING:	wft

*******************************************************************************/






#include <standard.h>
#include "Ntools.h"

pulsesequence()


{
 double dtau;	 /* dephase delay */
 char   fgrad[MAXSTR];	 /* grad flag*/
 char   fnloopsat[MAXSTR];	 /* flag for presat with nloopsat*/
 double grdelay;	 /* grad delay */
 char   gshape[MAXSTR];	 /* general gradient shape */
 double gt1;	 /* z-gradient pulse 1 */
 double gx1;	 /* x-gradient level 1 */
 double gy1;	 /* y-gradient level 1 */
 double gz1;	 /* z-gradient level 1 */
 int    nloopsat;	 /* number of presat loops  */
 double phcorr;	 /* phase corr. on water flip-back pulse*/
 double pw5;	 /* 1H pw*/
 char   satflg[MAXSTR];	 /* sat flag */
 double satfrq1;	 /* 1st presat freq witn nloopsat  */
 double satfrq2;	 /* 2nd presat freq witn nloopsat  */
 double wgdelay;	 /* WATERGATE  delay */

 d1	= getval("d1"); 
 dpwr	= getval("dpwr"); 
 dpwr2	= getval("dpwr2"); 
 dtau	= getval("dtau"); 
 (void) getstr("fgrad",fgrad); 
 (void) getstr("fnloopsat",fnloopsat); 
 grdelay	= getval("grdelay"); 
 (void) getstr("gshape",gshape); 
 gt1	= getval("gt1"); 
 gx1	= getval("gx1"); 
 gy1	= getval("gy1"); 
 gz1	= getval("gz1"); 
 nloopsat	= (int) (getval("nloopsat") + 0.5); 
 phcorr	= getval("phcorr"); 
 pw	= getval("pw"); 
 pw5	= getval("pw5"); 
 satdly	= getval("satdly"); 
 (void) getstr("satflg",satflg); 
 satfrq	= getval("satfrq"); 
 satfrq1	= getval("satfrq1"); 
 satfrq2	= getval("satfrq2"); 
 satpwr	= getval("satpwr"); 
 tof	= getval("tof"); 
 tpwr	= getval("tpwr"); 
 wgdelay	= getval("wgdelay"); 

 (void) loadtable("zgtest_mt8.tab"); 

{
 STEPSIZE( TDEV, 1.0 );
 initval((double) ((int) (phcorr+0.5)), v5);
 
  /* initialise hardloop nloops counter */
 initval((double) nloopsat, v10);

/* 
 CHECK( dm,     "nn*" ); 
 CHECK( dm2,    "nn*" ); 
 CHECK( dm3,    "nn*" ); 
 CHECK( satflg, "*nn" );
 */
 
setreceiver( t60 );


status( A );

LOCK_ON; 
 DECOUPLE_POWER( XDEV, dpwr, us(10) );
 DECOUPLE_POWER( YDEV, dpwr2, us(10) );
 
 if (satflg[A] == 'y')
	{
         POWER( TDEV,satpwr, us(10) );
	 delay( d1 - satdly );
	 
          if ( *fnloopsat == 'y' )
	     {
              starthardloop( v10 );

	         OFFSET( TDEV, satfrq1, us(5) );
	       rgpulse( 0.65*satdly/nloopsat, zero, rof1, 0.0);
	         OFFSET( TDEV, satfrq2, us(5) );
	       rgpulse( 0.35*satdly/nloopsat, zero, rof1, 0.0);

              endhardloop( );
	      OFFSET( TDEV, tof, us(5) );
	     }
	  else
	     {
	       OFFSET( TDEV, satfrq , us(5) );
    	      rgpulse( satdly, zero, rof1, 0.0);
	       OFFSET( TDEV, tof, us(5) );
	     }
	}
 else
	{delay( d1 );
	}


 LOCK_OFF;
  recoff();

 status( B );
 if ( *fgrad == 'y' )
   {
    GRAD( gshape, gt1, gx1, gy1, gz1 );
    delay(grdelay);  	
    }
  else
    delay(gt1 + grdelay);

    POWER( TDEV,tpwr, dtau - grdelay );
    delay(us(100));
  rgpulse(pw, t1, rof1, rof2);

      recon();
 status(C);

}
}
