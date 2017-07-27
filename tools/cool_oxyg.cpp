/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*CoolOxyg evaluate total cooling due to oxygen */
#include "cddefines.h"
#include "coolheavy.h"
#include "dense.h"
#include "taulines.h"
#include "h2.h"
#include "phycon.h"
#include "embesq.h"
#include "hmi.h"
#include "oxy.h"
#include "colden.h"
#include "mole.h"
#include "ligbar.h"
#include "thermal.h"
#include "lines_service.h"
#include "atoms.h"
#include "cooling.h"

#if defined (__ICC) && defined(__ia64) && __INTEL_COMPILER < 910
#pragma optimization_level 1
#endif
void CoolOxyg(void)
{
	double a21, 
	  a31, 
	  a32, 
	  a41, 
	  a42, 
	  a43, 
	  a51, 
	  a52, 
	  a53, 
	  a54, 
	  a6300, 
	  a6363, 
	  aeff,
	  cs2s2p, 
	  cs2s3p,
	  p[5],
	  r12 , 
	  r13;

	static double go2[5]={4.,6.,4.,4.,2.};
	static double exo2[4]={26808.4,21.0,13637.5,1.5};
	/* these will be used to update change in cs wrt temperature,
	 * and its affect on cooling derivative */
	static double oi_cs_tsave=-1. , oi_te_tsave=-1. , oi_dcdt_tsave=-1.;
	long int i;
	double rate_OH_dissoc;

	DEBUG_ENTRY( "CoolOxyg()" );

	{
		enum{DEBUG_LOC=false};
		if( DEBUG_LOC )
		{
			// simple unit test
			double TeTest = phycon.TEMP_LIMIT_LOW;
			double oi_a, oi_b, oi_c, oi_d, oi_e, oi_f;
			double oii_a,
				oii_b,
				oii_c,
				oii_d,
				oii_e,
				oii_f,
				oii_g,
				oii_h,
				oii_i,
				oii_j,
				oii_k;
			double oiii_a,
				oiii_b,
				oiii_c,
				oiii_d,
				oiii_e,
				oiii_f,
				oiii_g,
				oiii_h,
				oiii_i,
				oiii_j;
			double oiv_a,oiv_b;
			double ov_a,ov_b;
			double sii_a,
				sii_b,
				sii_c,
				sii_d,
				sii_e,
				sii_f,
				sii_g,
				sii_h,
				sii_i,
				sii_j,
				sii_k;
			double siii_a,
				siii_b,
				siii_c,
				siii_d,
				siii_e,
				siii_f,
				siii_g,
				siii_h,
				siii_i,
				siii_j,
				siii_k,
				siii_l;
			double siv_a;
			double sviii_a;
			double neiii_a,neiii_b,neiii_c,neiii_d,neiii_e;
			TempChange( TeTest , true );
			oi_cs(oi_a, oi_b, oi_c, oi_d, oi_e, oi_f);
			oii_cs(oii_a, oii_b, oii_c, oii_d, oii_e, oii_f, oii_g,
				oii_h, oii_i, oii_j, oii_k);
			oiii_cs(oiii_a, oiii_b, oiii_c, oiii_d, oiii_e, oiii_f,
				oiii_g, oiii_h, oiii_i, oiii_j);
			oiv_cs(oiv_a, oiv_b);
			ov_cs(ov_a, ov_b);
			sii_cs(sii_a, sii_b, sii_c, sii_d, sii_e, sii_f, sii_g,
				sii_h, sii_i, sii_j, sii_k);
			siii_cs(siii_a, siii_b, siii_c, siii_d, siii_e, siii_f,
				siii_g, siii_h, siii_i, siii_j, siii_k,siii_l);
			siv_cs(siv_a);
			sviii_cs(sviii_a);
			neiii_cs(neiii_a,neiii_b,neiii_c,neiii_d,neiii_e);
			double tinc = 0.05;
			for( TeTest=phycon.TEMP_LIMIT_LOW; TeTest<phycon.TEMP_LIMIT_HIGH;
				  TeTest *=(1+tinc) )
			{
				double oi_aold=oi_a,
					oi_bold=oi_b,
					oi_cold=oi_c,
					oi_dold=oi_d,
					oi_eold=oi_e,
					oi_fold=oi_f;
				double oii_aold=oii_a,
					oii_bold=oii_b,
					oii_cold=oii_c,
					oii_dold=oii_d,
					oii_eold=oii_e,
					oii_fold=oii_f,
					oii_gold=oii_g,
					oii_hold=oii_h,
					oii_iold=oii_i,
					oii_jold=oii_j,
					oii_kold=oii_k;
				double oiii_aold=oiii_a,
					oiii_bold=oiii_b,
					oiii_cold=oiii_c,
					oiii_dold=oiii_d,
					oiii_eold=oiii_e,
					oiii_fold=oiii_f,
					oiii_gold=oiii_g,
					oiii_hold=oiii_h,
					oiii_iold=oiii_i,
					oiii_jold=oiii_j;
				double oiv_aold=oiv_a,oiv_bold=oiv_b;
				double ov_aold=ov_a,ov_bold=ov_b;
				double sii_aold=sii_a,
					sii_bold=sii_b,
					sii_cold=sii_c,
					sii_dold=sii_d,
					sii_eold=sii_e,
					sii_fold=sii_f,
					sii_gold=sii_g,
					sii_hold=sii_h,
					sii_iold=sii_i,
					sii_jold=sii_j,
					sii_kold=sii_k;
				double siii_aold=siii_a,
					siii_bold=siii_b,
					siii_cold=siii_c,
					siii_dold=siii_d,
					siii_eold=siii_e,
					siii_fold=siii_f,
					siii_gold=siii_g,
					siii_hold=siii_h,
					siii_iold=siii_i,
					siii_jold=siii_j,
					siii_kold=siii_k,
					siii_lold=siii_l;
				double siv_aold=siv_a;
				double sviii_aold=sviii_a;
				double neiii_aold=neiii_a,
					neiii_bold=neiii_b,
					neiii_cold=neiii_c,
					neiii_dold=neiii_d,
					neiii_eold=neiii_e;
				TempChange( TeTest , true );
				oi_cs(oi_a, oi_b, oi_c, oi_d, oi_e, oi_f);
				ASSERT( fabs(oi_a-oi_aold) <= oi_a*3.*tinc);				
				ASSERT( fabs(oi_b-oi_bold) <= oi_b*3.*tinc);
				ASSERT( fabs(oi_c-oi_cold) <= oi_c*3.*tinc);
				ASSERT( fabs(oi_d-oi_dold) <= oi_d*3.*tinc);
				ASSERT( fabs(oi_e-oi_eold) <= oi_e*3.*tinc);
				ASSERT( fabs(oi_f-oi_fold) <= oi_f*3.*tinc);
				oii_cs(oii_a, oii_b, oii_c, oii_d, oii_e, oii_f,
					oii_g, oii_h, oii_i, oii_j, oii_k);
				ASSERT( fabs(oii_a-oii_aold) <= oii_a*3.*tinc);
				ASSERT( fabs(oii_b-oii_bold) <= oii_b*3.*tinc);
				ASSERT( fabs(oii_c-oii_cold) <= oii_c*3.*tinc);
				ASSERT( fabs(oii_d-oii_dold) <= oii_d*3.*tinc);
				ASSERT( fabs(oii_e-oii_eold) <= oii_e*3.*tinc);
				ASSERT( fabs(oii_f-oii_fold) <= oii_f*3.*tinc);
				ASSERT( fabs(oii_g-oii_gold) <= oii_g*3.*tinc);
				ASSERT( fabs(oii_h-oii_hold) <= oii_h*3.*tinc);
				ASSERT( fabs(oii_i-oii_iold) <= oii_i*3.*tinc);
				ASSERT( fabs(oii_j-oii_jold) <= oii_j*3.*tinc);
				ASSERT( fabs(oii_k-oii_kold) <= oii_k*3.*tinc);
				oiii_cs(oiii_a, oiii_b, oiii_c, oiii_d, oiii_e,
					oiii_f, oiii_g, oiii_h, oiii_i, oiii_j);
				ASSERT( fabs(oiii_a-oiii_aold) <= oiii_a*3.*tinc);
				ASSERT( fabs(oiii_b-oiii_bold) <= oiii_b*3.*tinc);
				ASSERT( fabs(oiii_c-oiii_cold) <= oiii_c*3.*tinc);
				ASSERT( fabs(oiii_d-oiii_dold) <= oiii_d*3.*tinc);
				ASSERT( fabs(oiii_e-oiii_eold) <= oiii_e*3.*tinc);
				ASSERT( fabs(oiii_f-oiii_fold) <= oiii_f*3.*tinc);
				ASSERT( fabs(oiii_g-oiii_gold) <= oiii_g*3.*tinc);
				ASSERT( fabs(oiii_h-oiii_hold) <= oiii_h*3.*tinc);
				ASSERT( fabs(oiii_i-oiii_iold) <= oiii_i*3.*tinc);
				ASSERT( fabs(oiii_j-oiii_jold) <= oiii_j*3.*tinc);
				oiv_cs(oiv_a, oiv_b);
				ASSERT( fabs(oiv_a-oiv_aold) <= oiv_a*3.*tinc);
				ASSERT( fabs(oiv_b-oiv_bold) <= oiv_b*3.*tinc);
				ov_cs(ov_a, ov_b);
				ASSERT( fabs(ov_a-ov_aold) <= ov_a*3.*tinc);
				ASSERT( fabs(ov_b-ov_bold) <= ov_b*3.*tinc);
				sii_cs(sii_a, sii_b, sii_c, sii_d, sii_e, sii_f,
					sii_g, sii_h, sii_i, sii_j, sii_k);
				ASSERT( fabs(sii_a-sii_aold) <= sii_a*3.*tinc);
				ASSERT( fabs(sii_b-sii_bold) <= sii_b*3.*tinc);
				ASSERT( fabs(sii_c-sii_cold) <= sii_c*3.*tinc);
				ASSERT( fabs(sii_d-sii_dold) <= sii_d*3.*tinc);
				ASSERT( fabs(sii_e-sii_eold) <= sii_e*3.*tinc);
				ASSERT( fabs(sii_f-sii_fold) <= sii_f*3.*tinc);
				ASSERT( fabs(sii_g-sii_gold) <= sii_g*3.*tinc);
				ASSERT( fabs(sii_h-sii_hold) <= sii_h*3.*tinc);
				ASSERT( fabs(sii_i-sii_iold) <= sii_i*3.*tinc);
				ASSERT( fabs(sii_j-sii_jold) <= sii_j*3.*tinc);
				ASSERT( fabs(sii_k-sii_kold) <= sii_k*3.*tinc);
				siii_cs(siii_a,
					siii_b,
					siii_c,
					siii_d,
					siii_e,
					siii_f,
					siii_g,
					siii_h,
					siii_i,
					siii_j,
					siii_k,
					siii_l);
				ASSERT( fabs(siii_a-siii_aold) <= siii_a*3.*tinc);
				ASSERT( fabs(siii_b-siii_bold) <= siii_b*3.*tinc);
				ASSERT( fabs(siii_c-siii_cold) <= siii_c*3.*tinc);
				ASSERT( fabs(siii_d-siii_dold) <= siii_d*3.*tinc);
				ASSERT( fabs(siii_e-siii_eold) <= siii_e*3.*tinc);
				ASSERT( fabs(siii_f-siii_fold) <= siii_f*3.*tinc);
				ASSERT( fabs(siii_g-siii_gold) <= siii_g*3.*tinc);
				ASSERT( fabs(siii_h-siii_hold) <= siii_h*3.*tinc);
				ASSERT( fabs(siii_i-siii_iold) <= siii_i*3.*tinc);
				ASSERT( fabs(siii_j-siii_jold) <= siii_j*3.*tinc);
				ASSERT( fabs(siii_k-siii_kold) <= siii_k*3.*tinc);
				ASSERT( fabs(siii_l-siii_lold) <= siii_l*3.*tinc);
				siv_cs(siv_a);
				ASSERT( fabs(siv_a-siv_aold) <= siv_a*3.*tinc);
				sviii_cs(sviii_a);
				ASSERT( fabs(sviii_a-sviii_aold) <= sviii_a*3.*tinc);
				neiii_cs(neiii_a, neiii_b, neiii_c, neiii_d, neiii_e);
				ASSERT( fabs(neiii_a-neiii_aold) <= neiii_a*3.*tinc);
				ASSERT( fabs(neiii_b-neiii_bold) <= neiii_b*3.*tinc);
				ASSERT( fabs(neiii_c-neiii_cold) <= neiii_c*3.*tinc);
				ASSERT( fabs(neiii_d-neiii_dold) <= neiii_d*3.*tinc);
				ASSERT( fabs(neiii_e-neiii_eold) <= neiii_e*3.*tinc);
			}
			cdEXIT(EXIT_SUCCESS);
		}
	}

	/***********************************************************************
	**************************************O I*******************************
	***********************************************************************/
	/* following does the OI Bowen Ly-bet pumping problem */
	atom_oi_calc(&CoolHeavy.coolOi);
	CoolAdd("o  1",8446,CoolHeavy.coolOi);
	double oi_cs3P23P1,
		oi_cs3P23P0,
		oi_cs3P13P0,
		oi_cs3P1D2,
		oi_cs1D21S0,
		oi_cs3P1S0;
	/*"oi_cs" calculates electron collision strengths for O I.
	  *>>refer	o1	cs	Bell, Berrington & Thomas 1998, MNRAS 293, L83
	 * cs variables are named based of the transition they represent (lower level first).
	 *  Transitions where triplet states are considered to be a single state
	 *  are denoted by omitting the J value from the triplet level
	 *  ex. (3P(J=2,1,0) -> 1D2 = oi_cs3P1D2
	 * written by Kirk Korista, 29 may 96
	 * adapted by Peter van Hoof, 30 march 99 (to include Bell et al. data)
	 * all data within the ground state 3P triplet, above 3000K, have been 
	 * adjusted down by a constant factor to make them line up with Bell et al. data.
	 *>>chng 10 mar 6 ML: combined oi3Pcs with the other cs calculations for O I.*/
	oi_cs(oi_cs3P23P1,
		oi_cs3P23P0,
		oi_cs3P13P0,
		oi_cs3P1D2,
		oi_cs1D21S0,
		oi_cs3P1S0);
	double csh01=-1.,
		cshe01=-1.,
		csh201=-1.,
		csh12=-1.,
		cshe12=-1.,
		csh212=-1.,
		csh02=-1.,
		cshe02=-1.,
		csh202 =-1.,
		csh2o01=-1.,
		csh2o02=-1.,
		csh2o12=-1.,
		csh2p01=-1.,
		csh2p02=-1.,
		csh2p12=-1.,
		csp01=-1.,
		csp02=-1.,
		csp12=-1.;
	/*"oi_othercs" calculates non-electron collision strengths for O I.
	 *The variables declared here have 2 digit indices at the end of the name.
	 *  The lower level is specified first, then the upper level.
	 * The levels are numbered in order of increasing energy. 
	 * Below is a small table to convert from the index used here (i) to J
	 * for the O I 3P triplet.
	 * i|J
	 * 0|2
	 * 1|1
	 * 2|0     */
	oi_othercs(csh01,cshe01,csh201,csh12,cshe12,csh212,csh02,cshe02,csh202,
	    csh2o01,csh2o02,csh2o12,csh2p01,csh2p02,csh2p12,csp01,csp02,csp12);

	oi_cs3P23P1 = oi_cs3P23P1+csp01+3.*(csh01*dense.xIonDense[ipHYDROGEN][0]
		+ cshe01*dense.xIonDense[ipHELIUM][0] + csh201*hmi.H2_total)/
		dense.cdsqte;
	oi_cs3P13P0 = oi_cs3P13P0+csp12+(csh12*dense.xIonDense[ipHYDROGEN][0] +
		cshe12*dense.xIonDense[ipHELIUM][0] + csh212*hmi.H2_total)/
		dense.cdsqte;
	oi_cs3P23P0 = oi_cs3P23P0+csp02+(csh02*dense.xIonDense[ipHYDROGEN][0] +
		cshe02*dense.xIonDense[ipHELIUM][0] + csh202*hmi.H2_total)/
		dense.cdsqte;

	/* O I 6300, 6363, A from
	* >>refer	all	all	Mendoza, C. 1982, in Planetary Nebulae, IAU Symp No. 103,
	* >>refercon ed by D.R. Flower, (D. Reidel: Holland), 143 */
	a6300 = TauLines[ipT6300].Emis().Aul()*TauLines[ipT6300].Emis().Pesc();
	TauLines[ipT6300].Emis().PopOpc() = (dense.xIonDense[ipOXYGEN][0]*5./5.);
	(*TauLines[ipT6300].Lo()).Pop() = (dense.xIonDense[ipOXYGEN][0]*5./5.);
	(*TauLines[ipT6300].Hi()).Pop() = 0.;
	TauLines[ipT6300].Coll().col_str() = (realnum)(oi_cs3P1D2*5./9.);
	TauLines[ipT6363].Emis().PopOpc() = (dense.xIonDense[ipOXYGEN][0]*5./3.);
	(*TauLines[ipT6363].Lo()).Pop() = (dense.xIonDense[ipOXYGEN][0]*5./3.);
	(*TauLines[ipT6363].Hi()).Pop() = 0.;
	TauLines[ipT6363].Coll().col_str() = (realnum)(oi_cs3P1D2*3./9.);
	a6363 = TauLines[ipT6363].Emis().Aul()*TauLines[ipT6363].Emis().Pesc();
	a21 = a6300 + a6363 + oxy.d6300;
	a32 = TauLines[ipT5577].Emis().Aul()*TauLines[ipT5577].Emis().Pesc();
	PutCS(oi_cs1D21S0,TauLines[ipT5577]);
	/* rate of new populations of O^0 due to dissociation of OH,
	 * co.rate_OH_dissoc is rate OH -> O + H [cm-3 s-1],
	 * must make it per unit O atom, so this rate is s-1 excitations per O atom */
	rate_OH_dissoc = mole.findrate("PHOTON,OH=>O,H");
	r12 = rate_OH_dissoc*0.55/SDIV( dense.xIonDense[ipOXYGEN][0] );
	r13 = rate_OH_dissoc*0.05/SDIV( dense.xIonDense[ipOXYGEN][0] );
	/* below is correction for fraction of excitations the produce emission */
	CoolHeavy.c6300_frac_emit = (a6300+a6363)/(a6300+a6363+oi_cs3P1D2*dense.cdsqte/5.);
	CoolHeavy.c5577_frac_emit = (a32)/(a32+oi_cs1D21S0*dense.cdsqte/3.);
	/* d6300 is the photoionization rate from the excited level
	 * was computed when ionization balance done */
	CoolHeavy.c5577 = atom_pop3(9.,5.,1.,oi_cs3P1D2,oi_cs3P1S0,oi_cs1D21S0,a21,7.56e-2,a32,
	  22590.,25775.7,&oxy.poiexc,dense.xIonDense[ipOXYGEN][0],0.,r12,r13)*a32*
	  3.57e-12;
	TauLines[ipT5577].Emis().PopOpc() = oxy.poiexc;
	(*TauLines[ipT5577].Lo()).Pop() = oxy.poiexc;
	(*TauLines[ipT5577].Hi()).Pop() = 0.;
	/* convert poiexc to relative abundances */
	/* >>chng 04 apr 20, include correction for fraction of 6300 due to OH pump */
	CoolHeavy.c5577 *= (1.-r13/(SDIV(atoms.c13)));
	CoolHeavy.c6300 = oxy.poiexc*a6300*TauLines[ipT6300].EnergyErg() *
		(1.-r12/(SDIV(atoms.c12)));
	CoolHeavy.c6363 = oxy.poiexc*a6363*TauLines[ipT6363].EnergyErg() *
		(1.-r12/(SDIV(atoms.c12)));
	/* must introduce correction for fraction of 6300 that is photo produced */
	thermal.dCooldT += CoolHeavy.c6300*(2.28e4*thermal.tsq1 + thermal.halfte) *
		/* note that atoms.c12 has all ra tes 1->2 */
		(1.-r12/(SDIV(atoms.c12)));
	oxy.poiexc /= (realnum)MAX2(1e-20,dense.xIonDense[ipOXYGEN][0]);
	CoolAdd("o  1",5577,CoolHeavy.c5577);
	CoolAdd("o  1",6300,CoolHeavy.c6300);
	CoolAdd("o  1",6363,CoolHeavy.c6363);	
	PutCS(oi_cs3P23P1,TauLines[ipT63]);
	PutCS(oi_cs3P13P0,TauLines[ipT146]);
	PutCS(oi_cs3P23P0,*TauDummy);
	atom_level3(TauLines[ipT63],TauLines[ipT146],*TauDummy);
	/* now save pops to add col den in radinc */
	for( i=0; i<3; ++i)
		colden.O1Pops[i] = (realnum)atoms.PopLevels[i];
	/* >>chng 02 jul 25, keep track of change in cs for 63 micron line */
	if( !fp_equal(phycon.te,oi_te_tsave) )
	{
		/* very first time we come through, previous values are -1 */
		if(oi_te_tsave>0. )
			oi_dcdt_tsave = (oi_cs3P23P1-oi_cs_tsave) /
				(phycon.te-oi_te_tsave);
		else
			oi_dcdt_tsave = 0.;
		oi_cs_tsave = oi_cs3P23P1;
		oi_te_tsave	= phycon.te;
		/* >>chng 03 jan 21, this factor could become very large - it should
		 * always be positive since neutral cs's are, and not much bigger than
		 * the usual derivative */
		/* can't be negative */
		oi_dcdt_tsave = MAX2( 0. , oi_dcdt_tsave);
		/* can't be bigger than several times normal dC/dT */
		oi_dcdt_tsave = MIN2( TauLines[ipT63].EnergyK()*thermal.tsq1*4.,oi_dcdt_tsave);
	}
	/* this is only derivative due to change in collision strength, which is capped to be
	 * less than 4x the thermal derivative just above */
	thermal.dCooldT += TauLines[ipT63].Coll().cool()*oi_dcdt_tsave;
	/* this is a bebug print statement for the numerical cs derivative */
	{
		enum{DEBUG_LOC=false};
		if( DEBUG_LOC )
		{
			fprintf(ioQQQ,"DEBUG OI\t%.2f\tte\t%.5e\tcool\t%.5e\tcs\t%.4e\told\t%.4e\tnew\t%.4e\n",
				fnzone,
				phycon.te, 
				TauLines[ipT63].Coll().cool() , 
				TauLines[ipT63].Coll().col_str() , 
				TauLines[ipT63].Coll().cool()*TauLines[ipT63].EnergyK()*thermal.tsq1, 
				TauLines[ipT63].Coll().cool()*oi_dcdt_tsave );
		}
	}

	/***********************************************************************
	**************************************O II******************************
	***********************************************************************/
	double oii_cs4S32D5,
		oii_cs4S32D3,
		oii_cs2D52D3,
		oii_cs4S32P3,
		oii_cs2D52P3,
		oii_cs2D32P3,
		oii_cs4S32P1,
		oii_cs2D52P1,
		oii_cs2D32P1,
		oii_cs2P32P1,
		oii_cs4S34P;
	/* >>chng 10 feb 14, update cs to
	 * >>referold	o2	cs	McLaughlin, B.M., & Bell, K.L. 1998, J Phys B 31, 4317 
	 * >>chng 02 mar 13, go back to older values as per Seaton/Osterbrock correspondence
	 * >>chng 04 nov 01, statistical weights were reversed, caught by Kevin Blagrave 
	 * >>refer	o2	as	Froese Fischer, C., & Tachiev, G. 2004, At. Data Nucl. Data Tables, 87, 1
	 */
	a21 = 4.12e-5;
	a31 = 1.63e-4;
	a32 = 1.24e-7;
	a41 = 5.65e-2;
	a42 = 1.11e-1;
	a43 = 5.87e-2;
	a51 = 2.27e-2;
	a52 = 5.82e-2;
	a53 = 9.67e-2;
	a54 = 3.15e-10;
	/*"oii_cs" calculates collision strengths for O II.
	 * cs variables are named based of the transition they represent (lower level first).
	 * Transitions where triplet states are considered to be a single state
	 * are denoted by omitting the J value from the triplet level
	 * ex. (3P(J=2,1,0) -> 1D2 = oi_cs3P1D2
	 * When the J is a half integer, only the numerator is given.
	 */
	oii_cs(oii_cs4S32D5,
		oii_cs4S32D3,
		oii_cs2D52D3,
		oii_cs4S32P3,
		oii_cs2D52P3,
		oii_cs2D32P3,
		oii_cs4S32P1,
		oii_cs2D52P1,
		oii_cs2D32P1,
		oii_cs2P32P1,
		oii_cs4S34P);

	double Cooling , CoolingDeriv;
	atom_pop5(go2,exo2,oii_cs4S32D5,oii_cs4S32D3,oii_cs4S32P3,oii_cs4S32P1,
		oii_cs2D52D3,oii_cs2D52P3,oii_cs2D52P1,oii_cs2D32P3,oii_cs2D32P1,
		oii_cs2P32P1,a21,a31,a41,a51,a32,a42,a52,a43,a53,a54,p,
		dense.xIonDense[ipOXYGEN][1], &Cooling , &CoolingDeriv, 0.,0.,0.,0.);

	CoolHeavy.O3730 = (realnum)(p[1]*a21*5.34e-12);
	CoolHeavy.O3726 = (realnum)(p[2]*a31*5.34e-12);
	CoolHeavy.O2471 = (realnum)((p[3]*a41 + p[4]*a51)*8.05e-12);
	CoolHeavy.O7323 = (realnum)((p[3]*a42 + p[4]*a52)*2.72e-12);
	CoolHeavy.O7332 = (realnum)((p[3]*a43 + p[4]*a53)*2.71e-12);
	CoolHeavy.c3727 = CoolHeavy.O3730 + CoolHeavy.O3726;
	CoolHeavy.c7325 = CoolHeavy.O7323 + CoolHeavy.O7332;

	// total cooling from lowest five levels of O II
	CoolAdd("O  2",3727,Cooling );
	thermal.dCooldT += CoolingDeriv;

	/* remember ratio of radiative to total decays, to use for estimating
	 * recombination contribution in lines_lv1_li_ne */
	if( (p[3] + p[4]) > SMALLFLOAT )
		CoolHeavy.O2_A3_tot = (p[3]*(a41+a42+a43) + p[4]*(a51+a52+a53) ) /
			( (p[3]*(a41+a42+a43) + p[4]*(a51+a52+a53) ) +
			( p[3]*(oii_cs4S32P3+oii_cs2D52P3+oii_cs2D32P3)/go2[3] +
			p[4]*(oii_cs4S32P1+oii_cs2D52P1+oii_cs2D32P1)/go2[4]) *
			dense.cdsqte );
	else
		CoolHeavy.O2_A3_tot = 0.;

	if( (p[1] + p[2]) > SMALLFLOAT )
		CoolHeavy.O2_A2_tot = (p[1]*a21 + p[2]*a31 ) /
			( (p[1]*a21 + p[2]*a31 ) +
			( p[1]*oii_cs4S32D5/go2[1] + p[2]*oii_cs4S32D3/go2[2]) *
			dense.cdsqte );
	else
		CoolHeavy.O2_A2_tot = 0.;

	/* O II 4S34P 833.9A, CS
	 * >>refer	o2	cs	McLaughlin, B.M., & Bell, K.L. 1993, ApJ, 408, 753 */
	/* >>chng 01 aug 10, turn this line back on - no reason given for turning it off. */
	PutCS(oii_cs4S34P,TauLines[ipT834]);
	atom_level2(TauLines[ipT834]);

	/***********************************************************************
	**************************************O III*****************************
	***********************************************************************/	
	double oiii_cs3P25S2,
		oiii_cs3P15S2,
		oiii_cs3P05S2,
		oiii_cs3P1D2,
		oiii_cs1D21S0,
		oiii_cs3P1S0,
		oiii_cs3P03P1,
		oiii_cs3P13P2,
		oiii_cs3P03P2,
		oiii_cs3P3D;
	/*"oiii_cs" calculates collision strengths for O III.
	* cs variables are named based of the transition they represent (lower level first).
	*Transitions where triplet states are considered to be a single state
	* are denoted by omitting the J value from the triplet level
	* ex. (3P(J=2,1,0) -> 1D2 = oi_cs3P1D2 */
	oiii_cs(oiii_cs3P25S2,
		oiii_cs3P15S2,
		oiii_cs3P05S2,
		oiii_cs3P1D2,
		oiii_cs1D21S0,
		oiii_cs3P1S0,
		oiii_cs3P03P1,
		oiii_cs3P13P2,
		oiii_cs3P03P2,
		oiii_cs3P3D);
	
	PutCS(oiii_cs3P25S2,TauLines[ipT1666]);
	PutCS(oiii_cs3P15S2,TauLines[ipT1661]);
	PutCS(oiii_cs3P05S2,*TauDummy);
	atom_level3(*TauDummy,TauLines[ipT1666],TauLines[ipT1661]);
	TauLines[ipT304].Emis().PopOpc() = dense.xIonDense[ipOXYGEN][2];
	(*TauLines[ipT304].Lo()).Pop() = dense.xIonDense[ipOXYGEN][2];
	(*TauLines[ipT304].Hi()).Pop() = 0.;

	/* o iii 5007+4959, As 96 NIST */
	/*The cs of the transitions 3P0,1,2 to 1D2 are added together to give oiii_cs3P1D2 */
	/*the cs of the transition 1D2-1S0 is mentioned as oiii_cs1D21S0*/
	/*The cs of the transitions 3P0,1,2 to 1S0 are added together to give oiii_cs3P1S0*/
	aeff = 0.027242 + oxy.d5007r;
	a21 = aeff - oxy.d5007r;   
	a31 = 0.2262;           
	a32 = 1.685;            
	oxy.o3ex23 = 32947.;
	oxy.o3br32 = (realnum)(a32/(a31 + a32));
	oxy.o3enro = (realnum)(4.56e-12/3.98e-12);
	// solve a 3 level system, collapsing ^3P ground term into single level
	/* POP3(G1,G2,G3,O12,O13,O23,A21,A31,A32,E12,E23,P2,ABUND,GAM2) */
	oxy.poiii3 = (realnum)(atom_pop3(9.,5.,1.,oiii_cs3P1D2,oiii_cs3P1S0,
		oiii_cs1D21S0,a21,a31,a32,28990.,oxy.o3ex23,&oxy.poiii2,
		dense.xIonDense[ipOXYGEN][2],oxy.d5007r,0.,0.));
	CoolHeavy.c4363 = oxy.poiii3*a32*4.56e-12;
	CoolHeavy.c5007 = oxy.poiii2*a21*3.98e-12;
	oxy.d5007t = (realnum)(CoolHeavy.c5007*oxy.d5007r/aeff);
	thermal.dCooldT += CoolHeavy.c5007*(2.88e4*thermal.tsq1 - thermal.halfte);
	thermal.dCooldT += CoolHeavy.c4363*6.20e4*thermal.tsq1;
	CoolAdd("O  3",5007,CoolHeavy.c5007);
	CoolAdd("O  3",4363,CoolHeavy.c4363);
	CoolAdd("O  3",2331,CoolHeavy.c4363*0.236);

	if( MAX2(oxy.poiii2,oxy.poiii3) > 0.f && dense.xIonDense[ipOXYGEN][2]>SMALLFLOAT)
	{
		oxy.poiii3 /= dense.xIonDense[ipOXYGEN][2];
		oxy.poiii2 /= dense.xIonDense[ipOXYGEN][2];
	}

	/* O III IR lines */
	PutCS(oiii_cs3P03P1,TauLines[ipTO88]);
	PutCS(oiii_cs3P13P2,TauLines[ipT52]);
	PutCS(oiii_cs3P03P2,*TauDummy);
	atom_level3(TauLines[ipTO88],TauLines[ipT52],*TauDummy);

	/* O III 3P to 3D triplets 835 angstroms */
	PutCS(oiii_cs3P3D,TauLines[ipT835]);
	atom_level2(TauLines[ipT835]);

	/***********************************************************************
	**************************************O IV******************************
	***********************************************************************/
	double oiv_cs2P2D,oiv_cs2P12P3;
	/*"oiv_cs" calculates collision strengths for O IV.
	 * cs variables are named based of the transition they represent (lower level first).
	 * Transitions where triplet states are considered to be a single state
	 * are denoted by omitting the J value from the triplet level
	 * ex. (3P(J=2,1,0) -> 1D2 = oi_cs3P1D2
	 * For half integer Js, only the numerator is given.  */
	oiv_cs(oiv_cs2P2D,oiv_cs2P12P3);
	/* O IV 789A,  2P2D CS from
	 * >>refer	o4	cs	Zhang, H.L., Graziani, M., Pradhan, A.K. 1994, A&A, 283, 319 */	
	PutCS(oiv_cs2P2D,TauLines[ipT789]);
	atom_level2(TauLines[ipT789]);
	/* O IV 26 micron, CS 
	 * >>referold	o4	cs	Blum, R.D., & Pradhan, A.K. 1992, ApJS 80, 425
	 * A=
	 * >>refer	o4	as	Brage, T., Judge, P.G., & Brekke, P. 1996, ApJ. 464, 1030 */
	/* >>chng 06 nov 08 - NPA.  Update collision strength to new data from:
	 * >>refer	o4	cs	Tayal, S. 2006, ApJS 166, 634 
	 * Equation derived by using TableCurve, and goes to zero as 
	 * T => 0 and T => infinity */	
	PutCS(oiv_cs2P12P3,TauLines[ipT26]);
	static vector< pair<TransitionList::iterator,double> > O4Pump;
	O4Pump.reserve(48);
	/* one time initialization if first call */
	if( O4Pump.empty() )
	{
		// set up level 1 pumping lines
		pair<TransitionList::iterator,double> pp( TauLines.begin()+ipT789, 1./6. ); 
		O4Pump.push_back( pp );
		// set up level 2 pumping lines
		for( i=0; i < nWindLine; ++i )
		{
			/* don't test on nelem==ipIRON since lines on physics, not C, scale */
			if( (*TauLine2[i].Hi()).nelem() == 8 && (*TauLine2[i].Hi()).IonStg() == 4 )
			{
#				if	0
				DumpLine( TauLine2.begin()+i );
#				endif
				double branch_ratio;
				// the branching ratios used here ignore cascades via intermediate levels
				// usually the latter are much slower, so this should be reasonable
				if( fp_equal( (*TauLine2[i].Hi()).g(), realnum(2.) ) )
					branch_ratio = 2./3.; // 2S upper level
				else if( fp_equal( (*TauLine2[i].Hi()).g(), realnum(6.) ) )
					branch_ratio = 1./2.; // 2P upper level
				else if( fp_equal( (*TauLine2[i].Hi()).g(), realnum(10.) ) )
					branch_ratio = 1./6.; // 2D upper level
				else
					TotalInsanity();
				pair<TransitionList::iterator,double> pp2( TauLine2.begin()+i, branch_ratio ); 
				O4Pump.push_back( pp2 );
			}
		}
	}

	/* now sum pump rates */
	double pump_rate = 0.;
	vector< pair<TransitionList::iterator,double> >::const_iterator o4p;
	for( o4p=O4Pump.begin(); o4p != O4Pump.end(); ++o4p )
	{
		const TransitionList::iterator t = o4p->first;
		double branch_ratio = o4p->second;
		pump_rate += (*t).Emis().pump()*branch_ratio;
#		if	0
		dprintf( ioQQQ, "O IV %.3e %.3e\n",
			 (*t).WLAng , (*t).Emis().pump()*branch_ratio );
#		endif
	}

	/*atom_level2(TauLines[ipT26]);*/
	/*AtomSeqBoron compute cooling from 5-level boron sequence model atom */
	/* >>refer	o4	cs	Blum, R.D., & Pradhan, A.K., 1992, ApJS 80, 425
	 * >>refer	o4	cs	Zhang, H.L., Graziani, M., Pradhan, A.K. 1994, A&A, 283, 319 */
	AtomSeqBoron(TauLines[ipT26], 
	  TauLines[ipO4_1400], 
	  TauLines[ipO4_1397], 
	  TauLines[ipO4_1407], 
	  TauLines[ipO4_1405], 
	  TauLines[ipO4_1401], 
	  0.1367 , 0.1560 , 1.1564 , 0.0457 , pump_rate,"O  4");
	/***********************************************************************
	**************************************O V*******************************
	***********************************************************************/

	double ov_cs1S01P1,ov_cs1S03P;
	/*"ov_cs" calculates collision strengths for O V.
	* cs variables are named based of the transition they represent (lower level first).
	*Transitions where triplet states are considered to be a single state
	*are denoted by omitting the J value from the triplet level
	*ex. (3P(J=2,1,0) -> 1D2 = oi_cs3P1D2  */
	ov_cs(ov_cs1S01P1,ov_cs1S03P);	
	PutCS(ov_cs1S01P1,TauLines[ipT630]);
	atom_level2(TauLines[ipT630]);	
	/* >>chng 01 sep 09, AtomSeqBeryllium will reset this to 1/3 so critical density correct */
	PutCS(ov_cs1S03P,TauLines[ipT1214]);
	/* c1214 = AtomSeqBeryllium( .87,1.05,3.32, t1214, .0216) * 1.64E-11
	 * AtomSeqBeryllium(CS23,CS24,CS34,tarray,A41)
	 * A's 
	 * >>refer	o5	as	Fleming, J., Bell, K.L, Hibbert, A., Vaeck, N., Godefroid, M.R.
	 * >>refercon	1996, MNRAS, 279, 1289 */
	AtomSeqBeryllium(.87,0.602,2.86,TauLines[ipT1214],.02198);
	embesq.em1218 = (realnum)(atoms.PopLevels[3]*0.0216*1.64e-11);

	/* O VI 1035 li seq
	 * generate collision strengths, then stuff them in
	 * >>refer	o6	vs	Cochrane, D.M., & McWhirter, R.W.P. 1983, PhyS, 28, 25 */
	ligbar(8,TauLines[ipT1032],TauLines[ipT150],&cs2s2p,&cs2s3p);
	PutCS(cs2s2p,TauLines[ipT1032]);
	PutCS(cs2s2p*0.5,TauLines[ipT1037]);
	/* no data for the 2-3 transition */
	PutCS(1.0,*TauDummy);
	/* solve the 3 level atom */
	atom_level3(TauLines[ipT1037],*TauDummy,TauLines[ipT1032]);

	PutCS(cs2s3p,TauLines[ipT150]);
	atom_level2(TauLines[ipT150]);
	return;
}

/*"oi_cs" calculates electron collision strengths for O I.
  *>>refer	o1	cs	Bell, Berrington & Thomas 1998, MNRAS 293, L83 
 * cs variables are named based of the transition they represent (lower level first).
 * Transitions where triplet states are considered to be a single state are
 * denoted by omitting the J value from the triplet level
 * ex. (3P(J=2,1,0) -> 1D2 = oi_cs3P1D2
 * written by Kirk Korista, 29 may 96
 * adapted by Peter van Hoof, 30 march 99 (to include Bell et al. data)
 * all data within the ground state 3P triplet, above 3000K, have been adjusted
 * down by a constant factor to make
 * them line up with Bell et al. data. 
 *>>chng 10 mar 6 ML: combined oi3Pcs with the other cs calculations for O I. */

void oi_cs(double& oi_cs3P23P1,double& oi_cs3P23P0,double& oi_cs3P13P0,
		   double& oi_cs3P1D2,double& oi_cs1D21S0,double& oi_cs3P1S0)
{
	double a, 
	  b, 
	  c, 
	  d;

	DEBUG_ENTRY( "oi_cs()" );
	/* local variables */

	/* 3P2 - 3P1 */
	if( phycon.te <= 3.0e3 )
	{
		oi_cs3P23P1 = 1.49e-4*phycon.sqrte/phycon.te02/phycon.te02;
	}
	else if( phycon.te <= 1.0e4 )
	{
		a = -5.5634127e-04;
		b = 8.3458102e-08;
		c = 2.3068232e-04;
		oi_cs3P23P1 = 0.228f*(a + b*phycon.te32 + c*phycon.sqrte);
	}
	else
	{
		oi_cs3P23P1 = 0.228*MIN2(0.222,5.563e-06*phycon.te*phycon.te05*
		  phycon.te02);
	}
	
	/* 3P2 - 3P0 */
	if( phycon.te <= 3.0e3 )
	{
		oi_cs3P23P0 = 4.98e-5*phycon.sqrte;
	}
	else if( phycon.te <= 1.0e4 )
	{
		a = -3.7178028e-04;
		b = 2.0569267e-08;
		c = 1.1898539e-04;
		oi_cs3P23P0 = 0.288*(a + b*phycon.te32 + c*phycon.sqrte);
	}
	else
	{
		oi_cs3P23P0 = 0.288*MIN2(0.0589,1.015e-05*phycon.te/phycon.te10/
		  phycon.te02/phycon.te005);
	}
	
	/* 3P1 - 3P0 */
	if( phycon.te <= 3.0e3 )
	{
		oi_cs3P13P0 = 1.83e-9*phycon.te*phycon.te30*phycon.te05;
	}
	else if( phycon.te <= 1.0e4 )
	{
		a = -5.9364373e-04;
		b = 0.02946867;
		c = 10768.675;
		d = 3871.9826;
		oi_cs3P13P0 = 0.0269*(a + b*exp(-0.5*POW2((phycon.te-c)/d)));
	}
	else
	{
		oi_cs3P13P0= 0.0269*MIN2(0.074,7.794e-08*phycon.te32/phycon.te10/
		  phycon.te01);
	}
	
	/* [OI] 6300, 6363, 5575, etc
	 * >>chng 06 oct 02, Humeshkar Nemala incorporate Barklem cs data for OI
	 * largest difference is 3P - 1D (6300+6363) which is now roughly 3x smaller */
	/* This is the transition from 3P(J=2,1,0;the levels are reversed) to 1D2
	 * The rate coefficients were converted to CS
	 *>>refer	oi	cs	Barklem,P.S.,2006,A&A (astroph 0609684)
	 * Data pts are avilable over 1000,3000,5000,8000,12000,20000 and 50000
	 * Fits between 1000K and 20000K are reliable;this is not the case
	 * between 20000K & 50000K*/

	if(phycon.te < 8E3)
		oi_cs3P1D2 = (4E-08)*(phycon.te*phycon.te70*phycon.te05);
	else if(phycon.te < 2E4)
		oi_cs3P1D2 = (4.630155E-05)*((phycon.te/phycon.te04)*phycon.te005*phycon.te0001);
	else
		oi_cs3P1D2 = (1.5394E-03)*(phycon.sqrte*phycon.te10*phycon.te01*phycon.te001*phycon.te0003);

	/* this block adds on collisional excitation by H0 */
	/* >>chng 06 aug 18, add atomic hydrogen collisional processes using rates from
	 * >>refer	O1	coll	Krems, R.V., Jameson, M.J. & Dalgarno, A. 2006, ApJ, 647, 1531
	 * their equation 10 - the deecxiation rate coefficient, cm3 s-1
	 * >>referold	oi	cs	Federman, S.R., & Shipsey, E.J. 1983, ApJ, 269, 791 */
	double te_scale = phycon.te / 6000.;
	double rate_H0 = (1.74*te_scale + 0.6)*1e-12*sexp(0.47*te_scale) / sqrt(te_scale ) *
		dense.xIonDense[ipHYDROGEN][0];
	oi_cs3P1D2 += ConvRate2CS( 5.f , (realnum)rate_H0 );

	if(phycon.te < 5E3)
		oi_cs1D21S0 = (7E-08)*(phycon.te*phycon.sqrte*phycon.te10*phycon.te007*phycon.te0001);
	else if(phycon.te<2E4)
		oi_cs1D21S0 = (1.98479e-04)*((phycon.te70/phycon.te03)*phycon.te003*phycon.te0007);
	else
		oi_cs1D21S0 = (9.31E-04)*(phycon.sqrte*phycon.te01*phycon.te007*phycon.te0005*phycon.te0001);

	/* 1,2,3 -> 5 transition */
	if(phycon.te < 2E4)
		oi_cs3P1S0 = (2E-05)*(phycon.sqrte*phycon.te30*phycon.te05*phycon.te01*(phycon.te004/phycon.te0002));
	else
		oi_cs3P1S0 = (1.054E-03)*(phycon.sqrte/phycon.te04)*phycon.te003*phycon.te0005;


	return;
}
/*"oi_othercs" calculates non-electron collision strengths for O I.
 *The variables declared here have 2 digit indices at the end of the name.
 *The lower level is specified first, then the upper level. The levels are
 *numbered in order of increasing energy. Below is a small table to convert from
 *the index used here (i) to J for the O I 3P triplet.
	 * i|J
	 * 0|2
	 * 1|1
	 * 2|0       */
void oi_othercs(double& csh01,double& cshe01,double& csh201,double& csh12,double& cshe12,
	double& csh212,double& csh02,double& cshe02,double& csh202,double& csh2o01,
	double& csh2o02,double& csh2o12,double& csh2p01,double& csh2p02,double& csh2p12,
	double& csp01,double& csp02,double& csp12)
{
	DEBUG_ENTRY( "oi_othercs()" );

	/* O I fine structure lines rad data from
	 * >>refer	all	cs	Mendoza, C. 1982, in Planetary Nebulae, IAU Symp No. 103,
	 * >>refercon ed by D.R. Flower, (D. Reidel: Holland), 143
	 * >>refer	o1	cs	Berrington, K.A. 1988, J.Phys. B, 21, 1083
	 * hydrogen collisions from
	 * >>refer	oi	cs	Tielens, A.G.G., & Hollenbach, D. 1985, ApJ, 291, 722
	 * factor in () is their rate coef
	 * assume H2 and H are same
	 * CDSQTE = 8.629E-6*EDEN/SQRTE
	 * cs01 = 9.8e-6*te + (4.2e-12*te70/te03) / cdsqte * 3. * hdcor
	 * cs12 = 2.6e-6*te  + (1.5e-10*sqrte/te03/te03)/cdsqte*hdcor
	 * cs02 = 2.9e-6*te  + (1.1e-12*te70*te10)/cdsqte*hdcor
	 * evaluate fits to OI electron rates, indices on var do not agree
	 * with Kirk's in sub, but are OK */
	/*==============================================================*/
	/* >>>chng 99 jun 01,
	 * following changed to be parallel to Peter van Hoof's changes
	 * in the Fortran C90.05 version
	 * following is collisions with electrons */

	/* remember the electron part of the cs */
	/* these were added by Peter van Hoof to update the collision
	 * data within the OI ground term */


	/* rate coefficients for collisional de-excitation of O^o(3P) with H^o(2S1/2)
	 * NOTE: due to lack of data these relations are extrapolated to higher Te !
	 * >>refer	o1	cs	Launay & Roueff 1977, AA 56, 289
	 * the first fit is for Te <= 300K, the second for higher temps
	 * these data are valid for 50K <= Te <= 1000K*/
	csh12 = MAX2(2.00e-11*phycon.te30*phycon.te05*phycon.te02,
		7.67e-12*phycon.sqrte*phycon.te03);

	/* these data are valid for 100K <= Te <= 1000K */
	csh01 = MIN2(3.12e-12*phycon.te70*phycon.te02*phycon.te02,
		7.51e-12*phycon.sqrte*phycon.te05*phycon.te03);

	/* these data are valid for 150K <= Te <= 1000K*/
	csh02 = MIN2(6.96e-13*phycon.te/phycon.te10/phycon.te02,
		1.49e-12*phycon.te70*phycon.te05);

	/*rate coefficients for collisional de-excitation of O^o(3P) with H^o(2S1/2)
	 * NOTE: due to lack of data these relations are extrapolated to higher Te !
	 * >>refer	o1	cs	Abrahamsson, E., Krems, R. V., & Dalgarno, A. 2007, ApJ 654
	 * fit by Morisset 2016/03 */
	csh12 = MAX2(MAX2(3.8e-10*phycon.te10*phycon.te05,1.7e-11*phycon.te70),8.3e-11*phycon.te20*phycon.te20);
	csh01 = 5.2e-11*phycon.te20*phycon.te20*phycon.te02;
	csh02 = 5.2e-11*phycon.te20*phycon.te20/phycon.te01;

	/* rate coefficients for collisional de-excitation of O^o(3P) with He^o(1S)
	 * NOTE: due to lack of data these relations are extrapolated to higher Te !
	 * >>refer	oi	cs	Monteiro & Flower 1987, MNRAS 228, 101
	 * the first fit is for Te <= 300K, the second for higher temps
	 * these data are valid for 100K <= Te <= 1000K */
	cshe12 = MIN2(8.09e-16*phycon.te32*phycon.te10*phycon.te03,
		9.72e-15*phycon.te*phycon.te20);

	cshe01 = 1.57e-12*phycon.te70/phycon.te01;

	cshe02 = MIN2(1.80e-12*phycon.te70*phycon.te05,
		4.45e-12*phycon.te70/phycon.te10);

	if( phycon.te<=1.5e3 )
	{
		/* >>chng 04 mar 15, use explicit ortho-para densities */
		double ortho_frac = h2.ortho_density/SDIV(hmi.H2_total);
		/* rate coefficients for collisional de-excitation of O^o(3P) with H2(J=1,0)
		 * >>refer	oi	cs	Jaquet et al. 1992, J.Phys.B 25, 285
		 * these data are valid for 40K <= Te <= 1500K
		 * the first entry is contribution from ortho H2, the second para H2.*/
		csh2o12 = 2.21e-14*phycon.te*phycon.te10/phycon.te01;
		csh2p12 = 9.45e-15*phycon.te*phycon.te20/phycon.te01;
		csh212 = ortho_frac*csh2o12 + (1.-ortho_frac)*csh2p12;

		csh2o01 = 2.37e-11*phycon.te30*phycon.te10/phycon.te02;
		csh2p01 = 3.40e-11*phycon.te30*phycon.te02;
		csh201 = ortho_frac*csh2o01 + (1.-ortho_frac)*csh2p01;

		csh2o02 = 4.86e-11*phycon.te30*phycon.te02*phycon.te02;
		csh2p02 = 6.82e-11*phycon.te30/phycon.te02;
		csh202 = ortho_frac*csh2o02 + (1.-ortho_frac)*csh2p02;
	}
	else
	{
		csh212 = 0.;
		csh201 = 0.;
		csh202 = 0.;
	}

	/* effective collision strength of O^o(3P) with p
	 * >>refer	oi	cs	Pequignot, D. 1990, A&A 231, 499
	 * original data:
	 * >>refer	oi	cs	Chambaud et al., 1980, J.Phys.B, 13, 4205 (upto 5000K)
	 * >>refer	oi	cs	Roueff, private communication (10,000K and 20,000K)*/
	if( phycon.te<=1000. )
	{
		csp01 = MAX2(2.22e-5*phycon.te/phycon.te10,
			/* >>chng 05 jul 05, eden to dense.EdenHCorr */
			/*2.69e-6*phycon.te*phycon.te30)*dense.xIonDense[ipHYDROGEN][1]/dense.eden;*/
			2.69e-6*phycon.te*phycon.te30)*dense.xIonDense[ipHYDROGEN][1]/dense.EdenHCorr;

		csp02 = MIN2(7.07e-8*phycon.te32*phycon.te10,2.46e-7*
			/* >>chng 05 jul 05, eden to dense.EdenHCorr */
			/*phycon.te32/phycon.te10)*dense.xIonDense[ipHYDROGEN][1]/dense.eden;*/
			phycon.te32/phycon.te10)*dense.xIonDense[ipHYDROGEN][1]/dense.EdenHCorr;
	}
	else
	{
		csp01 = MIN2(2.69e-6*phycon.te*phycon.te30,9.21e-5*phycon.te/phycon.te10/
			/* >>chng 05 jul 05, eden to dense.EdenHCorr */
			/*phycon.te03)*dense.xIonDense[ipHYDROGEN][1]/dense.eden;*/
			phycon.te03)*dense.xIonDense[ipHYDROGEN][1]/dense.EdenHCorr;

		csp02 = MIN2(2.46e-7*phycon.te32/phycon.te10,5.21e-5*phycon.te/
			/* >>chng 05 jul 05, eden to dense.EdenHCorr */
			/*phycon.te20)*dense.xIonDense[ipHYDROGEN][1]/dense.eden;*/
			phycon.te20)*dense.xIonDense[ipHYDROGEN][1]/dense.EdenHCorr;
	}

	csp12 = MIN2(2.35e-6*phycon.te*phycon.te05*phycon.te01,3.98e-5*
		/* >>chng 05 jul 05, eden to dense.EdenHCorr */
		/*phycon.te20)*dense.xIonDense[ipHYDROGEN][1]/dense.eden;*/
		/*phycon.te70/phycon.te01)*dense.xIonDense[ipHYDROGEN][1]/dense.eden;*/
		phycon.te70/phycon.te01)*dense.xIonDense[ipHYDROGEN][1]/dense.EdenHCorr;

	return;
}

/*"oii_cs" calculates collision strengths for O II.
 * cs variables are named based of the transition they represent (lower level first).
 * Transitions where triplet states are considered to be a single state are
 * denoted by omitting the J value from the triplet level
 * ex. (3P(J=2,1,0) -> 1D2 = oi_cs3P1D2
 * When the J is a half integer, only the numerator is given.
 */

void oii_cs(double& oii_cs4S32D5,
	    double& oii_cs4S32D3,
	    double& oii_cs2D52D3,
	    double& oii_cs4S32P3,
	    double& oii_cs2D52P3,
	    double& oii_cs2D32P3,
	    double& oii_cs4S32P1,
	    double& oii_cs2D52P1,
	    double& oii_cs2D32P1,
	    double& oii_cs2P32P1,
	    double& oii_cs4S34P)
{
	DEBUG_ENTRY( "oii_cs()" );
	/* >>refer	o2	cs	Kisielius, R., Storey, P. J., Ferland, G. J., & Keenan, F. P.,
	 * >>refercon	2009, MNRAS, 397, 903 */
	oii_cs4S32D5 = (realnum)(0.7776*(phycon.te007*phycon.te0005*phycon.te0001));
	oii_cs4S32D3 = (realnum)(0.5643/phycon.te002);
	oii_cs2D52D3 = (realnum)(2.2277/(phycon.te07/(phycon.te003*phycon.te0001)));
	oii_cs4S32P3 = (realnum)(0.2004*phycon.te02*(phycon.te007/phycon.te0004));
	oii_cs2D52P3 = (realnum)(0.6211*phycon.te03*phycon.te004*phycon.te0002);
	oii_cs2D32P3 = (realnum)(0.3159*phycon.te04*(phycon.te004/phycon.te0004));
	oii_cs4S32P1 = (realnum)(0.1112*(phycon.te02/(phycon.te001*phycon.te0004)));
	oii_cs2D52P1 = (realnum)(0.2337*phycon.te04*phycon.te0004);
	oii_cs2D32P1 = (realnum)(0.2226*phycon.te04*phycon.te003*phycon.te0001);
	oii_cs2P32P1 = (realnum)(0.1943*phycon.te04*(phycon.te002/phycon.te0004));

	/* O II 4S34P 833.9A, CS
	 * >>refer	o2	cs	McLaughlin, B.M., & Bell, K.L. 1993, ApJ, 408, 753 */
	/* >>chng 01 aug 10, turn this line back on - no reason given for turning it off. */
	oii_cs4S34P = 0.355*phycon.te10*phycon.te10*phycon.te003*phycon.te001*phycon.te001;

	return;
}

/*"oiii_cs" calculates collision strengths for O III.
 * cs variables are named based of the transition they represent (lower level first). 
 * Transitions where triplet states are considered to be a single state are
 * denoted by omitting the J value from the triplet level
 * ex. (3P(J=2,1,0) -> 1D2 = oi_cs3P1D2    */
void oiii_cs(double& oiii_cs3P25S2,
	    double& oiii_cs3P15S2,
	    double& oiii_cs3P05S2,
	    double& oiii_cs3P1D2,
	    double& oiii_cs1D21S0,
	    double& oiii_cs3P1S0,
	    double& oiii_cs3P03P1,
	    double& oiii_cs3P13P2,
	    double &oiii_cs3P03P2,
	    double& oiii_cs3P3D)
{
    DEBUG_ENTRY( "oiii_cs()" );

    /* O III 1666, 3P25S2 crit den=2.6+10, A from
	 * >>refer	all	cs	Mendoza, C. 1982, in Planetary Nebulae, IAU Symp No. 103,
	 * >>refercon ed by D.R. Flower, (D. Reidel: Holland), 143
	 * c.s.
	 * >>referold	o3	cs	Lennon, D.J. Burke, V.M. 1994, A&AS, 103, 273 */
	/* >>chng 06 jun 30- Humeshkar Nemala */
	/* >>refer	o3	cs	Aggarwal,K.M. & Keenan,F. P.1999,ApJS,123,311*/
	/*Data available in the temperature range 2500 K to 2E6 K*/
	/*fits at temperatures below 30000K and at temperatures above 30000K*/
	if(phycon.te < 30000)
		oiii_cs3P25S2 = (realnum)(0.2519*(phycon.te07*phycon.te02*phycon.te007*phycon.te001*phycon.te0002));
	else
		oiii_cs3P25S2 = (realnum)(6.166388*(1/(phycon.te20*phycon.te01*phycon.te002)));


	/* O III 1661, 3P15S2
	 *  >>refer	o3	cs	Aggarwal,K.M. & Keenan,F. P.1999,ApJS,123,311*/
	/*Data available in the temperature range 2500 K to 2E6 K*/
	/*fits at temperatures below 30000K and at temperatures above 30000K*/
	if(phycon.te < 30000)
		oiii_cs3P15S2 = (realnum)((0.1511)*(phycon.te07*phycon.te02*phycon.te007*phycon.te001*phycon.te0002));
	else
		oiii_cs3P15S2 = (realnum)((3.668474)*(1/(phycon.te20*phycon.te01*phycon.te001*phycon.te0002)));

	/* O III 3P0-5S^o2*/
	/*>>chng 06 jun 30- Humeshkar Nemala*/
	/*>>refer	o3	cs	Aggarwal,K.M. & Keenan,F. P.1999,ApJS,123,311*/
	/*Data available in the temperature range 2500 K to 2E6 K*/
	/*fits at temperatures below 30000K and at temperatures above 30000K*/
	if(phycon.te < 30000)
		oiii_cs3P05S2 = (realnum)(0.0504*((phycon.te10/phycon.te01)*phycon.te005*phycon.te003*phycon.te0002));
	else
		oiii_cs3P05S2 = (realnum)(1.223633/(phycon.te20*phycon.te01*phycon.te001*phycon.te0002));


	/* o iii 5007+4959, As 96 NIST */
	/*The cs of the transitions 3P0,1,2 to 1D2 are added together to give oiii_cs3P1D2 */
	/*the cs of the transition 1D2-1S0 is mentioned as oiii_cs1D21S0*/
	/*The cs of the transitions 3P0,1,2 to 1S0 are added together to give oiii_cs3P1S0*/
	/* >>chng 01 may 04, As updated to
	 * >>referold	o3	as	Storey, P.J., & Zeippen, C.J., 2000, MNRAS, 312, 813-816,
	 * changed by 10 percent! */
	/* >>chng 10 feb 11 aeff = 0.027205 + oxy.d5007r;
	 * >>refer	o3	as	Froese Fischer, C., & Tachiev, G. 2004, At. Data Nucl. Data Tables, 87, 1
	 * term oxy.d5007r is photoioniozation loss */
	/* following data used in routine that deduces OIII temp from spectrum
	 *
	 * >>refer	o3	cs	Lennon, D.J. Burke, V.M. 1994, A&AS, 103, 273
	 * IP paper Y(ki) differ significantly from those
	 * calculated by
	 * >>refer	o3	cs	Burke, V.M., Lennon, D.J., & Seaton, M.J. 1989, MNRAS, 236, 353
	 * especially for 1D2-1S0.
	 * BLS89 is adopted for 1D2-1S0 and LB94 for 3P2,1-1D2.
	 * NB!!  these cs's must be kept parallel with those in Peimbert analysis */

	/* >>refer	o3	cs	Aggarwal,K.M. & Keenan,F. P.1999,ApJS,123,311*/
	/* >>chng 06 jul 24- Humeshkar Nemala */
	/*cs are available over two temperature ranges: below 30000K and above 30000K*/
	/*The old values seemed to saturate at around 100,000K.The new values are around 
	10-15% different from the old values*/
	
	if(phycon.te < 3E4)
	{
		oiii_cs3P1D2 = (realnum)(0.9062*(phycon.te10/phycon.te002));
		oiii_cs1D21S0 = (realnum)(0.0995*phycon.te10*phycon.te07*(phycon.te007/phycon.te0002));
		oiii_cs3P1S0 = (realnum)(0.1237*phycon.te07*phycon.te02*phycon.te005*phycon.te0005);

	}
	else if(phycon.te < 6E4)
	{	
		oiii_cs3P1D2 = (realnum)(1.710073*(phycon.te04/phycon.te004)*phycon.te0004);
		oiii_cs3P1S0 = (realnum)(0.1963109*phycon.te05*phycon.te0007);
		oiii_cs1D21S0 = (realnum)(0.781266/(phycon.te02*phycon.te003*phycon.te0001));
	}
	else
	{
		oiii_cs3P1D2 = (realnum)(6.452638/((phycon.te10/phycon.te02)*phycon.te004*phycon.te0003));
		oiii_cs3P1S0 = (realnum)(0.841578/(phycon.te07*phycon.te01*(phycon.te002/phycon.te0004)));
		oiii_cs1D21S0 = (realnum)(0.781266/(phycon.te02*phycon.te003*phycon.te0001));
	}
	
	/* >>chng 10 feb 11  chng from a31 = 0.215634; a32 = 1.71;
	* >>refer	o3	as	Froese Fischer, C., & Tachiev, G. 2004, At. Data Nucl. Data Tables, 87, 1
	*/
	/* O III IR lines, col str from iron project,
	 * >>referold	o3	cs	Lennon, D.J. Burke, V.M. 1994, A&AS, 103, 273 */
	/*>>refer	o3	cs	Aggarwal,K.M. & Keenan,F. P.1999,ApJS,123,311*/
	/*88 microns refers to the 3P0-3P1 transition*/
	if(phycon.te < 3E4)
		oiii_cs3P03P1 = (realnum)(0.3993*(phycon.te03/phycon.te001)*phycon.te0001);
	else if(phycon.te < 1E5)
		oiii_cs3P03P1 = (realnum)(0.245712*phycon.te07*phycon.te005*phycon.te001*phycon.te0002);
	else
		oiii_cs3P03P1 = (realnum)(1.310467/((phycon.te07/phycon.te001)*phycon.te0002));

	/*O III 52 microns refers to the 3P1-3P2 transition*/
	if(phycon.te < 3E4)
		oiii_cs3P13P2 = (realnum)(0.7812*(phycon.te05/phycon.te0001));
	else if(phycon.te < 1.2E5)
		oiii_cs3P13P2 = (realnum)(0.516157*phycon.te07*phycon.te02*phycon.te0001);
	else
		oiii_cs3P13P2 = (realnum)(4.038402/(phycon.te05*phycon.te03*phycon.te005*phycon.te0007*phycon.te0001));

	/*O III TauDummy refers to the 3P0-3P2 transition*/
	if(phycon.te < 3E4)
		oiii_cs3P03P2 = (realnum)(0.1337*phycon.te07*phycon.te002*phycon.te0002);
	else if(phycon.te < 1.2E5)
		oiii_cs3P03P2 = (realnum)(0.086446*phycon.te10*phycon.te01*phycon.te004*phycon.te0005);
	else
		oiii_cs3P03P2 = (realnum)(0.82031/(phycon.te07*phycon.te007*phycon.te0007*phycon.te0002));

	/* O III 834A, 3P3D CS */
	 /* >>referold	    o3	    cs	    Aggarwal, K.M., 1985 A&A 146, 149. */
	/* >>chng 06 jul 22- Humeshkar Nemala */
	/*>>refer	o3	cs	Aggarwal,K.M. & Keenan,F. P.1999,ApJS,123,311*/
	/*the cs of OIII 834A was obtained by summing the cs of the nine transitions:
	3P(0,1,2)-3D^o(1,2,3)*/
	if(phycon.te < 3E4)
		oiii_cs3P3D = (realnum)(4.74*phycon.te04*phycon.te0002);
	else
		oiii_cs3P3D = (realnum)(0.533*phycon.te20*phycon.te05*phycon.te002*phycon.te0001);
    return;
}
/*"oiv_cs" calculates collision strengths for O IV.
 * cs variables are named based of the transition they represent (lower level first).
 * Transitions where triplet states are considered to be a single state are
 * denoted by omitting the J value from the triplet level
 * ex. (3P(J=2,1,0) -> 1D2 = oi_cs3P1D2
 * For half integer Js, only the numerator is given.
 */
void oiv_cs(double& oiv_cs2P2D,double& oiv_cs2P12P3)
{
    	double Te_bounded,Te_bounded_log;

	DEBUG_ENTRY( "oiv_cs()" );
    /* these cs data only extend over a modest Temp range */
	Te_bounded = MAX2(phycon.te,4500.);
	Te_bounded = MIN2(Te_bounded,450000.);
	Te_bounded_log = log(Te_bounded);
	/* O IV 789A,  2P2D CS from
	 * >>refer	o4	cs	Zhang, H.L., Graziani, M., Pradhan, A.K. 1994, A&A, 283, 319 */
	oiv_cs2P2D = -3.0102462 + 109.22973/Te_bounded_log - 18666.357/Te_bounded;
	/* O IV 26 micron, CS
	 * >>referold	o4	cs	Blum, R.D., & Pradhan, A.K. 1992, ApJS 80, 425
	 * A=
	 * >>refer	o4	as	Brage, T., Judge, P.G., & Brekke, P. 1996, ApJ. 464, 1030 */
	/* >>chng 06 nov 08 - NPA.  Update collision strength to new data from:
	 * >>refer	o4	cs	Tayal, S. 2006, ApJS 166, 634
	 * Equation derived by using TableCurve, and goes to zero as
	 * T => 0 and T => infinity */
	oiv_cs2P12P3 = (realnum)(exp(3.265723 - 0.00014686984*phycon.alnte*phycon.sqrte
		 -22.035066/phycon.alnte));
	/* cs goes to zero at very high T, which will cause an assert in the
	 * n-level solver - don't let this happen */
	oiv_cs2P12P3 = MAX2( oiv_cs2P12P3 , 3.25e-2);
	return;
}
/*"ov_cs" calculates collision strengths for O V.
 * cs variables are named based of the transition they represent (lower level first). 
 * Transitions where triplet states are considered to be a single state are
 * denoted by omitting the J value from the triplet level
 * ex. (3P(J=2,1,0) -> 1D2 = oi_cs3P1D2               */
void ov_cs(double& ov_cs1S01P1,double& ov_cs1S03P)
{
    DEBUG_ENTRY( "ov_cs()" );
    /* O V 630, 1S01P1 CS from
	 * >>refer	o5	cs	Berrington, K.A., Burke, P.G., Dufton, P.L., Kingston, A.E. 1985,
	 * >>refercon	At. Data Nucl. Data Tables, 33, 195 */
	ov_cs1S01P1 = MIN2(4.0,1.317*phycon.te10/phycon.te03);
    /* O V 1218; 1S03P coll data from
	 * >>refer	o5	cs	Berrington, K.A., Burke, P.G., Dufton, P.L., Kingston, A.E. 1985,
	 * >>refercon	At. Data Nucl. Data Tables, 33, 195 */
	if( phycon.te <= 3.16e4 )
	{
		ov_cs1S03P = 3.224/(phycon.te10*phycon.te03*phycon.te03*phycon.te003);
	}
	else
	{
		ov_cs1S03P = 10.549/(phycon.te10*phycon.te10*phycon.te10/phycon.te03);
	}

    return;
}

