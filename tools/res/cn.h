/**
 * cooper-nathans calculation
 * @author tweber
 * @date 2013-2016
 * @license GPLv2
 *
 * @desc This is a reimplementation in C++ of the file rc_cnmat.m of the
 *		rescal5 package by Zinkin, McMorrow, Tennant, Farhi, and Wildes:
 *		http://www.ill.eu/en/instruments-support/computing-for-science/cs-software/all-software/matlab-ill/rescal-for-matlab/
 * @desc see: [cn67] M. J. Cooper and R. Nathans, Acta Cryst. 23, 357 (1967),
 * 		[ch73] N. J. Chesser and J. D. Axe, Acta Cryst. A 29, 160 (1973)
 */

#ifndef __TAKIN_CN_H__
#define __TAKIN_CN_H__

#include "defs.h"
#include "tlibs/math/neutrons.hpp"

namespace units = boost::units;
namespace codata = boost::units::si::constants::codata;

struct CNParams
{
	// monochromator
	tl::t_length_si<t_real_reso> mono_d;
	tl::t_angle_si<t_real_reso> mono_mosaic;
	t_real_reso dmono_sense = -1.;

	// analyser
	tl::t_length_si<t_real_reso> ana_d;
	tl::t_angle_si<t_real_reso> ana_mosaic;
	t_real_reso dana_sense = -1.;

	// sample
	tl::t_angle_si<t_real_reso> sample_mosaic;
	tl::t_length_si<t_real_reso> sample_lattice[3];
	tl::t_angle_si<t_real_reso> sample_angles[3];
	t_real_reso dsample_sense = 1.;

	// collimators
	tl::t_angle_si<t_real_reso> coll_h_pre_mono;
	tl::t_angle_si<t_real_reso> coll_h_pre_sample;
	tl::t_angle_si<t_real_reso> coll_h_post_sample;
	tl::t_angle_si<t_real_reso> coll_h_post_ana;
	tl::t_angle_si<t_real_reso> coll_v_pre_mono;
	tl::t_angle_si<t_real_reso> coll_v_pre_sample;
	tl::t_angle_si<t_real_reso> coll_v_post_sample;
	tl::t_angle_si<t_real_reso> coll_v_post_ana;

	tl::t_wavenumber_si<t_real_reso> ki, kf, Q;
	tl::t_energy_si<t_real_reso> E;

	tl::t_angle_si<t_real_reso> thetaa, thetam;
	tl::t_angle_si<t_real_reso> twotheta;

	tl::t_angle_si<t_real_reso> angle_ki_Q;
	tl::t_angle_si<t_real_reso> angle_kf_Q;

	// resolution volume stuff
	t_real_reso dmono_refl;
	t_real_reso dana_effic;

	bool bCalcR0 = 1;
};

extern ResoResults calc_cn(const CNParams& cn);

#endif
