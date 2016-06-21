/*
 * loads reso settings
 * @author tweber
 * @date jul-2015
 * @license GPLv2
 */

#include "TASReso.h"
#include "tlibs/math/lattice.h"
#include "tlibs/file/prop.h"
#include "tlibs/log/log.h"
#include "tlibs/helper/thread.h"

#include <boost/units/io.hpp>


typedef t_real_reso t_real;

using t_vec = ublas::vector<t_real>;
using t_mat = ublas::matrix<t_real>;

static const auto angs = tl::get_one_angstrom<t_real>();
static const auto rads = tl::get_one_radian<t_real>();
static const auto meV = tl::get_one_meV<t_real>();
static const auto cm = tl::get_one_centimeter<t_real>();
static const auto meters = tl::get_one_meter<t_real>();
static const auto sec = tl::get_one_second<t_real>();

using wavenumber = tl::t_wavenumber_si<t_real>;


TASReso::TASReso()
{
	m_opts.bCenter = 0;
	m_opts.coords = McNeutronCoords::RLU;
}

TASReso::TASReso(const TASReso& res)
{
	operator=(res);
}

const TASReso& TASReso::operator=(const TASReso& res)
{
	this->m_algo = res.m_algo;
	this->m_foc = res.m_foc;
	this->m_opts = res.m_opts;
	this->m_reso = res.m_reso;
	this->m_tofreso = res.m_tofreso;
	this->m_res = res.m_res;
	this->m_bKiFix = res.m_bKiFix;
	this->m_dKFix = res.m_dKFix;

	return *this;
}


bool TASReso::LoadLattice(const char* pcXmlFile)
{
	const std::string strXmlRoot("taz/");

	tl::Prop<std::string> xml;
	if(!xml.Load(pcXmlFile, tl::PropType::XML))
	{
		tl::log_err("Cannot load crystal file \"", pcXmlFile, "\".");
		return false;
	}

	t_real a = xml.Query<t_real>((strXmlRoot + "sample/a").c_str(), 0.);
	t_real b = xml.Query<t_real>((strXmlRoot + "sample/b").c_str(), 0.);
	t_real c = xml.Query<t_real>((strXmlRoot + "sample/c").c_str(), 0.);
	t_real alpha = tl::d2r(xml.Query<t_real>((strXmlRoot + "sample/alpha").c_str(), 90.));
	t_real beta = tl::d2r(xml.Query<t_real>((strXmlRoot + "sample/beta").c_str(), 90.));
	t_real gamma = tl::d2r(xml.Query<t_real>((strXmlRoot + "sample/gamma").c_str(), 90.));

	t_real dPlaneX0 = xml.Query<t_real>((strXmlRoot + "plane/x0").c_str(), 1.);
	t_real dPlaneX1 = xml.Query<t_real>((strXmlRoot + "plane/x1").c_str(), 0.);
	t_real dPlaneX2 = xml.Query<t_real>((strXmlRoot + "plane/x2").c_str(), 0.);
	t_real dPlaneY0 = xml.Query<t_real>((strXmlRoot + "plane/y0").c_str(), 0.);
	t_real dPlaneY1 = xml.Query<t_real>((strXmlRoot + "plane/y1").c_str(), 1.);
	t_real dPlaneY2 = xml.Query<t_real>((strXmlRoot + "plane/y2").c_str(), 0.);

	t_vec vec1 = tl::make_vec({dPlaneX0, dPlaneX1, dPlaneX2});
	t_vec vec2 = tl::make_vec({dPlaneY0, dPlaneY1, dPlaneY2});

	if(!SetLattice(a, b, c, alpha, beta, gamma, vec1, vec2))
		return false;

	return true;
}

bool TASReso::LoadRes(const char* pcXmlFile)
{
	const std::string strXmlRoot("taz/");

	tl::Prop<std::string> xml;
	if(!xml.Load(pcXmlFile, tl::PropType::XML))
	{
		tl::log_err("Cannot load resolution file \"", pcXmlFile, "\".");
		return false;
	}

	// CN
	m_reso.mono_d = xml.Query<t_real>((strXmlRoot + "reso/mono_d").c_str(), 0.) * angs;
	m_reso.mono_mosaic = tl::m2r(xml.Query<t_real>((strXmlRoot + "reso/mono_mosaic").c_str(), 0.)) * rads;
	m_reso.ana_d = xml.Query<t_real>((strXmlRoot + "reso/ana_d").c_str(), 0.) * angs;
	m_reso.ana_mosaic = tl::m2r(xml.Query<t_real>((strXmlRoot + "reso/ana_mosaic").c_str(), 0.)) * rads;
	m_reso.sample_mosaic = tl::m2r(xml.Query<t_real>((strXmlRoot + "reso/sample_mosaic").c_str(), 0.)) * rads;

	m_reso.coll_h_pre_mono = tl::m2r(xml.Query<t_real>((strXmlRoot + "reso/h_coll_mono").c_str(), 0.)) * rads;
	m_reso.coll_h_pre_sample = tl::m2r(xml.Query<t_real>((strXmlRoot + "reso/h_coll_before_sample").c_str(), 0.)) * rads;
	m_reso.coll_h_post_sample = tl::m2r(xml.Query<t_real>((strXmlRoot + "reso/h_coll_after_sample").c_str(), 0.)) * rads;
	m_reso.coll_h_post_ana = tl::m2r(xml.Query<t_real>((strXmlRoot + "reso/h_coll_ana").c_str(), 0.)) * rads;

	m_reso.coll_v_pre_mono = tl::m2r(xml.Query<t_real>((strXmlRoot + "reso/v_coll_mono").c_str(), 0.)) * rads;
	m_reso.coll_v_pre_sample = tl::m2r(xml.Query<t_real>((strXmlRoot + "reso/v_coll_before_sample").c_str(), 0.)) * rads;
	m_reso.coll_v_post_sample = tl::m2r(xml.Query<t_real>((strXmlRoot + "reso/v_coll_after_sample").c_str(), 0.)) * rads;
	m_reso.coll_v_post_ana = tl::m2r(xml.Query<t_real>((strXmlRoot + "reso/v_coll_ana").c_str(), 0.)) * rads;

	m_reso.dmono_refl = xml.Query<t_real>((strXmlRoot + "reso/mono_refl").c_str(), 0.);
	m_reso.dana_effic = xml.Query<t_real>((strXmlRoot + "reso/ana_effic").c_str(), 0.);

	m_reso.dmono_sense = (xml.Query<int>((strXmlRoot+"reso/mono_scatter_sense").c_str(), 0) ? +1. : -1.);
	m_reso.dana_sense = (xml.Query<int>((strXmlRoot+"reso/ana_scatter_sense").c_str(), 0) ? +1. : -1.);
	m_reso.dsample_sense = (xml.Query<int>((strXmlRoot+"reso/sample_scatter_sense").c_str(), 1) ? +1. : -1.);


	// Pop
	m_reso.mono_w = xml.Query<t_real>((strXmlRoot + "reso/pop_mono_w").c_str(), 0.)*cm;
	m_reso.mono_h = xml.Query<t_real>((strXmlRoot + "reso/pop_mono_h").c_str(), 0.)*cm;
	m_reso.mono_thick = xml.Query<t_real>((strXmlRoot + "reso/pop_mono_thick").c_str(), 0.)*cm;
	m_reso.mono_curvh = xml.Query<t_real>((strXmlRoot + "reso/pop_mono_curvh").c_str(), 0.)*cm;
	m_reso.mono_curvv = xml.Query<t_real>((strXmlRoot + "reso/pop_mono_curvv").c_str(), 0.)*cm;
	m_reso.bMonoIsCurvedH = (xml.Query<int>((strXmlRoot + "reso/pop_mono_use_curvh").c_str(), 0) != 0);
	m_reso.bMonoIsCurvedV = (xml.Query<int>((strXmlRoot + "reso/pop_mono_use_curvv").c_str(), 0) != 0);
	m_reso.bMonoIsOptimallyCurvedH = (xml.Query<int>((strXmlRoot + "reso/pop_mono_use_curvh_opt").c_str(), 1) != 0);
	m_reso.bMonoIsOptimallyCurvedV = (xml.Query<int>((strXmlRoot + "reso/pop_mono_use_curvv_opt").c_str(), 1) != 0);

	m_reso.ana_w = xml.Query<t_real>((strXmlRoot + "reso/pop_ana_w").c_str(), 0.)*cm;
	m_reso.ana_h = xml.Query<t_real>((strXmlRoot + "reso/pop_ana_h").c_str(), 0.)*cm;
	m_reso.ana_thick = xml.Query<t_real>((strXmlRoot + "reso/pop_ana_thick").c_str(), 0.)*cm;
	m_reso.ana_curvh = xml.Query<t_real>((strXmlRoot + "reso/pop_ana_curvh").c_str(), 0.)*cm;
	m_reso.ana_curvv = xml.Query<t_real>((strXmlRoot + "reso/pop_ana_curvv").c_str(), 0.)*cm;
	m_reso.bAnaIsCurvedH = (xml.Query<int>((strXmlRoot + "reso/pop_ana_use_curvh").c_str(), 0) != 0);
	m_reso.bAnaIsCurvedV = (xml.Query<int>((strXmlRoot + "reso/pop_ana_use_curvv").c_str(), 0) != 0);
	m_reso.bAnaIsOptimallyCurvedH = (xml.Query<int>((strXmlRoot + "reso/pop_ana_use_curvh_opt").c_str(), 1) != 0);
	m_reso.bAnaIsOptimallyCurvedV = (xml.Query<int>((strXmlRoot + "reso/pop_ana_use_curvv_opt").c_str(), 1) != 0);

	m_reso.bSampleCub = (xml.Query<int>((strXmlRoot + "reso/pop_sample_cuboid").c_str(), 0) != 0);
	m_reso.sample_w_q = xml.Query<t_real>((strXmlRoot + "reso/pop_sample_wq").c_str(), 0.)*cm;
	m_reso.sample_w_perpq = xml.Query<t_real>((strXmlRoot + "reso/pop_sampe_wperpq").c_str(), 0.)*cm;
	m_reso.sample_h = xml.Query<t_real>((strXmlRoot + "reso/pop_sample_h").c_str(), 0.)*cm;

	m_reso.bSrcRect = (xml.Query<int>((strXmlRoot + "reso/pop_source_rect").c_str(), 0) != 0);
	m_reso.src_w = xml.Query<t_real>((strXmlRoot + "reso/pop_src_w").c_str(), 0.)*cm;
	m_reso.src_h = xml.Query<t_real>((strXmlRoot + "reso/pop_src_h").c_str(), 0.)*cm;

	m_reso.bDetRect = (xml.Query<int>((strXmlRoot + "reso/pop_det_rect").c_str(), 0) != 0);
	m_reso.det_w = xml.Query<t_real>((strXmlRoot + "reso/pop_det_w").c_str(), 0.)*cm;
	m_reso.det_h = xml.Query<t_real>((strXmlRoot + "reso/pop_det_h").c_str(), 0.)*cm;

	m_reso.bGuide = (xml.Query<int>((strXmlRoot + "reso/use_guide").c_str(), 0) != 0);
	m_reso.guide_div_h = tl::m2r(xml.Query<t_real>((strXmlRoot + "reso/pop_guide_divh").c_str(), 0.)) * rads;
	m_reso.guide_div_v = tl::m2r(xml.Query<t_real>((strXmlRoot + "reso/pop_guide_divv").c_str(), 0.)) * rads;

	m_reso.dist_mono_sample = xml.Query<t_real>((strXmlRoot + "reso/pop_dist_mono_sample").c_str(), 0.)*cm;
	m_reso.dist_sample_ana = xml.Query<t_real>((strXmlRoot + "reso/pop_dist_sample_ana").c_str(), 0.)*cm;
	m_reso.dist_ana_det = xml.Query<t_real>((strXmlRoot + "reso/pop_dist_ana_det").c_str(), 0.)*cm;
	m_reso.dist_src_mono = xml.Query<t_real>((strXmlRoot + "reso/pop_dist_src_mono").c_str(), 0.)*cm;


	// Eck
	m_reso.mono_mosaic_v = tl::m2r(xml.Query<t_real>((strXmlRoot + "reso/eck_mono_mosaic_v").c_str(), 0.)) * rads;
	m_reso.ana_mosaic_v = tl::m2r(xml.Query<t_real>((strXmlRoot + "reso/eck_ana_mosaic_v").c_str(), 0.)) * rads;
	m_reso.pos_x = xml.Query<t_real>((strXmlRoot + "reso/eck_sample_pos_x").c_str(), 0.)*cm;
	m_reso.pos_y = xml.Query<t_real>((strXmlRoot + "reso/eck_sample_pos_y").c_str(), 0.)*cm;
	m_reso.pos_z = xml.Query<t_real>((strXmlRoot + "reso/eck_sample_pos_z").c_str(), 0.)*cm;

	// TODO
	m_reso.mono_numtiles_h = 1;
	m_reso.mono_numtiles_v = 1;
	m_reso.ana_numtiles_h = 1;
	m_reso.ana_numtiles_v = 1;


	// TOF
	m_tofreso.len_pulse_mono = xml.Query<t_real>((strXmlRoot + "reso/viol_dist_pulse_mono").c_str(), 0.) * cm;
	m_tofreso.len_mono_sample = xml.Query<t_real>((strXmlRoot + "reso/viol_dist_mono_sample").c_str(), 0.) * cm;
	m_tofreso.len_sample_det = xml.Query<t_real>((strXmlRoot + "reso/viol_dist_sample_det").c_str(), 0.) * cm;

	m_tofreso.sig_len_pulse_mono = xml.Query<t_real>((strXmlRoot + "reso/viol_dist_pulse_mono_sig").c_str(), 0.) * cm;
	m_tofreso.sig_len_mono_sample = xml.Query<t_real>((strXmlRoot + "reso/viol_dist_mono_sample_sig").c_str(), 0.) * cm;
	m_tofreso.sig_len_sample_det = xml.Query<t_real>((strXmlRoot + "reso/viol_dist_sample_det_sig").c_str(), 0.) * cm;

	m_tofreso.sig_pulse = (xml.Query<t_real>((strXmlRoot + "reso/viol_time_pulse_sig").c_str(), 0.) * 1e-6) * sec;
	m_tofreso.sig_mono = (xml.Query<t_real>((strXmlRoot + "reso/viol_time_mono_sig").c_str(), 0.) * 1e-6) * sec;
	m_tofreso.sig_det = (xml.Query<t_real>((strXmlRoot + "reso/viol_time_det_sig").c_str(), 0.) * 1e-6) * sec;

	m_tofreso.twotheta_i = tl::d2r(xml.Query<t_real>((strXmlRoot + "reso/viol_angle_tt_i").c_str(), 0.)) * rads;
	m_tofreso.angle_outplane_i = tl::d2r(xml.Query<t_real>((strXmlRoot + "reso/viol_angle_ph_i").c_str(), 0.)) * rads;
	m_tofreso.angle_outplane_f = tl::d2r(xml.Query<t_real>((strXmlRoot + "reso/viol_angle_ph_f").c_str(), 0.)) * rads;

	m_tofreso.sig_twotheta_i = tl::d2r(xml.Query<t_real>((strXmlRoot + "reso/viol_angle_tt_i_sig").c_str(), 0.)) * rads;
	m_tofreso.sig_twotheta_f = tl::d2r(xml.Query<t_real>((strXmlRoot + "reso/viol_angle_tt_f_sig").c_str(), 0.)) * rads;
	m_tofreso.sig_outplane_i = tl::d2r(xml.Query<t_real>((strXmlRoot + "reso/viol_angle_ph_i_sig").c_str(), 0.)) * rads;
	m_tofreso.sig_outplane_f = tl::d2r(xml.Query<t_real>((strXmlRoot + "reso/viol_angle_ph_f_sig").c_str(), 0.)) * rads;


	m_algo = ResoAlgo(xml.Query<int>((strXmlRoot + "reso/algo").c_str(), 0));


	// preliminary position
	m_tofreso.ki = m_reso.ki = xml.Query<t_real>((strXmlRoot + "reso/ki").c_str(), 0.) / angs;
	m_tofreso.kf = m_reso.kf = xml.Query<t_real>((strXmlRoot + "reso/kf").c_str(), 0.) / angs;
	m_tofreso.E = m_reso.E = xml.Query<t_real>((strXmlRoot + "reso/E").c_str(), 0.) * meV;
	m_tofreso.Q = m_reso.Q = xml.Query<t_real>((strXmlRoot + "reso/Q").c_str(), 0.) / angs;

	m_dKFix = m_bKiFix ? m_reso.ki*angs : m_reso.kf*angs;
	return true;
}


bool TASReso::SetLattice(t_real a, t_real b, t_real c,
	t_real alpha, t_real beta, t_real gamma,
	const t_vec& vec1, const t_vec& vec2)
{
	tl::Lattice<t_real> latt(a, b, c, alpha, beta, gamma);

	t_mat matB = tl::get_B(latt, 1);
	t_mat matU = tl::get_U(vec1, vec2, &matB);
	t_mat matUB = ublas::prod(matU, matB);

	t_mat matBinv, matUinv;
	bool bHasB = tl::inverse(matB, matBinv);
	bool bHasU = tl::inverse(matU, matUinv);
	t_mat matUBinv = ublas::prod(matBinv, matUinv);

	bool bHasUB = bHasB && bHasU;

	if(!bHasUB)
	{
		tl::log_err("Cannot invert UB matrix");
		return false;
	}


	m_opts.matU = matU;
	m_opts.matB = matB;
	m_opts.matUB = matUB;
	m_opts.matUinv = matUinv;
	m_opts.matBinv = matBinv;
	m_opts.matUBinv = matUBinv;

	ublas::matrix<t_real>* pMats[] = {&m_opts.matU, &m_opts.matB, &m_opts.matUB, 
		&m_opts.matUinv, &m_opts.matBinv, &m_opts.matUBinv};

	for(ublas::matrix<t_real> *pMat : pMats)
	{
		pMat->resize(4,4,1);

		for(int i0=0; i0<3; ++i0)
			(*pMat)(i0,3) = (*pMat)(3,i0) = 0.;
		(*pMat)(3,3) = 1.;
	}

	return true;
}

bool TASReso::SetHKLE(t_real h, t_real k, t_real l, t_real E)
{
	if(m_opts.matUB.size1() < 3 || m_opts.matUB.size2() < 3)
	{
		const char* pcErr = "Invalid UB matrix.";
		tl::log_err(pcErr);
		m_res.strErr = pcErr;
		m_res.bOk = false;
		return false;
	}

	t_vec vecHKLE;
	if(m_opts.matUB.size1() == 3)
		vecHKLE = tl::make_vec({h, k, l});
	else
		vecHKLE = tl::make_vec({h, k, l, E});

	t_vec vecQ = ublas::prod(m_opts.matUB, vecHKLE);
	if(vecQ.size() > 3)
		vecQ.resize(3, true);

	m_tofreso.Q = m_reso.Q = ublas::norm_2(vecQ) / angs;
	m_tofreso.E = m_reso.E = E * meV;

	wavenumber kother = tl::get_other_k(m_reso.E, m_dKFix/angs, m_bKiFix);
	if(m_bKiFix)
	{
		m_tofreso.ki = m_reso.ki = m_dKFix / angs;
		m_tofreso.kf = m_reso.kf = kother;
	}
	else
	{
		m_tofreso.ki = m_reso.ki = kother;
		m_tofreso.kf = m_reso.kf = m_dKFix / angs;
	}

	m_reso.thetam = units::abs(tl::get_mono_twotheta(m_reso.ki, m_reso.mono_d, /*m_reso.dmono_sense>=0.*/1)*t_real(0.5));
	m_reso.thetaa = units::abs(tl::get_mono_twotheta(m_reso.kf, m_reso.ana_d, /*m_reso.dana_sense>=0.*/1)*t_real(0.5));
	m_tofreso.twotheta = m_reso.twotheta = units::abs(tl::get_sample_twotheta(m_reso.ki, m_reso.kf, m_reso.Q, 1));

	m_tofreso.angle_ki_Q = m_reso.angle_ki_Q = tl::get_angle_ki_Q(m_reso.ki, m_reso.kf, m_reso.Q, /*m_reso.dsample_sense>=0.*/1);
	m_tofreso.angle_kf_Q = m_reso.angle_kf_Q = tl::get_angle_kf_Q(m_reso.ki, m_reso.kf, m_reso.Q, /*m_reso.dsample_sense>=0.*/1);


	if(m_foc == ResoFocus::FOC_NONE)
	{
		m_reso.bMonoIsCurvedH = m_reso.bMonoIsCurvedV = 0;
		m_reso.bAnaIsCurvedH = m_reso.bAnaIsCurvedV = 0;
	}
	else
	{
		m_reso.bMonoIsCurvedH = m_reso.bMonoIsOptimallyCurvedH =
			(unsigned(m_foc) & unsigned(ResoFocus::FOC_MONO_H));
		m_reso.bMonoIsCurvedV = m_reso.bMonoIsOptimallyCurvedV =
			(unsigned(m_foc) & unsigned(ResoFocus::FOC_MONO_V));
		m_reso.bAnaIsCurvedH = m_reso.bAnaIsOptimallyCurvedH =
			(unsigned(m_foc) & unsigned(ResoFocus::FOC_ANA_H));
		m_reso.bAnaIsCurvedV = m_reso.bAnaIsOptimallyCurvedV =
			(unsigned(m_foc) & unsigned(ResoFocus::FOC_ANA_V));


		// remove collimators
		/*if(m_reso.bMonoIsCurvedH)
		{
			m_reso.coll_h_pre_mono = 99999. * rads;
			m_reso.coll_h_pre_sample = 99999. * rads;
		}
		if(m_reso.bMonoIsCurvedV)
		{
			m_reso.coll_v_pre_mono = 99999. * rads;
			m_reso.coll_v_pre_sample = 99999. * rads;
		}
		if(m_reso.bAnaIsCurvedH)
		{
			m_reso.coll_h_post_sample = 99999. * rads;
			m_reso.coll_h_post_ana = 99999. * rads;
		}
		if(m_reso.bAnaIsCurvedV)
		{
			m_reso.coll_v_post_sample = 99999. * rads;
			m_reso.coll_v_post_ana = 99999. * rads;
		}*/
	}


	if(std::fabs(vecQ[2]) > tl::get_plane_dist_tolerance<t_real>())
	{
		tl::log_err("Position Q = (", h, " ", k, " ", l, "),",
			" E = ", E, " meV not in scattering plane.");

		m_res.strErr = "Not in scattering plane.";
		m_res.bOk = false;
		return false;
	}

	vecQ.resize(2, true);
	m_opts.dAngleQVec0 = -tl::vec_angle(vecQ);

	// calculate resolution at (hkl) and E
	if(m_algo == ResoAlgo::CN)
	{
		m_reso.bCalcR0 = false;
		m_res = calc_cn(m_reso);
	}
	else if(m_algo == ResoAlgo::POP)
	{
		//m_reso.bCalcR0 = true;
		m_res = calc_pop(m_reso);
	}
	else if(m_algo == ResoAlgo::ECK)
	{
		m_reso.bCalcR0 = true;
		m_res = calc_eck(m_reso);
	}
	else if(m_algo == ResoAlgo::VIOL)
	{
		m_res = calc_viol(m_tofreso);
	}
	else
	{
		const char* pcErr = "Unknown algorithm selected.";
		tl::log_err(pcErr);
		m_res.strErr = pcErr;
		m_res.bOk = false;
		return false;
	}

	if(!m_res.bOk)
	{
		tl::log_err("Error calculating resolution: ", m_res.strErr);
		tl::log_debug("R0: ", m_res.dR0);
		tl::log_debug("res: ", m_res.reso);
	}

	return m_res.bOk;
}

Ellipsoid4d<t_real> TASReso::GenerateMC(std::size_t iNum, std::vector<t_vec>& vecNeutrons) const
{
	Ellipsoid4d<t_real> ell4d = calc_res_ellipsoid4d<t_real>(m_res.reso, m_res.Q_avg);
	if(vecNeutrons.size() != iNum)
		vecNeutrons.resize(iNum);

	unsigned int iNumThreads = std::thread::hardware_concurrency();
	std::size_t iNumPerThread = iNum / iNumThreads;
	std::size_t iRemaining = iNum % iNumThreads;

	tl::ThreadPool<void()> tp(iNumThreads);
	for(unsigned iThread=0; iThread<iNumThreads; ++iThread)
	{
		std::vector<t_vec>::iterator iterBegin = vecNeutrons.begin() + iNumPerThread*iThread;
		std::size_t iNumNeutr = iNumPerThread;
		if(iThread == iNumThreads-1)
			iNumNeutr = iNumPerThread + iRemaining;

		tp.AddTask([iterBegin, iNumNeutr, this, &ell4d]()
			{ mc_neutrons(ell4d, iNumNeutr, this->m_opts, iterBegin); });
	}

	tp.StartTasks();

	auto& lstFut = tp.GetFutures();
	for(auto& fut : lstFut)
		fut.get();

	return ell4d;
}
