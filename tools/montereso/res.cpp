/**
 * Monte-Carlo resolution calculation
 *
 * @author Tobias Weber
 * @date July 2012, Sep. 2014
 * @license GPLv2
 *
 * (based on Mcstas' mcresplot.pl perl program (www.mcstas.org)
 *  and the rescal5 matlab program)
 */

#include "res.h"

#include "tlibs/math/neutrons.hpp"
#include "tlibs/log/log.h"

#include <algorithm>
#include <boost/algorithm/minmax_element.hpp>
#include <fstream>
#include <ios>
#include <limits>
#include <memory>

using namespace ublas;


/*
 * this function tries to be a 1:1 C++ reimplementation of the Perl function
 * 'read_mcstas_res' of the McStas 'mcresplot' program
 */
Resolution calc_res(const std::vector<vector<t_real_reso>>& Q_vec,
	const vector<t_real_reso>& Q_avg, const std::vector<t_real_reso>* pp_vec)
{
	vector<t_real_reso> Q_dir = tl::make_vec({Q_avg[0], Q_avg[1], Q_avg[2]});
	Q_dir = Q_dir / norm_2(Q_dir);

	vector<t_real_reso> vecUp = tl::make_vec({0., 0., 1.});
	vector<t_real_reso> Q_perp = tl::cross_3(vecUp, Q_dir);
	vecUp = tl::cross_3(Q_dir, Q_perp);

	Q_perp[0]=Q_dir[1]; Q_perp[1]=Q_dir[0]; Q_perp[2]=Q_dir[2];


	matrix<t_real_reso> trafo(4, 4);
	trafo(0,0)=Q_dir[0];	trafo(1,0)=Q_dir[1];	trafo(2,0)=Q_dir[2];
	trafo(0,1)=Q_perp[0];	trafo(1,1)=Q_perp[1];	trafo(2,1)=Q_perp[2];
	trafo(0,2)=vecUp[0];	trafo(1,2)=vecUp[1];	trafo(2,2)=vecUp[2];

	trafo(3,0)=trafo(0,3)=trafo(3,1)=trafo(1,3)=trafo(3,2)=trafo(2,3) = 0.;
	trafo(3,3) = 1.;

	//std::cout << "trafo = " << trafo << std::endl;

	Resolution reso;
	reso.Q_avg = prod(trans(trafo), Q_avg);

	matrix<t_real_reso>& res = reso.res;
	matrix<t_real_reso>& cov = reso.cov;
	res.resize(4,4,0);
	cov.resize(4,4,0);

	cov = tl::covariance(Q_vec, pp_vec);
	cov = tl::transform<matrix<t_real_reso>>(cov, trafo, true);

	tl::log_info("Covariance matrix: ", cov);
	reso.bHasRes = tl::inverse(cov, res);

	if(reso.bHasRes)
	{
		reso.dQ.resize(4, 0);
		for(int iQ=0; iQ<4; ++iQ)
			reso.dQ[iQ] = tl::SIGMA2HWHM/sqrt(res(iQ,iQ));

		tl::log_info("Resolution matrix: ", res);

		const vector<t_real_reso>& dQ = reso.dQ;
		const vector<t_real_reso>& Q_avg = reso.Q_avg;

		std::ostringstream ostrVals;
		ostrVals << "Gaussian HWHM values: ";
		std::copy(dQ.begin(), dQ.end(), std::ostream_iterator<t_real_reso>(ostrVals, ", "));

		std::ostringstream ostrElli;
		ostrElli << "Ellipsoid offsets: ";
		std::copy(Q_avg.begin(), Q_avg.end(), std::ostream_iterator<t_real_reso>(ostrElli, ", "));

		tl::log_info(ostrVals.str());
		tl::log_info(ostrElli.str());
	}
	else
	{
		tl::log_err("Covariance matrix could not be inverted!");
	}

	return reso;
}


Resolution calc_res(unsigned int uiLen,
	const t_real_reso *_Q_x, const t_real_reso *_Q_y, const t_real_reso *_Q_z,
	const t_real_reso *_E)
{
	vector<t_real_reso> Q_avg(4);
	Q_avg[0] = Q_avg[1] = Q_avg[2] = Q_avg[3] = 0.;

	std::vector<vector<t_real_reso>> Q_vec;
	Q_vec.reserve(uiLen);

	for(unsigned int uiRow=0; uiRow<uiLen; ++uiRow)
	{
		vector<t_real_reso> Q(4);

		Q[0] = _Q_x[uiRow];
		Q[1] = _Q_y[uiRow];
		Q[2] = _Q_z[uiRow];
		Q[3] = _E[uiRow];

		Q_avg += Q;

		Q_vec.push_back(Q);
	}
	Q_avg /= t_real_reso(uiLen);
	tl::log_info("Average Q vector: ", Q_avg);

	return calc_res(Q_vec, Q_avg);
}


/*
 * this function tries to be a 1:1 C++ reimplementation of the Perl function
 * 'read_mcstas_res' of the McStas 'mcresplot' program
 */
Resolution calc_res(unsigned int uiLen,
	const t_real_reso *_ki_x, const t_real_reso *_ki_y, const t_real_reso *_ki_z,
	const t_real_reso *_kf_x, const t_real_reso *_kf_y, const t_real_reso *_kf_z,
	const t_real_reso *_p_i, const t_real_reso *_p_f)
{
	tl::log_info("Calculating resolution...");

	std::vector<vector<t_real_reso>> Q_vec;
	std::vector<t_real_reso> p_vec;
	Q_vec.reserve(uiLen);
	p_vec.reserve(uiLen);

	std::unique_ptr<t_real_reso[]> ptr_dE_vec(new t_real_reso[uiLen]);

	t_real_reso *dE_vec = ptr_dE_vec.get();


	const t_real_reso pi_max = _p_i ? *std::max_element(_p_i, _p_i+uiLen) : 1.;
	const t_real_reso pf_max = _p_f ? *std::max_element(_p_f, _p_f+uiLen) : 1.;
	const t_real_reso p_max = fabs(pi_max*pf_max);


	vector<t_real_reso> Q_avg(4);
	Q_avg[0] = Q_avg[1] = Q_avg[2] = Q_avg[3] = 0.;

	t_real_reso p_sum = 0.;

	for(unsigned int uiRow=0; uiRow<uiLen; ++uiRow)
	{
		vector<t_real_reso> Q(3);
		t_real_reso p;

		p = (_p_i && _p_f) ? fabs(_p_i[uiRow]*_p_f[uiRow]) : 1.;
		p /= p_max;		// normalize p to 0..1
		p_sum += p;

		t_real_reso &dE = dE_vec[uiRow];

		vector<t_real_reso> ki(3), kf(3);
		ki[0]=_ki_x[uiRow]; ki[1]=_ki_y[uiRow]; ki[2]=_ki_z[uiRow];
		kf[0]=_kf_x[uiRow]; kf[1]=_kf_y[uiRow]; kf[2]=_kf_z[uiRow];

		Q = ki - kf;
		t_real_reso Ei = tl::KSQ2E * inner_prod(ki, ki);
		t_real_reso Ef = tl::KSQ2E * inner_prod(kf, kf);
		dE = Ei - Ef;

		// insert the energy into the Q vector
		Q.resize(4, true);
		Q[3] = dE;

		Q_avg += Q*p;

		Q_vec.push_back(Q);
		p_vec.push_back(p);
	}
	Q_avg /= p_sum;
	tl::log_info("Average Q vector: ", Q_avg);


	return calc_res(Q_vec, Q_avg, &p_vec);
}
