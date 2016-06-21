/*
 * MIEZE formulas
 * @author Tobias Weber
 * @date May 2012, 29-may-2013
 * @license GPLv2 or GPLv3
 */

#ifndef __MIEZE_FORMULAS__
#define __MIEZE_FORMULAS__

#include "neutrons.hpp"
#include "linalg.h"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>


namespace tl {

namespace ublas = boost::numeric::ublas;


//------------------------------------------------------------------------------
// MIEZE time (eq. 117 from [Keller, Golub, Gähler, 2000])
// tau = hbar * omega * Ls / (m*v^3)
template<class Sys, class Y>
t_time<Sys,Y> mieze_tau(const t_freq<Sys,Y>& fm,
	const t_length<Sys,Y>& Ls, const t_length<Sys,Y>& lam)
{
	t_velocity<Sys,Y> v = lam2p(lam) / co::m_n;
	return 2.*get_pi<Y>() * fm * Ls * co::hbar / (co::m_n * v*v*v);
}

template<class Sys, class Y>
t_freq<Sys,Y> mieze_tau_fm(const t_time<Sys,Y>& tau,
	const t_length<Sys,Y>& Ls, const t_length<Sys,Y>& lam)
{
	t_velocity<Sys,Y> v = lam2p(lam) / co::m_n;
	return tau / (2.*get_pi<Y>() * Ls * co::hbar) * (co::m_n * v*v*v);
}

template<class Sys, class Y>
t_length<Sys,Y> mieze_tau_Ls(const t_time<Sys,Y>& tau,
	const t_freq<Sys,Y>& fm, const t_length<Sys,Y>& lam)
{
	t_velocity<Sys,Y> v = lam2p(lam) / co::m_n;
	return tau / (2.*get_pi<Y>() * fm * co::hbar) * (co::m_n * v*v*v);
}

template<class Sys, class Y>
t_length<Sys,Y> mieze_tau_lam(const t_time<Sys,Y>& tau,
	const t_freq<Sys,Y>& fm, const t_length<Sys,Y>& Ls)
{
	t_velocity<Sys,Y> v;
	v = boost::units::root<3>(2.*get_pi<Y>() * fm * Ls * co::hbar / (tau * co::m_n));

	t_momentum<Sys,Y> p = v * co::m_n;
	return p2lam(p);
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// Larmor precession

// gamma*B = omega
template<class Sys, class Y=double>
t_freq<Sys,Y> larmor_om(const t_flux<Sys,Y>& B)
{
	return co::gamma_n * B;
}

template<class Sys, class Y=double>
t_flux<Sys,Y> larmor_B(const t_freq<Sys,Y>& om)
{
	return om/co::gamma_n;
}

/* omega = -gamma*B
 * omega*t = -gamma*B*t
 * phi = - gamma * B * l/v
 * B = -phi*v / (gamma*l)
 * phi = -pi  =>  B = pi*v / (gamma*l)
 */
template<class Sys, class Y=double>
t_flux<Sys,Y> larmor_field(const t_length<Sys,Y>& lam,
	const t_length<Sys,Y>& len,
	const t_angle<Sys,Y>& phi)
{
	t_velocity<Sys,Y> v = lam2p(lam) / co::m_n;
	t_freq<Sys,Y> om = -Y(phi/radians)*v/len;
	return om/co::gamma_n;
}

//------------------------------------------------------------------------------



//------------------------------------------------------------------------------
// MIEZE contrast reduction due to detector geometry
template<class Sys, class Y>
Y mieze_reduction_det(const t_length<Sys,Y>& lx, const t_length<Sys,Y>& ly,
	const t_length<Sys,Y>& xpos, const t_length<Sys,Y>& ypos,
	const t_length<Sys,Y>& Ls,
	const t_time<Sys,Y>& tau,
	const t_length<Sys,Y>& lam,
	const t_angle<Sys,Y>& central_phase,
	std::size_t iXPixels=128, std::size_t iYPixels=128,
	ublas::matrix<Y>* pPhaseMatrix=0)
{
	using namespace units;
	using namespace co;

	const t_mass<Sys,Y> mn = co::m_n;
	t_velocity<Sys,Y> v0 = lam2p(lam)/mn;

	const t_freq<Sys,Y> fM = mieze_tau_fm(tau, Ls, lam);
	const t_freq<Sys,Y> omegaM = 2.*get_pi<Y>()*fM;

	t_length<Sys,Y> lx_inc = lx / Y(iXPixels);   // pixel size
	t_length<Sys,Y> ly_inc = ly / Y(iYPixels);   // pixel size

	if(pPhaseMatrix)
		pPhaseMatrix->resize(iXPixels, iYPixels, false);

	t_length<Sys,Y> dx, dy;
	t_area<Sys,Y> int_red = 0.*meters*meters;


	std::size_t iX, iY;
	// integrate over detector
	for(dx=-lx/2., iX=0; dx<lx/2. && iX<iXPixels; dx+=lx_inc, ++iX)
	{
		for(dy=-ly/2., iY=0; dy<ly/2. && iY<iYPixels; dy+=ly_inc, ++iY)
		{
			t_length<Sys,Y> dx_new = dx-xpos;
			t_length<Sys,Y> dy_new = dy-ypos;

			t_length<Sys,Y> path_diff = sqrt(dx_new*dx_new + dy_new*dy_new + Ls*Ls) - Ls;

			// additional time needed for the neutron
			t_time<Sys,Y> dt = path_diff / v0;

			// additional phase
			Y phase = -omegaM * dt + central_phase/radians;
			phase = fmod(phase, 2.*get_pi<Y>());

			int_red += cos(phase/2.)*lx_inc*ly_inc;

			if(pPhaseMatrix)
			{
				phase += 2.*get_pi<Y>();
				phase = fmod(phase, 2.*get_pi<Y>());
				(*pPhaseMatrix)(iX, iY) = phase;
			}
		}
	}

	Y dreduction = int_red / (lx*ly);
	dreduction = fabs(dreduction);

	return dreduction;
}

// reduction factor due to detector thickness
template<class Sys, class Y>
Y mieze_reduction_det_d(const t_length<Sys,Y>& d,
	const t_freq<Sys,Y>& fM, const t_length<Sys,Y>& lam)
{
	using namespace units;
	using namespace co;

	const t_mass<Sys,Y> mn = co::m_n;
	t_velocity<Sys,Y> v0 = lam2p(lam)/mn;

	const t_freq<Sys,Y> omegaM = 2.*get_pi<Y>()*fM;

	t_length<Sys,Y> int_red = 0.*meters;
	int_red = std::sin(-0.5*omegaM/v0 * d)/(-0.5*omegaM/v0);

	Y dreduction = int_red / d;
	dreduction = fabs(dreduction);

	return dreduction;
}

//------------------------------------------------------------------------------




//------------------------------------------------------------------------------
// MIEZE contrast reduction due to sample geometry

// numerical approximation to the R_sample integral of formula (9) in [Brandl 11]
template<class Sys, class Y>
Y mieze_reduction_sample_cuboid(const t_length<Sys,Y>& len_x,
	const t_length<Sys,Y>& len_y, const t_length<Sys,Y>& len_z,
	const t_freq<Sys,Y>& fM,
	const t_length<Sys,Y>& lam,
	const t_angle<Sys,Y>& twotheta, const t_angle<Sys,Y>& theta_s,
	std::size_t ITERS=100)
{
	using namespace units;
	using namespace co;

	const t_freq<Sys,Y> omegaM = 2.*get_pi<Y>()*fM;
	t_velocity<Sys,Y> v = lam2p(lam)/co::m_n;

	ublas::vector<Y> ki(3);
	ki[0] = 0.;
	ki[1] = 0.;
	ki[2] = 1.;

	ublas::vector<Y> kf(3);
	kf[0] = sin(twotheta);
	kf[1] = 0.;
	kf[2] = cos(twotheta);

	ublas::vector<Y> q_dir = ki-kf;

	t_length<Sys,Y> dX = len_x / Y(ITERS);
	t_length<Sys,Y> dY = len_y / Y(ITERS);
	t_length<Sys,Y> dZ = len_z / Y(ITERS);

	t_length<Sys,Y> x, y, z;

	t_volume<Sys,Y> integral = 0.*meters*meters*meters;
	t_volume<Sys,Y> vol = 0.*meters*meters*meters;

	const Y stheta_s = sin(theta_s);
	const Y ctheta_s = cos(theta_s);

	for(x=-len_x/2.; x<len_x/2.; x+=dX)
		for(y=-len_y/2.; y<len_y/2.; y+=dY)
			for(z=-len_z/2.; z<len_z/2.; z+=dZ)
			{
				t_length<Sys,Y> pos[3];
				// rotate sample
				pos[0] = ctheta_s*x + stheta_s*z;
				pos[1] = y;
				pos[2] = -stheta_s*x + ctheta_s*z;

				t_length<Sys,Y> path_diff = q_dir[0]*pos[0] + q_dir[1]*pos[1] + q_dir[2]*pos[2];
				Y phase = omegaM * path_diff / v;

				t_volume<Sys,Y> func_det = dX*dY*dZ;

				vol += func_det;
				integral += func_det * cos(phase);
			}

	vol = len_x*len_y*len_z;
	return integral / vol;
}

// Scattering with extinction
template<class Sys, class Y>
Y mieze_reduction_sample_cuboid_extinction(const t_length<Sys,Y>& len_x,
	const t_length<Sys,Y>& len_y, const t_length<Sys,Y>& len_z,
	const t_length_inverse<Sys,Y>& mu,
	const t_freq<Sys,Y>& fM,
	const t_length<Sys,Y>& lam,
	const t_angle<Sys,Y>& twotheta,
	std::size_t ITERS=100)
{
	const Y SUBDIVS = Y(ITERS);

	using namespace co;

	t_frequency<Sys,Y> omegaM = 2.*get_pi<Y>()*fM;
	t_velocity<Sys,Y> v = lam2p(lam)/co::m_n;

	ublas::vector<Y> ki(3);
	ki[0] = 0.;
	ki[1] = 0.;
	ki[2] = 1.;

	ublas::vector<Y> kf(3);
	kf[0] = sin(twotheta);
	kf[1] = 0.;
	kf[2] = cos(twotheta);

	ublas::vector<Y> q_dir = ki-kf;


	t_length<Sys,Y> x, y, z;
	t_volume<Sys,Y> integral = 0.*meters*meters*meters;
	t_volume<Sys,Y> vol = 0.*meters*meters*meters;


	t_angle<Sys,Y> theta_s = twotheta/2. - get_pi<Y>()/2.*radians;
	const Y stheta_s = sin(theta_s);
	const Y ctheta_s = cos(theta_s);


	t_length<Sys,Y> zpath = 1./(mu*2.);	// reflexive: path taken twice
	//zpath /= ctheta_s;
	//if(zpath > len_z)
		zpath = len_z;

	t_length<Sys,Y> dX = len_x / SUBDIVS;
	t_length<Sys,Y> dY = len_y / SUBDIVS;
	t_length<Sys,Y> dZ = zpath / SUBDIVS;

	const t_volume<Sys,Y> func_det = dX*dY*dZ;

	for(x=-len_x/2.; x<len_x/2.; x+=dX)
		for(y=-len_y/2.; y<len_y/2.; y+=dY)
			for(z=0.*meters; z<zpath; z+=dZ)
			{
				//z = len_z/2.;
				t_length<Sys,Y> pos[3];
				// rotate sample
				pos[0] = ctheta_s*x + stheta_s*z;
				pos[1] = y;
				pos[2] = -stheta_s*x + ctheta_s*z;

				t_length<Sys,Y> path_diff = q_dir[0]*pos[0] + q_dir[1]*pos[1] + q_dir[2]*pos[2];
				Y phase = omegaM * path_diff / v;

				t_length<Sys,Y> dist = z/*+len_z/2.*/;
				dist /= ctheta_s;

				Y extinction_factor = exp(-mu * dist * 2.); // reflexive: path taken twice

				//if(extinction_factor > 1 / std::exp(1))
				{
					vol += extinction_factor * func_det;
					integral += extinction_factor * func_det * cos(phase);
				}
			}

	return integral / vol;
}


// Bragg scattering with extinction
template<class Sys, class Y>
Y mieze_reduction_sample_cuboid_bragg(const t_length<Sys,Y>& len_x,
	const t_length<Sys,Y>& len_y,
	const t_length<Sys,Y>& len_z,
	const t_length_inverse<Sys, Y>& mu,
	const t_freq<Sys,Y>& fM,
	const t_length<Sys,Y>& lam,
	const t_length<Sys,Y>& d_ana,
	std::size_t ITERS=100)
{
	const t_angle<Sys,Y>
		twotheta = bragg_real_twotheta(d_ana, lam, 1.);

	return mieze_reduction_sample_cuboid_extinction(len_x, len_y, len_z, mu, fM, lam, twotheta, ITERS);
}

//------------------------------------------------------------------------------





//------------------------------------------------------------------------------

template<class Sys, class Y>
t_length<Sys,Y> mieze_condition_L2(const t_freq<Sys,Y>& f1,
	const t_freq<Sys,Y>& f2, const t_length<Sys,Y>& L1)
{
        return L1 / (f2/f1 - 1.);
}

template<class Sys, class Y>
t_length<Sys,Y> mieze_condition_L1(const t_freq<Sys,Y>& f1,
	const t_freq<Sys,Y>& f2, const t_length<Sys,Y>& L2)
{
        return L2 * (f2/f1 - 1.);
}

template<class Sys, class Y>
t_freq<Sys,Y> mieze_condition_f2(const t_freq<Sys,Y>& f1,
	const t_length<Sys,Y>& L1, const t_length<Sys,Y>& L2)
{
        return (L1/L2 + 1.)*f1;
}

template<class Sys, class Y>
t_freq<Sys,Y> mieze_condition_f1(const t_freq<Sys,Y>& f2,
	const t_length<Sys,Y>& L1, const t_length<Sys,Y>& L2)
{
        return f2 / (L1/L2 + 1.);
}



template<class Sys, class Y>
t_length<Sys,Y> mieze_condition_inel_Ls(const t_length<Sys,Y>& Ls0,
	const t_energy<Sys,Y>& dE, const t_length<Sys,Y>& lam)
{
	t_velocity<Sys,Y> v0 = lam2p(lam)/co::m_n;

	t_velocity<Sys,Y> dv = sqrt(v0*v0 + 2.*dE/co::m_n) - v0;
	t_velocity<Sys,Y> v1 = v0 + dv;
	return Ls0 * (v1*v1*v1/(v0*v0*v0));
}


//------------------------------------------------------------------------------

template<class Sys, class Y>
t_freq<Sys,Y> mieze_det_misaligned_df1(const t_length<Sys,Y>& L1,
	const t_length<Sys,Y>& L2, const t_length<Sys,Y>& dL,
	const t_freq<Sys,Y>& f1, const t_freq<Sys,Y>& f2)
{
	t_freq<Sys,Y> df;
	df = (f2-f1) - (L1*f2)/(L1+L2+dL);
	return df;
}

template<class Sys, class Y>
t_freq<Sys,Y> mieze_det_misaligned_df2(const t_length<Sys,Y>& L1,
	const t_length<Sys,Y>& L2, const t_length<Sys,Y>& dL,
	const t_freq<Sys,Y>& f1, const t_freq<Sys,Y>& f2)
{
	t_freq<Sys,Y> df;
	df = (L1 / (L2+dL) + 1.) * f1 - f2;
	return df;
}

//------------------------------------------------------------------------------
}

#endif
