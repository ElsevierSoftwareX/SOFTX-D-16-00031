/*
 * neutron formulas
 * @author Tobias Weber
 * @date 2012-2016
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_NEUTRONS__
#define __TLIBS_NEUTRONS__

#include "math.h"
#include "linalg.h"
#include "units.h"
#include "../helper/exception.h"

#include <boost/units/pow.hpp>
#include <cmath>


namespace tl {

// --------------------------------------------------------------------------------

template<typename T=double> T get_KSQ2E()
{
	const auto _A = get_one_angstrom<T>();
	const auto _meV = get_one_meV<T>();
	const auto _mn = get_m_n<T>();
	const auto _hbar = get_hbar<T>();

	// factor order important for small value types like "float"!
	return T(0.5) * _hbar/_A/_mn * _hbar/_A/_meV;
}

template<typename T=double> T get_E2KSQ()
{
	return T(1)/get_KSQ2E<T>();
}

#if __cplusplus >= 201402L
	template<class T=double> T t_KSQ2E = get_KSQ2E<T>();
	template<class T=double> T t_E2KSQ = T(1)/t_KSQ2E<T>;
#endif

static const double KSQ2E = (co::hbar*co::hbar / (2.*co::m_n)) / one_meV / (angstrom*angstrom);
static const double E2KSQ = 1./KSQ2E;

// --------------------------------------------------------------------------------


// --------------------------------------------------------------------------------
// de Broglie stuff
// lam = h/p
template<class Sys, class Y>
t_momentum<Sys,Y> lam2p(const t_length<Sys,Y>& lam)
{
	return get_h<Y>() / lam;
}

template<class Sys, class Y>
t_length<Sys,Y> p2lam(const t_momentum<Sys,Y>& p)
{
	return get_h<Y>() / p;
}


// lam = 2pi/k
template<class Sys, class Y>
t_length<Sys,Y> k2lam(const t_wavenumber<Sys,Y>& k)
{
	return Y(2.)*get_pi<Y>() / k;
}

template<class Sys, class Y>
t_wavenumber<Sys,Y> lam2k(const t_length<Sys,Y>& lam)
{
	return Y(2.)*get_pi<Y>() / lam;
}

template<class Sys, class Y>
t_momentum<Sys,Y> k2p(const t_wavenumber<Sys,Y>& k)
{
	return get_hbar<Y>()*k;
}

template<class Sys, class Y>
t_wavenumber<Sys,Y> p2k(const t_momentum<Sys,Y>& p)
{
	return p/get_hbar<Y>();
}

template<class Sys, class Y>
t_velocity<Sys,Y> k2v(const t_wavenumber<Sys,Y>& k)
{
	return k2p(k) / get_m_n<Y>();
}

template<class Sys, class Y>
t_wavenumber<Sys,Y> v2k(const t_velocity<Sys,Y>& v)
{
	return get_m_n<Y>()*v/get_hbar<Y>();
}
// --------------------------------------------------------------------------------


// --------------------------------------------------------------------------------
// E = hbar*omega
template<class Sys, class Y>
t_energy<Sys,Y> omega2E(const t_freq<Sys,Y>& omega)
{
	return get_hbar<Y>() * omega;
}

template<class Sys, class Y>
t_freq<Sys,Y> E2omega(const t_energy<Sys,Y>& en)
{
	return en / get_hbar<Y>();
}

template<class Sys, class Y>
t_energy<Sys,Y> k2E_direct(const t_wavenumber<Sys,Y>& k)
{
	t_momentum<Sys,Y> p = get_hbar<Y>()*k;
	t_energy<Sys,Y> E = p*p / (Y(2.)*get_m_n<Y>());
	return E;
}

template<class Sys, class Y>
t_wavenumber<Sys,Y> E2k_direct(const t_energy<Sys,Y>& _E, bool &bImag)
{
	bImag = (_E < Y(0.)*get_one_meV<Y>());
	t_energy<Sys,Y> E = bImag ? -_E : _E;

	auto pp = Y(2.) * get_m_n<Y>() * E;
	//t_momentum<Sys,Y> p = units::sqrt<typename decltype(pp)::unit_type, Y>(pp);
	t_momentum<Sys,Y> p = my_units_sqrt<t_momentum<Sys,Y>>(pp);
	t_wavenumber<Sys,Y> k = p / get_hbar<Y>();
	return k;
}
// --------------------------------------------------------------------------------


// ----------------------------------------------------------------------------
// indirect calculations using conversion factors for numerical stability
template<class Sys, class Y>
t_energy<Sys,Y> k2E(const t_wavenumber<Sys,Y>& k)
{
	Y dk = k*get_one_angstrom<Y>();
	Y dE = get_KSQ2E<Y>() * dk*dk;
	return dE * get_one_meV<Y>();
}

template<class Sys, class Y>
t_wavenumber<Sys,Y> E2k(const t_energy<Sys,Y>& _E, bool &bImag)
{
	bImag = (_E < Y(0.)*get_one_meV<Y>());
	t_energy<Sys,Y> E = bImag ? -_E : _E;
	const Y dE = E / get_one_meV<Y>();
	const Y dk = std::sqrt(get_E2KSQ<Y>() * dE);
	return dk / get_one_angstrom<Y>();
}
// ----------------------------------------------------------------------------



// --------------------------------------------------------------------------------
// Bragg
// real: n * lam = 2d * sin(twotheta/2)
template<class Sys, class Y>
t_length<Sys,Y> bragg_real_lam(const t_length<Sys,Y>& d,
	const t_angle<Sys,Y>& twotheta, Y n)
{
	return Y(2.)*d/n * units::sin(twotheta/Y(2.));
}

template<class Sys, class Y>
t_length<Sys,Y> bragg_real_d(const t_length<Sys,Y>& lam,
	const t_angle<Sys,Y>& twotheta, Y n)
{
	return n * lam / (Y(2.)* units::sin(twotheta/Y(2.)));
}

template<class Sys, class Y>
t_angle<Sys,Y> bragg_real_twotheta(const t_length<Sys,Y>& d,
	const t_length<Sys,Y>& lam, Y n)
{
	return units::asin(n*lam/(Y(2.)*d)) * Y(2.);
}

// reciprocal: Q * lam = 4pi * sin(twotheta/2)
template<class Sys, class Y>
t_angle<Sys,Y> bragg_recip_twotheta(const t_wavenumber<Sys,Y>& Q,
	const t_length<Sys,Y>& lam, Y n)
{
	return units::asin(Q*n*lam/(Y(4)*get_pi<Y>())) * Y(2);
}

template<class Sys, class Y>
t_wavenumber<Sys,Y> bragg_recip_Q(const t_length<Sys,Y>& lam,
	const t_angle<Sys,Y>& twotheta, Y n)
{
	return Y(4)*get_pi<Y>() / (n*lam) * units::sin(twotheta/Y(2));
}

template<class Sys, class Y>
t_length<Sys,Y> bragg_recip_lam(const t_wavenumber<Sys,Y>& Q,
	const t_angle<Sys,Y>& twotheta, Y n)
{
	return Y(4)*get_pi<Y>() / Q * units::sin(twotheta/Y(2)) / n;
}
// --------------------------------------------------------------------------------


// --------------------------------------------------------------------------------

// see e.g. ILL blue book sec. 2.6-2
template<class Sys, class Y>
t_wavenumber<Sys,Y> kinematic_plane(bool bFixedKi,
	const t_energy<Sys,Y>& EiEf, const t_energy<Sys,Y>& DeltaE,
	const t_angle<Sys,Y>& twotheta)
{
	const t_energy<Sys,Y> dE = DeltaE;
	if(bFixedKi)
		dE = -dE;

	t_wavenumber<Sys,Y> Q =
		units::sqrt(Y(2.)*get_m_n<Y>() / get_hbar<Y>()) *
		(Y(2.)*EiEf + dE - Y(2.)*units::cos(twotheta) *
		units::sqrt(EiEf*(EiEf + dE)));

	return Q;
}

template<class Sys, class Y>
t_energy<Sys,Y> kinematic_plane(bool bFixedKi, bool bBranch,
	const t_energy<Sys,Y>& EiEf, const t_wavenumber<Sys,Y>& Q,
	const t_angle<Sys,Y>& twotheta)
{
	auto c = Y(2.)*get_m_n<Y>() / (get_hbar<Y>()*get_hbar<Y>());
	Y ctt = units::cos(twotheta);
	Y c2tt = units::cos(Y(2.)*twotheta);

	Y dSign = Y(-1.);
	if(bBranch)
		dSign = Y(1.);

	Y dSignFixedKf = Y(1.);
	if(bFixedKi)
		dSignFixedKf = Y(-1.);

	using t_sqrt_rt = decltype(c*c*EiEf*ctt);
	using t_rt = decltype(t_sqrt_rt()*t_sqrt_rt());
	t_rt rt = c*c*c*c * (-EiEf*EiEf)*ctt*ctt
		+ c*c*c*c*EiEf*EiEf*ctt*ctt*c2tt
		+ Y(2.)*c*c*c*EiEf*Q*Q*ctt*ctt;

	t_energy<Sys,Y> E =
		Y(1.)/(c*c)*(dSignFixedKf*Y(2.)*c*c*EiEf*ctt*ctt
		- dSignFixedKf*Y(2.)*c*c*EiEf
		+ dSign*std::sqrt(Y(2.)) * my_units_sqrt<t_sqrt_rt>(rt)
		+ dSignFixedKf*c*Q*Q);

	return E;
}
// --------------------------------------------------------------------------------



// --------------------------------------------------------------------------------
// Debye-Waller factor

template<class Sys, class Y>
Y debye_waller_high_T(const t_temperature<Sys,Y>& T_D,
	const t_temperature<Sys,Y>& T, const t_mass<Sys,Y>& M,
	const t_wavenumber<Sys,Y>& Q, t_length_square<Sys,Y>* pZeta_sq=0)
{
	t_length_square<Sys,Y> zeta_sq;
	zeta_sq = Y(9.)*get_hbar<Y>()/get_kB<Y>() / (T_D * M) * T/T_D * get_hbar<Y>();
	Y dwf = units::exp(Y(-1./3.) * Q*Q * zeta_sq);

	if(pZeta_sq) *pZeta_sq = zeta_sq;
	return dwf;
}


template<class Sys, class Y>
Y debye_waller_low_T(const t_temperature<Sys,Y>& T_D,
	const t_temperature<Sys,Y>& T, const t_mass<Sys,Y>& M,
	const t_wavenumber<Sys,Y>& Q, t_length_square<Sys,Y>* pZeta_sq=0)
{
	t_length_square<Sys,Y> zeta_sq;
	zeta_sq = Y(9.)*get_hbar<Y>()/get_kB<Y>() / (Y(4.)*T_D*M) * get_hbar<Y>() *
		(Y(1.) + Y(2./3.) * get_pi<Y>()*get_pi<Y>() * (T/T_D)*(T/T_D));
	Y dwf = units::exp(Y(-1./3.) * Q*Q * zeta_sq);

	if(pZeta_sq) *pZeta_sq = zeta_sq;
	return dwf;
}

// --------------------------------------------------------------------------------


// --------------------------------------------------------------------------------
// scattering triangle / TAS stuff

// Q_vec = ki_vec - kf_vec
// kf_vec = ki_vec - Q_vec
// kf^2 = ki^2 + Q^2 - 2ki Q cos th
// cos th = (-kf^2 + ki^2 + Q^2) / (2kiQ)
template<class Sys, class Y>
t_angle<Sys,Y> get_angle_ki_Q(const t_wavenumber<Sys,Y>& ki,
	const t_wavenumber<Sys,Y>& kf,
	const t_wavenumber<Sys,Y>& Q,
	bool bPosSense=1, bool bAngleOutsideTriag=0)
{
	t_angle<Sys,Y> angle;

	if(Q*get_one_angstrom<Y>() == Y(0.))
		angle = get_pi<Y>()/Y(2) * units::si::radians;
	else
	{
		auto c = (ki*ki - kf*kf + Q*Q)/(Y(2.)*ki*Q);
		if(units::abs(c) > Y(1.))
			throw Err("Scattering triangle not closed.");

		angle = units::acos(c);
	}

	if(bAngleOutsideTriag) angle = get_pi<Y>()*units::si::radians - angle;
	if(!bPosSense) angle = -angle;

	return angle;
}

// Q_vec = ki_vec - kf_vec
// ki_vec = Q_vec + kf_vec
// ki^2 = Q^2 + kf^2 + 2Q kf cos th
// cos th = (ki^2 - Q^2 - kf^2) / (2Q kf)
template<class Sys, class Y>
t_angle<Sys,Y> get_angle_kf_Q(const t_wavenumber<Sys,Y>& ki,
	const t_wavenumber<Sys,Y>& kf,
	const t_wavenumber<Sys,Y>& Q,
	bool bPosSense=1, bool bAngleOutsideTriag=1)
{
	t_angle<Sys,Y> angle;

	if(Q*get_one_angstrom<Y>() == Y(0.))
		angle = get_pi<Y>()/Y(2) * units::si::radians;
	else
	{
		auto c = (-kf*kf + ki*ki - Q*Q)/(Y(2.)*kf*Q);
		if(units::abs(c) > Y(1.))
			throw Err("Scattering triangle not closed.");

		angle = units::acos(c);
	}

	if(!bAngleOutsideTriag) angle = get_pi<Y>()*units::si::radians - angle;
	if(!bPosSense) angle = -angle;

	return angle;
}


template<class Sys, class Y>
t_angle<Sys,Y> get_mono_twotheta(const t_wavenumber<Sys,Y>& k,
	const t_length<Sys,Y>& d, bool bPosSense=1)
{
	const Y dOrder = Y(1.);
	t_angle<Sys,Y> tt = bragg_real_twotheta(d, k2lam(k), dOrder);
	if(!bPosSense)
		tt = -tt;
	return tt;
}

template<class Sys, class Y>
t_wavenumber<Sys,Y> get_mono_k(const t_angle<Sys,Y>& _theta,
	const t_length<Sys,Y>& d, bool bPosSense=1)
{
	t_angle<Sys,Y> theta = _theta;
	if(!bPosSense)
		theta = -theta;

	const Y dOrder = Y(1.);
	return lam2k(bragg_real_lam(d, Y(2.)*theta, dOrder));
}


// Q_vec = ki_vec - kf_vec
// Q^2 = ki^2 + kf^2 - 2ki kf cos 2th
// cos 2th = (-Q^2 + ki^2 + kf^2) / (2ki kf)
template<class Sys, class Y>
t_angle<Sys,Y> get_sample_twotheta(const t_wavenumber<Sys,Y>& ki,
	const t_wavenumber<Sys,Y>& kf, const t_wavenumber<Sys,Y>& Q,
	bool bPosSense=1)
{
	t_dimensionless<Sys,Y> ttCos = (ki*ki + kf*kf - Q*Q)/(Y(2.)*ki*kf);
	if(units::abs(ttCos) > Y(1.))
		throw Err("Scattering triangle not closed.");

	t_angle<Sys,Y> tt;
	tt = units::acos(ttCos);

	if(!bPosSense) tt = -tt;
	return tt;
}


// again cos theorem:
// Q_vec = ki_vec - kf_vec
// Q^2 = ki^2 + kf^2 - 2ki kf cos 2th
// Q = sqrt(ki^2 + kf^2 - 2ki kf cos 2th)
template<class Sys, class Y>
const t_wavenumber<Sys,Y>
get_sample_Q(const t_wavenumber<Sys,Y>& ki,
	const t_wavenumber<Sys,Y>& kf, const t_angle<Sys,Y>& tt)
{
	t_dimensionless<Sys,Y> ctt = units::cos(tt);
	decltype(ki*ki) Qsq = ki*ki + kf*kf - Y(2.)*ki*kf*ctt;
	if(Y(Qsq*get_one_angstrom()*get_one_angstrom()) < Y(0.))
	{
		// TODO

		Qsq = -Qsq;
	}

	//t_wavenumber<Sys,Y> Q = units::sqrt(Qsq);
	t_wavenumber<Sys,Y> Q = my_units_sqrt<t_wavenumber<Sys,Y>>(Qsq);
	return Q;
}



template<class Sys, class Y>
t_energy<Sys,Y> get_energy_transfer(const t_wavenumber<Sys,Y>& ki,
	const t_wavenumber<Sys,Y>& kf)
{
	return k2E<Sys,Y>(ki) - k2E<Sys,Y>(kf);
}


// (hbar*ki)^2 / (2*mn)  -  (hbar*kf)^2 / (2mn)  =  E
// 1) ki^2  =  +E * 2*mn / hbar^2  +  kf^2
// 2) kf^2  =  -E * 2*mn / hbar^2  +  ki^2
template<class Sys, class Y>
t_wavenumber<Sys,Y> get_other_k(const t_energy<Sys,Y>& E,
	const t_wavenumber<Sys,Y>& kfix, bool bFixedKi)
{
	auto kE_sq = E*Y(2.)*(get_m_n<Y>()/get_hbar<Y>())/get_hbar<Y>();
	if(bFixedKi) kE_sq = -kE_sq;

	auto k_sq = kE_sq + kfix*kfix;
	if(k_sq*get_one_angstrom()*get_one_angstrom() < Y(0.))
		throw Err("Scattering triangle not closed.");

	//return units::sqrt(k_sq);
	return my_units_sqrt<t_wavenumber<Sys,Y>>(k_sq);
}

// --------------------------------------------------------------------------------



// --------------------------------------------------------------------------------
// kf^3 factor, see e.g. Shirane p. 125

template<class Sys, class Y>
Y ana_effic_factor(const t_wavenumber<Sys, Y>& kf, const t_length<Sys, Y>& d)
{
	t_angle<Sys, Y> theta = Y(0.5)*units::abs(get_mono_twotheta<Sys, Y>(kf, d, true));
	return kf*kf*kf / units::tan(theta) *
		get_one_angstrom()*get_one_angstrom()*get_one_angstrom();
}

// --------------------------------------------------------------------------------



// --------------------------------------------------------------------------------
// spurions

// Bragg tail -> see Shirane p. 152
template<class Sys, class Y>
t_energy<Sys,Y> get_bragg_tail(t_wavenumber<Sys,Y> k,
	t_wavenumber<Sys,Y> q, bool bConstEi=0)
{
	Y t = q / (Y(2.)*k);
	if(!bConstEi)
		t = -t;
	t_energy<Sys,Y> E = get_hbar<Y>()/get_m_n<Y>() * k*q*(Y(1.)+t) * get_hbar<Y>();
	return E;
}


// higher-order inelastic spurions -> Shirane pp. 146-148
template<class Sys, class Y>
t_energy<Sys,Y> get_inelastic_spurion(bool bConstEi, t_energy<Sys,Y> E,
	unsigned int iOrderMono, unsigned int iOrderAna)
{
	const Y dOrderMonoSq = Y(iOrderMono)*Y(iOrderMono);
	const Y dOrderAnaSq = Y(iOrderAna)*Y(iOrderAna);

	t_energy<Sys,Y> E_sp;

	// formulas from Shirane, p. 147
	if(bConstEi)
		E_sp = (Y(1.) - dOrderMonoSq/dOrderAnaSq) * E;
	else
		E_sp = (dOrderAnaSq/dOrderMonoSq - Y(1.)) * E;

	return E_sp;
}

template<class Y=double>
struct InelasticSpurion
{
	Y dE_meV = Y(0.);
	unsigned int iOrderMono = 1;
	unsigned int iOrderAna = 1;
};

template<class Sys, class Y>
std::vector<InelasticSpurion<Y>> check_inelastic_spurions(bool bConstEi,
	t_energy<Sys,Y> Ei, t_energy<Sys,Y> Ef,
	t_energy<Sys,Y> E, unsigned int iMaxOrder=5)
{
	const Y dESensitivity = Y(0.25);	// meV

	std::vector<InelasticSpurion<Y>> vecSpuris;

	for(unsigned int iOrder=1; iOrder<=iMaxOrder; ++iOrder)
	{
		InelasticSpurion<Y> spuri;
		t_energy<Sys,Y> EiEf;

		if(bConstEi)
		{
			spuri.iOrderAna = iOrder;
			EiEf = Ei;
		}
		else
		{
			spuri.iOrderMono = iOrder;
			EiEf = Ef;
		}

		spuri.dE_meV = get_inelastic_spurion(bConstEi, EiEf,
			spuri.iOrderMono, spuri.iOrderAna) / get_one_meV<Y>();

		//std::cout << spuri.dE_meV << " *** " << Y(E/get_one_meV<Y>()) << std::endl;
		if(spuri.dE_meV!=Y(0.) && float_equal<Y>(spuri.dE_meV, Y(E/get_one_meV<Y>()), dESensitivity))
			vecSpuris.push_back(spuri);
	}

	return vecSpuris;
}

struct ElasticSpurion
{
	bool bAType = 0;
	bool bMType = 0;

	bool bAKfSmallerKi = 0;
	bool bMKfSmallerKi = 0;
};

// accidental elastic (currat-axe) spurions -> Shirane pp. 150-155 (esp. fig. 6.2)
template<typename T=double>
ElasticSpurion check_elastic_spurion(const ublas::vector<T>& ki,
	const ublas::vector<T>& kf, const ublas::vector<T>& q)
{
	const T dKi = ublas::norm_2(ki);
	const T dKf = ublas::norm_2(kf);
	const T dq = ublas::norm_2(q);

	const T dAngleSensitivity = T(2.);
	const T dQSensitivity = std::max(dKi, dKf) / T(50.);


	ElasticSpurion result;

	ublas::vector<T> ki_norm = ki;	ki_norm /= dKi;
	ublas::vector<T> kf_norm = kf;	kf_norm /= dKf;

	// Q, q and G point in the opposite direction in Shirane!
	// Shirane: Q = kf - ki, E = Ei - Ef
	// here: Q = ki - kf, E = Ei - Ef
	ublas::vector<T> q_norm = -q;	q_norm /= dq;

	T dAngleKfq = std::acos(ublas::inner_prod(kf_norm, q_norm));
	T dAngleKiq = std::acos(ublas::inner_prod(ki_norm, q_norm));

	//std::cout << "angle ki q: " << dAngleKiq/M_PI*180. << std::endl;
	//std::cout << "angle kf q: " << dAngleKfq/M_PI*180. << std::endl;

	bool bKiqParallel = 0, bkiqAntiParallel = 0;
	bool bKfqParallel = 0, bKfqAntiParallel = 0;

	if(float_equal<T>(dAngleKiq, 0., tl::d2r(dAngleSensitivity)))
		bKiqParallel = 1;
	else if(float_equal<T>(dAngleKiq, get_pi<T>(), tl::d2r(dAngleSensitivity)))
		bkiqAntiParallel = 1;
	if(float_equal<T>(dAngleKfq, 0., tl::d2r(dAngleSensitivity)))
		bKfqParallel = 1;
	else if(float_equal<T>(dAngleKfq, get_pi<T>(), tl::d2r(dAngleSensitivity)))
		bKfqAntiParallel = 1;

	// type A: q || kf, kf > ki
	if(bKfqParallel)
	{
		T dApparentKf = dKf - dq;

		if(float_equal<T>(dApparentKf, dKi, dQSensitivity))
		{
			result.bAType = 1;
			result.bAKfSmallerKi = 0;
		}
	}
	// type A: q || kf, kf < ki
	else if(bKfqAntiParallel)
	{
		T dApparentKf = dKf + dq;

		if(float_equal<T>(dApparentKf, dKi, dQSensitivity))
		{
			result.bAType = 1;
			result.bAKfSmallerKi = 1;
		}
	}

	// type M: q || ki, kf > ki
	if(bKiqParallel)
	{
		T dApparentKi = dKi + dq;

		if(float_equal<T>(dApparentKi, dKf, dQSensitivity))
		{
			result.bMType = 1;
			result.bMKfSmallerKi = 0;
		}
	}
	// type M: q || ki, kf < ki
	else if(bkiqAntiParallel)
	{
		T dApparentKi = dKi - dq;

		if(float_equal<T>(dApparentKi, dKf, dQSensitivity))
		{
			result.bMType = 1;
			result.bMKfSmallerKi = 1;
		}
	}

	return result;
}


// --------------------------------------------------------------------------------

template<class t_real=double>
t_real bose(t_real E, t_real T)
{
	t_real kB = get_kB<t_real>() * units::si::kelvin/tl::get_one_meV();

	if(E >= 0.)
		return 1./(std::exp(std::abs(E)/(kB*T)) - 1.) + 1.;
	else
		return 1./(std::exp(std::abs(E)/(kB*T)) - 1.);
}

template<class Sys, class Y>
Y bose(const t_energy<Sys,Y>& E, const t_temperature<Sys,Y>& T)
{
	return bose<Y>(Y(E/tl::get_one_meV()), Y(T/tl::kelvin));
}

// see: B. Fak, B. Dorner, Physica B 234-236 (1997) pp. 1107-1108
template<class t_real=double>
t_real DHO_model(t_real E, t_real T, t_real E0, t_real hwhm, t_real amp, t_real offs)
{
	//if(E0*E0 - hwhm*hwhm < 0.) return 0.;
	return std::abs(bose<t_real>(E, T)*amp/(E0*get_pi<t_real>()) *
		(hwhm/((E-E0)*(E-E0) + hwhm*hwhm) - hwhm/((E+E0)*(E+E0) + hwhm*hwhm)))
		+ offs;
}


// --------------------------------------------------------------------------------


// get macroscopic from microscopic cross-section
template<class Sys, class Y=double>
t_length_inverse<Sys, Y> macro_xsect(const t_area<Sys, Y>& xsect,
	unsigned int iNumAtoms, const t_volume<Sys, Y>& volUC)
{
	return xsect * Y(iNumAtoms) / volUC;
}



// --------------------------------------------------------------------------------

// thin lens equation: 1/f = 1/lenB + 1/lenA
template<class Sys, class Y=double>
t_length<Sys, Y> focal_len(const t_length<Sys, Y>& lenBefore, const t_length<Sys, Y>& lenAfter)
{
	const t_length_inverse<Sys, Y> f_inv = Y(1)/lenBefore + Y(1)/lenAfter;
	return Y(1) / f_inv;
}

// optimal mono/ana curvature, see e.g. Monochromator_curved.comp in McStas or Shirane p. 66
template<class Sys, class Y=double>
t_length<Sys, Y> foc_curv(const t_length<Sys, Y>& lenBefore, const t_length<Sys, Y>& lenAfter,
	const t_angle<Sys, Y>& tt, bool bVert)
{
	const t_length<Sys, Y> f = focal_len<Sys, Y>(lenBefore, lenAfter);
	const Y s = Y(units::abs(units::sin(Y(0.5)*tt)));

	const t_length<Sys, Y> curv = bVert ? Y(2)*f*s : Y(2)*f/s;
	return curv;
}


}

#endif
