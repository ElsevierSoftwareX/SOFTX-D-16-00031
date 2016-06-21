/**
 * magnetic dispersion relations
 * @author tweber
 * @date 7-jul-15
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_MAGDISP_H__
#define __TLIBS_MAGDISP_H__

#include <initializer_list>
#include <vector>
#include <cmath>
#include <complex>
#include <cassert>

#include "linalg.h"
#include "atoms.h"
#include "nn.h"


namespace tl {
// ----------------------------------------------------------------------------

/**
 * Simple ferromagnetic dispersion
 * @param lstNeighbours list of distances to neighbour atoms and their coupling constants
 * @param vecq q position
 * @param tS spin
 * @return E(q)
 */
template<class t_vec = ublas::vector<double>,
	typename T = typename t_vec::value_type,
	template<class...> class t_cont = std::vector>
T ferromag(const t_cont<t_vec>& vecNeighbours, const t_cont<std::complex<T>>& vecJ,
	const ublas::vector<T>& vecq, T tS)
{
	std::complex<T> J(0., 0.), J0(0., 0.);
	J = structfact<T, std::complex<T>, t_vec, t_cont>
		(vecNeighbours, vecq, vecJ, &J0).real();
	return T(2)*tS*(J0 - J).real();
}

template<class t_vec = ublas::vector<double>,
	typename T = typename t_vec::value_type,
	typename t_cont = std::initializer_list<std::tuple<t_vec, std::complex<T>>>>
T ferromag(const t_cont& lstNeighbours, const ublas::vector<T>& vecq, T tS)
{
	return ferromag(vec_from_pairvec<0,std::vector,t_cont>()(lstNeighbours),
		vec_from_pairvec<1,std::vector,t_cont>()(lstNeighbours),
		vecq, tS);
}
// ----------------------------------------------------------------------------


// Magnetic form factors
// see: ILL Neutron Data Booklet sec. 2.5-1 (p. 60)
// also see: https://www.ill.eu/sites/ccsl/ffacts/

template<class T=double, template<class...> class t_vec=std::initializer_list>
T j0_avg(T Q, const t_vec<T>& A, const t_vec<T>& a)
{
	assert(A.size() == a.size()+1);

	T tJ = T(0);
	for(std::size_t i=0; i<a.size(); ++i)
		tJ += A[i] * std::exp(-a[i] * Q/(T(4)*get_pi<T>())*Q/(T(4)*get_pi<T>()));
	tJ += *A.rbegin();
	return tJ;
}

template<class T=double, template<class...> class t_vec=std::initializer_list>
T j2_avg(T Q, const t_vec<T>& A, const t_vec<T>& a)
{
	return j0_avg<T, t_vec>(Q, A, a) * Q/(T(4)*get_pi<T>()) * Q/(T(4)*get_pi<T>());
}

template<class T=double, template<class...> class t_vec=std::initializer_list>
T mag_formfact(T Q, T L, T S,
	const t_vec<T>& A0, const t_vec<T>& a0,
	const t_vec<T>& A2, const t_vec<T>& a2)
{
	return (L+T(2)*S) * j0_avg<T, t_vec>(Q, A0, a0) * L * j2_avg<T, t_vec>(Q, A2, a2);
}

// see: Squires, p. 139
template<class T=double, template<class...> class t_vec=std::initializer_list>
T mag_formfact(T Q, T L, T S, T J,
	const t_vec<T>& A0, const t_vec<T>& a0,
	const t_vec<T>& A2, const t_vec<T>& a2)
{
	T j0 = j0_avg<T, t_vec>(Q, A0, a0);
	T j2 = j2_avg<T, t_vec>(Q, A2, a2);

	T gL = T(0.5) + (L*(L+T(1)) - S*(S+T(1))) / (T(2)*J* (J+T(1)));
	T gS = T(1) + (S*(S+T(1)) - L*(L+T(1))) / (J * (J+T(1)));

	return (gS*j0 + gL*(j0+j2)) / (gL + gS);
}
// ----------------------------------------------------------------------------

}
#endif
