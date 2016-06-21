/*
 * random numbers
 * @author tweber
 * @date 16-aug-2013
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_RAND_H__
#define __TLIBS_RAND_H__

#include <random>
#include <vector>
#include <initializer_list>
#include <type_traits>
#include <future>

namespace tl {

extern std::mt19937/*_64*/ g_randeng;


extern void init_rand();
extern void init_rand_seed(unsigned int uiSeed);


extern unsigned int simple_rand(unsigned int iMax);


template<typename INT>
INT rand_int(int iMin, int iMax)
{
	std::uniform_int_distribution<INT> dist(iMin, iMax);
	return dist(g_randeng);
}

template<typename REAL>
REAL rand_real(REAL dMin, REAL dMax)
{
	std::uniform_real_distribution<REAL> dist(dMin, dMax);
	return dist(g_randeng);
}

template<typename REAL>
REAL rand01()
{
	return rand_real<REAL>(0., 1.);
}

template<typename REAL>
REAL rand_norm(REAL dMu, REAL dSigma)
{
	std::normal_distribution<REAL> dist(dMu, dSigma);
	return dist(g_randeng);
}

template<typename INT, typename REAL=double>
INT rand_poisson(REAL dMu)
{
	std::poisson_distribution<INT> dist(dMu);
	return dist(g_randeng);
}


template<class t_vec=std::vector<double>, class REAL=typename t_vec::value_type,
	class t_initlst=std::initializer_list<REAL>>
t_vec rand_norm_nd(const t_initlst& vecMu, const t_initlst& vecSigma)
{
	if(vecMu.size() != vecSigma.size())
		return t_vec();

	unsigned int iDim = vecMu.size();
	t_vec vecRet(iDim);

	std::vector<std::future<REAL>> vecFut;
	vecFut.reserve(iDim);

	using iter = typename t_initlst::const_iterator;
	iter iterMu = vecMu.begin();
	iter iterSig = vecSigma.begin();

	for(; iterMu!=vecMu.end(); ++iterMu, ++iterSig)
	{
		vecFut.emplace_back(
			std::async(std::launch::deferred | std::launch::async,
			std::function<REAL(REAL,REAL)>(rand_norm<REAL>), *iterMu, *iterSig));
	}

	for(unsigned int i=0; i<iDim; ++i)
		vecRet[i] = vecFut[i].get();

	return vecRet;
}
}

#endif
