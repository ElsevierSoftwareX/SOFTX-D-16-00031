/**
 * calculating term symbols
 * @author Tobias Weber
 * @date 2016
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_TERM_H__
#define __TLIBS_TERM_H__

#include <tuple>
#include <vector>
#include <numeric>
#include <cstdint>
#include <cmath>
#include <sstream>
#include <unordered_map>

#include "../string/string.h"
#include "../helper/exception.h"

namespace tl
{
	/**
	 * Hund's rules
	 *
	 * @return [S, L, J]
	 */
	template<class t_real = double>
	std::tuple<t_real, t_real, t_real>
	hund(std::uint16_t l, std::uint16_t iNumEs)
	{
		std::uint16_t iNumOrbitals = 2*l+1;
		if(iNumEs > iNumOrbitals*2)
			throw Err("Too many electrons.");

		std::vector<std::uint8_t> vecOrbitals;	// orbitals
		std::vector<std::int16_t> vec_ml;	// mag. q.number
		vecOrbitals.resize(iNumOrbitals);
		vec_ml.resize(iNumOrbitals);
		std::iota(vec_ml.rbegin(), vec_ml.rend(), -l);

		for(std::uint16_t iE=0; iE<iNumEs; ++iE)
			++vecOrbitals[iE%iNumOrbitals];

		t_real S=0, L=0, J=0;
		for(std::size_t iOrbital=0; iOrbital<vecOrbitals.size(); ++iOrbital)
		{
			std::uint8_t iEs = vecOrbitals[iOrbital];
			if(iEs==1)	// unpaired electron
				S += t_real(0.5);

			std::int16_t ml = vec_ml[iOrbital];
			L += t_real(std::int16_t(iEs)*ml);
		}

		if(iNumEs <= iNumOrbitals)
			J = std::abs(L-S);
		else
			J = L+S;

		return std::make_tuple(S,L,J);
	}


	template<class t_real=double, class t_str=std::string>
	t_str get_termsymbol(t_real S, t_real L, t_real J)
	{
		static const std::vector<t_str> vecL =
			{"S","P","D","F","G","H","I","K","L","M","N","O"};

		t_str strS = var_to_str<t_real, t_str>(t_real(2)*S+1);
		t_str strL = vecL[std::size_t(L)];
		t_str strJ = var_to_str<t_real, t_str>(J);

		return strS + strL + strJ;
	}


	/**
	* transforms e.g. 1s2 -> [1,0,2]
	* @return [n, l, #electrons]
	*/
	template<class t_str=std::string>
	std::tuple<uint16_t, uint16_t, uint16_t>
	get_orbital(const t_str& strOrbital)
	{
		using t_ch = typename t_str::value_type;
		static const std::unordered_map<t_ch, uint16_t> mapSubOrbitals =
		{
			{'s',0}, {'p',1}, {'d',2}, {'f',3},
			{'g',4}, {'h',5}, {'i',6}, {'k',7},
			{'l',8}, {'m',9}, {'n',10}, {'o',11},
		};

		std::istringstream istr(strOrbital);
		uint16_t n = 0, l = 0, iNumE = 0;
		t_ch cSub = 's';

		istr >> n >> cSub >> iNumE;
		auto iter = mapSubOrbitals.find(cSub);
		if(iter == mapSubOrbitals.end())
			throw Err("Invalid orbital.");
		l = iter->second;

		return std::make_tuple(n,l,iNumE);
	}


	/*
	 * gets term symbol from orbitals
	 * @return [S, L, J]
	 */
	template<class t_real=double, class t_str=std::string>
	std::tuple<t_real, t_real, t_real>
	hund(const t_str& strOrbitals)
	{
		std::tuple<t_real, t_real, t_real> tupTerm(0,0,0);
		std::vector<t_str> vecOrbitals;
		tl::get_tokens<t_str,t_str>(strOrbitals, " ,;", vecOrbitals);

		// all orbitals
		for(const t_str& strOrbital : vecOrbitals)
		{
			std::tuple<uint16_t, uint16_t, uint16_t> tup_nle =
				get_orbital(strOrbital);
			std::tuple<t_real, t_real, t_real> tup =
				hund(std::get<1>(tup_nle), std::get<2>(tup_nle));

			std::get<0>(tupTerm) += std::get<0>(tup);
			std::get<1>(tupTerm) += std::get<1>(tup);
			std::get<2>(tupTerm) += std::get<2>(tup);
		}

		return tupTerm;
	}
}

#endif
