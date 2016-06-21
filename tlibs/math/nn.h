/*
 * next neighbours
 * @author tweber
 * @date may-2016
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_NN_H__
#define __TLIBS_NN_H__

#include <initializer_list>
#include <vector>
#include <tuple>
#include <cmath>
#include <complex>

#include "linalg.h"
#include "../helper/misc.h"


namespace tl {

/**
 * get nearest neighbours
 * @return vector of neighbours
 */
template<class t_vec = ublas::vector<double>,
	template<class...> class t_cont = std::vector,
	class t_real = typename t_vec::value_type>
t_cont<t_cont<std::size_t>>
get_neighbours(const t_cont<t_vec>& vecAtoms, const t_vec& vecCentre,
	t_real epsShell = t_real(1e-3))
{
	t_cont<t_cont<std::size_t>> vecN;

	// generate lengths
	t_cont<t_real> vecLens;
	vecLens.reserve(vecAtoms.size());
	for(const t_vec& vec : vecAtoms)
		vecLens.push_back(ublas::norm_2(vec - vecCentre));

	// sort by lengths
	t_cont<std::size_t> vecIdx = sorted_idx(vecLens,
		[&vecAtoms, &vecLens](t_real len0, t_real len1) -> bool
	{
		return len0 < len1;
	});

	if(vecIdx.size() == 0) return vecN;

	// sort by nearest neighbour, next-nearest neighbour, etc.
	t_real distLast = vecLens[vecIdx[0]];
	t_cont<std::size_t> vecNearest;

	for(typename decltype(vecIdx)::const_iterator iter=vecIdx.begin(); iter!=vecIdx.end(); ++iter)
	{
		std::size_t iIdx = *iter;

		t_real distCur = vecLens[iIdx];
		if(float_equal(distCur, distLast, epsShell))
		{
			vecNearest.push_back(iIdx);
		}
		else
		{
			vecN.push_back(std::move(vecNearest));

			vecNearest.clear();
			vecNearest.push_back(iIdx);
			distLast = distCur;
		}

		if(iter+1 == vecIdx.end() && vecNearest.size())
			vecN.push_back(vecNearest);
	}
	return vecN;
}


// ----------------------------------------------------------------------------
// cubic systems

// atom positions
enum class UCType { SIMPLE, FCC, BCC, };

/**
 * Next neighbours
 * iDist == 0: nearest neighbours
 * iDist == 1: next-nearest neighbours
 */
template<typename T=double>
std::vector<ublas::vector<T>> get_neighbour_atoms(UCType crys, int iDist=0, T a=1.)
{
	std::vector<ublas::vector<T>> vecAtoms;

	if(crys == UCType::SIMPLE)
	{
		if(iDist == 0)
			vecAtoms = {
				make_vec({1., 0., 0.}),
				make_vec({0., 1., 0.}),
				make_vec({0., 0., 1.}),
				make_vec({-1., 0., 0.}),
				make_vec({0., -1., 0.}),
				make_vec({0., 0., -1.}) };
		else if(iDist == 1)
			vecAtoms = {
				make_vec({1., 1., 0.}),
				make_vec({1., -1., 0.}),
				make_vec({-1., 1., 0.}),
				make_vec({-1., -1., 0.}),
				make_vec({1., 0., 1.}),
				make_vec({1., 0., -1.}),
				make_vec({-1., 0., 1.}),
				make_vec({-1., 0., -1.}),
				make_vec({0., 1., 1.}),
				make_vec({0., 1., -1.}),
				make_vec({0., -1., 1.}),
				make_vec({0., -1., -1.}) };
	}
	else if(crys == UCType::FCC)
	{
		if(iDist == 0)
			vecAtoms = {
				make_vec({0.5, 0.5, 0.}),
				make_vec({0.5, -0.5, 0.}),
				make_vec({-0.5, 0.5, 0.}),
				make_vec({-0.5, -0.5, 0.}),
				make_vec({0.5, 0., 0.5}),
				make_vec({0.5, 0., -0.5}),
				make_vec({-0.5, 0., 0.5}),
				make_vec({-0.5, 0., -0.5}),
				make_vec({0., 0.5, 0.5}),
				make_vec({0., 0.5, -0.5}),
				make_vec({0., -0.5, 0.5}),
				make_vec({0., -0.5, -0.5}) };
		else if(iDist == 1)
			vecAtoms = {
				make_vec({1., 0., 0.}),
				make_vec({0., 1., 0.}),
				make_vec({0., 0., 1.}),
				make_vec({-1., 0., 0.}),
				make_vec({0., -1., 0.}),
				make_vec({0., 0., -1.}) };
	}
	else if(crys == UCType::BCC)
	{
		if(iDist == 0)
			vecAtoms = {
				make_vec({0.5, 0.5, 0.5}),
				make_vec({0.5, 0.5, -0.5}),
				make_vec({0.5, -0.5, 0.5}),
				make_vec({0.5, -0.5, -0.5}),
				make_vec({-0.5, 0.5, 0.5}),
				make_vec({-0.5, 0.5, -0.5}),
				make_vec({-0.5, -0.5, 0.5}),
				make_vec({-0.5, -0.5, -0.5}) };
		else if(iDist == 1)
			vecAtoms = {
				make_vec({1., 0., 0.}),
				make_vec({0., 1., 0.}),
				make_vec({0., 0., 1.}),
				make_vec({-1., 0., 0.}),
				make_vec({0., -1., 0.}),
				make_vec({0., 0., -1.}) };
	}

	if(!float_equal<T>(a, T(1)))
		for(ublas::vector<T>& vec : vecAtoms)
			vec *= a;

	return vecAtoms;
}
// ----------------------------------------------------------------------------

}
#endif
