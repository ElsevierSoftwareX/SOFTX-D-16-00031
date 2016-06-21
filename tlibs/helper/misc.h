/*
 * misc helper
 * @author tweber
 * @date 07-mar-2013
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_MISC_HELPER__
#define __TLIBS_MISC_HELPER__

#include <vector>
#include <list>
#include <map>
#include <initializer_list>
#include <algorithm>
#include <iterator>
#include <tuple>
#include <string>
#include "../math/math.h"
#include "array.h"


namespace tl {


template<typename T=double>
T max3(T t1, T t2, T t3)
{
	T tmax = t1;
	tmax = std::max(tmax, t2);
	tmax = std::max(tmax, t3);
	return tmax;
}

template<typename T=double>
T min3(T t1, T t2, T t3)
{
	T tmin = t1;
	tmin = std::min(tmin, t2);
	tmin = std::min(tmin, t3);
	return tmin;
}

template<typename T=double>
T min4(T t1, T t2, T t3, T t4)
{
	T tmin = t1;
	tmin = std::min(tmin, t2);
	tmin = std::min(tmin, t3);
	tmin = std::min(tmin, t4);
	return tmin;
}

template<typename T> T safe_log10(T t, T tInvalid=T(-10))
{
	if(t > T(0))
		return std::log10(t);

	return tInvalid;
}

// pixel -> val
template<typename T=double>
T tic_trafo(unsigned int iDim, T dMin, T dMax, bool bLog, T dPix)
{
	if(bLog)
	{
		dMin = safe_log10<T>(dMin);
		dMax = safe_log10<T>(dMax);
	}

	T dval = dMin + dPix/T(iDim) * (dMax-dMin);
	if(bLog)
		dval = pow(10., dval);

	return dval;
}
// val -> pixel
template<typename T=double>
T tic_trafo_inv(unsigned int iDim, T dMin, T dMax, bool bLog, T dVal)
{
	if(bLog)
	{
		dMin = safe_log10<T>(dMin);
		dMax = safe_log10<T>(dMax);

		dVal = safe_log10<T>(dVal);
	}

	T dpix = (dVal-dMin)/(dMax-dMin) * double(iDim);
	return dpix;
}

template<typename T1, typename T2>
void convert(T1* pDst, const T2* pSrc, std::size_t iSize)
{
	for(std::size_t i=0; i<iSize; ++i)
		pDst[i] = T1(pSrc[i]);
}

template<class vec_type>
typename vec_type::value_type sum_vec(const vec_type& vec)
{
	typename vec_type::value_type val = 0;
	for(const typename vec_type::value_type& v : vec)
		val += v;
	return val;
}

template<typename T>
void apply_fkt(const T* pIn, T* pOut, T(*fkt)(T), std::size_t iSize)
{
	for(std::size_t i=0; i<iSize; ++i)
		pOut[i] = (*fkt)(pIn[i]);
}

template<typename T=double>
unsigned int lerprgb(unsigned char r1, unsigned char g1, unsigned char b1,
	unsigned char r2, unsigned char g2, unsigned char b2,
	T dval)
{
	unsigned char r = lerp(r1, r2, dval);
	unsigned char g = lerp(g1, g2, dval);
	unsigned char b = lerp(b1, b2, dval);

	return (0xff<<24) | (r<<16) | (g<<8) | (b);
}

template<typename T=double>
unsigned int lerprgb(unsigned int col1, unsigned int col2, T dval)
{
	unsigned char r1 = (unsigned char)((col1&0x00ff0000) >> 16);
	unsigned char r2 = (unsigned char)((col2&0x00ff0000) >> 16);

	unsigned char g1 = (unsigned char)((col1&0x0000ff00) >> 8);
	unsigned char g2 = (unsigned char)((col2&0x0000ff00) >> 8);

	unsigned char b1 = (unsigned char)(col1&0x000000ff);
	unsigned char b2 = (unsigned char)(col2&0x000000ff);

	unsigned char r = lerp(r1, r2, dval);
	unsigned char g = lerp(g1, g2, dval);
	unsigned char b = lerp(b1, b2, dval);

	return (0xff<<24) | (r<<16) | (g<<8) | (b);
}

// -----------------------------------------------------------------------------

template<class T>
struct sort_obj
{
	std::vector<T> vec;
};

template<class T>
bool comp_fkt(sort_obj<T> t0, sort_obj<T> t1)
{ return t0.vec[0] < t1.vec[0]; }


// simultaneously sort two arrays
template<class Iter=double*>
void sort_2(Iter begin1, Iter end1, Iter begin2)
{
	typedef typename std::iterator_traits<Iter>::value_type T;

	const std::size_t N = end1-begin1;
	sort_obj<T> *pObj = new sort_obj<T>[N];
	for(std::size_t i=0; i<N; ++i)
	{
		pObj[i].vec.push_back(*(begin1+i));
		pObj[i].vec.push_back(*(begin2+i));
	}

	std::sort(pObj, pObj+N, comp_fkt<T>);
	for(std::size_t i=0; i<N; ++i)
	{
		*(begin1+i) = pObj[i].vec[0];
		*(begin2+i) = pObj[i].vec[1];
	}

	delete[] pObj;
}

// simultaneously sort three arrays
template<class Iter=double*>
void sort_3(Iter begin1, Iter end1, Iter begin2, Iter begin3)
{
	typedef typename std::iterator_traits<Iter>::value_type T;

	const std::size_t N = end1-begin1;
	sort_obj<T> *pObj = new sort_obj<T>[N];
	for(std::size_t i=0; i<N; ++i)
	{
		pObj[i].vec.push_back(*(begin1+i));
		pObj[i].vec.push_back(*(begin2+i));
		pObj[i].vec.push_back(*(begin3+i));
	}

	std::sort(pObj, pObj+N, comp_fkt<T>);
	for(std::size_t i=0; i<N; ++i)
	{
		*(begin1+i) = pObj[i].vec[0];
		*(begin2+i) = pObj[i].vec[1];
		*(begin3+i) = pObj[i].vec[2];
	}

	delete[] pObj;
}


/**
 * get the sorted indices of a container
 */
template<template<class...> class t_cont = std::vector,
	class t_val, class t_func>
t_cont<std::size_t> sorted_idx(const t_cont<t_val>& vec, t_func fkt)
{
	t_cont<std::size_t> vecIdx;
	vecIdx.reserve(vec.size());
	for(std::size_t i=0; i<vec.size(); ++i)
		vecIdx.push_back(i);

	std::sort(vecIdx.begin(), vecIdx.end(),
		[&vec, &fkt] (std::size_t iIdx0, std::size_t iIdx1) -> bool
	{
		return fkt(vec[iIdx0], vec[iIdx1]);
	});

	return vecIdx;
}

template<template<class...> class t_cont = std::vector, class t_val>
t_cont<std::size_t> sorted_idx(const t_cont<t_val>& vec)
{
	auto fkt = [](const t_val& v0, const t_val& v1) -> bool { return v0<v1; };
	return sorted_idx<t_cont, t_val, decltype(fkt)>(vec, fkt);
}

// -----------------------------------------------------------------------------


template<typename T1, typename T2>
void merge_map(std::map<T1, T2>& mapThis, const std::map<T1, T2>& mapOther)
{
	for(const std::pair<T1, T2>& thepair : mapOther)
		mapThis.insert(thepair);
}


template<class t_map>
bool has_map_all_keys(const t_map& map, const std::initializer_list<typename t_map::key_type>& lst)
{
	for(const typename t_map::key_type& key : lst)
		if(map.find(key) == map.end())
			return false;
	return true;
}


// -----------------------------------------------------------------------------

}
#endif
