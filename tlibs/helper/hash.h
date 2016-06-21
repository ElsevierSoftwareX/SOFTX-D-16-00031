/*
 * hashes
 * @author tweber
 * @date apr-2015
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIB_HASH_H__
#define __TLIB_HASH_H__

#include <vector>
#include <algorithm>
#include <boost/functional/hash.hpp>

namespace tl {

template<typename T>
std::size_t hash(const T& t)
{
	boost::hash<T> hsh;
	return hsh(t);
}

// order of elements in container matters
template<class t_cont>
std::size_t hash_ordered(const t_cont& cont)
{
	typedef typename t_cont::const_iterator t_iter;
	typedef typename t_cont::value_type T;

	std::size_t iseed = 0;
	boost::hash<T> hsh;

	for(t_iter iter=cont.begin(); iter!=cont.end(); ++iter)
	{
		T t = *iter;
		std::size_t iHsh = hsh(t);

		boost::hash_combine(iseed, iHsh);
	}

	return iseed;
}

// order of elements in container doesn't matter
template<class t_cont>
std::size_t hash_unordered(const t_cont& cont)
{
	typedef typename t_cont::const_iterator t_iter;
	typedef typename t_cont::value_type T;

	std::vector<T> vec;

	for(t_iter iter=cont.begin(); iter!=cont.end(); ++iter)
		vec.push_back(*iter);

	std::sort(vec.begin(), vec.end());
	return hash_ordered<std::vector<T>>(vec);
}

}
#endif
