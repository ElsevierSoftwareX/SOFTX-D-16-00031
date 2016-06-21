/**
 * array helpers
 *
 * @author: tweber
 * @date: nov-2015
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_ARRAYS_H__
#define __TLIBS_ARRAYS_H__

#include <cstddef>
#include <vector>
#include <list>
#include <array>
#include <tuple>
#include <algorithm>


namespace tl {


/**
 * Minimalistic wrapper for plain old arrays
 */
template<class T>
class wrapper_array
{
	public:
		using value_type = T;
		using size_type = std::size_t;

		using iterator = T*;
		using const_iterator = const T*;

		using reference = T&;
		using const_reference = const T&;


	protected:
		T* m_pt = nullptr;
		size_type m_len = 0;

	public:
		wrapper_array(T* p, size_type len)
			: m_pt(p), m_len(len)
		{}

		wrapper_array() = default;
		~wrapper_array() = default;

		size_type size() const { return m_len; }

		iterator begin() { return m_pt; }
		iterator end() { return m_pt+m_len; }
		const_iterator begin() const { return m_pt; }
		const_iterator end() const { return m_pt+m_len; }

		reference operator[](size_type i) { return m_pt[i]; }
		const_reference operator[](size_type i) const { return m_pt[i]; }
};



// ----------------------------------------------------------------------------
// conversions

/**
 * Converts a t_cont_in to a t_cont_out
 */
template<template<class...> class t_cont_out, class t_cont_in>
t_cont_out<typename t_cont_in::value_type> convert_containers(const t_cont_in& cont1)
{
	using T = typename t_cont_in::value_type;

	t_cont_out<T> cont2;

	for(const T& val : cont1)
		cont2.push_back(val);

	return cont2;
}

template<typename T>
std::list<T> vector_to_list(const std::vector<T>& vec)
{
	return convert_containers<std::list, std::vector<T>>(vec);
}

template<typename T>
T* vec_to_array(const std::vector<T>& vec)
{
	T* t_arr = new T[vec.size()];

	std::size_t i=0;
	for(const T& t : vec)
		t_arr[i++] = t;

	return t_arr;
}

// ----------------------------------------------------------------------------


// sort tuple-vector
template<const std::size_t isortidx,
	class... Ts,
	template<class...> class t_cont = std::vector>
void sorttuples(t_cont<std::tuple<Ts...> >& vec)
{
	std::sort(vec.begin(), vec.end(),
		[](const std::tuple<Ts...>& tup1, const std::tuple<Ts...>& tup2) -> bool
		{ return std::get<isortidx>(tup1) < std::get<isortidx>(tup2);});
}


// ----------------------------------------------------------------------------

template<class t_cont>
t_cont arrayunion(const std::initializer_list<t_cont>& lst)
{
	t_cont contRet;
	for(const t_cont& cont : lst)
		contRet.insert(contRet.end(), cont.begin(), cont.end());
	return contRet;
}


// ----------------------------------------------------------------------------


template<class t_to, class t_from, template<class...> class t_cont,
	bool bIsEqu = std::is_same<t_from, t_to>::value>
struct container_cast {};

template<class t_to, class t_from, template<class...> class t_cont>
struct container_cast<t_to, t_from, t_cont, 1>
{
	const t_cont<t_to>& operator()(const t_cont<t_from>& vec) const
	{ return vec; }
};

template<class t_to, class t_from, template<class...> class t_cont>
struct container_cast<t_to, t_from, t_cont, 0>
{
	t_cont<t_to> operator()(const t_cont<t_from>& vec) const
	{
		t_cont<t_to> vecTo;
		for(const t_from& t : vec)
			vecTo.push_back(t_to(t));
		return vecTo;
	}
};

// ----------------------------------------------------------------------------


template<std::size_t iWhichElem = 0,
	template<class...> class t_cont = std::vector,
	class t_pairvec = std::vector<std::pair<double,double>>>
struct vec_from_pairvec
{
	using t_pair = typename t_pairvec::value_type;
	using t_val = typename std::tuple_element<iWhichElem, t_pair>::type;

	t_cont<t_val> operator()(const t_pairvec& vec) const
	{
		t_cont<t_val> vecVal;
		vecVal.reserve(vec.size());

		for(const t_pair& elem : vec)
			vecVal.push_back(std::get<iWhichElem>(elem));

		return vecVal;
	}
};

// ----------------------------------------------------------------------------


}
#endif
