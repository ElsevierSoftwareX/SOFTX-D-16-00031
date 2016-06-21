/*
 * abstract function model base classes
 * @author tweber
 * @date 2013
 * @license GPLv2 or GPLv3
 */

#ifndef __FUNC_MOD_H__
#define __FUNC_MOD_H__

#include <boost/numeric/ublas/vector.hpp>

namespace tl {

// parametric function
template<class t_vec = boost::numeric::ublas::vector<double>, class T = typename t_vec::value_type>
class FunctionModel_param_gen
{
public:
	virtual ~FunctionModel_param_gen() = default;

	// t = 0..1
	virtual t_vec operator()(T t) const = 0;
	virtual const char* GetModelName() const = 0;
};

typedef class FunctionModel_param_gen<boost::numeric::ublas::vector<double>> FunctionModel_param;


// ----------------------------------------------------------------------------


// explicit function
template<class T = double> class FunctionModel_gen
{
public:
	virtual ~FunctionModel_gen() = default;

	virtual T operator()(T x) const = 0;
	virtual const char* GetModelName() const = 0;
};

// synonym
template<class T=double> using FunctionModel = class FunctionModel_gen<T>;


// ----------------------------------------------------------------------------


// explicit function with multiple internal parameter sets
template<class T = double> class FunctionModel_multi_gen : public FunctionModel_gen<T>
{
public:
	virtual ~FunctionModel_multi_gen() = default;

	virtual std::size_t GetParamSetCount() const = 0;
	virtual void SetParamSet(std::size_t iSet) = 0;
};

// synonym
template<class T=double> using FunctionModel_multi = class FunctionModel_multi_gen<T>;


// ----------------------------------------------------------------------------


// interface for n dimensional function
template<class T = double> class FunctionModel_nd_gen
{
protected:

public:
	virtual ~FunctionModel_nd_gen() = default;

	virtual T operator()(const T* px) const = 0;
	virtual const char* GetModelName() const = 0;
};

// synonyme
template<class T=double> using FunctionModel_nd = class FunctionModel_nd_gen<T>;



}

#endif
