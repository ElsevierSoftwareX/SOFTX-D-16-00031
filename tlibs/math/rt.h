/*
 * wrapper for boost r*trees
 * @author tweber
 * @date oct-2015
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_RT_H__
#define __TLIBS_RT_H__

#include <vector>
#include <list>
#include <algorithm>
#include <type_traits>
//#include <iostream>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>


namespace tl {

namespace geo = boost::geometry;


// ----------------------------------------------------------------------------
// static getter and setter loops for geo point
template<class t_point, class t_iter, std::size_t I, std::size_t MAX>
struct _rt_set_pt
{
	void operator()(t_point& pt, t_iter iter)
	{
		geo::set<I>(pt, *iter);

		_rt_set_pt<t_point, t_iter, I+1, MAX> s;
		s(pt, ++iter);
	}
};
template<class t_point, class t_iter, std::size_t MAX>
struct _rt_set_pt<t_point, t_iter, MAX, MAX>
{
	void operator()(t_point& pt, t_iter iter) {}
};

template<class t_point, class t_iter, std::size_t I, std::size_t MAX>
struct _rt_get_pt
{
	void operator()(const t_point& pt, t_iter iter)
	{
		*iter = geo::get<I>(pt);

		_rt_get_pt<t_point, t_iter, I+1, MAX> g;
		g(pt, ++iter);
	}
};
template<class t_point, class t_iter, std::size_t MAX>
struct _rt_get_pt<t_point, t_iter, MAX, MAX>
{
	void operator()(const t_point& pt, t_iter iter) {}
};
// ----------------------------------------------------------------------------


template<class T=double, std::size_t IDIM=3, std::size_t MAX_ELEMS=32>
class Rt
{
public:
	using t_point = geo::model::point<T, IDIM, geo::cs::cartesian>;
	using t_rest = std::vector<T>;
	using t_node = std::pair<t_point, t_rest>;

protected:
	std::vector<T> m_vecMin, m_vecMax;
	geo::index::rtree<t_node, geo::index::rstar<MAX_ELEMS>> *m_prt = nullptr;

public:
	void Unload()
	{
		m_vecMin.clear();
		m_vecMax.clear();

		if(m_prt) { delete m_prt; m_prt = nullptr; }
	}

	void Load(const std::list<std::vector<T>>& lstPoints)
	{
		Unload();

		// get min/max
		for(typename std::list<std::vector<T>>::const_iterator iter = lstPoints.begin();
			iter != lstPoints.end(); ++iter)
		{
			if(m_vecMin.size()==0)
			{
				m_vecMin.resize(IDIM);
				m_vecMax.resize(IDIM);

				for(std::size_t i0=0; i0<IDIM; ++i0)
					m_vecMin[i0] = m_vecMax[i0] = (*iter)[i0];
			}
			else
			{
				for(std::size_t i0=0; i0<IDIM; ++i0)
				{
					m_vecMin[i0] = std::min(m_vecMin[i0], (*iter)[i0]);
					m_vecMax[i0] = std::max(m_vecMax[i0], (*iter)[i0]);
				}
			}
		}


		m_prt = new typename std::remove_pointer<decltype(m_prt)>::type();

		for(const std::vector<T>& vec : lstPoints)
		{
			t_point pt;
			_rt_set_pt<t_point, typename std::vector<T>::const_iterator, std::size_t(0), IDIM> s;
			s(pt, vec.begin());

			t_rest vecRest;
			vecRest.reserve(vec.size() - IDIM);
			std::copy(vec.begin()+IDIM, vec.end(), std::back_inserter(vecRest));

			m_prt->insert(t_node(pt, vecRest));
		}
	}

	std::list<std::vector<T>> GetNearestNodes(const std::vector<T>& vec, std::size_t iNum=1) const
	{
		t_point pt;
		_rt_set_pt<t_point, typename std::vector<T>::const_iterator, std::size_t(0), IDIM> s;
		s(pt, vec.begin());

		std::vector<t_node> vecRes;
		m_prt->query(geo::index::nearest(pt, iNum), std::back_inserter(vecRes));


		std::list<std::vector<T>> lstRet;
		for(const t_node& nd : vecRes)
		{
			std::vector<T> vec;
			vec.reserve(IDIM + nd.second.size());

			_rt_get_pt<t_point, std::back_insert_iterator<std::vector<T>>, std::size_t(0), IDIM> g;
			g(nd.first, std::back_inserter(vec));

			std::copy(nd.second.begin(), nd.second.end(), std::back_inserter(vec));

			lstRet.push_back(vec);
		}

		return lstRet;
	}

	std::vector<T> GetNearestNode(const std::vector<T>& vec) const
	{
		return *GetNearestNodes(vec, 1).begin();
	}

	bool IsPointInGrid(const std::vector<T>& vec) const
	{
		for(std::size_t i=0; i<IDIM; ++i)
			if(vec[i] < m_vecMin[i] || vec[i] > m_vecMax[i])
				return false;
		return true;
	}

	Rt() = default;
	Rt(std::list<std::vector<T>>& lstPoints) { Load(lstPoints); }
	virtual ~Rt() { Unload(); }
};

}
#endif
