/*
 * data point interpolation
 *
 * @author: Tobias Weber
 * @date: 25-04-2013
 * @license GPLv2 or GPLv3
 */

#ifndef __INTERPOLATION_H__
#define __INTERPOLATION_H__

#include <cmath>
#include <boost/math/special_functions/binomial.hpp>
#include <vector>
#include <algorithm>
#include <limits>

#include "../math/math.h"
#include "../math/geo.h"
#include "../helper/misc.h"
#include "../log/log.h"
#include "funcmod.h"

namespace tl {

namespace ublas = boost::numeric::ublas;


// see: http://mathworld.wolfram.com/BernsteinPolynomial.html
template<typename T> T bernstein(int i, int n, T t)
{
	T bino = boost::math::binomial_coefficient<T>(n, i);
	return bino * pow(t, i) * pow(1-t, n-i);
}

// see: http://mathworld.wolfram.com/BezierCurve.html
template<typename T>
ublas::vector<T> bezier(const ublas::vector<T>* P, std::size_t N, T t)
{
	if(N==0) return ublas::vector<T>(0);
	const int n = N-1;

	ublas::vector<T> vec(P[0].size());
	for(std::size_t i=0; i<vec.size(); ++i) vec[i] = T(0);

	for(int i=0; i<=n; ++i)
		vec += P[i]*bernstein(i, n, t);

	return vec;
}


// see: http://mathworld.wolfram.com/B-Spline.html
template<typename T> 
T bspline_base(int i, int j, T t, const std::vector<T>& knots)
{
	if(j==0)
	{
		if((knots[i] <= t) && (t < knots[i+1]) && (knots[i]<knots[i+1]))
			return 1.;
		return 0.;
	}

	T val11 = (t - knots[i]) / (knots[i+j]-knots[i]);
	T val12 = bspline_base(i, j-1, t, knots);
	T val1 = val11 * val12;

	T val21 = (knots[i+j+1]-t) / (knots[i+j+1]-knots[i+1]);
	T val22 = bspline_base(i+1, j-1, t, knots);
	T val2 = val21 * val22;

	T val = val1 + val2;
	return val;
}


// see: http://mathworld.wolfram.com/B-Spline.html
template<typename T>
ublas::vector<T> bspline(const ublas::vector<T>* P, std::size_t N, T t, const std::vector<T>& knots)
{
	if(N==0) return ublas::vector<T>(0);
	const int n = N-1;
	const int m = knots.size()-1;
	const int degree = m-n-1;

	ublas::vector<T> vec(P[0].size());
	for(std::size_t i=0; i<vec.size(); ++i) vec[i] = T(0);

	for(int i=0; i<=n; ++i)
		vec += P[i]*bspline_base(i, degree, t, knots);

	return vec;
}


template<typename T=double>
class Bezier : public FunctionModel_param_gen<ublas::vector<T>>
{
	protected:
		ublas::vector<T> *m_pvecs;
		std::size_t m_iN;

	public:
		Bezier(std::size_t N, const T *px, const T *py);
		virtual ~Bezier();

		virtual ublas::vector<T> operator()(T t) const override;
		virtual const char* GetModelName() const override { return "bezier"; };
};


template<typename T=double>
class BSpline : public FunctionModel_param_gen<ublas::vector<T>>
{
	protected:
		ublas::vector<T> *m_pvecs;
		std::size_t m_iN, m_iDegree;
		std::vector<T> m_vecKnots;

	public:
		BSpline(std::size_t N, const T *px, const T *py, unsigned int iDegree=3);
		virtual ~BSpline();

		virtual ublas::vector<T> operator()(T t) const override;
		virtual const char* GetModelName() const override { return "bspline"; };
};


template<typename T>
void find_peaks(std::size_t iLen, const T* px, const T* py, unsigned int iOrder,
	std::vector<T>& vecMaximaX, std::vector<T>& vecMaximaSize, std::vector<T>& vecMaximaWidth)
{
	BSpline<T> spline(iLen, px, py, iOrder);
	const std::size_t iNumSpline = 512;    // TODO: in config

	T *pSplineX = new T[iNumSpline];
	T *pSplineY = new T[iNumSpline];
	T *pSplineDiff = new T[iNumSpline];
	T *pSplineDiff2 = new T[iNumSpline];

	const T* pdyMin = std::min_element(py, py+iLen);

	for(std::size_t iSpline=0; iSpline<iNumSpline; ++iSpline)
	{
		const T dT = T(iSpline) / T(iNumSpline-1);
		ublas::vector<T> vec = spline(dT);

		pSplineX[iSpline] = vec[0];
		pSplineY[iSpline] = vec[1];
	}

	tl::diff(iNumSpline, pSplineX, pSplineY, pSplineDiff);
	tl::diff(iNumSpline, pSplineX, pSplineDiff, pSplineDiff2);
	std::vector<std::size_t> vecZeroes = tl::find_zeroes<T>(iNumSpline, pSplineDiff);


	for(std::size_t iZeroIdx = 0; iZeroIdx<vecZeroes.size(); ++iZeroIdx)
	{
		const std::size_t iZero = vecZeroes[iZeroIdx];

		// minima / saddle points
		if(pSplineDiff2[iZero] >= 0.)
			continue;

		vecMaximaX.push_back(pSplineX[iZero]);

		int iMinIdxLeft = -1;
		int iMinIdxRight = -1;
		if(iZeroIdx > 0)
			iMinIdxLeft = vecZeroes[iZeroIdx-1];
		if(iZeroIdx+1 < vecZeroes.size())
			iMinIdxRight = vecZeroes[iZeroIdx+1];

		T dHeight = 0.;
		T dWidth = 0.;
		T dDiv = 0.;

		// minimum left of the peak
		if(iMinIdxLeft>=0)
		{
			dHeight += (pSplineY[iZero]-pSplineY[iMinIdxLeft]);
			dWidth += fabs((pSplineX[iZero]-pSplineX[iMinIdxLeft]));
			dDiv += 1.;
		}

		// minimum right of the peak
		if(iMinIdxRight>=0)
		{
			dHeight += (pSplineY[iZero]-pSplineY[iMinIdxRight]);
			dWidth += fabs((pSplineX[iZero]-pSplineX[iMinIdxRight]));
			dDiv += 1.;
		}

		// no adjacent minima...
		if(iMinIdxLeft<0 && iMinIdxRight<0)
		{
			dHeight = pSplineY[iZero]- *pdyMin;
			dWidth = (px[iLen-1] - px[0]) / 10.;	// guess something...
			dDiv = 1.;
		}

		if(dDiv != 0.)
		{
			dHeight /= dDiv;
			dWidth /= dDiv;
		}

		vecMaximaSize.push_back(dHeight);
		vecMaximaWidth.push_back(dWidth);
	}

	tl::sort_3<typename std::vector<T>::iterator>(vecMaximaSize.begin(), vecMaximaSize.end(),
						vecMaximaWidth.begin(), vecMaximaX.begin());
	std::reverse(vecMaximaSize.begin(), vecMaximaSize.end());
	std::reverse(vecMaximaWidth.begin(), vecMaximaWidth.end());
	std::reverse(vecMaximaX.begin(), vecMaximaX.end());



	std::ostringstream ostrDbg;
	ostrDbg << "Prefitter found peaks at: ";
	for(T dValX : vecMaximaX)
		ostrDbg << dValX << ", ";
	log_debug(ostrDbg.str());


	delete[] pSplineX;
	delete[] pSplineY;
	delete[] pSplineDiff;
	delete[] pSplineDiff2;
}



template<class T>
Bezier<T>::Bezier(std::size_t N, const T *px, const T *py)
	: m_pvecs(0), m_iN(N)
{
	m_pvecs = new ublas::vector<T>[m_iN];

	for(std::size_t i=0; i<N; ++i)
	{
		m_pvecs[i].resize(2);
		m_pvecs[i][0] = px[i];
		m_pvecs[i][1] = py[i];
	}

	//auto MinMax = boost::minmax_element(px, px+N);
	//m_dMin = *MinMax.first;
	//m_dMax = *MinMax.second;
}

template<class T>
Bezier<T>::~Bezier()
{
	if(m_pvecs) delete[] m_pvecs;
}

template<class T>
ublas::vector<T> Bezier<T>::operator()(T t) const
{
	return tl::bezier<T>(m_pvecs, m_iN, t);
}



template<class T>
BSpline<T>::BSpline(std::size_t N, const T *px, const T *py, unsigned int iDegree)
	: m_pvecs(0), m_iN(N), m_iDegree(iDegree)
{
	m_pvecs = new ublas::vector<T>[m_iN];

	for(std::size_t i=0; i<m_iN; ++i)
	{
		m_pvecs[i].resize(2);
		m_pvecs[i][0] = px[i];
		m_pvecs[i][1] = py[i];
	}

	std::size_t iM = m_iDegree + m_iN + 1;
	m_vecKnots.resize(iM);

	const T eps = std::numeric_limits<T>::epsilon();

	// set knots to uniform, nonperiodic B-Spline
	for(unsigned int i=0; i<m_iDegree+1; ++i)
		m_vecKnots[i] = 0.+i*eps;
	for(unsigned int i=iM-m_iDegree-1; i<iM; ++i)
		m_vecKnots[i] = 1.-i*eps;
	for(unsigned int i=m_iDegree+1; i<iM-m_iDegree-1; ++i)
		m_vecKnots[i] = T(i+1-m_iDegree-1) / T(iM-2*m_iDegree-2 + 1);

	//for(unsigned int i=0; i<iM; ++i)
	//	m_vecKnots[i] = T(i) / T(iM-1);
}

template<class T>
BSpline<T>::~BSpline()
{
	if(m_pvecs) delete[] m_pvecs;
}

template<class T>
ublas::vector<T> BSpline<T>::operator()(T t) const
{
	if(m_iN==0)
	{
		ublas::vector<T> vecNull(2);
		vecNull[0] = vecNull[1] = 0.;
		return vecNull;
	}

	ublas::vector<T> vec = tl::bspline<T>(m_pvecs, m_iN, t, m_vecKnots);

	// remove epsilon dependence
	if(t<=0.) vec = m_pvecs[0];
	if(t>=1.) vec = m_pvecs[m_iN-1];

	return vec;
}

}

#endif
