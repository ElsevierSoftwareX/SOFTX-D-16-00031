/*
 * Calculation of first Brillouin zone
 * @author Tobias Weber
 * @date jun-2014
 * @license GPLv2 or GPLv3
 */

#ifndef __BZ_H__
#define __BZ_H__

#include "../helper/exception.h"
#include "math.h"
#include "geo.h"
#include "../log/log.h"
#include <vector>

namespace tl {

template<typename T=double>
class Brillouin2D
{
	public:
		template<typename _T> using t_vec = ublas::vector<_T>;
		template<typename _T> using t_vecpair = std::pair<t_vec<_T>, t_vec<_T> >;
		template<typename _T> using t_vertices = std::vector<t_vecpair<_T> >;

	protected:
		t_vec<T> m_vecCentralReflex;
		std::vector<t_vec<T>> m_vecNeighbours;
		std::vector<t_vecpair<T> > m_vecVertices;
		bool m_bValid = 1;
		bool m_bHasCentralPeak = 0;

		const T eps = 0.001;

	protected:
		const t_vecpair<T>* GetNextVertexPair(const t_vertices<T>& vecPts, std::size_t *pIdx=0)
		{
			const t_vec<T>& vecLast = m_vecVertices.rbegin()->second;
			for(std::size_t iPt=0; iPt<vecPts.size(); ++iPt)
			{
				const t_vecpair<T>& vecp = vecPts[iPt];

				if(vec_equal(vecp.first, vecLast, eps))
				{
					if(pIdx) *pIdx = iPt;
					return &vecp;
				}
			}

			if(vecPts.size())
			{
				log_err("Invalid vertices in Brillouin zone.");
				m_bValid = 0;
			}
			return 0;
		}

	public:
		Brillouin2D() {}
		virtual ~Brillouin2D() {}

		const t_vertices<T>& GetVertices() const { return m_vecVertices; }
		const t_vec<T>& GetCentralReflex() const { return m_vecCentralReflex; }
		const std::vector<t_vec<T>>& GetNeighbours() const { return m_vecNeighbours; }

		void Clear()
		{
			m_vecCentralReflex.clear();
			m_vecNeighbours.clear();
			m_vecVertices.clear();
			m_bValid = 0;
		}

		bool IsValid() const { return m_bValid; }

		void SetCentralReflex(const t_vec<T>& vec)
		{
			if(vec.size() != 2)
				throw Err("Brillouin2D needs 2d vectors.");

			m_vecCentralReflex = vec;
			m_bHasCentralPeak = 1;
		}
		void AddReflex(const t_vec<T>& vec)
		{
			if(vec.size() != 2)
				throw Err("Brillouin2D needs 2d vectors.");

			m_vecNeighbours.push_back(vec);
		}

		void CalcBZ()
		{
			if(!m_bHasCentralPeak) return;
			m_bValid = 1;

			// calculate perpendicular lines
			std::vector<Line<T> > vecMiddlePerps;
			vecMiddlePerps.reserve(m_vecNeighbours.size());

			t_vertices<T> vecPts;

			for(const t_vec<T>& vecN : m_vecNeighbours)
			{
				Line<T> line(m_vecCentralReflex, vecN-m_vecCentralReflex);
				Line<T> lineperp;
				if(!line.GetMiddlePerp(lineperp))
					continue;

				vecMiddlePerps.push_back(lineperp);
			}


			// calculate intersections
			for(std::size_t iThisLine=0; iThisLine<vecMiddlePerps.size(); ++iThisLine)
			{
				const Line<T>& lineThis = vecMiddlePerps[iThisLine];
				T tPos = std::numeric_limits<T>::max();
				T tNeg = -tPos;

				for(std::size_t iOtherLine=0; iOtherLine<vecMiddlePerps.size(); ++iOtherLine)
				{
					if(iThisLine == iOtherLine)
						continue;

					T t;
					if(lineThis.IsParallel(vecMiddlePerps[iOtherLine], eps))
						continue;
					if(!lineThis.intersect(vecMiddlePerps[iOtherLine], t))
						continue;

					if(t>0.)
					{
						if(t < tPos) tPos = t;
					}
					else
					{
						if(t > tNeg) tNeg = t;
					}
				}

				t_vec<T> vecUpper = lineThis(tPos);
				t_vec<T> vecLower = lineThis(tNeg);
				if(vec_equal(vecUpper, vecLower, eps))
					continue;

				vecPts.push_back(t_vecpair<T>(vecUpper, vecLower));
			}


			// remove unnecessary vertices
			for(const Line<T>& line : vecMiddlePerps)
			{
				bool bSideReflex = line.GetSide(m_vecCentralReflex);

				for(std::size_t iPt=0; iPt<vecPts.size(); ++iPt)
				{
					t_vec<T>& vecUpper = vecPts[iPt].first;
					t_vec<T>& vecLower = vecPts[iPt].second;

					T tDistUpper=T(0), tDistLower=T(0);
					bool bSideUpper = line.GetSide(vecUpper, &tDistUpper);
					bool bSideLower = line.GetSide(vecLower, &tDistLower);

					if((bSideUpper!=bSideReflex && tDistUpper>eps) ||
						(bSideLower!=bSideReflex && tDistLower>eps))
					{
						vecPts.erase(vecPts.begin()+iPt);
						--iPt;
						continue;
					}
				}
			}

			if(vecPts.size() == 0)
				return;

			// sort vertices
			m_vecVertices.clear();
			m_vecVertices.reserve(vecPts.size());
			m_vecVertices.push_back(*vecPts.begin());


			vecPts.erase(vecPts.begin());

			std::size_t iIdx;
			while(const t_vecpair<T>* pPair = GetNextVertexPair(vecPts, &iIdx))
			{
				m_vecVertices.push_back(*pPair);
				vecPts.erase(vecPts.begin()+iIdx);
			}
		}
};

}

#endif
