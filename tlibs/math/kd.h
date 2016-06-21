/*
 * kd trees
 * @author tweber
 * @date jun-2015
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_KD_H__
#define __TLIBS_KD_H__

#include <vector>
#include <list>
#include <algorithm>
#include <iostream>

namespace tl {


template<class T=double>
struct KdNode
{
	unsigned int iAxis = 0;
	std::vector<T> vecMid;
	KdNode<T> *pParent = nullptr;
	KdNode<T> *pLeft = nullptr;
	KdNode<T> *pRight = nullptr;

	void print(std::ostream& ostr, unsigned int iLevel=0) const
	{
		for(unsigned int i=0; i<iLevel; ++i) ostr << "\t";
		ostr << "mid: ";
		for(T t : vecMid)
			ostr << t << ", ";

		if(pLeft)
		{
			ostr << "\n";
			for(unsigned int i=0; i<iLevel; ++i) ostr << "\t";
			ostr << "left:\n";
			pLeft->print(ostr, iLevel+1);
		}
		if(pRight)
		{
			ostr << "\n";
			for(unsigned int i=0; i<iLevel; ++i) ostr << "\t";
			ostr << "right:\n";
			pRight->print(ostr, iLevel+1);
		}

		ostr << "\n";
	}
};


template<class T=double>
class Kd
{
private:
	static KdNode<T>* make_kd(std::list<std::vector<T>>& lstPoints,
		int &iDim, unsigned int iLevel=0)
	{
		const unsigned int iSize = lstPoints.size();
		if(iSize == 0) return nullptr;

		if(iDim < 0)
			iDim = lstPoints.begin()->size();
		const unsigned int iAxis = iLevel % iDim;

		KdNode<T> *pNode = new KdNode<T>;

		if(iSize == 1)
		{
			pNode->vecMid = *lstPoints.begin();
			return pNode;
		}

		lstPoints.sort(
			[iAxis](const std::vector<T>& vec0, const std::vector<T>& vec1) -> bool
			{
				return vec0[iAxis] <= vec1[iAxis];
			});

		typename std::list<std::vector<T>>::iterator iterMid = std::next(lstPoints.begin(), iSize/2);
		pNode->vecMid = *iterMid;
		pNode->iAxis = iAxis;
		std::list<std::vector<T>> lstLeft(lstPoints.begin(), iterMid);
		std::list<std::vector<T>> lstRight(std::next(iterMid), lstPoints.end());

		pNode->pLeft = make_kd(lstLeft, iDim, iLevel+1);
		pNode->pRight = make_kd(lstRight, iDim, iLevel+1);

		if(pNode->pLeft) pNode->pLeft->pParent = pNode;
		if(pNode->pRight) pNode->pRight->pParent = pNode;

		return pNode;
	}

	static T get_radius_sq(const std::vector<T>& vec0, const std::vector<T>& vec1, unsigned int iDim)
	{
		T tRad = T(0);
		for(unsigned int i=0; i<iDim; ++i)
			tRad += (vec0[i]-vec1[i])*(vec0[i]-vec1[i]);
		return tRad;
	}

	static void get_best_match(const KdNode<T>* pNode, const std::vector<T>& vec,
		const KdNode<T>** ppBestNode, T* pRad, unsigned int iDim)
	{
		T tRad = get_radius_sq(pNode->vecMid, vec, iDim);
		if(tRad <= *pRad)
		{
			*pRad = tRad;
			*ppBestNode = pNode;
		}


		T tDistVecCut = vec[pNode->iAxis] - pNode->vecMid[pNode->iAxis];
		T tDistVecCutSq = tDistVecCut*tDistVecCut;

		if(tDistVecCutSq <= *pRad)							// intersects cut line?
		{
			if(pNode->pLeft)
				get_best_match(pNode->pLeft, vec, ppBestNode, pRad, iDim);
			if(tDistVecCutSq <= *pRad && pNode->pRight)		// still intersects cut line?
				get_best_match(pNode->pRight, vec, ppBestNode, pRad, iDim);
		}
		else
		{
			if(tDistVecCut <= 0.)
			{
				if(pNode->pLeft)
					get_best_match(pNode->pLeft, vec, ppBestNode, pRad, iDim);
			}
			else
			{
				if(pNode->pRight)
					get_best_match(pNode->pRight, vec, ppBestNode, pRad, iDim);
			}
		}
	}

	static void clear_kd(KdNode<T> *pNode)
	{
		if(!pNode) return;

		if(pNode->pLeft) clear_kd(pNode->pLeft);
		if(pNode->pRight) clear_kd(pNode->pRight);
		delete pNode;
	}

protected:
	KdNode<T> *m_pNode = nullptr;
	unsigned int m_iDim = 3;
	std::vector<T> m_vecMin, m_vecMax;

public:
	void Unload()
	{
		clear_kd(m_pNode);
		m_pNode = nullptr;
		m_vecMin.clear();
		m_vecMax.clear();
	}

	// alters lstPoints!
	void Load(std::list<std::vector<T>>& lstPoints, int iDim=-1)
	{
		Unload();

		// get min/max
		for(typename std::list<std::vector<T>>::const_iterator iter = lstPoints.begin();
			iter != lstPoints.end(); ++iter)
		{
			int iTheDim = iDim;
			if(iTheDim < 0)
				iTheDim = iter->size();

			if(m_vecMin.size()==0)
			{
				m_vecMin.resize(iTheDim);
				m_vecMax.resize(iTheDim);

				for(int i0=0; i0<iTheDim; ++i0)
					m_vecMin[i0] = m_vecMax[i0] = (*iter)[i0];
			}
			else
			{
				for(int i0=0; i0<iTheDim; ++i0)
				{
					m_vecMin[i0] = std::min(m_vecMin[i0], (*iter)[i0]);
					m_vecMax[i0] = std::max(m_vecMax[i0], (*iter)[i0]);
				}
			}
		}

		m_pNode = make_kd(lstPoints, iDim);
		m_iDim = unsigned(iDim);
	}

	const std::vector<T>& GetNearestNode(const std::vector<T>& vec) const
	{
		static const std::vector<T> vecNull;

		const KdNode<T> *pnodeBest = nullptr;
		const KdNode<T>** ppBestNode = const_cast<const KdNode<T>**>(&pnodeBest);

		T tRad = get_radius_sq(m_pNode->vecMid, vec, m_iDim);
		get_best_match(m_pNode, vec, ppBestNode, &tRad, m_iDim);

		if(pnodeBest)
			return pnodeBest->vecMid;
		return vecNull;
	}

	bool IsPointInGrid(const std::vector<T>& vec) const
	{
		for(unsigned int i=0; i<m_iDim; ++i)
			if(vec[i] < m_vecMin[i] || vec[i] > m_vecMax[i])
				return false;
		return true;
	}

	Kd() = default;
	Kd(std::list<std::vector<T>>& lstPoints) { Load(lstPoints); }
	virtual ~Kd() { Unload(); }

	const KdNode<T>* GetRootNode() const { return m_pNode; }
};

}
#endif
