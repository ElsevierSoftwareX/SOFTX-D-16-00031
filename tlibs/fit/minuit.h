/**
 * Minuit interface
 *
 * @author Tobias Weber
 * @date April 2012
 * @license GPLv2 or GPLv3
 *
 * @desc general fitter structure (i.e. function => chi^2 calculation => calling
 * 	minuit) originally based on the examples in the Minuit user's guide:
 * 	http://seal.cern.ch/documents/minuit/mnusersguide.pdf
 */

#ifndef __MINUIT_IFACE_H__
#define __MINUIT_IFACE_H__

#include <Minuit2/FCNBase.h>
#include <Minuit2/MnFcn.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>
#include <type_traits>
#include <limits>
#include <cmath>
#include <boost/numeric/ublas/vector.hpp>

#include "funcmod.h"
#include "../helper/misc.h"
#include "../log/log.h"


namespace tl {
//using t_real_min = double;
using t_real_min = typename std::result_of<
	decltype(&ROOT::Minuit2::MnFcn::Up)(ROOT::Minuit2::MnFcn) >::type;


class MinuitFuncModel : public FunctionModel<t_real_min>
{
public:
	virtual ~MinuitFuncModel() = default;

	virtual bool SetParams(const std::vector<t_real_min>& vecParams) = 0;
	virtual t_real_min operator()(t_real_min x) const override = 0;

	virtual MinuitFuncModel* copy() const = 0;
	virtual std::string print(bool bFillInSyms=true) const /*= 0;*/ { return ""; }

	virtual const char* GetModelName() const override /*= 0;*/ { return "MinuitFuncModel"; }
	virtual std::vector<std::string> GetParamNames() const = 0;
	virtual std::vector<t_real_min> GetParamValues() const = 0;
	virtual std::vector<t_real_min> GetParamErrors() const = 0;

	friend std::ostream& operator<<(std::ostream& ostr, const MinuitFuncModel& fkt)
	{
		ostr << fkt.print();
		return ostr;
	}
};


template<class t_real = t_real_min>
class MinuitMultiFuncModel : public FunctionModel_multi<t_real_min>
{
public:
	virtual ~MinuitMultiFuncModel() = default;

	virtual bool SetParams(const std::vector<t_real_min>& vecParams) = 0;
	virtual t_real_min operator()(t_real_min x) const override = 0;

	virtual MinuitMultiFuncModel* copy() const = 0;
	virtual std::string print(bool bFillInSyms=true) const /*= 0;*/ { return ""; }

	virtual const char* GetModelName() const override /*= 0;*/ { return "MinuitMultiFuncModel"; }
	virtual std::vector<std::string> GetParamNames() const = 0;
	virtual std::vector<t_real_min> GetParamValues() const = 0;
	virtual std::vector<t_real_min> GetParamErrors() const = 0;

	virtual void SetParamSet(std::size_t iSet) override /*= 0;*/ {}
	virtual std::size_t GetParamSetCount() const override /*= 0;*/ { return 1; }
	// optional intrinsic measured values for multi-parameter functions
	virtual std::size_t GetExpLen() const /*= 0;*/ { return 0; }
	virtual const t_real* GetExpX() const /*= 0;*/ { return nullptr; }
	virtual const t_real* GetExpY() const /*= 0*/ { return nullptr; }
	virtual const t_real* GetExpDY() const /*= 0*/ { return nullptr; }

	friend std::ostream& operator<<(std::ostream& ostr, const MinuitMultiFuncModel<t_real>& fkt)
	{
		ostr << fkt.print();
		return ostr;
	}
};


class MinuitFuncModel_nd : public FunctionModel_nd<t_real_min>
{
public:
	virtual ~MinuitFuncModel_nd() = default;

	virtual std::size_t GetDim() const = 0;

	virtual bool SetParams(const std::vector<t_real_min>& vecParams) = 0;
	virtual t_real_min operator()(const t_real_min* px) const = 0;

	virtual MinuitFuncModel_nd* copy() const = 0;
	virtual std::string print(bool bFillInSyms=true) const = 0;

	virtual const char* GetModelName() const = 0;


	friend std::ostream& operator<<(std::ostream& ostr, const MinuitFuncModel_nd& fkt)
	{
		ostr << fkt.print();
		return ostr;
	}
};


// ----------------------------------------------------------------------------


// generic chi^2 calculation
template<class t_real = t_real_min>
class Chi2Function_gen : public ROOT::Minuit2::FCNBase
{
protected:
	const MinuitFuncModel *m_pfkt;

	std::size_t m_uiLen;
	const t_real* m_px;
	const t_real* m_py;
	const t_real* m_pdy;

	t_real_min m_dSigma = 1.;
	bool m_bDebug = 0;

public:
	Chi2Function_gen(const MinuitFuncModel* fkt=0,
		std::size_t uiLen=0, const t_real *px=0,
		const t_real *py=0, const t_real *pdy=0)
		: m_pfkt(fkt), m_uiLen(uiLen), m_px(px), m_py(py), m_pdy(pdy)
	{}

	virtual ~Chi2Function_gen() {}

	/*
	 * chi^2 calculation
	 * based on the example in the Minuit user's guide:
	 * http://seal.cern.ch/documents/minuit/mnusersguide.pdf
	 */
	t_real_min chi2(const std::vector<t_real_min>& vecParams) const
	{
		// cannot operate on m_pfkt directly, because Minuit
		// uses more than one thread!
		std::unique_ptr<MinuitFuncModel> uptrFkt(m_pfkt->copy());
		MinuitFuncModel* pfkt = uptrFkt.get();

		/*bool bParamsOk = */pfkt->SetParams(vecParams);
		//if(!bParamsOk)
		//	return std::numeric_limits<t_real_min>::max();

		return tl::chi2<t_real_min, decltype(*pfkt), const t_real*>(*pfkt, m_uiLen, m_px, m_py, m_pdy);
	}

	virtual t_real_min Up() const override { return m_dSigma*m_dSigma; }

	virtual t_real_min operator()(const std::vector<t_real_min>& vecParams) const override
	{
		t_real_min dChi2 = chi2(vecParams);
		if(m_bDebug) tl::log_debug("Chi2 = ", dChi2);
		return dChi2;
	}

	void SetSigma(t_real_min dSig) { m_dSigma = dSig; }
	t_real_min GetSigma() const { return m_dSigma; }

	void SetDebug(bool b) { m_bDebug = b; }
};

// for most cases data type of measured values and internal data type is the same: t_real_min
using Chi2Function = Chi2Function_gen<t_real_min>;


// ----------------------------------------------------------------------------


/*
 * chi^2 calculation using multiple simultaneous functions where each
 * function can additionally have multiple parameter sets
 */
template<class t_real = t_real_min, template<class...> class t_cont=std::vector>
class Chi2Function_mult_gen : public ROOT::Minuit2::FCNBase
{
protected:
	t_cont<const MinuitMultiFuncModel<t_real>*> m_vecFkt;

	t_cont<std::size_t> m_vecLen;
	t_cont<const t_real*> m_vecpX;
	t_cont<const t_real*> m_vecpY;
	t_cont<const t_real*> m_vecpDY;

	t_real_min m_dSigma = 1.;
	bool m_bDebug = 0;

public:
	Chi2Function_mult_gen() = default;
	virtual ~Chi2Function_mult_gen() = default;

	void AddFunc(const MinuitMultiFuncModel<t_real>* pMod, std::size_t iNumDat,
		const t_real *pX, const t_real *pY, const t_real *pdY)
	{
		m_vecFkt.push_back(pMod);
		m_vecLen.push_back(iNumDat);
		m_vecpX.push_back(pX);
		m_vecpY.push_back(pY);
		m_vecpDY.push_back(pdY);
	}

	t_real_min chi2(std::size_t iFkt, const std::vector<t_real_min>& vecParams) const
	{
		std::unique_ptr<MinuitMultiFuncModel<t_real>> uptrFkt(m_vecFkt[iFkt]->copy());
		MinuitMultiFuncModel<t_real>* pfkt = uptrFkt.get();

		const std::size_t iNumParamSets = pfkt->GetParamSetCount();
		t_real_min dChi = t_real_min(0);
		for(std::size_t iParamSet=0; iParamSet<iNumParamSets; ++iParamSet)
		{
			pfkt->SetParamSet(iParamSet);
			std::size_t iLen = pfkt->GetExpLen();
			const t_real* pX = pfkt->GetExpX();
			const t_real* pY = pfkt->GetExpY();
			const t_real* pDY = pfkt->GetExpDY();

			// default experimental values if none are given in the model
			if(!pX || !pY || !pDY)
			{
				iLen = m_vecLen[iFkt];
				pX = m_vecpX[iFkt];
				pY = m_vecpY[iFkt];
				pDY = m_vecpDY[iFkt];
			}

			pfkt->SetParams(vecParams);

			t_real_min dSingleChi = tl::chi2<t_real_min, decltype(*pfkt), const t_real*>
				(*pfkt, iLen, pX, pY, pDY);
			dChi += dSingleChi;

			if(m_bDebug && iNumParamSets>1)
				tl::log_debug("Function ", iParamSet, " chi2 = ", dSingleChi);
		}
		dChi /= t_real_min(iNumParamSets);
		return dChi;
	}

	t_real_min chi2(const std::vector<t_real_min>& vecParams) const
	{
		t_real_min dChi = t_real_min(0);
		for(std::size_t iFkt=0; iFkt<m_vecFkt.size(); ++iFkt)
		{
			const t_real_min dSingleChi = chi2(iFkt, vecParams);
			dChi += dSingleChi;

			if(m_bDebug && m_vecFkt.size()>1)
				tl::log_debug("Function ", iFkt, " chi2 = ", dSingleChi);
		}
		dChi /= t_real_min(m_vecFkt.size());
		return dChi;
	}

	virtual t_real_min Up() const override { return m_dSigma*m_dSigma; }

	virtual t_real_min operator()(const std::vector<t_real_min>& vecParams) const override
	{
		t_real_min dChi2 = chi2(vecParams);
		if(m_bDebug) tl::log_debug("Total chi2 = ", dChi2);
		return dChi2;
	}

	void SetSigma(t_real_min dSig) { m_dSigma = dSig; }
	t_real_min GetSigma() const { return m_dSigma; }

	void SetDebug(bool b) { m_bDebug = b; }
};

// for most cases data type of measured values and internal data type is the same: t_real_min
using Chi2FunctionMult = Chi2Function_mult_gen<t_real_min, std::vector>;


// ----------------------------------------------------------------------------


// in n dimensions
class Chi2Function_nd : public ROOT::Minuit2::FCNBase
{
protected:
	const MinuitFuncModel_nd *m_pfkt;
	std::size_t m_uiDim;

	std::size_t m_uiLen;
	std::vector<const t_real_min*> m_vecpx;

	const t_real_min* m_py;
	const t_real_min* m_pdy;

	bool m_bDebug = 0;

public:
	Chi2Function_nd(const MinuitFuncModel_nd* fkt=0,
		std::size_t uiLen=0, const t_real_min **ppx=0,
		const t_real_min *py=0, const t_real_min *pdy=0)
		: m_pfkt(fkt), m_uiDim(fkt->GetDim()), m_uiLen(uiLen), m_py(py), m_pdy(pdy)
	{
		m_vecpx.resize(m_uiDim);

		for(std::size_t i=0; i<m_uiDim; ++i)
			m_vecpx[i] = ppx[i];
	}

	virtual ~Chi2Function_nd() {}
	t_real_min chi2(const std::vector<t_real_min>& vecParams) const
	{
		std::unique_ptr<MinuitFuncModel_nd> uptrFkt(m_pfkt->copy());
		MinuitFuncModel_nd* pfkt = uptrFkt.get();

		/*bool bParamsOk = */pfkt->SetParams(vecParams);

		std::unique_ptr<t_real_min[]> uptrX(new t_real_min[m_uiDim]);

		t_real_min dChi2 = 0.;
		for(std::size_t i=0; i<m_uiLen; ++i)
		{
			for(std::size_t iX=0; iX<m_uiDim; ++iX)
				uptrX[iX] = m_vecpx[iX][i];

			t_real_min d = (*pfkt)(uptrX.get()) - m_py[i];
			t_real_min dy = m_pdy ? m_pdy[i] : 0.1*d;	// assume 10% error if none given
			if(fabs(dy) < std::numeric_limits<t_real_min>::min())
				dy = std::numeric_limits<t_real_min>::min();

			d /= dy;
			dChi2 += d*d;
		}
		return dChi2;
	}


	virtual t_real_min Up() const override
	{
		// 1. for chi^2
		return 1.;
	}

	virtual t_real_min operator()(const std::vector<t_real_min>& vecParams) const override
	{
		t_real_min dChi2 = chi2(vecParams);
		if(m_bDebug) tl::log_debug("Chi2 = ", dChi2);
		return dChi2;
	}

	void SetDebug(bool b) { m_bDebug = b; }
};
}

#endif
