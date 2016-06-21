/*
 * external fit functions
 * @author tweber
 * @date jan 2014
 * @license GPLv2 or GPLv3
 */

#include "lang/types.h"
#include "string/string.h"
#include "log/log.h"
#include "fit/interpolation.h"
#include "calls_fit.h"
#include "lang/calls.h"
#include "lang/node.h"


template<typename T> using t_stdvec = std::vector<T>;

// --------------------------------------------------------------------------------
// fitting

#include "fit/minuit.h"
#include <algorithm>
#include <exception>

#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnPrint.h>
#include <Minuit2/MnMinos.h>

/*
static inline std::vector<std::string>
convert_string_vector(const std::vector<t_string>& vecStr)
{
	std::vector<std::string> vecRet;
	vecRet.reserve(vecStr.size());

	for(const t_string& str : vecStr)
		vecRet.push_back(WSTR_TO_STR(str));

	return vecRet;
}

static inline std::vector<t_string>
convert_string_vector(const std::vector<std::string>& vecStr)
{
	std::vector<t_string> vecRet;
	vecRet.reserve(vecStr.size());

	for(const std::string& str : vecStr)
		vecRet.push_back(STR_TO_WSTR(str));

	return vecRet;
}*/

class GenericModel : public tl::MinuitFuncModel
{
protected:
	const NodeFunction *m_pFkt = 0;

	ParseInfo *m_pinfo = 0;
	const SymbolTable *m_pCallerSymTab = 0;
	SymbolTable *m_pTable = 0;

	std::string m_strFreeParam;
	std::vector<std::string> m_vecParamNames;
	std::vector<Symbol*> m_vecSyms;

	bool m_bUseVecParams = false;

public:
	GenericModel(const GenericModel& mod)
		: m_pFkt((NodeFunction*)mod.m_pFkt/*->clone()->optimize()*/),
		  m_pinfo(mod.m_pinfo),
		  m_pCallerSymTab(mod.m_pCallerSymTab),
		  m_pTable(new SymbolTable()),
		  m_strFreeParam(mod.m_strFreeParam),
		  m_vecParamNames(mod.m_vecParamNames),
		  m_bUseVecParams(mod.m_bUseVecParams)
	{
		if(m_bUseVecParams)
		{
			m_vecSyms.resize(2);
			m_vecSyms[0] = mod.m_vecSyms[0]->clone();
			m_vecSyms[1] = mod.m_vecSyms[1]->clone();
		}
		else
		{
			m_vecSyms.reserve(m_vecParamNames.size()+1);
			for(unsigned int i=0; i<m_vecParamNames.size()+1; ++i)
				m_vecSyms.push_back(new SymbolReal(mod.m_vecSyms[i]->GetValDouble()));
		}
	}

	GenericModel(const NodeFunction *pFkt, ParseInfo& info,
		 SymbolTable *pCallerSymTab, const std::vector<t_string>* pvecParamNames=0)
		: m_pFkt((NodeFunction*)pFkt/*->clone()->optimize()*/),
		  m_pinfo(&info),
		  m_pCallerSymTab(pCallerSymTab),
		  m_pTable(new SymbolTable())
	{
		std::vector<std::string> vecParams = /*convert_string_vector*/(m_pFkt->GetParamNames());
		m_strFreeParam = vecParams[0];

		if(pvecParamNames)
		{
			m_bUseVecParams = 1;
			m_vecParamNames = *pvecParamNames;

			SymbolArray* pSymParams = new SymbolArray();
			for(unsigned int i=0; i<m_vecParamNames.size(); ++i)
				pSymParams->GetArr().push_back(new SymbolReal(0.));

			m_vecSyms.resize(2);
			m_vecSyms[0] = new SymbolReal(0.);	// free parameter
			m_vecSyms[1] = pSymParams;		// parameter vector

			pSymParams->UpdateIndices();
		}
		else
		{
			m_bUseVecParams = 0;
			m_vecParamNames.resize(vecParams.size()-1);
			std::copy(vecParams.begin()+1, vecParams.end(), m_vecParamNames.begin());

			m_vecSyms.reserve(m_vecParamNames.size()+1);
			for(unsigned int i=0; i<m_vecParamNames.size()+1; ++i)
				m_vecSyms.push_back(new SymbolReal(0.));
		}
	}

	virtual ~GenericModel()
	{
		//if(m_pFkt) { delete m_pFkt; m_pFkt=0; }
		if(m_pTable) { delete m_pTable; m_pTable=0; }
	}

	virtual bool SetParams(const std::vector<tl::t_real_min>& vecParams) override
	{
		if(vecParams.size() != m_vecParamNames.size())
		{
			tl::log_err("Parameter array length mismatch.");
			return 0;
		}

		for(unsigned int iParam=0; iParam<vecParams.size(); ++iParam)
		{
			//const std::string& strName = m_vecParamNames[iParam];
			tl::t_real_min dVal = vecParams[iParam];

			if(m_bUseVecParams)
				((SymbolReal*)((SymbolArray*)m_vecSyms[1])->GetArr()[iParam])->SetVal(dVal);
			else
				((SymbolReal*)m_vecSyms[iParam+1])->SetVal(dVal);
		}

		return 1;
	}

	virtual tl::t_real_min operator()(tl::t_real_min x) const override
	{
		((SymbolReal*)m_vecSyms[0])->SetVal(t_real(x));

		SymbolArray arrArgs;
		arrArgs.SetDontDel(1);
		arrArgs.GetArr() = m_vecSyms;
		arrArgs.UpdateIndices();
		//m_pFkt->SetArgSyms(&m_vecSyms);

		RuntimeInfo runinfo;
		m_pTable->InsertSymbol(T_STR"<args>", &arrArgs);
		Symbol *pSymRet = m_pFkt->eval(*m_pinfo, runinfo, m_pTable);
		m_pTable->RemoveSymbolNoDelete(T_STR"<args>");

		t_real dRetVal = 0.;
		if(pSymRet)
		{
			dRetVal = pSymRet->GetValDouble();
			//if(std::isnan(dRetVal) || std::isinf(dRetVal))
			//	dRetVal = std::numeric_limits<double>::max();
		}
		safe_delete(pSymRet, m_pCallerSymTab, m_pinfo);

		return tl::t_real_min(dRetVal);
	}

	virtual GenericModel* copy() const override
	{ return new GenericModel(*this); }
	virtual std::string print(bool bFillInSyms=true) const override
	{ return "<not implemented>"; }
	virtual const char* GetModelName() const override
	{ return "generic fitter model"; }
	virtual std::vector<std::string> GetParamNames() const override
	{ return m_vecParamNames; }

	virtual std::vector<tl::t_real_min> GetParamValues() const override
	{ throw tl::Err("Called invalid function in generic fitter model"); }
	virtual std::vector<tl::t_real_min> GetParamErrors() const override
	{ throw tl::Err("Called invalid function in generic fitter model"); }
};


static void get_values(const std::vector<t_string>& vecParamNames,
	const Symbol* pSym, std::vector<t_real>& vec, std::vector<bool>& vecActive)
{
	vecActive.resize(vecParamNames.size());

	if(pSym->GetType() == SYMBOL_ARRAY)
	{
		vec = sym_to_vec<t_stdvec>(pSym);

		for(unsigned int i=0; i<vecActive.size(); ++i)
			vecActive[i] = 1;
	}
	else if(pSym->GetType() == SYMBOL_MAP)
	{
		vec.resize(vecParamNames.size());
		std::map<t_string, t_real> mymap = sym_to_map<t_string, t_real>(pSym);

		for(unsigned int iParam=0; iParam<vecParamNames.size(); ++iParam)
		{
			const t_string& strKey = vecParamNames[iParam];

			std::map<t_string, t_real>::iterator iter = mymap.find(strKey);
			if(iter == mymap.end())
			{
				vecActive[iParam] = 0;
			}
			else
			{
				vecActive[iParam] = 1;
				vec[iParam] = iter->second;
			}
		}
	}
}

// fit("function", x, y, yerr, params)
static Symbol* fkt_fit(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	int iDebug = 0;
	bool bDoMinos = 0;
	t_real dSigma = 1.;

	if(vecSyms.size()<4 || !is_vec(vecSyms[1]) || !is_vec(vecSyms[2]) || !is_vec(vecSyms[3]))
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "Invalid arguments for fit.";
		throw tl::Err(ostrErr.str(),0);
	}

	if(vecSyms[0]->GetType() != SYMBOL_STRING)
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "Need a fit function name.";
		throw tl::Err(ostrErr.str(),0);
	}

	const t_string& strFkt = ((SymbolString*)vecSyms[0])->GetVal();
	NodeFunction *pFkt = info.GetFunction(strFkt);
	if(!pFkt)
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "Invalid function \"" << strFkt << "\".";
		throw tl::Err(ostrErr.str(),0);
	}


	std::vector<t_string> vecParamNames;

	bool bUseParamVec = 0, bHasParamMap = 0;
	if(vecSyms.size()==5 && vecSyms[4]->GetType()==SYMBOL_MAP)
		bHasParamMap = 1;

	if(bHasParamMap)
	{
		SymbolMap::t_map& mapSym = ((SymbolMap*)vecSyms[4])->GetMap();

		SymbolMap::t_map::iterator iterParamNames = mapSym.find(SymbolMapKey(T_STR"use_param_vec"));
		if(iterParamNames != mapSym.end())
		{
			vecParamNames = sym_to_vec<t_stdvec, t_string>(iterParamNames->second);
			bUseParamVec = 1;
		}
	}

	GenericModel mod(pFkt, info, pSymTab, bUseParamVec?&vecParamNames:0);

	if(!bUseParamVec)
		vecParamNames = /*convert_string_vector*/(mod.GetParamNames());
	const unsigned int iParamSize = vecParamNames.size();


	std::vector<t_real> vecHints, vecHintsErr;
	std::vector<bool> vecHintsActive, vecHintsErrActive;

	std::vector<t_real> vecLimMin, vecLimMax;
	std::vector<bool> vecLimMinActive, vecLimMaxActive;

	std::vector<t_string> vecFittingSteps;
	std::vector<t_string> vecFixedParams;

	vecLimMinActive.resize(iParamSize);
	vecLimMaxActive.resize(iParamSize);

	// parameter map
	if(bHasParamMap)
	{
		SymbolMap::t_map& mapSym = ((SymbolMap*)vecSyms[4])->GetMap();

		SymbolMap::t_map::iterator iterHints = mapSym.find(SymbolMapKey(T_STR"hints"));
		SymbolMap::t_map::iterator iterHintsErr = mapSym.find(SymbolMapKey(T_STR"hints_errors"));

		SymbolMap::t_map::iterator iterLimitsMin = mapSym.find(SymbolMapKey(T_STR"lower_limits"));
		SymbolMap::t_map::iterator iterLimitsMax = mapSym.find(SymbolMapKey(T_STR"upper_limits"));

		SymbolMap::t_map::iterator iterFixed = mapSym.find(SymbolMapKey(T_STR"fixed"));

		SymbolMap::t_map::iterator iterSteps = mapSym.find(SymbolMapKey(T_STR"steps"));

		SymbolMap::t_map::iterator iterDebug = mapSym.find(SymbolMapKey(T_STR"debug"));

		SymbolMap::t_map::iterator iterSigma = mapSym.find(SymbolMapKey(T_STR"sigma"));
		SymbolMap::t_map::iterator iterErrAnalysis = mapSym.find(SymbolMapKey(T_STR"error_analysis"));



		if(iterHints != mapSym.end())
			get_values(vecParamNames, iterHints->second, vecHints, vecHintsActive);
		if(iterHintsErr != mapSym.end())
			get_values(vecParamNames, iterHintsErr->second, vecHintsErr, vecHintsErrActive);

		if(iterLimitsMin != mapSym.end())
			get_values(vecParamNames, iterLimitsMin->second, vecLimMin, vecLimMinActive);
		if(iterLimitsMax != mapSym.end())
			get_values(vecParamNames, iterLimitsMax->second, vecLimMax, vecLimMaxActive);

		if(iterSteps != mapSym.end())
		{
			vecFittingSteps = sym_to_vec<t_stdvec, t_string>(iterSteps->second);
			for(t_string& strStep : vecFittingSteps)
				strStep = tl::remove_char(strStep, ' ');
		}

		if(iterFixed != mapSym.end())
			vecFixedParams = sym_to_vec<t_stdvec, t_string>(iterFixed->second);

		if(iterDebug != mapSym.end())
			iDebug = iterDebug->second->GetValInt();

		if(iterSigma != mapSym.end())
			dSigma = iterSigma->second->GetValDouble();
		if(iterErrAnalysis != mapSym.end())
			bDoMinos = (iterErrAnalysis->second->GetValInt()!=0);
	}

	bool bFitterDebug = (iDebug>0);


	std::vector<t_real> vecX = sym_to_vec<t_stdvec>(vecSyms[1]);
	std::vector<t_real> vecY = sym_to_vec<t_stdvec>(vecSyms[2]);
	std::vector<t_real> vecYErr = sym_to_vec<t_stdvec>(vecSyms[3]);

	unsigned int iSize = std::min<unsigned int>(vecX.size(), vecY.size());
	iSize = std::min<unsigned int>(iSize, vecYErr.size());

	tl::Chi2Function_gen<t_real> chi2fkt(&mod, iSize, vecX.data(), vecY.data(), vecYErr.data());
	chi2fkt.SetSigma(dSigma);


	ROOT::Minuit2::MnUserParameters params;
	for(unsigned int iParam=0; iParam<iParamSize; ++iParam)
	{
		const t_string& wstrParam = vecParamNames[iParam];
		std::string strParam = WSTR_TO_STR(wstrParam);

		t_real dHint = 0.;
		t_real dErr = 0.;

		if(iParam < vecHints.size())
			dHint = vecHints[iParam];
		if(iParam < vecHintsErr.size())
			dErr = vecHintsErr[iParam];

		params.Add(strParam, dHint, dErr);



		t_real dLimMin = 0.;
		t_real dLimMax = 0.;

		if(iParam < vecLimMin.size())
			dLimMin = vecLimMin[iParam];
		if(iParam < vecLimMax.size())
			dLimMax = vecLimMax[iParam];

		if(vecLimMinActive[iParam] && vecLimMaxActive[iParam])
		{
			params.SetLimits(strParam, dLimMin, dLimMax);
		}
		else if(vecLimMinActive[iParam] && vecLimMaxActive[iParam]==0)
			params.SetLowerLimit(strParam, dLimMin);
		else if(vecLimMinActive[iParam]==0 && vecLimMaxActive[iParam])
			params.SetUpperLimit(strParam, dLimMax);

		if(std::find(vecFixedParams.begin(), vecFixedParams.end(), /*w*/strParam) != vecFixedParams.end())
			params.Fix(strParam);
	}


	bool bValidFit = 1;
	std::vector<ROOT::Minuit2::FunctionMinimum> minis;

	if(vecFittingSteps.size() == 0)
	{
		if(bFitterDebug)
			tl::log_info("Using default fitting steps.");

		minis.reserve(2);

		{
			// step 1: free fit (limited)

			ROOT::Minuit2::MnMigrad migrad(chi2fkt, params, 2);
			ROOT::Minuit2::FunctionMinimum mini = migrad();
			bValidFit = mini.IsValid() && mini.HasValidParameters() && mini.UserState().IsValid();

			for(const t_string& strSym : vecParamNames)
			{
				std::string _strSym = WSTR_TO_STR(strSym);

				params.SetValue(_strSym, mini.UserState().Value(_strSym));
				params.SetError(_strSym, mini.UserState().Error(_strSym));
			}

			minis.push_back(mini);
		}


		{
			// step 2: free fit (unlimited)

			for(const t_string& strSym : vecParamNames)
				params.RemoveLimits(WSTR_TO_STR(strSym));

			ROOT::Minuit2::MnMigrad migrad(chi2fkt, params, 2);
			ROOT::Minuit2::FunctionMinimum mini = migrad();
			bValidFit = mini.IsValid() && mini.HasValidParameters() && mini.UserState().IsValid();

			minis.push_back(mini);
		}
	}
	else	// custom fitting steps
	{
		if(bFitterDebug)
			tl::log_info("Using ", vecFittingSteps.size(), " custom fitting steps.");

		minis.reserve(vecFittingSteps.size());

		// one strip per fitting step
		unsigned int iFitStep = 0;
		for(const t_string& strStep : vecFittingSteps)
		{
			if(bFitterDebug)
				tl::log_info("Performing fitting step ", iFitStep+1, ": ", strStep);

			// one character for each parameter
			for(unsigned int iStepParam = 0; iStepParam<strStep.length(); ++iStepParam)
			{
				t_char cOp = strStep[iStepParam];
				const t_string& strParam = vecParamNames[iStepParam];

				switch(cOp)
				{
					case 'f':			// fix param
						params.Fix(WSTR_TO_STR(strParam));
						break;
					case 'r': 			// release param
						params.Release(WSTR_TO_STR(strParam));
						break;
					case 'u': 			// remove limits
						params.RemoveLimits(WSTR_TO_STR(strParam));
						break;
					case 'x': 			// release param and remove limits
						params.Release(WSTR_TO_STR(strParam));
						params.RemoveLimits(WSTR_TO_STR(strParam));
						break;
					case 'n':			// nop
						break;
					default:
						tl::log_err("Unknow fitting step operation \'", cOp, "\'.");
						break;
				}
			} // ops


			ROOT::Minuit2::MnMigrad migrad(chi2fkt, params, 2);
			ROOT::Minuit2::FunctionMinimum mini = migrad();
			bValidFit = mini.IsValid() && mini.HasValidParameters() && mini.UserState().IsValid();

			for(unsigned int iStepParam = 0; iStepParam<strStep.length(); ++iStepParam)
			{
				char cOp = strStep[iStepParam];
				const t_string& strParam = vecParamNames[iStepParam];
				std::string _strSym = WSTR_TO_STR(strParam);

				// copy fitted values only for non-fixed params
				if(cOp != 'f')
				{
					params.SetValue(_strSym, mini.UserState().Value(_strSym));
					params.SetError(_strSym, mini.UserState().Error(_strSym));
				}
			}

			minis.push_back(mini);
			++iFitStep;
		} // steps
	}




	const ROOT::Minuit2::FunctionMinimum& lastmini = *minis.rbegin();


	std::vector<tl::t_real_min> vecLastParams;
	vecLastParams.reserve(vecParamNames.size());

	SymbolMap *pSymMap = new SymbolMap();
	for(const t_string& strSym : vecParamNames)
	{
		std::string _strSym = WSTR_TO_STR(strSym);

		t_real dVal = lastmini.UserState().Value(_strSym);
		t_real dErr = lastmini.UserState().Error(_strSym);
		dErr = fabs(dErr);

		vecLastParams.push_back(dVal);

		SymbolArray* pArr = new SymbolArray();
		pArr->GetArr().push_back(new SymbolReal(dVal));
		pArr->GetArr().push_back(new SymbolReal(dErr));
		pArr->UpdateIndices();

		pSymMap->GetMap().insert(SymbolMap::t_map::value_type(strSym, pArr));
	}
	pSymMap->UpdateIndices();


	std::vector<std::pair<tl::t_real_min, tl::t_real_min>> vecMinosErrs;
	vecMinosErrs.reserve(iParamSize);

	if(bDoMinos)
	{
		// error analysis (negative (first) and positive (second) errors)
		// TODO: make available for scripts (not just debug output)
		ROOT::Minuit2::MnMinos minos(chi2fkt, lastmini);
		for(unsigned int iParam=0; iParam<iParamSize; ++iParam)
		{
			//const std::string& strCurParam = vecParamNames[iParam];
			std::pair<tl::t_real_min, tl::t_real_min> err = minos(iParam);
			vecMinosErrs.push_back(err);
		}
	}


	tl::t_real_min dChi2 = chi2fkt(vecLastParams);
	tl::t_real_min dChi2red = dChi2 / tl::t_real_min(iSize-(vecLastParams.size()-vecFixedParams.size()));

	if(bFitterDebug)
	{
		unsigned int uiMini=0;
		for(const auto& mini : minis)
		{
			std::ostringstream ostrMini;
			ostrMini << mini;

			tl::log_info("result of user-defined fit step ", (++uiMini));
			tl::log_info(STR_TO_WSTR(ostrMini.str()));
		}

		if(bDoMinos)
		{
			tl::log_info("Minos error analysis");
			for(unsigned int iParam=0; iParam<iParamSize; ++iParam)
			{
				tl::log_info(vecParamNames[iParam], ": lower error: ",
					vecMinosErrs[iParam].first, ", upper error: ",
					vecMinosErrs[iParam].second);
			}
		}

		tl::log_info("Chi^2 = ", dChi2);
		tl::log_info("Chi^2 / ndf = ", dChi2red);
	}

	pSymMap->GetMap()[SymbolMapKey("<model>")] = new SymbolString(strFkt);
	pSymMap->GetMap()[SymbolMapKey("<valid>")] = new SymbolInt(bValidFit);
	pSymMap->GetMap()[SymbolMapKey("<chi2_red>")] = new SymbolReal(dChi2red);
	return pSymMap;
}
// --------------------------------------------------------------------------------



// --------------------------------------------------------------------------------
// interpolation

namespace ublas = boost::numeric::ublas;

enum FktParam
{
	FKT_SPLINE,
	FKT_BEZIER
};


// bezier(x, y, 128)
// spline(x, y, 128, degree)
static Symbol* _fkt_param(FktParam whichfkt, const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(vecSyms.size() < 2 || !is_vec(vecSyms[0]) || !is_vec(vecSyms[1]))
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "Function needs x and y vector arguments." << std::endl;
		throw tl::Err(ostrErr.str(),0);
	}

	unsigned int ideg = 3;
	unsigned int N = 128;

	if(vecSyms.size() > 2)
		N = vecSyms[2]->GetValInt();
	if(vecSyms.size() > 3)
		ideg = vecSyms[3]->GetValInt();

	std::vector<t_real> vecX = sym_to_vec<t_stdvec>(vecSyms[0]);
	std::vector<t_real> vecY = sym_to_vec<t_stdvec>(vecSyms[1]);

	unsigned int iSize = std::min(vecX.size(), vecY.size());

	tl::FunctionModel_param_gen<ublas::vector<t_real>>* pfkt = 0;
	if(whichfkt == FKT_SPLINE)
		pfkt = new tl::BSpline<t_real>(iSize, vecX.data(), vecY.data(), ideg);
	else if(whichfkt == FKT_BEZIER)
		pfkt = new tl::Bezier<t_real>(iSize, vecX.data(), vecY.data());
	else
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "Unknown parametric function selected." << std::endl;
		throw tl::Err(ostrErr.str(),0);
	}

	SymbolArray* pArrX = new SymbolArray();
	SymbolArray* pArrY = new SymbolArray();
	pArrX->GetArr().reserve(N);
	pArrY->GetArr().reserve(N);

	for(unsigned int i=0; i<N; ++i)
	{
		t_real t = t_real(i)/t_real(N-1);
		ublas::vector<t_real> vecSpl = (*pfkt)(t);

		if(vecSpl.size() < 2)
			continue;

		pArrX->GetArr().push_back(new SymbolReal(vecSpl[0]));
		pArrY->GetArr().push_back(new SymbolReal(vecSpl[1]));
	}

	delete pfkt;

	SymbolArray *pArr = new SymbolArray();
	pArr->GetArr().push_back(pArrX);
	pArr->GetArr().push_back(pArrY);

	pArrX->UpdateIndices();
	pArrY->UpdateIndices();
	pArr->UpdateIndices();

	return pArr;
}

static Symbol* fkt_bezier(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	return _fkt_param(FKT_BEZIER, vecSyms, info, runinfo, pSymTab);
}

static Symbol* fkt_spline(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	return _fkt_param(FKT_SPLINE, vecSyms, info, runinfo, pSymTab);
}

static Symbol* fkt_find_peaks(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(vecSyms.size() < 2)
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "find_peaks needs x and y arrays." << std::endl;
		throw tl::Err(ostrErr.str(),0);
	}

	const unsigned int iOrder = 5;

	std::vector<t_real> vecX = sym_to_vec<t_stdvec>(vecSyms[0]);
	std::vector<t_real> vecY = sym_to_vec<t_stdvec>(vecSyms[1]);
	const unsigned int iLen = std::min(vecX.size(), vecY.size());

	std::vector<t_real> vecMaximaX, vecMaximaSize, vecMaximaWidth;
	tl::find_peaks<t_real>(iLen, vecX.data(), vecY.data(), iOrder,
						vecMaximaX, vecMaximaSize, vecMaximaWidth);

	SymbolArray* pArrX = (SymbolArray*)vec_to_sym<t_stdvec>(vecMaximaX);
	SymbolArray* pArrSizes = (SymbolArray*)vec_to_sym<t_stdvec>(vecMaximaSize);
	SymbolArray* pArrWidths = (SymbolArray*)vec_to_sym<t_stdvec>(vecMaximaWidth);

	SymbolArray *pArr = new SymbolArray();
	pArr->GetArr().push_back(pArrX);
	pArr->GetArr().push_back(pArrSizes);
	pArr->GetArr().push_back(pArrWidths);
	pArr->UpdateIndices();

	return pArr;
}
// --------------------------------------------------------------------------------


// ["a" : [val0, val1]]  =>  ["a" : val0]
static Symbol* fkt_map_vec_to_val(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(vecSyms.size()<1 || vecSyms[0]->GetType()!=SYMBOL_MAP)
	{
		tl::log_err(linenr(runinfo), "Need a map of vectors.");
		return 0;
	}


	// index parameter
	t_int iIdx = 0;
	if(vecSyms.size()>1)
	{
		iIdx = vecSyms[1]->GetValInt();
		if(iIdx < 0)
		{
			tl::log_warn(linenr(runinfo), "Ignoring negative index.");
			iIdx = 0;
		}
	}


	SymbolMap* pMapRet = new SymbolMap();

	for(const SymbolMap::t_map::value_type& pair : ((SymbolMap*)vecSyms[0])->GetMap())
	{
		const t_string& strKey = pair.first.strKey;
		const Symbol* pSym = pair.second;

		if(pSym->GetType() != SYMBOL_ARRAY)
		{
			// this also ignores the "<valid>" fit variable automatically

			//G_CERR << linenr(T_STR"Warning", info)
			//			<< "Ignoring non-vector variable " 
			//			<< "\"" << strKey << "\""
			//			<< " in map." << std::endl;
			continue;
		}

		const std::vector<Symbol*>& arr = ((SymbolArray*)pSym)->GetArr();
		if(iIdx >= int(arr.size()))
		{
			tl::log_warn(linenr(runinfo), "Ignoring invalid index.");
			continue;
		}

		const Symbol* pElem = arr[iIdx];

		pMapRet->GetMap().insert(SymbolMap::t_map::value_type(strKey, pElem->clone()));
	}

	pMapRet->UpdateIndices();
	return pMapRet;
}

extern void init_ext_fit_calls()
{
	t_mapFkts mapFkts =
	{
		t_mapFkts::value_type(T_STR"fit", fkt_fit),

		t_mapFkts::value_type(T_STR"map_vec_to_val", fkt_map_vec_to_val),

		t_mapFkts::value_type(T_STR"spline", fkt_spline),
		t_mapFkts::value_type(T_STR"bezier", fkt_bezier),

		t_mapFkts::value_type(T_STR"find_peaks", fkt_find_peaks),
	};

	add_ext_calls(mapFkts);
}
