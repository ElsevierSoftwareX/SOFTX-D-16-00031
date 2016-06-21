/*
 * external plot functions
 * @author tweber
 * @date dec 2013
 * @license GPLv2 or GPLv3
 */

#include "lang/types.h"
#include "calls_plot.h"
#include "lang/calls.h"
#include "gfx/gnuplot.h"
#include "string/string.h"


// --------------------------------------------------------------------------------
// plotting

//#define DEFAULT_TERM T_STR"qt";
#define DEFAULT_TERM T_STR"x11";
static tl::GnuPlot_gen<t_real> g_plot;

static inline bool is_array_of_arrays(const Symbol* pSym)
{
	if(!pSym) return 0;
	if(pSym->GetType()!=SYMBOL_ARRAY) return 0;

	SymbolArray* pSymArr = (SymbolArray*)pSym;
	if(pSymArr->GetArr().size()==0) return 0;

	Symbol* pSymInArr = pSymArr->GetArr()[0];
	return (pSymInArr->GetType()==SYMBOL_ARRAY);
}

static inline bool is_array_of_array_of_arrays(const Symbol* pSym)
{
	if(!pSym) return 0;
	if(pSym->GetType() != SYMBOL_ARRAY)
		return 0;

	if(((SymbolArray*)pSym)->GetArr().size() == 0)
		return 0;

	return is_array_of_arrays(((SymbolArray*)pSym)->GetArr()[0]);
}

struct XYLimits
{
	bool bHasX, bHasY, bHasCB, bCBCyclic;
	t_real dMinX, dMaxX;
	t_real dMinY, dMaxY;
	t_real dMinCB, dMaxCB;

	XYLimits() : bHasX(0), bHasY(0), bHasCB(0), bCBCyclic(0)
	{}
};

static XYLimits get_plot_limits(const SymbolMap* pParamMap)
{
	XYLimits lim;

	bool bHasVal=0;
	t_string strVal = pParamMap->GetStringVal(T_STR"xylimits", &bHasVal);
	if(bHasVal)
	{
		t_istringstream istr(strVal);
		istr >> lim.dMinX >> lim.dMaxX >> lim.dMinY >> lim.dMaxY;

		lim.bHasX = 1;
		lim.bHasY = 1;
	}
	else
	{
		t_string strX = pParamMap->GetStringVal(T_STR"xlimits", &lim.bHasX);
		t_string strY = pParamMap->GetStringVal(T_STR"ylimits", &lim.bHasY);

		if(lim.bHasX)
		{
			t_istringstream istrX(strX);
			istrX >> lim.dMinX >> lim.dMaxX;
		}

		if(lim.bHasY)
		{
			t_istringstream istrY(strY);
			istrY >> lim.dMinY >> lim.dMaxY;
		}
	}


	t_string strValCB = pParamMap->GetStringVal(T_STR"cblimits", &lim.bHasCB);
	t_istringstream istrCB(strValCB);
	istrCB >> lim.dMinCB >> lim.dMaxCB;

	bool bHasCyc = 0;
	t_string strValCBCyc = pParamMap->GetStringVal(T_STR"cbcyclic", &bHasCyc);
	if(bHasCyc)
	{
		t_istringstream istrCyc(strValCBCyc);
		istrCyc >> lim.bCBCyclic;
	}

	return lim;
}


static void set_plot_params(tl::GnuPlot_gen<t_real>& plot, const SymbolMap* pParamMap, 
	tl::PlotObj_gen<t_real>* pCurPlotObj=0, XYLimits* pLimits=0)
{
	bool bHasVal = 0;
	t_string strTitle = pParamMap->GetStringVal(T_STR"title", &bHasVal);
	if(bHasVal) plot.SetTitle(WSTR_TO_STR(strTitle).c_str());
	t_string strXLab = pParamMap->GetStringVal(T_STR"xlabel", &bHasVal);
	if(bHasVal) plot.SetXLabel(WSTR_TO_STR(strXLab).c_str());
	t_string strYLab = pParamMap->GetStringVal(T_STR"ylabel", &bHasVal);
	if(bHasVal) plot.SetYLabel(WSTR_TO_STR(strYLab).c_str());

	t_real dLogX = pParamMap->GetRealVal(T_STR"xlog", &bHasVal);
	if(bHasVal) plot.SetLogX(dLogX);
	t_real dLogY = pParamMap->GetRealVal(T_STR"ylog", &bHasVal);
	if(bHasVal) plot.SetLogY(dLogY);
	//tl::log_debug("Log bases: ", dLogX, ", ", dLogY);

	if(pCurPlotObj)
	{
		t_string strStyle = pParamMap->GetStringVal(T_STR"style", &bHasVal);
		if(bHasVal)
		{
			if(strStyle == T_STR"line" || strStyle == T_STR"lines")
				pCurPlotObj->linestyle = tl::STYLE_LINES_SOLID;
			else if(strStyle == T_STR"line_dashed" || strStyle == T_STR"lines_dashed")
				pCurPlotObj->linestyle = tl::STYLE_LINES_DASHED;
			else if(strStyle == T_STR"point" || strStyle == T_STR"points")
				pCurPlotObj->linestyle = tl::STYLE_POINTS;
		}

		t_string strLegend = pParamMap->GetStringVal(T_STR"legend", &bHasVal);
		if(bHasVal) pCurPlotObj->strLegend = WSTR_TO_STR(strLegend);

		bool bHasSize=0, bHasColor=0;
		t_real dSize = pParamMap->GetRealVal(T_STR"size", &bHasSize);
		if(bHasSize) pCurPlotObj->odSize = dSize;

		unsigned int iColor = pParamMap->GetIntVal(T_STR"color", &bHasColor);
		if(bHasColor) pCurPlotObj->oiColor = iColor;
	}

	XYLimits lim = get_plot_limits(pParamMap);
	if(lim.bHasX) plot.SetXRange(lim.dMinX, lim.dMaxX);
	if(lim.bHasY) plot.SetYRange(lim.dMinY, lim.dMaxY);


	// for 2D plot
	if(pLimits)
	{
		*pLimits = lim;
		if(pLimits->bHasCB)
			plot.SetColorBarRange(pLimits->dMinCB, pLimits->dMaxCB, pLimits->bCBCyclic);
	}


	// terminal
	t_string strTerm = DEFAULT_TERM;
	t_string strUserTerm = pParamMap->GetStringVal(T_STR"term", &bHasVal);
	if(bHasVal) strTerm = strUserTerm;


	int iGrid = pParamMap->GetIntVal(T_STR"grid", &bHasVal);
	if(bHasVal) plot.SetGrid(iGrid!=0);


	t_string strCmdFile = pParamMap->GetStringVal(T_STR"cmdfile", &bHasVal);


	t_string strArr = pParamMap->GetStringVal(T_STR"arrow", &bHasVal);
	if(bHasVal)
	{
		std::vector<t_string> vecArrToks;
		tl::get_tokens<t_string, t_string>(strArr, ",", vecArrToks);

		for(const t_string& strArrTok : vecArrToks)
		{
			std::istringstream istrArr(strArrTok);
			t_real dX0, dY0, dX1, dY1;
			bool bHead = 0;
			istrArr >> dX0 >> dY0 >> dX1 >> dY1 >> bHead;

			plot.AddArrow(dX0, dY0, dX1, dY1, bHead);
		}
	}


	// legend options
	t_string strLegendOpts = pParamMap->GetStringVal(T_STR"legend_opts", &bHasVal);
	if(bHasVal) plot.SetLegendOpts(WSTR_TO_STR(strLegendOpts));
	t_string strLegendPlace = pParamMap->GetStringVal(T_STR"legend_place", &bHasVal);
	if(bHasVal) plot.SetLegendPlace(WSTR_TO_STR(strLegendPlace));


	int iPlotWnd = 0;
	int iUserPlotWnd = atoi(WSTR_TO_STR(pParamMap->GetStringVal(T_STR"window", &bHasVal)).c_str());
	if(bHasVal) iPlotWnd = iUserPlotWnd;

	plot.SetTerminal(iPlotWnd, WSTR_TO_STR(strTerm).c_str());
	plot.SetCmdFileOutput(WSTR_TO_STR(strCmdFile).c_str());
}

static Symbol* fkt_fileplot(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab);

static Symbol* fkt_plot(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	g_plot.Init();
	unsigned int iNumSyms = vecSyms.size();
	/*const*/ SymbolMap *pPlotParams = 0;

	// plot([[x, y, yerr, xerr, mapParams], ...]);
	if(iNumSyms==1 && is_array_of_array_of_arrays(vecSyms[0]))
		return fkt_plot(((SymbolArray*)vecSyms[0])->GetArr(), info, runinfo, pSymTab);
	// plot([x, y, yerr, xerr, mapParams], [x2, y2, yerr2, xerr2, mapParams2], ...)
	else if(iNumSyms>=1 && is_array_of_arrays(vecSyms[0]))
	{
		g_plot.StartPlot();
		for(Symbol *pArr : vecSyms)
		{
			SymbolType symType = pArr->GetType();

			// ignore non-array arguments
			if(symType == SYMBOL_ARRAY)
				fkt_plot(((SymbolArray*)pArr)->GetArr(), info, runinfo, pSymTab);
			else if(symType == SYMBOL_MAP)
			{
				pPlotParams = (SymbolMap*)pArr;
				set_plot_params(g_plot, pPlotParams);
			}
		}
		g_plot.FinishPlot();
	}
	// plot(x, y, yerr, xerr, mapParams)
	else if(iNumSyms >= 2 && vecSyms[0]->GetType()==SYMBOL_ARRAY && vecSyms[1]->GetType()==SYMBOL_ARRAY)
	{
		std::vector<t_real> vecX = ((SymbolArray*)vecSyms[0])->ToDoubleArray();
		std::vector<t_real> vecY = ((SymbolArray*)vecSyms[1])->ToDoubleArray();
		std::vector<t_real> vecYErr, vecXErr;

		if(iNumSyms >= 3 && vecSyms[2]->GetType()==SYMBOL_ARRAY)
			vecYErr = ((SymbolArray*)vecSyms[2])->ToDoubleArray();
		if(iNumSyms >= 4 && vecSyms[3]->GetType()==SYMBOL_ARRAY)
			vecXErr = ((SymbolArray*)vecSyms[3])->ToDoubleArray();


		tl::PlotObj_gen<t_real> obj;
		obj.vecX = vecX;
		obj.vecY = vecY;
		obj.vecErrX = vecXErr;
		obj.vecErrY = vecYErr;

		// parameter map given as last argument
		if(vecSyms[iNumSyms-1]->GetType()==SYMBOL_MAP)
		{
			pPlotParams = (SymbolMap*)vecSyms[iNumSyms-1];
			set_plot_params(g_plot, pPlotParams, &obj);
		}

		g_plot.StartPlot();
		g_plot.AddLine(obj);
		g_plot.FinishPlot();
	}
	else
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "Invalid call to plot." << std::endl;
		throw tl::Err(ostrErr.str(), 0);
	}

	if(pPlotParams)
	{
		// also plot to file
		if(pPlotParams->GetMap().find(SymbolMapKey("outfile")) != pPlotParams->GetMap().end()
			&& pPlotParams->GetMap().find(SymbolMapKey("<outfile_written>")) == pPlotParams->GetMap().end())
		{
			// hack to prevent infinite recursion
			pPlotParams->GetMap().insert(
				SymbolMap::t_map::value_type(SymbolMapKey("<outfile_written>"), new SymbolInt(1)));

			Symbol* pSymFileName = pPlotParams->GetMap()[SymbolMapKey("outfile")];
			std::vector<Symbol*> vecNewSyms;
			vecNewSyms.reserve(vecSyms.size()+1);
			vecNewSyms.push_back(pSymFileName);

			for(unsigned int iSym=0; iSym<vecSyms.size(); ++iSym)
				vecNewSyms.push_back(vecSyms[iSym]);

			fkt_fileplot(vecNewSyms, info, runinfo, pSymTab);
		}
	}
	return 0;
}

static Symbol* fkt_plot2d(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	g_plot.Init();

	unsigned int iNumSyms = vecSyms.size();

	// e.g. plot2d([[1,2],[3,4]], params)
	if(iNumSyms >= 1 && is_array_of_arrays(vecSyms[0]))
	{
		SymbolArray* pArr = (SymbolArray*)vecSyms[0];
		SymbolMap* pMapParam = 0;

		// has parameter map
		if(iNumSyms == 2 && vecSyms[1]->GetType()==SYMBOL_MAP)
			pMapParam = (SymbolMap*)vecSyms[1];

		std::vector<std::vector<t_real> > vecXY;
		vecXY.reserve(pArr->GetArr().size());

		for(const Symbol* _pX : pArr->GetArr())
		{
			SymbolArray* pX = (SymbolArray*)_pX;
			std::vector<t_real> vecX = pX->ToDoubleArray();

			vecXY.push_back(vecX);
		}

		t_real dRMinX=1., dRMaxX=-1., dRMinY=1., dRMaxY=-1.;
		if(pMapParam)
		{
			XYLimits lim;
			set_plot_params(g_plot, (SymbolMap*)pMapParam, 0, &lim);
        	        if(lim.bHasX) { dRMinX = lim.dMinX; dRMaxX = lim.dMaxX; }
	                if(lim.bHasY) { dRMinY = lim.dMinY; dRMaxY = lim.dMaxY; }
		}

		g_plot.SimplePlot2d(vecXY, dRMinX, dRMaxX, dRMinY, dRMaxY);
	}
	else
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "Invalid call to plot2d." << std::endl;
		throw tl::Err(ostrErr.str(),0);
	}

	return 0;
}


static Symbol* _fkt_fileplot(const std::vector<Symbol*>& vecSyms, 
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab,
	Symbol* (*pPltFkt)(const::std::vector<Symbol*>&, ParseInfo&, RuntimeInfo&, SymbolTable*))
{
	g_plot.Init();
	if(vecSyms.size() < 1 || vecSyms[0]->GetType()!=SYMBOL_STRING)
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo)
			<< "First argument to fileplot has to be the file name." << std::endl;
		throw tl::Err(ostrErr.str(),0);
	}

	const t_string& strFile = ((SymbolString*)vecSyms[0])->GetVal();
	g_plot.SetFileTerminal(WSTR_TO_STR(strFile).c_str());
	g_plot.LockTerminal();

	std::vector<Symbol*> vecPlot;
	vecPlot.reserve(vecSyms.size()-1);
	for(unsigned int iSym=1; iSym<vecSyms.size(); ++iSym)
		vecPlot.push_back(vecSyms[iSym]);

	Symbol* pSymRet = pPltFkt(vecPlot, info, runinfo, pSymTab);

	g_plot.UnlockTerminal();
	return pSymRet;
}

static Symbol* fkt_fileplot(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	return _fkt_fileplot(vecSyms, info, runinfo, pSymTab, fkt_plot);
}

static Symbol* fkt_fileplot2d(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	return _fkt_fileplot(vecSyms, info, runinfo, pSymTab, fkt_plot2d);
}


// --------------------------------------------------------------------------------



extern void init_ext_plot_calls()
{
	t_mapFkts mapFkts =
	{
		t_mapFkts::value_type(T_STR"plot", fkt_plot),
		t_mapFkts::value_type(T_STR"plot2d", fkt_plot2d),

		t_mapFkts::value_type(T_STR"fileplot", fkt_fileplot),
		t_mapFkts::value_type(T_STR"fileplot2d", fkt_fileplot2d),
	};

	add_ext_calls(mapFkts);
}
