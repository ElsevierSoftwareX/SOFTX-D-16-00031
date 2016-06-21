/*
 * invoke gnuplot
 * @autor tweber
 * @date 24-dec-2013
 * @license GPLv2 or GPLv3
 */

#ifndef __M_PLOTTER_GPL__
#define __M_PLOTTER_GPL__

#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/optional.hpp>
#include <iostream>
#include <string>
#include <vector>
#include <memory>

namespace tl {
enum LineStyle
{
	STYLE_POINTS,
	STYLE_LINES_SOLID,
	STYLE_LINES_DASHED
};

template<class t_real = double>
struct PlotObj_gen
{
	std::vector<t_real> vecX, vecY;
	std::vector<t_real> vecErrX, vecErrY;

	std::string strLegend;

	LineStyle linestyle = STYLE_POINTS;

	boost::optional<t_real> odSize;
	boost::optional<unsigned int> oiColor;

	PlotObj_gen() = default;
	~PlotObj_gen() = default;
};


template<class t_real = double>
class GnuPlot_gen
{
protected:
	FILE *m_pipe = nullptr;
	std::unique_ptr<boost::iostreams::file_descriptor_sink> m_pfds;
	std::unique_ptr<boost::iostreams::stream_buffer<boost::iostreams::file_descriptor_sink>> m_psbuf;
	std::unique_ptr<std::ostream> m_postr;

	std::vector<PlotObj_gen<t_real>> m_vecObjs;
	// has to be 0 to show plot
	int m_iStartCounter = 0;

	bool m_bTermLocked = false;
	bool m_bHasLegend = false;

	std::string m_strLegendOpts;
	std::string m_strLegendPlacement = std::string("default");

	std::string m_strCmdFileOutput;

protected:
	std::string BuildCmd();
	std::string BuildTable(const std::vector<t_real>& vecX, const std::vector<t_real>& vecY,
		const std::vector<t_real>& vecYErr, const std::vector<t_real>& vecXErr);
	void RefreshVars();

public:
	GnuPlot_gen() = default;
	virtual ~GnuPlot_gen() { DeInit(); }

	void Init();
	void DeInit();

	bool IsReady() const;
	std::ostream& GetStream();

	void SetTerminal(int iWnd=0, const char* pcBackend="x11");
	void SetFileTerminal(const char* pcFile);
	void SetCmdFileOutput(const char* pcFile);

	void StartPlot();
	void AddLine(const PlotObj_gen<t_real>& obj);
	void FinishPlot();

	void SimplePlot(const std::vector<t_real>& vecX, const std::vector<t_real>& vecY,
		const std::vector<t_real>& vecYErr, const std::vector<t_real>& vecXErr,
		LineStyle style=STYLE_POINTS);
	void SimplePlot2d(const std::vector<std::vector<t_real> >& vec,
		t_real dMinX=0., t_real dMaxX=-1., t_real dMinY=0., t_real dMaxY=-1.);

	void SetXLabel(const char* pcLab);
	void SetYLabel(const char* pcLab);
	void SetTitle(const char* pcTitle);
	void SetGrid(bool bOn);
	void AddArrow(t_real dX0, t_real dY0, t_real dX1, t_real dY1, bool bHead=1);

	void SetXRange(t_real dMin, t_real dMax);
	void SetYRange(t_real dMin, t_real dMax);
	void SetLogX(t_real tBase);
	void SetLogY(t_real tBase);

	void SetColorBarRange(t_real dMin, t_real dMax, bool bCyclic=0);

	void LockTerminal() { m_bTermLocked = 1; }
	void UnlockTerminal() { m_bTermLocked = 0; }

	void SetLegendOpts(const std::string& strOpts) { m_strLegendOpts = strOpts; }
	void SetLegendPlace(const std::string& strPlace) { m_strLegendPlacement = strPlace; }
};


	using PlotObj = PlotObj_gen<>;
	using GnuPlot = GnuPlot_gen<>;
}


#ifdef TLIBS_INC_HDR_IMPLS
        #include "gnuplot_impl.h"
#endif

#endif
