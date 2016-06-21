/*
 * invoke gnuplot
 * @autor tweber
 * @date 24-dec-2013
 * @license GPLv2 or GPLv3
 */

#ifndef __GPL_PLOTTER_IMPL_H__
#define __GPL_PLOTTER_IMPL_H__

#include "../helper/flags.h"
#include "../helper/misc.h"
#include "../string/string.h"
#include "../log/log.h"

#include <cstdio>
#include <sstream>
#include <fstream>
#include <iomanip>

#include "gnuplot.h"


namespace tl {

namespace ios = boost::iostreams;

template<class t_real>
void GnuPlot_gen<t_real>::DeInit()
{
	m_postr.reset();
	m_psbuf.reset();
	m_pfds.reset();

	if(m_pipe) ::pclose(m_pipe);
	m_pipe = 0;

	m_iStartCounter = 0;
}

template<class t_real>
void GnuPlot_gen<t_real>::Init()
{
	if(IsReady()) return;
	DeInit();

	m_pipe = (FILE*)::/*my_*/popen("gnuplot -p 2>/dev/null 1>/dev/null", "w");
	if(!m_pipe)
	{
		log_err("Could not load gnuplot.");
		return;
	}

	m_pfds.reset(new ios::file_descriptor_sink(fileno(m_pipe), ios::close_handle /*ios::never_close_handle*/));
	m_psbuf.reset(new ios::stream_buffer<ios::file_descriptor_sink>(*m_pfds));
	m_postr.reset(new std::ostream(m_psbuf.get()));

	(*m_postr) << "set grid\n";
	(*m_postr) << "set nokey\n";
	//(*m_postr) << "set noborder\n";
	(*m_postr) << "set size 1,1\n";
	(*m_postr) << "set palette rgbformulae 33,13,10\n";
}

template<class t_real>
void GnuPlot_gen<t_real>::SetTerminal(int iWnd, const char* pcBackend)
{
	if(m_bTermLocked) return;

	(*m_postr) << "set output\n";
	(*m_postr) << "set obj 1 rectangle behind fillcolor rgbcolor \"white\" from screen 0,0 to screen 1,1\n";

	(*m_postr) << "set term " << pcBackend <<  " " << iWnd << " "
		<< "size 640,480 "
		<< "enhanced "
		<< "font \"NimbusSanL-Regu,12\" "
//		<< "title \"" << "Plot " << (iWnd+1) << "\" " 
		<< "persist "
		<< "dashed "
		<<  "\n";
}

template<class t_real>
void GnuPlot_gen<t_real>::SetFileTerminal(const char* pcFile)
{
	if(m_bTermLocked) return;

	std::string strFile = pcFile;
	std::string strExt = get_fileext(strFile);

	if(str_is_equal(strExt, std::string("pdf"), 0))
	{
		(*m_postr) << "set term pdf enhanced color font \"NimbusSanL-Regu,16\" \n";
	}
	else if(str_is_equal(strExt, std::string("ps"), 0))
	{
		(*m_postr) << "set term postscript eps enhanced color font \"NimbusSanL-Regu,16\" \n";
	}
	else
	{
		log_err("Unknown file extension \"", strExt, "\" for output terminal.");
		return;
	}

	(*m_postr) << "set output \"" << strFile << "\"\n";
}

template<class t_real>
void GnuPlot_gen<t_real>::RefreshVars()
{
	if(m_bHasLegend)
		(*m_postr) << "set key on " << m_strLegendPlacement
			<< " box 1 " << m_strLegendOpts << "\n";
	else
		(*m_postr) << "set nokey\n";
}

template<class t_real>
void GnuPlot_gen<t_real>::SimplePlot(const std::vector<t_real>& vecX, const std::vector<t_real>& vecY,
	const std::vector<t_real>& vecYErr, const std::vector<t_real>& vecXErr,
	LineStyle style)
{
	if(!IsReady()) return;
	(*m_postr) << "plot \"-\" ";

	switch(style)
	{
		case STYLE_LINES_SOLID:
			(*m_postr) << "with lines linetype 1 linewidth 1";
			break;
		case STYLE_LINES_DASHED:
			(*m_postr) << "with lines linetype 2 linewidth 1";
			break;
		default:
		case STYLE_POINTS:
			break;
	}
	(*m_postr) << "\n";

	std::string strTable = BuildTable(vecX, vecY, vecYErr, vecXErr);
	(*m_postr) << strTable;
	m_postr->flush();
}


template<class t_real>
void GnuPlot_gen<t_real>::SimplePlot2d(const std::vector<std::vector<t_real> >& vec,
	t_real dMinX, t_real dMaxX, t_real dMinY, t_real dMaxY)
{
	if(!IsReady()) return;

	std::vector<unsigned int> vecSizes;
	vecSizes.reserve(vec.size());

	for(unsigned int iY=0; iY<vec.size(); ++iY)
		vecSizes.push_back(vec[iY].size());

	std::vector<unsigned int>::iterator iterMin =
		std::min_element(vecSizes.begin(), vecSizes.end());
	unsigned int iXCntMin = *iterMin;


	unsigned int iYDim = vec.size();
	unsigned int iXDim = iXCntMin;

	// invalid values select image dimensions
	if(dMinX > dMaxX)
	{
		dMinX = 0.;
		dMaxX = iXDim-1;
	}
	if(dMinY > dMaxY)
	{
		dMinY = 0.;
		dMaxY = iYDim-1;
	}

	// ----------------------------------------
	// ranges
	(*m_postr) << "set tics out scale 0.8\n";

	t_real dRangeMinX = tic_trafo<t_real>(iXDim, dMinX, dMaxX, 0, -0.5);
	t_real dRangeMaxX = tic_trafo<t_real>(iXDim, dMinX, dMaxX, 0, t_real(iXDim)-0.5);
	t_real dRangeMinY = tic_trafo<t_real>(iYDim, dMinY, dMaxY, 0, -0.5);
	t_real dRangeMaxY = tic_trafo<t_real>(iYDim, dMinY, dMaxY, 0, t_real(iYDim)-0.5);

	(*m_postr) << "set xrange [" << dRangeMinX << ":" << dRangeMaxX << "]\n";
	(*m_postr) << "set yrange [" << dRangeMinY << ":" << dRangeMaxY << "]\n";
	// ----------------------------------------

	// ----------------------------------------
	// tics
	std::ostringstream ostrTicsX, ostrTicsY;
	ostrTicsX << "(" << dMinX << " + " << "($1)/" << iXDim
			<< " * (" << dMaxX << "-" << dMinX << "))";
	ostrTicsY << "(" << dMinY << " + " << "($2)/" << iYDim
			<< " * (" << dMaxY << "-" << dMinY << "))";

	std::string strTics = "using " + ostrTicsX.str() + ":" + ostrTicsY.str() + ":3";
	// ----------------------------------------

	(*m_postr) << "plot \"-\" " << strTics << " matrix with image\n";


	for(unsigned int iY=0; iY<vec.size(); ++iY)
	{
		for(unsigned int iX=0; iX<iXCntMin; ++iX)
			(*m_postr) << vec[iY][iX] << " ";
		(*m_postr) << "\n";
	}

	(*m_postr) << "end\nend\n";


	m_postr->flush();
}

template<class t_real>
void GnuPlot_gen<t_real>::AddLine(const PlotObj_gen<t_real>& obj)
{
	m_vecObjs.push_back(obj);
}

template<class t_real>
void GnuPlot_gen<t_real>::StartPlot()
{
	if(!IsReady()) return;

	if(m_iStartCounter == 0)
		m_vecObjs.clear();

	++m_iStartCounter;
}

template<class t_real>
void GnuPlot_gen<t_real>::SetCmdFileOutput(const char* pcFile)
{
	m_strCmdFileOutput = pcFile;
}

template<class t_real>
void GnuPlot_gen<t_real>::FinishPlot()
{
	if(!IsReady()) return;

	if(--m_iStartCounter == 0)
	{
		std::string strCmd = BuildCmd();
		RefreshVars();

		if(m_strCmdFileOutput != "")
		{
			std::ofstream ofCmd(m_strCmdFileOutput);
			if(!!ofCmd)
			{
				ofCmd << "#!/usr/bin/gnuplot -p\n\n";
				ofCmd << strCmd << "\n";
				ofCmd.close();
				m_strCmdFileOutput = "";
			}
		}

		//std::cout << "Plot cmd: " << strCmd << std::endl;
		(*m_postr) << strCmd;
		//(*m_postr) << "replot\n";
		m_postr->flush();
		m_vecObjs.clear();
	}
}

template<class t_real>
std::string GnuPlot_gen<t_real>::BuildCmd()
{
	m_bHasLegend = 0;

	std::ostringstream ostr;
	ostr << "plot ";

	for(const PlotObj_gen<t_real>& obj : m_vecObjs)
	{
		const bool bConnectLines = (obj.linestyle != STYLE_POINTS);
		const bool bHasXErr = (obj.vecErrX.size() != 0);
		const bool bHasYErr = (obj.vecErrY.size() != 0);
		t_real dSize = boost::get_optional_value_or(obj.odSize, bConnectLines ? 1. : 1.25);

		std::ostringstream ostrTmp;
		std::string strPointStyle;

		if(bHasXErr && bHasYErr)
			ostrTmp << "with xyerrorbars";
		else if(bHasXErr && !bHasYErr)
			ostrTmp << "with xerrorbars";
		else if(!bHasXErr && bHasYErr)
			ostrTmp << "with yerrorbars";
		else if(!bHasXErr && !bHasYErr)
			ostrTmp << "with points";

		ostrTmp << " pointtype 7 pointsize " << dSize;
		strPointStyle = ostrTmp.str();


		ostr << "\"-\" ";
		switch(obj.linestyle)
		{
			case STYLE_LINES_SOLID:
				ostr << "with lines linetype 1 linewidth " << dSize;
				break;
			case STYLE_LINES_DASHED:
				ostr << "with lines linetype 2 linewidth " << dSize;
				break;
			default:
				log_warn("Unknown line style.");
			case STYLE_POINTS:
				ostr << strPointStyle;
				break;
		}

		if(obj.oiColor)
		{
			unsigned int iColor = *obj.oiColor;
			char chFill = ostr.fill();
			ostr << std::setfill('0');
			ostr << "linecolor rgb \"#" << std::hex
				<< std::setw(2) << ((iColor & 0xff0000) >> 16)
				<< std::setw(2) << ((iColor & 0x00ff00) >> 8)
				<< std::setw(2) << ((iColor & 0x0000ff))
				<< "\" " << std::dec;
			ostr << std::setfill(chFill);
		}

		if(obj.strLegend != "")
		{
			m_bHasLegend = 1;
			ostr << " title \"" << obj.strLegend << "\" ";
		}

		if(&obj != &(*m_vecObjs.rbegin()))
			ostr << ", ";
	}
	ostr << "\n";

	for(const PlotObj_gen<t_real>& obj : m_vecObjs)
	{
		std::string strTab = BuildTable(obj.vecX, obj.vecY, obj.vecErrY, obj.vecErrX);
		ostr << strTab;
	}

	return ostr.str();
}

template<class t_real>
std::string GnuPlot_gen<t_real>::BuildTable(const std::vector<t_real>& vecX, const std::vector<t_real>& vecY,
	const std::vector<t_real>& vecYErr, const std::vector<t_real>& vecXErr)
{
	std::ostringstream ostr;

	const unsigned int iSize = std::min(vecX.size(), vecY.size());
	const bool bHasXErr = (vecXErr.size() != 0);
	const bool bHasYErr = (vecYErr.size() != 0);

	for(unsigned int iDat=0; iDat<iSize; ++iDat)
	{
		ostr << vecX[iDat] << " " << vecY[iDat];

		if(bHasXErr) ostr << " " << vecXErr[iDat];
		if(bHasYErr) ostr << " " << vecYErr[iDat];

		ostr << "\n";
	}
	ostr << "end\n";

	return ostr.str();
}

template<class t_real>
void GnuPlot_gen<t_real>::SetXLabel(const char* pcLab)
{
	if(!IsReady()) return;
	(*m_postr) << "set xlabel \"" << pcLab << "\"\n";
	//m_postr->flush();
}

template<class t_real>
void GnuPlot_gen<t_real>::SetYLabel(const char* pcLab)
{
	if(!IsReady()) return;
	(*m_postr) << "set ylabel \"" << pcLab << "\"\n";
	//m_postr->flush();
}

template<class t_real>
void GnuPlot_gen<t_real>::SetTitle(const char* pcTitle)
{
	if(!IsReady()) return;
	(*m_postr) << "set title \"" << pcTitle << "\"\n";
	//m_postr->flush();
}

template<class t_real>
void GnuPlot_gen<t_real>::SetXRange(t_real dMin, t_real dMax)
{
	if(!IsReady()) return;
	//std::cout << "xmin: "  << dMin << ", xmax: " << dMax << std::endl;
	(*m_postr) << "set xrange [" << dMin << ":" << dMax << "]\n";
	m_postr->flush();
}

template<class t_real>
void GnuPlot_gen<t_real>::SetYRange(t_real dMin, t_real dMax)
{
	if(!IsReady()) return;
	(*m_postr) << "set yrange [" << dMin << ":" << dMax << "]\n";
	m_postr->flush();
}

template<class t_real>
void GnuPlot_gen<t_real>::SetLogX(t_real tBase)
{
	if(!IsReady()) return;
	if(tBase >= 0.)
		(*m_postr) << "set logscale x " << tBase << "\n";
	else
		(*m_postr) << "unset logscale x\n";
	m_postr->flush();
}

template<class t_real>
void GnuPlot_gen<t_real>::SetLogY(t_real tBase)
{
	if(!IsReady()) return;
	if(tBase >= 0.)
		(*m_postr) << "set logscale y " << tBase << "\n";
	else
		(*m_postr) << "unset logscale y\n";
	m_postr->flush();
}

template<class t_real>
void GnuPlot_gen<t_real>::SetGrid(bool bOn)
{
	if(!IsReady()) return;

	if(bOn)
		(*m_postr) << "set grid\n";
	else
		(*m_postr) << "unset grid\n";

	m_postr->flush();
}

template<class t_real>
void GnuPlot_gen<t_real>::AddArrow(t_real dX0, t_real dY0, t_real dX1, t_real dY1, bool bHead)
{
	if(!IsReady()) return;

	(*m_postr) << "set arrow from " << dX0 << "," << dY0
			<< " to " << dX1 << "," << dY1;
	if(!bHead)
		(*m_postr) << " nohead";

	(*m_postr) << "\n";
	m_postr->flush();
}

template<class t_real>
void GnuPlot_gen<t_real>::SetColorBarRange(t_real dMin, t_real dMax, bool bCyclic)
{
	if(!IsReady()) return;

	(*m_postr) << "set cbrange [" << dMin << ":" << dMax << "]\n";

	if(bCyclic)
		(*m_postr) << "set palette defined (0 \"blue\", 0.25 \"cyan\", 0.5 \"yellow\", 0.75 \"red\", 1 \"blue\")\n";
	else
		(*m_postr) << "set palette defined (0 \"blue\", 0.3333 \"cyan\", 0.6666 \"yellow\", 1 \"red\")\n";

	m_postr->flush();
}

template<class t_real>
bool GnuPlot_gen<t_real>::IsReady() const { return m_postr!=0; }

template<class t_real>
std::ostream& GnuPlot_gen<t_real>::GetStream() { return *m_postr; }

}

#endif
