/*
 * invoke gnuplot
 * @autor tweber
 * @date 24-dec-2013
 * @license GPLv2 or GPLv3
 */

#include "gnuplot.h"
#include "gnuplot_impl.h"

namespace tl
{
	template struct PlotObj_gen<double>;
	//template struct PlotObj_gen<float>;

	template class GnuPlot_gen<double>;
	//template class GnuPlot_gen<float>;
}
