/*
 * Py helpers
 * @author tweber
 * @date 27-aug-2014, 17-feb-2015
 * @copyright GPLv2 or GPLv3
 */

#ifndef __PY_HELPERS_H__
#define __PY_HELPERS_H__

#include "../string/string.h"

namespace tl {

template<class t_str=std::string, class t_cont=std::vector<double>>
t_cont get_py_array(const t_str& str)
{
	typedef typename t_cont::value_type t_elems;
	t_cont vecArr;

	std::size_t iStart = str.find('[');
	std::size_t iEnd = str.find(']');

	// search for list instead
	if(iStart==t_str::npos || iEnd==t_str::npos)
	{
		iStart = str.find('(');
		iEnd = str.find(')');
	}

	// invalid array
	if(iStart==t_str::npos || iEnd==t_str::npos || iEnd<iStart)
		return vecArr;

	t_str strArr = str.substr(iStart+1, iEnd-iStart-1);
	tl::get_tokens<t_elems, t_str>(strArr, ",", vecArr);

	return vecArr;
}

template<class t_str=std::string>
t_str get_py_string(const t_str& str)
{
	std::size_t iStart = str.find_first_of("\'\"");
	std::size_t iEnd = str.find_last_of("\'\"");

	// invalid string
	if(iStart==t_str::npos || iEnd==t_str::npos || iEnd<iStart)
		return "";

	return str.substr(iStart+1, iEnd-iStart-1);
}

}

#endif
