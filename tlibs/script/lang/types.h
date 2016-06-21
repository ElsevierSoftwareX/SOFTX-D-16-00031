/*
 * Data types
 * @author tweber
 * @date 08-mar-2014
 * @license GPLv2 or GPLv3
 */

#ifndef __SCRIPT_TYPES_H__
#define __SCRIPT_TYPES_H__

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>
#include "helper/traits.h"


//typedef float t_real;
typedef double t_real;
typedef int t_int;
typedef std::complex<t_real> t_complex;


template<typename tyEnum>
struct EnumDirectHash
{
	std::size_t operator()(tyEnum ty) const
	{
		return static_cast<std::size_t>(ty);
	}
};

//#define USE_WCHAR


#ifdef USE_WCHAR
	typedef wchar_t t_char;
	typedef unsigned wchar_t t_uchar;
	typedef std::wstring t_string;

	typedef std::wostream t_ostream;
	typedef std::wistream t_istream;

	typedef std::wofstream t_ofstream;
	typedef std::wifstream t_ifstream;

	typedef std::wostringstream t_ostringstream;
	typedef std::wistringstream t_istringstream;

	#define G_CIN std::wcin
	#define G_COUT std::wcout
	#define G_CERR std::wcerr

	#define T_STR L

	#define STR_TO_WSTR str_to_wstr
	#define WSTR_TO_STR wstr_to_str
#else

	typedef char t_char;
	typedef unsigned char t_uchar;
	typedef std::string t_string;

	typedef std::ostream t_ostream;
	typedef std::istream t_istream;

	typedef std::ofstream t_ofstream;
	typedef std::ifstream t_ifstream;

	typedef std::ostringstream t_ostringstream;
	typedef std::istringstream t_istringstream;

	#define G_CIN std::cin
	#define G_COUT std::cout
	#define G_CERR std::cerr

	#define T_STR

        #define STR_TO_WSTR 
        #define WSTR_TO_STR 
#endif

#endif
