/*
 * Special chars
 * @author tweber
 * @date 09-mar-14
 * @license GPLv2 or GPLv3
 */

#include "spec_char.h"
#include <iostream>

namespace tl {

static t_mapSpecChars *g_pmapSpecChars = 0;

void init_spec_chars()
{
	deinit_spec_chars();
	g_pmapSpecChars = new t_mapSpecChars;

	(*g_pmapSpecChars)["AA"] = SpecChar("\xe2\x84\xab", L"\x212b");
	(*g_pmapSpecChars)["deg"] = SpecChar("\xc2\xb0", L"\x00b0");
	(*g_pmapSpecChars)["pm"] = SpecChar("\xc2\xb1", L"\x00b1");
	(*g_pmapSpecChars)["mp"] = SpecChar("\xe2\x88\x93", L"\x2213");
	(*g_pmapSpecChars)["cross"] = SpecChar("\xc3\x97", L"\x00d7");
	(*g_pmapSpecChars)["cdot"] = SpecChar("\xe2\x8b\x85", L"\x22c5");
	(*g_pmapSpecChars)["neq"] = SpecChar("\xe2\x89\xa0", L"\x2260");
	(*g_pmapSpecChars)["approx"] = SpecChar("\xe2\x89\x88", L"\x2248");

	(*g_pmapSpecChars)["hbar"] = SpecChar("\xe2\x84\x8f", L"\x210f");

	(*g_pmapSpecChars)["sub0"] = SpecChar("\xe2\x82\x80", L"\x2080");
	(*g_pmapSpecChars)["sub1"] = SpecChar("\xe2\x82\x81", L"\x2081");
	(*g_pmapSpecChars)["sub2"] = SpecChar("\xe2\x82\x82", L"\x2082");
	(*g_pmapSpecChars)["sub3"] = SpecChar("\xe2\x82\x83", L"\x2083");
	(*g_pmapSpecChars)["sub4"] = SpecChar("\xe2\x82\x84", L"\x2084");
	(*g_pmapSpecChars)["sub5"] = SpecChar("\xe2\x82\x85", L"\x2085");
	(*g_pmapSpecChars)["sub6"] = SpecChar("\xe2\x82\x86", L"\x2086");
	(*g_pmapSpecChars)["sub7"] = SpecChar("\xe2\x82\x87", L"\x2087");
	(*g_pmapSpecChars)["sub8"] = SpecChar("\xe2\x82\x88", L"\x2088");
	(*g_pmapSpecChars)["sub9"] = SpecChar("\xe2\x82\x89", L"\x2089");

	(*g_pmapSpecChars)["sup0"] = SpecChar("\xe2\x81\xb0", L"\x2070");
	(*g_pmapSpecChars)["sup1"] = SpecChar("\xc2\xb9", L"\x00b9");
	(*g_pmapSpecChars)["sup2"] = SpecChar("\xc2\xb2", L"\x00b2");
	(*g_pmapSpecChars)["sup3"] = SpecChar("\xc2\xb3", L"\x00b3");
	(*g_pmapSpecChars)["sup4"] = SpecChar("\xe2\x81\xb4", L"\x2074");
	(*g_pmapSpecChars)["sup5"] = SpecChar("\xe2\x81\xb5", L"\x2075");
	(*g_pmapSpecChars)["sup6"] = SpecChar("\xe2\x81\xb6", L"\x2076");
	(*g_pmapSpecChars)["sup7"] = SpecChar("\xe2\x81\xb7", L"\x2077");
	(*g_pmapSpecChars)["sup8"] = SpecChar("\xe2\x81\xb8", L"\x2078");
	(*g_pmapSpecChars)["sup9"] = SpecChar("\xe2\x81\xb9", L"\x2079");
	(*g_pmapSpecChars)["sup-"] = SpecChar("\xe2\x81\xbb", L"\x207b");
	(*g_pmapSpecChars)["sup+"] = SpecChar("\xe2\x81\xba", L"\x207a");

	(*g_pmapSpecChars)["alpha"] = SpecChar("\xce\xb1", L"\x03b1");
	(*g_pmapSpecChars)["beta"] = SpecChar("\xce\xb2", L"\x03b2");
	(*g_pmapSpecChars)["gamma"] = SpecChar("\xce\xb3", L"\x03b3");
	(*g_pmapSpecChars)["delta"] = SpecChar("\xce\xb4", L"\x03b4");
	(*g_pmapSpecChars)["epsilon"] = SpecChar("\xce\xb5", L"\x03b5");
	(*g_pmapSpecChars)["zeta"] = SpecChar("\xce\xb6", L"\x03b6");
	(*g_pmapSpecChars)["eta"] = SpecChar("\xce\xb7", L"\x03b7");
	(*g_pmapSpecChars)["theta"] = SpecChar("\xce\xb8", L"\x03b8");
	(*g_pmapSpecChars)["iota"] = SpecChar("\xce\xb9", L"\x03b9");
	(*g_pmapSpecChars)["kappa"] = SpecChar("\xce\xba", L"\x03ba");
	(*g_pmapSpecChars)["lambda"] = SpecChar("\xce\xbb", L"\x03bb");
	(*g_pmapSpecChars)["mu"] = SpecChar("\xce\xbc", L"\x03bc");
	(*g_pmapSpecChars)["nu"] = SpecChar("\xce\xbd", L"\x03bd");
	(*g_pmapSpecChars)["xi"] = SpecChar("\xce\xbe", L"\x03be");
	(*g_pmapSpecChars)["omicron"] = SpecChar("\xce\xbf", L"\x03bf");
	(*g_pmapSpecChars)["pi"] = SpecChar("\xcf\x80", L"\x03c0");
	(*g_pmapSpecChars)["rho"] = SpecChar("\xcf\x81", L"\x03c1");
	(*g_pmapSpecChars)["sigma"] = SpecChar("\xcf\x83", L"\x03c3");
	(*g_pmapSpecChars)["tau"] = SpecChar("\xcf\x84", L"\x03c4");
	(*g_pmapSpecChars)["upsilon"] = SpecChar("\xcf\x85", L"\x03c5");
	(*g_pmapSpecChars)["phi"] = SpecChar("\xcf\x86", L"\x03c6");
	(*g_pmapSpecChars)["chi"] = SpecChar("\xcf\x87", L"\x03c7");
	(*g_pmapSpecChars)["psi"] = SpecChar("\xcf\x88", L"\x03c8");
	(*g_pmapSpecChars)["omega"] = SpecChar("\xcf\x89", L"\x03c9");

	(*g_pmapSpecChars)["Alpha"] = SpecChar("\xce\x91", L"\x0391");
	(*g_pmapSpecChars)["Beta"] = SpecChar("\xce\x92", L"\x0392");
	(*g_pmapSpecChars)["Gamma"] = SpecChar("\xce\x93", L"\x0393");
	(*g_pmapSpecChars)["Delta"] = SpecChar("\xce\x94", L"\x0394");
	(*g_pmapSpecChars)["Epsilon"] = SpecChar("\xce\x95", L"\x0395");
	(*g_pmapSpecChars)["Zeta"] = SpecChar("\xce\x96", L"\x0396");
	(*g_pmapSpecChars)["Eta"] = SpecChar("\xce\x97", L"\x0397");
	(*g_pmapSpecChars)["Theta"] = SpecChar("\xce\x98", L"\x0398");
	(*g_pmapSpecChars)["Iota"] = SpecChar("\xce\x99", L"\x0399");
	(*g_pmapSpecChars)["Kappa"] = SpecChar("\xce\x9a", L"\x039a");
	(*g_pmapSpecChars)["Lambda"] = SpecChar("\xce\x9b", L"\x039b");
	(*g_pmapSpecChars)["Mu"] = SpecChar("\xce\x9c", L"\x039c");
	(*g_pmapSpecChars)["Nu"] = SpecChar("\xce\x9d", L"\x039d");
	(*g_pmapSpecChars)["Xi"] = SpecChar("\xce\x9e", L"\x039e");
	(*g_pmapSpecChars)["Omicron"] = SpecChar("\xce\x9f", L"\x039f");
	(*g_pmapSpecChars)["Pi"] = SpecChar("\xce\xa0", L"\x03a0");
	(*g_pmapSpecChars)["Rho"] = SpecChar("\xce\xa1", L"\x03a1");
	(*g_pmapSpecChars)["Sigma"] = SpecChar("\xce\xa3", L"\x03a3");
	(*g_pmapSpecChars)["Tau"] = SpecChar("\xce\xa4", L"\x03a4");
	(*g_pmapSpecChars)["Upsilon"] = SpecChar("\xce\xa5", L"\x03a5");
	(*g_pmapSpecChars)["Phi"] = SpecChar("\xce\xa6", L"\x03a6");
	(*g_pmapSpecChars)["Chi"] = SpecChar("\xce\xa7", L"\x03a7");
	(*g_pmapSpecChars)["Psi"] = SpecChar("\xce\xa8", L"\x03a8");
	(*g_pmapSpecChars)["Omega"] = SpecChar("\xce\xa9", L"\x03a9");

	(*g_pmapSpecChars)["leftarrow"] = SpecChar("\xe2\x86\x90", L"\x2190");
	(*g_pmapSpecChars)["rightarrow"] = SpecChar("\xe2\x86\x92", L"\x2192");
	(*g_pmapSpecChars)["uparrow"] = SpecChar("\xe2\x86\x91", L"\x2191");
	(*g_pmapSpecChars)["downarrow"] = SpecChar("\xe2\x86\x93", L"\x2193");
	//(*g_pmapSpecChars)["Leftarrow"] = SpecChar("\xe2\x87\x90", L"\x21d0");
	//(*g_pmapSpecChars)["Rightarrow"] = SpecChar("\xe2\x87\x92", L"\x21d2");
	//(*g_pmapSpecChars)["Uparrow"] = SpecChar("\xe2\x87\x91", L"\x21d1");
	//(*g_pmapSpecChars)["Downarrow"] = SpecChar("\xe2\x87\x93", L"\x21d3");

	(*g_pmapSpecChars)["return"] = SpecChar("\xe2\x8f\x8e", L"\x23ce");
	(*g_pmapSpecChars)["bullet"] = SpecChar("\xe2\x80\xa2", L"\x2022");
}

void deinit_spec_chars()
{
	if(g_pmapSpecChars)
	{
		delete g_pmapSpecChars;
		g_pmapSpecChars = 0;
	}
}

const std::string& get_spec_char_utf8(const std::string& strChar)
{
	static const std::string strNull;
	if(!g_pmapSpecChars)
		init_spec_chars();

	t_mapSpecChars::const_iterator iter = g_pmapSpecChars->find(strChar);
	if(iter == g_pmapSpecChars->end())
		return strNull;

	return iter->second.strUTF8;
}

const std::wstring& get_spec_char_utf16(const std::string& strChar)
{
	static const std::wstring strNull;
	if(!g_pmapSpecChars)
		init_spec_chars();

	t_mapSpecChars::const_iterator iter = g_pmapSpecChars->find(strChar);
	if(iter == g_pmapSpecChars->end())
		return strNull;

	return iter->second.strUTF16;
}

const t_mapSpecChars& get_spec_chars()
{
	if(!g_pmapSpecChars)
		init_spec_chars();

	return *g_pmapSpecChars;
}

}

