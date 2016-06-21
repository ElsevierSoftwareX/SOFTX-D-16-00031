/*
 * external functions
 * @author tweber
 * @date 2013-2014
 * @license GPLv2 or GPLv3
 */

#include "helper/flags.h"
#include "calls.h"
#include "info.h"
#include "script_helper.h"


// external functions
static t_mapFkts g_mapFkts;


extern t_string linenr(const RuntimeInfo &info)
{
	if(info.pCurCaller)
		return info.pCurCaller->linenr(info);
	return T_STR"";
}


extern const t_string& get_type_name(SymbolType ty)
{
	static const t_string strInvalid = "invalid";

	static const std::unordered_map<SymbolType, t_string, EnumDirectHash<SymbolType>> mapTypes =
		{
			{SYMBOL_DOUBLE, "real"},
			{SYMBOL_INT, "int"},
			{SYMBOL_COMPLEX, "complex"},
			{SYMBOL_STRING, "string"},
			{SYMBOL_ARRAY, "vector"},
			{SYMBOL_MAP, "map"},

			{SYMBOL_SCALAR, "scalar"},
			{SYMBOL_CONTAINER, "container"},

			{SYMBOL_ANY, "any"},
		};

	auto iter = mapTypes.find(ty);
	if(iter == mapTypes.end())
		return strInvalid;

	return iter->second;
}
extern const t_string& get_type_name(unsigned int ty)
{
	return get_type_name(SymbolType(ty));
}

extern bool check_args(RuntimeInfo& runinfo, 
	const std::vector<Symbol*>& vecSyms, 
	const std::initializer_list<unsigned int>& lstTypes, 
	const std::initializer_list<bool> &lstOptional, 
	const char* pcFkt, const char* pcErr)
{
	const std::size_t iSyms = vecSyms.size();
	const std::size_t iTypes = lstTypes.size();
	const std::size_t iTotalOpt = lstOptional.size();
	std::size_t iCompulsory = 0;
	for(bool bOpt : lstOptional)
		if(!bOpt) ++iCompulsory;

	if(iSyms < iCompulsory)
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) 
				<< "Function \"" << pcFkt << "\""
				<< " requires " << iCompulsory
				<< " arguments, but only "
				<< iSyms << " were given.";
		if(pcErr) ostrErr << " " << pcErr;
		ostrErr << std::endl;
		throw tl::Err(ostrErr.str(),0);
	}


	std::vector<Symbol*>::const_iterator iterSym = vecSyms.begin();
	std::initializer_list<unsigned int>::iterator iterTypes = lstTypes.begin();
	std::initializer_list<bool>::iterator iterOptional = lstOptional.begin();

	std::size_t iCurSym = 0;
	for(; iterSym!=vecSyms.end(); ++iterSym, ++iterTypes, ++iterOptional, ++iCurSym)
	{
		// ignore remaining symbols
		if(iCurSym >= iTypes || iCurSym >= iTotalOpt)
			break;

		if(!*iterSym)
		{
			std::ostringstream ostrErr;
			ostrErr << linenr(runinfo)
				<< "Argument " << (iCurSym+1) 
				<< " of function \"" << pcFkt << "\""
				<< " is invalid.";
			if(pcErr) ostrErr << " " << pcErr;
			ostrErr << std::endl;

			throw tl::Err(ostrErr.str(),0);
		}

		if(!(*iterTypes & (*iterSym)->GetType()))
		{
			std::ostringstream ostrErr;
			ostrErr << linenr(runinfo)
				<< "Argument " << (iCurSym+1)
				<< " of function \"" << pcFkt << "\""
				<< " has wrong type. "
				<< "Expected " << get_type_name(*iterTypes)
				<< ", received " << get_type_name((*iterSym)->GetType()) 
				<< ".";
			if(pcErr) ostrErr << " " << pcErr;
			ostrErr << std::endl;

			throw tl::Err(ostrErr.str(),0);
		}
	}

	return 1;
}


extern Symbol* ext_call(const t_string& strFkt,
	const std::vector<Symbol*>& vecSyms, ParseInfo& info,
	RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	t_mapFkts::iterator iter = g_mapFkts.find(strFkt);
	if(iter == g_mapFkts.end())
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo)
			<< "Tried to call unknown function \""
			<< strFkt << "\"."
			<< std::endl;
		throw tl::Err(ostrErr.str(),0);
	}

	return (*iter).second(vecSyms, info, runinfo, pSymTab);
}

// --------------------------------------------------------------------------------

extern bool add_ext_call(const t_string& strFkt, t_extcall pExtCall)
{
	bool bInserted = g_mapFkts.insert(t_mapFkts::value_type(strFkt, pExtCall)).second;
	return bInserted;
}

extern void add_ext_calls(t_mapFkts& mapFkt)
{
	g_mapFkts.insert(mapFkt.begin(), mapFkt.end());
}

extern bool has_ext_call(const t_string& strFkt)
{
	return g_mapFkts.find(strFkt) != g_mapFkts.end();
}

extern const t_mapFkts* get_ext_calls()
{
	return &g_mapFkts;
}
