/*
 * external functions
 * @author tweber
 * @date 2013-2014
 * @license GPLv2 or GPLv3
 */
#ifndef __EXT_CALLS__
#define __EXT_CALLS__

#include "types.h"

#include <vector>
#include <unordered_map>
#include "symbol.h"
#include "node.h"
#include "info.h"

#include <initializer_list>

// helper functions
extern t_string linenr(const RuntimeInfo &info);
extern const std::string& get_type_name(SymbolType ty);
extern const std::string& get_type_name(unsigned int ty);

// errors
extern bool check_args(RuntimeInfo& info, 
	const std::vector<Symbol*>& vecSyms, 
	const std::initializer_list<unsigned int>& lstTypes, 
	const std::initializer_list<bool> &lstOptional, 
	const char* pcFkt, const char* pcErr=0);

// typedefs
typedef Symbol*(*t_extcall)(const std::vector<Symbol*>&, ParseInfo&, RuntimeInfo&, SymbolTable*);
typedef std::unordered_map<t_string, t_extcall> t_mapFkts;


// adding new external calls
extern bool add_ext_call(const t_string& strFkt, t_extcall pExtCall);
extern void add_ext_calls(t_mapFkts&);

extern bool has_ext_call(const t_string& strFkt);
extern const t_mapFkts* get_ext_calls();

// calling external functions from interpreter
extern Symbol* ext_call(const t_string& strFkt,
	const std::vector<Symbol*>& vecSyms, ParseInfo &info,
	RuntimeInfo &runinfo, SymbolTable* pSymTab);

#endif
