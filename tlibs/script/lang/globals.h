/*
 * global symbols
 * @author tweber
 * @date 2013-2014
 * @license GPLv2 or GPLv3
 */

#ifndef __SCRIPT_GLOBALS_H__
#define __SCRIPT_GLOBALS_H__

#include "types.h"
#include "symbol.h"

extern "C" int yydebug;

extern const t_char* g_pcVersion;
extern void init_global_syms(SymbolTable *pSymTab);

#endif
