/*
 * Parse & Runtime Info
 * @author tweber
 * @date 12-oct-14
 * @license GPLv2 or GPLv3
 */

#ifndef __SCRIPT_INFOS_H__
#define __SCRIPT_INFOS_H__

#include "node.h"
#include "handles.h"
#include "symbol.h"
#include "types.h"
#include "lexer.h"


// stuff that can change during execution
struct RuntimeInfo
{
	// function to execute, e.g. "main" (with external local symbol table)
	t_string strExecFkt;
	t_string strInitScrFile;
	SymbolTable *pLocalSymsOverride = nullptr;

	// currently active function
	const NodeFunction *pCurFunction = nullptr;
	const NodeCall *pCurCaller = nullptr;
	bool bWantReturn = 0;

	const Node* pCurLoop = nullptr;
	bool bWantBreak = 0;
	bool bWantContinue = 0;


	// implicitely return last symbol in function
	bool bImplicitRet = 0;


	bool IsExecDisabled() const
	{
		return bWantReturn || bWantBreak || bWantContinue;
	}
	void EnableExec()
	{
		bWantReturn = bWantBreak = bWantContinue = 0;
	}


	RuntimeInfo() = default;
	RuntimeInfo(const RuntimeInfo&) = delete;
	~RuntimeInfo() = default;
};


// stuff that is (more or less) fixed after parsing
struct ParseInfo
{
	// external imported modules
	typedef std::unordered_map<t_string, Node*> t_mods;
	t_mods *pmapModules = nullptr;

	// all functions from all modules
	typedef std::vector<NodeFunction*> t_funcs;
	t_funcs vecFuncs;

	// global symbol table
	SymbolTable *pGlobalSyms = nullptr;
	std::mutex *pmutexGlobalSyms = nullptr;

	HandleManager *phandles = nullptr;


	// mutex for script if no explicit mutex given
	std::mutex *pmutexGlobal = nullptr;


	bool bEnableDebug = 0;
	std::mutex *pmutexTraceback = nullptr;
	typedef std::deque<std::string> t_oneTraceback;
	typedef std::unordered_map<std::thread::id, t_oneTraceback> t_stckTraceback;
	t_stckTraceback stckTraceback;

	void PushTraceback(std::string&& strTrace);
	void PopTraceback();



	ParseInfo();
	~ParseInfo();

	NodeFunction* GetFunction(const t_string& strName);
};


struct ParseObj
{
	Lexer* pLexer;
	Node* pRoot;

	// only used during parsing/lexing for yyerror(), NOT during exec
	unsigned int iCurLine;
	t_string strCurFile;
};

#endif
