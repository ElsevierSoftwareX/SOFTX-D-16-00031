/*
 * Hermelin script interpreter
 * @author tweber
 * @date 2013
 * @license GPLv2 or GPLv3
 */

#include "types.h"

#include "helper/flags.h"
#include "string/string.h"
#include "string/spec_char.h"
#include "log/log.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <array>
#include <algorithm>
#include <boost/version.hpp>

#include "script_helper.h"
#include "globals.h"
#include "calls.h"
#include "info.h"

extern int yyparse(void*);
static bool g_bShowTiming = 0;


static inline void usage(const char* pcProg)
{
	G_COUT << "This is the " << g_pcVersion << ".\n";
	G_COUT << "Written by Tobias Weber, 2013-2014.\n";
	G_COUT << "Built on " << __DATE__ << ", " << __TIME__;
	//G_COUT << " with CC version " << __VERSION__ << ".\n\n";
	G_COUT << " with " << BOOST_COMPILER << ".\n";
	G_COUT << "Internal data type sizes: int = " << sizeof(t_int)*8 << " bits, "
		<< "real = " << sizeof(t_real)*8 << " bits.\n\n";
	G_COUT << "Usage: " << pcProg << " [arguments to hermelin]"
		<< " <script file> [arguments to script]"
		<< "\n";

	G_COUT << "\n" <<
	R"RAW(Arguments to hermelin:
	-h, --help            This message.
	-i, --interactive     Interactive mode.
	-t, --timing          Show timing information.
	-s, --symbols         Show symbol tables.
	-d[0-4]               Verbosity (0=none, 1=errors, 2=warnings, 3=infos, 4=debug).
	)RAW";

	G_COUT << std::endl;
}



// execute interactive commands
static inline int interactive(bool bShowSymbols=0, unsigned int uiDebugLevel=3)
{
	static const t_char* pcCmdFunc = T_STR"__cmd__";

	ParseInfo info;
	RuntimeInfo runinfo;

	std::unique_ptr<SymbolTable> ptrLocalSym(new SymbolTable);
	runinfo.pLocalSymsOverride = ptrLocalSym.get();
	runinfo.bImplicitRet = 1;
	runinfo.strExecFkt = pcCmdFunc;
	runinfo.strInitScrFile = "<interactive>";

	init_global_syms(info.pGlobalSyms);

	std::function<void()> remove_cmdfunc = [&info]()
	{
		ParseInfo::t_funcs::iterator iterNewEnd =
			std::remove_if(info.vecFuncs.begin(), info.vecFuncs.end(),
			[](const NodeFunction* pFunc)->bool
			{ return pFunc->GetName() == pcCmdFunc; } );
		info.vecFuncs.resize(iterNewEnd-info.vecFuncs.begin());
	};

	while(1)
	{
		try
		{
			std::cout << "> ";
			std::string strLine;
			if(!std::getline(std::cin, strLine))
				break;

			runinfo.EnableExec();

			t_string strInput = t_string(pcCmdFunc) + "() { " + strLine + " }";
			Lexer lex(strInput);

			if(!lex.IsOk())
			{
				tl::log_err("Lexer returned with errors.");
				continue;
			}


			ParseObj par;
			info.bEnableDebug = (uiDebugLevel>=4);
			par.pLexer = &lex;

			int iParseRet = yyparse(&par);
			if(iParseRet != 0)
			{
				tl::log_err("Parser returned with error code ", iParseRet, ".");
				remove_cmdfunc();
				continue;
			}

			par.pRoot = par.pRoot->optimize();
			Symbol *pSymRet = par.pRoot->eval(info, runinfo, 0);

			if(bShowSymbols)
				runinfo.pLocalSymsOverride->print();

			if(pSymRet)
			{
				std::cout << pSymRet->print() << std::endl;
				safe_delete(pSymRet, runinfo.pLocalSymsOverride, &info);
			}

			remove_cmdfunc();

			if(par.pRoot)
			{
				delete par.pRoot;
				par.pRoot = 0;
			}
		}
		catch(const std::exception& ex)
		{
			tl::log_crit(ex.what());
			remove_cmdfunc();
		}
	}

	return 0;
}



// execute script
static inline int script_main(int argc, char** argv)
{
	if(argc<=1)
	{
		usage(argv[0]);
		return -1;
	}

	bool bShowSymbols = 0;
	bool bInteractive = 0;
	unsigned int uiDebugLevel = 3;
#ifndef NDEBUG
	uiDebugLevel = 4;
#endif
	unsigned int iStartArg = 1;
	for(iStartArg=1; iStartArg<unsigned(argc); ++iStartArg)
	{
		t_string strArg = STR_TO_WSTR(argv[iStartArg]);
		tl::trim(strArg);

		// end of arguments to hermelin
		if(strArg[0] != T_STR'-')
			break;

		if(strArg=="-s" || strArg == "--symbols")
			bShowSymbols = 1;
		else if(strArg=="-i" || strArg == "--interactive")
			bInteractive = 1;
		else if(strArg=="-h" || strArg == "--help")
			{ usage(argv[0]); return 0; }

		else if(strArg=="-d0") uiDebugLevel = 0;
		else if(strArg=="-d1") uiDebugLevel = 1;
		else if(strArg=="-d2") uiDebugLevel = 2;
		else if(strArg=="-d3") uiDebugLevel = 3;
		else if(strArg=="-d4") uiDebugLevel = 4;
	}

	const std::array<tl::Log*, 5> arrLogs{{&tl::log_crit, &tl::log_err, &tl::log_warn, &tl::log_info, &tl::log_debug}};
	for(unsigned int iLog=0; iLog<arrLogs.size(); ++iLog)
		arrLogs[iLog]->SetEnabled(uiDebugLevel>=iLog);

	// debug in script.yy needs to be set
	yydebug = (uiDebugLevel>=4);

	if(bInteractive)
		return interactive(bShowSymbols, uiDebugLevel);


	if(iStartArg >= unsigned(argc))
	{
		tl::log_err("No input file given.");
		return -1;
	}



	// loading of input file
	const char* pcFile = argv[iStartArg];
	t_string strFile = STR_TO_WSTR(pcFile);

	t_char* pcInput = load_file(pcFile);
	if(!pcInput)
		return -2;


	ParseObj par;
	ParseInfo info;
	RuntimeInfo runinfo;

	info.bEnableDebug = (uiDebugLevel>=4);


	// lexing
	par.strCurFile = strFile;
	par.pLexer = new Lexer(pcInput, strFile.c_str());

	delete[] pcInput;
	pcInput = 0;

	if(!par.pLexer->IsOk())
	{
		tl::log_err("Lexer returned with errors.");
		return -3;
	}

	init_global_syms(info.pGlobalSyms);


	// parsing
	int iParseRet = yyparse(&par);

	delete par.pLexer;
	par.pLexer = 0;

	if(iParseRet != 0)
	{
		tl::log_err("Parser returned with error code ", iParseRet, ".");
		return -4;
	}


	// optimizing
	par.pRoot = par.pRoot->optimize();



	// executing
	SymbolArray *parrMainArgs = new SymbolArray();
	for(int iArg=iStartArg; iArg<argc; ++iArg)
	{
		SymbolString *pSymArg = new SymbolString();
		pSymArg->SetVal(STR_TO_WSTR(argv[iArg]));
		parrMainArgs->GetArr().push_back(pSymArg);
	}

	SymbolTable *pTableSup = new SymbolTable();

	info.pmapModules->insert(ParseInfo::t_mods::value_type(strFile, par.pRoot));
	runinfo.strExecFkt = T_STR"main";
	runinfo.strInitScrFile = strFile;

	SymbolArray arrMainArgs;
	arrMainArgs.GetArr().push_back(parrMainArgs);
	pTableSup->InsertSymbol(T_STR"<args>", &arrMainArgs);
	par.pRoot->eval(info, runinfo, pTableSup);
	pTableSup->RemoveSymbolNoDelete(T_STR"<args>");
	delete pTableSup;


	if(bShowSymbols)
	{
		tl::log_info("================================================================================");
		tl::log_info("Global symbols:");
		info.pGlobalSyms->print();

		std::ostringstream ostrFkts;
		for(const NodeFunction* pFunc : info.vecFuncs)
			ostrFkts << pFunc->GetName() << ", ";
		tl::log_info("Script functions: ", ostrFkts.str());


		const t_mapFkts* pExtFkts = get_ext_calls();

		std::ostringstream ostrSysFkts;
		for(const auto& fktpair : *pExtFkts)
			ostrSysFkts << fktpair.first << ", ";
		tl::log_info("System functions: ", ostrSysFkts.str());
		tl::log_info("================================================================================");
	}

	return 0;
}



#include "time/stopwatch.h"

int main(int argc, char** argv)
{
	// global flags
	for(int iStartArg=1; iStartArg<argc; ++iStartArg)
	{
		t_string strArg = STR_TO_WSTR(argv[iStartArg]);
		tl::trim(strArg);

		// end of arguments to hermelin
		if(strArg[0] != T_STR'-')
			break;

		if(strArg=="-t" || strArg == "--timing")
			g_bShowTiming = 1;
	}


	const std::array<tl::Log*, 5> arrLogs{{&tl::log_crit, &tl::log_err, &tl::log_warn, &tl::log_info, &tl::log_debug}};
	for(tl::Log* pLog : arrLogs)
	{
		pLog->SetShowDate(0);
		pLog->SetShowThread(0);
	}

	int iRet = -99;

	tl::Stopwatch<t_real> watch;

	try
	{
		tl::init_spec_chars();
		if(g_bShowTiming) watch.start();
		iRet = script_main(argc, argv);
		if(g_bShowTiming) watch.stop();
	}
	catch(const std::exception& ex)
	{
		tl::log_crit(ex.what());
	}

	if(g_bShowTiming && iRet==0)
	{
		tl::log_info("================================================================================");
		tl::log_info("Script start time:     ", watch.GetStartTimeStr());
		tl::log_info("Script stop time:      ", watch.GetStopTimeStr());
		tl::log_info("Script execution time: ", watch.GetDur(), " s");
		tl::log_info("================================================================================");
	}

	return iRet;
}
