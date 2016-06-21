/*
 * basic external functions
 * @author tweber
 * @date dec-2013
 * @license GPLv2 or GPLv3
 */

#include "helper/flags.h"
#include "log/log.h"
#include "calls_basic.h"
#include "lang/types.h"
#include "lang/calls.h"

#include "lang/info.h"
#include "lang/script_helper.h"
#include "string/string.h"

#include <sstream>
#include <iostream>
#include <unordered_map>
#include <map>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <utility>

#include <ratio>
#include <chrono>
#include <thread>

#ifndef USE_BOOST_REX
	#include <regex>
	namespace rex = ::std;
#else
	#include <boost/tr1/regex.hpp>
	namespace rex = ::boost;
#endif


static Symbol* fkt_version(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	extern const t_char* g_pcVersion;
	return new SymbolString(g_pcVersion);
}

static Symbol* fkt_traceback(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!info.bEnableDebug)
		tl::log_warn(linenr(runinfo), "Debugging is disabled.");

	ParseInfo::t_stckTraceback& mapstck = info.stckTraceback;
	typedef typename ParseInfo::t_stckTraceback::value_type::second_type t_stck;
	t_stck& stck = mapstck[std::this_thread::get_id()];

	SymbolArray *pArr = new SymbolArray();
	unsigned int iCur = 0;
	for(t_stck::const_iterator iter = stck.begin(); iter!=stck.end(); ++iter)
	{
		const t_string& str = *iter;

		SymbolString *pStr = new SymbolString(std::to_string(iCur) + ": " + str);
		pArr->GetArr().push_back(pStr);

		++iCur;
	}

	pArr->UpdateIndices();
	return pArr;
}

static Symbol* fkt_is_valid(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	for(const Symbol* pSym : vecSyms)
	{
		if(!pSym)
			return new SymbolInt(0);
	}
	return new SymbolInt(1);
}

static Symbol* fkt_is_rval(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_ANY}, {0}, "is_rval", "Symbol needed."))
		return 0;

	return new SymbolInt(vecSyms[0]->IsRval());
}

static Symbol* fkt_null(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	return 0;
}

static Symbol* fkt_getargnames(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_STRING}, {0}, "get_argnames", "Function name needed."))
		return 0;

	const t_string& strFkt = ((SymbolString*)vecSyms[0])->GetVal();
	NodeFunction *pFkt = info.GetFunction(strFkt);
	if(!pFkt)
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "Unknown user function \"" << strFkt << "\"." << std::endl;
		throw tl::Err(ostrErr.str(),0);
	}


	SymbolArray *pSymArr = new SymbolArray();

	std::vector<t_string> vecParams = pFkt->GetParamNames();
	for(const t_string& strParam : vecParams)
		pSymArr->GetArr().push_back(new SymbolString(strParam));
	pSymArr->UpdateIndices();

	return pSymArr;
}

// call("fkt", [arg1, arg2, ...]);
static Symbol* fkt_call(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_STRING, SYMBOL_ARRAY|SYMBOL_MAP}, {0,1}, "call", "Function name needed."))
		return 0;

	const t_string& strFkt = ((SymbolString*)vecSyms[0])->GetVal();
	NodeFunction *pFkt = info.GetFunction(strFkt);

	SymbolArray *pArr = 0;
	if(vecSyms.size()>=2 && vecSyms[1]->GetType()==SYMBOL_ARRAY)
		pArr = (SymbolArray*)vecSyms[1];

	Symbol *pSymRet = 0;
	if(pFkt)	// user function
	{
		// map to function args
		if(vecSyms.size()>=2 && vecSyms[1]->GetType() == SYMBOL_MAP)
		{
			if(!pArr) pArr = new SymbolArray();

			SymbolMap::t_map& mapSyms = ((SymbolMap*)vecSyms[1])->GetMap();
			const std::vector<t_string> vecParamNames = pFkt->GetParamNames();

			for(unsigned int iArg=0; iArg<vecParamNames.size(); ++iArg)
			{
				const t_string& strParamName = vecParamNames[iArg];
				SymbolMap::t_map::iterator iter = mapSyms.find(SymbolMapKey(strParamName));
				if(iter == mapSyms.end())
				{
					tl::log_err(linenr(runinfo), "Parameter \"", strParamName, "\" not in argument map. Using 0.");
					pArr->GetArr().push_back(new SymbolReal(0.));
					continue;
				}

				pArr->GetArr().push_back(iter->second->clone());
			}
			pArr->UpdateIndices();
		}

		pSymTab->InsertSymbol(T_STR"<args>", pArr);
		pSymRet = pFkt->eval(info, runinfo, pSymTab);
		runinfo.bWantReturn = 0;
		pSymTab->RemoveSymbolNoDelete(T_STR"<args>");
	}
	else		// system function
	{
		pSymRet = ext_call(strFkt, pArr ? pArr->GetArr() : std::vector<Symbol*>(),
			  info, runinfo, pSymTab);
	}

	return pSymRet;
}

static Symbol* fkt_int(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(vecSyms.size() == 0)
		return new SymbolInt(0);

	return vecSyms[0]->ToType(SYMBOL_INT);
}

static Symbol* fkt_cplx_real(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab);

static Symbol* fkt_double(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(vecSyms.size() == 0)
		return new SymbolReal(0.);

	// real([1,2,3])
	else if(vecSyms[0]->GetType() == SYMBOL_ARRAY)
	{
		SymbolArray *pArr = new SymbolArray();

		for(Symbol* pSym : ((SymbolArray*)vecSyms[0])->GetArr())
		{
			std::vector<Symbol*> vecTmp{pSym};
			Symbol* pSymElem = fkt_double(vecTmp, info, runinfo, pSymTab);

			pArr->GetArr().push_back(pSymElem);
		}

		pArr->UpdateIndices();
		return pArr;
	}

	// for complex args delegate to other "real" function
	else if(vecSyms[0]->GetType() == SYMBOL_COMPLEX)
		return fkt_cplx_real(vecSyms, info, runinfo, pSymTab);

	// cast arg to real value
	return vecSyms[0]->ToType(SYMBOL_DOUBLE);
}

static Symbol* fkt_double_vec(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_ARRAY}, {0}, "real_vec"))
		return 0;

	SymbolArray* pSymRet = new SymbolArray();
	pSymRet->GetArr().reserve(vecSyms.size());

	for(Symbol *pSym : ((SymbolArray*)vecSyms[0])->GetArr())
	{
		Symbol *pSymCast = 0;

		if(pSym->GetType() == SYMBOL_ARRAY)
			pSymCast = fkt_double_vec(((SymbolArray*)pSym)->GetArr(), info, runinfo, pSymTab);
		else
		{
			std::vector<Symbol*> vecSymTmp{pSym};
			pSymCast = fkt_double(vecSymTmp, info, runinfo, pSymTab);
		}

		if(pSymCast)
			pSymRet->GetArr().push_back(pSymCast);
	}

	pSymRet->UpdateIndices();
	return pSymRet;
}

static Symbol* fkt_str(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	SymbolString *pSymRet = new SymbolString;
	for(const Symbol *pSym : vecSyms)
	{
		if(!pSym)
			continue;
		pSymRet->GetVal() += pSym->print();
	}
	return pSymRet;
}

static Symbol* fkt_output(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	t_ostream& ostr = G_COUT;

	for(Symbol *pSym : vecSyms)
		if(pSym)
			ostr << pSym->print();

	ostr.flush();
	return 0;
}

static Symbol* fkt_print(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	fkt_output(vecSyms, info, runinfo, pSymTab);
	G_COUT << std::endl;
	return 0;
}

static Symbol* fkt_input(const std::vector<Symbol*>& vecSyms,
		ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	fkt_output(vecSyms, info, runinfo, pSymTab);

	SymbolString* pSymStr = new SymbolString();
	std::getline(G_CIN, pSymStr->GetVal());
	return pSymStr;
}

static Symbol* fkt_sleep(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_DOUBLE|SYMBOL_INT}, {0}, "sleep", "Number of ms needed."))
		return 0;

	t_real dMillis = vecSyms[0]->GetValDouble();

	typedef std::chrono::duration<t_real, std::milli> t_millis;
	t_millis delay(dMillis);
	std::this_thread::sleep_for(delay);
	return 0;
}


extern int yyparse(void*);
static bool _import_file(const t_string& strFile, ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	ParseInfo::t_mods::iterator iterMod = info.pmapModules->find(strFile);
	if(iterMod != info.pmapModules->end())
	{
		//G_CERR << "Warning: Module \"" << strFile << "\" already loaded." << std::endl;
		return 0;
	}

	t_string _strFile = WSTR_TO_STR(strFile);
	t_char* pcInput = load_file(_strFile.c_str());
	if(!pcInput)
		return 0;

	ParseObj par;
	par.strCurFile = strFile;
	par.pLexer = new Lexer(pcInput, strFile.c_str());

	delete[] pcInput;
	pcInput = 0;

	if(!par.pLexer->IsOk())
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "Lexer returned with errors"
			<< " for file \"" << _strFile << "\"."<< std::endl;
		throw tl::Err(ostrErr.str(),0);
	}

	int iParseRet = yyparse(&par);

	delete par.pLexer;
	par.pLexer = 0;

	if(iParseRet != 0)
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "Parser returned with error code "
			<< iParseRet
			<< " for file \"" << _strFile << "\"." << std::endl;
		throw tl::Err(ostrErr.str(),0);
	}

	Node *pRoot = par.pRoot;
	if(pRoot) pRoot = pRoot->optimize();


	info.pmapModules->insert(ParseInfo::t_mods::value_type(strFile, pRoot));

	t_string strExecOld = runinfo.strExecFkt;
	t_string strFileOld = runinfo.strInitScrFile;
	runinfo.strExecFkt = T_STR"";
	runinfo.strInitScrFile = strFile;
	pRoot->eval(info, runinfo);
	runinfo.strExecFkt = strExecOld;
	runinfo.strInitScrFile = strFileOld;

	return 1;
}

static Symbol* fkt_import(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo& runinfo, SymbolTable* pSymTab)
{
	for(Symbol *pSym : vecSyms)
		if(pSym && pSym->GetType()==SYMBOL_STRING)
		{
			const t_string& strFile = ((SymbolString*)pSym)->GetVal();
			bool bOk = _import_file(strFile, info, runinfo, pSymTab);

#ifdef INSTALL_PREFIX
			if(!bOk)
			{
				t_string strInstPrefix = t_string(INSTALL_PREFIX) + "/share/hermelin/";
				tl::log_debug("Importing from ", strInstPrefix);

				bOk = _import_file(strInstPrefix + strFile, info, runinfo, pSymTab);
			}
#endif
		}

	return 0;
}

static Symbol* fkt_has_var(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_STRING}, {0}, "has_var", "Symbol name needed."))
		return 0;

	const t_string& strVar = ((SymbolString*)vecSyms[0])->GetVal();
	bool bHasVar = 0;

	// check local variables
	if(pSymTab->GetSymbol(strVar))
		bHasVar = 1;

	{
		// check global variables
		std::lock_guard<std::mutex> _lck(*info.pmutexGlobalSyms);

		if(info.pGlobalSyms->GetSymbol(strVar))
			bHasVar = 1;
	}

	return new SymbolInt(bHasVar);
}

// register a variable in the symbol table
static Symbol* fkt_register_var(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	bool bUseGlobal = 0;

	if(vecSyms.size() == 2)
	{
		if(!check_args(runinfo, vecSyms, {SYMBOL_STRING, SYMBOL_ANY}, {0,0}, "register_var", "Symbol name and symbol needed."))
			return 0;

		const t_string& strVar = ((SymbolString*)vecSyms[0])->GetVal();
		Symbol* pVar = vecSyms[1]->clone();

		if(bUseGlobal)
		{
			std::lock_guard<std::mutex> _lck(*info.pmutexGlobalSyms);
			info.pGlobalSyms->InsertSymbol(strVar, pVar);
		}
		else
		{
			pSymTab->InsertSymbol(strVar, pVar);
		}
	}
	else if(vecSyms.size() == 1)
	{
		if(!check_args(runinfo, vecSyms, {SYMBOL_MAP}, {0}, "register_var", "Map needed."))
			return 0;

		SymbolMap::t_map& varmap = ((SymbolMap*)vecSyms[0])->GetMap();
		for(SymbolMap::t_map::value_type& val : varmap)
		{
			const t_string& strKey = val.first.strKey;
			if(strKey == "")
				continue;

			SymbolString symKey(strKey);
			Symbol *pSymVal = val.second;

			std::vector<Symbol*> vecDummy;
			vecDummy.resize(2);

			vecDummy[0] = &symKey;
			vecDummy[1] = pSymVal;

			fkt_register_var(vecDummy, info, runinfo, pSymTab);
		}
	}
	else
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "Invalid call to register_var." << std::endl;
		throw tl::Err(ostrErr.str(),0);
	}
	return 0;
}

static Symbol* fkt_typeof(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_ANY}, {0}, "typeof"))
		return 0;

	Symbol *pSymbol = vecSyms[0];
	SymbolString *pType = new SymbolString(pSymbol->GetTypeName().c_str());
	return pType;
}

static Symbol* fkt_setprec(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_SCALAR}, {0}, "set_prec"))
		return 0;

	Symbol *pSymbol = vecSyms[0];
	int iPrec = pSymbol->GetValInt();

	if(iPrec >= 0)
		SymbolReal::SetPrec(iPrec);
	else
		SymbolReal::SetPrec(SymbolReal::GetDefPrec());

	return 0;
}


// --------------------------------------------------------------------------------
// map

static Symbol* fkt_map(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(vecSyms.size()>=1 && vecSyms[0]->GetType()==SYMBOL_MAP)
		return vecSyms[0]->clone();

	return new SymbolMap();
}
// --------------------------------------------------------------------------------



// --------------------------------------------------------------------------------
// complex
static Symbol* fkt_complex(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(vecSyms.size() == 0)
		return new SymbolComplex(0., 0.);
	// complex(z)
	else if(vecSyms.size() == 1 && vecSyms[0]->GetType()==SYMBOL_COMPLEX)
		return vecSyms[0]->clone();
	// complex([1,2,3])
	else if(vecSyms.size() >= 1 && vecSyms[0]->GetType()==SYMBOL_ARRAY)
	{
		SymbolArray *pSymArr = new SymbolArray();
		for(const Symbol* pSym : ((SymbolArray*)vecSyms[0])->GetArr())
		{
			std::vector<Symbol*> vecTmp{const_cast<Symbol*>(pSym)};
			for(unsigned int iRest=1; iRest<vecSyms.size(); ++iRest)
				vecTmp.push_back(vecSyms[iRest]);

			Symbol *pSymNew = fkt_complex(vecTmp, info, runinfo, pSymTab);
			pSymArr->GetArr().push_back(pSymNew);
		}
		pSymArr->UpdateIndices();

		return pSymArr;
	}
	// complex(1)
	else if(vecSyms.size() == 1)
		return new SymbolComplex(vecSyms[0]->GetValDouble(), 0.);
	// complex(1,2)
	else if(vecSyms.size() > 1)
		return new SymbolComplex(vecSyms[0]->GetValDouble(), vecSyms[1]->GetValDouble());

	return 0;
}

static Symbol* fkt_complex_polar(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	t_real tRho = 0.;
	t_real tTheta = 0.;

	if(vecSyms.size() > 0 && vecSyms[0])
		tRho = vecSyms[0]->GetValDouble();
	if(vecSyms.size() > 1 && vecSyms[1])
		tTheta = vecSyms[1]->GetValDouble();

	return new SymbolComplex(std::polar<t_real>(tRho, tTheta));
}

static Symbol* fkt_cplx_real(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_COMPLEX}, {0}, "real"))
		return 0;

	SymbolComplex *pComplex = (SymbolComplex*)vecSyms[0];
	return new SymbolReal(pComplex->GetValReal());
}

static Symbol* fkt_cplx_imag(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_COMPLEX}, {0}, "imag"))
		return 0;

	SymbolComplex *pComplex = (SymbolComplex*)vecSyms[0];
	return new SymbolReal(pComplex->GetValImag());
}

static Symbol* fkt_cplx_arg(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_COMPLEX}, {0}, "carg"))
		return 0;

	SymbolComplex *pComplex = (SymbolComplex*)vecSyms[0];
	return new SymbolReal(std::arg<t_real>(pComplex->GetVal()));
}

// --------------------------------------------------------------------------------



// --------------------------------------------------------------------------------
// array

static Symbol* fkt_array(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(vecSyms.size()<1)
		return new SymbolArray();

	if(vecSyms.size()>=1 && vecSyms[0]->GetType()==SYMBOL_ARRAY)
		return vecSyms[0]->clone();


	Symbol *pSymSize = vecSyms[0];
	if(pSymSize->GetType() != SYMBOL_INT)
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "\"num\" in vec(num, val=0) has to be integer." << std::endl;
		throw tl::Err(ostrErr.str(),0);
	}

	t_int iVal = ((SymbolInt*)pSymSize)->GetVal();
	if(iVal < 0) iVal = 0;



	bool bOwnVal = 0;
	Symbol *pSymVal = 0;
	if(vecSyms.size()>1)
	{
		pSymVal = vecSyms[1];
	}
	else
	{
		pSymVal = new SymbolReal(0.);
		bOwnVal = 1;
	}


	SymbolArray* pSymRet = new SymbolArray;
	pSymRet->GetArr().reserve(iVal);
	for(t_int i=0; i<iVal; ++i)
		pSymRet->GetArr().push_back(pSymVal->clone());

	if(bOwnVal)
		delete pSymVal;

	pSymRet->UpdateIndices();
	return pSymRet;
}

static Symbol* fkt_array_size(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_CONTAINER}, {0}, "vec_size"))
		return 0;

	Symbol *pSymArr = vecSyms[0];
	SymbolInt *pSymRet = new SymbolInt(0);

	if(pSymArr->GetType() == SYMBOL_ARRAY)
		pSymRet->SetVal(((SymbolArray*)pSymArr)->GetArr().size());
	else if(pSymArr->GetType() == SYMBOL_STRING)
		pSymRet->SetVal(((SymbolString*)pSymArr)->GetVal().length());
	else if(pSymArr->GetType() == SYMBOL_MAP)
		pSymRet->SetVal(((SymbolMap*)pSymArr)->GetMap().size());
	else
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "vec_size needs a vector type argument." << std::endl;
		throw tl::Err(ostrErr.str(),0);
	}

	return pSymRet;
}

static int pos_in_string(const SymbolString* pSymStr, const SymbolString *pSym)
{
	const t_string& str = pSymStr->GetVal();
	const t_string& strToFind = pSym->GetVal();

	std::size_t pos = str.find(strToFind);
	if(pos == t_string::npos)
		return -1;

	return int(pos);
}

static int pos_in_array(const SymbolArray* pSymArr, const Symbol *pSym)
{
	unsigned int iPos;
	bool bFound = 0;

	for(iPos=0; iPos<pSymArr->GetArr().size(); ++iPos)
	{
		const Symbol& sym = *pSymArr->GetArr()[iPos];
		bool bSym0Scalar = sym.IsScalar() || sym.GetType()==SYMBOL_STRING;
		bool bSym1Scalar = pSym->IsScalar() || pSym->GetType()==SYMBOL_STRING;

		if(bSym0Scalar != bSym1Scalar)
			continue;
		else if(bSym0Scalar && bSym1Scalar)
		{
			Symbol *pEqu = Node::Op(&sym, pSym, NODE_LOG_EQ);
			if(!pEqu) continue;
			int iEqu = pEqu->GetValInt();
			delete pEqu;

			if(iEqu)
			{
				bFound = 1;
				break;
			}
		}
		else
		{
			std::ostringstream ostrErr;
			ostrErr << "Error: Array compare not yet implemented." << std::endl;
			throw(ostrErr.str(),0);

			// TODO: array/map compare
		}
	}

	if(!bFound) return -1;
	return int(iPos);
}

static Symbol* fkt_find(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_CONTAINER, SYMBOL_ANY}, {0,0}, "find"))
		return 0;

	const Symbol *pContainer = vecSyms[0];
	const Symbol *pContainee = vecSyms[1];

	if(pContainer->GetType() == SYMBOL_ARRAY)
	{
		int iPos = pos_in_array((SymbolArray*)pContainer, pContainee);
		return new SymbolInt(iPos);
	}
	else if(pContainer->GetType() == SYMBOL_MAP)
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "Find not yet implemented for map." << std::endl;
		throw tl::Err(ostrErr.str(),0);

		// TODO: Implement and also adapt fkt_contains
	}
	else if(pContainer->GetType() == SYMBOL_STRING)
	{
		if(pContainee->GetType() != SYMBOL_STRING)
		{
			std::ostringstream ostrErr;
			ostrErr << linenr(runinfo)
				<< "Second argument to find has to be of string type."
				<< std::endl;
			throw tl::Err(ostrErr.str(),0);
		}

		int iPos = pos_in_string((SymbolString*)pContainer, (SymbolString*)pContainee);
		return new SymbolInt(iPos);
	}

	return 0;
}

static Symbol* fkt_contains(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	int iIdx = -1;
	Symbol *pFind = fkt_find(vecSyms, info, runinfo, pSymTab);
	if(pFind)
	{
		iIdx = pFind->GetValInt();
		delete pFind;
	}

	return new SymbolInt(iIdx>=0);
}

static Symbol* fkt_cur_iter(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_ANY}, {0}, "cur_iter"))
		return 0;

	const t_string& strIdent = vecSyms[0]->GetIdent();
	if(strIdent == T_STR"")
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "No identifier given for cur_iter." << std::endl;
		throw tl::Err(ostrErr.str(),0);
	}


	t_string strIter = T_STR"<cur_iter_" + strIdent + T_STR">";

	Symbol* pSymIter = pSymTab->GetSymbol(strIter);
	if(!pSymIter || pSymIter->GetType()!=SYMBOL_INT)
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "cur_iter could not determine iteration index \""
					<< strIter << "\"." << std::endl;
		throw tl::Err(ostrErr.str(),0);
	}

	return pSymIter;
}


static Symbol* fkt_splice(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(vecSyms.size() == 0)
	{
		tl::log_err(linenr(runinfo), "Need argument for splice.");
		return nullptr;
	}
	if(vecSyms[0] == nullptr)
	{
		tl::log_err(linenr(runinfo), "Need valid argument for splice.");
		return nullptr;
	}

	Symbol *pSymRet = nullptr;
	const SymbolType tyFirstSym = vecSyms[0]->GetType();

	if(tyFirstSym == SYMBOL_ARRAY)
		pSymRet = new SymbolArray();
	else if(tyFirstSym == SYMBOL_MAP)
		pSymRet = new SymbolMap();
	else
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "Cannot splice symbols of type "
			<< vecSyms[0]->GetTypeName() << "." << std::endl;
		throw tl::Err(ostrErr.str(),0);
	}


	for(const Symbol* pSym : vecSyms)
	{
		if(!pSym || pSym->GetType() != tyFirstSym)
		{
			tl::log_warn(linenr(runinfo), "Type mismatch for symbol \"", pSym->GetName(), "\" in splice. Ignoring symbol.");
			continue;
		}

		if(tyFirstSym == SYMBOL_ARRAY)
		{
			for(const Symbol* pSymInArr : ((SymbolArray*)pSym)->GetArr())
				((SymbolArray*)pSymRet)->GetArr().push_back(pSymInArr->clone());
		}
		else if(tyFirstSym == SYMBOL_MAP)
		{
			typedef SymbolMap::t_map::value_type t_pair;

			for(const t_pair& pair : ((SymbolMap*)pSym)->GetMap())
				((SymbolMap*)pSymRet)->GetMap().insert(t_pair(pair.first, pair.second->clone()));
		}
	}

	if(tyFirstSym == SYMBOL_ARRAY)
		((SymbolArray*)pSymRet)->UpdateIndices();
	else if(tyFirstSym == SYMBOL_MAP)
		((SymbolMap*)pSymRet)->UpdateIndices();
	return pSymRet;
}

static Symbol* fkt_append(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_ARRAY, SYMBOL_ANY}, {0,0}, "append"))
		return 0;
	//if(vecSyms[1].GetType() == SYMBOL_ARRAY)
	//	return fkt_splice(vecSyms, info, runinfo, pSymTab);

	((SymbolArray*)vecSyms[0])->GetArr().push_back(vecSyms[1]->clone());
	((SymbolArray*)vecSyms[0])->UpdateIndices();

	return vecSyms[0];
}



typedef std::tuple<const Symbol*, unsigned int> t_symtup;

static void _sortarr(std::vector<t_symtup>& vec, bool bReverse=0)
{
	auto comp = [bReverse](const t_symtup& tup1, const t_symtup& tup2) -> bool
	{
		const Symbol* pSym0 = std::get<0>(tup1);
		const Symbol* pSym1 = std::get<0>(tup2);

		bool bLess = pSym0->IsLessThan(*pSym1);
		if(bReverse) bLess = !bLess;

		return bLess;
	};

	std::sort(vec.begin(), vec.end(), comp);
}

static void _rearrangearr(std::vector<Symbol*>& vecSyms, const std::vector<unsigned int>& vecIdx)
{
	std::vector<Symbol*> vecTmp = vecSyms;
	for(unsigned int i=0; i<vecSyms.size(); ++i)
		vecSyms[i] = vecTmp[vecIdx[i]];
}

static Symbol* _fkt_sort(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab,
	bool bReverse=0)
{
	if(vecSyms.size()<1 || vecSyms[0]->GetType()!=SYMBOL_ARRAY)
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo)
			<< "Arguments to sort have to be arrays."
			<< std::endl;
		throw tl::Err(ostrErr.str(),0);
	}

	const SymbolArray* pArr = (SymbolArray*)vecSyms[0];
	const unsigned int iArrSize = pArr->GetArr().size();

	std::vector<t_symtup> vecTups;
	vecTups.reserve(iArrSize);
	for(unsigned int iElem=0; iElem<iArrSize; ++iElem)
		vecTups.push_back(t_symtup(pArr->GetArr()[iElem], iElem));

	_sortarr(vecTups, bReverse);

	SymbolArray* pArrRet = new SymbolArray();
	pArrRet->GetArr().reserve(iArrSize);

	std::vector<unsigned int> vecSortedIndices;
	vecSortedIndices.reserve(iArrSize);

	for(unsigned int iElem=0; iElem<iArrSize; ++iElem)
	{
		pArrRet->GetArr().push_back(std::get<0>(vecTups[iElem])->clone());
		vecSortedIndices.push_back(std::get<1>(vecTups[iElem]));
	}

	// no other arguments to sort
	if(vecSyms.size() == 1)
		return pArrRet;



	// sort other arrays in the same way
	SymbolArray *pArrArr = new SymbolArray();
	pArrArr->GetArr().reserve(vecSyms.size());
	pArrArr->GetArr().push_back(pArrRet);

	for(unsigned int iElem=1; iElem<vecSyms.size(); ++iElem)
	{
		const Symbol* pSym = vecSyms[iElem];
		if(pSym->GetType() != SYMBOL_ARRAY)
		{
			tl::log_err(linenr(runinfo), "Arguments to sort have to be arrays. Ignoring.");
			continue;
		}

		if(((SymbolArray*)pSym)->GetArr().size() != vecSortedIndices.size())
		{
			tl::log_err(linenr(runinfo), "Array size mismatch in sort. Ignoring.");
			continue;
		}

		SymbolArray *pNextArr = (SymbolArray*)pSym->clone();
		_rearrangearr(pNextArr->GetArr(), vecSortedIndices);

		pArrArr->GetArr().push_back(pNextArr);
	}

	pArrArr->UpdateIndices();
	return pArrArr;
}

static Symbol* fkt_sort(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	return _fkt_sort(vecSyms, info, runinfo, pSymTab, 0);
}

static Symbol* fkt_sort_rev(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	return _fkt_sort(vecSyms, info, runinfo, pSymTab, 1);
}

static Symbol* fkt_zip(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	SymbolArray* pArrRet = new SymbolArray();
	std::vector<Symbol*>* pVec = &pArrRet->GetArr();

	bool bFirst = 1;
	unsigned int iSize = 0;
	for(Symbol* pSym : vecSyms)
	{
		if(pSym->GetType() != SYMBOL_ARRAY)
		{
			tl::log_err(linenr(runinfo), "Symbol \"", pSym->GetName(),
				"\" is of type ", pSym->GetTypeName(),
				", but zip needs vectors. Ignoring.");
			continue;
		}

		std::vector<Symbol*>& curSym = ((SymbolArray*)pSym)->GetArr();
		if(bFirst)
		{
			iSize = curSym.size();
			pVec->reserve(iSize);
			for(unsigned int iCurSym=0; iCurSym<iSize; ++iCurSym)
				pVec->push_back(new SymbolArray());

			bFirst = 0;
		}

		if(curSym.size() < iSize)
		{
			for(unsigned int i=0; i<iSize-curSym.size(); ++i)
			{
				delete *pVec->rbegin();
				pVec->pop_back();
			}

			iSize = curSym.size();
		}

		for(unsigned int i=0; i<iSize; ++i)
			((SymbolArray*)(*pVec)[i])->GetArr().push_back(curSym[i]->clone());
	}

	pArrRet->UpdateIndices();
	return pArrRet;
}
// --------------------------------------------------------------------------------



// --------------------------------------------------------------------------------
// string operations
static Symbol* fkt_trim(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_ANY}, {0}, "trim"))
		return 0;

	if(vecSyms[0]->GetType() == SYMBOL_ARRAY)
	{
		SymbolArray *pArr = new SymbolArray();
		pArr->GetArr().reserve(((SymbolArray*)vecSyms[0])->GetArr().size());

		for(Symbol* pSymArr : ((SymbolArray*)vecSyms[0])->GetArr())
		{
			std::vector<Symbol*> vecDummy = { pSymArr };
			pArr->GetArr().push_back(fkt_trim(vecDummy, info, runinfo, pSymTab));
		}

		pArr->UpdateIndices();
		return pArr;
	}
	else if(vecSyms[0]->GetType() == SYMBOL_STRING)
	{
		t_string str = ((SymbolString*)vecSyms[0])->GetVal();
		tl::trim(str);
		return new SymbolString(str);
	}

	// simply copy non-string arguments
	//G_CERR << linenr(T_STR"Warning", info)
	//		<< "Called trim with invalid argument." << std::endl;
	return vecSyms[0]->clone();
}

// split(string, delim)
static Symbol* fkt_split(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_STRING, SYMBOL_STRING}, {0, 1}, "split"))
		return 0;

	static const std::string strDefaultDelim = " ";

	const t_string& str = ((SymbolString*)vecSyms[0])->GetVal();
	const t_string* pstrDelim;

	if(vecSyms.size() >= 2)
		pstrDelim = &((SymbolString*)vecSyms[1])->GetVal();
	else
		pstrDelim = &strDefaultDelim;

	std::size_t iPos = str.find(*pstrDelim);

	t_string str0 = str.substr(0, iPos);
	t_string str1;

	if(iPos != t_string::npos)
		str1 = str.substr(iPos + pstrDelim->length(), t_string::npos);

	SymbolArray* pSymRet = new SymbolArray();
	pSymRet->GetArr().push_back(new SymbolString(str0));
	pSymRet->GetArr().push_back(new SymbolString(str1));

	pSymRet->UpdateIndices();
	return pSymRet;
}

// tokens(string, delim)
static Symbol* fkt_tokens(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	t_string *pstrInput = 0;
	t_string strDelim = T_STR" \t\n";

	// split("Test 123")
	if(vecSyms.size() >= 1 && vecSyms[0]->GetType()==SYMBOL_STRING)
		pstrInput = &((SymbolString*)vecSyms[0])->GetVal();

	// split("Test 123", " \t\n")
	if(vecSyms.size() >= 2 && vecSyms[1]->GetType()==SYMBOL_STRING)
		strDelim = ((SymbolString*)vecSyms[1])->GetVal();

	if(!pstrInput)
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "Called tokens with invalid arguments." << std::endl;
		throw tl::Err(ostrErr.str(),0);
	}

	std::vector<t_string> vecTokens;
	tl::get_tokens<t_string>(*pstrInput, strDelim, vecTokens);

	SymbolArray* pArr = new SymbolArray;
	for(const t_string& strTok : vecTokens)
		pArr->GetArr().push_back(new SymbolString(strTok));

	pArr->UpdateIndices();
	return pArr;
}


static Symbol* fkt_replace(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_STRING, SYMBOL_STRING, SYMBOL_STRING}, {0,0,0}, "replace"))
		return 0;

	t_string str = ((SymbolString*)vecSyms[0])->GetVal();
	const t_string& strOld = ((SymbolString*)vecSyms[1])->GetVal();
	const t_string& strNew = ((SymbolString*)vecSyms[2])->GetVal();

	tl::find_all_and_replace(str, strOld, strNew);

	return new SymbolString(str);
}


// replace_regex(str, regex, repl);
static Symbol* fkt_replace_regex(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_STRING, SYMBOL_STRING, SYMBOL_STRING}, {0,0,0}, "regex_replace"))
		return 0;

	const t_string& str = ((SymbolString*)vecSyms[0])->GetVal();
	const t_string& strRegex = ((SymbolString*)vecSyms[1])->GetVal();
	const t_string& strRepl = ((SymbolString*)vecSyms[2])->GetVal();

	t_string strRet;
	try
	{
		rex::basic_regex<t_char> rex(strRegex, rex::regex::ECMAScript);
		strRet = rex::regex_replace(str, rex, strRepl);
	}
	catch(const std::exception& ex)
	{
		tl::log_err(linenr(runinfo), "Regex evaluation failed with error: ", ex.what());
		return 0;
	}

	return new SymbolString(std::move(strRet));
}

static Symbol* fkt_find_regex(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_STRING, SYMBOL_STRING}, {0,0}, "regex_find"))
		return 0;

	const t_string& strOrg = ((SymbolString*)vecSyms[0])->GetVal();
	const t_string& strRegex = ((SymbolString*)vecSyms[1])->GetVal();

	struct Match
	{
		t_string strMatch;
		int iPos;
		int iLen;
	};
	std::vector<Match> vecMatches;

	try
	{
		t_string str = strOrg;
		rex::basic_regex<t_char> rex(strRegex, rex::regex::ECMAScript);
		rex::smatch mRes;

		int iPos = 0;
		while(1)
		{
			if(!rex::regex_search(str, mRes, rex) || mRes.empty())
				break;
			if(!mRes.length())
				break;

			Match match;
			match.strMatch = mRes[0];
			match.iPos = mRes.position() + iPos;
			match.iLen = mRes.length();

			iPos = match.iPos + match.iLen;
			vecMatches.push_back(std::move(match));

			str = mRes.suffix().str();
		}
	}
	catch(const std::exception& ex)
	{
		tl::log_err(linenr(runinfo), "Regex evaluation failed with error: ", ex.what());
		return 0;
	}


	SymbolArray *pRet = new SymbolArray();
	pRet->GetArr().reserve(vecMatches.size());

	for(const Match& match : vecMatches)
	{
		SymbolArray *pMatch = new SymbolArray();
		pRet->GetArr().push_back(pMatch);

		pMatch->GetArr().reserve(3);
		pMatch->GetArr().push_back(new SymbolString(match.strMatch));
		pMatch->GetArr().push_back(new SymbolInt(match.iPos));
		pMatch->GetArr().push_back(new SymbolInt(match.iLen));
	}

	pRet->UpdateIndices();
	return pRet;
}

static Symbol* fkt_match_regex(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_STRING, SYMBOL_STRING}, {0,0}, "regex_match"))
		return 0;

	const t_string& str = ((SymbolString*)vecSyms[0])->GetVal();
	const t_string& strRegex = ((SymbolString*)vecSyms[1])->GetVal();
	bool bMatch = 0;

	try
	{
		rex::basic_regex<t_char> rex(strRegex, rex::regex::ECMAScript);
		rex::smatch mRes;
		bMatch = rex::regex_match(str, mRes, rex);
	}
	catch(const std::exception& ex)
	{
		tl::log_err(linenr(runinfo), "Regex evaluation failed with error: ", ex.what());
		return 0;
	}

	return new SymbolInt(bMatch);
}

static Symbol* fkt_subfind_regex(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_STRING, SYMBOL_STRING, SYMBOL_INT}, {0,0,1}, "regex_subfind"))
		return 0;

	const t_string& str = ((SymbolString*)vecSyms[0])->GetVal();
	const t_string& strRegex = ((SymbolString*)vecSyms[1])->GetVal();

	// search or match regex?
	bool bSearch = 1;

	SymbolArray *pSymRet = new SymbolArray();

	try
	{
		rex::basic_regex<t_char> rex(strRegex, rex::regex::ECMAScript);
		rex::smatch mRes;

		if(bSearch)
			rex::regex_search(str, mRes, rex);
		else
			rex::regex_match(str, mRes, rex);

		for(const rex::ssub_match& strMatch : mRes)
		{
			pSymRet->GetArr().push_back(new SymbolString(strMatch));
		}
	}
	catch(const std::exception& ex)
	{
		tl::log_err(linenr(runinfo), "Regex evaluation failed with error: ", ex.what());
		return 0;
	}

	pSymRet->UpdateIndices();
	return pSymRet;
}
// --------------------------------------------------------------------------------




// --------------------------------------------------------------------------------
// map operations

static Symbol* fkt_has_key(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_MAP, SYMBOL_ANY}, {0,0}, "has_key"))
		return 0;

	const SymbolMap* pMap = (SymbolMap*)vecSyms[0];

	int bHasKey = (pMap->GetMap().find(vecSyms[1]) != pMap->GetMap().end());
	return new SymbolInt(bHasKey);
}

static Symbol* fkt_rm_key(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_MAP, SYMBOL_ANY}, {0,0}, "rm_key"))
		return 0;

	SymbolMap* pMap = (SymbolMap*)vecSyms[0];
	SymbolMap::t_map::iterator iter = pMap->GetMap().find(vecSyms[1]);

	Symbol *pErasedSym = nullptr;
	bool bHasKey = (iter != pMap->GetMap().end());
	if(bHasKey)
	{
		pErasedSym = iter->second;
		pMap->GetMap().erase(iter);
	}

	return pErasedSym;
}

// --------------------------------------------------------------------------------



extern void init_ext_basic_calls()
{
	t_mapFkts mapFkts =
	{
		// basic stuff
		t_mapFkts::value_type(T_STR"interp_ver", fkt_version),
		t_mapFkts::value_type(T_STR"register_var", fkt_register_var),
		t_mapFkts::value_type(T_STR"is_valid", fkt_is_valid),
		t_mapFkts::value_type(T_STR"null", fkt_null),
		t_mapFkts::value_type(T_STR"call", fkt_call),
		t_mapFkts::value_type(T_STR"get_argnames", fkt_getargnames),
		t_mapFkts::value_type(T_STR"is_rval", fkt_is_rval),
		t_mapFkts::value_type(T_STR"traceback", fkt_traceback),

		// input/output
		t_mapFkts::value_type(T_STR"input", fkt_input),
		t_mapFkts::value_type(T_STR"output", fkt_output),
		t_mapFkts::value_type(T_STR"print", fkt_print),	// output with "\n" at the end
		t_mapFkts::value_type(T_STR"sleep", fkt_sleep),

		// modules
		t_mapFkts::value_type(T_STR"import", fkt_import),

		// symbols & casts
		t_mapFkts::value_type(T_STR"int", fkt_int),
		t_mapFkts::value_type(T_STR"real", fkt_double),
		t_mapFkts::value_type(T_STR"real_vec", fkt_double_vec),
		t_mapFkts::value_type(T_STR"str", fkt_str),
		t_mapFkts::value_type(T_STR"map", fkt_map),
		t_mapFkts::value_type(T_STR"vec", fkt_array),
		t_mapFkts::value_type(T_STR"complex", fkt_complex),
		t_mapFkts::value_type(T_STR"complex_polar", fkt_complex_polar),
		t_mapFkts::value_type(T_STR"imag", fkt_cplx_imag),
		t_mapFkts::value_type(T_STR"carg", fkt_cplx_arg),
		t_mapFkts::value_type(T_STR"has_var", fkt_has_var),
		t_mapFkts::value_type(T_STR"typeof", fkt_typeof),
		t_mapFkts::value_type(T_STR"set_prec", fkt_setprec),

		// string operations
		t_mapFkts::value_type(T_STR"trim", fkt_trim),
		t_mapFkts::value_type(T_STR"split", fkt_split),
		t_mapFkts::value_type(T_STR"tokens", fkt_tokens),
		t_mapFkts::value_type(T_STR"replace", fkt_replace),
		t_mapFkts::value_type(T_STR"length", fkt_array_size),

		t_mapFkts::value_type(T_STR"regex_replace", fkt_replace_regex),
		t_mapFkts::value_type(T_STR"regex_find", fkt_find_regex),
		t_mapFkts::value_type(T_STR"regex_match", fkt_match_regex),
		t_mapFkts::value_type(T_STR"regex_subfind", fkt_subfind_regex),


		// array operations
		t_mapFkts::value_type(T_STR"vec_size", fkt_array_size),	// deprecated, use "length" instead
		t_mapFkts::value_type(T_STR"cur_iter", fkt_cur_iter),
		t_mapFkts::value_type(T_STR"zip", fkt_zip),
		t_mapFkts::value_type(T_STR"sort", fkt_sort),
		t_mapFkts::value_type(T_STR"sort_rev", fkt_sort_rev),
		t_mapFkts::value_type(T_STR"splice", fkt_splice),
		t_mapFkts::value_type(T_STR"append", fkt_append),

		// map/array operations
		t_mapFkts::value_type(T_STR"contains", fkt_contains),
		t_mapFkts::value_type(T_STR"has_key", fkt_has_key),
		t_mapFkts::value_type(T_STR"rm_key", fkt_rm_key),
		t_mapFkts::value_type(T_STR"find", fkt_find),
	};

	add_ext_calls(mapFkts);
}
