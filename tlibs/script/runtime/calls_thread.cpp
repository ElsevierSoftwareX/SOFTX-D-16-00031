/*
 * external thread functions
 * @author tweber
 * @date dec 2013
 * @license GPLv2 or GPLv3
 */

#include "lang/types.h"
#include "helper/flags.h"
#include "string/string.h"
#include "log/log.h"
#include "calls_thread.h"
#include "lang/calls.h"
#include <thread>
#include <future>
#include <wait.h>
#include <cstdlib>

static inline Symbol* fkt_exec(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	t_string strExec;

	for(Symbol *pSym : vecSyms)
		if(pSym)
		{
			strExec += pSym->print();
			strExec += T_STR" ";
		}

	bool bOk = 0;

	std::string _strExec = WSTR_TO_STR(strExec);
	FILE *pPipe = (FILE*)tl::my_popen(_strExec.c_str(), "w");

	if(pPipe)
	{
		bOk = 1;
		int iRet = tl::my_pclose(pPipe);
		if(iRet == -1)
		{
			bOk = 0;
		}
		else
		{
			int iExitCode = int(char(WEXITSTATUS(iRet)));
			bOk = (iExitCode==0);
		}
	}

	return new SymbolInt(bOk);
}

static inline Symbol* fkt_exit(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	int iStatus = 0;
	if(vecSyms.size()>=1 && vecSyms[0])
		iStatus = vecSyms[0]->GetValInt();

	std::exit(iStatus);		// TODO: change to less brutal way to exit
	return 0;
}

// --------------------------------------------------------------------------------
// thread

std::vector<Symbol*>* clone_symbols(const std::vector<Symbol*>* pvecSyms,
	unsigned int iBegin=0)
{
	if(!pvecSyms)
		return 0;

	std::vector<Symbol*> *pvec = new std::vector<Symbol*>;
	pvec->reserve(pvecSyms->size());

	for(unsigned int i=iBegin; i<pvecSyms->size(); ++i)
	{
		Symbol *pSym = (*pvecSyms)[i];
		if(pSym) pSym = pSym->clone();
		pvec->push_back(pSym);
	}

	return pvec;
}

void delete_symbols(std::vector<Symbol*>* pvecSyms)
{
	for(Symbol *pSym : *pvecSyms)
		delete pSym;
	delete pvecSyms;
}

static Symbol* task_proc(NodeFunction* pFunc, ParseInfo* pinfo, std::vector<Symbol*>* pvecSyms)
{
	if(!pFunc || !pinfo) return 0;

	SymbolTable *pTable = new SymbolTable();
	SymbolArray arrArgs;
	arrArgs.SetDontDel(1);

	const NodeFunction *pThreadFunc = 0;
	RuntimeInfo *pruninfo2 = new RuntimeInfo();
	pThreadFunc = (NodeFunction*)pFunc/*->clone()->optimize()*/;

	if(pvecSyms) arrArgs.GetArr() = *pvecSyms;
	arrArgs.UpdateIndices();

	pTable->InsertSymbol(T_STR"<args>", &arrArgs);
	Symbol* pRet = pThreadFunc->eval(*pinfo, *pruninfo2, pTable);
	pTable->RemoveSymbolNoDelete(T_STR"<ret>");
	pTable->RemoveSymbolNoDelete(T_STR"<args>");

	if(pTable) delete pTable;
	if(pvecSyms) delete_symbols(pvecSyms);
	//if(pThreadFunc) delete pThreadFunc;
	if(pruninfo2) delete pruninfo2;

	return pRet;
}

static void thread_proc(NodeFunction* pFunc, ParseInfo* pinfo, std::vector<Symbol*>* pvecSyms)
{
	// ignore return value
	Symbol *pRet = task_proc(pFunc, pinfo, pvecSyms);
	safe_delete(pRet, 0, pinfo);
}

static Symbol* fkt_thread_task(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info,
						RuntimeInfo &runinfo, 
						SymbolTable* pSymTab,
						bool bTask=0)
{
	if(vecSyms.size()<1)
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "Need thread proc identifier." << std::endl;
		throw tl::Err(ostrErr.str(), 0);
	}

	Symbol* _pSymIdent = vecSyms[0];
	if(_pSymIdent->GetType() != SYMBOL_STRING)
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "Thread proc identifier needs to be a string." << std::endl;
		throw tl::Err(ostrErr.str(), 0);
	}

	SymbolString *pSymIdent = (SymbolString*)_pSymIdent;
	const t_string& strIdent = pSymIdent->GetVal();


	NodeFunction* pFunc = info.GetFunction(strIdent);
	if(pFunc == 0)
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "Thread proc \"" << strIdent << "\" not defined." << std::endl;
		throw tl::Err(ostrErr.str(), 0);
	}

	std::vector<Symbol*>* vecThreadSymsClone = clone_symbols(&vecSyms, 1);

	t_int iHandle = -1;
	if(bTask)
	{
		bool bIsThread = 1;
		std::launch policy = std::launch::async /*| std::launch::deferred*/;
		unsigned int iNumThreads = info.phandles->CountAllThreads();

		// start deferred
		if(iNumThreads >= std::thread::hardware_concurrency())
		{
			bIsThread = 0;
			// let system decide
			policy |= std::launch::deferred;
		}

		std::future<Symbol*> *pFuture = new std::future<Symbol*>(
			std::async(policy, ::task_proc, pFunc, &info, vecThreadSymsClone));
		iHandle = info.phandles->AddHandle(new HandleTask(pFuture, bIsThread));
	}
	else
	{
		std::thread* pThread = new std::thread(::thread_proc, pFunc, &info, vecThreadSymsClone);
		iHandle = info.phandles->AddHandle(new HandleThread(pThread));
	}

	return new SymbolInt(iHandle);
}

static Symbol* fkt_thread(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	return fkt_thread_task(vecSyms, info, runinfo, pSymTab, 0);
}

static Symbol* fkt_task(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	return fkt_thread_task(vecSyms, info, runinfo, pSymTab, 1);
}

static Symbol* fkt_thread_hwcount(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	unsigned int iNumThreads = std::thread::hardware_concurrency();
	if(iNumThreads == 0)
		iNumThreads = 1;

	return new SymbolInt(t_int(iNumThreads));
}

static Symbol* fkt_mutex(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	std::mutex* pMutex = new std::mutex;
	t_int iHandle = info.phandles->AddHandle(new HandleMutex(pMutex));

	return new SymbolInt(iHandle);
}

static Symbol* fkt_begin_critical(const std::vector<Symbol*>& vecSyms,
								ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	// no argument given: lock global mutex
	if(vecSyms.size() == 0)
	{
		info.pmutexGlobal->lock();
	}
	else
	{
		for(Symbol* pSym : vecSyms)
		{
			if(pSym->GetType() == SYMBOL_ARRAY)
				fkt_begin_critical(((SymbolArray*)pSym)->GetArr(), info, runinfo, pSymTab);
			else if(pSym->GetType() == SYMBOL_INT)
			{
				//std::lock_guard<std::mutex> lck(*info.pmutexCrit);

				int iHandle = pSym->GetValInt();
				Handle *pHandle = info.phandles->GetHandle(iHandle);

				if(pHandle==0 || pHandle->GetType()!=HANDLE_MUTEX)
				{
					std::ostringstream ostrErr;
					ostrErr << linenr(runinfo) << "Handle (" 
						<< iHandle << ") does not exist"
						<< " or is not a mutex handle. Ignoring." 
						<< std::endl;
					throw tl::Err(ostrErr.str(), 0);
					//continue;
				}

				std::mutex *pMutex = ((HandleMutex*)pHandle)->GetInternalHandle();
				pMutex->lock();
			}
			else
			{
				std::ostringstream ostrErr;
				ostrErr << linenr(runinfo) << "Invalid mutex handle: "
						<< pSym->print() << " Ignoring." << std::endl;
				throw tl::Err(ostrErr.str(), 0);
			}
		}
	}

	return 0;
}

static Symbol* fkt_end_critical(const std::vector<Symbol*>& vecSyms,
								ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	// no argument given: unlock global mutex
	if(vecSyms.size() == 0)
	{
		info.pmutexGlobal->unlock();
	}
	else
	{
		for(Symbol* pSym : vecSyms)
		{
			if(pSym->GetType() == SYMBOL_ARRAY)
				fkt_begin_critical(((SymbolArray*)pSym)->GetArr(), info, runinfo, pSymTab);
			else if(pSym->GetType() == SYMBOL_INT)
			{
				//std::lock_guard<std::mutex> lck(*info.pmutexCrit);

				int iHandle = pSym->GetValInt();
				Handle *pHandle = info.phandles->GetHandle(iHandle);

				if(pHandle==0 || pHandle->GetType()!=HANDLE_MUTEX)
				{
					std::ostringstream ostrErr;
					ostrErr << linenr(runinfo) << "Handle (" 
						<< iHandle << ") does not exist"
						<< " or is not a mutex handle." << std::endl;
					throw tl::Err(ostrErr.str(), 0);
					//continue;
				}

				std::mutex *pMutex = ((HandleMutex*)pHandle)->GetInternalHandle();
				pMutex->unlock();
			}
			else
			{
				std::ostringstream ostrErr;
				ostrErr << linenr(runinfo) << "Invalid mutex handle: "
						<< pSym->print() << std::endl;
				throw tl::Err(ostrErr.str(), 0);
			}
		}
	}

	return 0;
}


static Symbol* fkt_thread_join(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(vecSyms.size()<1)
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "join needs at least one argument." << std::endl;
		throw tl::Err(ostrErr.str(),0);
	}

	Symbol *pRet = 0;
	for(Symbol* pSym : vecSyms)
	{
		if(pSym == 0) continue;

		if(pSym->GetType() == SYMBOL_ARRAY)
		{
			return fkt_thread_join(((SymbolArray*)pSym)->GetArr(), info, runinfo, pSymTab);
		}

		if(pSym->GetType() != SYMBOL_INT)
		{
			std::ostringstream ostrErr;
			ostrErr << linenr(runinfo) << "join needs thread handles." << std::endl;
			throw tl::Err(ostrErr.str(), 0);
			//continue;
		}

		t_int iHandle = ((SymbolInt*)pSym)->GetVal();
		Handle *pHandle = info.phandles->GetHandle(iHandle);

		if(pHandle==0)
		{
			std::ostringstream ostrErr;
			ostrErr << linenr(runinfo) << "Handle (" 
				<< iHandle << ") does not exist." << std::endl;
			throw tl::Err(ostrErr.str(), 0);
			//continue;
		}

		if(pHandle->GetType() == HANDLE_THREAD)
		{
			HandleThread *pThreadHandle = (HandleThread*)pHandle;
			std::thread *pThread = pThreadHandle->GetInternalHandle();

			pThread->join();
		}
		else if(pHandle->GetType() == HANDLE_TASK)
		{
			HandleTask *pTaskHandle = (HandleTask*)pHandle;
			std::future<Symbol*> *pFut = pTaskHandle->GetInternalHandle();

			// TODO: array task joins
			if(pRet != 0)
				delete pRet;
			pRet = pFut->get();
		}
		else
		{
			std::ostringstream ostrErr;
			ostrErr << linenr(runinfo) << "Handle (" 
				<< iHandle << ") is invalid." << std::endl;
			throw tl::Err(ostrErr.str(), 0);
		}
		
		info.phandles->CloseHandle(iHandle);
	}

	return pRet;
}


// --------------------------------------------------------------------------------
// nthread

// nthread(iNumThreads, strFunc, vecArgs, ...)
static Symbol* fkt_nthread(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(vecSyms.size()<3)
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) 
			<< "nthread needs at least 3 arguments: N, func, arg." 
			<< std::endl;
		throw tl::Err(ostrErr.str(),0);
	}

	Symbol* _pSymN = vecSyms[0];
	if(_pSymN->GetType() != SYMBOL_INT)
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "Number of threads has to be integer." << std::endl;
		throw tl::Err(ostrErr.str(),0);
	}

	SymbolInt *pSymN = (SymbolInt*)_pSymN;
	t_int iNumThreads = pSymN->GetVal();



	Symbol* _pSymIdent = vecSyms[1];
	if(_pSymIdent->GetType() != SYMBOL_STRING)
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "Thread proc identifier needs to be a string." << std::endl;
		throw tl::Err(ostrErr.str(), 0);
	}

	SymbolString *pSymIdent = (SymbolString*)_pSymIdent;
	const t_string& strIdent = pSymIdent->GetVal();



	Symbol* _pSymArr = vecSyms[2];
	if(_pSymArr->GetType() != SYMBOL_ARRAY)
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "Thread arg has to be an array." << std::endl;
		throw tl::Err(ostrErr.str(), 0);
	}

	SymbolArray *pSymArr = (SymbolArray*)_pSymArr;
	const std::vector<Symbol*>& vecArr = pSymArr->GetArr();



	NodeFunction* pFunc = info.GetFunction(strIdent);
	if(pFunc == 0)
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "Thread proc \"" << strIdent 
			<< "\" not defined." << std::endl;
		throw tl::Err(ostrErr.str(), 0);
	}





	if(iNumThreads > int(vecArr.size()))
	{
		iNumThreads = vecArr.size();
		tl::log_warn(linenr(runinfo), "More threads requested in nthread than necessary, ",
					"reducing to array size (", iNumThreads, ").");
	}


	std::vector<SymbolArray*> vecSymArrays;
	vecSymArrays.resize(iNumThreads);

	t_int iCurTh = 0;
	for(Symbol* pThisSym : vecArr)
	{
		if(!vecSymArrays[iCurTh])
			vecSymArrays[iCurTh] = new SymbolArray();

		vecSymArrays[iCurTh]->GetArr().push_back(pThisSym->clone());
		vecSymArrays[iCurTh]->UpdateLastNIndices(1);

		++iCurTh;
		if(iCurTh == iNumThreads)
			iCurTh = 0;
	}



	std::vector<std::thread*> vecThreads;
	vecThreads.reserve(iNumThreads);

	for(iCurTh=0; iCurTh<iNumThreads; ++iCurTh)
	{
		std::vector<Symbol*>* vecThreadSyms = new std::vector<Symbol*>;
		vecThreadSyms->reserve(vecSyms.size()-3+1);

		vecThreadSyms->push_back(vecSymArrays[iCurTh]);

		for(unsigned int iSym=3; iSym<vecSyms.size(); ++iSym)
			vecThreadSyms->push_back(vecSyms[iSym]->clone());

		std::thread *pth = new std::thread(::thread_proc, pFunc, &info, vecThreadSyms);
		vecThreads.push_back(pth);
	}

	/*
	// automatically join
	for(iCurTh=0; iCurTh<iNumThreads; ++iCurTh)
	{
		vecThreads[iCurTh]->join();
		delete vecThreads[iCurTh];
		vecThreads[iCurTh] = 0;
	}*/


	SymbolArray* pArrThreads = new SymbolArray();

	for(iCurTh=0; iCurTh<iNumThreads; ++iCurTh)
	{
		std::thread* pCurThread = vecThreads[iCurTh];
		t_int iHandle = info.phandles->AddHandle(new HandleThread(pCurThread));
		SymbolInt *pSymThreadHandle = new SymbolInt(iHandle);

		pArrThreads->GetArr().push_back(pSymThreadHandle);
	}

	pArrThreads->UpdateIndices();
	return pArrThreads;
}

// --------------------------------------------------------------------------------


extern void init_ext_thread_calls()
{
	t_mapFkts mapFkts =
	{
		// threads & tasks
		t_mapFkts::value_type(T_STR"thread", fkt_thread),
		t_mapFkts::value_type(T_STR"nthread", fkt_nthread),
		t_mapFkts::value_type(T_STR"task", fkt_task),

		t_mapFkts::value_type(T_STR"thread_hwcount", fkt_thread_hwcount),
		t_mapFkts::value_type(T_STR"join", fkt_thread_join),
		t_mapFkts::value_type(T_STR"mutex", fkt_mutex),
		t_mapFkts::value_type(T_STR"begin_critical", fkt_begin_critical),
		t_mapFkts::value_type(T_STR"end_critical", fkt_end_critical),

		// processes
		t_mapFkts::value_type(T_STR"exec", fkt_exec),
		t_mapFkts::value_type(T_STR"exit", fkt_exit),
	};

	add_ext_calls(mapFkts);
}
