/*
 * Script interpreter
 * Node evaluation
 * @author tweber
 * @date 10 oct 2013
 * @license GPLv2 or GPLv3
 */

#include "node.h"
#include "info.h"
#include "calls.h"
#include "log/log.h"

Symbol* NodeReturn::eval(ParseInfo &info, RuntimeInfo& runinfo, SymbolTable *pSym) const
{
	if(runinfo.IsExecDisabled()) return 0;
	Symbol *pRet = 0;

	if(m_pExpr)
	{
		Symbol *pEval = m_pExpr->eval(info, runinfo, pSym);
		if(clone_if_needed(pEval, pRet))
			safe_delete(pEval, pSym, &info);
	}
	pSym->InsertSymbol(T_STR"<ret>", pRet ? pRet : 0);

	runinfo.bWantReturn = 1;
	return pRet;
}

Symbol* NodeBreak::eval(ParseInfo &info, RuntimeInfo& runinfo, SymbolTable *pSym) const
{
	if(runinfo.IsExecDisabled()) return 0;

	if(runinfo.pCurLoop==0)
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo)
				<< "Cannot use break outside loop."
				<< std::endl;
		throw tl::Err(ostrErr.str(),0);
	}

	runinfo.bWantBreak = 1;
	return 0;
}

Symbol* NodeContinue::eval(ParseInfo &info, RuntimeInfo& runinfo, SymbolTable *pSym) const
{
	if(runinfo.IsExecDisabled()) return 0;

	if(runinfo.pCurLoop==0)
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo)
				<< "Cannot use continue outside loop."
				<< std::endl;
		throw tl::Err(ostrErr.str(),0);
	}

	runinfo.bWantContinue = 1;
	return 0;
}

Symbol* NodeIdent::eval(ParseInfo &info, RuntimeInfo& runinfo, SymbolTable *pSym) const
{
	if(runinfo.IsExecDisabled()) return 0;

	// local symbol
	Symbol *pSymbol = 0;
	if(pSym)
		pSymbol = pSym->GetSymbol(m_strIdent);

	// global symbol
	Symbol *pSymbolGlob = 0;
	if(info.pGlobalSyms)
	{
		std::lock_guard<std::mutex> _lck(*info.pmutexGlobalSyms);
		pSymbolGlob = info.pGlobalSyms->GetSymbol(m_strIdent);
	}

	if(pSymbol && pSymbolGlob)
	{
		tl::log_warn(linenr(runinfo), "Symbol \"", m_strIdent,
			"\" exists in local and global scope, using local one.");
	}

	// if no local symbol is available, use global symbol instead
	if(!pSymbol)
		pSymbol = pSymbolGlob;

	if(!pSymbol)
	{
		// no throw, so null symbol or symbol validity query can be used in script
		tl::log_err(linenr(runinfo), "Symbol \"", m_strIdent,
			"\" not in symbol table.");
		return 0;
	}


	if(pSymbol == pSymbolGlob)
		info.pmutexGlobalSyms->lock();

	pSymbol->SetIdent(m_strIdent);

	if(pSymbol == pSymbolGlob)
		info.pmutexGlobalSyms->unlock();
	return pSymbol;
}

Symbol* NodeCall::eval(ParseInfo &info, RuntimeInfo& runinfo, SymbolTable *pSym) const
{
	if(runinfo.IsExecDisabled()) return 0;
	runinfo.pCurCaller = this;

	if(m_pIdent->GetType() != NODE_IDENT)
		return 0;
	//if(m_pArgs->GetType() != NODE_ARGS)
	//	return 0;

	const NodeIdent* pIdent = (NodeIdent*) m_pIdent;
	const t_string& strFkt = pIdent->GetIdent();


	bool bCallUserFkt = 0;
	// user-defined function
	NodeFunction *pFkt = info.GetFunction(strFkt);;
	if(pFkt)
		bCallUserFkt = 1;


	SymbolArray arrArgs;
	arrArgs.SetDontDel(1);
	std::vector<Symbol*> &vecArgSyms = arrArgs.GetArr();
	for(const Node* pNode : m_vecArgs)
	{
		if(pNode->GetType() == NODE_UNPACK)
		{
			Node *pChild = ((NodeUnaryOp*)pNode)->GetChild();
			if(!pChild)
			{
				std::ostringstream ostrErr;
				ostrErr << linenr(runinfo)
					<< "Invalid symbol to unpack." << std::endl;
				throw tl::Err(ostrErr.str(),0);
			}
			Symbol *pSymChild = pChild->eval(info, runinfo, pSym);

			if(pSymChild->GetType() == SYMBOL_MAP)
			{
				if(!pFkt)
				{
					std::ostringstream ostrErr;
					ostrErr << linenr(runinfo)
						<< "Tried to unpack map for (external) function with unnamed parameters." << std::endl;
					throw tl::Err(ostrErr.str(),0);
				}

				SymbolMap::t_map& mapSyms = ((SymbolMap*)pSymChild)->GetMap();
				std::vector<t_string> vecParamNames;
				if(pFkt) vecParamNames = pFkt->GetParamNames();

				for(unsigned int iParam=vecArgSyms.size(); iParam<vecParamNames.size(); ++iParam)
				{
					const t_string& strParamName = vecParamNames[iParam];
					SymbolMap::t_map::iterator iter = mapSyms.find(SymbolMapKey(strParamName));
					if(iter == mapSyms.end())
					{
						tl::log_err(linenr(runinfo), "Parameter \"", strParamName,
							"\" not in argument map. Using 0.");

						vecArgSyms.push_back(new SymbolReal(0.));
						continue;
					}

					Symbol *pSymClone = nullptr;
					clone_if_needed(iter->second, pSymClone);
					vecArgSyms.push_back(pSymClone);
				}
			}
			else if(pSymChild->GetType() == SYMBOL_ARRAY)
			{
				SymbolArray::t_arr& arrSyms = ((SymbolArray*)pSymChild)->GetArr();
				std::vector<t_string> vecParamNames;
				if(pFkt) vecParamNames = pFkt->GetParamNames();

				for(unsigned int iParam=0; iParam<arrSyms.size(); ++iParam)
				{
					if(pFkt && vecArgSyms.size()>=vecParamNames.size())
					{
						tl::log_err(linenr(runinfo), "Exceeded call parameter size in vector-unpack.");
						break;
					}

					Symbol *pSymClone = nullptr;
					clone_if_needed(arrSyms[iParam], pSymClone);
					vecArgSyms.push_back(pSymClone);
				}
			}
			else
			{
				std::ostringstream ostrErr;
				ostrErr << linenr(runinfo)
					<< "Unpack only valid for maps or vectors." << std::endl;
				throw tl::Err(ostrErr.str(),0);
			}

			safe_delete(pSymChild, pSym, &info);
		}
		else
		{
			//tl::log_debug("type: ", pNode->GetType());
			Symbol *pSymbol = pNode->eval(info, runinfo, pSym);
			vecArgSyms.push_back(pSymbol);
		}
	}

	arrArgs.UpdateIndices(false);
	Symbol* pFktRet = 0;
	if(bCallUserFkt)	// call user-defined function
	{
		//pFkt->SetArgSyms(&vecArgSyms);
		pSym->InsertSymbol(T_STR"<args>", &arrArgs);
		pFktRet = pFkt->eval(info, runinfo, pSym);
		pSym->RemoveSymbolNoDelete(T_STR"<args>");
		runinfo.bWantReturn = 0;
	}
	else if(has_ext_call(strFkt))		// call system function
	{
		if(info.bEnableDebug)
		{
			std::string strTrace = "syscall: " + strFkt + ", " 
					+ std::to_string(vecArgSyms.size()) + " args";
			info.PushTraceback(std::move(strTrace));
		}

		pFktRet = ext_call(strFkt, vecArgSyms, info, runinfo, pSym);

		if(info.bEnableDebug)
		{
			info.PopTraceback();
		}
	}
	else
	{
		tl::log_err(linenr(runinfo), "Tried to call undefined function \"", strFkt, "\".");
	}

	arrArgs.ClearIndices();

	for(Symbol *pArgSym : vecArgSyms)
		safe_delete(pArgSym, pSym, &info);
	return pFktRet;
}

Symbol* NodeReal::eval(ParseInfo &info, RuntimeInfo& runinfo, SymbolTable *pSym) const
{
	if(runinfo.IsExecDisabled()) return 0;
	return m_pSymbol;
}

Symbol* NodeInt::eval(ParseInfo &info, RuntimeInfo& runinfo, SymbolTable *pSym) const
{
	if(runinfo.IsExecDisabled()) return 0;
	return m_pSymbol;
}

Symbol* NodeString::eval(ParseInfo &info, RuntimeInfo& runinfo, SymbolTable *pSym) const
{
	if(runinfo.IsExecDisabled()) return 0;
	return m_pSymbol;
}

Symbol* NodePair::eval(ParseInfo &info, RuntimeInfo& runinfo, SymbolTable *pSym) const
{
	if(runinfo.IsExecDisabled()) return 0;

	std::ostringstream ostrErr;
	ostrErr << linenr(runinfo)
		<< "Pairs should not be evaluated directly." 
		<< std::endl;
	throw tl::Err(ostrErr.str(),0);
	return 0;
}

Symbol* NodeMap::eval(ParseInfo &info, RuntimeInfo& runinfo, SymbolTable *pSym) const
{
	if(runinfo.IsExecDisabled()) return 0;

	NodeBinaryOp *pMap = (NodeBinaryOp*)m_pMap;
	std::vector<Node*> vecNodes;
	if(pMap) vecNodes = pMap->flatten(NODE_ARGS);

	SymbolMap *pSymMap = new SymbolMap;

	for(Node* pNode : vecNodes)
	{
		if(!pNode) continue;
		if(pNode->GetType() != NODE_PAIR)
		{
			std::ostringstream ostrErr;
			ostrErr << linenr(runinfo)
				<< "Maps have to consist of key-value pairs." 
				<< std::endl;
			throw tl::Err(ostrErr.str(),0);
			//continue;
		}

		Symbol* pSymFirst = 0;
		Symbol* pSymSecond = 0;

		if(((NodePair*)pNode)->GetFirst())
			pSymFirst = ((NodePair*)pNode)->GetFirst()->eval(info, runinfo, pSym);
		if(((NodePair*)pNode)->GetSecond())
			pSymSecond = ((NodePair*)pNode)->GetSecond()->eval(info, runinfo, pSym);

		bool bSecondCloned = 0;
		if(pSymFirst && pSymSecond)
		{
			Symbol *pSymClone;
			bSecondCloned = clone_if_needed(pSymSecond, pSymClone);

			pSymMap->GetMap().insert(
				SymbolMap::t_map::value_type(SymbolMapKey(pSymFirst),
				pSymClone));
		}

		safe_delete(pSymFirst, pSym, &info);
		if(bSecondCloned)
			safe_delete(pSymSecond, pSym, &info);
	}

	pSymMap->UpdateIndices();
	return pSymMap;
}

Symbol* NodeArray::eval(ParseInfo &info, RuntimeInfo& runinfo, SymbolTable *pSym) const
{
	if(runinfo.IsExecDisabled()) return 0;

	NodeBinaryOp *pArr = (NodeBinaryOp*)m_pArr;
	std::vector<Node*> vecNodes;
	if(pArr) vecNodes = pArr->flatten(NODE_ARGS);

	SymbolArray *pSymArr = new SymbolArray;
	pSymArr->GetArr().reserve(vecNodes.size());

	for(Node* pNode : vecNodes)
	{
		if(!pNode) continue;

		Symbol *pSymbol = pNode->eval(info, runinfo, pSym);
		Symbol *pClone;
		if(clone_if_needed(pSymbol, pClone))
			safe_delete(pSymbol, pSym, &info);

		pSymArr->GetArr().push_back(pClone);
		pSymArr->UpdateLastNIndices(1);
	}
	return pSymArr;
}

Symbol* NodeRange::eval(ParseInfo &info, RuntimeInfo& runinfo, SymbolTable *pSym) const
{
	if(runinfo.IsExecDisabled()) return 0;

	std::ostringstream ostrErr;
	ostrErr << linenr(runinfo) << "Range should not be evaluated directly." << std::endl;
	throw tl::Err(ostrErr.str(), 0);
}

void NodeRange::GetRangeIndices(ParseInfo &info, RuntimeInfo& runinfo, SymbolTable *pSym,
					t_int iMaxLen, t_int& iBeginIdx, t_int& iEndIdx)
{
	if(m_rangetype == RANGE_FULL)
	{
		iBeginIdx = 0;
		iEndIdx = iMaxLen;
	}
	else if(m_rangetype == RANGE_BEGINEND)
	{
		if(m_pBegin==0 || m_pEnd==0)
		{
			std::ostringstream ostrErr;
			ostrErr << linenr(runinfo) << "Invalid range." << std::endl;
			throw tl::Err(ostrErr.str(),0);
		}

		Symbol *pSymBeg = m_pBegin->eval(info, runinfo, pSym);
		Symbol *pSymEnd = m_pEnd->eval(info, runinfo, pSym);

		iBeginIdx = pSymBeg->GetValInt();
		iEndIdx = pSymEnd->GetValInt();

		// convert negative indices
		if(iBeginIdx < 0) iBeginIdx = iMaxLen + iBeginIdx;
		if(iEndIdx < 0) iEndIdx = iMaxLen + iEndIdx;


		if(iBeginIdx<0 || iBeginIdx>=iMaxLen)
		{
			tl::log_err(linenr(runinfo), "Lower array index out of bounds: ",
				iBeginIdx, ". Adjusting to lower limit.");
			iBeginIdx = 0;
		}
		if(iEndIdx<-1 || iEndIdx>iMaxLen)
		{
			tl::log_err(linenr(runinfo), "Upper array index out of bounds: ",
				iEndIdx, ". Adjusting to upper limit.");
			iEndIdx = iMaxLen;
		}

		safe_delete(pSymBeg, pSym, &info);
		safe_delete(pSymEnd, pSym, &info);
	}
	else
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "Invalid range operation." << std::endl;
		throw tl::Err(ostrErr.str(),0);
	}
}

unsigned int get_max_cols(SymbolArray* pArr, std::vector<unsigned int>* pvecCols=0)
{
	unsigned int iCols = 0;
	for(Symbol* pSym : pArr->GetArr())
	{
		if(!pSym)
		{
			if(pvecCols) pvecCols->push_back(0);
			continue;
		}
		if(pSym->GetType() == SYMBOL_ARRAY)
			iCols = std::max<unsigned int>(iCols, ((SymbolArray*)pSym)->GetArr().size());
		else
			iCols = std::max<unsigned int>(iCols, 1);

		if(pvecCols) pvecCols->push_back(iCols);
	}

	return iCols;
}

Symbol* get_mat_elem(SymbolArray* pArr, unsigned int iLine, unsigned int iCol)
{
	if(iLine >= pArr->GetArr().size())
		return 0;

	Symbol* pLine = pArr->GetArr()[iLine];
	if(pLine->GetType() == SYMBOL_ARRAY)
	{
		SymbolArray* pArrLine = (SymbolArray*)pLine;
		if(iCol >= pArrLine->GetArr().size())
			return 0;

		return pArrLine->GetArr()[iCol]->clone();
	}
	else	// single scalar element
	{
		if(iCol == 0)
			return pLine;
		else
			return 0;
	}
}

SymbolArray* transpose(SymbolArray* pArr, std::vector<unsigned int>* pvecCols=0)
{
	unsigned int iLines = pArr->GetArr().size();
	unsigned int iCols = get_max_cols(pArr, pvecCols);

	SymbolArray *pArrNew = new SymbolArray();
	pArrNew->GetArr().reserve(iCols);

	for(unsigned int iCol=0; iCol<iCols; ++iCol)
	{
		SymbolArray* pArrCol = new SymbolArray();
		pArrCol->GetArr().reserve(iLines);

		for(unsigned int iLine=0; iLine<iLines; ++iLine)
		{
			Symbol *pElem = get_mat_elem(pArr, iLine, iCol);
			if(!pElem) pElem = new SymbolInt(0);
			pArrCol->GetArr().push_back(pElem);
		}

		pArrCol->UpdateIndices();
		pArrNew->GetArr().push_back(pArrCol);
	}

	pArrNew->UpdateIndices();
	return pArrNew;
}

Symbol* NodeArrayAccess::eval(ParseInfo &info, RuntimeInfo& runinfo, SymbolTable *pSym) const
{
	if(runinfo.IsExecDisabled()) return 0;

	if(!m_pIdent /*|| m_pIdent->GetType() != NODE_IDENT*/)
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "Tried to access non-array." << std::endl;
		throw tl::Err(ostrErr.str(),0);
	}

	t_string strIdent;
	Symbol* pSymbol = m_pIdent->eval(info, runinfo, pSym);

	if(!pSymbol)
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) 
			<< "Symbol for array not found." << std::endl;
		throw tl::Err(ostrErr.str(),0);
	}

	if(m_pIdent->GetType() == NODE_IDENT)
		strIdent = pSymbol->GetIdent();
	else
		strIdent = T_STR"<tmp_sym>";


	if(pSymbol->GetType() == SYMBOL_ARRAY)
	{
		bool bCreatedSym = 0;
		bool bEvenIndex = 0;
		bool bAlreadyTransposed = 0;

		for(const Node *pIndices : m_vecIndices)
		{
			SymbolArray *pArr = (SymbolArray*)pSymbol;

			if(bCreatedSym /*&& bEvenIndex*/)
			{
				pArr = transpose(pArr);
				safe_delete(pSymbol, pSym, &info);
				pSymbol = pArr;

				bAlreadyTransposed = 1;
			}

			// TODO: assigment for ranged access
			if(pIndices->GetType() == NODE_RANGE)	// range index
			{
				NodeRange *pRange = (NodeRange*)pIndices;
				t_int iBeginIdx = 0, iEndIdx = 0;
				pRange->GetRangeIndices(info, runinfo, pSym, pArr->GetArr().size(), iBeginIdx, iEndIdx);

				t_int iStep = 1;
				t_int iSize = iEndIdx - iBeginIdx;
				if(iEndIdx < iBeginIdx)
				{
					iSize = -iSize;
					iStep = -1;
				}

				SymbolArray *pSubArr = new SymbolArray();
				pSubArr->GetArr().reserve(iSize);

				for(t_int iIdx=iBeginIdx, iNewIdx=0; iIdx!=iEndIdx && iNewIdx<iSize; iIdx+=iStep, ++iNewIdx)
				{
					Symbol *pElemClone = pArr->GetArr()[iIdx]->clone();
					pSubArr->GetArr().push_back(pElemClone);
					pSubArr->UpdateIndex(iNewIdx);
				}

				if(bAlreadyTransposed)
				{
					SymbolArray *pOrgSubArr = pSubArr;
					pSubArr = transpose(pSubArr);
					delete pOrgSubArr;

					bAlreadyTransposed = 0;
				}

				if(bCreatedSym)
				{
					delete pSymbol;
					pSymbol = 0;
					bCreatedSym = 0;
				}

				safe_delete(pSymbol, pSym, &info);
				pSymbol = pSubArr;
				bCreatedSym = 1;
			}
			else								// integer index
			{
				Symbol *pSymExpr = pIndices->eval(info, runinfo, pSym);
				if(pSymExpr==0 || pSymExpr->GetType()!=SYMBOL_INT)
				{
					std::ostringstream ostrErr;
					ostrErr << linenr(runinfo)
							<< "Array index has to be of integer type."
							<< std::endl;
					throw tl::Err(ostrErr.str(),0);
				}

				t_int iIdx = pSymExpr->GetValInt();
				safe_delete(pSymExpr, pSym, &info);

				// convert negative indices
				if(iIdx < 0)
					iIdx = pArr->GetArr().size()  + iIdx;

				if(iIdx < 0)
				{
					std::ostringstream ostrErr;
					ostrErr << linenr(runinfo)
						<< "Invalid array index."
						<< std::endl;
					throw tl::Err(ostrErr.str(),0);
				}

				// index too high -> fill up with zeroes
				if(iIdx>=int(pArr->GetArr().size()))
				{
					unsigned int iOldSize = pArr->GetArr().size();
					for(unsigned int iRem=0; iRem<iIdx+1-iOldSize; ++iRem)
					{
						SymbolReal *pNewSym = new SymbolReal(0.);
						pNewSym->SetConst(1);
						pArr->GetArr().push_back(pNewSym);
						//G_COUT << "Inserting: " << iRem << std::endl;
					}
				}

				if((void*)pSymbol != (void*)pArr)
					safe_delete(pSymbol, pSym, &info);

				pSymbol = pArr->GetArr()[iIdx];
				pArr->UpdateIndex(iIdx);
			}

			bEvenIndex = !bEvenIndex;
		}
	}
	else if(pSymbol->GetType() == SYMBOL_MAP)
	{
		if(m_vecIndices.size()==0)
		{
			std::ostringstream ostrErr;
			ostrErr << linenr(runinfo) << "No key given for map." << std::endl;
			throw tl::Err(ostrErr.str(),0);
		}
		else if(m_vecIndices.size()>1)
		{
			tl::log_warn(linenr(runinfo), "Multiple keys given for map, using first one.");
		}

		const Node *pNodeKey = m_vecIndices[0];
		Symbol *pSymExpr = pNodeKey->eval(info, runinfo, pSym);
		if(pSymExpr==0)
		{
			std::ostringstream ostrErr;
			ostrErr << linenr(runinfo) << "Map key is invalid."
						<< std::endl;
			throw tl::Err(ostrErr.str(),0);
		}

		SymbolMap *pMap = (SymbolMap*)pSymbol;
		SymbolMapKey key = SymbolMapKey(pSymExpr);
		SymbolMap::t_map::iterator iterMap = pMap->GetMap().find(key);

		std::ostringstream ostrHash;
		ostrHash << std::hex << "0x" << key.key;

		// key not yet in map -> insert it
		if(iterMap == pMap->GetMap().end())
		{
			SymbolString *pNewSym = new SymbolString();
			pNewSym->SetConst(1);
			iterMap = pMap->GetMap().insert(SymbolMap::t_map::value_type(key, pNewSym)).first;
		}

		if((void*)pSymbol != (void*)iterMap->second)
			safe_delete(pSymbol, pSym, &info);

		pSymbol = iterMap->second;
		pMap->UpdateIndex(key);
		safe_delete(pSymExpr, pSym, &info);
	}
	else if(pSymbol->GetType() == SYMBOL_STRING)
	{
		SymbolString *pStr = (SymbolString*)pSymbol;
		t_int iStrLen = pStr->GetVal().length();

		if(m_vecIndices.size()!=1)
		{
			std::ostringstream ostrErr;
			ostrErr << linenr(runinfo)
					<< "Need exactly one string index."
					<< std::endl;
			throw(tl::Err(ostrErr.str()),0);
		}

		const Node *pIndices = m_vecIndices[0];

		if(pIndices->GetType() == NODE_RANGE)	// range index
		{
			NodeRange *pRange = (NodeRange*)pIndices;
			t_int iBeginIdx = 0, iEndIdx = 0;
			pRange->GetRangeIndices(info, runinfo, pSym, iStrLen, iBeginIdx, iEndIdx);

			t_int iStep = 1;
			t_int iSize = iEndIdx - iBeginIdx;
			if(iEndIdx < iBeginIdx)
			{
				iSize = -iSize;
				iStep = -1;
			}

			t_string strNew;
			strNew.resize(iSize);

			for(t_int iIdx=iBeginIdx, iNewIdx=0; iIdx!=iEndIdx && iNewIdx<iSize; iIdx+=iStep, ++iNewIdx)
				strNew[iNewIdx] = pStr->GetVal()[iIdx];

			SymbolString *pSubStr = new SymbolString(strNew);
			safe_delete(pSymbol, pSym, &info);
			pSymbol = pSubStr;
		}
		else								// integer index
		{
			Symbol *pSymExpr = pIndices->eval(info, runinfo, pSym);
			if(pSymExpr==0 || pSymExpr->GetType()!=SYMBOL_INT)
			{
				std::ostringstream ostrErr;
				ostrErr << linenr(runinfo)
						<< "String index has to be of integer type."
						<< std::endl;
				throw tl::Err(ostrErr.str(),0);
			}

			t_int iIdx = pSymExpr->GetValInt();
			safe_delete(pSymExpr, pSym, &info);

			// convert negative indices
			if(iIdx < 0)
				iIdx = iStrLen  + iIdx;

			if(iIdx < 0 || iIdx >= iStrLen)
			{
				std::ostringstream ostrErr;
				ostrErr << linenr(runinfo) << "String index out of bounds."
					<< std::endl;

				throw tl::Err(ostrErr.str(),0);
			}

			t_string strNew;
			strNew += pStr->GetVal()[iIdx];
			safe_delete(pSymbol, pSym, &info);
			pSymbol = new SymbolString(strNew);
		}
	}
	else
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "Symbol \"" << strIdent
				<< "\" is neither an array nor a map." << std::endl;
		throw tl::Err(ostrErr.str(),0);
	}

	return pSymbol;
}

static void uminus_inplace(Symbol* pSym, ParseInfo& info, RuntimeInfo& runinfo)
{
	if(!pSym) return;

	if(pSym->GetType() == SYMBOL_DOUBLE)
		((SymbolReal*)pSym)->SetVal(-((SymbolReal*)pSym)->GetVal());
	else if(pSym->GetType() == SYMBOL_INT)
		((SymbolInt*)pSym)->SetVal(-((SymbolInt*)pSym)->GetVal());
	else if(pSym->GetType() == SYMBOL_COMPLEX)
		((SymbolComplex*)pSym)->SetVal(-((SymbolComplex*)pSym)->GetVal());
	else if(pSym->GetType() == SYMBOL_ARRAY)
	{
		for(Symbol* pElem : ((SymbolArray*)pSym)->GetArr())
			uminus_inplace(pElem, info, runinfo);
	}
	/*else if(pSym->GetType() == SYMBOL_MAP)
	{
		for(SymbolMap::t_map::value_type& pair : ((SymbolMap*)pSym)->GetMap())
			uminus_inplace(pair.second, info, runinfo);
	}*/
	else
	{
		tl::log_err(linenr(runinfo), "Unary minus not defined for ", pSym->GetTypeName(), ".");
	}
}

Symbol* NodeUnaryOp::eval(ParseInfo &info, RuntimeInfo& runinfo, SymbolTable *pSym) const
{
	if(runinfo.IsExecDisabled()) return 0;

	switch(GetType())
	{
		case NODE_UMINUS:
		{
			Symbol *pSymbolEval = m_pChild->eval(info, runinfo, pSym);
			Symbol *pClone;
			if(clone_if_needed(pSymbolEval, pClone))
				safe_delete(pSymbolEval, pSym, &info);

			uminus_inplace(pClone, info, runinfo);
			return pClone;
		}

		case NODE_LOG_NOT:
		{
			Symbol *pSymbolEval = m_pChild->eval(info, runinfo, pSym);
			SymbolInt *pSymbolInt = new SymbolInt();

			if(pSymbolEval->GetType() == SYMBOL_DOUBLE)
				pSymbolInt->SetVal(!((SymbolReal*)pSymbolEval)->GetVal());
			else if(pSymbolEval->GetType() == SYMBOL_INT)
				pSymbolInt->SetVal(!((SymbolInt*)pSymbolEval)->GetVal());

			safe_delete(pSymbolEval, pSym, &info);
			return pSymbolInt;
		}

		case NODE_STMTS:
		{
			if(m_pChild)
			{
				Symbol *pSymbol = m_pChild->eval(info, runinfo, pSym);
				safe_delete(pSymbol, pSym, &info);
			}
			return 0;
		}

		case NODE_UNPACK:
		{
			tl::log_warn(linenr(runinfo), "Unpack operation only allowed in function call. Ignoring.");

			if(m_pChild)
				return m_pChild->eval(info, runinfo, pSym);
		}
		default:
			tl::log_warn(linenr(runinfo), "Unknown node type: ", GetType());
			break;
	}

	return 0;
}

Symbol* NodeBinaryOp::eval_assign(ParseInfo &info, RuntimeInfo& runinfo, SymbolTable *pSym,
 				Node *pLeft, Node *pRight, Symbol *pSymRightAlt,
				const bool *pbGlob) const
{
	if(pLeft==0) pLeft = m_pLeft;
	if(pRight==0) pRight = m_pRight;
	if(pbGlob==0) pbGlob = &m_bGlobal;

	if(pLeft==0 || pRight==0)
	{
		tl::log_err(linenr(runinfo), "NULL assignment.");
		return 0;
	}

	Symbol *pSymbolOrg = 0;
	if(pSymRightAlt)	// use RHS symbol if given instead of RHS node
		pSymbolOrg = pSymRightAlt;
	else
		pSymbolOrg = pRight->eval(info, runinfo, pSym);

	if(!pSymbolOrg)
	{
		tl::log_err(linenr(runinfo), "Invalid rhs expression in assignment.");
		return 0;
	}

	Symbol *pSymbol;
	if(clone_if_needed(pSymbolOrg, pSymbol))
		safe_delete(pSymbolOrg, pSym, &info);
	pSymbol->SetRval(0);

	if(pLeft->GetType() == NODE_IDENT)		// single variable
	{
		const t_string& strIdent = ((NodeIdent*)pLeft)->GetIdent();
		//tl::log_debug("Assigning ", strIdent, " = ", pSymbol);

		Symbol* pSymGlob = 0;
		if(info.pGlobalSyms)
		{
			std::lock_guard<std::mutex> _lck(*info.pmutexGlobalSyms);
			pSymGlob = info.pGlobalSyms->GetSymbol(strIdent);
		}

		Symbol* pSymLoc = 0;
		if(!*pbGlob && pSym)
			pSymLoc = pSym->GetSymbol(strIdent);

		if(pSymLoc && pSymGlob)
		{
			tl::log_warn(linenr(runinfo), "Symbol \"", strIdent,
				"\" exists in local and global scope, using local one.");
		}

		if(pSymGlob && !pSymLoc && !*pbGlob && info.pGlobalSyms)
		{
			tl::log_warn(linenr(runinfo), "Overwriting global symbol \"", strIdent, "\".");

			std::lock_guard<std::mutex> _lck(*info.pmutexGlobalSyms);
			info.pGlobalSyms->InsertSymbol(strIdent, pSymbol);
		}
		else
		{
			if(*pbGlob && info.pGlobalSyms)
			{
				std::lock_guard<std::mutex> _lck(*info.pmutexGlobalSyms);
				info.pGlobalSyms->InsertSymbol(strIdent, pSymbol);
			}
			else if(pSym)
			{
				pSym->InsertSymbol(strIdent, pSymbol);
			}
		}

		return pSymbol;
	}
	else if(pLeft->GetType() == NODE_ARRAY)				// e.g. [a,b] = [1,2];
	{	// TODO: check that LHS does not want eval!
		if(pSymbol->GetType() != SYMBOL_ARRAY)
		{
			tl::log_err(linenr(runinfo), "Assignment needs an array on the right-hand side.");
			return pSymbol;
		}
		SymbolArray *pArrRight = (SymbolArray*)pSymbol;

		std::vector<Node*> vecLeftArgs = ((NodeBinaryOp*)((NodeArray*)pLeft)->GetArr())->flatten(NODE_ARGS);
		if(vecLeftArgs.size() != pArrRight->GetArr().size())
		{
                        tl::log_warn(linenr(runinfo),
				"Size mismatch between assigned and returned array: ",
				vecLeftArgs.size(), " != ", pArrRight->GetArr().size(), ".");
		}

		for(unsigned int iArr=0; iArr<vecLeftArgs.size(); ++iArr)
		{
			Node *pNodeLeft = vecLeftArgs[iArr];
			Symbol *pSymRight = pArrRight->GetArr()[iArr];

			//std::cout << ((NodeIdent*)pNode)->GetIdent() << std::endl;
			eval_assign(info, runinfo, pSym, pNodeLeft, 0, pSymRight, pbGlob);
		}

		return pSymbol;
	}
	else								// array or map
	{
		Symbol *pSymLeft = pLeft->eval(info, runinfo, pSym);
		if(!pSymLeft)
		{
			std::ostringstream ostrErr;
			ostrErr << linenr(runinfo)
				<< "No array element found." << std::endl;
			throw tl::Err(ostrErr.str(),0);
		}

		if(pSymLeft->GetType() == pSymbol->GetType())
		{
			pSymLeft->assign(pSymbol);
		}
		else if(pSymLeft->GetArrPtr())		// array
		{
			t_int iArrIdx = pSymLeft->GetArrIdx();
			SymbolArray* pArr = pSymLeft->GetArrPtr();

			if(int(pArr->GetArr().size()) <= iArrIdx)
			{
				unsigned int iOldSize = pArr->GetArr().size();
				for(unsigned int iRem=0; iRem<iArrIdx+1-iOldSize; ++iRem)
				{
					SymbolReal *pNewSym = new SymbolReal(0.);
					pNewSym->SetConst(1);
					pArr->GetArr().push_back(pNewSym);
					pArr->UpdateLastNIndices(1);
				}
			}


			Symbol* pSymOld = pArr->GetArr()[iArrIdx];
			if((void*)pSymOld != (void*)pSymLeft)
			{
				std::ostringstream ostrErr;
				ostrErr << linenr(runinfo)
					<< "Array member mismatch." << std::endl;
				throw tl::Err(ostrErr.str(),0);
			}


			pArr->GetArr()[iArrIdx] = pSymbol;
			pSymbol->SetArrPtr(pArr);
			pSymbol->SetArrIdx(iArrIdx);

			pSymOld->SetArrPtr(0);
			safe_delete(pSymOld, pSym, &info);
		}
		else if(pSymLeft->GetMapPtr())		// map
		{
			const SymbolMapKey& MapKey = pSymLeft->GetMapKey();
			SymbolMap* pMap = pSymLeft->GetMapPtr();

			SymbolMap::t_map::iterator iterMap = pMap->GetMap().find(MapKey);

			// symbol not in map -> insert a zero
			if(iterMap == pMap->GetMap().end())
			{
				SymbolReal *pNewSym = new SymbolReal(0.);
				pNewSym->SetConst(1);
				iterMap = pMap->GetMap().insert(SymbolMap::t_map::value_type(MapKey, pNewSym)).first;
			}

			Symbol* pSymOld = iterMap->second;
			if((void*)pSymOld != (void*)pSymLeft)
			{
				std::ostringstream ostrErr;
				ostrErr << linenr(runinfo) << "Map member mismatch." << std::endl;
				throw tl::Err(ostrErr.str(),0);
			}

			pSymbol->SetMapPtr(pMap);
			pSymbol->SetMapKey(MapKey);
			iterMap->second = pSymbol;

			pSymOld->SetMapPtr(0);
			safe_delete(pSymOld, pSym, &info);
		}
		else
		{
			std::ostringstream ostrErr;
			ostrErr << linenr(runinfo)
				<< "Trying to access array/map member with no associated array/map."
				<< std::endl;
			throw tl::Err(ostrErr.str(),0);
		}

		//safe_delete(pSymLeft, pSym, &info);
	}
	return pSymbol;
}

Symbol* NodeBinaryOp::eval_funcinit(ParseInfo &info, RuntimeInfo& runinfo, SymbolTable *pSym) const
{
	std::vector<NodeFunction*>& vecFuncs = info.vecFuncs;
	NodeFunction *pToRun = 0;

	//tl::log_debug("# funcs: ", m_vecNodesFlat.size());
	for(Node *_pNodeFunc : m_vecNodesFlat)
	{
		if(!_pNodeFunc || _pNodeFunc->GetType()!=NODE_FUNC)
			continue;

		NodeFunction *pNodeFunc = (NodeFunction*)_pNodeFunc;
		pNodeFunc->SetScrFile(runinfo.strInitScrFile);

		const t_string& strFktName = pNodeFunc->GetName();

		if(strFktName == T_STR"__init__" || strFktName == T_STR"module_init")
		{
			pToRun = pNodeFunc;
		}
		else if(info.GetFunction(strFktName))
		{
			if(strFktName != "main")
				tl::log_warn(linenr(runinfo),
					"Function \"", strFktName,
					"\" redefined in \"", runinfo.strInitScrFile, "\".",
					" Ignoring.");
		}
		else
		{
			if(has_ext_call(strFktName))
			{
				tl::log_warn(linenr(runinfo),
					"Function \"", strFktName,
					"\" in \"", runinfo.strInitScrFile,
					"\" overwrites a system function.");
			}
			vecFuncs.push_back(pNodeFunc);

			if(strFktName != T_STR"__cmd__")
				tl::log_debug("Imported function \"", strFktName, 
					"\" from module \"", runinfo.strInitScrFile, "\".");
		}
	}

	// execute general entry point function
	if(pToRun)
	{
		t_string strExecBck = runinfo.strExecFkt;

		Symbol *pSymInitRet = pToRun->eval(info, runinfo, /*pSym*/0);
		safe_delete(pSymInitRet, pSym, &info);

		runinfo.strExecFkt = strExecBck;
	}

	// execute named entry point function
	if(runinfo.strExecFkt != T_STR"")
	{
		for(NodeFunction* pFkt : vecFuncs)
		{
			if(pFkt->GetName() == runinfo.strExecFkt)
			{
				if(pSym && pFkt->GetArgVec().size() == 0)
					pSym->RemoveSymbolNoDelete(T_STR"<args>");

				//tl::log_info("Executing ", runinfo.strExecFkt);
				Symbol *pSymRet = pFkt->eval(info, runinfo, pSym);
				if(pSym)
					pSym->RemoveSymbolNoDelete(T_STR"<args>");

				return pSymRet;
			}
		}

		tl::log_err(linenr(runinfo), "Function \"", runinfo.strExecFkt, "\" not defined.");
	}
	return 0;
}

Symbol* NodeBinaryOp::eval_recursive(ParseInfo &info, RuntimeInfo& runinfo, SymbolTable *pSym) const
{
	// prefer eval_sequential
	tl::log_warn(linenr(runinfo), "Non-optimal statement list.");

	if(m_pLeft)
	{
		if(runinfo.IsExecDisabled()) return 0;

		Symbol *pSymbol = m_pLeft->eval(info, runinfo, pSym);
		safe_delete(pSymbol, pSym, &info);
	}
	if(m_pRight)
	{
		if(runinfo.IsExecDisabled()) return 0;

		Symbol *pSymbol = m_pRight->eval(info, runinfo, pSym);
		safe_delete(pSymbol, pSym, &info);
	}
	return 0;
}

Symbol* NodeBinaryOp::eval_sequential(ParseInfo &info, RuntimeInfo& runinfo, SymbolTable *pSym) const
{
	Symbol *pSymRet = 0;

	for(std::size_t iNode=0; iNode<m_vecNodesFlat.size(); ++iNode)
	{
		if(runinfo.IsExecDisabled()) break;

		Node *pNode = m_vecNodesFlat[iNode];
		Symbol *pSymbol = pNode->eval(info, runinfo, pSym);
		if(runinfo.bImplicitRet && iNode==m_vecNodesFlat.size()-1)
		{
			pSymRet = pSymbol /*? pSymbol->clone() : 0*/;
			break;
		}

		safe_delete(pSymbol, pSym, &info);
	}

	return pSymRet;
}

Symbol* NodeBinaryOp::eval(ParseInfo &info, RuntimeInfo& runinfo, SymbolTable *pSym) const
{
	if(runinfo.IsExecDisabled()) return 0;
	const bool bIsFunctionRoot = (GetType()==NODE_FUNC || GetType()==NODE_FUNCS);

	if(!bIsFunctionRoot && m_vecNodesFlat.size())
		return eval_sequential(info, runinfo, pSym);

	switch(GetType())
	{
		case NODE_STMTS:
		//case NODE_ARGS:
			return eval_recursive(info, runinfo, pSym);

		case NODE_ASSIGN: return eval_assign(info, runinfo, pSym);

		// should only be called once per module
		case NODE_FUNCS: return eval_funcinit(info, runinfo, pSym);

		default: break;
	};

	Symbol *pSymbolLeft = m_pLeft->eval(info, runinfo, pSym);
	Symbol *pSymbolRight = 0;
	Symbol *pSymbol = 0;

	// optimisation: 0 && x == 0
	if(GetType() == NODE_LOG_AND && pSymbolLeft->GetValInt()==0)
		pSymbol = new SymbolInt(0);
	// optimisation: 1 || x == 1
	else if(GetType() == NODE_LOG_OR && pSymbolLeft->GetValInt()==1)
		pSymbol = new SymbolInt(1);
	else
	{
		pSymbolRight = m_pRight->eval(info, runinfo, pSym);
		pSymbol = Op(pSymbolLeft, pSymbolRight, GetType());
	}

	if(pSymbol!=pSymbolLeft) safe_delete(pSymbolLeft, pSym, &info);
	if(pSymbol!=pSymbolRight) safe_delete(pSymbolRight, pSym, &info);

	return pSymbol;
}


Symbol* NodeFunction::eval(ParseInfo &info, RuntimeInfo& runinfo, SymbolTable* pTableSup) const
{
	if(runinfo.IsExecDisabled()) return 0;
	runinfo.pCurFunction = this;

	const t_string& strName = GetName();

	bool bOverrideSymTab = 0;
	if(runinfo.pLocalSymsOverride && strName == runinfo.strExecFkt)
		bOverrideSymTab = 1;

	std::unique_ptr<SymbolTable> ptrLocalSym(bOverrideSymTab ? 0 : new SymbolTable);
	SymbolTable *pLocalSym = ptrLocalSym.get();

	// externally supplied symbol table
	if(bOverrideSymTab)
		pLocalSym = runinfo.pLocalSymsOverride;

	SymbolArray* pArgs = 0;
	if(pTableSup)
		pArgs = (SymbolArray*)pTableSup->GetSymbol(T_STR"<args>");

	unsigned int iArgSize = 0;
	if(pArgs)
	{
		const std::vector<Symbol*> *pVecArgSyms = &pArgs->GetArr();

		if(m_vecArgs.size() != pVecArgSyms->size())
		{
			tl::log_warn(linenr(runinfo), "Function \"",
					strName, "\"", " takes ",
					m_vecArgs.size(), " arguments, but ",
					pVecArgSyms->size(), " given.");
		}

		iArgSize = /*std::min(*/m_vecArgs.size()/*, pVecArgSyms->size())*/;
		for(unsigned int iArg=0; iArg<iArgSize; ++iArg)
		{
			const NodeIdent* pIdent = (NodeIdent*)m_vecArgs[iArg];
			const Node* pDefArg = pIdent->GetDefArg();

			Symbol *pSymbol = 0;

			if(iArg < pVecArgSyms->size())		// argument given by caller
				pSymbol = (*pVecArgSyms)[iArg];
			if(pSymbol==0 && pDefArg)		// default argument
			{
				pSymbol = pDefArg->eval(info, runinfo, pTableSup);

				if(iArg < pVecArgSyms->size() && pSymbol)
					tl::log_warn(linenr(runinfo),
						"Given argument \"", pIdent->GetIdent(), 
						"\" for function \"", 
						strName, "\" not valid. ", 
						"Using default argument.");
			}
			if(pSymbol==0)
			{
				tl::log_err(linenr(runinfo), "Argument \"",
					pIdent->GetIdent(), "\" for function \"",
					strName, "\" not given. Ignoring.");
				continue;
			}

			//G_COUT << "arg: " << pIdent->GetIdent() << std::endl;
			Symbol* pSymToInsert;
			if(clone_if_needed(pSymbol, pSymToInsert))
				safe_delete(pSymbol, pLocalSym, &info);
			pSymToInsert->SetRval(0);
			pLocalSym->InsertSymbol(pIdent->GetIdent(), pSymToInsert);
		}
	}


	if(info.bEnableDebug)
	{
		std::string strTrace = "call: " + GetName() + ", " 
					+ std::to_string(iArgSize) + " args";
		info.PushTraceback(std::move(strTrace));
	}

	Symbol *pRet = 0;
	if(m_pStmts)
	{
		pRet = m_pStmts->eval(info, runinfo, pLocalSym);
		if(!pRet)
		{
			pRet = pLocalSym->GetSymbol(T_STR"<ret>");
			pLocalSym->RemoveSymbolNoDelete(T_STR"<ret>");
		}
	}

	if(info.bEnableDebug)
	{
		info.PopTraceback();
	}

	//tl::log_debug("Local symbols for \"", strName, "\":");
	//pLocalSym->print();

	return pRet;
}


Symbol* NodeIf::eval(ParseInfo &info, RuntimeInfo& runinfo, SymbolTable *pSym) const
{
	if(runinfo.IsExecDisabled()) return 0;

	Symbol *pSymExpr = 0;
	Symbol *pSymRet = 0;
	if(m_pExpr)
		pSymExpr = m_pExpr->eval(info, runinfo, pSym);

	if(pSymExpr && pSymExpr->IsNotZero())
		pSymRet = (m_pIf ? m_pIf->eval(info, runinfo, pSym) : 0);
	else
		pSymRet = (m_pElse ? m_pElse->eval(info, runinfo, pSym) : 0);

	safe_delete(pSymExpr, pSym, &info);
	safe_delete(pSymRet, pSym, &info);

	return 0;
}


Symbol* NodeWhile::eval(ParseInfo &info, RuntimeInfo& runinfo, SymbolTable *pSym) const
{
	if(runinfo.IsExecDisabled()) return 0;

	if(!m_pExpr) return 0;
	if(!m_pStmt) return 0;

	const Node* pLastLoop = runinfo.pCurLoop;
	runinfo.pCurLoop = this;
	while(1)
	{
		Symbol *pSymRet = 0;
		Symbol *pSymExpr = m_pExpr->eval(info, runinfo, pSym);

		if(pSymExpr && pSymExpr->IsNotZero())
			pSymRet = m_pStmt->eval(info, runinfo, pSym);
		else
			break;

		safe_delete(pSymRet, pSym, &info);
		safe_delete(pSymExpr, pSym, &info);

		if(runinfo.bWantBreak)
		{
			runinfo.bWantBreak = 0;
			break;
		}
		if(runinfo.bWantContinue)
		{
			runinfo.bWantContinue = 0;
			continue;
		}
	}
	runinfo.pCurLoop = pLastLoop;

	return 0;
}

Symbol* NodeFor::eval(ParseInfo &info, RuntimeInfo& runinfo, SymbolTable *pSym) const
{
	if(runinfo.IsExecDisabled()) return 0;

	if(!m_pExprCond) return 0;
	if(!m_pStmt) return 0;

	const Node* pLastLoop = runinfo.pCurLoop;
	runinfo.pCurLoop = this;

	if(m_pExprInit)
	{
		Symbol *pSymInit = m_pExprInit->eval(info, runinfo, pSym);
		safe_delete(pSymInit, pSym, &info);
	}

	while(1)
	{
		Symbol *pSymRet = 0;
		Symbol *pSymExpr = m_pExprCond->eval(info, runinfo, pSym);

		if(pSymExpr && pSymExpr->IsNotZero())
			pSymRet = m_pStmt->eval(info, runinfo, pSym);
		else
			break;

		safe_delete(pSymRet, pSym, &info);
		safe_delete(pSymExpr, pSym, &info);

		if(runinfo.bWantBreak)
		{
			runinfo.bWantBreak = 0;
			break;
		}

		if(runinfo.bWantContinue)
		{
			runinfo.bWantContinue = 0;
			//continue;
		}

		if(m_pExprEnd)
		{
			Symbol *pSymEnd = m_pExprEnd->eval(info, runinfo, pSym);
			safe_delete(pSymEnd, pSym, &info);
		}
	}
	runinfo.pCurLoop = pLastLoop;

	return 0;
}

Symbol* NodeRangedFor::eval(ParseInfo &info, RuntimeInfo& runinfo, SymbolTable *pSym) const
{
	if(runinfo.IsExecDisabled()) return 0;
	if(!m_pIdent || !m_pExpr || !m_pStmt) return 0;

	if(m_pIdent->GetType() != NODE_IDENT)
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) 
			<< "Range-based for loop needs identifier." << std::endl;
		throw tl::Err(ostrErr.str(),0);
	}

	Symbol *_pArr = m_pExpr->eval(info, runinfo, pSym);
	if(!_pArr)
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo)
			<< "Invalid array for loop." << std::endl;
		return 0;
	}

	if(_pArr->GetType() != SYMBOL_ARRAY)
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo)
				<< "Range-based for loop needs array." << std::endl;
		safe_delete(_pArr, pSym, &info);

		throw tl::Err(ostrErr.str(),0);
	}

	SymbolArray *pArr = (SymbolArray*)_pArr;


	const t_string& strIdent = ((NodeIdent*)m_pIdent)->GetIdent();

	SymbolInt *pSymIter = new SymbolInt(0);
	t_string strIter = T_STR"<cur_iter_" + strIdent + T_STR">";
	pSym->InsertSymbol(strIter, pSymIter);

	const Node* pLastLoop = runinfo.pCurLoop;
	runinfo.pCurLoop = this;
	for(unsigned int iArr=0; iArr<pArr->GetArr().size(); ++iArr)
	{
		Symbol *pSymInArr = pArr->GetArr()[iArr];
		pSym->InsertSymbol(strIdent, pSymInArr);

		Symbol *pBodyRet = m_pStmt->eval(info, runinfo, pSym);
		safe_delete(pBodyRet, pSym, &info);


		// write back symbol in case an assignment has taken place
		Symbol *pNewSym = pSym->GetSymbol(strIdent);
		if(pSymInArr != pNewSym)
		{
			pArr->GetArr()[iArr] = pNewSym;
			//delete pSymInArr;
			pSymInArr = pNewSym;
			pArr->UpdateIndices();
		}
		pSym->RemoveSymbolNoDelete(strIdent);


		++pSymIter->GetVal();

		if(runinfo.bWantBreak)
		{
			runinfo.bWantBreak = 0;
			break;
		}
		if(runinfo.bWantContinue)
		{
			runinfo.bWantContinue = 0;
			continue;
		}
	}
	runinfo.pCurLoop = pLastLoop;

	pSym->RemoveSymbol(strIter);

	safe_delete(_pArr, pSym, &info);
	return 0;
}
