/*
 * external file functions
 * @author tweber
 * @date dec-2013
 * @license GPLv2 or GPLv3
 */

#include "lang/types.h"
#include "string/string.h"
#include "file/file.h"
#include "log/log.h"
#include "calls_file.h"
#include "lang/calls.h"
#include "file/loaddat.h"
#include "file/loadinstr.h"
#include <sstream>
#include <fstream>
#include <iomanip>
#include <algorithm>

// --------------------------------------------------------------------------------
// file operations
static Symbol* fkt_file_exists(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_STRING}, {0}, "file_exists"))
		return 0;

	const t_string& strFile = ((SymbolString*)vecSyms[0])->GetVal();
	t_ifstream ifstr(strFile);

	bool bFileExists = ifstr.is_open();
	return new SymbolInt(bFileExists);
}

static Symbol* fkt_read_file(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_STRING}, {0}, "read_file"))
		return 0;

	const t_string& strFile = ((SymbolString*)vecSyms[0])->GetVal();

	t_ifstream ifstr(strFile);
	if(!ifstr.is_open())
		return 0;

	t_ostringstream ostr;

	std::copy(std::istreambuf_iterator<t_char>(ifstr),
		std::istreambuf_iterator<t_char>(),
		std::ostreambuf_iterator<t_char>(ostr));

	return new SymbolString(ostr.str());
}

static Symbol* fkt_write_file(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_STRING, SYMBOL_ANY}, {0,0}, "write_file"))
		return 0;

	const t_string& strFile = ((SymbolString*)vecSyms[0])->GetVal();
	t_ofstream ofstr(strFile);
	if(!ofstr.is_open())
		return new SymbolInt(0);


	t_string* pStr = 0;
	bool bAllocatedStr = 0;
	if(vecSyms[1]->GetType() == SYMBOL_STRING)
		pStr = &((SymbolString*)vecSyms[1])->GetVal();
	else
	{
		pStr = new std::string;
		bAllocatedStr = 1;
		*pStr = vecSyms[1]->print();
	}

	std::copy(pStr->begin(), pStr->end(),
		std::ostreambuf_iterator<t_char>(ofstr));

	if(bAllocatedStr)
		delete pStr;
	return new SymbolInt(1);
}
// --------------------------------------------------------------------------------



// --------------------------------------------------------------------------------
// loading and saving of .dat files

static Symbol* fkt_loadtxt(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_STRING}, {0}, "loadtxt"))
		return 0;

	const t_string& strFile = ((SymbolString*)vecSyms[0])->GetVal();
	SymbolArray *pArr = new SymbolArray();

	tl::DatFile<t_real, typename t_string::value_type> dat;
	bool bLoaded = dat.Load(strFile);
	if(!bLoaded)
	{
		tl::log_err(linenr(runinfo), "loadtxt could not open \"", strFile, "\".");
		return pArr;
	}

	pArr->GetArr().reserve(dat.GetColumnCount());
	for(std::size_t iCol=0; iCol<dat.GetColumnCount(); ++iCol)
	{
		const std::vector<t_real>& vecCol = dat.GetColumn(iCol);
		const std::size_t iColLen = vecCol.size();

		SymbolArray *pArrCol = new SymbolArray();
		pArrCol->GetArr().reserve(iColLen);

		for(std::size_t iRow=0; iRow<iColLen; ++iRow)
		{
			SymbolReal* pSymD = new SymbolReal();
			pSymD->SetVal(vecCol[iRow]);

			pArrCol->GetArr().push_back(pSymD);
		}

		pArrCol->UpdateIndices();
		pArr->GetArr().push_back(pArrCol);
	}

	// load the parameter map
	SymbolMap *pSymMap = new SymbolMap();
	for(const typename decltype(dat)::t_map::value_type &val : dat.GetHeader())
	{
		t_string strKey = STR_TO_WSTR(val.first);

		SymbolString *pSymStrVal = new SymbolString();
		pSymStrVal->SetVal(STR_TO_WSTR(val.second));

		pSymMap->GetMap().insert(SymbolMap::t_map::value_type(strKey, pSymStrVal));
	}
	pArr->GetArr().push_back(pSymMap);

	pArr->UpdateIndices();
	return pArr;
}

static Symbol* fkt_loadinstr(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	//if(!check_args(runinfo, vecSyms, {SYMBOL_STRING}, {0}, "loadinstr"))
	//	return 0;

	if(vecSyms.size() != 1)
	{
		tl::log_err(linenr(runinfo), "loadinstr needs one argument of string or array type.");
		return 0;
	}

	std::vector<tl::FileInstrBase<t_real>*> vecToMerge;

	t_string strMainFile;
	if(vecSyms[0]->GetType() == SYMBOL_STRING)
		strMainFile = ((SymbolString*)vecSyms[0])->GetVal();
	else if(vecSyms[0]->GetType() == SYMBOL_ARRAY)
	{
		const std::vector<Symbol*>& vecArr = ((SymbolArray*)vecSyms[0])->GetArr();
		if(vecArr.size() == 0 || vecArr[0]->GetType() != SYMBOL_STRING)
		{
			tl::log_err(linenr(runinfo), "loadinstr got no valid file name.");
			return 0;
		}

		strMainFile = ((SymbolString*)vecArr[0])->GetVal();

		for(std::size_t iFile=1; iFile<vecArr.size(); ++iFile)
		{
			Symbol *pSymAux = vecArr[iFile];
			if(!pSymAux || pSymAux->GetType() != SYMBOL_STRING)
			{
				tl::log_warn(linenr(runinfo), "loadinstr cannot load file ", iFile, " in array.");
				continue;
			}

			const t_string& strFileAux = ((SymbolString*)pSymAux)->GetVal();

			tl::FileInstrBase<t_real> *pToMerge = tl::FileInstrBase<t_real>::LoadInstr(strFileAux.c_str());
			if(!pToMerge)
			{
				tl::log_warn(linenr(runinfo), "loadinstr could not open \"", strFileAux, "\" for merging.");
				continue;
			}
			else
				vecToMerge.push_back(pToMerge);
		}
	}


	SymbolMap* pmapRet = new SymbolMap();
	SymbolInt* pOk = new SymbolInt(0);
	pmapRet->GetMap().insert(SymbolMap::t_map::value_type(t_string("ok"), pOk));

	const t_string& strFile = strMainFile;
	tl::FileInstrBase<t_real>* pInstr = tl::FileInstrBase<t_real>::LoadInstr(strFile.c_str());
	if(!pInstr)
	{
		tl::log_err(linenr(runinfo), "loadinstr could not open \"", strFile, "\".");
		return pmapRet;
	}

	for(tl::FileInstrBase<t_real>* pToMerge : vecToMerge)
	{
		pInstr->MergeWith(pToMerge);
		delete pToMerge;
	}
	vecToMerge.clear();


	const tl::FileInstrBase<t_real>::t_vecDat& vecDat = pInstr->GetData();
	const tl::FileInstrBase<t_real>::t_vecColNames& vecColNames = pInstr->GetColNames();
	const tl::FileInstrBase<t_real>::t_mapParams& mapParams = pInstr->GetAllParams();
	std::vector<std::string> vecScanVars = pInstr->GetScannedVars();


	// data
	SymbolArray* pArrDat = new SymbolArray();
	pmapRet->GetMap().insert(SymbolMap::t_map::value_type(t_string("data"), pArrDat));

	pArrDat->GetArr().reserve(vecDat.size());
	for(const tl::FileInstrBase<t_real>::t_vecDat::value_type& vecRow : vecDat)
	{
		SymbolArray *pArrRow = new SymbolArray();
		pArrDat->GetArr().push_back(pArrRow);
		pArrRow->GetArr().reserve(vecRow.size());

		for(const tl::FileInstrBase<t_real>::t_vecDat::value_type::value_type& val : vecRow)
			pArrRow->GetArr().push_back(new SymbolReal(val));
	}


	// column labels
	SymbolArray* pArrLab = new SymbolArray();
	pmapRet->GetMap().insert(SymbolMap::t_map::value_type(t_string("labels"), pArrLab));

	pArrLab->GetArr().reserve(vecColNames.size());
	for(const tl::FileInstrBase<t_real>::t_vecColNames::value_type& strLab : vecColNames)
		pArrLab->GetArr().push_back(new SymbolString(strLab));


	// param map
	SymbolMap *pmapParams = new SymbolMap();
	pmapRet->GetMap().insert(SymbolMap::t_map::value_type(t_string("params"), pmapParams));
	for(const tl::FileInstrBase<t_real>::t_mapParams::value_type& pair : mapParams)
	{
		const t_string& strKey = pair.first;
		SymbolString *pVal = new SymbolString(pair.second);

		//std::cout << "Inserting \"" << strKey << "\" = \"" << pair.second << "\"" << std::endl;
		pmapParams->GetMap().insert(SymbolMap::t_map::value_type(strKey, pVal));
	}


	// scan vars
	SymbolArray* pArrSc = new SymbolArray();
	pmapRet->GetMap().insert(SymbolMap::t_map::value_type(t_string("scanvars"), pArrSc));

	pArrSc->GetArr().reserve(vecScanVars.size());
	for(const std::string& strVar : vecScanVars)
		pArrSc->GetArr().push_back(new SymbolString(strVar));



	// misc
	pmapRet->GetMap().insert(SymbolMap::t_map::value_type(t_string("cki"), new SymbolInt(pInstr->IsKiFixed())));
	pmapRet->GetMap().insert(SymbolMap::t_map::value_type(t_string("kfix"), new SymbolReal(pInstr->GetKFix())));
	pmapRet->GetMap().insert(SymbolMap::t_map::value_type(t_string("countvar"), new SymbolString(pInstr->GetCountVar())));
	pmapRet->GetMap().insert(SymbolMap::t_map::value_type(t_string("monvar"), new SymbolString(pInstr->GetMonVar())));


	delete pInstr;
	pOk->SetVal(1);
	return pmapRet;
}

static void get_2darr_size(const SymbolArray* pArr,
	unsigned int& iColLen, unsigned int& iRowLen)
{
	iColLen = pArr->GetArr().size();
	iRowLen = 0;

	if(iColLen)
	{
		// look for first real array (not the parameter map)
		for(unsigned int iCol=0; iCol<iColLen; ++iCol)
		{
			Symbol* pSym = pArr->GetArr()[iCol];
			if(pSym->GetType() == SYMBOL_ARRAY)
			{
				iRowLen = ((SymbolArray*)pSym)->GetArr().size();
				break;
			}
		}
	}

	unsigned int iNonArray = 0;
	// don't count the parameter map
	for(unsigned int iCol=0; iCol<iColLen; ++iCol)
	{
		Symbol* pSym = pArr->GetArr()[iCol];
		if(pSym->GetType() != SYMBOL_ARRAY)
			++iNonArray;
	}

	iColLen -= iNonArray;
}

static std::string get_2darr_strval(const SymbolArray* pArr,
	unsigned int iCol, unsigned int iRow)
{
	unsigned int iColLen = pArr->GetArr().size();
	if(iCol >= iColLen)
		return "0";

	bool bFoundCol = 0;
	unsigned int iColRealArray = 0;
	for(unsigned int iCurCol=0; iCurCol<iColLen; ++iCurCol)
	{
		Symbol *pSym = pArr->GetArr()[iCurCol];
		if(pSym->GetType() == SYMBOL_ARRAY)
		{
			if(iColRealArray == iCol)
			{
				bFoundCol = 1;
				iCol = iCurCol;
				break;
			}

			++iColRealArray;
		}
	}

	if(!bFoundCol)
	{
		tl::log_err("Invalid column index: ", iCol, ".");
		return "0";
	}


	Symbol *pSym = pArr->GetArr()[iCol];

	const std::vector<Symbol*>& veccol = ((SymbolArray*)pSym)->GetArr();
	if(iRow >= veccol.size())
		return "0";

	return veccol[iRow]->print();
}

static Symbol* fkt_savetxt(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_STRING, SYMBOL_ARRAY}, {0,0}, "savetxt"))
		return 0;

	const t_string& strFile = ((SymbolString*)vecSyms[0])->GetVal();
	SymbolArray* pArr = (SymbolArray*)vecSyms[1];

	t_ofstream ofstr(WSTR_TO_STR(strFile).c_str());
	if(!ofstr.is_open())
	{
		tl::log_err(linenr(runinfo), "Cannot open \"", strFile, "\".");
		return 0;
	}

	// save parameter map
	for(unsigned int iCol=0; iCol<pArr->GetArr().size(); ++iCol)
	{
		Symbol *_pSymMap = pArr->GetArr()[iCol];
		if(_pSymMap && _pSymMap->GetType() == SYMBOL_MAP)
		{
			SymbolMap *pSymMap = (SymbolMap*)_pSymMap;

			for(const SymbolMap::t_map::value_type& val : pSymMap->GetMap())
			{
				ofstr << "# " << val.first.strKey << " : "
					<< (val.second?val.second->print():T_STR"") << "\n";
			}
		}
	}

	unsigned int iColLen=0, iRowLen=0;
	get_2darr_size(pArr, iColLen, iRowLen);
	//G_COUT << "col len: " << iColLen << ", row len: " << iRowLen << std::endl;

	for(unsigned int iRow=0; iRow<iRowLen; ++iRow)
	{
		for(unsigned int iCol=0; iCol<iColLen; ++iCol)
			ofstr << std::setw(20) << std::left << get_2darr_strval(pArr, iCol, iRow) << " ";
		ofstr << "\n";
	}

	ofstr.flush();
	ofstr.close();
	return 0;
}

// --------------------------------------------------------------------------------



extern void init_ext_file_calls()
{
	t_mapFkts mapFkts =
	{
		// general files
		t_mapFkts::value_type(T_STR"read_file", fkt_read_file),
		t_mapFkts::value_type(T_STR"write_file", fkt_write_file),

		t_mapFkts::value_type(T_STR"file_exists", fkt_file_exists),

		// dat files
		t_mapFkts::value_type(T_STR"loadtxt", fkt_loadtxt),
		t_mapFkts::value_type(T_STR"savetxt", fkt_savetxt),

		// special intrument dat files
		t_mapFkts::value_type(T_STR"loadinstr", fkt_loadinstr),
	};

	add_ext_calls(mapFkts);
}
