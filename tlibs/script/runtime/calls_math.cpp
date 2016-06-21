/*
 * external math functions
 * @author tweber
 * @date 2013
 * @license GPLv2 or GPLv3
 */

#include "lang/types.h"
#include "calls_math.h"
#include "lang/calls.h"
#include "math/math.h"
#include "math/fourier.h"
#include "math/linalg.h"
#include "math/linalg2.h"
#include "math/rand.h"
#include "log/log.h"
#include <boost/math/special_functions/erf.hpp>

namespace ublas = boost::numeric::ublas;


static inline Symbol* _fkt_linlogspace(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab, bool bLog)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_SCALAR, SYMBOL_SCALAR, SYMBOL_INT}, {0,0,0}, "linspace"))
		return 0;

	t_int iNum = vecSyms[2]->GetValInt();
	t_real dmin = vecSyms[0]->GetValDouble();
	t_real dmax = vecSyms[1]->GetValDouble();

	SymbolArray* pSymRet = new SymbolArray;
	pSymRet->GetArr().reserve(iNum);
	for(t_int i=0; i<iNum; ++i)
	{
		SymbolReal *pSymD = new SymbolReal();
		t_real dDiv = (iNum!=1 ? t_real(iNum-1) : 1);
		pSymD->SetVal(t_real(i)*(dmax-dmin)/dDiv + dmin);
		if(bLog)
		{
			const t_real dBase = 10.;
			pSymD->SetVal(std::pow(dBase, pSymD->GetVal()));
		}
		pSymRet->GetArr().push_back(pSymD);
	}

	pSymRet->UpdateIndices();
	return pSymRet;
}

static Symbol* fkt_linspace(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	return _fkt_linlogspace(vecSyms, info, runinfo, pSymTab, false);
}

static Symbol* fkt_logspace(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	return _fkt_linlogspace(vecSyms, info, runinfo, pSymTab, true);
}



// --------------------------------------------------------------------------------
// math
enum MathFkts
{
	MATH_MAX,
	MATH_MIN
};

template<MathFkts fkt, typename T>
const T& math_fkt(const T& t1, const T& t2)
{
	if(fkt == MATH_MAX)
		return std::max<T>(t1, t2);
	else if(fkt == MATH_MIN)
		return std::min<T>(t1, t2);

	std::ostringstream ostrErr;
	ostrErr << "Error: Invalid function selected in math_fkt." << std::endl;
	throw tl::Err(ostrErr.str(),0);
}

template<MathFkts fkt>
static Symbol* fkt_math_for_every(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(vecSyms.size() < 1)
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "fkt_math_for_every needs at least one argument" << std::endl;
		throw tl::Err(ostrErr.str(),0);
	}

	t_real dRes = 0.;
	t_int iRes = 0;

	bool bHadInt = 0,
		bHadDouble = 0;

	for(Symbol* pSym : vecSyms)
	{
		Symbol *pThisSym = pSym;
		bool bCleanSym = 0;
		if(pSym->GetType() == SYMBOL_ARRAY)
		{
			pThisSym = fkt_math_for_every<fkt>(
					((SymbolArray*)pSym)->GetArr(),
					info, runinfo, pSymTab);

			bCleanSym = 1;
		}

		if(pThisSym->GetType() == SYMBOL_INT)
		{
			if(!bHadInt)
				iRes = ((SymbolInt*)pThisSym)->GetVal();
			else
				iRes = math_fkt<fkt, t_int>(iRes, ((SymbolInt*)pThisSym)->GetVal());

			bHadInt = 1;
		}
		else if(pThisSym->GetType() == SYMBOL_DOUBLE)
		{
			if(!bHadDouble)
				dRes = ((SymbolReal*)pThisSym)->GetVal();
			else
				dRes = math_fkt<fkt, t_real>(dRes, ((SymbolReal*)pThisSym)->GetVal());

			bHadDouble = 1;
		}

		if(bCleanSym)
			delete pThisSym;
	}

	if(bHadInt && !bHadDouble)
		return new SymbolInt(iRes);
	else if(bHadInt && bHadDouble)
	{
		dRes = math_fkt<fkt, t_real>(dRes, t_real(iRes));
		return new SymbolReal(dRes);
	}
	else if(!bHadInt && bHadDouble)
		return new SymbolReal(dRes);


	std::ostringstream ostrErr;
	ostrErr << linenr(runinfo) << "No valid arguments given for fkt_math_for_every." << std::endl;
	throw tl::Err(ostrErr.str(), 0);
	//return 0;
}


template<t_real (*FKT)(t_real), t_complex (*FKT_C)(const t_complex&)=nullptr>
static Symbol* fkt_math_1arg(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_ANY}, {0}, "fkt_math_1arg"))
		return 0;

	if(vecSyms[0]->GetType() == SYMBOL_ARRAY)
	{
		SymbolArray* pArrRet = new SymbolArray();

		SymbolArray* pSymArr = (SymbolArray*)vecSyms[0];
		for(Symbol* pArrElem : pSymArr->GetArr())
		{
			std::vector<Symbol*> vecDummy;
			vecDummy.push_back(pArrElem);

			pArrRet->GetArr().push_back(fkt_math_1arg<FKT, FKT_C>(vecDummy, info, runinfo, pSymTab));
		}

		pArrRet->UpdateIndices();
		return pArrRet;
	}
	else
	{
		if(vecSyms[0]->GetType() == SYMBOL_COMPLEX || (FKT_C && !FKT))
		{
			if(FKT_C == nullptr)
			{
				tl::log_err(linenr(runinfo), "Undefined complex function.");
				return 0;
			}

			t_complex cResult = FKT_C(vecSyms[0]->GetValComplex());
			return new SymbolComplex(std::move(cResult));
		}
		else
		{
			if(FKT == nullptr)
			{
				tl::log_err(linenr(runinfo), "Undefined real function.");
				return 0;
			}


			t_real tVal = vecSyms[0]->GetValDouble();

			if(((FKT==static_cast<t_real(*)(t_real)>(std::sqrt) && tVal<0.) ||
			   (FKT==static_cast<t_real(*)(t_real)>(std::log) && tVal<0.) ||
			   (FKT==static_cast<t_real(*)(t_real)>(std::log10) && tVal<0.))
				&& FKT_C!=nullptr)
			{
				t_complex tCompVal = t_complex(tVal, 0.);
				return new SymbolComplex(FKT_C(tCompVal));
			}

			t_real dResult = FKT(tVal);
			return new SymbolReal(dResult);
		}
	}

	return 0;
}

// TODO: for arrays
template<t_real (*FKT)(t_real, t_real), t_complex (*FKT_C)(const t_complex&, const t_complex&)=nullptr>
static Symbol* fkt_math_2args(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_ANY, SYMBOL_ANY}, {0,0}, "fkt_math_2args"))
		return 0;

	Symbol* pFirst = vecSyms[0];
	Symbol* pSecond = vecSyms[1];
	Symbol* pToDel = 0;

	// complex
	bool bFirstComplex = (pFirst->GetType()==SYMBOL_COMPLEX);
	bool bSecondComplex = (pSecond->GetType()==SYMBOL_COMPLEX);

	if(bFirstComplex && !bSecondComplex)
		pSecond = pToDel = pSecond->ToType(SYMBOL_COMPLEX);
	else if(!bFirstComplex && bSecondComplex)
		pFirst = pToDel = pFirst->ToType(SYMBOL_COMPLEX);

	if(bFirstComplex || bSecondComplex)
	{
		if(FKT_C == nullptr)
		{
			tl::log_err(linenr(runinfo), "Undefined complex function.");
			return 0;
		}

		t_complex cResult = FKT_C(((SymbolComplex*)pFirst)->GetVal(), 
					((SymbolComplex*)pSecond)->GetVal());

		if(pToDel) delete pToDel;
		return new SymbolComplex(std::move(cResult));
	}


	// real
	if(FKT == nullptr)
	{
		tl::log_err(linenr(runinfo), "Undefined real function.");
		return 0;
	}

	t_real dResult = FKT(pFirst->GetValDouble(), pSecond->GetValDouble());
	return new SymbolReal(dResult);
}


template<typename T> T myabs(T t)
{
	if(t < T(0))
		return -t;
	return t;
}

template<> t_real myabs(t_real t) { return ::fabs(t); }


// TODO: integrate with fkt_math_1arg
static Symbol* fkt_math_abs(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_ANY}, {0}, "abs"))
		return 0;

	if(vecSyms[0]->GetType() == SYMBOL_ARRAY)
	{
		SymbolArray* pArrRet = new SymbolArray();

		SymbolArray* pSymArr = (SymbolArray*)vecSyms[0];
		for(Symbol* pArrElem : pSymArr->GetArr())
		{
			std::vector<Symbol*> vecDummy;
			vecDummy.push_back(pArrElem);

			pArrRet->GetArr().push_back(fkt_math_abs(vecDummy, info, runinfo, pSymTab));
		}

		pArrRet->UpdateIndices();
		return pArrRet;
	}
	else if(vecSyms[0]->GetType() == SYMBOL_INT)
	{
		SymbolInt* pSymInt = (SymbolInt*)vecSyms[0];
		return new SymbolInt(myabs(pSymInt->GetVal()));
	}
	else if(vecSyms[0]->GetType() == SYMBOL_DOUBLE)
	{
		SymbolReal* pSymD = (SymbolReal*)vecSyms[0];
		return new SymbolReal(myabs(pSymD->GetVal()));
	}
	else if(vecSyms[0]->GetType() == SYMBOL_COMPLEX)
	{
		SymbolComplex* pSymC = (SymbolComplex*)vecSyms[0];
		return new SymbolReal(std::abs<t_real>(pSymC->GetVal()));
	}


	std::ostringstream ostrErr;
	ostrErr << linenr(runinfo) << "abs received unsupported symbol type." << std::endl;
	throw tl::Err(ostrErr.str(),0);
	//return 0;
}


// complex norm = abs^2
static Symbol* fkt_math_cnorm(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_ANY}, {0}, "cnorm"))
		return 0;

	t_complex cVal(0., 0.);

	if(vecSyms[0]->GetType() == SYMBOL_COMPLEX)
		cVal = ((SymbolComplex*)vecSyms[0])->GetVal();
	else
		cVal.real(vecSyms[0]->GetValDouble());

	t_real dVal = std::norm<t_real>(cVal);
	return new SymbolReal(dVal);
}

// functions with boolean return value
template<bool (*FKT)(t_real)>
static Symbol* fkt_math_1arg_bret(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_ANY}, {0}, "fkt_math_1arg_bret"))
		return 0;

	if(vecSyms[0]->GetType() == SYMBOL_COMPLEX)
	{
		tl::log_err(linenr(runinfo), "Function not implemented for complex values.");
		return 0;
	}
	else if(vecSyms[0]->GetType() == SYMBOL_ARRAY)
	{
		SymbolArray* pArrRet = new SymbolArray();

		SymbolArray* pSymArr = (SymbolArray*)vecSyms[0];
		for(Symbol* pArrElem : pSymArr->GetArr())
		{
			std::vector<Symbol*> vecDummy;
			vecDummy.push_back(pArrElem);

			pArrRet->GetArr().push_back(fkt_math_1arg_bret<FKT>(vecDummy, info, runinfo, pSymTab));
		}

		pArrRet->UpdateIndices();
		return pArrRet;
	}
	else
	{
		t_int iResult = FKT(vecSyms[0]->GetValDouble());
		return new SymbolInt(iResult);
	}

	return 0;
}


// --------------------------------------------------------------------------------
// FFT

static Symbol* _fkt_fft(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab,
	bool bInv)
{
	void (tl::Fourier_gen<t_real>::*pFkt)(const t_real*, const t_real*, t_real*, t_real*)
		= (bInv ? &tl::Fourier_gen<t_real>::ifft : &tl::Fourier_gen<t_real>::fft);

	bool bArgsOk=1;
	std::vector<t_real> vecRealIn, vecImagIn;

	// real and imag part as separate arguments
	if(vecSyms.size()==2 && vecSyms[0]->GetType()==SYMBOL_ARRAY
			&& vecSyms[1]->GetType()==SYMBOL_ARRAY)
	{
		vecRealIn = ((SymbolArray*)vecSyms[0])->ToDoubleArray();
		vecImagIn = ((SymbolArray*)vecSyms[1])->ToDoubleArray();
	}
	// arrays in one array
	else if(vecSyms.size()==1 && vecSyms[0]->GetType()==SYMBOL_ARRAY)
	{
		SymbolArray* pSymArr = (SymbolArray*)vecSyms[0];
		unsigned int iSymArrSize = pSymArr->GetArr().size();

		if(iSymArrSize==0)
			bArgsOk = 0;
		if(iSymArrSize>=1 && pSymArr->GetArr()[0]->GetType()==SYMBOL_ARRAY)
			vecRealIn = ((SymbolArray*)pSymArr->GetArr()[0])->ToDoubleArray();
		if(iSymArrSize>=2 && pSymArr->GetArr()[1]->GetType()==SYMBOL_ARRAY)
			vecImagIn = ((SymbolArray*)pSymArr->GetArr()[1])->ToDoubleArray();

		// simple array containing real data
		if(pSymArr->GetArr()[0]->GetType()!=SYMBOL_ARRAY)
			vecRealIn = pSymArr->ToDoubleArray();
	}

	if(!bArgsOk)
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "fft received invalid arguments." << std::endl;
		throw tl::Err(ostrErr.str(),0);
	}

	if(vecRealIn.size() != vecImagIn.size())
	{
		unsigned int iSize = std::max(vecRealIn.size(), vecImagIn.size());
		vecRealIn.resize(iSize);
		vecImagIn.resize(iSize);
	}

	std::vector<t_real> vecRealOut, vecImagOut;
	vecRealOut.resize(vecRealIn.size());
	vecImagOut.resize(vecImagIn.size());

	tl::Fourier_gen<t_real> fourier(vecRealIn.size());
	(fourier.*pFkt)(vecRealIn.data(), vecImagIn.data(),
		vecRealOut.data(), vecImagOut.data());

	SymbolArray* pArrReal = new SymbolArray();
	SymbolArray* pArrImag = new SymbolArray();

	pArrReal->FromDoubleArray(vecRealOut);
	pArrImag->FromDoubleArray(vecImagOut);

	SymbolArray* pRet = new SymbolArray();
	pRet->GetArr().push_back(pArrReal);
	pRet->GetArr().push_back(pArrImag);

	return pRet;
}

static Symbol* fkt_fft(const std::vector<Symbol*>& vecSyms, ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{ return _fkt_fft(vecSyms, info, runinfo, pSymTab, false); }
static Symbol* fkt_ifft(const std::vector<Symbol*>& vecSyms, ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{ return _fkt_fft(vecSyms, info, runinfo, pSymTab, true); }

// --------------------------------------------------------------------------------


// --------------------------------------------------------------------------------
// linalg stuff

template<typename T=t_real> using t_vec = ublas::vector<T>;
template<typename T=t_real> using t_mat = ublas::matrix<T>;

template<typename T=t_real> using t_stdvec = std::vector<T>;


static Symbol* fkt_length(const std::vector<Symbol*>& vecSyms, ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_ARRAY}, {0}, "len"))
		return 0;

	if(!is_vec(vecSyms[0]))
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "len needs a vector argument." << std::endl;
		throw tl::Err(ostrErr.str(),0);
	}

	t_vec<t_real> vec = sym_to_vec<t_vec>(vecSyms[0]);
	t_real dLen = std::sqrt(ublas::inner_prod(vec, vec));
	return new SymbolReal(dLen);
}

static Symbol* fkt_modf(const std::vector<Symbol*>& vecSyms, ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_DOUBLE}, {0}, "modf"))
		return 0;

	t_real d = vecSyms[0]->GetValDouble();
	t_real dIntPart = 0.;
	t_real dFracPart = std::modf(d, &dIntPart);

	return new SymbolArray({new SymbolReal(dIntPart), new SymbolReal(dFracPart)});
}

static Symbol* fkt_mean(const std::vector<Symbol*>& vecSyms, ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_ARRAY}, {0}, "mean"))
		return 0;

	if(!is_vec(vecSyms[0]))
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "Mean value needs vector arguments." << std::endl;
		throw tl::Err(ostrErr.str(),0);
	}

	t_vec<t_real> vecLeft = sym_to_vec<t_vec>(vecSyms[0]);
	t_real dMean = tl::mean_value(vecLeft);

	return new SymbolReal(dMean);
}

static Symbol* fkt_stddev(const std::vector<Symbol*>& vecSyms, ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_ARRAY}, {0}, "stddev"))
		return 0;

	if(!is_vec(vecSyms[0]))
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "Standard deviation needs vector arguments." << std::endl;
		throw tl::Err(ostrErr.str(),0);
	}

	t_vec<t_real> vecLeft = sym_to_vec<t_vec>(vecSyms[0]);
	t_real dStd = tl::std_dev(vecLeft);

	return new SymbolReal(dStd);
}

static Symbol* fkt_cross(const std::vector<Symbol*>& vecSyms, ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_ARRAY, SYMBOL_ARRAY}, {0,0}, "cross"))
		return 0;

	if(!is_vec(vecSyms[0]) || !is_vec(vecSyms[1]))
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "Cross product needs vector arguments." << std::endl;
		throw tl::Err(ostrErr.str(),0);
	}

	t_vec<t_real> vecLeft = sym_to_vec<t_vec>(vecSyms[0]);
	t_vec<t_real> vecRight = sym_to_vec<t_vec>(vecSyms[1]);

	if(vecLeft.size()!=3 || vecRight.size()!=3)
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "Cross product needs 3-vectors." << std::endl;
		throw tl::Err(ostrErr.str(), 0);
	}

	t_vec<t_real> vecCross = tl::cross_3(vecLeft, vecRight);
	return vec_to_sym<t_vec>(vecCross);
}

// matrix(rows, cols)
// matrix(dim)
static Symbol* fkt_matrix(const std::vector<Symbol*>& vecSyms, ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(vecSyms.size()<1)
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "Need size of matrix." << std::endl;
		throw tl::Err(ostrErr.str(), 0);
	}

	t_int iRows = vecSyms[0]->GetValInt();
	t_int iCols = iRows;

	// cols also given
	if(vecSyms.size() >= 2)
		iCols = vecSyms[1]->GetValInt();

	t_mat<t_real> mat = ublas::zero_matrix<t_real>(iRows, iCols);
	return mat_to_sym<t_mat>(mat);
}

static Symbol* fkt_transpose(const std::vector<Symbol*>& vecSyms, 
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_ARRAY}, {0}, "trans"))
		return 0;

	bool bIsMat = 0;
	t_mat<t_real> mat = sym_to_mat<t_mat, t_vec>(vecSyms[0], &bIsMat);
	if(!bIsMat)
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "Transpose needs a matrix." << std::endl;
		throw tl::Err(ostrErr.str(),0);
	}

	t_mat<t_real> mat_trans = ublas::trans(mat);
	return mat_to_sym<t_mat>(mat_trans);
}

static Symbol* fkt_inverse(const std::vector<Symbol*>& vecSyms, 
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_ARRAY}, {0}, "inv"))
		return 0;

	bool bIsMat = 0;
	t_mat<t_real> mat = sym_to_mat<t_mat, t_vec>(vecSyms[0], &bIsMat);
	if(!bIsMat)
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "Inverse needs a matrix." << std::endl;
		throw tl::Err(ostrErr.str(),0);
	}

	t_mat<t_real> mat_inv;
	if(!tl::inverse(mat, mat_inv))
	{
		tl::log_warn(linenr(runinfo), "Matrix inversion failed.");
		return 0;
	}

	return mat_to_sym<t_mat>(mat_inv);
}

static Symbol* fkt_determinant(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_ARRAY}, {0}, "det"))
		return 0;

	bool bIsMat = 0;
	t_mat<t_real> mat = sym_to_mat<t_mat, t_vec>(vecSyms[0], &bIsMat);
	if(!bIsMat || mat.size1()!=mat.size2())
	{
		tl::log_err(linenr(runinfo), "Determinant needs a square matrix.");
		return 0;
	}

	t_real dDet = tl::determinant<t_mat<t_real>>(mat);
	return new SymbolReal(dDet);
}

static Symbol* fkt_unitmatrix(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_INT}, {0}, "unity"))
		return 0;

	t_int iSize = vecSyms[0]->GetValInt();
	t_mat<t_real> mat = tl::unit_matrix<t_mat<t_real>>(iSize);

	return mat_to_sym<t_mat>(mat);
}

static Symbol* fkt_outerproduct(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_ARRAY, SYMBOL_ARRAY}, {0,0}, "outer_prod"))
		return 0;

	if(!is_vec(vecSyms[0]) || !is_vec(vecSyms[1]))
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "Outer product needs two vector arguments." << std::endl;
		throw tl::Err(ostrErr.str(),0);
	}

	t_vec<t_real> vec1 = sym_to_vec<t_vec>(vecSyms[0]);
	t_vec<t_real> vec2 = sym_to_vec<t_vec>(vecSyms[1]);

	t_mat<t_real> mat = ublas::outer_prod(vec1, vec2);
	return mat_to_sym<t_mat>(mat);
}

static Symbol* fkt_dot(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_ARRAY, SYMBOL_ARRAY}, {0,0}, "dot"))
		return 0;

	if(!is_vec(vecSyms[0]) || !is_vec(vecSyms[1]))
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "Inner product needs two vector arguments." << std::endl;
		throw tl::Err(ostrErr.str(),0);
	}

	t_vec<t_real> vec1 = sym_to_vec<t_vec>(vecSyms[0]);
	t_vec<t_real> vec2 = sym_to_vec<t_vec>(vecSyms[1]);

	return new SymbolReal(ublas::inner_prod(vec1, vec2));
}

static Symbol* fkt_product(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_ANY, SYMBOL_ANY}, {0,0}, "prod"))
		return 0;

	Symbol* pRet = 0;

	bool bFirstIsVec = is_vec(vecSyms[0]);
	bool bSecondIsVec = is_vec(vecSyms[1]);

	// dot product
	if(bFirstIsVec && bSecondIsVec)
	{
		pRet = fkt_dot(vecSyms, info, runinfo, pSymTab);
	}
	else
	{
		unsigned int iCols1, iRows1, iCols2, iRows2;
		bool bFirstIsMat = is_mat(vecSyms[0], &iCols1, &iRows1);
		bool bSecondIsMat = is_mat(vecSyms[1], &iCols2, &iRows2);

		// matrix product
		if(bFirstIsMat && bSecondIsMat)
		{
			if(iCols1!=iRows2 /*|| iCols1!=iRows2*/)
			{
				tl::log_err(linenr(runinfo), "Row and column counts of matrices do not match: ",
							"Rows: ", iRows1, ", ", iRows2, 
							", columns: ", iCols1, ", ", iCols2, ".");
				return 0;
			}

			t_mat<t_real> mat1 = sym_to_mat<t_mat, t_vec>(vecSyms[0]);
			t_mat<t_real> mat2 = sym_to_mat<t_mat, t_vec>(vecSyms[1]);

			t_mat<t_real> matProd = ublas::prod(mat1, mat2);
			pRet = mat_to_sym<t_mat>(matProd);
		}
		else if(bFirstIsMat && bSecondIsVec)
		{
			t_mat<t_real> mat = sym_to_mat<t_mat, t_vec>(vecSyms[0]);
			t_vec<t_real> vec = sym_to_vec<t_vec>(vecSyms[1]);

			t_vec<t_real> vecProd = ublas::prod(mat, vec);
			pRet = vec_to_sym<t_vec>(vecProd);
		}
		else if(bFirstIsVec && bSecondIsMat)
		{
			t_vec<t_real> vec = sym_to_vec<t_vec>(vecSyms[0]);
			t_mat<t_real> mat = sym_to_mat<t_mat, t_vec>(vecSyms[1]);

			t_vec<t_real> vecProd = ublas::prod(vec, mat);
			pRet = vec_to_sym<t_vec>(vecProd);
		}
	}

	if(!pRet)
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "Invalid call to prod." << std::endl;
		throw tl::Err(ostrErr.str(),0);
	}
	return pRet;
}

// --------------------------------------------------------------------------------



// --------------------------------------------------------------------------------
// advanced linalg stuff

#ifdef USE_LAPACK
static Symbol* fkt_eigenvecs(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_ARRAY}, {0}, "eigenvecs"))
		return 0;

	if(!is_mat(vecSyms[0]))
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "Invalid call to eigenvecs." << std::endl;
		throw tl::Err(ostrErr.str(),0);
	}

	t_mat<t_real> mat = sym_to_mat<t_mat, t_vec>(vecSyms[0]);

	std::vector<t_vec<t_real>> evecs_real, evecs_imag;
	std::vector<t_real> evals_real, evals_imag;

	bool bOk = tl::eigenvec<t_real>(mat, evecs_real, evecs_imag, evals_real, evals_imag);
	if(!bOk) return 0;

	SymbolArray* pSymRet = new SymbolArray();
	SymbolArray* pSymEvecs_real = new SymbolArray();
	SymbolArray* pSymEvecs_imag = new SymbolArray();
	SymbolArray* pSymEvals_real = new SymbolArray();
	SymbolArray* pSymEvals_imag = new SymbolArray();

	pSymRet->GetArr().push_back(pSymEvecs_real);
	pSymRet->GetArr().push_back(pSymEvecs_imag);
	pSymRet->GetArr().push_back(pSymEvals_real);
	pSymRet->GetArr().push_back(pSymEvals_imag);

	for(const ublas::vector<t_real>& evec : evecs_real)
	{
		Symbol* pEvec = vec_to_sym<t_vec>(evec);
		pSymEvecs_real->GetArr().push_back(pEvec);
	}
	for(const ublas::vector<t_real>& evec : evecs_imag)
	{
		Symbol* pEvec = vec_to_sym<t_vec>(evec);
		pSymEvecs_imag->GetArr().push_back(pEvec);
	}

	for(const t_real dEval : evals_real)
	{
		Symbol *pSymReal = new SymbolReal(dEval);
		pSymEvals_real->GetArr().push_back(pSymReal);
	}
	for(const t_real dEval : evals_imag)
	{
		Symbol *pSymReal = new SymbolReal(dEval);
		pSymEvals_imag->GetArr().push_back(pSymReal);
	}

	pSymRet->UpdateIndices();
	return pSymRet;
}
#else
static Symbol* fkt_eigenvecs(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	std::ostringstream ostrErr;
	ostrErr << linenr(runinfo) << "Function eigenvecs not linked." << std::endl;
	throw tl::Err(ostrErr.str(),0);
}
#endif

static Symbol* fkt_eigenvecs_sym(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_ARRAY}, {0}, "eigenvecs_sym"))
		return 0;

	if(!is_mat(vecSyms[0]))
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "Invalid call to eigenvecs_sym." << std::endl;
		throw tl::Err(ostrErr.str(),0);
	}

	t_mat<t_real> mat = sym_to_mat<t_mat, t_vec>(vecSyms[0]);

	std::vector<t_vec<t_real>> evecs_real;
	std::vector<t_real> evals_real;

	bool bOk = tl::eigenvec_sym<t_real>(mat, evecs_real, evals_real);
	if(!bOk) return 0;

	SymbolArray* pSymRet = new SymbolArray();
	SymbolArray* pSymEvecs_real = new SymbolArray();
	SymbolArray* pSymEvals_real = new SymbolArray();

	pSymRet->GetArr().push_back(pSymEvecs_real);
	pSymRet->GetArr().push_back(pSymEvals_real);

	for(const ublas::vector<t_real>& evec : evecs_real)
	{
		Symbol* pEvec = vec_to_sym<t_vec>(evec);
		pSymEvecs_real->GetArr().push_back(pEvec);
	}

	for(const t_real dEval : evals_real)
	{
		Symbol *pSymReal = new SymbolReal(dEval);
		pSymEvals_real->GetArr().push_back(pSymReal);
	}

	pSymRet->UpdateIndices();
	return pSymRet;
}

static Symbol* fkt_qr(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_ARRAY}, {0}, "qr"))
		return 0;

	if(!is_mat(vecSyms[0]))
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "Invalid call to qr." << std::endl;
		throw tl::Err(ostrErr.str(),0);
	}

	t_mat<t_real> mat = sym_to_mat<t_mat, t_vec>(vecSyms[0]);

	ublas::matrix<t_real> Q, R;
	bool bOk = tl::qr(mat, Q, R);
	if(!bOk) return 0;

	Symbol *pQ = mat_to_sym<t_mat>(Q);
	Symbol *pR = mat_to_sym<t_mat>(R);

	SymbolArray* pSymRet = new SymbolArray();
	pSymRet->GetArr().push_back(pQ);
	pSymRet->GetArr().push_back(pR);

	pSymRet->UpdateIndices();
	return pSymRet;
}

// --------------------------------------------------------------------------------



// --------------------------------------------------------------------------------
// rand stuff


static Symbol* fkt_rand01(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	t_real dRand = tl::rand01<t_real>();
	return new SymbolReal(dRand);
}

static Symbol* fkt_rand_real(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_SCALAR, SYMBOL_SCALAR}, {0,0}, "rand_real"))
		return 0;

	t_real dMin = vecSyms[0]->GetValDouble();
	t_real dMax = vecSyms[1]->GetValDouble();

	t_real dRand = tl::rand_real<t_real>(dMin, dMax);
	return new SymbolReal(dRand);
}

static Symbol* fkt_rand_int(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_SCALAR, SYMBOL_SCALAR}, {0,0}, "rand_int"))
		return 0;

	t_int iMin = vecSyms[0]->GetValInt();
	t_int iMax = vecSyms[1]->GetValInt();

	t_int iRand = tl::rand_int<t_int>(iMin, iMax);
	return new SymbolInt(iRand);
}

static Symbol* fkt_rand_norm(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_SCALAR, SYMBOL_SCALAR}, {0,0}, "rand_norm"))
		return 0;

	t_real dMu = 0.;
	t_real dSigma = 1.;

	if(vecSyms.size() >= 1)
		dMu = vecSyms[0]->GetValDouble();
	if(vecSyms.size() >= 2)
		dSigma = vecSyms[1]->GetValDouble();

	t_real dRand = tl::rand_norm<t_real>(dMu, dSigma);
	return new SymbolReal(dRand);
}

static Symbol* fkt_rand_norm_nd(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_ARRAY, SYMBOL_ARRAY}, {0,0}, "rand_norm_nd"))
		return 0;

	t_stdvec<t_real> vecMu = sym_to_vec<t_stdvec>(vecSyms[0]);
	t_stdvec<t_real> vecSigma = sym_to_vec<t_stdvec>(vecSyms[1]);

	if(vecMu.size() != vecSigma.size())
	{
		tl::log_err(linenr(runinfo), "Mu and sigma arrays have different sizes.");
		return 0;
	}

	t_stdvec<t_real> vecResult = tl::rand_norm_nd<t_stdvec<t_real>, t_real, t_stdvec<t_real>>(vecMu, vecSigma);
	return vec_to_sym<t_stdvec>(vecResult);
}
// --------------------------------------------------------------------------------




// --------------------------------------------------------------------------------
// minmax

static Symbol* fkt_minmax_elem(const std::vector<Symbol*>& vecSyms, 
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_ARRAY}, {0}, "minmax_elem"))
		return 0;

	const Symbol *pSymMin=0, *pSymMax=0;
	unsigned int iIdxMin=0, iIdxMax=0;
	unsigned int iCurIter=0;
	for(const Symbol* pSym : ((SymbolArray*)vecSyms[0])->GetArr())
	{
		if(pSymMin==0) { pSymMin = pSym; iIdxMin=iCurIter; }
		if(pSymMax==0) { pSymMax = pSym; iIdxMax=iCurIter; }

		if(pSym->IsLessThan(*pSymMin))
		{
			iIdxMin = iCurIter;
			pSymMin = pSym;
		}
		if(pSym->IsGreaterThan(*pSymMax))
		{
			iIdxMax = iCurIter;
			pSymMax = pSym;
		}

		++iCurIter;
	}

	return new SymbolArray({new SymbolInt(iIdxMin), new SymbolInt(iIdxMax)});
}

static Symbol* fkt_minmax(const std::vector<Symbol*>& vecSyms,
	ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	SymbolArray* pElem = (SymbolArray*)fkt_minmax_elem(vecSyms, info, runinfo, pSymTab);
	if(!pElem || pElem->GetArr().size()!=2)
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(runinfo) << "Invalid input for minmax." << std::endl;
		throw tl::Err(ostrErr.str(),0);
	}

	int iIdxMin = pElem->GetArr()[0]->GetValInt();
	int iIdxMax = pElem->GetArr()[1]->GetValInt();
	delete pElem;

	//std::cout << iIdxMin << ", " << iIdxMax << std::endl;

	SymbolArray* pArr = (SymbolArray*)vecSyms[0];
	return new SymbolArray({pArr->GetArr()[iIdxMin]->clone(), pArr->GetArr()[iIdxMax]->clone()});
}


// --------------------------------------------------------------------------------

template<class t_ret, class t_arg=t_ret> t_ret cot(t_arg t);
template<> t_real cot(t_real t) { return std::cos(t)/std::sin(t); }
template<> t_complex cot(const t_complex& t) { return std::cos<t_real>(t)/std::sin<t_real>(t); }

// --------------------------------------------------------------------------------



extern void init_ext_math_calls()
{
	tl::init_rand();

	t_mapFkts mapFkts =
	{
		// math stuff
		t_mapFkts::value_type(T_STR"sign", fkt_math_1arg< tl::sign<t_real> >),
		t_mapFkts::value_type(T_STR"sqrt", fkt_math_1arg< std::sqrt, std::sqrt<t_real> >),
		t_mapFkts::value_type(T_STR"cbrt", fkt_math_1arg< std::cbrt >),
		t_mapFkts::value_type(T_STR"exp", fkt_math_1arg< std::exp, std::exp<t_real> >),
		t_mapFkts::value_type(T_STR"exp2", fkt_math_1arg< std::exp2 >),
		t_mapFkts::value_type(T_STR"expm1", fkt_math_1arg< std::expm1 >),
		t_mapFkts::value_type(T_STR"log", fkt_math_1arg< std::log, std::log<t_real> >),
		t_mapFkts::value_type(T_STR"log1p", fkt_math_1arg< std::log1p >),
		t_mapFkts::value_type(T_STR"log10", fkt_math_1arg< std::log10, std::log10<t_real> >),
		t_mapFkts::value_type(T_STR"log2", fkt_math_1arg< std::log2 >),
		t_mapFkts::value_type(T_STR"logb", fkt_math_1arg< std::logb >),
		t_mapFkts::value_type(T_STR"pow", fkt_math_2args< std::pow, std::pow<t_real> >),

		t_mapFkts::value_type(T_STR"sin", fkt_math_1arg< std::sin, std::sin<t_real> >),
		t_mapFkts::value_type(T_STR"cos", fkt_math_1arg< std::cos, std::cos<t_real> >),
		t_mapFkts::value_type(T_STR"tan", fkt_math_1arg< std::tan, std::tan<t_real> >),
		t_mapFkts::value_type(T_STR"cot", fkt_math_1arg< cot<t_real>, cot<t_complex> >),
		t_mapFkts::value_type(T_STR"asin", fkt_math_1arg< std::asin, std::asin<t_real> >),
		t_mapFkts::value_type(T_STR"acos", fkt_math_1arg< std::acos, std::acos<t_real> >),
		t_mapFkts::value_type(T_STR"atan", fkt_math_1arg< std::atan, std::atan<t_real> >),
		t_mapFkts::value_type(T_STR"atan2", fkt_math_2args< std::atan2 >),
		t_mapFkts::value_type(T_STR"hypot", fkt_math_2args< std::hypot >),

		t_mapFkts::value_type(T_STR"sinh", fkt_math_1arg< std::sinh, std::sinh<t_real> >),
		t_mapFkts::value_type(T_STR"cosh", fkt_math_1arg< std::cosh, std::cosh<t_real> >),
		t_mapFkts::value_type(T_STR"tanh", fkt_math_1arg< std::tanh, std::tanh<t_real> >),
		t_mapFkts::value_type(T_STR"asinh", fkt_math_1arg< std::asinh, std::asinh<t_real> >),
		t_mapFkts::value_type(T_STR"acosh", fkt_math_1arg< std::acosh, std::acosh<t_real> >),
		t_mapFkts::value_type(T_STR"atanh", fkt_math_1arg< std::atanh, std::atanh<t_real> >),

		t_mapFkts::value_type(T_STR"erf", fkt_math_1arg< std::erf >),
		t_mapFkts::value_type(T_STR"erfc", fkt_math_1arg< std::erfc >),
		t_mapFkts::value_type(T_STR"erf_inv", fkt_math_1arg< boost::math::erf_inv >),
		t_mapFkts::value_type(T_STR"tgamma", fkt_math_1arg< std::tgamma >),
		t_mapFkts::value_type(T_STR"lgamma", fkt_math_1arg< std::lgamma >),

#ifdef HAS_COMPLEX_ERF
		t_mapFkts::value_type(T_STR"erf_cplx", fkt_math_1arg< nullptr, tl::erf<t_real> >),
		t_mapFkts::value_type(T_STR"erfc_cplx", fkt_math_1arg< nullptr, tl::erfc<t_real> >),
		t_mapFkts::value_type(T_STR"faddeeva", fkt_math_1arg< nullptr, tl::faddeeva<t_real> >),
#endif

		t_mapFkts::value_type(T_STR"round", fkt_math_1arg< std::round >),
		t_mapFkts::value_type(T_STR"trunc", fkt_math_1arg< std::trunc >),
		t_mapFkts::value_type(T_STR"rint", fkt_math_1arg< std::rint >),
		t_mapFkts::value_type(T_STR"nearbyint", fkt_math_1arg< std::nearbyint >),
		t_mapFkts::value_type(T_STR"fmod", fkt_math_2args< std::fmod >),
		t_mapFkts::value_type(T_STR"nextafter", fkt_math_2args< std::nextafter >),
		//t_mapFkts::value_type(T_STR"nexttoward", fkt_math_2args< std::nexttoward >),
		t_mapFkts::value_type(T_STR"ceil", fkt_math_1arg< std::ceil >),
		t_mapFkts::value_type(T_STR"floor", fkt_math_1arg< std::floor >),
		t_mapFkts::value_type(T_STR"abs", fkt_math_abs),
		t_mapFkts::value_type(T_STR"max", fkt_math_for_every<MATH_MAX>),
		t_mapFkts::value_type(T_STR"min", fkt_math_for_every<MATH_MIN>),
		t_mapFkts::value_type(T_STR"fdim", fkt_math_2args< std::fdim >),
		t_mapFkts::value_type(T_STR"remainder", fkt_math_2args< std::remainder >),

		t_mapFkts::value_type(T_STR"conj", fkt_math_1arg<nullptr, std::conj<t_real> >),
		t_mapFkts::value_type(T_STR"cnorm", fkt_math_cnorm),

		t_mapFkts::value_type(T_STR"modf", fkt_modf),

		// statistical stuff
		t_mapFkts::value_type(T_STR"mean", fkt_mean),
		t_mapFkts::value_type(T_STR"stddev", fkt_stddev),

		// minmax
		t_mapFkts::value_type(T_STR"minmax", fkt_minmax),
		t_mapFkts::value_type(T_STR"minmax_elem", fkt_minmax_elem),


		// fft
		t_mapFkts::value_type(T_STR"fft", fkt_fft),
		t_mapFkts::value_type(T_STR"ifft", fkt_ifft),

		// arrays
		t_mapFkts::value_type(T_STR"linspace", fkt_linspace),
		t_mapFkts::value_type(T_STR"logspace", fkt_logspace),

		// vector operations
		t_mapFkts::value_type(T_STR"dot", fkt_dot),
		t_mapFkts::value_type(T_STR"cross", fkt_cross),
		t_mapFkts::value_type(T_STR"outer_prod", fkt_outerproduct),

		t_mapFkts::value_type(T_STR"len", fkt_length),

		// matrix operations
		t_mapFkts::value_type(T_STR"mat", fkt_matrix),
		t_mapFkts::value_type(T_STR"unity", fkt_unitmatrix),
		t_mapFkts::value_type(T_STR"trans", fkt_transpose),
		t_mapFkts::value_type(T_STR"inv", fkt_inverse),
		t_mapFkts::value_type(T_STR"det", fkt_determinant),

		// matrix-vector operations
		t_mapFkts::value_type(T_STR"prod", fkt_product),

		// advanced linalg
		t_mapFkts::value_type(T_STR"eigenvecs", fkt_eigenvecs),
		t_mapFkts::value_type(T_STR"eigenvecs_sym", fkt_eigenvecs_sym),
		t_mapFkts::value_type(T_STR"qr", fkt_qr),


		// random numbers
		t_mapFkts::value_type(T_STR"rand01", fkt_rand01),
		t_mapFkts::value_type(T_STR"rand_real", fkt_rand_real),
		t_mapFkts::value_type(T_STR"rand_int", fkt_rand_int),
		t_mapFkts::value_type(T_STR"rand_norm", fkt_rand_norm),
		t_mapFkts::value_type(T_STR"rand_norm_nd", fkt_rand_norm_nd),

		// float classification
		t_mapFkts::value_type(T_STR"isnan", fkt_math_1arg_bret<std::isnan>),
		t_mapFkts::value_type(T_STR"isinf", fkt_math_1arg_bret<std::isinf>),
		t_mapFkts::value_type(T_STR"isfinite", fkt_math_1arg_bret<std::isfinite>),
	};

	add_ext_calls(mapFkts);
}
