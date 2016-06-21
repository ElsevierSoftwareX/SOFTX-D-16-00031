/*
 * Symbol Table
 * @author tweber
 * @date 2013
 * @license GPLv2 or GPLv3
 */

#ifndef __HERMELIN_SYM__
#define __HERMELIN_SYM__

#include "types.h"
//#include "info.h"
struct ParseInfo;

#include "helper/exception.h"
#include "log/log.h"

#include <map>
#include <unordered_map>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <initializer_list>
#include <functional>
#include <boost/functional/hash.hpp>



enum SymbolType : unsigned int
{
	SYMBOL_DOUBLE = 1<<0,
	SYMBOL_INT = 1<<1,
	SYMBOL_STRING = 1<<2,

	SYMBOL_ARRAY = 1<<3,
	SYMBOL_MAP = 1<<4,

	SYMBOL_COMPLEX = 1<<5,

	SYMBOL_SCALAR = SYMBOL_DOUBLE | SYMBOL_INT,
	SYMBOL_CONTAINER = SYMBOL_ARRAY | SYMBOL_MAP | SYMBOL_STRING,
	SYMBOL_ANY = 0xffffffff
};


class Symbol;
class SymbolArray;
class SymbolMap;
class SymbolString;
class SymbolInt;
class SymbolReal;
class SymbolComplex;


struct SymbolMapKey
{
	std::size_t key;	// the actual hash
	t_string strKey;	// key as string
	SymbolType tyKey;	// key type

	SymbolMapKey(const t_string& _str);
	SymbolMapKey(t_string&& _str);

	SymbolMapKey(const Symbol* pSym);

	SymbolMapKey() = default;
};



class Symbol
{
protected:
	t_string m_strName;
	t_string m_strIdent;			// last seen identifier
	bool m_bIsRval = 1;			// is the symbol an rvalue?
	bool m_bIsConst = 0;			// is the symbol a constant? -> forces clone

	unsigned int m_iArrIdx = 0;		// if symbol is contained in an array
	SymbolArray *m_pArr = 0;

	SymbolMapKey m_MapKey;			// if symbol is contained in a map
	SymbolMap *m_pMap = 0;

public:
	Symbol() = default;
	virtual ~Symbol() = 0;

	virtual SymbolType GetType() const = 0;
	virtual t_string GetTypeName() const = 0;

	bool IsRval() const { return m_bIsRval; }
	void SetRval(bool bRval) { m_bIsRval = bRval; }

	bool IsConst() const { return m_bIsConst; }
	void SetConst(bool bConst) { m_bIsConst = bConst; }

	virtual std::size_t hash() const = 0;

	virtual const t_string& GetName() const { return m_strName; }
	virtual const t_string& GetIdent() const { return m_strIdent; }
	virtual void SetName(const t_string& strName) { m_strName = strName; }
	virtual void SetIdent(const t_string& strIdent) { m_strIdent = strIdent; }

	// cast and clone symbol
	virtual Symbol* ToType(SymbolType stype) const = 0;

	virtual t_string print() const = 0;
	virtual Symbol* clone() const = 0;
	virtual Symbol* alloc() const = 0;		// alloc new symbol of same type
	virtual void assign(Symbol *pSym) = 0;

	//virtual bool equals(Symbol *pSym) const = 0;
	virtual bool IsNotZero() const = 0;
	virtual bool IsLessThan(const Symbol&) const { return 0; }
	virtual bool IsGreaterThan(const Symbol&) const { return 0; }

	// cast value without converting or cloning symbol
	virtual t_int GetValInt() const { return 0; }
	virtual t_real GetValDouble() const { return 0.; }
	virtual t_complex GetValComplex() const { return t_complex(0.,0.); }

	virtual bool IsScalar() const { return 0; }

public:
	SymbolMap* GetMapPtr() const { return m_pMap; }
	void SetMapPtr(SymbolMap* pMap) { m_pMap = pMap; }

	const SymbolMapKey& GetMapKey() const { return m_MapKey; }
	void SetMapKey(const SymbolMapKey& key) { m_MapKey = key; }
	void SetMapKey(SymbolMapKey&& key) { m_MapKey = key; }


	SymbolArray* GetArrPtr() const { return m_pArr; }
	void SetArrPtr(SymbolArray* pArr) { m_pArr = pArr; }

	unsigned int GetArrIdx() const { return m_iArrIdx; }
	void SetArrIdx(unsigned int iIdx) { m_iArrIdx = iIdx; }

	void ClearIndices()
	{
		m_pMap = 0;
		m_pArr = 0;
	}
};

class SymbolReal : public Symbol
{
protected:
	t_real m_dVal;
	static const int m_defprec;
	static int m_prec;

	static /*std::*/boost::hash<t_real> s_hsh;

public:
	SymbolReal() : Symbol(), m_dVal(0.) {}
	SymbolReal(t_real dVal) : m_dVal(dVal) {}
	SymbolReal(const t_string&) { throw tl::Err("Invalid SymbolReal constructor."); }

	virtual SymbolType GetType() const override { return SYMBOL_DOUBLE; }
	virtual t_string GetTypeName() const override { return T_STR"real"; }
	virtual Symbol* ToType(SymbolType stype) const override;

	virtual std::size_t hash() const override { return s_hsh(m_dVal); }

	virtual t_string print() const override;
	virtual Symbol* clone() const override;
	virtual Symbol* alloc() const override { return new SymbolReal(); }
	virtual void assign(Symbol *pSym) override;

	//virtual bool equals(Symbol *pSym) const;
	virtual bool IsLessThan(const Symbol&) const override;
	virtual bool IsGreaterThan(const Symbol&) const override;
	virtual bool IsNotZero() const override { return m_dVal != 0.; }

	virtual t_int GetValInt() const override { return t_int(m_dVal); }
	virtual t_real GetValDouble() const override { return m_dVal; }
	virtual t_complex GetValComplex() const override { return t_complex(m_dVal, 0.); }

	static const int GetDefPrec() { return m_defprec; }
	static const int GetPrec() { return m_prec; }
	static void SetPrec(int iPrec) { m_prec = iPrec; }

	void SetVal(t_real dVal) { m_dVal = dVal; }
	const t_real& GetVal() const { return m_dVal; }
	t_real& GetVal() { return m_dVal; }

	virtual bool IsScalar() const override { return 1; }
};

class SymbolInt : public Symbol
{
protected:
	t_int m_iVal;

	static /*std::*/boost::hash<t_int> s_hsh;

public:
	SymbolInt() : Symbol(), m_iVal(0) {}
	SymbolInt(t_int iVal) : m_iVal(iVal) {}

	virtual SymbolType GetType() const override { return SYMBOL_INT; }
	virtual t_string GetTypeName() const override { return T_STR"int"; }
	virtual Symbol* ToType(SymbolType stype) const override;

	virtual std::size_t hash() const override { return s_hsh(m_iVal); }

	virtual t_string print() const override;
	virtual Symbol* clone() const override;
	virtual Symbol* alloc() const override { return new SymbolInt(); }
	virtual void assign(Symbol *pSym) override;

	//virtual bool equals(Symbol *pSym) const;
	virtual bool IsLessThan(const Symbol&) const override;
	virtual bool IsGreaterThan(const Symbol&) const override;
	virtual bool IsNotZero() const override { return m_iVal != 0; }

	virtual t_int GetValInt() const override { return m_iVal; }
	virtual t_real GetValDouble() const override { return t_real(m_iVal); }
	virtual t_complex GetValComplex() const override { return t_complex(t_real(m_iVal), 0.); }

	void SetVal(int iVal) { m_iVal = iVal; }
	const t_int& GetVal() const { return m_iVal; }
	t_int& GetVal() { return m_iVal; }

	virtual bool IsScalar() const override { return 1; }
};

class SymbolString : public Symbol
{
protected:
	t_string m_strVal;

	static /*std::*/boost::hash<std::string> s_hsh;

public:
	SymbolString() : Symbol() {}
	SymbolString(const t_char* pcStr) : m_strVal(pcStr) {}
	template<typename _t_string = t_string>
	SymbolString(_t_string&& str) : m_strVal(std::forward<_t_string>(str)) {}
	SymbolString(t_real dVal) { throw tl::Err("Invalid SymbolString constructor."); }

	virtual SymbolType GetType() const override { return SYMBOL_STRING; }
	virtual t_string GetTypeName() const override { return T_STR"string"; }
	virtual Symbol* ToType(SymbolType stype) const override;

	virtual std::size_t hash() const override { return s_hsh(m_strVal); }
	static std::size_t hash(const t_string& str) { return s_hsh(str); }

	//virtual bool equals(Symbol *pSym) const;
	virtual bool IsLessThan(const Symbol&) const override;
	virtual bool IsGreaterThan(const Symbol&) const override;

	virtual t_string print() const override;
	virtual Symbol* clone() const override;
	virtual Symbol* alloc() const override { return new SymbolString(); }
	virtual void assign(Symbol *pSym) override;

	void SetVal(const t_string& str) { m_strVal = str; }
	const t_string& GetVal() const { return m_strVal; }
	t_string& GetVal() { return m_strVal; }

	virtual bool IsNotZero() const override { return 0; }
};

class SymbolComplex : public Symbol
{
protected:
	t_complex m_val;

	static const int m_defprec;
	static int m_prec;

	static /*std::*/boost::hash<t_complex> s_hsh;

public:
	SymbolComplex() : Symbol(), m_val(0., 0.) {}
	SymbolComplex(t_real dReal, t_real dImag) : m_val(dReal, dImag) {}
	SymbolComplex(const t_complex& val) : m_val(val) {}
	SymbolComplex(t_complex&& val) : m_val(val) {}
	SymbolComplex(const t_string&) { throw tl::Err("Invalid SymbolComplex constructor."); }

	virtual SymbolType GetType() const override { return SYMBOL_COMPLEX; }
	virtual t_string GetTypeName() const override { return T_STR"complex"; }
	virtual Symbol* ToType(SymbolType stype) const override;

	virtual std::size_t hash() const override { return s_hsh(m_val); }

	virtual t_string print() const override;
	virtual Symbol* clone() const override;
	virtual Symbol* alloc() const override { return new SymbolComplex(); }
	virtual void assign(Symbol *pSym) override;

	//virtual bool equals(Symbol *pSym) const;
	virtual bool IsLessThan(const Symbol&) const override;
	virtual bool IsGreaterThan(const Symbol&) const override;
	virtual bool IsNotZero() const override { return m_val != t_complex(0.,0.); }

	virtual t_int GetValInt() const override { throw tl::Err("Cannot convert complex to scalar type."); /*return t_int(GetValReal());*/ }
	virtual t_real GetValDouble() const override { throw tl::Err("Cannot convert complex to scalar type."); /*return GetValReal();*/ }
	virtual t_complex GetValComplex() const override { return m_val; }

	static const int GetDefPrec() { return m_defprec; }
	static const int GetPrec() { return m_prec; }
	static void SetPrec(int iPrec) { m_prec = iPrec; }

	void SetVal(t_real dReal, t_real dImag) { m_val=t_complex(dReal, dImag); }
	void SetVal(const t_complex& val) { m_val = val; }
	const t_complex& GetVal() const { return m_val; }

	void SetValReal(t_real dVal) { m_val.real(dVal); }
	void SetValImag(t_real dVal) { m_val.imag(dVal); }
	t_real GetValReal() const { return m_val.real(); }
	t_real GetValImag() const { return m_val.imag(); }

	virtual bool IsScalar() const override { return 0; }
};


class SymbolArray : public Symbol
{
public:
	using t_arr = std::vector<Symbol*>;

protected:
	bool m_bDontDel;
	t_arr m_arr;

public:
	SymbolArray() : Symbol(), m_bDontDel(0) { /*std::cout << "symarr -> new" << std::endl;*/ }
	SymbolArray(const std::initializer_list<Symbol*>& lst);
	virtual ~SymbolArray();

	virtual SymbolType GetType() const override { return SYMBOL_ARRAY; }
	virtual t_string GetTypeName() const override { return T_STR"vector"; }
	virtual Symbol* ToType(SymbolType stype) const override;

	virtual std::size_t hash() const override;

	std::vector<t_real> ToDoubleArray() const;
	void FromDoubleArray(const std::vector<t_real>& vec);

	virtual t_string print() const override;
	virtual Symbol* clone() const override;
	virtual Symbol* alloc() const override { return new SymbolArray(); }
	virtual void assign(Symbol *pSym) override;
	//virtual bool equals(Symbol *pSym) const;

	virtual bool IsNotZero() const override { return 0; }

	void UpdateIndex(unsigned int, bool bOverwrite=true);
	void UpdateLastNIndices(unsigned int N, bool bOverwrite=true);
	void UpdateIndices(bool bOverwrite=true);
	void ClearIndices();

	const std::vector<Symbol*>& GetArr() const { return m_arr; }
	std::vector<Symbol*>& GetArr() { return m_arr; }

	virtual bool IsScalar() const override { return 0; }


	bool GetDontDel() const { return m_bDontDel; }
	void SetDontDel(bool b) { m_bDontDel = b; }
};



struct SymbolMapKeyHash
{
	std::size_t operator()(const SymbolMapKey& key) const { return key.key; }
};

struct SymbolMapKeyEqual
{
	bool operator()(const SymbolMapKey& key0, const SymbolMapKey& key1) const
	{ return key0.key==key1.key && key0.strKey==key1.strKey && key0.tyKey==key1.tyKey; }
};

class SymbolMap : public Symbol
{
public:
	typedef std::unordered_map<SymbolMapKey, Symbol*, SymbolMapKeyHash, SymbolMapKeyEqual> t_map;
protected:
	t_map m_map;

public:
	SymbolMap() : Symbol() {}
	virtual ~SymbolMap();

	// no copy/move constructors -> use clone
	SymbolMap(const SymbolMap& map) = delete;
	SymbolMap(SymbolMap&& map) = delete;

	virtual SymbolType GetType() const override { return SYMBOL_MAP; }
	virtual t_string GetTypeName() const override { return T_STR"map"; }
	virtual Symbol* ToType(SymbolType stype) const override;

	virtual std::size_t hash() const override;

	virtual t_string print() const override;
	virtual Symbol* clone() const override;
	virtual Symbol* alloc() const override { return new SymbolMap(); }
	virtual void assign(Symbol *pSym) override;
	//virtual bool equals(Symbol *pSym) const;

	virtual bool IsNotZero() const override { return 0; }

	void UpdateIndex(const t_map::key_type& strKey);
	void UpdateIndices();

	t_string GetStringVal(const SymbolMapKey& key, bool *pbHasVal=0) const;
	t_int GetIntVal(const SymbolMapKey& key, bool *pbHasVal=0) const;
	t_real GetRealVal(const SymbolMapKey& key, bool *pbHasVal=0) const;

	t_string GetStringVal(const t_string& strKey, bool *pbHasVal=0) const
	{ return GetStringVal(SymbolMapKey(strKey), pbHasVal); }
	t_int GetIntVal(const t_string& strKey, bool *pbHasVal=0) const
	{ return GetIntVal(SymbolMapKey(strKey), pbHasVal); }
	t_real GetRealVal(const t_string& strKey, bool *pbHasVal=0) const
	{ return GetRealVal(SymbolMapKey(strKey), pbHasVal); }

	virtual bool IsScalar() const override { return 0; }

	const t_map& GetMap() const { return m_map; }
	t_map& GetMap() { return m_map; }
};



// --------------------------------------------------------------------------------


class SymbolTable
{
protected:
	typedef std::unordered_map<t_string, Symbol*> t_syms;
	t_syms m_syms;

public:
	SymbolTable();
	virtual ~SymbolTable();

	void print() const;

	Symbol* GetSymbol(const t_string& strKey);
	void InsertSymbol(const t_string& strKey, Symbol *pSym);
	void RemoveSymbol(const t_string& strKey);
	void RemoveSymbolNoDelete(const t_string& strKey);
	bool IsPtrInMap(const Symbol* pSym) const;
};



// --------------------------------------------------------------------------------
// conversions

#include <type_traits>

template<typename T> static T convert_symbol(const Symbol* pSym)
{ std::cerr << "Error: Invalid symbol conversion." << std::endl; return T(); }

template<> t_string convert_symbol<t_string>(const Symbol* pSym)
{ return pSym->print(); }
template<>  t_real convert_symbol<t_real>(const Symbol* pSym)
{ return pSym->GetValDouble(); }
template<> t_int convert_symbol<t_int>(const Symbol* pSym)
{ return pSym->GetValInt(); }


template<typename T> static Symbol* create_symbol(const T& t)
{ G_CERR << "Error: Invalid symbol creation." << std::endl; return 0; }

template<> Symbol* create_symbol<t_real>(const t_real& t)
{ return new SymbolReal(t); }
template<> Symbol* create_symbol<t_int>(const t_int& t)
{ return new SymbolInt(t); }
template<> Symbol* create_symbol<t_string>(const t_string& t)
{ return new SymbolString(t); }


template<typename T1=t_string, typename T2=t_real>
static std::map<T1, T2> sym_to_map(const Symbol* pSym)
{
	if(pSym->GetType() != SYMBOL_MAP)
		return std::map<T1, T2>();

	SymbolMap* pSymMap = (SymbolMap*)pSym;
	std::map<T1, T2> _map;

	for(const typename SymbolMap::t_map::value_type& pair : pSymMap->GetMap())
		_map[pair.first.strKey] = convert_symbol<T2>(pair.second);

	return _map;
}

extern bool is_vec(const Symbol* pSym);
extern bool is_mat(const Symbol* pSym, unsigned int *piNumCols=0, unsigned int *piNumRows=0);

template<template<class> class t_vec, typename T=t_real>
static t_vec<T> sym_to_vec(const Symbol* pSym)
{
	if(!pSym || pSym->GetType() != SYMBOL_ARRAY)
		return t_vec<T>();

	SymbolArray* pSymArr = (SymbolArray*)pSym;
	t_vec<T> vec(pSymArr->GetArr().size());

	unsigned int iIdx = 0;
	for(const Symbol* pSymInArr : pSymArr->GetArr())
	{
		vec[iIdx] = convert_symbol<T>(pSymInArr);
		++iIdx;
	}

	return vec;
}

template<template<class> class t_vec, typename T=t_real>
static Symbol* vec_to_sym(const t_vec<T>& vec)
{
	SymbolArray* pSym = new SymbolArray();
	pSym->GetArr().reserve(vec.size());

	for(const T& t : vec)
		pSym->GetArr().push_back(create_symbol<T>(t));

	pSym->UpdateIndices();
	return pSym;
}

template<template<class> class t_mat, template<class> class t_vec, typename T=t_real>
static t_mat<T> sym_to_mat(const Symbol* pSym, bool* pbIsMat=0)
{
	unsigned int iNumCols=0, iNumRows=0;
	if(!is_mat(pSym, &iNumCols, &iNumRows))
	{
		if(pbIsMat) *pbIsMat = 0;
		return t_mat<T>();
	}
	if(pbIsMat) *pbIsMat = 1;

	t_mat<T> mat(iNumRows, iNumCols);
	const SymbolArray* pSymArr = (SymbolArray*)pSym;

	unsigned int iRow=0;
	for(const Symbol* pSymInArr : pSymArr->GetArr())
	{
		t_vec<T> vecRow = sym_to_vec<t_vec>(pSymInArr);
		unsigned int iNumActCols = std::min<unsigned int>(vecRow.size(), iNumCols);

		for(unsigned int iCol=0; iCol<iNumActCols; ++iCol)
			mat(iRow, iCol) = vecRow[iCol];

		// fill rest with 0
		for(unsigned int iCol=iNumActCols; iCol<iNumCols; ++iCol)
			mat(iRow, iCol) = 0.;

		++iRow;
	}

	return mat;
}

template<template<class> class t_mat, typename T=t_real>
static Symbol* mat_to_sym(const t_mat<T>& mat)
{
	unsigned int iNumRows = mat.size1();
	unsigned int iNumCols = mat.size2();

	SymbolArray* pSym = new SymbolArray();
	pSym->GetArr().reserve(iNumRows);

	for(unsigned int iRow=0; iRow<iNumRows; ++iRow)
	{
		SymbolArray* pRow = new SymbolArray();
		pRow->GetArr().reserve(iNumCols);

		for(unsigned int iCol=0; iCol<iNumCols; ++iCol)
		{
			Symbol *pSymVal = create_symbol<T>(mat(iRow, iCol));
			pRow->GetArr().push_back(pSymVal);
		}

		pRow->UpdateIndices();
		pSym->GetArr().push_back(pRow);
	}

	pSym->UpdateIndices();
	return pSym;
}


extern void safe_delete(Symbol *&pSym, const SymbolTable* pSymTab, ParseInfo* pParseInfo=0);
extern bool is_tmp_sym(const Symbol* pSym);
extern bool clone_if_needed(Symbol *pSym, Symbol*& pClone);
extern Symbol* recycle_or_alloc(const std::initializer_list<const Symbol*>& lstSyms, bool bAlwaysAlloc=0);

#endif
