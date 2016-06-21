/*
 * Lexer
 * @author tweber
 * @date 2013
 * @license GPLv2 or GPLv3
 */

#ifndef __MIEZE_LEXER__
#define __MIEZE_LEXER__

#include "types.h"

#include <istream>
#include <string>
#include <vector>

enum TokenType : unsigned int
{
	LEX_TOKEN_INVALID,
	LEX_TOKEN_END,

	LEX_TOKEN_DOUBLE,
	LEX_TOKEN_INT,
	LEX_TOKEN_STRING,
	LEX_TOKEN_IDENT,
	LEX_TOKEN_CHAROP,

	LEX_TOKEN_IF,
	LEX_TOKEN_ELSE,
	LEX_TOKEN_FOR,
	LEX_TOKEN_WHILE,

	LEX_TOKEN_RETURN,
	LEX_TOKEN_BREAK,
	LEX_TOKEN_CONTINUE,

	LEX_TOKEN_LOG_AND,
	LEX_TOKEN_LOG_OR,
	LEX_TOKEN_LOG_NOT,
	LEX_TOKEN_LOG_EQ,
	LEX_TOKEN_LOG_NEQ,
	LEX_TOKEN_LOG_LESS,
	LEX_TOKEN_LOG_GREATER,
	LEX_TOKEN_LOG_LEQ,
	LEX_TOKEN_LOG_GEQ,

	LEX_TOKEN_GLOBAL
};

struct Token
{
	TokenType type;

	t_char cOp;
	t_int iVal;
	t_real dVal;
	t_string strVal;

	unsigned int iLine;

	Token() : type(LEX_TOKEN_INVALID), cOp(0), dVal(0), iLine(0)
	{}
};

class Lexer
{
protected:
	bool m_bOk;
	t_string m_strWhitespace, m_strSep;

	unsigned int m_iLexPos;
	unsigned int m_iNumToks;
	std::vector<Token> m_vecToks;
	Token m_tokEnd;

	t_string m_strFile;

	void FixTokens();

public:
	Lexer();
	Lexer(const t_string& str, const t_char* pcFile=0);
	virtual ~Lexer();

	void load(const t_string& strInput);
	void print();
	const Token& lex();

	unsigned int GetNumTokens() const { return m_vecToks.size(); }
	const Token& GetToken(unsigned int i) const { return m_vecToks[i]; }

	static t_string RemoveComments(const t_string& strInput);
	static std::vector<t_string> GetStringTable(const t_string& strInput);
	static void ReplaceEscapes(t_string& str);

	bool IsOk() const { return m_bOk; }
};

#endif
