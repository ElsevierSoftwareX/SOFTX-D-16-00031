/*
 * Simple Script
 * @author tweber
 * @date 2013
 * @license GPLv2 or GPLv3
 */

#include "types.h"
#include "log/log.h"

#include <iostream>
#include <string.h>

#include "tokens.h"
#include "yylexer.h"


extern "C" int yylex(void* _yylval, void* _pParseObj)
{
	if(_yylval==0 || _pParseObj==0)
	{
		tl::log_err("Invalid lval or parseobj in lexer.");
		return 0;
	}

	YYSTYPE& yylval = *(YYSTYPE*)_yylval;
	ParseObj* pParseObj = (ParseObj*)_pParseObj;
	Lexer* pLexer = pParseObj->pLexer;
	//pLexer->print();

	const Token& tok = pLexer->lex();
	pParseObj->iCurLine = tok.iLine;

	if(tok.type == LEX_TOKEN_END)
		return 0;
	else if(tok.type == LEX_TOKEN_INVALID)
		return 1;
	else if(tok.type == LEX_TOKEN_DOUBLE)
	{
		yylval.pNode = new NodeReal(tok.dVal);
		return TOK_DOUBLE;
	}
	else if(tok.type == LEX_TOKEN_INT)
	{
		yylval.pNode = new NodeInt(tok.iVal);
		return TOK_INT;
	}
	else if(tok.type == LEX_TOKEN_STRING)
	{
		yylval.pNode = new NodeString(tok.strVal);
		return TOK_STRING;
	}
	else if(tok.type == LEX_TOKEN_CHAROP)
		return tok.cOp;
	else if(tok.type == LEX_TOKEN_IDENT)
	{
		yylval.pNode = new NodeIdent(tok.strVal);
		return TOK_IDENT;
	}
	else
	{
		switch(tok.type)
		{
			case LEX_TOKEN_IF: 		return TOK_IF;
			case LEX_TOKEN_ELSE:		return TOK_ELSE;
			case LEX_TOKEN_FOR:		return TOK_FOR;
			case LEX_TOKEN_WHILE:		return TOK_WHILE;
			case LEX_TOKEN_RETURN:		return TOK_RETURN;
			case LEX_TOKEN_BREAK:		return TOK_BREAK;
			case LEX_TOKEN_CONTINUE:	return TOK_CONTINUE;
			case LEX_TOKEN_LOG_AND:		return TOK_LOG_AND;
			case LEX_TOKEN_LOG_OR:		return TOK_LOG_OR;
			case LEX_TOKEN_LOG_NOT:		return TOK_LOG_NOT;
			case LEX_TOKEN_LOG_EQ:		return TOK_LOG_EQ;
			case LEX_TOKEN_LOG_NEQ:		return TOK_LOG_NEQ;
			case LEX_TOKEN_LOG_LESS:	return TOK_LOG_LESS;
			case LEX_TOKEN_LOG_GREATER:	return TOK_LOG_GREATER;
			case LEX_TOKEN_LOG_LEQ:		return TOK_LOG_LEQ;
			case LEX_TOKEN_LOG_GEQ:		return TOK_LOG_GEQ;
			case LEX_TOKEN_GLOBAL:		return TOK_GLOBAL;
			default: break;
		}
	}

	tl::log_err("line ", pParseObj->iCurLine,
		" in \"", pParseObj->strCurFile,
		"\": Invalid token: ", tok.type, 
		".");
	return 0;
}


extern "C" void yyerror(void *_pParseObj, const char* pc)
{
	ParseObj *pParseObj = (ParseObj*)_pParseObj;

	t_string strLine;
	if(pParseObj)
	{
		t_ostringstream ostrLine;
		ostrLine << "line " << pParseObj->iCurLine;
		ostrLine << " in \"" << pParseObj->strCurFile << "\"";
		strLine = ostrLine.str();
	}

	tl::log_err(strLine, ": ", pc, ".");
}
