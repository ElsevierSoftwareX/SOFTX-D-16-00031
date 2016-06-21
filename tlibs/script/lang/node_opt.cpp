/*
 * Optimizer
 * @author tweber
 * @date 27-sep-2014
 * @license GPLv2 or GPLv3
 */

#include "node.h"
#include "log/log.h"

static bool is_symbol_node(const Node *pNode)
{
	NodeType ty = pNode->GetType();

	if(ty==NODE_DOUBLE || ty==NODE_INT || ty==NODE_STRING)
		return 1;

	return 0;
}

static Symbol* get_symbol(Node *pNode)
{
	if(pNode->GetType() == NODE_INT)
		return ((NodeInt*)pNode)->GetSym();
	else if(pNode->GetType() == NODE_DOUBLE)
		return ((NodeReal*)pNode)->GetSym();
	else if(pNode->GetType() == NODE_STRING)
		return ((NodeString*)pNode)->GetSym();

	return 0;
}

static Node* new_node_for_sym(Symbol* pSym)
{
	if(pSym->GetType() == SYMBOL_DOUBLE)
		return new NodeReal((SymbolReal*)pSym);
	if(pSym->GetType() == SYMBOL_INT)
		return new NodeInt((SymbolInt*)pSym);
	if(pSym->GetType() == SYMBOL_STRING)
		return new NodeString((SymbolString*)pSym);

	return 0;
}

// is optimizable operator?
static bool is_operator(NodeType ty)
{
	if(ty==NODE_PLUS) return true;
	if(ty==NODE_MINUS) return true;
	if(ty==NODE_MULT) return true;
	if(ty==NODE_DIV) return true;
	if(ty==NODE_POW) return true;
	if(ty==NODE_MOD) return true;

	return false;
}



Node* NodeBreak::optimize() { return this; }
Node* NodeContinue::optimize() { return this; }
Node* NodeReal::optimize() { return this; }
Node* NodeInt::optimize() { return this; }
Node* NodeString::optimize() { return this; }


Node* NodeReturn::optimize()
{
	if(m_pExpr) m_pExpr = m_pExpr->optimize();

	return this;
}

Node* NodeIf::optimize()
{
	if(m_pExpr) m_pExpr = m_pExpr->optimize();
	if(m_pIf) m_pIf = m_pIf->optimize();
	if(m_pElse) m_pElse = m_pElse->optimize();

	return this;
}

Node* NodeWhile::optimize()
{
	if(m_pExpr) m_pExpr = m_pExpr->optimize();
	if(m_pStmt) m_pStmt = m_pStmt->optimize();

	return this;
}

Node* NodeFor::optimize()
{
	if(m_pExprInit) m_pExprInit = m_pExprInit->optimize();
	if(m_pExprCond) m_pExprCond = m_pExprCond->optimize();
	if(m_pExprEnd) m_pExprEnd = m_pExprEnd->optimize();
	if(m_pStmt) m_pStmt = m_pStmt->optimize();

	return this;
}

Node* NodeRangedFor::optimize()
{
	if(m_pIdent) m_pIdent = m_pIdent->optimize();
	if(m_pExpr) m_pExpr = m_pExpr->optimize();
	if(m_pStmt) m_pStmt = m_pStmt->optimize();

	return this;
}



Node* NodeIdent::optimize()
{
	if(m_pDefArg) m_pDefArg = m_pDefArg->optimize();

	return this;
}



Node* NodeArray::optimize()
{
	if(m_pArr) m_pArr = m_pArr->optimize();

	return this;
}

Node* NodeMap::optimize()
{
	if(m_pMap) m_pMap = m_pMap->optimize();

	return this;
}



Node* NodePair::optimize()
{
	if(m_pFirst) m_pFirst = m_pFirst->optimize();
	if(m_pSecond) m_pSecond = m_pSecond->optimize();

	return this;
}

Node* NodeRange::optimize()
{
	if(m_pBegin) m_pBegin = m_pBegin->optimize();
	if(m_pEnd) m_pEnd = m_pEnd->optimize();

	return this;
}



Node* NodeArrayAccess::optimize()
{
	if(m_pIdent) m_pIdent = m_pIdent->optimize();
	if(m_pExpr)
	{
		m_pExpr = m_pExpr->optimize();

		NodeBinaryOp* pIndices = (NodeBinaryOp*) m_pExpr;
		m_vecIndices = pIndices->flatten(NODE_ARGS);
	}

	return this;
}



Node* NodeUnaryOp::optimize()
{
	if(!m_pChild) return this;
	m_pChild = m_pChild->optimize();

	if(m_type==NODE_UMINUS)
	{
		Node *pChild = m_pChild;
		const NodeType ty = pChild->GetType();
		bool bSubstNode = 0;

		if(ty == NODE_DOUBLE)
		{
			SymbolReal *pSym = ((NodeReal*)pChild)->GetSym();
			pSym->GetVal() = -pSym->GetVal();
			bSubstNode = 1;
		}
		else if(ty == NODE_INT)
		{
			SymbolInt *pSym = ((NodeInt*)pChild)->GetSym();
			pSym->GetVal() = -pSym->GetVal();
			bSubstNode = 1;
		}


		if(bSubstNode)
		{
			this->m_pChild = 0;
			delete this;
			return pChild;
		}
	}

	return this;
}

Node* NodeBinaryOp::optimize()
{
	if(m_pLeft) m_pLeft = m_pLeft->optimize();
	if(m_pRight) m_pRight = m_pRight->optimize();

	if(GetType()==NODE_STMTS || GetType()==NODE_FUNCS)
		FlattenNodes(GetType());

	if(!m_pLeft || !m_pRight) return this;

	if(is_symbol_node(m_pLeft) && is_symbol_node(m_pRight) && is_operator(GetType()))
	{
		Symbol *pSymLeft = get_symbol(m_pLeft);
		Symbol *pSymRight = get_symbol(m_pRight);

		Symbol *pNewSym = Op(pSymLeft, pSymRight, GetType(), 0);
		pNewSym->SetConst(1);

		delete this;
		return new_node_for_sym(pNewSym);
	}

	return this;
}



Node* NodeCall::optimize()
{
	if(m_pIdent) m_pIdent = m_pIdent->optimize();
	if(m_pArgs)
	{
		m_pArgs = m_pArgs->optimize();
		m_vecArgs = ((NodeBinaryOp*)m_pArgs)->flatten(NODE_ARGS);
	}

	return this;
}

Node* NodeFunction::optimize()
{
	if(m_pIdent) m_pIdent = m_pIdent->optimize();

	if(m_pArgs)
	{
		m_pArgs = m_pArgs->optimize();

		if(m_pArgs->GetType()==NODE_IDENTS || m_pArgs->GetType()==NODE_IDENT)
		{
			NodeBinaryOp* pArgs = (NodeBinaryOp*)m_pArgs;
			m_vecArgs = pArgs->flatten(NODE_IDENTS);
		}
	}

	if(m_pStmts) m_pStmts = m_pStmts->optimize();

	return this;
}


//--------------------------------------------------------------------------------

void NodeBinaryOp::FlattenNodes(NodeType ntype)
{
	if(GetType() == ntype)
	{
		m_vecNodesFlat = this->flatten(ntype);
	}
}

// TODO: Fix: Gets called for non casted NodeBinaryOps which are of other type!
std::vector<Node*> NodeBinaryOp::flatten(NodeType ntype) const
{
	//log_debug(GetType(), ", ", ntype);
	std::vector<Node*> vecNodes;

	if(GetType() == ntype)
	{
		NodeBinaryOp *pLeft = (NodeBinaryOp*) m_pLeft;
		NodeBinaryOp *pRight = (NodeBinaryOp*) m_pRight;

		if(m_pLeft)
		{
			if(m_pLeft->GetType() == ntype)
			{
				std::vector<Node*> vecLeft = pLeft->flatten(ntype);
				vecNodes.insert(vecNodes.begin(), vecLeft.begin(), vecLeft.end());
			}
			else
			{
				vecNodes.push_back(m_pLeft);
			}
		}

		if(m_pRight)
		{
			if(m_pRight->GetType() == ntype)
			{
				std::vector<Node*> vecRight = pRight->flatten(ntype);
				vecNodes.insert(vecNodes.end(), vecRight.begin(), vecRight.end());
			}
			else
			{
				vecNodes.push_back(m_pRight);
			}
		}
	}
	else
	{
		vecNodes.push_back(const_cast<NodeBinaryOp*>(this));
	}

	return vecNodes;
}
