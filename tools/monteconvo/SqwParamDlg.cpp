/*
 * S(q,w) parameters dialog
 * @author tweber
 * @date aug-2015
 * @license GPLv2
 */

#include "SqwParamDlg.h"


SqwParamDlg::SqwParamDlg(QWidget* pParent, QSettings* pSett)
	: QDialog(pParent), m_pSett(pSett)
{
	setupUi(this);
	if(m_pSett)
	{
		QFont font;
		if(m_pSett->contains("main/font_gen") && font.fromString(m_pSett->value("main/font_gen", "").toString()))
			setFont(font);
	}

	tableParams->verticalHeader()->setDefaultSectionSize(tableParams->verticalHeader()->minimumSectionSize()+2);

	connect(buttonBox, SIGNAL(clicked(QAbstractButton*)), this, SLOT(ButtonBoxClicked(QAbstractButton*)));


	if(m_pSett && m_pSett->contains("monteconvo/param_geo"))
		restoreGeometry(m_pSett->value("monteconvo/param_geo").toByteArray());
}

SqwParamDlg::~SqwParamDlg()
{}


void SqwParamDlg::SqwLoaded(const std::vector<SqwBase::t_var>& vecVars)
{
	const bool bSortTable = tableParams->isSortingEnabled();
	tableParams->setSortingEnabled(0);

	tableParams->setRowCount(vecVars.size());
	tableParams->setColumnWidth(0, 100);
	tableParams->setColumnWidth(1, 75);
	tableParams->setColumnWidth(2, 175);

	int iRow=0;
	for(const SqwBase::t_var& var : vecVars)
	{
		const std::string& strName = std::get<0>(var);
		const std::string& strType = std::get<1>(var);
		const std::string& strVal = std::get<2>(var);


		QTableWidgetItem *pItemName = tableParams->item(iRow, 0);
		if(!pItemName)
		{
			pItemName = new QTableWidgetItem();
			tableParams->setItem(iRow, 0, pItemName);
		}
		pItemName->setFlags(pItemName->flags() & ~Qt::ItemIsEditable);
		pItemName->setText(strName.c_str());


		QTableWidgetItem *pItemType = tableParams->item(iRow, 1);
		if(!pItemType)
		{
			pItemType = new QTableWidgetItem();
			tableParams->setItem(iRow, 1, pItemType);
		}
		pItemType->setFlags(pItemType->flags() & ~Qt::ItemIsEditable);
		pItemType->setText(strType.c_str());


		QTableWidgetItem *pItemVal = tableParams->item(iRow, 2);
		if(!pItemVal)
		{
			pItemVal = new QTableWidgetItem();
			tableParams->setItem(iRow, 2, pItemVal);
		}
		pItemVal->setFlags(pItemVal->flags() | Qt::ItemIsEditable);
		pItemVal->setText(strVal.c_str());

		++iRow;
	}

	tableParams->setSortingEnabled(bSortTable);
}

void SqwParamDlg::SaveSqwParams()
{
	std::vector<SqwBase::t_var> vecVars;

	for(int iRow=0; iRow<tableParams->rowCount(); ++iRow)
	{
		SqwBase::t_var var;
		std::get<0>(var) = tableParams->item(iRow, 0)->text().toStdString();
		std::get<2>(var) = tableParams->item(iRow, 2)->text().toStdString();

		vecVars.push_back(std::move(var));
	}

	emit SqwParamsChanged(vecVars);
}



void SqwParamDlg::showEvent(QShowEvent *pEvt)
{
	QDialog::showEvent(pEvt);
}

void SqwParamDlg::ButtonBoxClicked(QAbstractButton *pBtn)
{
	if(buttonBox->buttonRole(pBtn) == QDialogButtonBox::ApplyRole ||
	   buttonBox->buttonRole(pBtn) == QDialogButtonBox::AcceptRole)
	{
		SaveSqwParams();
	}
	else if(buttonBox->buttonRole(pBtn) == QDialogButtonBox::RejectRole)
		QDialog::reject();

	if(buttonBox->buttonRole(pBtn) == QDialogButtonBox::AcceptRole)
	{
		if(m_pSett)
			m_pSett->setValue("monteconvo/param_geo", saveGeometry());

		QDialog::accept();
	}
}


#include "SqwParamDlg.moc"
