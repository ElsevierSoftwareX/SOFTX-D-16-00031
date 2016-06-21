/**
 * Atom Positions Dialog
 * @author Tobias Weber
 * @date nov-2015
 * @license GPLv2
 */

#include "AtomsDlg.h"
#include "tlibs/string/string.h"
#include "tlibs/math/linalg.h"

using t_real = t_real_glob;


enum class AtInfo : int
{
	NAME = 0,
	POS_X = 1,
	POS_Y = 2,
	POS_Z = 3,

	J = 4
};

AtomsDlg::AtomsDlg(QWidget* pParent, QSettings *pSettings, bool bEnableJ)
	: QDialog(pParent), m_pSettings(pSettings), m_bEnableJ(bEnableJ)
{
	setupUi(this);
	if(m_pSettings)
	{
		QFont font;
		if(m_pSettings->contains("main/font_gen") && font.fromString(m_pSettings->value("main/font_gen", "").toString()))
			setFont(font);
	}

	if(m_bEnableJ)
	{
		tableAtoms->setColumnCount(5);
		tableAtoms->setHorizontalHeaderItem(static_cast<int>(AtInfo::J), new QTableWidgetItem("J (meV/K)"));
	}

	tableAtoms->setColumnWidth(0, 75);
	btnAdd->setIcon(load_icon("res/list-add.svg"));
	btnDel->setIcon(load_icon("res/list-remove.svg"));

#if QT_VER >= 5
	QObject::connect(btnAdd, &QAbstractButton::clicked, this, &AtomsDlg::AddAtom);
	QObject::connect(btnDel, &QAbstractButton::clicked, this, &AtomsDlg::RemoveAtom);
	QObject::connect(buttonBox, &QDialogButtonBox::clicked, this, &AtomsDlg::ButtonBoxClicked);
#else
	QObject::connect(btnAdd, SIGNAL(clicked(bool)), this, SLOT(AddAtom()));
	QObject::connect(btnDel, SIGNAL(clicked(bool)), this, SLOT(RemoveAtom()));
	QObject::connect(buttonBox, SIGNAL(clicked(QAbstractButton*)), this,
		SLOT(ButtonBoxClicked(QAbstractButton*)));
#endif

	if(m_pSettings && m_pSettings->contains("atoms/geo"))
		restoreGeometry(m_pSettings->value("atoms/geo").toByteArray());
}

AtomsDlg::~AtomsDlg() {}


void AtomsDlg::RemoveAtom()
{
	const bool bSort = tableAtoms->isSortingEnabled();
	tableAtoms->setSortingEnabled(0);

	bool bNothingRemoved = 1;

	// remove selected rows
	QList<QTableWidgetSelectionRange> lstSel = tableAtoms->selectedRanges();
	for(QTableWidgetSelectionRange& range : lstSel)
	{
		for(int iRow=range.bottomRow(); iRow>=range.topRow(); --iRow)
		{
			tableAtoms->removeRow(iRow);
			bNothingRemoved = 0;
		}
	}

	// remove last row if nothing is selected
	if(bNothingRemoved)
		tableAtoms->removeRow(tableAtoms->rowCount()-1);

	tableAtoms->setSortingEnabled(bSort);
}

void AtomsDlg::AddAtom()
{
	const bool bSort = tableAtoms->isSortingEnabled();
	tableAtoms->setSortingEnabled(0);

	int iRow = tableAtoms->rowCount();
	tableAtoms->insertRow(iRow);
	tableAtoms->setItem(iRow, 0, new QTableWidgetItem("H"));
	for(unsigned int i=0; i<3; ++i)
		tableAtoms->setItem(iRow, static_cast<int>(AtInfo::POS_X)+i, new QTableWidgetItem("0"));
	if(m_bEnableJ)
		tableAtoms->setItem(iRow, static_cast<int>(AtInfo::J), new QTableWidgetItem("0"));

	tableAtoms->setSortingEnabled(bSort);
}


void AtomsDlg::SetAtoms(const std::vector<AtomPos<t_real>>& vecAtoms)
{
	const bool bSort = tableAtoms->isSortingEnabled();
	tableAtoms->setSortingEnabled(0);

	tableAtoms->setRowCount(vecAtoms.size());

	for(std::size_t iRow=0; iRow<vecAtoms.size(); ++iRow)
	{
		// add missing items
		for(int iCol=0; iCol<tableAtoms->columnCount(); ++iCol)
			if(!tableAtoms->item(iRow, iCol))
				tableAtoms->setItem(iRow, iCol, new QTableWidgetItem(""));

		const AtomPos<t_real>& atom = vecAtoms[iRow];
		tableAtoms->item(iRow, 0)->setText(atom.strAtomName.c_str());
		for(unsigned int i=0; i<3; ++i)
			tableAtoms->item(iRow, static_cast<int>(AtInfo::POS_X)+i)->setText(tl::var_to_str(atom.vecPos[i]).c_str());

		if(m_bEnableJ)
			tableAtoms->item(iRow, static_cast<int>(AtInfo::J))->setText(tl::var_to_str(atom.J).c_str());
	}

	tableAtoms->setSortingEnabled(bSort);
}

void AtomsDlg::SendApplyAtoms()
{
	std::vector<AtomPos<t_real>> vecAtoms;
	vecAtoms.reserve(tableAtoms->rowCount());

	for(int iRow=0; iRow<tableAtoms->rowCount(); ++iRow)
	{
		AtomPos<t_real> atom;
		atom.strAtomName = tableAtoms->item(iRow, static_cast<int>(AtInfo::NAME))->text().toStdString();
		tl::trim(atom.strAtomName);
		t_real dX = tl::str_to_var<t_real>(tableAtoms->item(iRow, static_cast<int>(AtInfo::POS_X))->text().toStdString());
		t_real dY = tl::str_to_var<t_real>(tableAtoms->item(iRow, static_cast<int>(AtInfo::POS_Y))->text().toStdString());
		t_real dZ = tl::str_to_var<t_real>(tableAtoms->item(iRow, static_cast<int>(AtInfo::POS_Z))->text().toStdString());
		atom.vecPos = tl::make_vec({dX, dY, dZ});

		if(m_bEnableJ)
			atom.J = tl::str_to_var<t_real>(tableAtoms->item(iRow, static_cast<int>(AtInfo::J))->text().toStdString());

		vecAtoms.push_back(std::move(atom));
	}

	emit ApplyAtoms(vecAtoms);
}

void AtomsDlg::ButtonBoxClicked(QAbstractButton* pBtn)
{
	if(buttonBox->buttonRole(pBtn) == QDialogButtonBox::ApplyRole ||
		buttonBox->buttonRole(pBtn) == QDialogButtonBox::AcceptRole)
	{
		SendApplyAtoms();
	}
	else if(buttonBox->buttonRole(pBtn) == QDialogButtonBox::RejectRole)
	{
		reject();
	}

	if(buttonBox->buttonRole(pBtn) == QDialogButtonBox::AcceptRole)
	{
		if(m_pSettings)
			m_pSettings->setValue("atoms/geo", saveGeometry());

		QDialog::accept();
	}
}

void AtomsDlg::closeEvent(QCloseEvent* pEvt)
{
	QDialog::closeEvent(pEvt);
}


#include "AtomsDlg.moc"
