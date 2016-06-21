/*
 * Scan viewer
 * @author tweber
 * @date mar-2015
 * @license GPLv2
 */

#include "scanviewer.h"
#include <QFileDialog>
#include <QTableWidget>
#include <QTableWidgetItem>
//#include <QProcess>
#include <iostream>
#include <set>
#include <string>
#include <algorithm>
#include <iterator>
#include <boost/filesystem.hpp>

#include "tlibs/string/string.h"
#include "tlibs/log/log.h"

using t_real = t_real_glob;
namespace fs = boost::filesystem;


#ifndef QWT_VER
	#define QWT_VER 6
#endif


ScanViewerDlg::ScanViewerDlg(QWidget* pParent)
	: QDialog(pParent, Qt::WindowTitleHint|Qt::WindowCloseButtonHint|Qt::WindowMinMaxButtonsHint),
		m_settings("tobis_stuff", "scanviewer"),
		m_vecExts({	".dat", ".DAT", ".scn", ".SCN", ".ng0", ".NG0", ".log", ".LOG" })
{
	this->setupUi(this);
	QFont font;
	if(m_settings.contains("main/font_gen") && font.fromString(m_settings.value("main/font_gen", "").toString()))
		setFont(font);

	splitter->setStretchFactor(0, 1);
	splitter->setStretchFactor(1, 2);


	// -------------------------------------------------------------------------
	// plot stuff
	QColor colorBck(240, 240, 240, 255);
	plot->setCanvasBackground(colorBck);

	m_plotwrap.reset(new QwtPlotWrapper(plot, 2, true));

	QPen penCurve;
	penCurve.setColor(QColor(0,0,0x99));
	penCurve.setWidth(2);
	m_plotwrap->GetCurve(0)->setPen(penCurve);
	m_plotwrap->GetCurve(0)->setStyle(QwtPlotCurve::CurveStyle::Lines);
	m_plotwrap->GetCurve(0)->setTitle("Scan Curve");

	QPen penPoints;
	penPoints.setColor(QColor(0xff,0,0));
	penPoints.setWidth(4);
	m_plotwrap->GetCurve(1)->setPen(penPoints);
	m_plotwrap->GetCurve(1)->setStyle(QwtPlotCurve::CurveStyle::Dots);
	m_plotwrap->GetCurve(1)->setTitle("Scan Points");
	// -------------------------------------------------------------------------


	// -------------------------------------------------------------------------
	// property map stuff
	tableProps->setColumnCount(2);
	tableProps->setColumnWidth(0, 150);
	tableProps->setColumnWidth(1, 350);
	//tableProps->sortByColumn(0);

	tableProps->setHorizontalHeaderItem(0, new QTableWidgetItem("Property"));
	tableProps->setHorizontalHeaderItem(1, new QTableWidgetItem("Value"));

	tableProps->verticalHeader()->setVisible(false);
	tableProps->verticalHeader()->setDefaultSectionSize(tableProps->verticalHeader()->minimumSectionSize()+4);
	// -------------------------------------------------------------------------

#if QT_VER>=5
	ScanViewerDlg *pThis = this;
	QObject::connect(editPath, &QLineEdit::textEdited, pThis, &ScanViewerDlg::ChangedPath);
	QObject::connect(listFiles, &QListWidget::currentItemChanged, pThis, &ScanViewerDlg::FileSelected);
	QObject::connect(btnBrowse, &QToolButton::clicked, pThis, &ScanViewerDlg::SelectDir);
	QObject::connect(comboX, static_cast<void (QComboBox::*)(const QString&)>(&QComboBox::currentIndexChanged), pThis, &ScanViewerDlg::XAxisSelected);
	QObject::connect(comboY, static_cast<void (QComboBox::*)(const QString&)>(&QComboBox::currentIndexChanged), pThis, &ScanViewerDlg::YAxisSelected);
	QObject::connect(tableProps, &QTableWidget::currentItemChanged, pThis, &ScanViewerDlg::PropSelected);
	QObject::connect(comboExport, static_cast<void (QComboBox::*)(int)>(&QComboBox::currentIndexChanged), pThis, &ScanViewerDlg::GenerateExternal);
#else
	QObject::connect(editPath, SIGNAL(textEdited(const QString&)),
		this, SLOT(ChangedPath()));
	//QObject::connect(listFiles, SIGNAL(itemSelectionChanged()),
	//	this, SLOT(FileSelected()));
	QObject::connect(listFiles, SIGNAL(currentItemChanged(QListWidgetItem*, QListWidgetItem*)),
		this, SLOT(FileSelected(QListWidgetItem*, QListWidgetItem*)));
	QObject::connect(btnBrowse, SIGNAL(clicked(bool)),
		this, SLOT(SelectDir()));
	QObject::connect(comboX, SIGNAL(currentIndexChanged(const QString&)),
		this, SLOT(XAxisSelected(const QString&)));
	QObject::connect(comboY, SIGNAL(currentIndexChanged(const QString&)),
		this, SLOT(YAxisSelected(const QString&)));
	QObject::connect(tableProps, SIGNAL(currentItemChanged(QTableWidgetItem*, QTableWidgetItem*)),
		this, SLOT(PropSelected(QTableWidgetItem*, QTableWidgetItem*)));
	QObject::connect(comboExport, SIGNAL(currentIndexChanged(int)),
		this, SLOT(GenerateExternal(int)));
#endif

	QString strDir = m_settings.value("last_dir", tl::wstr_to_str(fs::current_path().native()).c_str()).toString();
	editPath->setText(strDir);

	m_bDoUpdate = 1;
	ChangedPath();


	if(m_settings.contains("geo"))
		restoreGeometry(m_settings.value("geo").toByteArray());
}

ScanViewerDlg::~ScanViewerDlg()
{
	ClearPlot();
	tableProps->setRowCount(0);
}

void ScanViewerDlg::closeEvent(QCloseEvent* pEvt)
{
	m_settings.setValue("geo", saveGeometry());
	QDialog::closeEvent(pEvt);
}

void ScanViewerDlg::ClearPlot()
{
	if(m_pInstr)
	{
		delete m_pInstr;
		m_pInstr = nullptr;
	}

	m_vecX.clear();
	m_vecY.clear();

	set_qwt_data<t_real>()(*m_plotwrap, m_vecX, m_vecY, 0, 0);
	set_qwt_data<t_real>()(*m_plotwrap, m_vecX, m_vecY, 1, 0);

	m_strX = m_strY = m_strCmd = "";
	plot->setAxisTitle(QwtPlot::xBottom, "");
	plot->setAxisTitle(QwtPlot::yLeft, "");
	plot->setTitle("");

	auto edits = { editA, editB, editC,
		editAlpha, editBeta, editGamma,
		editPlaneX0, editPlaneX1, editPlaneX2,
		editPlaneY0, editPlaneY1, editPlaneY2,
		editTitle, editSample,
		editUser, editContact,
		editKfix, editTimestamp };
	for(auto* pEdit : edits)
		pEdit->setText("");

	comboX->clear();
	comboY->clear();
	textRoot->clear();

	m_plotwrap->GetPlot()->replot();
}

void ScanViewerDlg::SelectDir()
{
	QFileDialog::Option fileopt = QFileDialog::Option(0);
	if(!m_settings.value("main/native_dialogs", 1).toBool())
		fileopt = QFileDialog::DontUseNativeDialog;

	QString strCurDir = (m_strCurDir==""?".":m_strCurDir.c_str());
	QString strDir = QFileDialog::getExistingDirectory(this, "Select directory",
		strCurDir, QFileDialog::ShowDirsOnly | fileopt);
	if(strDir != "")
	{
		editPath->setText(strDir);
		ChangedPath();
	}
}

void ScanViewerDlg::XAxisSelected(const QString& strLab) { PlotScan(); }
void ScanViewerDlg::YAxisSelected(const QString& strLab) { PlotScan(); }

void ScanViewerDlg::FileSelected(QListWidgetItem *pItem, QListWidgetItem *pItemPrev)
{
	if(!pItem) return;

	m_strCurFile = pItem->text().toStdString();


	ClearPlot();
	std::string strFile = m_strCurDir + m_strCurFile;
	m_pInstr = tl::FileInstrBase<t_real>::LoadInstr(strFile.c_str());
	if(!m_pInstr) return;

	std::vector<std::string> vecScanVars = m_pInstr->GetScannedVars();
	std::string strCntVar = m_pInstr->GetCountVar();

	m_bDoUpdate = 0;
	int iIdxX=-1, iIdxY=-1, iCurIdx=0;
	const tl::FileInstrBase<t_real>::t_vecColNames& vecColNames = m_pInstr->GetColNames();
	for(const tl::FileInstrBase<t_real>::t_vecColNames::value_type& strCol : vecColNames)
	{
		comboX->addItem(strCol.c_str());
		comboY->addItem(strCol.c_str());

		if(vecScanVars.size() && vecScanVars[0]==strCol)
			iIdxX = iCurIdx;
		if(strCntVar==strCol)
			iIdxY = iCurIdx;

		++iCurIdx;
	}

	comboX->setCurrentIndex(iIdxX);
	comboY->setCurrentIndex(iIdxY);

	m_bDoUpdate = 1;

	ShowProps();
	PlotScan();
}

void ScanViewerDlg::PlotScan()
{
	if(m_pInstr==nullptr || !m_bDoUpdate)
		return;

	m_strX = comboX->currentText().toStdString();
	m_strY = comboY->currentText().toStdString();
	std::string strTitle = m_pInstr->GetTitle();
	m_strCmd = m_pInstr->GetScanCommand();

	m_vecX = m_pInstr->GetCol(m_strX.c_str());
	m_vecY = m_pInstr->GetCol(m_strY.c_str());

	std::array<t_real, 3> arrLatt = m_pInstr->GetSampleLattice();
	std::array<t_real, 3> arrAng = m_pInstr->GetSampleAngles();
	std::array<t_real, 3> arrPlaneX = m_pInstr->GetScatterPlane0();
	std::array<t_real, 3> arrPlaneY = m_pInstr->GetScatterPlane1();

	editA->setText(tl::var_to_str(arrLatt[0]).c_str());
	editB->setText(tl::var_to_str(arrLatt[1]).c_str());
	editC->setText(tl::var_to_str(arrLatt[2]).c_str());
	editAlpha->setText(tl::var_to_str(tl::r2d(arrAng[0])).c_str());
	editBeta->setText(tl::var_to_str(tl::r2d(arrAng[1])).c_str());
	editGamma->setText(tl::var_to_str(tl::r2d(arrAng[2])).c_str());

	editPlaneX0->setText(tl::var_to_str(arrPlaneX[0]).c_str());
	editPlaneX1->setText(tl::var_to_str(arrPlaneX[1]).c_str());
	editPlaneX2->setText(tl::var_to_str(arrPlaneX[2]).c_str());
	editPlaneY0->setText(tl::var_to_str(arrPlaneY[0]).c_str());
	editPlaneY1->setText(tl::var_to_str(arrPlaneY[1]).c_str());
	editPlaneY2->setText(tl::var_to_str(arrPlaneY[2]).c_str());

	labelKfix->setText(m_pInstr->IsKiFixed()
		? QString::fromWCharArray(L"ki (1/\x212b):")
		: QString::fromWCharArray(L"kf (1/\x212b):"));
	editKfix->setText(tl::var_to_str(m_pInstr->GetKFix()).c_str());

	editTitle->setText(strTitle.c_str());
	editSample->setText(m_pInstr->GetSampleName().c_str());
	editUser->setText(m_pInstr->GetUser().c_str());
	editContact->setText(m_pInstr->GetLocalContact().c_str());
	editTimestamp->setText(m_pInstr->GetTimestamp().c_str());


	plot->setAxisTitle(QwtPlot::xBottom, m_strX.c_str());
	plot->setAxisTitle(QwtPlot::yLeft, m_strY.c_str());
	plot->setTitle(m_strCmd.c_str());


	if(m_vecX.size()==0 || m_vecY.size()==0)
		return;

	set_qwt_data<t_real>()(*m_plotwrap, m_vecX, m_vecY, 0, 0);
	set_qwt_data<t_real>()(*m_plotwrap, m_vecX, m_vecY, 1, 1);

	GenerateExternal(comboExport->currentIndex());
}

void ScanViewerDlg::GenerateExternal(int iLang)
{
	textRoot->clear();
	if(!m_vecX.size() || !m_vecY.size())
		return;

	if(iLang == 0)
		GenerateForRoot();
	else if(iLang == 1)
		GenerateForGnuplot();
	else if(iLang == 2)
		GenerateForPython();
	else if(iLang == 3)
		GenerateForHermelin();
	else
		tl::log_err("Unknown external language.");
}

void ScanViewerDlg::GenerateForGnuplot()
{
	const std::string& strTitle = m_strCmd;
	const std::string& strLabelX = m_strX;
	const std::string& strLabelY = m_strY;

	std::string strPySrc =
R"RAWSTR(set term wxt

set xlabel "%%LABELX%%"
set ylabel "%%LABELY%%"
set title "%%TITLE%%"
set grid

set xrange [%%MINX%%:%%MAXX%%]
set yrange [%%MINY%%:%%MAXY%%]

plot "-" using 1:2:3 pt 7 with yerrorbars title "Data"
%%POINTS%%
end)RAWSTR";


	std::vector<t_real> vecYErr = m_vecY;
	std::for_each(vecYErr.begin(), vecYErr.end(), [](t_real& d) { d = std::sqrt(d); });

	auto minmaxX = std::minmax_element(m_vecX.begin(), m_vecX.end());
	auto minmaxY = std::minmax_element(m_vecY.begin(), m_vecY.end());
	t_real dMaxErrY = *std::max_element(vecYErr.begin(), vecYErr.end());

	std::ostringstream ostrPoints;

	for(std::size_t i=0; i<std::min(m_vecX.size(), m_vecY.size()); ++i)
	{
		ostrPoints << m_vecX[i]
			<< " " << m_vecY[i]
			<< " " << std::sqrt(m_vecY[i])
			<< "\n";
	}

	tl::find_and_replace<std::string>(strPySrc, "%%MINX%%", tl::var_to_str(*minmaxX.first));
	tl::find_and_replace<std::string>(strPySrc, "%%MAXX%%", tl::var_to_str(*minmaxX.second));
	tl::find_and_replace<std::string>(strPySrc, "%%MINY%%", tl::var_to_str(*minmaxY.first-dMaxErrY));
	tl::find_and_replace<std::string>(strPySrc, "%%MAXY%%", tl::var_to_str(*minmaxY.second+dMaxErrY));
	tl::find_and_replace<std::string>(strPySrc, "%%TITLE%%", strTitle);
	tl::find_and_replace<std::string>(strPySrc, "%%LABELX%%", strLabelX);
	tl::find_and_replace<std::string>(strPySrc, "%%LABELY%%", strLabelY);
	tl::find_and_replace<std::string>(strPySrc, "%%POINTS%%", ostrPoints.str());

	textRoot->setText(strPySrc.c_str());
}

void ScanViewerDlg::GenerateForPython()
{
	const std::string& strTitle = m_strCmd;
	const std::string& strLabelX = m_strX;
	const std::string& strLabelY = m_strY;

	std::string strPySrc =
R"RAWSTR(import numpy as np

x = np.array([ %%VECX%% ])
y = np.array([ %%VECY%% ])
yerr = np.array([ %%VECYERR%% ])

min = np.array([ %%MINX%%, %%MINY%% ])
max = np.array([ %%MAXX%%, %%MAXY%% ])
range = max-min
mid = min + range*0.5

yerr = [a if a!=0. else 0.001*range[1] for a in yerr]



import scipy.optimize as opt

def gauss_model(x, x0, sigma, amp, offs):
        return amp * np.exp(-0.5 * ((x-x0) / sigma)**2.) + offs

hints = [mid[0], range[0]*0.5, range[1]*0.5, min[1]]
popt, pcov = opt.curve_fit(gauss_model, x, y, sigma=yerr, absolute_sigma=True, p0=hints)

x_fine = np.linspace(min[0], max[0], 128)
y_fit = gauss_model(x_fine, *popt)



import matplotlib.pyplot as plt

plt.figure()

plt.xlim(min[0], max[0])
plt.ylim(min[1], max[1])

plt.title("%%TITLE%%")
plt.xlabel("%%LABELX%%")
plt.ylabel("%%LABELY%%")

plt.grid(True)
plt.errorbar(x,y,yerr, fmt="o")
plt.plot(x_fine, y_fit)
plt.show())RAWSTR";


	std::vector<t_real> vecYErr = m_vecY;
	std::for_each(vecYErr.begin(), vecYErr.end(), [](t_real& d) { d = std::sqrt(d); });

	auto minmaxX = std::minmax_element(m_vecX.begin(), m_vecX.end());
	auto minmaxY = std::minmax_element(m_vecY.begin(), m_vecY.end());
	t_real dMaxErrY = *std::max_element(vecYErr.begin(), vecYErr.end());

	std::ostringstream ostrX, ostrY, ostrYErr;

	for(std::size_t i=0; i<std::min(m_vecX.size(), m_vecY.size()); ++i)
	{
		ostrX << m_vecX[i] << ", ";
		ostrY << m_vecY[i] << ", ";
		ostrYErr << vecYErr[i] << ", ";
	}

	tl::find_and_replace<std::string>(strPySrc, "%%MINX%%", tl::var_to_str(*minmaxX.first));
	tl::find_and_replace<std::string>(strPySrc, "%%MAXX%%", tl::var_to_str(*minmaxX.second));
	tl::find_and_replace<std::string>(strPySrc, "%%MINY%%", tl::var_to_str(*minmaxY.first-dMaxErrY));
	tl::find_and_replace<std::string>(strPySrc, "%%MAXY%%", tl::var_to_str(*minmaxY.second+dMaxErrY));
	tl::find_and_replace<std::string>(strPySrc, "%%TITLE%%", strTitle);
	tl::find_and_replace<std::string>(strPySrc, "%%LABELX%%", strLabelX);
	tl::find_and_replace<std::string>(strPySrc, "%%LABELY%%", strLabelY);
	tl::find_and_replace<std::string>(strPySrc, "%%VECX%%", ostrX.str());
	tl::find_and_replace<std::string>(strPySrc, "%%VECY%%", ostrY.str());
	tl::find_and_replace<std::string>(strPySrc, "%%VECYERR%%", ostrYErr.str());

	textRoot->setText(strPySrc.c_str());
}

void ScanViewerDlg::GenerateForHermelin()
{
    std::string strStoatSrc =
R"RAWSTR(#!./hermelin -t

module_init()
{
	import("apps/instr.scr");

	global fit_dbg = 1;
	global theterm = "wxt";
	global norm_to_mon = 1;
}

scan_plot()
{
	scanfile = "%%FILE%%";
	[instr, datx, daty, datyerr, xlab, ylab] = load_instr(scanfile, norm_to_mon);
	title = "\"" + scanfile + "\"";

	maxx = max(datx); minx = min(datx);
	maxy = max(daty); miny = min(daty);
	rangex = maxx-minx; rangey = maxy-miny;
	midx = minx + rangex*0.5;


	gauss_pos = [midx];
	gauss_amp = [rangey*0.5];
	gauss_sig = [rangex*0.5];
	bckgrd = miny;
	fitsteps = ["xxx x"];

	outfile = "";
	#outfile = "scan_plot.pdf";

	thefit = fit_gauss_manual_singlestep(datx, daty, datyerr, gauss_pos, gauss_amp, gauss_sig, bckgrd, fitsteps);


	plotmap = map();
	plotmap = ["xlimits" : minx + " " + maxx];
	#plotmap += ["ylimits" : "0 0.5"];

	plot_gausses(1, thefit, [datx, daty, datyerr], title, xlab, ylab, outfile, plotmap);
}

main(args)
{
	scan_plot();
})RAWSTR";

	const std::string strFile = m_strCurDir + m_strCurFile;
	tl::find_and_replace<std::string>(strStoatSrc, "%%FILE%%", strFile);

	textRoot->setText(strStoatSrc.c_str());
}

void ScanViewerDlg::GenerateForRoot()
{
	const std::string& strTitle = m_strCmd;
	const std::string& strLabelX = m_strX;
	const std::string& strLabelY = m_strY;

	std::string strRootSrc =
R"RAWSTR(void scan_plot()
{
	const Double_t vecX[] = { %%VECX%% };
	const Double_t vecY[] = { %%VECY%% };
	const Double_t vecYErr[] = { %%VECYERR%% };

	const Double_t dMin[] = { %%MINX%%, %%MINY%% };
	const Double_t dMax[] = { %%MAXX%%, %%MAXY%% };
	const Int_t iSize = sizeof(vecX)/sizeof(*vecX);

	gStyle->SetOptFit(1);
	TCanvas *pCanvas = new TCanvas("canvas0", "Root Canvas", 800, 600);
	pCanvas->SetGrid(1,1);
	pCanvas->SetTicks(1,1);
	//pCanvas->SetLogy();

	TH1F *pFrame = pCanvas->DrawFrame(dMin[0], dMin[1], dMax[0], dMax[1], "");
	pFrame->SetTitle("%%TITLE%%");
	pFrame->SetXTitle("%%LABELX%%");
	pFrame->SetYTitle("%%LABELY%%");

	TGraphErrors *pGraph = new TGraphErrors(iSize, vecX, vecY, 0, vecYErr);
	pGraph->SetMarkerStyle(20);
	pGraph->Draw("P");
})RAWSTR";


	std::vector<t_real> vecYErr = m_vecY;
	std::for_each(vecYErr.begin(), vecYErr.end(), [](t_real& d) { d = std::sqrt(d); });

	auto minmaxX = std::minmax_element(m_vecX.begin(), m_vecX.end());
	auto minmaxY = std::minmax_element(m_vecY.begin(), m_vecY.end());
	t_real dMaxErrY = *std::max_element(vecYErr.begin(), vecYErr.end());

	std::ostringstream ostrX, ostrY, ostrYErr;

	for(std::size_t i=0; i<std::min(m_vecX.size(), m_vecY.size()); ++i)
	{
		ostrX << m_vecX[i] << ", ";
		ostrY << m_vecY[i] << ", ";
		ostrYErr << vecYErr[i] << ", ";
	}

	tl::find_and_replace<std::string>(strRootSrc, "%%MINX%%", tl::var_to_str(*minmaxX.first));
	tl::find_and_replace<std::string>(strRootSrc, "%%MAXX%%", tl::var_to_str(*minmaxX.second));
	tl::find_and_replace<std::string>(strRootSrc, "%%MINY%%", tl::var_to_str(*minmaxY.first-dMaxErrY));
	tl::find_and_replace<std::string>(strRootSrc, "%%MAXY%%", tl::var_to_str(*minmaxY.second+dMaxErrY));
	tl::find_and_replace<std::string>(strRootSrc, "%%TITLE%%", strTitle);
	tl::find_and_replace<std::string>(strRootSrc, "%%LABELX%%", strLabelX);
	tl::find_and_replace<std::string>(strRootSrc, "%%LABELY%%", strLabelY);
	tl::find_and_replace<std::string>(strRootSrc, "%%VECX%%", ostrX.str());
	tl::find_and_replace<std::string>(strRootSrc, "%%VECY%%", ostrY.str());
	tl::find_and_replace<std::string>(strRootSrc, "%%VECYERR%%", ostrYErr.str());

	textRoot->setText(strRootSrc.c_str());
}

void ScanViewerDlg::PropSelected(QTableWidgetItem *pItem, QTableWidgetItem *pItemPrev)
{
	if(!pItem)
		m_strSelectedKey = "";

	for(int iItem=0; iItem<tableProps->rowCount(); ++iItem)
	{
		const QTableWidgetItem *pKey = tableProps->item(iItem, 0);
		const QTableWidgetItem *pVal = tableProps->item(iItem, 1);

		if(pKey==pItem || pVal==pItem)
		{
			m_strSelectedKey = pKey->text().toStdString();
			break;
		}
	}
}

void ScanViewerDlg::ShowProps()
{
	if(m_pInstr==nullptr || !m_bDoUpdate)
		return;

	const tl::FileInstrBase<t_real>::t_mapParams& params = m_pInstr->GetAllParams();
	tableProps->setRowCount(params.size());

	const bool bSort = tableProps->isSortingEnabled();
	tableProps->setSortingEnabled(0);
	unsigned int iItem = 0;
	for(const tl::FileInstrBase<t_real>::t_mapParams::value_type& pair : params)
	{
		QTableWidgetItem *pItemKey = tableProps->item(iItem, 0);
		if(!pItemKey)
		{
			pItemKey = new QTableWidgetItem();
			tableProps->setItem(iItem, 0, pItemKey);
		}

		QTableWidgetItem* pItemVal = tableProps->item(iItem, 1);
		if(!pItemVal)
		{
			pItemVal = new QTableWidgetItem();
			tableProps->setItem(iItem, 1, pItemVal);
		}

		pItemKey->setText(pair.first.c_str());
		pItemVal->setText(pair.second.c_str());

		++iItem;
	}

	tableProps->setSortingEnabled(bSort);
	//tableProps->sortItems(0, Qt::AscendingOrder);


	// retain previous selection
	bool bHasSelection = 0;
	for(int iItem=0; iItem<tableProps->rowCount(); ++iItem)
	{
		const QTableWidgetItem *pItem = tableProps->item(iItem, 0);
		if(!pItem) continue;

		if(pItem->text().toStdString() == m_strSelectedKey)
		{
			tableProps->selectRow(iItem);
			bHasSelection = 1;
			break;
		}
	}

	if(!bHasSelection)
		tableProps->selectRow(0);
}

void ScanViewerDlg::ChangedPath()
{
	listFiles->clear();
	ClearPlot();
	tableProps->setRowCount(0);

	std::string strPath = editPath->text().toStdString();
	fs::path dir(strPath);
	if(fs::exists(dir) && fs::is_directory(dir))
	{
		m_strCurDir = tl::wstr_to_str(dir.native());
		tl::trim(m_strCurDir);
		if(*(m_strCurDir.begin()+m_strCurDir.length()-1) != fs::path::preferred_separator)
			m_strCurDir += fs::path::preferred_separator;
		UpdateFileList();

		m_settings.setValue("last_dir", QString(m_strCurDir.c_str()));
	}
}

void ScanViewerDlg::UpdateFileList()
{
	listFiles->clear();

	try
	{
		fs::path dir(m_strCurDir);
		fs::directory_iterator dir_begin(dir), dir_end;

		std::set<fs::path> lst;
		std::copy_if(dir_begin, dir_end, std::insert_iterator<decltype(lst)>(lst, lst.end()),
			[this](const fs::path& p) -> bool
			{
				std::string strExt = tl::wstr_to_str(p.extension().native());
				if(strExt == ".bz2" || strExt == ".gz" || strExt == ".z")
					strExt = "." + tl::wstr_to_str(tl::get_fileext2(p.filename().native()));

				if(this->m_vecExts.size() == 0)
					return true;
				return std::find(this->m_vecExts.begin(), this->m_vecExts.end(),
					strExt) != this->m_vecExts.end();
			});

		for(const fs::path& d : lst)
			listFiles->addItem(tl::wstr_to_str(d.filename().native()).c_str());
	}
	catch(const std::exception& ex)
	{}
}

#include "scanviewer.moc"
