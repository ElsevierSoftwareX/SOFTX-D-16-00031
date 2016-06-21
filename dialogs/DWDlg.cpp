/*
 * Scattering factors dialog (e.g. Debye-Waller factor)
 * @author tweber
 * @date 2013, jan-2015
 * @license GPLv2
 */

#include "DWDlg.h"

#include "tlibs/string/string.h"
#include "tlibs/math/neutrons.hpp"
#include "tlibs/math/atoms.h"
#include "tlibs/helper/array.h"

#include <boost/units/io.hpp>

using t_real = t_real_glob;
static const tl::t_length_si<t_real> angs = tl::get_one_angstrom<t_real>();
static const tl::t_energy_si<t_real> meV = tl::get_one_meV<t_real>();
static const tl::t_angle_si<t_real> rads = tl::get_one_radian<t_real>();
static const tl::t_temperature_si<t_real> kelvin = tl::get_one_kelvin<t_real>();


DWDlg::DWDlg(QWidget* pParent, QSettings *pSettings)
	: QDialog(pParent), m_pSettings(pSettings)
{
	this->setupUi(this);
	if(m_pSettings)
	{
		QFont font;
		if(m_pSettings->contains("main/font_gen") && font.fromString(m_pSettings->value("main/font_gen", "").toString()))
			setFont(font);
	}


	// -------------------------------------------------------------------------
	// Bose Factor stuff
	std::vector<QDoubleSpinBox*> vecSpinBoxesBose = {spinBoseT, spinBoseEMin, spinBoseEMax};
	for(QDoubleSpinBox* pSpin : vecSpinBoxesBose)
		QObject::connect(pSpin, SIGNAL(valueChanged(double)), this, SLOT(CalcBose()));

	m_plotwrapBose.reset(new QwtPlotWrapper(plotBose, 2));

	// positive Bose factor
	QPen penCurveBosePos;
	penCurveBosePos.setColor(QColor(0,0,0x99));
	penCurveBosePos.setWidth(2);
	m_plotwrapBose->GetCurve(0)->setTitle("Boson Creation");
	m_plotwrapBose->GetCurve(0)->setPen(penCurveBosePos);

	// negative Bose factor
	QPen penCurveBoseNeg;
	penCurveBoseNeg.setColor(QColor(0x99,0,0));
	penCurveBoseNeg.setWidth(2);
	m_plotwrapBose->GetCurve(1)->setTitle("Boson Annihilation");
	m_plotwrapBose->GetCurve(1)->setPen(penCurveBoseNeg);

	if(m_plotwrapBose->HasTrackerSignal())
		connect(m_plotwrapBose->GetPicker(), SIGNAL(moved(const QPointF&)), this, SLOT(cursorMoved(const QPointF&)));

	m_pLegendBose = new QwtLegend();
	m_plotwrapBose->GetPlot()->insertLegend(m_pLegendBose, QwtPlot::TopLegend);

	m_plotwrapBose->GetPlot()->setAxisTitle(QwtPlot::xBottom, "E (meV)");
	m_plotwrapBose->GetPlot()->setAxisTitle(QwtPlot::yLeft, "Bose Factor");

	CalcBose();


	// -------------------------------------------------------------------------
	// DW Factor stuff
	std::vector<QDoubleSpinBox*> vecSpinBoxes = {spinAMU_deb, spinTD_deb, spinT_deb, spinMinQ_deb, spinMaxQ_deb};
	for(QDoubleSpinBox* pSpin : vecSpinBoxes)
		QObject::connect(pSpin, SIGNAL(valueChanged(double)), this, SLOT(CalcDW()));

	m_plotwrapDW.reset(new QwtPlotWrapper(plot, 1));
	m_plotwrapDW->GetCurve(0)->setTitle("Debye-Waller Factor");

	if(m_plotwrapDW->HasTrackerSignal())
		connect(m_plotwrapDW->GetPicker(), SIGNAL(moved(const QPointF&)), this, SLOT(cursorMoved(const QPointF&)));

	m_plotwrapDW->GetPlot()->setAxisTitle(QwtPlot::xBottom, "Q (1/A)");
	m_plotwrapDW->GetPlot()->setAxisTitle(QwtPlot::yLeft, "DW Factor");

	CalcDW();


	// -------------------------------------------------------------------------
	// Ana Factor stuff
	std::vector<QDoubleSpinBox*> vecSpinBoxesAna = {spinAnad, spinMinkf, spinMaxkf};
	for(QDoubleSpinBox* pSpin : vecSpinBoxesAna)
		QObject::connect(pSpin, SIGNAL(valueChanged(double)), this, SLOT(CalcAna()));

	m_plotwrapAna.reset(new QwtPlotWrapper(plotAna, 1));
	m_plotwrapAna->GetCurve(0)->setTitle("Analyser Factor");

	if(m_plotwrapAna->HasTrackerSignal())
		connect(m_plotwrapAna->GetPicker(), SIGNAL(moved(const QPointF&)), this, SLOT(cursorMoved(const QPointF&)));

	m_plotwrapAna->GetPlot()->setAxisTitle(QwtPlot::xBottom, "kf (1/A)");
	m_plotwrapAna->GetPlot()->setAxisTitle(QwtPlot::yLeft, "Intensity (a.u.)");

	CalcAna();


	// -------------------------------------------------------------------------
	// Lorentz Factor stuff
	std::vector<QDoubleSpinBox*> vecSpinBoxesLor = {spinMin2Th, spinMax2Th};
	for(QDoubleSpinBox* pSpin : vecSpinBoxesLor)
		QObject::connect(pSpin, SIGNAL(valueChanged(double)), this, SLOT(CalcLorentz()));
	QObject::connect(checkPol, SIGNAL(toggled(bool)), this, SLOT(CalcLorentz()));

	m_plotwrapLor.reset(new QwtPlotWrapper(plotLorentz, 1));
	m_plotwrapLor->GetCurve(0)->setTitle("Lorentz Factor");

	if(m_plotwrapLor->HasTrackerSignal())
		connect(m_plotwrapLor->GetPicker(), SIGNAL(moved(const QPointF&)), this, SLOT(cursorMoved(const QPointF&)));

	m_plotwrapLor->GetPlot()->setAxisTitle(QwtPlot::xBottom, "Scattering Angle (deg)");
	m_plotwrapLor->GetPlot()->setAxisTitle(QwtPlot::yLeft, "Lorentz Factor");

	CalcLorentz();


	if(m_pSettings && m_pSettings->contains("dw/geo"))
		restoreGeometry(m_pSettings->value("dw/geo").toByteArray());
}


DWDlg::~DWDlg()
{}

void DWDlg::cursorMoved(const QPointF& pt)
{
	std::string strX = tl::var_to_str(pt.x(), g_iPrecGfx);
	std::string strY = tl::var_to_str(pt.y(), g_iPrecGfx);

	std::ostringstream ostr;
	ostr << "(" << strX << ", " << strY << ")";

	this->labelStatus->setText(ostr.str().c_str());
}


void DWDlg::CalcBose()
{
	const unsigned int NUM_POINTS = 512;

	const t_real dMinE = spinBoseEMin->value();
	const t_real dMaxE = spinBoseEMax->value();

	const tl::t_temperature_si<t_real> T = t_real(spinBoseT->value()) * kelvin;

	m_vecBoseE.clear();
	m_vecBoseIntPos.clear();
	m_vecBoseIntNeg.clear();

	m_vecBoseE.reserve(NUM_POINTS);
	m_vecBoseIntPos.reserve(NUM_POINTS);
	m_vecBoseIntNeg.reserve(NUM_POINTS);

	for(unsigned int iPt=0; iPt<NUM_POINTS; ++iPt)
	{
		tl::t_energy_si<t_real> E = (dMinE + (dMaxE - dMinE)/t_real(NUM_POINTS)*t_real(iPt)) * meV;
		m_vecBoseE.push_back(E / meV);

		m_vecBoseIntPos.push_back(tl::bose(E, T));
		m_vecBoseIntNeg.push_back(tl::bose(-E, T));
	}

	set_qwt_data<t_real>()(*m_plotwrapBose, m_vecBoseE, m_vecBoseIntPos, 0, false);
	set_qwt_data<t_real>()(*m_plotwrapBose, m_vecBoseE, m_vecBoseIntNeg, 1, false);

	std::vector<t_real> vecMerged;
	vecMerged.resize(m_vecBoseIntPos.size() + m_vecBoseIntNeg.size());
	std::merge(m_vecBoseIntPos.begin(), m_vecBoseIntPos.end(),
		m_vecBoseIntNeg.begin(), m_vecBoseIntNeg.end(), vecMerged.begin());
	set_zoomer_base(m_plotwrapBose->GetZoomer(),
		tl::container_cast<t_real_qwt, t_real, std::vector>()(m_vecBoseE),
		tl::container_cast<t_real_qwt, t_real, std::vector>()(vecMerged));

	m_plotwrapBose->GetPlot()->replot();
}


void DWDlg::CalcLorentz()
{
	const unsigned int NUM_POINTS = 512;

	const t_real dMin2th = tl::d2r(spinMin2Th->value());
	const t_real dMax2th = tl::d2r(spinMax2Th->value());

	const bool bPol = checkPol->isChecked();

	m_vecLor2th.clear();
	m_vecLor.clear();

	m_vecLor2th.reserve(NUM_POINTS);
	m_vecLor.reserve(NUM_POINTS);

	for(unsigned int iPt=0; iPt<NUM_POINTS; ++iPt)
	{
		t_real d2th = (dMin2th + (dMax2th - dMin2th)/t_real(NUM_POINTS)*t_real(iPt));
		t_real dLor = tl::lorentz_factor(d2th);
		if(bPol)
			dLor *= tl::lorentz_pol_factor(d2th);

		m_vecLor2th.push_back(tl::r2d(d2th));
		m_vecLor.push_back(dLor);
	}

	set_qwt_data<t_real>()(*m_plotwrapLor, m_vecLor2th, m_vecLor, 0, true);
}


void DWDlg::CalcDW()
{
	const unsigned int NUM_POINTS = 512;

	t_real dMinQ = spinMinQ_deb->value();
	t_real dMaxQ = spinMaxQ_deb->value();

	tl::t_temperature_si<t_real> T = t_real(spinT_deb->value()) * kelvin;
	tl::t_temperature_si<t_real> T_D = t_real(spinTD_deb->value()) * kelvin;
	tl::t_mass_si<t_real> M = t_real(spinAMU_deb->value()) * tl::get_amu<t_real>();

	m_vecQ.clear();
	m_vecDeb.clear();

	m_vecQ.reserve(NUM_POINTS);
	m_vecDeb.reserve(NUM_POINTS);

	bool bHasZetaSq = 0;
	for(unsigned int iPt=0; iPt<NUM_POINTS; ++iPt)
	{
		tl::t_wavenumber_si<t_real> Q = (dMinQ + (dMaxQ - dMinQ)/t_real(NUM_POINTS)*t_real(iPt)) / angs;
		t_real dDWF = 0.;
		auto zetasq = angs*angs;

		if(T <= T_D)
			dDWF = tl::debye_waller_low_T(T_D, T, M, Q, &zetasq);
		else
			dDWF = tl::debye_waller_high_T(T_D, T, M, Q, &zetasq);

		m_vecQ.push_back(Q * angs);

		if(!bHasZetaSq)
		{
			std::string strZetaSq = tl::var_to_str(t_real(tl::my_units_sqrt<tl::t_length_si<t_real>>(zetasq) / angs));
			editZetaSq->setText(strZetaSq.c_str());

			bHasZetaSq = 1;
		}

		m_vecDeb.push_back(dDWF);
	}

	set_qwt_data<t_real>()(*m_plotwrapDW, m_vecQ, m_vecDeb, 0, true);
}


void DWDlg::CalcAna()
{
	const unsigned int NUM_POINTS = 512;

	const tl::t_length_si<t_real> d = t_real(spinAnad->value()) * angs;
	const t_real dMinKf = spinMinkf->value();
	const t_real dMaxKf = spinMaxkf->value();

	t_real dAngMax = 0.5*tl::r2d(tl::get_mono_twotheta(dMinKf/angs, d, 1) / rads);
	t_real dAngMin = 0.5*tl::r2d(tl::get_mono_twotheta(dMaxKf/angs, d, 1) / rads);

	editAngMin->setText(tl::var_to_str(dAngMin).c_str());
	editAngMax->setText(tl::var_to_str(dAngMax).c_str());

	m_veckf.clear();
	m_vecInt.clear();

	m_veckf.reserve(NUM_POINTS);
	m_vecInt.reserve(NUM_POINTS);

	for(unsigned int iPt=0; iPt<NUM_POINTS; ++iPt)
	{
		tl::t_wavenumber_si<t_real> kf = (dMinKf + (dMaxKf - dMinKf)/t_real(NUM_POINTS)*t_real(iPt)) / angs;
		t_real dEffic = tl::ana_effic_factor(kf, d);

		m_veckf.push_back(kf * angs);
		m_vecInt.push_back(dEffic);
	}

	set_qwt_data<t_real>()(*m_plotwrapAna, m_veckf, m_vecInt, 0, true);
}


void DWDlg::showEvent(QShowEvent *pEvt)
{
	QDialog::showEvent(pEvt);
}

void DWDlg::accept()
{
	if(m_pSettings)
		m_pSettings->setValue("dw/geo", saveGeometry());

	QDialog::accept();
}


#include "DWDlg.moc"

