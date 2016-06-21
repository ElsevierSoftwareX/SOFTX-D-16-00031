/*
 * Ellipse Dialog
 * @author Tobias Weber
 * @date 2013 - 2016
 * @license GPLv2
 */

#include "EllipseDlg.h"
#include "tlibs/string/string.h"
#include "tlibs/string/spec_char.h"
#include "tlibs/helper/flags.h"

#include <future>


EllipseDlg::EllipseDlg(QWidget* pParent, QSettings* pSett)
	: QDialog(pParent), m_pSettings(pSett)
{
	setupUi(this);
	setWindowTitle(m_pcTitle);

	if(m_pSettings)
	{
		QFont font;
		if(m_pSettings->contains("main/font_gen") && font.fromString(m_pSettings->value("main/font_gen", "").toString()))
			setFont(font);
	}

	m_vecplotwrap.reserve(4);
	m_elliProj.resize(4);
	m_elliSlice.resize(4);
	m_vecXCurvePoints.resize(8);
	m_vecYCurvePoints.resize(8);

	QwtPlot* pPlots[] = {plot1, plot2, plot3, plot4};
	for(unsigned int i=0; i<4; ++i)
	{
		m_vecplotwrap.push_back(std::unique_ptr<QwtPlotWrapper>(new QwtPlotWrapper(pPlots[i], 2)));
		m_vecplotwrap[i]->GetPlot()->setMinimumSize(200,200);

		m_vecplotwrap[i]->GetCurve(0)->setTitle("Projected Ellipse");
		m_vecplotwrap[i]->GetCurve(1)->setTitle("Sliced Ellipse");

		QPen penProj, penSlice;
		penProj.setColor(QColor(0, 0x99,0));
		penSlice.setColor(QColor(0,0,0x99));
		penProj.setWidth(2);
		penSlice.setWidth(2);

		m_vecplotwrap[i]->GetCurve(0)->setPen(penProj);
		m_vecplotwrap[i]->GetCurve(1)->setPen(penSlice);

		if(m_vecplotwrap[i]->HasTrackerSignal())
		{
#if QT_VER >= 5
			connect(m_vecplotwrap[i]->GetPicker(), &QwtPlotPicker::moved,
				this, &EllipseDlg::cursorMoved);
#else
			connect(m_vecplotwrap[i]->GetPicker(), SIGNAL(moved(const QPointF&)),
				this, SLOT(cursorMoved(const QPointF&)));
#endif
		}
	}

#if QT_VER >= 5
	QObject::connect(comboCoord, static_cast<void(QComboBox::*)(int)>(&QComboBox::currentIndexChanged),
		this, &EllipseDlg::Calc);
#else
	QObject::connect(comboCoord, SIGNAL(currentIndexChanged(int)), this, SLOT(Calc()));
#endif

	if(m_pSettings && m_pSettings->contains("reso/ellipse_geo"))
		restoreGeometry(m_pSettings->value("reso/ellipse_geo").toByteArray());
	m_bReady = 1;
}

EllipseDlg::~EllipseDlg()
{
	m_vecplotwrap.clear();
}

void EllipseDlg::SetTitle(const char* pcTitle)
{
	QString strTitle = m_pcTitle;
	strTitle += " (";
	strTitle += pcTitle;
	strTitle += ")";
	this->setWindowTitle(strTitle);
}


void EllipseDlg::cursorMoved(const QPointF& pt)
{
	std::string strX = tl::var_to_str(pt.x(), g_iPrecGfx);
	std::string strY = tl::var_to_str(pt.y(), g_iPrecGfx);

	std::ostringstream ostr;
	ostr << "(" << strX << ", " << strY << ")";

	this->labelStatus->setText(ostr.str().c_str());
}


void EllipseDlg::Calc()
{
	if(!m_bReady) return;
	const EllipseCoordSys coord = static_cast<EllipseCoordSys>(comboCoord->currentIndex());

	const ublas::matrix<t_real_reso> *pReso = nullptr;
	const ublas::vector<t_real_reso> *pQavg = nullptr;

	switch(coord)
	{
		case EllipseCoordSys::Q_AVG:		// Q|| Qperp system in 1/A
			pReso = &m_reso; pQavg = &m_Q_avg;
			break;
		case EllipseCoordSys::RLU:		// rlu system
			pReso = &m_resoHKL; pQavg = &m_Q_avgHKL;
			break;
		case EllipseCoordSys::RLU_ORIENT:	// rlu system
			pReso = &m_resoOrient; pQavg = &m_Q_avgOrient;
			break;
		default:
			tl::log_err("Unknown coordinate system selected."); return;
	}


	const ublas::matrix<t_real_reso>& reso = *pReso;
	const ublas::vector<t_real_reso>& _Q_avg = *pQavg;

	try
	{
		int iParams[2][4][5] =
		{
			{	// projected
				{0, 3, 1, 2, -1},
				{1, 3, 0, 2, -1},
				{2, 3, 0, 1, -1},
				{0, 1, 3, 2, -1}
			},
			{	// sliced
				{0, 3, -1, 2, 1},
				{1, 3, -1, 2, 0},
				{2, 3, -1, 1, 0},
				{0, 1, -1, 2, 3}
			}
		};

		static const std::string strDeg = tl::get_spec_char_utf8("deg");


		bool bCenterOn0 = 1;
		if(m_pSettings)
			bCenterOn0 = m_pSettings->value("reso/center_around_origin", 1).toInt() != 0;

		ublas::vector<t_real_reso> Q_avg = _Q_avg;
		if(bCenterOn0)
			Q_avg = ublas::zero_vector<t_real_reso>(Q_avg.size());


		std::vector<std::future<Ellipse2d<t_real_reso>>> tasks_ell_proj, tasks_ell_slice;

		for(unsigned int iEll=0; iEll<4; ++iEll)
		{
			if(m_pSettings)
			{
				std::ostringstream ostrCfg;
				ostrCfg << "reso/ellipse_2d_" << iEll;
				std::string strProjPath = ostrCfg.str() + "_proj";
				std::string strSlicePath = ostrCfg.str() + "_slice";

				std::string strProj = m_pSettings->value(strProjPath.c_str(), "").toString().toStdString();
				std::string strSlice = m_pSettings->value(strSlicePath.c_str(), "").toString().toStdString();

				const std::string* strEllis[] = {&strProj, &strSlice};

				for(unsigned int iWhichEll=0; iWhichEll<2; ++iWhichEll)
				{
					if(*strEllis[iWhichEll] != "")
					{
						std::vector<int> vecIdx;
						tl::get_tokens<int>(*strEllis[iWhichEll], std::string(","), vecIdx);

						if(vecIdx.size() == 5)
						{
							for(unsigned int iParam=0; iParam<5; ++iParam)
								iParams[iWhichEll][iEll][iParam] = vecIdx[iParam];
						}
						else
						{
							tl::log_err("Error in res.conf: Wrong size of parameters for ",
								strProj,  ".");
						}
					}
				}
			}

			const int *iP = iParams[0][iEll];
			const int *iS = iParams[1][iEll];

			std::future<Ellipse2d<t_real_reso>> ell_proj =
				std::async(std::launch::deferred|std::launch::async,
				[=, &reso, &Q_avg]()
				{ return ::calc_res_ellipse<t_real_reso>(reso, Q_avg, iP[0], iP[1], iP[2], iP[3], iP[4]); });
			std::future<Ellipse2d<t_real_reso>> ell_slice =
				std::async(std::launch::deferred|std::launch::async,
				[=, &reso, &Q_avg]()
				{ return ::calc_res_ellipse<t_real_reso>(reso, Q_avg, iS[0], iS[1], iS[2], iS[3], iS[4]); });

			tasks_ell_proj.push_back(std::move(ell_proj));
			tasks_ell_slice.push_back(std::move(ell_slice));
		}

		for(unsigned int iEll=0; iEll<4; ++iEll)
		{
			m_elliProj[iEll] = tasks_ell_proj[iEll].get();
			m_elliSlice[iEll] = tasks_ell_slice[iEll].get();

			std::vector<t_real_reso>& vecXProj = m_vecXCurvePoints[iEll*2+0];
			std::vector<t_real_reso>& vecYProj = m_vecYCurvePoints[iEll*2+0];
			std::vector<t_real_reso>& vecXSlice = m_vecXCurvePoints[iEll*2+1];
			std::vector<t_real_reso>& vecYSlice = m_vecYCurvePoints[iEll*2+1];

			t_real_reso dBBProj[4], dBBSlice[4];
			m_elliProj[iEll].GetCurvePoints(vecXProj, vecYProj, 512, dBBProj);
			m_elliSlice[iEll].GetCurvePoints(vecXSlice, vecYSlice, 512, dBBSlice);

			set_qwt_data<t_real_reso>()(*m_vecplotwrap[iEll], vecXProj, vecYProj, 0, false);
			set_qwt_data<t_real_reso>()(*m_vecplotwrap[iEll], vecXSlice, vecYSlice, 1, false);


			std::ostringstream ostrSlope;
			ostrSlope.precision(4);
			ostrSlope << "Projected ellipse (green):\n";
			ostrSlope << "\tSlope: " << m_elliProj[iEll].slope << "\n";
			ostrSlope << "\tAngle: " << tl::r2d(m_elliProj[iEll].phi) << strDeg << "\n";
			ostrSlope << "\tArea " << m_elliProj[iEll].area << "\n";
			ostrSlope << "Sliced ellipse (blue):\n";
			ostrSlope << "\tSlope: " << m_elliSlice[iEll].slope << "\n";
			ostrSlope << "\tAngle: " << tl::r2d(m_elliSlice[iEll].phi) << strDeg << "\n";
			ostrSlope << "\tArea " << m_elliSlice[iEll].area;
			m_vecplotwrap[iEll]->GetPlot()->setToolTip(QString::fromUtf8(ostrSlope.str().c_str()));

			const std::string& strLabX = ellipse_labels(iParams[0][iEll][0], coord);
			const std::string& strLabY = ellipse_labels(iParams[0][iEll][1], coord);
			m_vecplotwrap[iEll]->GetPlot()->setAxisTitle(QwtPlot::xBottom, strLabX.c_str());
			m_vecplotwrap[iEll]->GetPlot()->setAxisTitle(QwtPlot::yLeft, strLabY.c_str());

			m_vecplotwrap[iEll]->GetPlot()->replot();

			QRectF rect;
			rect.setLeft(std::min(dBBProj[0], dBBSlice[0]));
			rect.setRight(std::max(dBBProj[1], dBBSlice[1]));
			rect.setTop(std::max(dBBProj[2], dBBSlice[2]));
			rect.setBottom(std::min(dBBProj[3], dBBSlice[3]));
			if(m_vecplotwrap[iEll]->GetZoomer())
				m_vecplotwrap[iEll]->GetZoomer()->setZoomBase(rect);

			switch(m_algo)
			{
				case ResoAlgo::CN: SetTitle("Cooper-Nathans Algorithm (TAS)"); break;
				case ResoAlgo::POP: SetTitle("Popovici Algorithm (TAS)"); break;
				case ResoAlgo::ECK: SetTitle("Eckold-Sobolev Algorithm (TAS)"); break;
				case ResoAlgo::VIOL: SetTitle("Violini Algorithm (TOF)"); break;
				default: SetTitle("Unknown Resolution Algorithm"); break;
			}
		}
	}
	catch(const std::exception& ex)
	{
		tl::log_err("Cannot calculate ellipses.");
		SetTitle("Error");
	}
}

void EllipseDlg::SetParams(const ublas::matrix<t_real_reso>& reso, const ublas::vector<t_real_reso>& Q_avg,
	const ublas::matrix<t_real_reso>& resoHKL, const ublas::vector<t_real_reso>& Q_avgHKL,
	const ublas::matrix<t_real_reso>& resoOrient, const ublas::vector<t_real_reso>& Q_avgOrient,
	ResoAlgo algo)
{
	m_reso = reso;
	m_resoHKL = resoHKL;
	m_resoOrient = resoOrient;
	m_Q_avg = Q_avg;
	m_Q_avgHKL = Q_avgHKL;
	m_Q_avgOrient = Q_avgOrient;
	m_algo = algo;

	Calc();
}

void EllipseDlg::accept()
{
	if(m_pSettings)
		m_pSettings->setValue("reso/ellipse_geo", saveGeometry());

	QDialog::accept();
}

void EllipseDlg::showEvent(QShowEvent *pEvt)
{
	QDialog::showEvent(pEvt);
}


#include "EllipseDlg.moc"
