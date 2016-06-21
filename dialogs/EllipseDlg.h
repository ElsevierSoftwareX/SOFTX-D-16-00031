/*
 * Ellipse Dialog
 * @author Tobias Weber
 * @date 2013 - 2016
 * @license GPLv2
 */

#ifndef __RESO_ELLI_DLG__
#define __RESO_ELLI_DLG__

#include <QDialog>
#include <QSettings>
#if QT_VER>=5
	#include <QtWidgets>
#endif

#include <vector>
#include <memory>

#include "tlibs/math/linalg.h"
#include "tools/res/ellipse.h"
#include "tools/res/pop.h"
#include "ui/ui_ellipses.h"
#include "libs/qthelper.h"


class EllipseDlg : public QDialog, Ui::EllipseDlg
{ Q_OBJECT
	private:
		const char* m_pcTitle = "Resolution Ellipses";

	protected:
		bool m_bReady = 0;
		std::vector<std::unique_ptr<QwtPlotWrapper>> m_vecplotwrap;

		std::vector<struct Ellipse2d<t_real_reso>> m_elliProj;
		std::vector<struct Ellipse2d<t_real_reso>> m_elliSlice;

		std::vector<std::vector<t_real_reso> > m_vecXCurvePoints;
		std::vector<std::vector<t_real_reso> > m_vecYCurvePoints;

		QSettings *m_pSettings = 0;

	protected:
		ublas::matrix<t_real_reso> m_reso, m_resoHKL, m_resoOrient;
		ublas::vector<t_real_reso> m_Q_avg, m_Q_avgHKL, m_Q_avgOrient;
		ResoAlgo m_algo = ResoAlgo::UNKNOWN;

	public:
		EllipseDlg(QWidget* pParent=0, QSettings* pSett=0);
		virtual ~EllipseDlg();

	protected:
		virtual void showEvent(QShowEvent *pEvt) override;
		virtual void accept() override;

	protected slots:
		void cursorMoved(const QPointF& pt);

	public slots:
		void SetParams(const ublas::matrix<t_real_reso>& reso, const ublas::vector<t_real_reso>& Q_avg,
			const ublas::matrix<t_real_reso>& resoHKL, const ublas::vector<t_real_reso>& Q_avgHKL,
			const ublas::matrix<t_real_reso>& resoOrient, const ublas::vector<t_real_reso>& Q_avgOrient,
			ResoAlgo algo);
		void Calc();

	public:
		void SetTitle(const char* pcTitle);
};

#endif
