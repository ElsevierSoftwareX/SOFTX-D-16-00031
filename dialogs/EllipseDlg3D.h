/*
 * 3D Ellipsoid Dialog
 * @author Tobias Weber
 * @date may-2013, 29-apr-2014
 * @license GPLv2
 */

#ifndef __RESO_ELLI_DLG_3D__
#define __RESO_ELLI_DLG_3D__

#include <QDialog>
#include <QSettings>
#include <QComboBox>

#include <vector>

#include "libs/plotgl.h"
#include "tlibs/math/linalg.h"
#include "tools/res/ellipse.h"
#include "tools/res/defs.h"


class EllipseDlg3D : public QDialog
{Q_OBJECT
	protected:
		std::vector<PlotGl*> m_pPlots;
		std::vector<Ellipsoid3d<t_real_reso>> m_elliProj;
		std::vector<Ellipsoid3d<t_real_reso>> m_elliSlice;

		QComboBox *m_pComboCoord = nullptr;
		QSettings *m_pSettings = nullptr;

		ublas::matrix<t_real_reso> m_reso, m_resoHKL, m_resoOrient;
		ublas::vector<t_real_reso> m_Q_avg, m_Q_avgHKL, m_Q_avgOrient;
		ResoAlgo m_algo = ResoAlgo::UNKNOWN;

	protected:
		ublas::vector<t_real_reso>
		ProjRotatedVec(const ublas::matrix<t_real_reso>& rot,
			const ublas::vector<t_real_reso>& vec);

	public:
		EllipseDlg3D(QWidget* pParent, QSettings *pSett=0);
		virtual ~EllipseDlg3D();

	protected:
		void hideEvent(QHideEvent *event);
		void showEvent(QShowEvent *event);

	public slots:
		void SetParams(const ublas::matrix<t_real_reso>& reso, const ublas::vector<t_real_reso>& Q_avg,
			const ublas::matrix<t_real_reso>& resoHKL, const ublas::vector<t_real_reso>& Q_avgHKL,
			const ublas::matrix<t_real_reso>& resoOrient, const ublas::vector<t_real_reso>& Q_avgOrient,
			ResoAlgo algo);
		void Calc();
};

#endif
