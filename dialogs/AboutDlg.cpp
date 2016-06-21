/*
 * About Dialog
 * @author Tobias Weber
 * @date nov-2015
 * @license GPLv2
 */

#include "AboutDlg.h"

#include <boost/config.hpp>
#include <boost/version.hpp>
#include <qwt_global.h>
#include "tlibs/version.h"

#include "tlibs/string/string.h"
#include "libs/formfactors/formfact.h"
#include "libs/spacegroups/spacegroup.h"
#include "libs/globals.h"
#include "libs/version.h"
#include <sstream>


AboutDlg::AboutDlg(QWidget* pParent, QSettings *pSett)
	: QDialog(pParent), m_pSettings(pSett)
{
	setupUi(this);
	if(m_pSettings)
	{
		QFont font;
		if(m_pSettings->contains("main/font_gen") && font.fromString(m_pSettings->value("main/font_gen", "").toString()))
			setFont(font);

		if(m_pSettings->contains("about/geo"))
			restoreGeometry(m_pSettings->value("about/geo").toByteArray());
	}

	labelVersion->setText("Version " TAKIN_VER);
	labelWritten->setText("Written by Tobias Weber <tobias.weber@tum.de>");
	labelYears->setText("2014 - 2016");

	std::string strCC = "Built using " + std::string(BOOST_COMPILER);
#ifdef __cplusplus
	strCC += " (standard: " + tl::var_to_str(__cplusplus) + ")";
#endif
#ifdef BOOST_STDLIB
		strCC += " with " + std::string(BOOST_STDLIB);
#endif
#ifdef BOOST_PLATFORM
		strCC += " for " + std::string(BOOST_PLATFORM);
#endif
	strCC += ".";
	labelCC->setText(strCC.c_str());
	labelBuildDate->setText(QString("Build date: ") +
		QString(__DATE__) + ", " + QString(__TIME__));


	// -------------------------------------------------------------------------


	std::ostringstream ostrLibs;
	ostrLibs << "<html><body>";

	ostrLibs << "<dl>";

	ostrLibs << "<dt>Uses Qt version " << QT_VERSION_STR << "</dt>";
	ostrLibs << "<dd><a href=\"http://qt-project.org\">http://qt-project.org</a></dd>";

	ostrLibs << "<dt>Uses Qwt version " << QWT_VERSION_STR << "</dt>";
	ostrLibs << "<dd><a href=\"http://qwt.sourceforge.net\">http://qwt.sourceforge.net</a></dd>";

	std::string strBoost = BOOST_LIB_VERSION;
	tl::find_all_and_replace<std::string>(strBoost, "_", ".");
	ostrLibs << "<dt>Uses Boost version " << strBoost << "</dt>";
	ostrLibs << "<dd><a href=\"http://www.boost.org\">http://www.boost.org</a></dd>";

#ifndef NO_LAPACK
	ostrLibs << "<dt>Uses Lapack/e version 3</dt>";
	ostrLibs << "<dd><a href=\"http://www.netlib.org/lapack\">http://www.netlib.org/lapack</a></dd>";
#endif

	ostrLibs << "<dt>Uses tLibs version " << TLIBS_VERSION << "</dt>";

	ostrLibs << "<dt>Uses Clipper crystallography library</dt>";
	ostrLibs << "<dd><a href=\"http://www.ysbl.york.ac.uk/~cowtan/clipper\">http://www.ysbl.york.ac.uk/~cowtan/clipper</a></dd>";

	ostrLibs << "<dt>Uses resolution algorithms ported from Rescal version 5</dt>";
	ostrLibs << "<dd><a href=\"http://www.ill.eu/en/html/instruments-support/computing-for-science/cs-software/all-software/matlab-ill/rescal-for-matlab\">http://www.ill.eu/en/html/instruments-support/computing-for-science/cs-software/all-software/matlab-ill/rescal-for-matlab</a></dd>";

	ostrLibs << "<dt>Uses Tango icons</dt>";
	ostrLibs << "<dd><a href=\"http://tango.freedesktop.org\">http://tango.freedesktop.org</a></dd>";

	ostrLibs << "</dl>";
	ostrLibs << "</body></html>";

	labelLibraries->setText(ostrLibs.str().c_str());


	// -------------------------------------------------------------------------


	std::ostringstream ostrConst;
	ostrConst << "<html><body>";

	ostrConst << "<dl>";

	ostrConst << "<dt>Physical constants from Boost Units</dt>";
	ostrConst << "<dd><a href=\"http://www.boost.org/doc/libs/release/libs/units/\">http://www.boost.org/doc/libs/release/libs/units/</a></dd>";

	std::shared_ptr<const SpaceGroups<t_real_glob>> sgs = SpaceGroups<t_real_glob>::GetInstance();
	ostrConst << "<dt>" << sgs->get_sgsource(0) <<"</dt>";
	ostrConst << "<dd><a href=\"" << sgs->get_sgsource(1) << "\">" << sgs->get_sgsource(1) << "</a></dd>";

	std::shared_ptr<const FormfactList<t_real_glob>> ff = FormfactList<t_real_glob>::GetInstance();
	std::shared_ptr<const MagFormfactList<t_real_glob>> mff = MagFormfactList<t_real_glob>::GetInstance();
	std::shared_ptr<const ScatlenList<t_real_glob>> sl = ScatlenList<t_real_glob>::GetInstance();

	if(g_bHasFormfacts)
	{
		ostrConst << "<dt>" << ff->GetSource() << "</dt>";
		ostrConst << "<dd><a href=\"" << ff->GetSourceUrl() << "\">" << ff->GetSourceUrl() << "</a></dd>";
	}
	if(g_bHasMagFormfacts)
	{
		ostrConst << "<dt>" << mff->GetSource() << "</dt>";
		ostrConst << "<dd><a href=\"" << mff->GetSourceUrl() << "\">" << mff->GetSourceUrl() << "</a></dd>";
	}
	if(g_bHasScatlens)
	{
		ostrConst << "<dt>" << sl->GetSource() << "</dt>";
		ostrConst << "<dd><a href=\"" << sl->GetSourceUrl() << "\">" << sl->GetSourceUrl() << "</a></dd>";
	}

	ostrConst << "</dl>";
	ostrConst << "</body></html>";
	labelConst->setText(ostrConst.str().c_str());



	std::string strLicensesFile = find_resource("LICENSES");
	std::ifstream ifstrLicenses(strLicensesFile);
	std::string strLicenses;
	while(ifstrLicenses)
	{
		std::string strLic;
		std::getline(ifstrLicenses, strLic);
		strLicenses += strLic + "\n";
	}
	editAllLicenses->setPlainText(strLicenses.c_str());
}


void AboutDlg::accept()
{
	if(m_pSettings)
		m_pSettings->setValue("about/geo", saveGeometry());

	QDialog::accept();
}


#include "AboutDlg.moc"
