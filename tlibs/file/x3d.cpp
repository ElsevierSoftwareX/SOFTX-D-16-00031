/**
 * simplified x3d file handling
 * @author tweber
 * @date nov-2015
 * @license GPLv2 or GPLv3
 */

#include "x3d.h"
#include <fstream>

namespace tl{


X3dElem::~X3dElem()
{
	for(X3dElem* pElem : m_vecChildren)
		delete pElem;
	m_vecChildren.clear();
}

void X3dElem::AddChild(X3dElem* pElem)
{
	if(pElem)
		m_vecChildren.push_back(pElem);
}

void X3dElem::Write(std::ostream& ostr) const
{
	for(const X3dElem *pElem : m_vecChildren)
		pElem->Write(ostr);
}

void X3dElem::WriteColor(std::ostream& ostr, const t_vec& vecColor)
{
	if(vecColor.size() >= 3)
	{
		ostr << "<Appearance>\n";
			ostr << "<Material diffuseColor=\""
				<< vecColor[0] << " " << vecColor[1] << " " << vecColor[2]
				<< "\" ";
		if(vecColor.size() >= 4)
			ostr << "transparency=\"" << vecColor[3] << "\" ";
		ostr << "/>\n";
		ostr << "</Appearance>\n";
	}

}

// -----------------------------------------------------------------------------

void X3dScene::Write(std::ostream& ostr) const
{
	ostr << "<Scene>\n";
	X3dElem::Write(ostr);
	ostr << "</Scene>\n";
}

// -----------------------------------------------------------------------------


void X3dTrafo::Write(std::ostream& ostr) const
{
	ostr << "<Transform ";

	if(m_vecTrans.size() >= 3)
	{
		ostr << "translation=\""
			<< m_vecTrans[0] << " " << m_vecTrans[1] << " " << m_vecTrans[2]
			<< "\" ";
	}
	if(m_vecScale.size() >= 3)
	{
		ostr << "scale=\""
			<< m_vecScale[0] << " " << m_vecScale[1] << " " << m_vecScale[2]
			<< "\" ";
	}
	if(m_bHasRot)
	{
		ostr << "rotation=\""
			<< m_quatRot.R_component_1() << " " << m_quatRot.R_component_2() << " "
			<< m_quatRot.R_component_3() << " " << m_quatRot.R_component_4()
			<< "\" ";
	}

	ostr << ">\n";

	X3dElem::Write(ostr);

	ostr << "</Transform>\n";
}

// -----------------------------------------------------------------------------

void X3dSphere::Write(std::ostream& ostr) const
{
	ostr << "<Shape>\n";
		ostr << "<Sphere radius=\"" << m_dRadius << "\" />\n";

	X3dElem::WriteColor(ostr, m_vecColor);
	ostr << "</Shape>\n";

	// TODO: child elements?
	//X3dElem::Write(ostr);
}

void X3dCube::Write(std::ostream& ostr) const
{
	ostr << "<Shape>\n";
		ostr << "<Box size=\""
			<< m_vecLength[0] << " "
			<< m_vecLength[1] << " "
			<< m_vecLength[2] << "\" />\n";

	X3dElem::WriteColor(ostr, m_vecColor);
	ostr << "</Shape>\n";

	// TODO: child elements?
	//X3dElem::Write(ostr);
}

void X3dCylinder::Write(std::ostream& ostr) const
{
	ostr << "<Shape>\n";
		ostr << "<Cylinder height=\"" << m_dHeight << "\""
			<< " radius=\"" << m_dRadius << "\" />\n";

	X3dElem::WriteColor(ostr, m_vecColor);
	ostr << "</Shape>\n";

	// TODO: child elements?
	//X3dElem::Write(ostr);
}

// -----------------------------------------------------------------------------


X3d::~X3d()
{}

void X3d::Write(std::ostream& ostr) const
{
	if(m_strComment.size())
		ostr << "<!--\n" << m_strComment << "\n-->\n";

	ostr << "<X3D>\n";
	m_scene.Write(ostr);
	ostr << "</X3D>\n";
}

bool X3d::Save(const char* pcFile) const
{
	std::ofstream ofstr(pcFile);
	if(!ofstr.is_open())
		return false;

	Write(ofstr);
	return true;
}


}
