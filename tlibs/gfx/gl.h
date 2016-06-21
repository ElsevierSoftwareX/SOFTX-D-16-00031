/*
 * gfx stuff
 * @author tweber
 * @date 22-dec-2014
 * @license GPLv2 or GPLv3
 */

#ifndef __GFX_STUFF_H__
#define __GFX_STUFF_H__

#include <string>
#include <ostream>
#include <unordered_map>
#include <GL/gl.h>

#include "../math/linalg.h"
#include "../math/geo.h"

namespace tl {

template<class T=double> 
using t_mat4_gen = ublas::matrix<T, ublas::row_major, ublas::bounded_array<T,4*4>>;
template<class T=double> 
using t_mat3_gen = ublas::matrix<T, ublas::row_major, ublas::bounded_array<T,3*3>>;
template<class T=double> 
using t_vec4_gen = ublas::vector<T, ublas::bounded_array<T,4>>;
template<class T=double> 
using t_vec3_gen = ublas::vector<T, ublas::bounded_array<T,3>>;

typedef t_mat4_gen<double> t_mat4;
typedef t_mat3_gen<double> t_mat3;
typedef t_vec4_gen<double> t_vec4;
typedef t_vec3_gen<double> t_vec3;


template<typename T=double, typename... Args>
void to_gl_array(const ublas::matrix<T, Args...>& mat, T* glmat)
{
	glmat[0]=mat(0,0);  glmat[1]=mat(1,0);  glmat[2]=mat(2,0);
	glmat[4]=mat(0,1);  glmat[5]=mat(1,1);  glmat[6]=mat(2,1);
	glmat[8]=mat(0,2);  glmat[9]=mat(1,2);  glmat[10]=mat(2,2);

	if(mat.size1()>=4 && mat.size2()>=4)
	{
		glmat[3]=mat(3,0); glmat[7]=mat(3,1); glmat[11]=mat(3,2);
		glmat[12]=mat(0,3); glmat[13]=mat(1,3); glmat[14]=mat(2,3); glmat[15]=mat(3,3);
	}
	else
	{
		glmat[3]=0; glmat[7]=0; glmat[11]=0;
		glmat[12]=0; glmat[13]=0; glmat[14]=0; glmat[15]=1;
	}
}

template<typename t_mat = t_mat4>
t_mat from_gl_array(const typename t_mat::value_type* glmat)
{
	t_mat mat(4,4);

	for(short j=0; j<4; ++j)
		for(short i=0; i<4; ++i)
			mat(i,j)=glmat[i + j*4];

	return mat;
}


template<typename t_mat = t_mat4, typename t_vec = t_vec4,
	typename T = typename t_mat::value_type>
void proj_pt(T dX, T dY, T dZ, const t_mat& matProj, const t_mat& matMV,
	T& dXProj, T& dYProj)
{
	t_mat mat = ublas::prod(matProj, matMV);
	t_vec vec = ublas::prod(mat, make_vec<t_vec>({dX, dY, dZ, 1.}));
	vec /= vec[3];

	//std::cout << vec << std::endl;
	dXProj = vec[0];
	dYProj = vec[1];
}

template<typename t_mat = t_mat4, typename t_vec = t_vec4,
	typename T = typename t_mat::value_type>
void gl_proj_pt(T dX, T dY, T dZ, T& dXProj, T& dYProj)
{
	GLdouble dMatMV[16], dMatProj[16];
	glGetDoublev(GL_PROJECTION_MATRIX, dMatProj);
	glGetDoublev(GL_MODELVIEW_MATRIX, dMatMV);

	t_mat matProj = from_gl_array<t_mat>(dMatProj);
	t_mat matMV = from_gl_array<t_mat>(dMatMV);

	proj_pt<t_mat, t_vec, T>(dX, dY, dZ, matProj, matMV, dXProj, dYProj);
}


template<typename t_mat = t_mat4, typename t_vec = t_vec4,
	typename T = typename t_mat::value_type>
void gl_mv_pt(const t_vec& vec, t_vec& vecOut)
{
	GLdouble dMatMV[16];
	glGetDoublev(GL_MODELVIEW_MATRIX, dMatMV);
	t_mat matMV = from_gl_array<t_mat>(dMatMV);
	t_mat matMV_inv;
	tl::inverse(matMV, matMV_inv);

	vecOut = ublas::prod(matMV_inv, vec);
}

// distance to the object defined by the current modelview matrix
template<typename t_mat = t_mat4, typename t_vec = t_vec4,
	typename T = typename t_mat::value_type>
T gl_dist_mv()
{
	t_vec vecPos;
	gl_mv_pt(make_vec<t_vec>({0.,0.,0.,1.}), vecPos);
	vecPos /= vecPos[3];
	vecPos[3] = 0.;
	T dDist = ublas::norm_2(vecPos);
	return dDist;
}

// size of projected sphere using fov angle
template<typename t_mat = t_mat4, typename t_vec = t_vec4,
	typename T = typename t_mat::value_type>
T gl_proj_sphere_size(T dFOV, T dRadius)
{
	return 2.*dRadius / (2.*gl_dist_mv<t_mat, t_vec, T>() * std::tan(0.5*dFOV));
}

// size of projected sphere using projection matrix
template<typename t_mat = t_mat4, typename t_vec = t_vec4,
	typename T = typename t_mat::value_type>
T gl_proj_sphere_size(T dRadius)
{
	GLdouble dMatProj[16];
	glGetDoublev(GL_PROJECTION_MATRIX, dMatProj);
	t_mat matProj = from_gl_array<t_mat>(dMatProj);

	T dDist = gl_dist_mv<t_mat, t_vec, T>();
	t_vec vec1 = make_vec<t_vec>({0., dRadius, dDist, 1.});
	t_vec vec2 = make_vec<t_vec>({0., -dRadius, dDist, 1.});

	t_vec vecProj1 = ublas::prod(matProj, vec1); vecProj1 /= vecProj1[3];
	t_vec vecProj2 = ublas::prod(matProj, vec2); vecProj2 /= vecProj2[3];
	t_vec vecProj = vecProj2 - vecProj1;

	vecProj[3] = vecProj[2] = 0.;
	return ublas::norm_2(vecProj);
}


// ray through screen coordinates (dX, dY)
// similar to: https://www.opengl.org/sdk/docs/man2/xhtml/gluUnProject.xml
template<typename t_mat = t_mat4, typename t_vec = t_vec4,
	typename T = typename t_mat::value_type>
Line<T> screen_ray(T dX, T dY, const t_mat& matProj, const t_mat& matMV)
{
	t_mat mat = ublas::prod(matProj, matMV);
	t_mat matInv;
	inverse(mat, matInv);

	t_vec vecNear = make_vec<t_vec>({dX, dY, -1., 1.});
	t_vec vecFar = make_vec<t_vec>({dX, dY, 1., 1.});

	vecNear = ublas::prod(matInv, vecNear);
	vecFar = ublas::prod(matInv, vecFar);

	vecNear /= vecNear[3];
	vecFar /= vecFar[3];

	ublas::vector<T> vecPos = vecNear; vecPos.resize(3, 1);
	ublas::vector<T> vecDir = vecFar-vecNear; vecDir.resize(3, 1);
	Line<T> line(vecPos, vecDir);
	return line;
}

template<typename t_mat = t_mat4, typename t_vec = t_vec4,
	typename T = typename t_mat::value_type>
Line<T> gl_screen_ray(T dX, T dY)
{
	GLdouble dMatMV[16], dMatProj[16];
	glGetDoublev(GL_PROJECTION_MATRIX, dMatProj);
	glGetDoublev(GL_MODELVIEW_MATRIX, dMatMV);

	t_mat matProj = from_gl_array<t_mat>(dMatProj);
	t_mat matMV = from_gl_array<t_mat>(dMatMV);

	return screen_ray<t_mat, t_vec, T>(dX, dY, matProj, matMV);
}

// --------------------------------------------------------------------------------


template<typename t_mat = ublas::matrix<double>, typename t_vec = ublas::vector<double>,
	typename T = typename t_mat::value_type>
class Cam
{
protected:
	t_vec m_vecDir = make_vec<t_vec>({0.,0.,1.});
	t_vec m_vecUp = make_vec<t_vec>({0.,1.,0.});
	t_vec m_vecRight = make_vec<t_vec>({1.,0.,0.});
	t_vec m_vecPos = make_vec<t_vec>({0.,0.,0.});

public:
	Cam() {}
	Cam(const t_vec& vecDir, const t_vec& vecUp)
		: m_vecDir(vecDir), m_vecUp(vecUp)
	{
		m_vecDir = ublas::norm_2(m_vecDir);
		m_vecUp = ublas::norm_2(m_vecUp);
		m_vecRight = cross_3(m_vecUp, m_vecDir);
		m_vecUp = cross_3(m_vecDir, m_vecRight);
	}

	virtual ~Cam() {}

	// only rotation matrix
	t_mat GetRotMatrix() const
	{
		return row_matrix({m_vecRight, m_vecUp, m_vecDir});
	}

	// homogeneous coordinates
	template<typename t_mat_h>
	t_mat_h GetHomMatrix() const
	{
		t_mat_h mat = GetRotMatrix();
		mat.resize(4,4,1);
		for(int i=0; i<3; ++i)
			mat(3,i) = mat(i,3) = 0.;
		mat(3,3) = 1.;

		t_mat_h matTrans = make_mat<t_mat_h>(
			{{ 1, 0, 0, m_vecPos[0]},
			{  0, 1, 0, m_vecPos[1]},
			{  0, 0, 1, m_vecPos[2]},
			{  0, 0, 0,          1 }});

		mat = ublas::prod(mat, matTrans);
		return mat;
	}

	void RotateAroundRight(T tAngle)
	{
		t_mat matRot = rotation_matrix<t_mat>(m_vecRight, tAngle);
		m_vecUp = ublas::prod(m_vecUp, matRot);
		m_vecDir = ublas::prod(m_vecDir, matRot);
	}
	void RotateAroundUp(T tAngle)
	{
		t_mat matRot = rotation_matrix<t_mat>(m_vecUp, tAngle);
		m_vecRight = ublas::prod(m_vecRight, matRot);
		m_vecDir = ublas::prod(m_vecDir, matRot);
	}
	void RotateAroundDir(T tAngle)
	{
		t_mat matRot = rotation_matrix<t_mat>(m_vecDir, tAngle);
		m_vecRight = ublas::prod(m_vecRight, matRot);
		m_vecUp = ublas::prod(m_vecUp, matRot);
	}

	void MoveForward(T t)
	{
		m_vecPos += m_vecDir * t;
	}

	const t_vec& GetPos() const { return m_vecPos; }
};

// --------------------------------------------------------------------------------
}


/*
 * Freetype rendering under OpenGL, inspired by:
 * https://en.wikibooks.org/wiki/OpenGL_Programming/Modern_OpenGL_Tutorial_Text_Rendering_01
 * https://en.wikibooks.org/wiki/OpenGL_Programming/Modern_OpenGL_Tutorial_Text_Rendering_02
 */

#define DEF_FONT "/usr/share/fonts/dejavu/DejaVuSansMono.ttf"
#define DEF_FONT_SIZE 12

#include <ft2build.h>
#include FT_FREETYPE_H

namespace tl {
class FontMap
{
	protected:
		bool m_bOk = 0;
		FT_Library m_ftLib = 0;
		FT_Face m_ftFace = 0;

		int m_iCharsPerLine=0, m_iLines=0;
		int m_iTileH=0, m_iTileW=0;
		int m_iPadH=DEF_FONT_SIZE/4, m_iPadW=DEF_FONT_SIZE/4;
		int m_iLargeW=0, m_iLargeH=0;
		unsigned char *m_pcLarge = nullptr;

		static const std::string m_strCharset;
		using t_offsmap = std::unordered_map<
			typename std::string::value_type,
			std::pair<int, int>>;
		t_offsmap m_mapOffs;

	protected:
		void UnloadFont();

	static void draw_tile(unsigned char* pcBuf,
		unsigned int iW, unsigned int iH,
		unsigned int iTileW, unsigned int iTileH,
		unsigned int iPosX, unsigned int iPosY,
		const unsigned char* pcGlyph,
		unsigned int iGlyphXOffs, unsigned int iGlyphYOffs,
		unsigned int iGlyphW, unsigned int iGlyphH);

	public:
		FontMap(const char* pcFont/*=DEF_FONT*/, int iSize/*=DEF_FONT_SIZE*/);
		FontMap();
		virtual ~FontMap();

		bool LoadFont(const char* pcFont, int iSize);
		bool LoadFont(FT_Face ftFace);

		void dump(std::ostream& ostr) const;
		void dump(std::ostream& ostr, const std::pair<int,int>& pair) const;
		virtual bool IsOk() const { return m_bOk; }

		const unsigned char* GetBuffer() const { return m_pcLarge; }
		std::pair<int, int> GetOffset(char ch) const;
};


class GlFontMap : public FontMap
{
	protected:
		bool m_bOk = 0;
		GLuint m_tex = 0;

	protected:
		bool CreateFontTexture();

	public:
		GlFontMap() = delete;
		GlFontMap(const char* pcFont=DEF_FONT, int iSize=DEF_FONT_SIZE);
		GlFontMap(FT_Face ftFace);
		virtual ~GlFontMap();

		virtual bool IsOk() const override { return m_bOk && FontMap::IsOk(); }

		void BindTexture();
		void DrawText(double dX, double dY, const std::string& str, bool bCenter=1);
		void DrawText(double dX, double dY, double dZ, const std::string& str, bool bCenter=1);
};

// --------------------------------------------------------------------------------


}

//
//#ifdef TLIBS_INC_HDR_IMPLS
//	#include "gl_impl.h"
//#endif

#endif
