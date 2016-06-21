/**
 * atoms and structural calculations
 * @author Tobias Weber
 * @date 2015-2016
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_ATOMS_H__
#define __TLIBS_ATOMS_H__


#include "linalg.h"
#include "linalg_ops.h"
#include "lattice.h"
#include <tuple>


namespace tl{

/**
 * Maps atom position back to units cell
 */
template<class t_vec>
void restrict_to_uc(t_vec& vec,
	typename t_vec::value_type tMin=0, typename t_vec::value_type tMax=1)
{
	using T = typename t_vec::value_type;

	for(std::size_t i=0; i<vec.size(); ++i)
	{
		vec[i] = std::fmod(vec[i], T(1));

		while(vec[i] < tMin) vec[i] += T(1);
		while(vec[i] >= tMax) vec[i] -= T(1);
	}
}


/**
 * Generates atom positions using trafo matrices
 */
template<class t_mat, class t_vec, template<class...> class t_cont>
t_cont<t_vec> generate_atoms(const t_cont<t_mat>& trafos, const t_vec& vecAtom,
	typename t_vec::value_type tUCMin=0, typename t_vec::value_type tUCMax=1,
	typename t_vec::value_type eps = std::numeric_limits<typename t_vec::value_type>::epsilon())
{
	//typedef typename t_vec::value_type t_real;
	t_cont<t_vec> vecvecRes;

	for(const t_mat& mat : trafos)
	{
		t_vec vecRes = mat * vecAtom;
		restrict_to_uc<t_vec>(vecRes, tUCMin, tUCMax);

		bool bPushBack = 1;
		// already have pos?
		for(const t_vec& vecOld : vecvecRes)
		{
			if(vec_equal(vecOld, vecRes, eps))
			{
				bPushBack = 0;
				break;
			}
		}

		if(bPushBack)
			vecvecRes.push_back(std::move(vecRes));
	}

	return vecvecRes;
}


/**
 * Generates atom positions using trafo matrices for all atoms in unit cell
 * @return tuple of (names, positions, positions in rlu, atom types)
 */
template<class t_mat, class t_vec, template<class...> class t_cont,
	class t_str=std::string, class t_real = typename t_mat::value_type>
std::tuple<t_cont<t_str>, t_cont<t_vec>, t_cont<t_vec>, t_cont<std::size_t>>
generate_all_atoms(const t_cont<t_mat>& trafos,
	const t_cont<t_vec>& vecAtoms, const t_cont<t_str>* pvecNames,
	const t_mat& matA, t_real tUCMin=0, t_real tUCMax=1,
	t_real eps = std::numeric_limits<t_real>::epsilon())
{
	t_cont<t_vec> vecAllAtoms, vecAllAtomsFrac;
	t_cont<t_str> vecAllNames;
	t_cont<std::size_t> vecAllAtomTypes;

	for(std::size_t iAtom=0; iAtom<vecAtoms.size(); ++iAtom)
	{
		t_vec vecAtom = vecAtoms[iAtom];
		t_str strNone;

		// homogeneous coordinates
		vecAtom.resize(4,1); vecAtom[3] = t_real(1);
		const t_str& strElem = pvecNames ? (*pvecNames)[iAtom] : strNone;

		t_cont<t_vec> vecOtherAtoms = vecAtoms;
		vecOtherAtoms.erase(vecOtherAtoms.begin() + iAtom);

		t_cont<t_vec> vecSymPos =
			tl::generate_atoms<t_mat, t_vec, t_cont>
				(trafos, vecAtom, tUCMin, tUCMax, eps);


		std::size_t iGeneratedAtoms = 0;
		for(t_vec vecThisAtom : vecSymPos)
		{
			vecThisAtom.resize(3,1);

			// is the atom position in the unit cell still free?
			if(std::find_if(vecAllAtomsFrac.begin(), vecAllAtomsFrac.end(),
				[&vecThisAtom, eps](const t_vec& _v) -> bool
				{ return tl::vec_equal(_v, vecThisAtom, eps); }) == vecAllAtomsFrac.end()
				&& // and is it not at a given initial atom position?
				std::find_if(vecOtherAtoms.begin(), vecOtherAtoms.end(),
				[&vecThisAtom, eps](const t_vec& _v) -> bool
				{ return tl::vec_equal(_v, vecThisAtom, eps); }) == vecOtherAtoms.end())
			{
				vecAllAtomsFrac.push_back(vecThisAtom);

				// converts from fractional coordinates
				vecThisAtom = matA * vecThisAtom;
				vecAllAtoms.push_back(std::move(vecThisAtom));
				vecAllNames.push_back(strElem);
				vecAllAtomTypes.push_back(iAtom);

				++iGeneratedAtoms;
			}
			else
			{
				tl::log_warn("Position ", vecThisAtom, " is already occupied,",
					" skipping current ", strElem, " atom.");
			}
		}
	}

	return std::make_tuple(vecAllNames,
		vecAllAtoms, vecAllAtomsFrac, vecAllAtomTypes);
}


/**
 * Generates supercell
 * @return tuple of positions, (user-defined) factors and indices
 */
template<class t_vec = ublas::vector<double>,
	template<class...> class t_cont = std::vector,
	class t_real = typename t_vec::value_type>
std::tuple<t_cont<t_vec>, t_cont<std::complex<t_real>>, t_cont<std::size_t>>
generate_supercell(const Lattice<t_real>& latt,
	const t_cont<t_vec>& vecAtomsUC,
	const t_cont<std::complex<t_real>>& vecFactsUC,
	std::ptrdiff_t N)
{
	using t_cplx = std::complex<t_real>;

	t_cont<t_vec> vecAllAtoms;
	t_cont<t_cplx> vecAllFacts;
	t_cont<std::size_t> vecAllIdx;

	for(std::ptrdiff_t h=-N+1; h<N; ++h)
	for(std::ptrdiff_t k=-N+1; k<N; ++k)
	for(std::ptrdiff_t l=-N+1; l<N; ++l)
	{
		t_vec vecPos = latt.GetPos(h,k,l);

		for(std::size_t iAtom=0; iAtom<vecAtomsUC.size(); ++iAtom)
		{
			vecAllIdx.push_back(iAtom);
			const t_vec& vecAtom = vecAtomsUC[iAtom];
			t_cplx cFact;
			if(vecFactsUC.size() == vecAtomsUC.size())
				cFact = vecFactsUC[iAtom];
			else if(vecFactsUC.size() == 1)		// use the same for all atoms
				cFact = vecFactsUC[0];

			vecAllAtoms.push_back(vecPos + vecAtom);
			if(vecFactsUC.size() != 0)
				vecAllFacts.push_back(cFact);
		}
	}

	return std::make_tuple(vecAllAtoms, vecAllFacts, vecAllIdx);
}
// ----------------------------------------------------------------------------


/**
 * calculates atomic form factors
 * @param G Length of lattice vector
 * @param vecA "a" coefficients
 * @param vecB "b" coefficients
 * @param c "c" coefficient
 * @return form factor
 * @desc: see Waasmaier and Kirfel, Acta Cryst. A51, 416-431 (1995)
 */
template<class T=double, template<class...> class t_cont>
T formfact(T G, const t_cont<T>& vecA, const t_cont<T>& vecB, T c)
{
	T ff = T(0);
	T s = G / T(4.*get_pi<T>());

	typename t_cont<T>::const_iterator iterA = vecA.begin();
	typename t_cont<T>::const_iterator iterB = vecB.begin();

	for(; iterA!=vecA.end() && iterB!=vecB.end(); ++iterA, ++iterB)
		ff += (*iterA)*std::exp(-(*iterB)*s*s);
	ff += c;

	return ff;
}



/**
 * calculates the structure factor F
 * @param lstAtoms List of atom positions
 * @param vecG Lattice vector
 * @param lstf G-dependent Atomic form factors (x-rays) or coherent scattering length (neutrons)
 * @param pF0 optional total form factor.
 * @return structure factor
 */
template<typename T = double, typename t_ff = std::complex<T>,
	class t_vec = ublas::vector<T>,
	template<class ...> class t_cont=std::initializer_list>
std::complex<T> structfact(const t_cont<t_vec>& lstAtoms, const t_vec& vecG,
	const t_cont<t_ff>& lstf = t_cont<t_ff>(),
	t_ff *pF0 = nullptr)
{
	constexpr std::complex<T> i(0., 1.);
	std::complex<T> F(0., 0.);

	using t_iter_atoms = typename t_cont<t_vec>::const_iterator;
	using t_iter_ffact = typename t_cont<t_ff>::const_iterator;

	t_iter_atoms iterAtom = lstAtoms.begin();
	t_iter_ffact iterFFact = lstf.begin();

	if(pF0) *pF0 = t_ff(0);

	for(; iterAtom != lstAtoms.end(); ++iterAtom)
	{
		// only use form factors or scattering lengths when available
		t_ff tFF = T(1);
		if(iterFFact != lstf.end())
			tFF = *iterFFact;

		F += tFF * std::exp(i * (vecG * *iterAtom));
		if(pF0) *pF0 += tFF;

		// if there is only one form factor in the list, use it for all positions
		if(iterFFact!=lstf.end() && std::next(iterFFact)!=lstf.end())
			++iterFFact;
	}
	return F;
}


/**
 * Lorentz factor
 * @param twotheta Scattering angle in rad
 */
template<typename T=double>
T lorentz_factor(T twotheta)
{
	T theta = 0.5*twotheta;
	return T(0.25) / (std::sin(theta)*std::sin(theta) * std::cos(theta));
}

/**
 * Lorentz polarisation factor (only for x-rays)
 * @param twotheta Scattering angle in rad
 */
template<typename T=double>
T lorentz_pol_factor(T twotheta)
{
	return T(0.5) + T(0.5)*std::cos(twotheta)*std::cos(twotheta);
}

}
#endif
