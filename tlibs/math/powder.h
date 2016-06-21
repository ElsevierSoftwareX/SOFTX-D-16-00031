/*
 * Powder peaks
 * @author Tobias Weber
 * @date apr-2015
 * @license GPLv2 or GPLv3
 */

#ifndef __POWDER_H__
#define __POWDER_H__

#include <tuple>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <initializer_list>
#include "lattice.h"
#include "../string/string.h"
#include "../helper/hash.h"
//#include <iostream>

namespace tl {

template<class t_int=int, class t_real=double>
class Powder
{
	public:
		// hkl G F
		typedef std::tuple<t_int, t_int, t_int, t_real, t_real> t_peak;
		typedef std::string t_str;

	private:
		static t_str to_str(t_real t)
		{
			static const int iPrec = 8;
			return tl::var_to_str<t_real, t_str>(t, iPrec);
		}

		static bool is_eq(t_real t0, t_real t1)
		{
			t_str str0 = to_str(t0);
			t_str str1 = to_str(t1);
			return str0 == str1;
		}

		static size_t hash_peak(const t_peak& peak)
		{
			return tl::hash_ordered<std::initializer_list<t_int>>
			({
				std::get<0>(peak),
				std::get<1>(peak),
				std::get<2>(peak)
			});
		}
		static size_t hash_peak_unique(const t_peak& peak)
		{
			t_str strG = to_str(std::get<3>(peak));
			t_str strF = to_str(std::get<4>(peak));

			return tl::hash_ordered<std::initializer_list<std::string>>({strG, strF});
		}

		static bool equ_peak(const t_peak& peak0, const t_peak& peak1)
		{
			return std::get<0>(peak0) == std::get<0>(peak1) &&
				std::get<1>(peak0) == std::get<1>(peak1) &&
				std::get<2>(peak0) == std::get<2>(peak1);
		}
		static bool equ_peak_unique(const t_peak& peak0, const t_peak& peak1)
		{
			bool bGEq = is_eq(std::get<3>(peak0), std::get<3>(peak1));
			bool bFEq = is_eq(std::get<4>(peak0), std::get<4>(peak1));

			return bGEq && bFEq;
		}

	public:
		typedef std::unordered_set<t_peak, decltype(&hash_peak), decltype(&equ_peak)> t_peaks;
		typedef std::unordered_set<t_peak, decltype(&hash_peak_unique), decltype(&equ_peak_unique)> t_peaks_unique;

	protected:
		t_peaks m_peaks;			// hashes & compares hkl
		t_peaks_unique m_peaks_unique;		// hashes & compares F & G

		// associated reciprocal lattice
		const Lattice<t_real> *m_pLatticeRecip = nullptr;

	public:
		Powder()
			: m_peaks(10, &hash_peak, &equ_peak),
			  m_peaks_unique(10, &hash_peak_unique, &equ_peak_unique)
		{}

		ublas::vector<t_real> GetRecipLatticePos(t_real dh, t_real dk, t_real dl) const
		{
			if(m_pLatticeRecip)
				return m_pLatticeRecip->GetPos(dh, dk, dl);

			return ublas::vector<t_real>();
		}

		t_real GetG(t_real dh, t_real dk, t_real dl) const
		{
			ublas::vector<t_real> vecG = GetRecipLatticePos(dh, dk, dl);
			if(vecG.size())
				return ublas::norm_2(vecG);
			return 0.;
		}

		void AddPeak(t_int h, t_int k, t_int l, t_real F=0.)
		{
			t_peak peak(h,k,l, GetG(h,k,l), F);

			m_peaks.insert(peak);
			m_peaks_unique.insert(peak);
		}

		void SetRecipLattice(const Lattice<t_real>* pLatt)
		{
			m_pLatticeRecip = pLatt;
		}

		const t_peaks& GetPeaks() const { return m_peaks; }
		const t_peaks& GetUniquePeaks() const { return m_peaks_unique; }

		bool HasPeak(int h, int k, int l) const
		{
			t_peak peak(h,k,l);
			return m_peaks.find(peak) != m_peaks.end();
		}

		bool HasUniquePeak(int h, int k, int l, t_real F) const
		{
			const t_real G = GetG(h,k,l);

			for(const t_peak& pk : m_peaks_unique)
			{
				if(is_eq(G, std::get<3>(pk)) && is_eq(F, std::get<4>(pk)))
					return 1;
			}
			return 0;
		}

		std::size_t GetMultiplicity(t_int h, t_int k, t_int l) const
		{
			t_real G = GetG(h,k,l);

			std::size_t iMult = 0;
			for(const t_peak& pk : m_peaks)
			{
				if(is_eq(G, std::get<3>(pk)))
					++iMult;
			}
			return iMult;
		}

		std::size_t GetMultiplicity(t_int h, t_int k, t_int l, t_real F) const
		{
			t_real G = GetG(h,k,l);

			std::size_t iMult = 0;
			for(const t_peak& pk : m_peaks)
			{
				if(is_eq(G, std::get<3>(pk)) && is_eq(F, std::get<4>(pk)))
					++iMult;
			}
			return iMult;
		}

		void clear()
		{
			m_peaks.clear();
			m_peaks_unique.clear();
			m_pLatticeRecip = nullptr;
		}
};

}

#include <ostream>

template<class t_int=int, class t_real=double>
std::ostream& operator<<(std::ostream& ostr, const tl::Powder<t_int,t_real>& powder)
{
	const typename tl::Powder<t_int,t_real>::t_peaks& peaks = powder.GetPeaks();
	const typename tl::Powder<t_int,t_real>::t_peaks_unique& peaks_unique = powder.GetUniquePeaks();

	ostr << "Peaks:\n";
	for(const typename tl::Powder<t_int,t_real>::t_peak& pk : peaks)
	{
		t_int h = std::get<0>(pk);
		t_int k = std::get<1>(pk);
		t_int l = std::get<2>(pk);

		ostr << "\t(" << h << k << l << ")\n";
	}

	ostr << "Unique Peaks:\n";
	for(const typename tl::Powder<t_int,t_real>::t_peak& pk : peaks_unique)
	{
		t_int h = std::get<0>(pk);
		t_int k = std::get<1>(pk);
		t_int l = std::get<2>(pk);

		ostr << "\t(" << h << k << l << ")";
		ostr << ", multiplicity: " << powder.GetMultiplicity(h,k,l) << "\n";
	}

	return ostr;
}

#endif
