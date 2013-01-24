///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010
//
//    This file is part of rootpwa
//
//    rootpwa is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    rootpwa is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with rootpwa. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
//
// Description:
//      sum accumulators with reduced round-off errors for large sums
//      of floating-point numbers
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef SUMACCUMULATORS_HPP
#define SUMACCUMULATORS_HPP


#include <vector>

#include <boost/accumulators/framework/accumulator_base.hpp>
#include <boost/accumulators/framework/parameters/sample.hpp>
#include <boost/accumulators/framework/depends_on.hpp>
#include <boost/accumulators/framework/extractor.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/sum.hpp>


namespace boost {
	namespace accumulators {


		///////////////////////////////////////////////////////////////////////////////
		// named parameter for cascaded sum
		BOOST_PARAMETER_NESTED_KEYWORD(tag, cascadedSumNmbElements, nmbElements)


		////////////////////////////////////////////////////////////////////////////
		// accumulator implementations
		namespace impl {

			// cascaded sum reduces (random) round-off errors form O[eps * sqrt(n)] (for naive accumulation)
			// to O[eps * sqrt(log n)] with very little overhead
			// see http://en.wikipedia.org/wiki/Pairwise_summation
			// the implementation is somewhat simplistic in that it requires buffer memory of
			// 3/4 of the size of the original data
			template<typename sampleT>
			struct cascadedSumAccumulator : accumulator_base	{

				typedef sampleT result_type;

				// constructor takes an Boost.Paramter argument pack there might be an initial value in
				// the argument pack; 'sample' is defined in sample.hpp
				template<typename argsT>
				cascadedSumAccumulator(const argsT& args)
					: _sumNextElement  (false),
					  _partialSumsIndex(0)
				{
					std::size_t nmbValues = args[cascadedSumNmbElements | 1];
					std::size_t nmbSums   = nmbValues / 2;
					if (nmbValues % 2 > 0)
						++nmbSums;
					if (nmbSums > 0) {
						_partialSums.resize(nmbSums, sampleT());
						_partialSums[0] = args[sample | sampleT()];
					}
				}

				template<typename argsT>
				void
				operator ()(const argsT& args)  // accumulate function also accepts an argument pack
				{
					// sum pairs and store result in _partialSums
					if (_sumNextElement) {
						_partialSums[_partialSumsIndex] += args[sample];
						++_partialSumsIndex;
					} else {
						if (_partialSumsIndex >= _partialSums.size())
							_partialSums.resize(_partialSumsIndex, sampleT());
						_partialSums[_partialSumsIndex] = args[sample];
					}
					_sumNextElement = !_sumNextElement;
				}

				result_type
				result(dont_care) const  // passed argument pack is not used
				{
					// sum up partial sums
					// treat special cases
					const std::size_t size = _partialSums.size();
					if (size == 0)
						return 0;
					else if (size == 1)
						return _partialSums[0];
					else if (size == 2)
						return _partialSums[0] + _partialSums[1];
					// allocate buffer for intermediate sums
					const std::size_t halfSize = size / 2;
					sampleT sum[halfSize + 1];
					// perform first summation pass
					std::size_t iNewSum, iOldSum;
					for (iNewSum = iOldSum = 0; iNewSum < halfSize; iNewSum++, iOldSum += 2)
						sum[iNewSum] = _partialSums[iOldSum] + _partialSums[iOldSum + 1];
					if (halfSize > 0)
						sum[halfSize] = _partialSums[size - 1];
					// sum pairs until only one pair remains
					std::size_t newSumSize = (size + 1) / 2;
					while (newSumSize > 2) {
						for (iNewSum = iOldSum = 0; iNewSum < newSumSize / 2; iNewSum++, iOldSum += 2)
							sum[iNewSum] = sum[iOldSum] + sum[iOldSum + 1];
						if (newSumSize % 2 > 0)
							sum[newSumSize / 2] = sum[newSumSize - 1];
						newSumSize = (newSumSize + 1) / 2;
					}
					return sum[0] + sum[1];
				}

			private:

				bool                 _sumNextElement;
				std::size_t          _partialSumsIndex;
				std::vector<sampleT> _partialSums;

			};


			// compensated (or Kahan) sum reduces (random) round-off errors form O[eps * sqrt(n)] (for
			// naive accumulation) to O[eps] at the cost of 4 summations (instead of 1) for each entry;
			// however, no additional buffer memory is required
			// see http://en.wikipedia.org/wiki/Kahan_summation
			template<typename sampleT>
			struct compensatedSumAccumulator : accumulator_base	{

				typedef sampleT result_type;

				// constructor takes an Boost.Paramter argument pack there might be an initial value in
				// the argument pack; 'sample' is defined in sample.hpp
				template<typename argsT>
				compensatedSumAccumulator(const argsT& args)
					: _sum(args[sample | sampleT()]),
					  _compensation(sampleT())
				{	}

				template<typename argsT>
				void
				operator ()(const argsT& args)  // accumulate function also accepts an argument pack
				{
					// at the beginning _compensation is zero
					const sampleT correctedTerm = args[sample] - _compensation;
					// _sum is big, correctedTerm small, so low-order digits of correctedTerm are lost
					const sampleT newSum        = _sum + correctedTerm;
					// (newSum - _sum) recovers the high-order part of correctedTerm;
					// subtracting correctedTerm recovers -(the lost low part of correctedTerm)
					_compensation = (newSum - _sum) - correctedTerm;
					_sum          = newSum;  // beware of eagerly optimising compilers
					// in the next call the lost low part stored in _compensation will be added to sum
				}

				result_type
				result(dont_care) const  // passed argument pack is not used
				{
					return _sum;
				}


			private:

				sampleT _sum;
				sampleT _compensation;  // running compensation for lost low-order bits

			};

		}  // namespace impl


		////////////////////////////////////////////////////////////////////////////
		// bind features to accumulator implementations
		namespace tag {

			struct cascadedSum : depends_on<>,
			                    tag::cascadedSumNmbElements
			{
				typedef accumulators::impl::cascadedSumAccumulator<mpl::_1> impl;
#ifdef BOOST_ACCUMULATORS_DOXYGEN_INVOKED
				static const parameter::keyword<tag::cascadedSumNmbElements> nmbElements;
#endif
			};


			struct compensatedSum : depends_on<>
			{
				typedef accumulators::impl::compensatedSumAccumulator<mpl::_1> impl;
			};

		}


		////////////////////////////////////////////////////////////////////////////
		// extractors for features
		namespace extract {

			const extractor<tag::cascadedSum> cascadedSum = {};
			BOOST_ACCUMULATORS_IGNORE_GLOBAL(cascadedSum)


			const extractor<tag::compensatedSum> compensatedSum = {};
			BOOST_ACCUMULATORS_IGNORE_GLOBAL(compensatedSum)

		}
		using extract::cascadedSum;
		using extract::compensatedSum;


		////////////////////////////////////////////////////////////////////////////
		// register feature variants

		// sum(cascaded) -> cascadedSum
		struct cascaded {};
		template<>
		struct as_feature<tag::sum(cascaded)>
		{
			typedef tag::cascadedSum type;
		};
    // for feature-based dependency resolution: provides the same feature as sum
		template<>
		struct feature_of<tag::cascadedSum> : feature_of<tag::sum>
    { };


		// sum(compensated) -> compensatedSum
		struct compensated {};
		template<>
		struct as_feature<tag::sum(compensated)>
		{
			typedef tag::compensatedSum type;
		};
    // for feature-based dependency resolution: provides the same feature as sum
		template<>
		struct feature_of<tag::compensatedSum> : feature_of<tag::sum>
    { };


	}  // namespace accumulators
}  // namespace boost


#endif  // SUMACCUMULATORS_HPP
