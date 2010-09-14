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
// File and Version Information:
// $Rev::                             $: revision of last commit
// $Author::                          $: author of last commit
// $Date::                            $: date of last commit
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

#include "boost/accumulators/framework/accumulator_base.hpp"
#include "boost/accumulators/framework/parameters/sample.hpp"
#include "boost/accumulators/framework/depends_on.hpp"
#include "boost/accumulators/framework/extractor.hpp"


namespace boost {
	namespace accumulators {
		
		
		///////////////////////////////////////////////////////////////////////////////
		// named parameter for cascaded sum
		BOOST_PARAMETER_NESTED_KEYWORD(tag, cascadeSumNmbElements,          nmbElements)
		BOOST_PARAMETER_NESTED_KEYWORD(tag, cascadeSumNmbOfSumsAtEachStage, nmbOfSumsAtEachStage)


		////////////////////////////////////////////////////////////////////////////
		// accumulator implementations
		namespace impl {

			// cascaded sum reduces (random) round-off errors form O[eps * sqrt(n)] (for naive accumulation)
			// to O[eps * sqrt(log n)] (for nmbOfSumsAtEachStage = 2) with very little overhead
			// see http://en.wikipedia.org/wiki/Pairwise_summation
			// the implementation is somewhat simplistic in that it requires buffer memory of
			// 1 / nmbOfSumsAtEachStage of the size of the original data
			template<typename sampleT>
			struct cascadeSumAccumulator : accumulator_base	{

				typedef sampleT result_type;
    
				// constructor takes an Boost.Paramter argument pack there might be an initial value in
				// the argument pack; 'sample' is defined in sample.hpp
				template<typename argsT>
				cascadeSumAccumulator(const argsT& args)
					: _countElements       (0),
					  _nmbOfSumsAtEachStage(args[cascadeSumNmbOfSumsAtEachStage | 2]),
					  _partialSumsIndex    (0)
				{
					std::size_t nmbValues = args[cascadeSumNmbElements | 1];
					std::size_t nmbSums   = nmbValues / _nmbOfSumsAtEachStage;
					if (nmbValues % _nmbOfSumsAtEachStage > 0)
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
					// sum _nmbOfSumsAtEachStage elements and store result in _partialSums
					if (_countElements < _nmbOfSumsAtEachStage)
						_partialSums[_partialSumsIndex] += args[sample];
					else {
						_countElements = 0;
						++_partialSumsIndex;
						if (_partialSumsIndex >= _partialSums.size())
							_partialSums.resize(_partialSumsIndex, sampleT());
						_partialSums[_partialSumsIndex] = args[sample];
					}
					++_countElements;
				}
    
				result_type
				result(dont_care) const  // passed argument pack is not used
				{
					// sum up partial sums
					const std::vector<sampleT>* sumsPrev       = &_partialSums;
					std::vector<sampleT>*       sumsNext       = 0;
					bool                        firstIteration = true;
					do {
						std::size_t nmbSumsNext = sumsPrev->size() / _nmbOfSumsAtEachStage;
						if (sumsPrev->size() % _nmbOfSumsAtEachStage > 0)
							++nmbSumsNext;
						sumsNext = new std::vector<sampleT>(nmbSumsNext, sampleT());
						for (std::size_t i = 0; i < nmbSumsNext; ++i) {
							for (std::size_t j = i * _nmbOfSumsAtEachStage;
							     j < min((i + 1) * _nmbOfSumsAtEachStage, sumsPrev->size()); ++j) {
								(*sumsNext)[i] += (*sumsPrev)[j];
							}
						}
						// cleanup
						if (firstIteration)
							firstIteration = false;
						else
							delete sumsPrev;
						// prepare for next iteration
						sumsPrev = sumsNext;
					} while (sumsPrev->size() > 1);
					return (*sumsPrev)[0];
				}
				

			private:

				size_t               _countElements;
				size_t               _nmbOfSumsAtEachStage;
				size_t               _partialSumsIndex;
				std::vector<sampleT> _partialSums;
				
			};


			// Kahan sum reduces (random) round-off errors form O[eps * sqrt(n)] (for naive accumulation)
			// to O[eps] at the cost of 4 summations (instead of 1) for each entry; however, no additional
			// buffer memory is required
			// see http://en.wikipedia.org/wiki/Kahan_summation
			template<typename sampleT>
			struct kahanSumAccumulator : accumulator_base	{

				typedef sampleT result_type;
    
				// constructor takes an Boost.Paramter argument pack there might be an initial value in
				// the argument pack; 'sample' is defined in sample.hpp
				template<typename argsT>
				kahanSumAccumulator(const argsT& args)
					: _sum(args[sample | sampleT()]),
					  _compensation(sampleT())
				{	}
    
				template<typename argsT>
				void
				operator ()(const argsT& args)  // accumulate function also accepts an argument pack
				{
					// at the beginning _compensation is zero
					const sampleT x = args[sample] - _compensation;
					// _sum is big, x small, so low-order digits of x are lost
					const sampleT y = _sum + x;
					// (y - _sum) recovers the high-order part of x;
					// subtracting x recovers -(the lost low part of x)
					_compensation = (y - _sum) - x;
					_sum = y;  // beware eagerly optimising compilers!
					// in the next call the lost low part stored in _compensation will be added to x
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

			struct cascadeSum : depends_on<>,
			                    tag::cascadeSumNmbElements,
			                    tag::cascadeSumNmbOfSumsAtEachStage
			{
				typedef accumulators::impl::cascadeSumAccumulator<mpl::_1> impl;
#ifdef BOOST_ACCUMULATORS_DOXYGEN_INVOKED
				static const parameter::keyword<tag::cascadeSumNmbElements>          nmbElements;
				static const parameter::keyword<tag::cascadeSumNmbOfSumsAtEachStage> nmbOfSumsAtEachStage;
#endif
			};


			struct kahanSum : depends_on<>
			{
				typedef accumulators::impl::kahanSumAccumulator<mpl::_1> impl;
			};

		}


		////////////////////////////////////////////////////////////////////////////
		// extractors for features
		namespace extract {

			const extractor<tag::cascadeSum> cascadeSum = {};
			BOOST_ACCUMULATORS_IGNORE_GLOBAL(cascadeSum)


			const extractor<tag::kahanSum> kahanSum = {};
			BOOST_ACCUMULATORS_IGNORE_GLOBAL(kahanSum)

		}
		using extract::cascadeSum;
		using extract::kahanSum;


		////////////////////////////////////////////////////////////////////////////
		// register feature variants

		// sum(cascade) -> cascadeSum
		struct cascade {};
		template<>
		struct as_feature<tag::sum(cascade)>
		{
			typedef tag::cascadeSum type;
		};
    // for feature-based dependency resolution: provides the same feature as sum
		template<>
		struct feature_of<tag::cascadeSum> : feature_of<tag::sum>
    { };

		
		// sum(kahan) -> kahanSum
		struct kahan {};
		template<>
		struct as_feature<tag::sum(kahan)>
		{
			typedef tag::kahanSum type;
		};
    // for feature-based dependency resolution: provides the same feature as sum
		template<>
		struct feature_of<tag::kahanSum> : feature_of<tag::sum>
    { };

		
	}  // namespace accumulators
}  // namespace boost


#endif  // SUMACCUMULATORS_HPP


