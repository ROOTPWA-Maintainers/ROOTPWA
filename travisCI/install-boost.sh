#!/bin/bash
set -exv

# fallback to version 1.63.0 in case no version is given as an argument
BOOST_VERSION="1.63.0"
if [ $# -gt 0 ]
then
	BOOST_VERSION=$1
fi
echo "Using BOOST version ${BOOST_VERSION}"

cd ${TRAVIS_BUILD_DIR}/deps/

if [ -d boost-${BOOST_VERSION} ] ; then
	echo "BOOST installation found in 'boost-${BOOST_VERSION}', using that."
else
	echo "No BOOST installation found, installing a fresh one."

	LIBS="accumulators algorithm align any array assert assign atomic bimap bind chrono concept_check config container conversion core date_time detail dynamic_bitset exception foreach function function_types functional fusion graph graph_parallel integer intrusive io iterator lexical_cast math move mpi mpl multi_array multi_index optional parameter predef preprocessor property_map proto python range ratio rational regex serialization smart_ptr spirit static_assert system test thread throw_exception timer tokenizer tti tuple type_index type_traits typeof unordered utility xpressive"
	TOOLS="build inspect"
	EXTRA="interval numeric_conversion ublas"

	# download required files
	wget https://github.com/boostorg/boost/archive/boost-${BOOST_VERSION}.tar.gz -O boost.tar.gz
	for i in ${LIBS} ${TOOLS} ${EXTRA}
	do
		wget https://github.com/boostorg/${i}/archive/boost-${BOOST_VERSION}.tar.gz -O ${i}.tar.gz
	done

	# extract required libraries and tools
	mkdir boost-${BOOST_VERSION}
	tar -xzf boost.tar.gz -C boost-${BOOST_VERSION} --strip-components=1
	for i in ${LIBS}
	do
		tar -xzf ${i}.tar.gz -C boost-${BOOST_VERSION}/libs/${i}/ --strip-components=1
	done
	for i in ${TOOLS}
	do
		tar -xzf ${i}.tar.gz -C boost-${BOOST_VERSION}/tools/${i}/ --strip-components=1
	done
	tar -xzf interval.tar.gz -C boost-${BOOST_VERSION}/libs/numeric/interval/ --strip-components=1
	tar -xzf numeric_conversion.tar.gz -C boost-${BOOST_VERSION}/libs/numeric/conversion/ --strip-components=1
	tar -xzf ublas.tar.gz -C boost-${BOOST_VERSION}/libs/numeric/ublas/ --strip-components=1

	# compile and install BOOST
	BOOST_ROOT=boost-${BOOST_VERSION} ${TRAVIS_BUILD_DIR}/compileBoostLibraries.sh
fi
ln -sfn boost-${BOOST_VERSION} boost
