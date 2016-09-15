#!/bin/bash
set -exv

# fallback to version 1.61.0 in case no version is given as an argument
BOOST_VERSION="1.61.0"
if [ $# -gt 0 ]
then
	BOOST_VERSION=$1
fi

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
tar -xzf boost.tar.gz
for i in ${LIBS}
do
	tar -xzf ${i}.tar.gz -C boost-boost-${BOOST_VERSION}/libs/${i}/ --strip-components=1
done
for i in ${TOOLS}
do
	tar -xzf ${i}.tar.gz -C boost-boost-${BOOST_VERSION}/tools/${i}/ --strip-components=1
done
tar -xzf interval.tar.gz -C boost-boost-${BOOST_VERSION}/libs/numeric/interval/ --strip-components=1
tar -xzf numeric_conversion.tar.gz -C boost-boost-${BOOST_VERSION}/libs/numeric/conversion/ --strip-components=1
tar -xzf ublas.tar.gz -C boost-boost-${BOOST_VERSION}/libs/numeric/ublas/ --strip-components=1

# compile and install BOOST
BOOST_ROOT=${PWD}/boost-boost-${BOOST_VERSION} ./compileBoostLibraries.sh
