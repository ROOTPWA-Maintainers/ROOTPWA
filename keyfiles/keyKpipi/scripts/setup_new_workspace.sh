#!/bin/bash

# setting up a new workspace
# author: P.Jasinski Promme@web.de , jasinski@kph.uni-mainz.de

# clear
if [ $# -gt 0 ] # in case of input parameters run the script else give the options
then
	echo " ************** warning ***************** "
	echo " * this script is ment to be run once   * "
	echo " * for a new workspace. Continue (Y/N)? * "
	echo " **************************************** "
	read CONTINUE
	if [ $CONTINUE == "Y" ]
	then
		echo " cleaning up $1 "
		rm -rf $1
		mkdir $1
		echo " compiling tree reader "
		g++ `root-config --glibs` `root-config --cflags` -O0 -o readtree Tree_to_BNL_event.C
		let BINWIDTH=($4-$3)/$2
		echo " creating $2 folder in ${BINWIDTH} (MeV) steps "
		for (( I=$3; $I < $4; I+=$BINWIDTH ))
		do
			BINLOW=$I
			let BINHIGH=$I+$BINWIDTH
			FOLDERNAME=${BINLOW}.${BINHIGH}
        		echo "creating folder ${FOLDERNAME}"
			mkdir $1/${FOLDERNAME}
			cd $1/${FOLDERNAME}
			mkdir "ACCAMPS"
			mkdir "AMPS"
			mkdir "MC"
			mkdir "PSPAMPS" 
			echo " filling with MC pasespace events "
			genpw -n 10000 -M ${BINLOW} -B ${BINWIDTH} -c -r $ROOTPWA/generators/KpipiDiffractive2008.conf
			# BINLOWGEV=`echo "scale=4; ${BINLOW}/1000" | bc`
			# BINHIGHGEV=`echo "scale=4; ${BINHIGH}/1000" | bc`
			# root -l -q -b ${ROOTPWA}/src/pwafitTest/genPhaseSpaceData.C+\(${BINLOWGEV},${BINHIGHGEV},\"\",\"${ROOTPWA}/src/pwafitTest/hTheta.root\",10000,true\)
			cd -
		done
		echo " filling with events "
		./readtree $5 $6 $1 $2 $3 $4 $7
		echo " cleaning up "
		rm -f readtree
		echo " done "
	fi
else
	echo " script to set up a new PWA workspace "
	echo " parameters: <workspace path> <number of bins> <bin min incl. (MeV)> <bin max excl. (MeV)> "
	echo " <path to rootfile> <path to Tree in the rootfile> <number of particles per event>"
fi
