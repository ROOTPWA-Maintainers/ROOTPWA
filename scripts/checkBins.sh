#!/bin/bash

# script checks pwa directory structure for consitency
# Checks: 
#         Equal Number of waves in every bin
#         Equal Number of events in every wave

#loop over mass-bin directories


TXT_BLD=$(tput bold)
TXT_RED=$(tput setaf 1)
TXT_GREEN=$(tput setaf 2)
TXT_YLW=$(tput setaf 3)
TXT_BLUE=$(tput setaf 4)
TXT_PURPLE=$(tput setaf 5)
TXT_CYAN=$(tput setaf 6)
TXT_WHITE=$(tput setaf 7)
TXT_RESET=$(tput sgr0)

rm /tmp/unfinAmps
rm /tmp/unfinPAmps
rm /tmp/unfinAAmps

rm /tmp/zeroAmps
rm /tmp/zeroPAmps
rm /tmp/zeroAAmps

#cd $1;

echo -e "Mass-Bin \t Num \t NumPSP  NumAcc \t FileSize \t PSPFileSize \t AccFileSize \tt NAN"

    


for i in $1*; do

    export NUMAMP=`ls $i/AMPS | grep .amp | wc -l`
    export NUMPSPAMP=`ls $i/PSPAMPS | grep .amp | wc -l`
    export NUMACCAMP=`ls $i/ACCAMPS | grep .amp | wc -l`

    MAXN=0;
    MINN=1000000;
    MAXNMC=0;
    MINNMC=1000000;
    MAXNACC=0;
    MINNACC=1000000;

    if [ $NUMAMP -gt 0 ]; then
    for j in $i/AMPS/*.amp; do
	N=$(stat -c%s "$j"); # use filesize as proxy
	#cat $j | vamp | wc
	#echo $j $N;
        if [ $N -gt $MAXN ]; then MAXN=$N; fi;
	if [ $N -lt $MINN ]; then MINN=$N; fi;
	if [ $N -lt $MAXN -a $N -gt 0 ] ; then ls $j >> /tmp/unfinAmps; fi
	if [ $N -lt $MAXN -a $N -eq 0 ] ; then ls $j >> /tmp/zeroAmps; fi
    done
    fi
    if [ $NUMPSPAMP -gt 0 ]; then
    for j in $i/PSPAMPS/*.amp; do
	NMC=$(stat -c%s "$j"); # use filesize as proxy
	#echo $j $N;
        if [ $NMC -gt $MAXNMC ]; then MAXNMC=$NMC; fi;
	if [ $NMC -lt $MINNMC ]; then MINNMC=$NMC; fi;
	if [ $NMC -lt $MAXNMC -a $NMC -gt 0 ]; then ls $j >> /tmp/unfinPAmps; fi
	if [ $NMC -lt $MAXNMC -a $NMC -eq 0 ]; then ls $j >> /tmp/zeroPAmps; fi
    done
    fi
    if [ $NUMACCAMP -gt 0 ]; then
    for j in $i/ACCAMPS/*.amp; do
	NACC=$(stat -c%s "$j"); # use filesize as proxy
	#echo $j $N;
        if [ $NACC -gt $MAXNACC ]; then MAXNACC=$NACC; fi;
	if [ $NACC -lt $MINNACC ]; then MINNACC=$NACC; fi;
	if [ $NACC -lt $MAXNACC -a $NACC -gt 0 ]; then ls $j >> /tmp/unfinAAmps; fi
	if [ $NACC -lt $MAXNACC -a $NACC -eq 0 ]; then ls $j >> /tmp/zeroAAmps; fi
    done
fi

    if [[ $NUMAMP != $NUMPSPAMP ]]; then 
       
	echo -e "${TXT_RED}$i \t $NUMAMP \t $NUMPSPAMP \t $NUMACCAMP \t\t $MINN..$MAXN \t $MINNMC..$MAXNMC \t $MINNACC..$MAXNACC ${TXT_RESET}\c";
    elif [[ $MAXNMC != $MINNMC ]]; then 
       
	echo -e "${TXT_RED}$i${TXT_RESET} \t $NUMAMP \t $NUMPSPAMP \t $NUMACCAMP \t\t $MINN..$MAXN \t ${TXT_RED}$MINNMC..$MAXNMC ${TXT_RESET} \t $MINNACC..$MAXNACC\c";
    elif [[ $MAXN != $MINN ]]; then 
	echo -e "${TXT_RED}$i${TXT_RESET} \t $NUMAMP \t $NUMPSPAMP \t $NUMACCAMP \t\t ${TXT_RED}$MINN..$MAXN${TXT_RESET} \t $MINNMC..$MAXNMC \c";
    else echo -e "$i \t $NUMAMP \t $NUMPSPAMP \t $NUMACCAMP \t\t $MINN..$MAXN \t $MINNMC..$MAXNMC \t $MINNACC..$MAXNACC \c"
    fi;	
    
   

    if [[ -e $i/AMPS/norm.int ]]; then
	if grep -q nan $i/AMPS/norm.int ; then echo "${TXT_RED}!!${TXT_RESET}";
	else echo -e "\t ok";
	fi
    else echo -e "\t no norm";
    fi

    

done; # and loop over bins

echo "Unfinshed Files (remove if you want to recalculate them!)" 
cat /tmp/unfinAmps
cat /tmp/unfinPAmps
cat /tmp/unfinAAmps
   
echo "Zero size Files" 
cat /tmp/zeroAmps
cat /tmp/zeroPAmps
cat /tmp/zeroAAmps
   

cd -;