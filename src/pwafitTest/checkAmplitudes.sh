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

RESULT=0

echo -e "Directory \t Num \t FileSize \t NAN"

for i in $1; do

    export NUMAMP=`ls $i | grep .amp | wc -l`
 
    MAXN=0;
    MINN=1000000;
    MAXNMC=0;
    MINNMC=1000000;
    MAXNACC=0;
    MINNACC=1000000;

    if [ $NUMAMP -gt 0 ]; then
    for j in $i/*.amp; do
	N=$(stat -c%s "$j"); # use filesize as proxy
	#echo $j $N;
        if [ $N -gt $MAXN ]; then MAXN=$N; fi;
	if [ $N -lt $MINN ]; then MINN=$N; fi;
	if [ $N -lt $MAXN -a $N -gt 0 ] ; then ls $j >> /tmp/unfinAmps; fi
    done
    fi
    if [[ $MAXN != $MINN ]]; then 
	echo -e "${TXT_RED}$i${TXT_RESET} \t $NUMAMP \t ${TXT_RED}$MINN..$MAXN${TXT_RESET}\c";
	RESULT=1;
    else echo -e "$i \t $NUMAMP \t $MINN..$MAXN \c"
    fi;	
    
   

    if [[ -e $i/norm.int ]]; then
	if grep -q nan $i/norm.int ; then echo "${TXT_RED}!!${TXT_RESET}";
	else echo -e "\t ok";
	fi
    else echo -e "\t no norm";
    fi

    

done; # and loop over bins

   
exit $RESULT