#!/bin/bash


source ./batchUtils.sh


MAX_NMB_JOBS=4
NMB_BINS=50


echo ">>> started ${0} on $(date)"


for i in $(seq 1 ${NMB_BINS})
do
		echo ">>> generating amplitudes for bin ${i}"
		CMD="./calcAmplitudesForMassBin.sh ${i} &> ${PWA_LOGS_DIR}/calcAmplitudes.${i}.log"
    if (( "${MAX_NMB_JOBS}" > 1 ))
    then
        echo "${CMD} &"
        eval "${CMD} &"
        waitForJobs "calcAmplitudesForMassBin.sh" ${MAX_NMB_JOBS}
    else
        echo ${CMD}
        eval ${CMD}
    fi
		
done


echo ">>> finished ${0} on $(date)"
exit 0
