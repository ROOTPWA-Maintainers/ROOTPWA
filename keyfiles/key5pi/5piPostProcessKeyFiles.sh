#!/bin/bash
#-------------------------------------------------------------------------
# File and Version Information:
# $Rev::                             $: revision of last commit
# $Author::                          $: author of last commit
# $Date::                            $: date of last commit
#
# Description:
#      sed script that adds keys necessary for amplitude calculation
#
#
# Author List:
#      Boris Grube          TUM            (original author)
#
#
#-------------------------------------------------------------------------


if [[ "${1}" = *sigma0*sigma0*sigma0* ]]
then

    # special treatment for sigma -> sigma sigma decays
		sed --in-place=.bak 's/\(^[[:blank:]]*\)name = "sigma0";/\1name = "sigma0";\n\1massDep : { name = "piPiSWaveAuMorganPenningtonKachaev"; };/' "${1}"
		# use Breit-Wigner for 4-body decay mode
		sed --in-place=.bak 's/^        massDep : { name = "piPiSWaveAuMorganPenningtonKachaev"; };/        massDep : { name = "BreitWigner"; };/' ${1}

else

    # set special mass dependencies for some particles
    # see http://unix.stackexchange.com/questions/14887/the-way-to-use-usr-bin-env-sed-f-in-shebang
		exec sed --in-place=.bak "$(cat <<'EOF'

# sigma
s/\(^[[:blank:]]*\)name = "sigma0";/\1name = "sigma0";\n\1massDep : { name = "piPiSWaveAuMorganPenningtonKachaev"; };/

# rho'
s/\(^[[:blank:]]*\)name = "rhoPrime0";/\1name = "rhoPrime0";\n\1massDep : { name = "rhoPrime"; };/

EOF
)" "${1}"

fi
