#!/bin/bash
#-------------------------------------------------------------------------
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


if [[ "${1}" = *sigma*sigma*sigma* ]]
then

    # special treatment for sigma -> sigma sigma decays
		sed --in-place=.bak 's/\(^[[:blank:]]*\)name = "sigma0";/\1name = "sigma0";\n\1massDep : { name = "piPiSWaveAuMorganPenningtonKachaev"; };/' "${1}"
		# use Breit-Wigner for 4-body decay mode
		sed --in-place 's/^        massDep : { name = "piPiSWaveAuMorganPenningtonKachaev"; };/        massDep : { name = "BreitWigner"; };/' "${1}"

else

    # set special mass dependencies for some particles
		exec sed --in-place=.bak "$(cat <<'EOF'

# sigma
s/\(^[[:blank:]]*\)name = "sigma0";/\1name = "sigma0";\n\1massDep : { name = "piPiSWaveAuMorganPenningtonKachaev"; };/

# rho'
s/\(^[[:blank:]]*\)name = "rhoPrime0";/\1name = "rhoPrime0";\n\1massDep : { name = "rhoPrime"; };/

# f_0(980)
s/\(^[[:blank:]]*\)name = "f0(980)0";/\1name = "f0(980)0";\n\1massDep : { name = "f_0(980)"; };/

EOF
)" "${1}"

fi
