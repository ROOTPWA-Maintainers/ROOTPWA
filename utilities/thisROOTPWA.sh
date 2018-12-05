# source this file to set up ROOTPWA


removePath(){
	oldPath="$1"
	toRemove="$2"
	cleanedPath="$(echo "$oldPath" | sed -e "s;^${toRemove}:;;g" -e "s;^${toRemove}$;;g" -e "s;:${toRemove}$;;g" -e "s;:${toRemove}:;:;g")"
	echo "${cleanedPath}"
}


old_ROOTPWA_BUILD="${ROOTPWA_BUILD}"

if ! realpath_loc="$(type -p "realpath")" || [[ -z $realpath_loc ]]; then
	realpath(){
		echo "$1"
	}
fi


ROOTPWA_BUILD="$(dirname "$(dirname "$(realpath "${BASH_SOURCE[0]}")")")"
export ROOTPWA_BUILD

# usually, but not neccessarily
ROOTPWA="$(dirname "${ROOTPWA_BUILD}")"
export ROOTPWA

if [ -n "${old_ROOTPWA_BUILD}" ]; then
	PATH="$(removePath "${PATH}" "${old_ROOTPWA_BUILD}/bin")"
	LD_LIBRARY_PATH="$(removePath "${LD_LIBRARY_PATH}" "${old_ROOTPWA_BUILD}/lib")"
	PYTHONPATH="$(removePath "${PYTHONPATH}" "${old_ROOTPWA_BUILD}/pyLib")"
fi

if [ -z "${PATH}" ]; then
	PATH="${ROOTPWA_BUILD}/bin"
else
	PATH="${ROOTPWA_BUILD}/bin:${PATH}"
fi
export PATH

if [ -z "${LD_LIBRARY_PATH}" ]; then
	LD_LIBRARY_PATH="${ROOTPWA_BUILD}/lib"
else
	LD_LIBRARY_PATH="${ROOTPWA_BUILD}/lib:${LD_LIBRARY_PATH}"
fi
export LD_LIBRARY_PATH




if [ -z "${PYTHONPATH}" ]; then
	PYTHONPATH="${ROOTPWA_BUILD}/pyLib"
else
	PYTHONPATH="${ROOTPWA_BUILD}/pyLib:${PYTHONPATH}"
fi
export PYTHONPATH



unset old_ROOTPWA_BUILD
unset -f removePath
