#!/bin/bash


function copyFiles {
		local _FILE_LIST="${1}"
		local _SOURCE="${2}"
		local _DEST="${3}"
		local _FILES=( $(cat "${_FILE_LIST}") )

		for (( IDX=0; IDX<"${#_FILES[@]}"; IDX++ ))
		do
				cp -v "${_SOURCE}/${_FILES[IDX]}" "${_DEST}"
		done
}


function findFiles {
		local _FILE_LIST="${1}"
		local _PATH="${2}"
		local _DEST="${3}"
		local _FILES=( $(cat "${_FILE_LIST}") )

		for (( IDX=0; IDX<"${#_FILES[@]}"; IDX++ ))
		do
				echo "looking for '${_FILES[IDX]}'"
				FILE=$(find "${_PATH}" -name "${_FILES[IDX]}")
				echo "${FILE}"
				cp -v "${FILE}" "${_DEST}"
				echo
		done
}


#copyFiles charly_sym.list ./charly ./charly/sym
#copyFiles charly_nosym.list ./charly ./charly/nosym

findFiles sebastian_sym.list ~/compass/pwa/rootpwa/keyfiles/key5pi ./sebastian/sym
#findFiles sebastian_nosym.list ~/compass/pwa/rootpwa/keyfiles/key5pi


exit 0
