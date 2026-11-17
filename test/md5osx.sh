#!/bin/bash

while read storedhash file
do
    #skip blank lines and comments, like md5sum -c does
    [[ -z $storedhash || $storedhash == \#* ]] && continue
    generatedhash=$(md5 -q "$file")
    if [[ $generatedhash != $storedhash ]]; then
	echo "Hash for file '$file 'does not match"
	exit 1
    fi
done <$2

