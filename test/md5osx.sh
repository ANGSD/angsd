#!/bin/bash

while read storedhash file
do
    generatedhash=$(md5 -q "$file")
    if [[ $generatedhash != $storedhash ]]; then
	echo "Hash for file '$file 'does not match"
	exit 1
    fi
done <$2

