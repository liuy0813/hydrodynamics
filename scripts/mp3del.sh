#!/bin/bash
echo "This will delete ALL MP3s in all subdirectories rooted at $(pwd)."
echo "Are you sure (y/N)?: "
read a
if [[ $a == 'y' || $a == 'Y' || $a == 'yes' || $a == 'YES' ]]; then
    if [ ! -f "list.txt" ] ; then
        find . -iname "*.mp3" > list.txt
    fi
    cat list.txt | while read file ; do
        echo "Removing $file"
        rm "$file" ;
    done
fi
if [ -f "list.txt" ] ; then rm list.txt ; fi