#!/bin/bash
# This script converts all MP3s in the directory tree routed 'here'
# to OGG vorbis format at 192kbps. It won't delete the original MP3s,
# so that will have to be done manually after it has finished. The
# script will generate its own list of files to process so that it can
# be resumed at any later point if it is stopped.

# Generate a list of all MP3s in the directory tree, or use an already
# existing one (when resuming).
if [ ! -f "list.txt" ] ; then
    find . -iname "*.mp3" > list.txt
fi

# Iterate over every file in the list, converting it to OGG
# This requires first decoding the MP3 with mpg321 into a WAV
# file, and then encoding to OGG (at 192 kbps) with oggenc.
# Finally, when the conversion is done, remove the file from the list.
cat list.txt | while read file ; do
    tempfile=`echo $file | sed 's/mp3/wav/'` ;
    oggfile=`echo $file | sed 's/mp3/ogg/'` ;
    mpg321 "$file" -w "$tempfile" ;
    oggenc "$tempfile" -o "$oggfile" -q 6 ;
    rm "$tempfile" ;
    tail --lines=+2 list.txt > list.tmp ;
    mv list.tmp list.txt ;
done

# Clean up
if [ -f "list.temp" ] ; then rm list.temp ; fi
if [ -f "list.txt" ] ; then
    if [ ! -s "list.txt" ] ; then rm list.txt ; fi
fi

echo "ALL DONE!"