#!/bin/bash
if [ ! -n "$1" ] && [ -t 0 ]; then
    echo "USAGE: seq2logo 20 1000 <reads.fastq >logo.ps # for length 20 in first 1000 lines" >&2
    echo "OR: seq2logo reads.fastq > logo.ps            # automatic length" >&2
    exit
fi

l=20
n="cat - "
f=""
if [ -n "$1" ]; then
    l="$1"
    [ -f "$1" ] && l=$(grep -v '^[>\+@]' $1 | head | wc -L | awk '{print $1}') && f="$1"
    [ -n "$2" ] && n="head -n$(($2 * 2))"
fi

awk '
/^[>@]/ {
    fastq=sub(/@/, ">");
    print $0;
    getline;
    if(length($0) > '$l') {
        $0 = substr($0, 0, '$l');
        print $0
    } else {
        printf $0;
        for (i=1; i+length($0) <= '$l'; i++)
            printf "-";
        printf "\n"
    } if(fastq==1) {
        getline;
        getline
    }
}' $f | $n | weblogo -c classic
