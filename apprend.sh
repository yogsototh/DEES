#! /bin/bash

if [ $# -lt 3 ]; then
    echo "utilisation : $0 fichier_automate nombre_d_essais nombre_de_base [nb_repetitions]"
    exit -1
fi

options="-v -p .5 --maxmots 10"
# options="-v -p .2"

if [ $# -gt 3 ]; then
    nombre_repetitions=$4
else
    nombre_repetitions=1
fi

cible=$1
name=`basename $cible`
svgrep="test/svg"

dist=$svgrep/dist-$name
cat /dev/null > $dist
k=1
while [ $k -le $2 ]; do
    i=$(( $k * $3 ))
    realTaille=$(( $i * $2 ))
    for j in `seq $nombre_repetitions`; do
        id=$name-$realTaille-$j
        sample=$svgrep/sample-$id
        log=$svgrep/log-$id
        ma=$svgrep/ma-$id.ma
        ./dees --sample $cible $realTaille $sample
        ./dees --deesha $options -o $ma $sample > $log
        ./dees --dist $cible $ma >> $dist 2>/dev/null
    done
    k=$(( $k + 1 ))
done
