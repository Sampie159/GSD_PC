#!/bin/bash

programa=./gsd
tamanhos=(small medium large)
threads=(1 2 4 8 16 32)

for tamanho in ${tamanhos[@]}; do
    for thread in ${threads[@]}; do
        echo "Executando $programa com $thread threads e $tamanho tamanho"
        echo "Executando $programa com $thread threads e $tamanho tamanho" >> saidas.txt
        $programa -d $tamanho -t $thread -S 5 >> saidas.txt
    done
done
