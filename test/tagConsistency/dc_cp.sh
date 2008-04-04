#!/bin/bash
for file in `ls $1`; do
  dccp $file $2;
  #echo $file $2;
done
