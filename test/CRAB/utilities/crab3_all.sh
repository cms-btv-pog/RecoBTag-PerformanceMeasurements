#!/bin/bash

if [ "$#" -ne 1 ]; then exit; fi;

for dir in crab_*; do

  if [ ! -d "${dir}" ]; then continue; fi;

  crab "$1" "${dir}"

done

unset -v dir
