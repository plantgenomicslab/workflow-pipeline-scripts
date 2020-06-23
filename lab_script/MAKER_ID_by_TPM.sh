#!/usr/bin/env bash

# BcaC09g46190.t1
# 123456789012345
# 012345678901234

# set -x
set -e

input=$1

curr_geneid=""
no=1
while read f1 f2; do
  geneid=${f1:0:12}
  if [ "$geneid" != "$curr_geneid" ]; then
    no=1
    curr_geneid=$geneid
  fi
  echo $f1 $f2 $curr_geneid.p$no
  ((no+=1))
done < <(sort -k1.1,1.12 -k2nr $input)
