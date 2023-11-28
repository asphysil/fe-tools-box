#!/bin/bash
cd `pwd`
for i in 0 1 2 3 4 5 6 7 8 9 10 # 11 12 13 14 15 16 17
#0.0 0.005 0.01 0.015 0.02 0.025 0.03 0.035 0.04 0.045 0.05
#-0.03 -0.02 -0.01 0.0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1
#10 40 80 120 160 200 240 280 320
do
#c0=4.00
#cunit=$(echo "($c0-($d*($c0)))" | bc -l)
#echo $cunit
dir=(4.403 4.436 4.393 4.385 4.374 4.361 4.346 4.323 4.304 4.300 4.288)

a=(3.8364822200000006 3.8344557000000004 3.8443259800000007 3.84892599 3.85411067 3.8597421899999995 3.8656302699999996 3.873680170000001 3.8788753 3.87978016 3.8831147699999997)
b=(3.83641394 3.8399166199999994 3.8444273300000007 3.84921914 3.8541075199999995 3.859770100000001 3.8655818200000005 3.8733455299999995 3.8790866700000004 3.8798922999999994 3.88316942)
c=(4.40363076 4.436717989999999 4.393410009999999 4.38506716 4.374744709999999 4.361817919999999 4.34640727 4.323252010000001 4.304168740000001 4.30042522 4.288226879999999)

#echo ${a[${i}]}


sed -i "14s/.*/REAL(dp_real):: a11=${a[${i}]}, a22=${b[${i}]}, a33=${c[${i}]}  ! a=b=c=3.86 PbTiO3 latice parameters/g" ABO3_dtype.f90
make clean
make
./domain.x
mv data.PTO ../../${dir[${i}]}/
done 