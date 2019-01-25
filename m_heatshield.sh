#!/bin/bash

if [ -e interim_mheatshield.dat ] ; then
rm interim_mheatshield.dat
fi
m_heatshield=$(awk 'NR==1 {print $3}' mass.out)
echo $m_heatshield > interim_mheatshield.dat
mv interim_mheatshield.dat INPUT/

