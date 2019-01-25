#!/bin/bash
m0=0;
m0t=0;
iter1=$(awk 'NR==7 {print $1}' c_mass_loop.dat)
i=1
index=0
thmax=1
mfuel=1
conv_line=$(awk 'NR==8 {print $1}' c_mass_loop.dat)
if [ -e interim_mheatshield.dat ] ; then
rm interim_mheatshield.dat
fi
m_heatshield=$(echo print 0.2*$m0 | perl)
echo $m_heatshield > interim_mheatshield.dat
mv interim_mheatshield.dat INPUT/
#---------------------------------------------------------------------------------------------------------------------------
./traj	
#---------------------------------------------------------------------------------------------------------------------------
while [ $index -eq 0 ] ; do
#echo "--------------- Iteration " $i " ------------------------" 
m0minone=$m0
./FindThrustSetting
#.......................................................................
if [ -e interim.dat ] ; then
rm interim.dat
fi 
if [ -e interim2.dat ] ; then
rm interim2.dat 
fi
if [ -e mass.out ] ; then
rm mass.out
fi
#.......................................................................
awk   '{print $4,$5,$26,$30}' thrust.dat >> interim.dat
./find_landing_val.exe >> interim2.dat
#------------------------------------------------------------
./mass_code.exe >> mass.out
m0=$(awk 'NR==1 {print $1}' mass.out)
./m_heatshield.sh
rm mass.out
#echo 
thmax=$(awk 'NR==3 {print $1}' interim2.dat)
mfuel=$(awk 'NR==4 {print $1}' interim2.dat)
#echo " Thrust:     " $thmax
#echo " Mfuel:      " $mfuel
#m0t=$(echo $m0t / 1000 | bc )  
echo $i " " $m0
#echo 
awk -v var="$m0" '{if(NR==1)for(i=1;i<=NF;i++) $1=var; print $0}' INPUT/c_ship.dat > c_ship.dat
mv c_ship.dat INPUT/
#rm interim.dat
#rm interim2.dat
i=$[ $i + 1]
#------------------------------------------------------
#        Convergence Control
#------------------------------------------------------
zero=0
diff=$(echo "$m0 - $m0minone"|bc)
abscheck=$(echo $diff'<'$zero|bc -l)
if [ $abscheck -eq 1 ] ; then
diff=$(echo "$zero - $diff"|bc)
fi
prozdiff=$(echo print $diff/$m0 | perl)
diffcheck=$(echo $prozdiff'<'$conv_line|bc -l)
if [ $diffcheck -eq 1 ] ; then
index=1
fi
#------------------------------------------------------
done
#---------------------------------------------------------------------------------------------------------------------------
