#!/bin/bash
#
#These commands set up the Grid Environment for your job:
#PBS -N Amber
#PBS -l nodes=2:ppn=32,walltime=48:00:00

model=ns5dal-3.pdb
numres=904

#cd $PBS_O_WORKDIR
#!/bin/bash
export AMBERHOME="/opt/amber12-at13"
#export AMBERHOME=/opt/biosoft/amber11
export PATH=${AMBERHOME}/bin:/home/sferanchuk/soft:${PATH}
export LD_LIBRARY_PATH=/opt/cuda-5.5/lib64
#cd /home/sferanchuk/work/screening/pbsa
NP="16"
water_size="9.0"
cut_size="9.0"

sed 's/XXX/MOL/g' gtp.mol2 >lig-s.mol2
antechamber -i lig-s.mol2 -fi mol2 -o ligand.mol2 -fo mol2 
parmchk -i ligand.mol2 -f mol2 -o ligand.frcmod

cat << _EOF1_ >leap.cmd
source leaprc.gaff
MOL = loadmol2 ligand.mol2
loadamberparams ligand.frcmod
saveamberparm MOL lig.prmtop lig.inpcrd
saveoff MOL lig.lib
savepdb MOL l0.pdb
quit
_EOF1_

tleap -s -f leap.cmd

cat $model l0.pdb >mod0.pdb

cat << _EOF2_ >leap2.cmd
source leaprc.ff03.r1
source leaprc.gaff
loadamberparams ligand.frcmod
loadoff lig.lib
mod = loadpdb mod0.pdb
addAtomTypes {
  { "ZN" "Zn" "sp3" }
}
addAtomTypes {
  { "MG" "Mg" "sp3" }
}
# load the frcmd file for Zn2+ after the definition done by "addAtomTypes"
FRCMOD = loadamberparams frcmod.zincmg
loadoff ZN2.lib
loadoff MG0.lib
set mod.903.1 type ZN
set mod.903.1 charge 2.0
set mod.904.1 type ZN
set mod.904.1 charge 2.0
set mod.905.1 type MG
set mod.905.1 charge 2.0
saveamberparm mod model.prmtop model.inpcrd
savepdb mod cmodel.pdb
charge mod
solvatebox mod TIP3PBOX 9.4
saveamberparm mod cmodel.prmtop cmodel.inpcrd
quit
_EOF2_

tleap -s -f leap2.cmd


#NP="8"
water_size="9.0"
cut_size="9.0"
#res_count=`/usr/sbin/pdbconvert count receptor.pdb`


cat << _EOF3_ > sander_en_min1.in
minimise ras-raf
 &cntrl
 imin=1,maxcyc=200,ncyc=100,
 cut=8.0,ntb=1,
 ntc=1,ntf=1,
 ntpr=100,
 ntr=1, restraintmask=':1-$numres',
 restraint_wt=2.0
/
_EOF3_

cmd="mpirun -np $NP $AMBERHOME/bin/sander.MPI -O -i sander_en_min1.in -o cmodel-1.out -c cmodel.inpcrd -ref cmodel.inpcrd -p cmodel.prmtop -r  cmodel_min1.rst -inf md.info -x min1_md.crd "
echo "[$cmd]"
$cmd

cat << _EOF4_ > sander_en_min2.in
minimise ras-raf
 &cntrl
 imin=1,maxcyc=1000,ncyc=500,
 cut=8.0,ntb=1,
 ntc=2,ntf=2,
 ntpr=20,
/
_EOF4_

cmd="mpirun -np $NP $AMBERHOME/bin/sander.MPI -O -i sander_en_min2.in -o cmodel-2.out -c cmodel_min1.rst -p cmodel.prmtop -r  cmodel_min2.rst -inf min2_md.info -x min2_md.crd "
echo "[$cmd]"
$cmd
#cp cmodel_min1.rst cmodel_min2.rst

cat << _EOF5_ > sander_heat.in
heat ras-raf
 &cntrl
  imin=0,irest=0,ntx=1,
  nstlim=25000,dt=0.002,
  ntc=2,ntf=2,
  cut=8.0, ntb=1,
  ntpr=500, ntwx=500,
  ntt=3, gamma_ln=2.0,
  tempi=0.0, temp0=300.0,
  nmropt=1, 
/
 &wt TYPE='TEMP0', istep1=0, istep2=25000,
  value1=0.0, value2=300.0, /
 &wt TYPE='END' /
Restraints for PMEMD
  2.0
  RES 1 $numres
END
END
_EOF5_
# &wt TYPE='TEMP0', istep1=0, istep2=25000,
#  value1=0.0, value2=300.0, /
# &wt TYPE='END' /
#Restraints for PMEMD
#  2.0
#  RES 1 $res_count
#END
#END

#
# run sander energy minimization
#
cmd="mpirun -np $NP $AMBERHOME/bin/sander.MPI -O -i sander_heat.in -o cmodel-3.out -c cmodel_min2.rst -ref cmodel_min2.rst -p cmodel.prmtop -r  cmodel_heat.rst -inf heat_md.info -x heat_md.crd "
echo "[$cmd]"
$cmd

cat << _EOF6_ > sander_dens.in
heat ras-raf
 &cntrl
 imin=0,irest=1,ntx=5,
 nstlim=25000,dt=0.002,
 ntc=2,ntf=2,
 cut=8.0, ntb=2, ntp=1, taup=1.0,
 ntpr=500, ntwx=500,
 ntt=3, gamma_ln=2.0,
 temp0=300.0,
 ntr=1,
/
Restraints for PMEMD
  2.0
  RES 1 $numres
END
END
_EOF6_

cmd="mpirun -np $NP $AMBERHOME/bin/sander.MPI -O -i sander_dens.in -o cmodel-4.out -c cmodel_heat.rst -ref cmodel_heat.rst -p cmodel.prmtop -r  cmodel_dens.rst -inf dens_md.info -x dens_md.crd "
echo [$cmd]
$cmd
