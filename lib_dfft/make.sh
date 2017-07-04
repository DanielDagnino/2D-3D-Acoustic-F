#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#$ -e make.error
#$ -o make.output
#$ -q new_ifb.q

##------------------------------------------------------------------------------#
#module purge   # treu els mòduls precarregats
#module load intel-2016 openmpi-1.10.3_ifb-intel-2016 

###module avail
##module purge
##module load rocks-openmpi_ib   #
###module load rocks-openmpi      #
###module load openmpi-x86_64     #

#------------------------------------------------------------------------------#
echo 'MPI module'
mpif90 --showme
echo

#------------------------------------------------------------------------------#

make
