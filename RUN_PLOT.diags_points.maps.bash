#!/bin/bash
##############################################################
# Main script to compute statistics and generate figures for the air-sea interface
#
# 	In this script, the user chooses and assigns the following criteria for collocation: experiments, dates, HPC partition
#
# 	To run this script on the command line: ./NameOfThisScript.bash
#
# CONTRIBUTORS: Katherine E. Lukens             NOAA/NWS/NCEP/EMC, CISESS at U. of Maryland
#               Guillaume Vernieres             NOAA/NWS/NCEP/EMC
#               Kayo Ide                        U. of Maryland
#               Mindo Choi                      NOAA/NWS/NCEP/EMC, SAIC
#  
# OUTPUT: 
#	- Innovation comparisons between 2 experiments
#		- Maps
#
# NOTES: 
#	1. This script runs all jobs in the background. 
#
# HISTORY:
#	2024-05-10	K.E. Lukens	Created
#
##############################################################

#set -x

#==========================================================================================
# BEGIN USER INPUT
#==========================================================================================

#----------------------------------------------
# Set date range over which to run the program
#	For computational efficiency, it is recommended to loop through one month at a time. 
#	Choose one year/month combination and an array of days.
#
# yyyy 	= year
# mm	= month
# ddarr = day array
#----------------------------------------------

	#`````````````````````````````````````
	# Start Date
yyyyS="2021"
mmS="07"
ddS="01"
hhS="00"

	#`````````````````````````````````````
	# End Date
yyyyE="2021"
mmE="07"
ddE="04"
hhE="18"

#----------------------------------------------
# Set paths to experiments
#----------------------------------------------

	# directory where this script is located
dir_home="/scratch1/NCEPDEV/da/Katherine.Lukens/NSST/tools/marineDA/"

	# paths to experiments
	# 	expt0 = control
	#	expt1 = test to compare against control
#dir_expt="/scratch1/NCEPDEV/stmp2/Katherine.Lukens/expts/cp0.no-ocn-da/COMROOT/cp0.no-ocn-da/saved/"
#exptname="cp0.no-ocn-da"
#dir_expt="/scratch2/NCEPDEV/stmp1/Katherine.Lukens/expts/cp0.ocn-da.sst/COMROOT/cp0.ocn-da.sst/saved/"
#exptname="cp0.ocn-da.sst"

#exptname="ctl_Gbranch"
#dir_expt="/scratch2/NCEPDEV/stmp1/Katherine.Lukens/expts/ctl_Gbranch/COMROOT/ctl_Gbranch/saved/"
#exptname="skint_w_ctl_Gbranch"
#dir_expt="/scratch2/NCEPDEV/stmp1/Katherine.Lukens/expts/skint_w_ctl_Gbranch/COMROOT/skint_w_ctl_Gbranch/saved/"
#exptname="sstavg_w_ctl_Gbranch"
#dir_expt="/scratch2/NCEPDEV/stmp1/Katherine.Lukens/expts/sstavg_w_ctl_Gbranch/COMROOT/sstavg_w_ctl_Gbranch/saved/"
#exptname="GEOerrors_inflated"
#dir_expt="/scratch2/NCEPDEV/stmp1/Katherine.Lukens/expts/GEOerrors_inflated/COMROOT/GEOerrors_inflated/saved/"
exptname="monitorGEOobs"
dir_expt="/scratch2/NCEPDEV/stmp1/Katherine.Lukens/expts/monitorGEOobs/COMROOT/monitorGEOobs/saved/"

#exptname="S2Smodel_S2Sda_C03"
#dir_expt="/scratch1/NCEPDEV/climate/role.ufscpara/ModelOutput/CycledTests/C03/"

#exptname="coolskin_test.v1.1"
#dir_expt="/scratch1/NCEPDEV/da/Katherine.Lukens/NSST/expts/coolskin_test.v1/coolskin_test.v1.1/COMROOT/coolskin_test.v1.1/"
#exptname="coolskin_test.v1.2"
#dir_expt="/scratch1/NCEPDEV/da/Katherine.Lukens/NSST/expts/coolskin_test.v1/coolskin_test.v1.2/COMROOT/coolskin_test.v1.2/"

	# variable to display
var="sst"
#var="seaice"

	# analysis cycle hours to display
#hour=("00" "06" "12" "18")
hour=("00")

	# observing instruments to display
	#	"all" = all available instruments for $var are used
#obstype=("all" "abi_g16" "abi_g17" "ahi_h08" "avhrr_ma" "avhrr_mb" "viirs_n20" "viirs_npp")
#obstype=("all")
obstype=("abi_g16" "abi_g17" "ahi_h08" "avhrr_ma" "avhrr_mb" "viirs_n20" "viirs_npp")
#obstype=("viirs_npp")

#----------------------------------------------
# Set HPC account, partition, and runtime limit for jobs
#	All variables are strings (use double quotes)
#----------------------------------------------

account="da-cpu"

	# Hera
qos="batch"		# max timelimit = 8 hours | default
#qos="debug"		# max timelimit = 30 minutes | highest priority | only 2 jobs can run at same time
#qos="windfall"		# max timelimit = 8 hours | lowest priority | won't affect FairShare value

partition="hera"	# default for batch, debug, windfall

timelimit="00:30:00"
#timelimit="01:00:00"
#timelimit="02:00:00"
#timelimit="04:00:00"
#timelimit="08:00:00"

ntasks=1

#==========================================================================================
# END USER INPUT
#==========================================================================================
#==========================================================================================
#==========================================================================================
# Set up and run SLURM commands
#==========================================================================================

#----------------------------------------------
# Set working paths
#       Always end path names with a slash "/"

dir_src=${dir_home}"src/"           #location of collocation source code
echo 'SOURCE CODE DIRECTORY = '${dir_src}

#----------------------------------------------
# Create output directory (where figures go) if it doesn't already exist

dir_out=${dir_home}"/output/"

if [[ ! -d ${dir_out} ]] ; then
  mkdir -p ${dir_out}
fi

#----------------------------------------------
# Count number of items in arrays

nhour=${#hour[@]}
nobstype=${#obstype[@]}

#----------------------------------------------

echo "-- Start Date: $yyyyS $mmS $ddS"
echo "-- End Date:   $yyyyE $mmE $ddE"
dateSTART=$yyyyS$mmS$ddS$hhS
dateEND=$yyyyE$mmE$ddE$hhE

cd ${dir_home}

#----------------------------------------------
# Set name of run job script (to pass as argument to run job script)
#       This script actually runs the python code

run_python_code=run.diags_points.maps.job

#----------------------------------------------
# Input arguments for SLURM (sbatch) commands
#
#       Arguments (arg) to customize sbatch command

arg1=${dateSTART}
arg2=${dateEND}
arg3=${dir_out}
arg4=${exptname}
arg5=${dir_expt}
#arg6=${exptname1}
#arg7=${dir_expt1}
arg6=${var}

#----------------------------------------------
# Set up sbatch command(s) and run through loop of computation scripts
#
#       Run command format for each job:
#               sbatch $output_logfile_name $job_name $input_directory $colloc_script $arg1 /
#               $arg2 $arg3 $arg4 $arg5 $arg6 $arg7 $arg8 $partition $timelimit $run_job_script

ihour=0
while [[ $ihour -lt $nhour ]]
do
  iobs=0
  while [[ $iobs -lt $nobstype ]]
  do

    #+++++++++++++++++++++++++++++++++
    # Innovation Statistics
    #+++++++++++++++++++++++++++++++++

    arg7=${hour[$ihour]}
    arg8=${obstype[$iobs]}

    	#---------------------------------
	# Diags
    jname=PLOT.diags_points
    python_code=MAIN_PLOT.diags_points.maps.py

    log=LOG_${jname}_${arg4}_${arg6}_${arg7}_${arg8}_${dateSTART}_${dateEND}            # log name containing output log info
    jlog=${jname}_${arg4}_${arg6}_${arg7}_${arg8}_${dateSTART}_${dateEND} 	      # log name for job in queue

    rm $log                                             # remove old log file

    sbatch --mem=0 --output=${log} --job-name=${jlog} --account=${account} --partition=${partition} --qos=${qos} --time=${timelimit} --ntasks=${ntasks} --export=INDIR=${dir_src},SCRIPT=${python_code},ARG1=${arg1},ARG2=${arg2},ARG3=${arg3},ARG4=${arg4},ARG5=${arg5},ARG6=${arg6},ARG7=${arg7},ARG8=${arg8} ${run_python_code}

    	#---------------------------------
        # Diags Diurnal Cycle
#!!! DOESN"T WORK
#   jname=PLOT.diags_points.diurnal
#    python_code=MAIN_PLOT.diags_points.diurnal_cycle.maps.py
#
#    log=LOG_${jname}_${arg4}_${arg6}_${arg7}_${arg8}_${dateSTART}_${dateEND}            # log name containing output log info
#    jlog=${jname}_${arg4}_${arg6}_${arg7}_${arg8}_${dateSTART}_${dateEND}             # log name for job in queue
#
#    rm $log                                             # remove old log file
#
#    sbatch --mem=0 --output=${log} --job-name=${jlog} --account=${account} --partition=${partition} --qos=${qos} --time=${timelimit} --ntasks=${ntasks} --export=INDIR=${dir_src},SCRIPT=${python_code},ARG1=${arg1},ARG2=${arg2},ARG3=${arg3},ARG4=${arg4},ARG5=${arg5},ARG6=${arg6},ARG7=${arg7},ARG8=${arg8} ${run_python_code}

    #+++++++++++++++++++++++++++++++++

    let iobs=iobs+1
  done
  let ihour=ihour+1
done

#----------------------------------------------
#==============================================

##############################################################
# END
##############################################################
