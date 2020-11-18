#!/bin/bash
#********************************************************************************
#*                                                                              *
#*             THERMINATOR 2: THERMal heavy-IoN generATOR 2                     *
#*                                                                              *
#* Version:                                                                     *
#*      Release, 2.0.3, 1 February 2011                                         *
#*                                                                              *
#* Authors:                                                                     *
#*      Mikolaj Chojnacki   (Mikolaj.Chojnacki@ifj.edu.pl)                      *
#*      Adam Kisiel         (kisiel@if.pw.edu.pl)                               *
#*      Wojciech Broniowski (Wojciech.Broniowski@ifj.edu.pl)                    *
#*      Wojciech Florkowski (Wojciech.Florkowski@ifj.edu.pl)                    *
#*                                                                              *
#* Project homepage:                                                            *
#*      http://therminator2.ifj.edu.pl/                                         *
#*                                                                              *
#* For the detailed description of the program and further references           *
#* to the description of the model please refer to                              *
#* http://arxiv.org/abs/1102.0273                                               *
#*                                                                              *
#* This code can be freely used and redistributed. However if you decide to     *
#* make modifications to the code, please, inform the authors.                  *
#* Any publication of results obtained using this code must include the         *
#* reference to arXiv:1102.0273 and the published version of it, when           *
#* available.                                                                   *
#*                                                                              *
#********************************************************************************/

function f_script_usage()
{
  echo "Script runall.sh is the general script for running THERMINATOR2 subscripts:"
  echo " - runevents.sh: generate events"
  echo " - runfigure.sh: create figures form event*.root files"
  echo " - runhbt.sh:    calculate two-particle correlation function, perform HBT fit"
  echo "This script can be executed on a cluster. At this time we prepared it to work with TORQUE."
  echo " - runtorque.sh: configured to execute this script at our cluster at IFJ PAN."
  echo
  echo "Usage:"
  echo "runall.sh [EVENTS_INI] [FEMTO_INI] [HBTFIT_INI]"
  echo "runall.sh [OPTION]"
  echo "  [EVENTS_INI]		settings file for therm2_events; default: \"\""
  echo "  [FEMTO_INI] 		settings file for therm2_femto;  default: \"\""
  echo "  [HBTFIT_INI]		settings file for therm2_hbtfit; default: \"\""
  echo "  [OPTION]"
  echo "    -h | --help		this screen"
  echo
  exit 0
}

function f_intro()
{
  INIT_SEC_=$(date +"%s")
  DATE_=$(date +"%F %T %Z")

  if [ -f /proc/cpuinfo ]; then
    NO_OF_CPU_=$(grep "processor" /proc/cpuinfo | wc -l)
    CPU_MODEL_=$(grep -m 1 "model name" /proc/cpuinfo | sed 's/[model name,:, \t]*//')
  else
    NO_OF_CPU_="N/A"
    CPU_MODEL_="N/A"
  fi
 
  if [ -f /proc/meminfo ]; then
    MEMORY_=$(grep "MemTotal:" /proc/meminfo | sed 's/[MemTotal:, \t]*//')
  else
    MEMORY_="N/A"
  fi

  if which lsb_release &> /dev/null; then
    DISTRO_=$(lsb_release -d | sed 's/[Description:, \t]*//')
  else
    DISTRO_="N/A"
  fi

  if which uname &> /dev/null; then
    KERNEL_=$(uname -r)
  else
    KERNEL_="N/A"
  fi

  echo "-=[ runall.sh ]=--------------------------------------------------"
  echo "- Nr of CPUs:     $NO_OF_CPU_"
  echo "- CPU(s) model:   $CPU_MODEL_"
  echo "- Memory:         $MEMORY_"
  echo "- Distribution:   $DISTRO_"
  echo "- Kernel:         $KERNEL_"
  echo "- User:           $USER@$HOSTNAME"
  echo "- Script PID:     $TH2_PPID_"
  echo "- Started at:     $DATE_"
  echo "- [EVENTS_INI]    $TH2_EVENTS_INI_"
  echo "- [FEMTO_INI]     $TH2_FEMTO_INI_"
  echo "- [HBTFIT_INI]    $TH2_HBTFIT_INI_"
  echo "------------------------------------------------------------------"
}

function f_error_message()
{
  case "$1" in
    1) ERRORMESS_="ERROR #$1: runevents.sh reported error on exit";;
    2) ERRORMESS_="ERROR #$1: temporary file ($TH2_EVENTS_TMP_) not found.";;
    3) ERRORMESS_="ERROR #$1: runfigure.sh reported error on exit";;
    4) ERRORMESS_="ERROR #$1: runhbt.sh reported error on exit";;
    *) ERRORMESS_="ERROR #$1: unknown type";;
  esac
  echo "-=[ runall.sh ]=--------------------------------------------------"
  echo "- $ERRORMESS_"
  echo "------------------------------------------------------------------"
  exit $1
}

function f_global_variables()
{
  TH2_PPID_=$$
  export TH2_PPID_

  if [ -z $TH2_EVENTS_INI_ ]; then
    TH2_EVENTS_INI_="$1"
  fi
  export TH2_EVENTS_INI_

  if [ -z $TH2_FEMTO_INI_ ]; then
    TH2_FEMTO_INI_="$2"
  fi
  export TH2_FEMTO_INI_

  if [ -z $TH2_HBTFIT_INI_ ]; then
    TH2_HBTFIT_INI_="$3"
  fi
  export TH2_HBTFIT_INI_
}

function f_events_temp_file()
{
  TH2_EVENT_TMP_="./events_$TH2_PPID_.tmp"

  if [ -e $TH2_EVENT_TMP_ ]; then
    TH2_EVENTS_DIR_=$(sed -n '1p' $TH2_EVENT_TMP_)
    TH2_EVENTS_COUNT_=$(sed -n '2p' $TH2_EVENT_TMP_)
  else
    f_error_message 2
  fi
# TORQUE: qsub -d path
  if [ ! -z $PBS_O_INITDIR ]; then
    TH2_EVENTS_DIR_=$PBS_O_INITDIR/$TH2_EVENTS_DIR_
  fi

  export TH2_EVENTS_DIR_
  export TH2_EVENTS_COUNT_
}

function f_outro()
{
  DATE_=$(date +"%F %T %Z")
  SEC_=$[$(date +"%s") - $INIT_SEC_]
  DAY_=$[$SEC_ / 86400];	SEC_=$[$SEC_ % 86400]
  HOR_=$[$SEC_ / 3600];		SEC_=$[$SEC_ % 3600]
  MIN_=$[$SEC_ / 60];		SEC_=$[$SEC_ % 60]
  echo "-=[ runall.sh ]=--------------------------------------------------"
  echo "- Finished at:    $DATE_"
  echo "- Total evaluation time: $DAY_ day(s) $HOR_ h $MIN_ m $SEC_ s"
  echo "------------------------------------------------------------------"
}

################################################################################
# MAIN
################################################################################

if [ "$1" == "-h" -o "$1" == "--help" ]; then
  f_script_usage
fi
f_global_variables $1 $2 $3

f_intro
if ! ./runevents.sh; then
  f_error_message 1
fi
f_events_temp_file
if ! ./runfigure.sh; then
  f_error_message 3
fi
if ! ./runhbt.sh; then
 f_error_message 4
fi
f_outro
exit 0
