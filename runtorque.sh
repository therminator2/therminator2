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
  echo "Script runtorque.sh executes the \"runall.sh\" script on a cluster."
  echo "This script is configured to run on the IFJ PAN cluster only."
  echo
  echo "Usage:"
  echo "runtorque.sh [JOB_NAME] [QUE_NAME] [EVENTS_INI] [FEMTO_INI] [HBTFIT_INI]"
  echo "runtorque.sh [OPTION]"
  echo "  [JOB_NAME]		job name;			 default: \"therminator2\""
  echo "  [QUE_NAME]		que name;			 default: \"qunlimitted\""
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
  echo "-=[ runtorque.sh ]=-----------------------------------------------"
  echo "- TORQUE at kasztan.ifj.edu.pl"
  echo "- Dir:            $TH2_DIR_"
  echo "- Varlist:        $VAR_LIST_"
  echo "- Started at:     $(date +"%F %T %Z")"
  echo "- [JOB_NAME]      $TH2_JOB_NAME_"
  echo "- [QUE_NAME]      $TH2_QUE_NAME_"
  echo "- [EVENTS_INI]    $TH2_EVENTS_INI_"
  echo "- [FEMTO_INI]     $TH2_FEMTO_INI_"
  echo "- [HBTFIT_INI]    $TH2_HBTFIT_INI_"
  echo "------------------------------------------------------------------"
}

function f_global_variables()
{
  TH2_DIR_=$PWD
  if [ -z $1 ]; then
    TH2_JOB_NAME_="therminator2"
  else
    TH2_JOB_NAME_=$1
  fi
  if [ -z $1 ]; then
    TH2_QUE_NAME_="qunlimitted"
  else
    TH2_QUE_NAME_=$2
  fi
  if [ -n "${ROOTSYS}" ]; then
    VAR_LIST_=ROOTSYS=$ROOTSYS,PATH=$PATH,LD_LIBRARY_PATH=$LD_LIBRARY_PATH
  fi
  VAR_LIST_=$VAR_LIST_,TH2_EVENTS_INI_="$3",TH2_FEMTO_INI_="$4",TH2_HBTFIT_INI_="$5"
}

################################################################################
# MAIN
################################################################################

if [ "$1" == "-h" -o "$1" == "--help" ]; then
  f_script_usage
fi
f_global_variables $1 $2

f_intro
/usr/local/opt/torque/bin/qsub -d $TH2_DIR_ -N $TH2_JOB_NAME_ -q $TH2_QUE_NAME_ -v $VAR_LIST_ ./runall.sh

exit 0
