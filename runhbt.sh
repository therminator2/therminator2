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
  echo "This script generates two particle corelation with the \"therm2_femto\""
  echo " and makes a HBT radii fit with \"therm2_hbtfit\""
  echo "Usage:"
  echo "runhbt.sh <EVENT_DIR> <EVENT_COUNT> [FEMTO_INI] [HBTFIT_INI] [PPID]"
  echo "runhbt.sh [OPTION]"
  echo "  <EVENT_DIR>		directory with the event*.root files to process"
  echo "  <EVENT_COUNT>		number of event*.root files to process"
  echo "  [FEMTO_INI]		settings file for therm2_hbtfit; default: \"\""
  echo "  [HBTFIT_INI]		settings file for therm2_events; default: \"\""
  echo "  [PPID]		parent's system process ID; default: 0"
  echo "  [OPTION]"
  echo "    -h | --help		this screen"
  echo "This script makes use of global variables TH2_EVENTS_DIR_, TH2_EVENTS_COUNT_,"
  echo "TH2_FEMTO_INI_, TH2_HBTFIT_INI_ and TH2_PPID_"
  echo
  exit 0
}

function f_intro
{
  echo "-=[ runhbt.sh ]=-----------------------------------------------"
  echo "- Started at:     $(date +"%F %T %Z")"
  echo "- <EVENT_DIR>     $TH2_EVENTS_DIR_"
  echo "- <EVENT_COUNT>   $TH2_EVENTS_COUNT_"
  echo "- [FEMTO_INI]     $TH2_FEMTO_INI_"
  echo "- [HBTFIT_INI]    $TH2_HBTFIT_INI_"
  echo "- [PPID]          $TH2_PPID_"
  echo "---------------------------------------------------------------"
}

function f_error_message()
{
  case "$1" in
    1) ERRORMESS_="ERROR #$1: therm2_femto reported error.";;
    2) ERRORMESS_="ERROR #$1: temporary file ($TH2_FEMTO_TMP_) not found..";;
    3) ERRORMESS_="ERROR #$1: therm2_hbtfit reported error.";;
    *) ERRORMESS_="ERROR #$1: unknown type";;
  esac
  echo "-=[ runhbt.sh ]=-----------------------------------------------"
  echo "- $ERRORMESS_"
  echo "------------------------------------------------------------------"
}

function f_global_variables()
{
  if [ -z $TH2_EVENTS_DIR_ ]; then
    if [ -z $1 ]; then
      f_error_message 1
    else
      TH2_EVENTS_DIR_=$1
    fi
  fi

  if [ -z $TH2_EVENTS_COUNT_ ]; then
    TH2_EVENTS_COUNT_=1
    if ! [ -z $2 ]; then
      TH2_EVENTS_COUNT_=$2
    fi
  fi

  if [ -z $TH2_FEMTO_INI_ ]; then
    if ! [ -z $3 ]; then
      TH2_FEMTO_INI_=$3
    fi
  fi

  if [ -z $TH2_HBTFIT_INI_ ]; then
    if ! [ -z $4 ]; then
      TH2_HBTFIT_INI_=$4
    fi
  fi

  if [ -z $TH2_PPID_ ]; then
    TH2_PPID_=0
    if ! [ -z $5 ]; then
      TH2_PPID_=$5
    fi
  fi

  TH2_FEMTO_TMP_="./femto_$TH2_PPID_.tmp"
  TH2_HBTFIT_TMP_="./hbtfit_$TH2_PPID_.tmp"

  if [ -e $TH2_FEMTO_TMP_ ]; then
    rm $TH2_FEMTO_TMP_
  fi 
  if [ -e $TH2_HBTFIT_TMP_ ]; then
    rm $TH2_HBTFIT_TMP_
  fi
}

################################################################################
# MAIN
################################################################################
if [ "$1" == "-h" -o "$1" == "--help" ]; then
  f_script_usage
fi
f_global_variables $1 $2 $3 $4 $5

f_intro

KTBIN_=0
while [ $KTBIN_ -le 3 ]; do
  echo "-=[ runhbt.sh ]=-----------------------------------------------"
  echo "- kT bin:      $KTBIN_"
  echo "---------------------------------------------------------------"
  if ! ./therm2_femto $KTBIN_ $TH2_EVENTS_DIR_ $TH2_EVENTS_COUNT_ $TH2_FEMTO_INI_ $TH2_PPID_; then
    f_error_message 1
  fi
  if ! [ -e $TH2_FEMTO_TMP_ ]; then
    f_error_message 2
  else
    TH2_FEMTO_ROOT_=$(tail -n 1 $TH2_FEMTO_TMP_)
    if ! ./therm2_hbtfit $TH2_FEMTO_ROOT_ $TH2_HBTFIT_INI_ $TH2_PPID_; then
      f_error_message 3
    fi
    TH2_HBTFIT_ROOT_=$(tail -n 1 $TH2_HBTFIT_TMP_)
    root -b -l -q -x "./macro/figure_fitfemto.C(\"$TH2_HBTFIT_ROOT_\")"
  fi
  KTBIN_=$[$KTBIN_+1]
done

root -b -l -q -x "./macro/figure_hbtradii.C(\"$TH2_EVENTS_DIR_\")"

exit 0
