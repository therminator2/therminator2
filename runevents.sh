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
  echo "This script generates events with the \"therm2_events\" program"
  echo "Usage:"
  echo "runevents.sh [EVENTS_INI] [PPID]"
  echo "runevents.sh [OPTION]"
  echo "  [EVENTS_INI]		settings file for therm2_events; default: \"\""
  echo "  [PPID]		parent's system process ID; default: 0"
  echo "  [OPTION]"
  echo "    -h | --help		this screen"
  echo "This script makes use of global variables TH2_EVENTS_INI_ and TH2_PPID_"
  echo
  exit 0
}

function f_intro()
{
  echo "-=[ runevents.sh ]=-----------------------------------------------"
  echo "- Started at:     $(date +"%F %T %Z")"
  echo "- [EVENTS_INI]    $TH2_EVENTS_INI_"
  echo "- [PPID]          $TH2_PPID_"
  echo "------------------------------------------------------------------"
}

function f_error_message()
{
  case "$1" in
    1) ERRORMESS_="ERROR #$1: therm2_events reported error.";;
    *) ERRORMESS_="ERROR #$1: unknown type";;
  esac
  echo "-=[ runevents.sh ]=-----------------------------------------------"
  echo "- $ERRORMESS_"
  echo "------------------------------------------------------------------"
  exit $1
}

function f_global_variables()
{
  if [ -z $TH2_EVENTS_INI_ ]; then
    if ! [ -z $1 ]; then
      TH2_EVENTS_INI_=$1
    fi
  fi

  if [ -z $TH2_PPID_ ]; then
    if ! [ -z $2 ]; then
      TH2_PPID_=$2
    fi
  fi
}

################################################################################
# MAIN
################################################################################
if [ "$1" == "-h" -o "$1" == "--help" ]; then
  f_script_usage
fi
f_global_variables $1 $2

f_intro
if ! ./therm2_events $TH2_EVENTS_INI_ $TH2_PPID_; then
  f_error_message 1
fi

exit 0
