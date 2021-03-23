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
  echo "This script creates figures from event files by executing the ROOT macros \"./macro/figure_*.C\"."
  echo "Usage:"
  echo "runfigure.sh <EVENT_DIR> <EVENT_COUNT>"
  echo "runfigure.sh [OPTION]"
  echo "  <EVENT_DIR>		directory with the event*.root files to process"
  echo "  <EVENT_COUNT>		number of event*.root files to process"
  echo "  [OPTION]"
  echo "    -h | --help		this screen"
  echo "This script makes use of global variables TH2_EVENTS_DIR_ and TH2_EVENTS_COUNT_"
  echo
  exit 0
}

function f_intro()
{
  echo "-=[ runfigure.sh ]=-----------------------------------------------"
  echo "- Started at:     $(date +"%F %T %Z")"
  echo "- <EVENT_DIR>     $TH2_EVENTS_DIR_"
  echo "- <EVENT_COUNT>   $TH2_EVENTS_COUNT_"
  echo "------------------------------------------------------------------"
}

function f_error_message()
{
  case "$1" in
    1) ERRORMESS_="ERROR #$1: events directory not defined";;
    *) ERRORMESS_="ERROR #$1: unknown type";;
  esac
  echo "-=[ runfigure.sh ]=-----------------------------------------------"
  echo "- $ERRORMESS_"
  echo "------------------------------------------------------------------"
  exit $1
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
}

################################################################################
# MAIN
################################################################################
if [ "$1" == "-h" -o "$1" == "--help" ]; then
  f_script_usage
fi
f_global_variables $1 $2

f_intro
#root -b -l -q -x "./macro/figure_distpt.C(\"$TH2_EVENTS_DIR_\", $TH2_EVENTS_COUNT_)"
#root -b -l -q -x "./macro/figure_distpt_exotic.C(\"$TH2_EVENTS_DIR_\", $TH2_EVENTS_COUNT_)"
#root -b -l -q -x "./macro/figure_distpt_pion.C(\"$TH2_EVENTS_DIR_\", $TH2_EVENTS_COUNT_)"
#root -b -l -q -x "./macro/figure_distpt_proton.C(\"$TH2_EVENTS_DIR_\", $TH2_EVENTS_COUNT_)"
root -b -l -q -x "./macro/figure_distmt.C(\"$TH2_EVENTS_DIR_\", $TH2_EVENTS_COUNT_)"
#root -b -l -q -x "./macro/figure_v2pt.C(\"$TH2_EVENTS_DIR_\", $TH2_EVENTS_COUNT_)"
#root -b -l -q -x "./macro/figure_v2rap.C(\"$TH2_EVENTS_DIR_\", $TH2_EVENTS_COUNT_)"
#root -b -l -q -x "./macro/figure_distxyt_pion.C(\"$TH2_EVENTS_DIR_\", $TH2_EVENTS_COUNT_)"
#root -b -l -q -x "./macro/figure_disteta.C(\"$TH2_EVENTS_DIR_\", $TH2_EVENTS_COUNT_)"
root -b -l -q -x "./macro/figure_distrap.C(\"$TH2_EVENTS_DIR_\", $TH2_EVENTS_COUNT_)"
root -b -l -q -x "./macro/figure_distrap_prim.C(\"$TH2_EVENTS_DIR_\", $TH2_EVENTS_COUNT_)"
#root -b -l -q -x "./macro/figure_xpx.C(\"$TH2_EVENTS_DIR_\", $TH2_EVENTS_COUNT_)"
#root -b -l -q -x "./macro/figure_hsxyt.C(\"$TH2_EVENTS_DIR_\", $TH2_EVENTS_COUNT_)"
#root -b -l -q -x "./macro/figure_hsret.C(\"$TH2_EVENTS_DIR_\", $TH2_EVENTS_COUNT_)"
exit 0
