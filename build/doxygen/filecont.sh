#!/bin/bash

function f_dox_content() {
  echo "/*! @page file_content_$2 therminator2/$1"
  echo "@verbatim"
  cat -n ./$1
  echo "@endverbatim"
  echo " */"
}

if [ ! -e $1file_makefile.dox ]; then
  f_dox_content $2Makefile makefile			>> $1file_makefile.dox
fi
if [ ! -e $1file_runall_sh.dox ]; then
f_dox_content $2runall.sh runall_sh			>> $1file_runall_sh.dox
fi
if [ ! -e $1file_runevents_sh.dox ]; then
f_dox_content $2runevents.sh runevents_sh		>> $1file_runevents_sh.dox
fi
if [ ! -e $1file_runfigures_sh.dox ]; then
f_dox_content $2runfigures.sh runfigures_sh		>> $1file_runfigures_sh.dox
fi
if [ ! -e $1file_runhbt_sh.dox ]; then
f_dox_content $2runhbt.sh runhbt_sh			>> $1file_runhbt_sh.dox
fi
if [ ! -e $1file_runtorque.dox ]; then
f_dox_content $2runtorque.sh runtorque_sh		>> $1file_runtorque_sh.dox
fi
if [ ! -e $1file_runevents.dox ]; then
f_dox_content $2runevents.sh runevents_sh		>> $1file_runevents_sh.dox
fi
if [ ! -e $1file_events_ini.dox ]; then
f_dox_content $2events.ini events_ini			>> $1file_events_ini.dox
fi
if [ ! -e $1file_femto_ini.dox ]; then
f_dox_content $2femto.ini femto_ini			>> $1file_femto_ini.dox
fi
if [ ! -e $1file_hbtfit_ini.dox ]; then
f_dox_content $2hbtfit.ini hbtfit_ini			>> $1file_hbtfit_ini.dox
fi
if [ ! -e $1file_blastwave_ini.dox ]; then
f_dox_content $2fomodel/blastwave.ini blastwave_ini	>> $1file_blastwave_ini.dox
fi
if [ ! -e $1file_bwa_ini.dox ]; then
f_dox_content $2fomodel/bwa.ini bwa_ini			>> $1file_bwa_ini.dox
fi
if [ ! -e $1file_krakow_ini.dox ]; then
f_dox_content $2fomodel/krakow.ini krakow_ini		>> $1file_krakow_ini.dox
fi
if [ ! -e $1file_lhyquid2dbi_ini.dox ]; then
f_dox_content $2fomodel/lhyquid2dbi.ini lhyquid2dbi_ini	>> $1file_lhyquid2dbi_ini.dox
fi
if [ ! -e $1file_lhyquid3d_ini.dox ]; then
f_dox_content $2fomodel/lhyquid3d.ini lhyquid3d_ini	>> $1file_lhyquid3d_ini.dox
fi
if [ ! -e $1file_decays_data.dox ]; then
f_dox_content $2share/decays.data decays_data		>> $1file_decays_data.dox
fi
if [ ! -e $1file_particles_data.dox ]; then
f_dox_content $2share/particles.data particles_data	>> $1file_particles_data.dox
fi
