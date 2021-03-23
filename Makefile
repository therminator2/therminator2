#################################################################################
# Makefile THERMINATOR 2                                       ver. 2011-02-01  #
# Author: Mikolaj Chojnacki (Mikolaj.Chojnacki@ifj.edu.pl)                      #
#################################################################################
#
#helpUsage: (g)make [LABEL] [PREP_DEF]
#help  [LABEL]
#help    [all]     compile THERMINATOR2
#help    doc       generate documentation with Doxygen
#help    package   compress project to tar-ball file
#help    clean     remove object files and binaries
#help    cleandoc  remove Doxygen documentation
#help    help      this screen
#help
#help  [PREP_DEF]  list of preprocesor definitions
#help    DEBUG                     debugging level;                            default: DEBUG=0
#help    RESCALE_CHANNELS          rescale decay channels;                     default: not defined
#help    DISABLE_THREE_BODY_DECAYS disable three body decays;                  default: not defined
#help    DISABLE_TWO_BODY_DECAYS   disable two body decays;                    default: not defined
#help    BACK_FLOW                 particles emitted back to the hydro region; default: not defined
#

#################################################################################
# DEFINITIONS                                                                   #
#################################################################################

# version information form build/include/THGlobal.h
TH2_VERSION = $(shell grep "\#define _THERMINATOR2_VERSION_" $(DIR_H)THGlobal.h | sed 's|[\#,_,2_,a-z,A-Z, ]*||')

# directory structure
DIR_MAIN    = ./
DIR_ADDONS  = $(DIR_MAIN)addons/
DIR_BUILD   = $(DIR_MAIN)build/
DIR_DOXYGEN = $(DIR_BUILD)doxygen/
DIR_H       = $(DIR_BUILD)include/
DIR_OBJ     = $(DIR_BUILD)obj/
DIR_CXX     = $(DIR_BUILD)src/
DIR_DOC     = $(DIR_MAIN)$(shell grep "OUTPUT_DIRECTORY *=" Doxyfile | sed 's|[A-Z,_, \t,=,./]*||')
DIR_DOCHTML = $(DIR_DOC)html/
DIR_DOCTEX  = $(DIR_DOC)latex/
DIR_EVENTS  = $(DIR_MAIN)events/
DIR_FOMODEL = $(DIR_MAIN)fomodel/
DIR_MACRO   = $(DIR_MAIN)macro/
DIR_SHARE   = $(DIR_MAIN)share/
# DIR_UNIGEN  = ${HOME}/unigen/v2.1/
$(info DIR_UNIGEN is $(DIR_UNIGEN))

# search paths
vpath %.h   $(DIR_H)
vpath %.cxx $(DIR_CXX)
vpath %.o   $(DIR_OBJ)

# distribution files
F_INCLUDE   = $(DIR_H)*.h
F_SOURCE    = $(DIR_CXX)*.cxx
F_DOXYGEN   = $(DIR_MAIN)Doxyfile $(DIR_DOXYGEN)*.*
F_MACRO     = $(DIR_MACRO)*.C
F_FOMODEL   = $(DIR_FOMODEL)*.ini $(DIR_FOMODEL)lhyquid2dbi/*.xml
F_SHARE     = $(DIR_SHARE)*.data
F_INI       = $(DIR_MAIN)*.ini
F_BASH      = $(DIR_MAIN)*.sh
F_ADDONS    = $(DIR_ADDONS)*.cxx $(DIR_ADDONS)*.mk
F_PACK      = therminator2_$(TH2_VERSION).tar.gz

# file lists
# THERM2_EVENTS
BIN_EVENTS  = therm2_events
HSRC_EVENTS = Parser.cxx Configurator.cxx ParticleDB.cxx ParticleType.cxx DecayTable.cxx DecayChannel.cxx EventGenerator.cxx Event.cxx Particle.cxx ParticleCoor.cxx Integrator.cxx ParticleDecayer.cxx Model.cxx Model_SR.cxx Model_BlastWave.cxx Model_BWA.cxx Model_KrakowSFO.cxx Model_Lhyquid2DBI.cxx Model_Lhyquid3D.cxx Hypersurface_Lhyquid2D.cxx Hypersurface_Lhyquid3D.cxx Thermodynamics.cxx Chemistry.cxx Energy.cxx Entropy.cxx Pressure.cxx SoundVelocity.cxx Temperature.cxx Viscosity.cxx Hypersurface_Library.cxx Crc32.cxx Vector3D.cxx AbstractEventSaver.cxx RootEventSaver.cxx TextEventSaver.cxx CollectionEventSaver.cxx
ifdef DIR_UNIGEN
HSRC_EVENTS += UnigenEventSaver.cxx
endif
SRC_EVENTS  = $(HSRC_EVENTS:%=$(DIR_CXX)%) $(BIN_EVENTS:%=$(DIR_CXX)%.cxx)
OBJ_EVENTS  = $(SRC_EVENTS:$(DIR_CXX)%.cxx=$(DIR_OBJ)%.o)
# THERM2_FEMTO
BIN_FEMTO   = therm2_femto 
HSRC_FEMTO  = Parser.cxx Configurator.cxx ParticleDB.cxx ParticleType.cxx DecayTable.cxx DecayChannel.cxx
SRC_FEMTO   = $(HSRC_FEMTO:%=$(DIR_CXX)%) $(BIN_FEMTO:%=$(DIR_CXX)%.cxx)
OBJ_FEMTO   = $(SRC_FEMTO:$(DIR_CXX)%.cxx=$(DIR_OBJ)%.o)
# THERM2_HBTFIT
BIN_HBTFIT  = therm2_hbtfit
HSRC_HBTFIT = Parser.cxx Configurator.cxx ParticleDB.cxx ParticleType.cxx DecayTable.cxx DecayChannel.cxx
SRC_HBTFIT  = $(HSRC_HBTFIT:%=$(DIR_CXX)%) $(BIN_HBTFIT:%=$(DIR_CXX)%.cxx)
OBJ_HBTFIT  = $(SRC_HBTFIT:$(DIR_CXX)%.cxx=$(DIR_OBJ)%.o)

# preprocessor
PREPROCESS  = -D_CXX_VER_="\"$(shell $(CXX) --version | grep $(CXX))\"" -D_ROOT_VER_="\"$(shell root-config --version)\""
ifdef DEBUG
  PREPROCESS  := $(PREPROCESS) -D_DEBUG_LEVEL_=$(DEBUG)
else
  PREPROCESS  := $(PREPROCESS) -D_DEBUG_LEVEL_=0
endif
ifdef BACK_FLOW
  PREPROCESS  := $(PREPROCESS) -D_MODEL_LHYQUID_ONLY_BACK_FLOW_=$(BACK_FLOW)
endif
ifdef RESCALE_CHANNELS
  PREPROCESS  := $(PREPROCESS) -D_PARTICLE_DECAYER_RESCALE_CHANNELS_=$(RESCALE_CHANNELS)
endif
ifdef DISABLE_THREE_BODY_DECAYS
  PREPROCESS  := $(PREPROCESS) -D_PARTICLE_DECAYER_DISABLE_THREE_BODY_DECAYS_=$(DISABLE_THREE_BODY_DECAYS)
endif
ifdef DISABLE_TWO_BODY_DECAYS
  PREPROCESS  := $(PREPROCESS) -D_PARTICLE_DECAYER_DISABLE_TWO_BODY_DECAYS_=$(DISABLE_TWO_BODY_DECAYS)
endif

# compilation
CXX         = g++
LD          = g++
CXXFLAGS    = -std=c++11 -O0 -g -Wno-deprecated -Wno-return-type -I $(DIR_H) $(PREPROCESS) `root-config --cflags`
LFLAGS      = -lm -g `root-config --libs`
ifdef DIR_UNIGEN
CXXFLAGS    += -DUSE_UNIGEN -I $(DIR_UNIGEN)/base/include
LFLAGS      += -L $(DIR_UNIGEN)/lib -lUniGen
endif
#################################################################################
# RULES                                                                         #
#################################################################################
all: $(BIN_EVENTS:%=$(DIR_OBJ)%) $(BIN_FEMTO:%=$(DIR_OBJ)%) $(BIN_HBTFIT:%=$(DIR_OBJ)%)
	cp $^ $(DIR_MAIN)
	echo
	echo "Type \"./therm2_events\" to generate events,"
	echo "Type \"./therm2_femto\" to generate two-particle corelation function,"
	echo "Type \"./therm2_hbtfit\" to fit and extract HBT radii."
	echo

$(DIR_OBJ)therm2_events: $(OBJ_EVENTS)
	echo "Linking:   $@ ($(LD))"
	$(LD) $^ -o $@ $(LFLAGS)

$(DIR_OBJ)therm2_femto: $(OBJ_FEMTO)
	echo "Linking:   $@ ($(LD))"
	$(LD) $^ -o $@ $(LFLAGS)

$(DIR_OBJ)therm2_hbtfit: $(OBJ_HBTFIT)
	echo "Linking:   $@ ($(LD))"
	$(LD) $^ -o $@ $(LFLAGS)

$(DIR_OBJ)%.o: %.cxx
	@[ -d $(DIR_OBJ) ] || mkdir -p $(DIR_OBJ)
	echo "Compiling: $< ($(CXX))"
	$(CXX) -c $< -o $@ $(CXXFLAGS)

doc:
	$(DIR_DOXYGEN)filecont.sh $(DIR_DOXYGEN)
	doxygen
	rm -f $(DIR_DOXYGEN)file_*.dox
	cp $(DIR_DOXYGEN)dirstruct.svg $(DIR_DOCHTML).
	cp $(DIR_DOXYGEN)figure.xml $(DIR_DOCHTML).
	cp $(DIR_FOMODEL)lhyquid2dbi/RHICAuAu200c2030Ti455ti025Tf145.xml $(DIR_DOCHTML)rhic.xml
	echo "<!DOCTYPE HTML PUBLIC '-//W3C//DTD HTML 4.0 Transitional//EN'><HTML><HEAD><META http-equiv='REFRESH' content='0;url=$(DIR_DOCHTML)index.html'></HEAD><BODY></BODY></HTML>" > manual.html
	$(MAKE) -C $(DIR_DOCTEX)
	ln -s $(DIR_DOCTEX)refman.pdf manual.pdf
	echo
	echo "THERMINATOR 2 reference manual generated with Doxygen:"
	echo "- manual.html"
	echo "- manual.pdf"

package: $(F_INCLUDE) $(F_SOURCE) $(F_MACRO) $(F_FOMODEL) $(F_SHARE) $(F_DOXYGEN) $(F_INI) $(F_BASH) $(F_ADDONS) Makefile
	echo "$(TH2_VERSION)" > version
	tar zcvf $(F_PACK) $^ version
	echo "Package '$(F_PACK)' created."

clean:
	rm -f $(DIR_OBJ)*.o
	rm -f $(DIR_OBJ)$(BIN_EVENTS) $(DIR_MAIN)$(BIN_EVENTS)
	rm -f $(DIR_OBJ)$(BIN_FEMTO)  $(DIR_MAIN)$(BIN_FEMTO)
	rm -f $(DIR_OBJ)$(BIN_HBTFIT) $(DIR_MAIN)$(BIN_HBTFIT)
	echo "*.o and binary files removed."

cleandoc:
	rm -rf $(DIR_DOC)
	rm -f $(DIR_MAIN)manual.html $(DIR_MAIN)manual.pdf
	echo "Doxygen documentation removed."

help:
	@grep -h "^#help" $(MAKEFILE_LIST) | sed 's|\#help||'

.SILENT :
