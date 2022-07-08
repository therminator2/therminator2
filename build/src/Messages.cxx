#include <iostream>
#include "Messages.h"
#include "THGlobal.h"

using namespace std;

void Messages::Intro()
{
    PRINT_MESSAGE("  ***********************************************************************"	);
    PRINT_MESSAGE("  *\t\tTHERMINATOR 2 EVENTS version "<<_THERMINATOR2_VERSION_<<"\t\t\t*"	);
    PRINT_MESSAGE("  *\t\t\t\t\t\t\t\t\t*"							);
    PRINT_MESSAGE("  * authors: M.Chojnacki, A.Kisiel, W.Florkowski, W.Broniowski\t\t*"		);
    PRINT_MESSAGE("  * cite as: arXiv:1102.0273\t\t\t\t\t\t*"					);
    PRINT_MESSAGE("  * webpage: http://therminator2.ifj.edu.pl/\t\t\t\t*"				);
    PRINT_MESSAGE("  ***********************************************************************"	);
}

void Messages::HelpHbt()
{
    PRINT_MESSAGE("Usage: ");
    PRINT_MESSAGE("thermiCorr <FOLDER_PATH> <PAIR_TYPE> <KT_MIN> <KT_MAX> [HBTFIT_INI]");
    PRINT_MESSAGE("thermiCorr [OPTION]");
    PRINT_MESSAGE("  <FOLDER_PATH>\t\tpath to folder with femto ROOT files created with therm2_femto");
    PRINT_MESSAGE("  <PAIR_TYPE>\t\tparticle type created using therm2_femto;\t\toptions: pion-pion or pionM-pionM");
    PRINT_MESSAGE("  <KT_MIN>\t\tstarting kT bin number (from therm2_femto)");
    PRINT_MESSAGE("  <KT_MAX>\t\tfinishing kT bin number (from therm2_femto)");
    PRINT_MESSAGE("  [HBTFIT_INI]\t\tmain settings file;\t\tdefault: hbtfit.ini");
    PRINT_MESSAGE("  [OPTION]");
    PRINT_MESSAGE("    -h | --help\t\tthis screen");
    PRINT_MESSAGE("    -v | --version\tversion information");
}

void Messages::HelpTst()
{
    PRINT_MESSAGE("Usage: ");
    PRINT_MESSAGE("Q2test <FILE_PATH> [Q2TEST_INI]");
    PRINT_MESSAGE("Q2test [OPTION]");
    PRINT_MESSAGE("  <FILE_PATH>\t\tpath to therminacorr output file containing fit parameters");
    PRINT_MESSAGE("  [Q2TEST_INI]\t\tmain settings file;\t\tdefault: statTest.ini");
    PRINT_MESSAGE("  [OPTION]");
    PRINT_MESSAGE("    -h | --help\t\tthis screen");
    PRINT_MESSAGE("    -v | --version\tversion information");
}

void Messages::HelpEvnts()
{
    PRINT_MESSAGE("Usage:");
    PRINT_MESSAGE("therm2_events [EVENTS_INI] [PPID] [HYPER_XML]");
    PRINT_MESSAGE("therm2_events [OPTION]");
    PRINT_MESSAGE("  [EVENTS_INI]\t\tmain settings file;\t\tdefault: events.ini");
    PRINT_MESSAGE("  [PPID]\t\tparent's system process ID;\tdefault: 0");
    PRINT_MESSAGE("  [HYPER_XML]\tlocation of the hypersurface XML file;\tdefault:");
    PRINT_MESSAGE("  [OPTION]");
    PRINT_MESSAGE("    -h | --help\t\tthis screen");
    PRINT_MESSAGE("    -v | --version\tversion information");
}

void Messages::HelpFmt()
{
    PRINT_MESSAGE("Usage: ");
    PRINT_MESSAGE("therm2_femto <KTBIN> <EVENT_DIR> <EVENT_FILES> [FEMTO_INI] [PPID]");
    PRINT_MESSAGE("therm2_femto [OPTION]");
    PRINT_MESSAGE("  <KTBIN>\t\tnumber of the kT bin");
    PRINT_MESSAGE("  <EVENT_DIR>\t\tdirectory with event*.root files");
    PRINT_MESSAGE("  <EVENT_FILES>\t\tnumber of event*.root files to include");
    PRINT_MESSAGE("  [FEMTO_INI]\t\tmain settings file;\t\tdefault: femto.ini");
    PRINT_MESSAGE("  [PPID]\t\tparent's system process ID;\tdefault: 0");
    PRINT_MESSAGE("  [OPTION]");
    PRINT_MESSAGE("    -h | --help\t\tthis screen");
    PRINT_MESSAGE("    -v | --version\tversion information");
}

void Messages::Summary(double low, double high, vector<TString> files, vector<Summ> su)
{
    PRINT_MESSAGE("===============================================================================");

    if(files.empty() && su.empty())
    {
        PRINT_MESSAGE("Everything was successful");
    }
    else
    {
        if(!files.empty())
        {
            PRINT_MESSAGE("Files that were not found: ");
            for(auto it = begin(files); it != end(files); it++)
            {
                PRINT_MESSAGE(*it);
            }
            if(!su.empty())
            {
                PRINT_MESSAGE("");
                PRINT_MESSAGE("");
            }
        }
        if(!su.empty())
        {
            char inv, osl;
            
            PRINT_MESSAGE("Fits outside bounds, where lambda is not in (" << low << ", " << high << "):");
            PRINT_MESSAGE("                   " << "kT bin" << "            " << "Lambda Inv" << "        " << "Lambda OSL");
            for(auto it = begin(su); it != end(su); it++)
            {
                inv = ' ';
                osl = ' ';
                if(it->Linv)
                    inv = '+';
                if(it->Losl)
                    osl = '+';
                PRINT_MESSAGE("                     " << it->kT << "             " << inv << "                   " << osl);
            }
        }
    }

    PRINT_MESSAGE("===============================================================================");
}

void Messages::Version()
{
    PRINT_MESSAGE("version:\t"<<_THERMINATOR2_VERSION_);
    PRINT_MESSAGE("compiled with:");
    PRINT_MESSAGE("\t\t"<<_CXX_VER_);
    PRINT_MESSAGE("\t\tROOT("<<_ROOT_VER_<<")");
    #ifdef _DEBUG_LEVEL_
        PRINT_MESSAGE("DEBUG=\t"<<_DEBUG_LEVEL_);
    #endif
    #ifdef _MODEL_LHYQUID_ONLY_BACK_FLOW_
    std::cout << "BACK_FLOW="<<_MODEL_LHYQUID_ONLY_BACK_FLOW_<<" ";
    #endif
    #ifdef _PARTICLE_DECAYER_RESCALE_CHANNELS_
    std::cout << "RESCALE_CHANNELS="<<_PARTICLE_DECAYER_RESCALE_CHANNELS_<<" ";
    #endif
    #ifdef _PARTICLE_DECAYER_DISABLE_THREE_BODY_DECAYS_
    std::cout << "DISABLE_THREE_BODY_DECAYS="<<_PARTICLE_DECAYER_DISABLE_THREE_BODY_DECAYS_<<" ";
    #endif
    #ifdef _PARTICLE_DECAYER_DISABLE_TWO_BODY_DECAYS_
    std::cout << "DISABLE_TWO_BODY_DECAYS="<<_PARTICLE_DECAYER_DISABLE_TWO_BODY_DECAYS_<<" ";
    #endif
    std::cout << std::endl;
}