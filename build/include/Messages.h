#ifndef _TH2_MESSAGES_H_
    #define _TH2_MESSAGES_H_

    #include <vector>
    #include <TString.h>

    using namespace std;

    class Messages
    {
        public:
            struct Summ
            {
                int kT;
                bool Linv;
                bool Losl;
            };   
            static void Intro();
            static void HelpHbt();
            static void HelpTst();
            static void HelpEvnts();
            static void HelpFmt();
            static void Summary(double low, double high, vector<TString> files, vector<Summ> su);
            static void Version();
            
    };

#endif