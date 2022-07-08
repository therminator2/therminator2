#ifndef _TH2_STORAGE_H_
    #define _TH2_STORAGE_H_

    #include <TString.h>
    #include <fstream>

    using namespace std;

    class Storage
    {
        public:
            Storage();
            Storage(TString sFileNam,TString ovrrd = "false");
            ~Storage();

            void appendToTxt(TString sValName, double parVal, double parErr);
            void appendToTxt(TString sValName, double parVal);

        private:
            ofstream *ofStr;
    };

#endif