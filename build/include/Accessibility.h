#ifndef _TH2_ACCESSIBILITY_H_
    #define _TH2_ACCESSIBILITY_H_

    #include <TString.h>
    #include "Configurator.h"

    class Accessibility
    {
        public:
            Accessibility();
            ~Accessibility();

            void GetParameter(Configurator* aINI,TString aKeyword, Double_t *parval, Int_t *isfixed, Double_t *parmin, Double_t *parmax);
            void CopyINIFile(TString sMainINI, TString sEventDir);
            static TString getEventDir(TString sFilePath);
    };

#endif
