#ifndef _TH2_HBTFIT_H_
    #define _TH2_HBTFIT_H_

    #include <TH3D.h>
    #include <TF3.h>

    class HBTFit 
    {
        public:
            HBTFit();
            ~HBTFit();

            TH3D *getfitprojc(TH3D *expden, TF3 *funqk);
            void preparepad(TH1D *hExp, TH1D* hFit);
            void preparehist(TH1D *hist, int projType, int wType, TString type);
            TH1D *getproj(TH3D *numq, TH3D *denq, int nproj, int wbin, double norm);
            Double_t fungek(Double_t *x, Double_t *par);
            Double_t fungek1D(Double_t *x, Double_t *par);
    };

#endif