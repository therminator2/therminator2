#include "TMath.h"

double integrate_mt(TH1 *h, float from, float to, bool useNextWhenEmpty = true) {
    float sum = 0;
    int N = 1;
    int bin1 = h->FindBin(from);
    int bin2 = h->FindBin(to);
    for (int b = bin1; b <= bin2; ++b) {
        float mt = h->GetBinCenter(b) + HPhysicsConstants::mass("pi+");
        float dmt = h->GetBinWidth(b);
        float dn = h->GetBinContent(b);
        sum += N * dn * mt * mt * dmt;
        if (useNextWhenEmpty && (dn <= 0. || !TMath::Finite(dn))) {
            ++N;
        } else {
            N = 1;
        }
    }
    return sum;
}
