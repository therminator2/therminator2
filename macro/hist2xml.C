/********************************************************************************
 *                                                                              *
 *             THERMINATOR 2: THERMal heavy-IoN generATOR 2                     *
 *                                                                              *
 * Version:                                                                     *
 *      Release, 2.0.3, 1 February 2011                                         *
 *                                                                              *
 * Authors:                                                                     *
 *      Mikolaj Chojnacki   (Mikolaj.Chojnacki@ifj.edu.pl)                      *
 *      Adam Kisiel         (kisiel@if.pw.edu.pl)                               *
 *      Wojciech Broniowski (Wojciech.Broniowski@ifj.edu.pl)                    *
 *      Wojciech Florkowski (Wojciech.Florkowski@ifj.edu.pl)                    *
 *                                                                              *
 * Project homepage:                                                            *
 *      http://therminator2.ifj.edu.pl/                                         *
 *                                                                              *
 * For the detailed description of the program and further references           *
 * to the description of the model please refer to                              *
 * http://arxiv.org/abs/1102.0273                                               *
 *                                                                              *
 * This code can be freely used and redistributed. However if you decide to     *
 * make modifications to the code, please, inform the authors.                  *
 * Any publication of results obtained using this code must include the         *
 * reference to arXiv:1102.0273 and the published version of it, when           *
 * available.                                                                   *
 *                                                                              *
 ********************************************************************************/

#include <fstream>
#include <iostream>
#include <TH1D.h>
#define _FIGURE_XSD_FILE_ "http://therminator2.ifj.edu.pl/figure.xsd"

using namespace std;

void hist2xml(const char* aFile, TH1D** aHist, int aNHis = 1, int aEvents = 1, int aEntries = 1)
{
  ofstream ofs;

  ofs.open(aFile);
  if((ofs) && (ofs.is_open())) {
    ofs << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << endl;
    ofs << "<FIGURE xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:noNamespaceSchemaLocation=\""<<_FIGURE_XSD_FILE_<<"\">" << endl;
    ofs << "  <LABEL>"<<aHist[0]->GetYaxis()->GetTitle()<<"</LABEL>" << endl;
    ofs << "  <EVENTS>"<<aEvents<<"</EVENTS>" << endl;
    ofs << "  <ENTRIES>"<<aEntries<<"</ENTRIES>" << endl;
    for(Int_t j=0; j<aNHis; j++) {
      ofs << "  <HISTOGRAM>" << endl;
      ofs << "    <PARTICLE entries=\""<<static_cast<int>(aHist[j]->GetEntries())<<"\">";
      ofs << aHist[j]->GetTitle();
      ofs << "</PARTICLE>" << endl;
      ofs << "    <BINS number=\"6\">" << "idx\t";
      ofs << aHist[j]->GetXaxis()->GetTitle() <<'\t'<< aHist[j]->GetXaxis()->GetTitle() <<"_[MIN]\t"<< aHist[j]->GetXaxis()->GetTitle()<<"_[MAX]\t";
      ofs << aHist[j]->GetYaxis()->GetTitle() <<'\t'<< aHist[j]->GetYaxis()->GetTitle() <<"_[ERROR]";
      ofs << "</BINS>" << endl;
      ofs << "    <DATA points=\""<<aHist[j]->GetNbinsX()<<"\">" << endl;
      for(Int_t i=1; i<=aHist[j]->GetNbinsX(); i++) {
        ofs << i <<'\t';
        ofs << aHist[j]->GetBinCenter(i)  <<'\t'<< aHist[j]->GetBinLowEdge(i) <<'\t'<< aHist[j]->GetBinLowEdge(i) + aHist[j]->GetBinWidth(i) <<'\t';
        ofs << aHist[j]->GetBinContent(i) <<'\t'<< aHist[j]->GetBinError(i)   << endl;
      }
      ofs << "    </DATA>" << endl;
      ofs << "  </HISTOGRAM>" << endl;
    }
    ofs << "</FIGURE>" << endl;
    ofs.close();
    cerr << "Info: XML file: " <<aFile<< " has been created" << endl;
  } else {
    cerr << "Unable to open file "<<aFile<< endl;
  }
}

/*! @file hist2xml.C
 * @brief ROOT macro with hist2xml() function.
 * 
 * @def _FIGURE_XSD_FILE_
 * @brief location of the XSD (XML Schema Definition) file.
 */
/*! @fn void hist2xml(const char* aFile, TH1D** aHist, int aNHis = 1, int aEvents = 1, int aEntries = 1)
 *@brief Saves the contents of TH1 type histograms to a given file in the XML format
 @param [in] aFile name of the XML file
 @param [in] aHist list of pointers to TH1 histograms
 @param [in] aNHis number of histograms in the list
 @param [in] aEvents number of events analyzed
 @param [in] aEntries total number of ParticleCoor entries analyzed
 */
