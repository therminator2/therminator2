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

#include <cstring>
#include <fstream>
#include <sstream>
#include <TString.h>
#include "THGlobal.h"
#include "Hypersurface_Library.h"

using namespace std;

Hypersurface_Library::Hypersurface_Library()
: mXMLFileName(""), mTags(0), mCurrTag(0)
{
}

Hypersurface_Library::~Hypersurface_Library()
{
  if(mTags) {
    mCurrTag = mTags;
    while(mTags->child) { // main XML tag only has children, next=0
      if(mCurrTag->next == 0)
        if(mCurrTag->child == 0)
	  if(mCurrTag->prev == 0) { //last child
            mCurrTag = mCurrTag->father;
	    delete (mCurrTag->child);
	    mCurrTag->child = 0;
	  } else {
	    mCurrTag = mCurrTag->prev;
	    delete (mCurrTag->next);
	    mCurrTag->next = 0;
	  }
        else
	  mCurrTag = mCurrTag->child;
      else
        mCurrTag = mCurrTag->next;
    }
    delete mTags;
  }
}

Hypersurface_Library* Hypersurface_Library::GetXMLTag(const char* aTagName, const char* aAttName, const char* aAttVal) noexcept(false)
{
  XMLtag* Temp;
  
  mCurrTag = mTags;
  if(mCurrTag == 0) {
    PRINT_DEBUG_1("<Hypersurface_Library::GetXMLTag>\tHypersurface_Library is empty!");
    throw _ERROR_LIBRARY_EMPTY_;
  }
  if((mCurrTag->child == 0) && (mCurrTag->name != aTagName) && (!XMLAttributeMatch(aAttName,aAttVal))) {
    PRINT_DEBUG_1("<Hypersurface_Library::GetXMLTag>\tHypersurface_Library with only one entry, unmatched");
    throw _ERROR_LIBRARY_TAG_NOT_FOUND_;
  }

  Temp     = mTags;
  mCurrTag = mTags->child;
  while(mCurrTag != mTags) {
    if((mCurrTag->name == aTagName) && (mCurrTag->attributes.empty() || XMLAttributeMatch(aAttName,aAttVal))) {
      return(this); // mCurrTag is set to the requested TAG. GetXMLAttribute(), GetXMLContant() or GetXMLVector3D() calls available.
    } else 
      if(mCurrTag == Temp)	// moving level up
	if(mCurrTag->next == 0) {
	  mCurrTag = mCurrTag->father;
	  Temp = mCurrTag;
	} else
	  mCurrTag = mCurrTag->next;
      else			// moving level down
        if(mCurrTag->child == 0)
          if(mCurrTag->next == 0) {
	    mCurrTag = mCurrTag->father;
	    Temp = mCurrTag;
          } else
	    mCurrTag = mCurrTag->next;
        else {
          Temp = mCurrTag;
          mCurrTag = mCurrTag->child;
        }
  }
  PRINT_DEBUG_1("<Hypersurface_Library::GetXMLTag>\t<"<<aTagName<<" "<<aAttName<<"=\""<<aAttVal<<"\"> not found.");
  throw _ERROR_LIBRARY_TAG_NOT_FOUND_;
}

TString Hypersurface_Library::GetXMLAttribute(const char* aAttName) noexcept(false)
{
  list<XMLattribute>::iterator Iter;
  
  if(mCurrTag) {
    for(Iter=mCurrTag->attributes.begin(); Iter !=mCurrTag->attributes.end(); Iter++)
      if(Iter->name == aAttName)
        return Iter->value;
    return "";
  } else {
    PRINT_MESSAGE("<Hypersurface_Library::GetXMLAttribute>\tTag attribute "<<aAttName<<" not found");
    throw _ERROR_LIBRARY_TAG_ATTRIB_NOT_FOUND_;
  }
}

TString Hypersurface_Library::GetXMLContent()
{
  TString  buff;
  ifstream ifs;
  char     chr;
  
  buff = "";
  if(mCurrTag) {
    ifs.open(mXMLFileName);
    if((ifs) && ifs.is_open()) {
      ifs.seekg(mCurrTag->begin);
      while(static_cast<long>(ifs.tellg()) != mCurrTag->end) {
	ifs.get(chr);
	buff += chr;
      }
      ifs.close();
    } else {
      PRINT_MESSAGE("<Hypersurface_Library::GetXMLContent>\tFile " << mXMLFileName<< " not opened.");
      exit(_ERROR_GENERAL_FILE_NOT_FOUND_);
    }
  } 
  return buff;
}

Vector3D* Hypersurface_Library::GetXMLVector3D()
{
  ifstream  ifs;
  Vector3D* tVec;
  TString   vName;
  TString   vAxis[3];
  double    vMin[3];
  double    vMax[3];
  int       vPts[3];

  if((mCurrTag) && (mCurrTag->name == "VECTOR3D") ) { // make sure mCurrTag points to a VECTOR3D type Tag
    vName = GetXMLAttribute("name"); mCurrTag = mCurrTag->child; // goto first AXIS
    // extract AXIS Tags
    for(int i=0; i<3; i++) {
      vAxis[i] = GetXMLAttribute("name"); mCurrTag = mCurrTag->child; // goto <DETAIL name="min">
      vMin[i]  = GetXMLContent().Atof(); mCurrTag = mCurrTag->next;  // goto <DETAIL name="max">
      vMax[i]  = GetXMLContent().Atof(); mCurrTag = mCurrTag->next;  // goto <DETAIL name="pts">
      vPts[i]  = GetXMLContent().Atoi(); mCurrTag = mCurrTag->father->next; // goto next AXIS or DATA
    }
    // extract ARRAY Tag -> Vector3D object
    tVec = new Vector3D(vName.Data(),vMin[0],vMax[0],vPts[0],vMin[1],vMax[1],vPts[1],vMin[2],vMax[2],vPts[2]);
    ifs.open(mXMLFileName);
    if((ifs) && ifs.is_open()) {
      ifs.seekg(mCurrTag->begin);
      ifs.precision(16);
      for(int i=0; i<vPts[0]; i++)
        for(int j=0; j<vPts[1]; j++)
	  for(int k=0; k<vPts[2]; k++)
	    ifs >> (*tVec)(i,j,k);
      ifs.close();
      PRINT_DEBUG_1("<Hypersurface_Library::GetXMLVector3D>\tExtracted vector "<<tVec->GetName()<<"("<<vAxis[0]<<"["<<tVec->GetXMin()<<","<<tVec->GetXMax()<<"], "<<vAxis[1]<<"["<<tVec->GetYMin()<<","<<tVec->GetYMax()<<"], "<<vAxis[2]<<"["<<tVec->GetZMin()<<","<<tVec->GetZMax()<<"])");
      return tVec;
    } else {
      PRINT_MESSAGE("<Hypersurface_Library::GetXMLVector3D>\tFile " << mXMLFileName<< " not opened.");
      exit(_ERROR_GENERAL_FILE_NOT_FOUND_);
    }
  } else {
    PRINT_MESSAGE("<Hypersurface_Library::GetXMLVector3D>\tUnable to extract vector. Use GetXMLTag(\"VECTOR3D\",\"name\",\"your_vector_name\").");
    exit(_ERROR_LIBRARY_VECTOR3D_NOT_FOUND_);
  }
}

TString Hypersurface_Library::GetXMLTagDescription()
{
  TString tDesc;
  
  if(mCurrTag == NULL) {
    tDesc = "<NULL />";
    return tDesc;
  }
  
  tDesc = "<" + mCurrTag->name;
  for(list<XMLattribute>::iterator Iter = mCurrTag->attributes.begin(); Iter != mCurrTag->attributes.end(); Iter++)
    tDesc += " " + Iter->name + "=\"" + Iter->value + "\"";
  tDesc += ">\n\t@:";
  tDesc += mCurrTag->begin;
  tDesc += " - ";
  tDesc += mCurrTag->end;
  if(mCurrTag->father != NULL)
    tDesc += "\n\tF:" + mCurrTag->father->name;
  else
    tDesc += "\n\tF:<NULL />";
  if(mCurrTag->child != NULL)
    tDesc += "\n\tC:" + mCurrTag->child->name;
  else
    tDesc += "\n\tC:<NULL />";
  if(mCurrTag->prev != NULL)
    tDesc += "\n\tP:" + mCurrTag->prev->name;
  else
    tDesc += "\n\tP:<NULL />";
  if(mCurrTag->next != NULL)
    tDesc += "\n\tN:" + mCurrTag->next->name;
  else
    tDesc += "\n\tN:<NULL />";
  return tDesc;
}

int Hypersurface_Library::XMLAttributeMatch(const char* aAttName, const char* aAttVal)
{
  for(list<XMLattribute>::iterator Iter = mCurrTag->attributes.begin(); Iter != mCurrTag->attributes.end(); Iter++)
    if(Iter->name.Contains(aAttName) && Iter->value.Contains(aAttVal))
      return 1;
  return 0;
}
