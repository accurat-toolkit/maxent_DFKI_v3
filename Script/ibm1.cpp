#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <cstdlib>




using namespace std;

class Lexicon{

private:
  int id;
  static const double zero0=1E-5;
public:
  map < pair < string , string > , double > entry;
  set < string > fInLex;

  Lexicon() : id(-1) {}
  Lexicon(char* l) {if(!read(l)) id=(-1);}


  //read GIZA output .ti.actual.final

  bool read(char* pLexTab ){

  fstream pLexFile(pLexTab);

  string puffer;

  while (getline(pLexFile, puffer)){
    istringstream istrLexTab(puffer);

    string tabsrc, tabtgt, tabprob;

    istrLexTab>>tabtgt>>tabsrc>>tabprob;

	if (atof(tabprob.c_str()) > 0.1){
      entry[make_pair(tabtgt, tabsrc)]=atof(tabprob.c_str());
      fInLex.insert(tabsrc);
	}

  }

  pLexFile.close();
  return 1;
  }


  double getProb(string sWord, string tWord){

    map < pair <string, string> , double >::iterator position = entry.find(make_pair(tWord,sWord));
    if (position!=entry.end()){
      if (position->second!=0){
      return position->second;
      }else{
return zero0;
      }
    }else{
      return zero0;
    }
  }
};




