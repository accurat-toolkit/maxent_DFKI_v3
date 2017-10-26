#ifndef _main_
#define _main_ 1


#include <string.h>
#include <stdio.h>
#include <list>
#include <cstdio>
#include "maxent.h"
#include <iostream>
#include <fstream>
#include <istream>
#include <vector>
#include <map>
#include <cctype>
#include <math.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <stdlib.h>
#include <ctype.h>

#include "ibm1.cpp"


using namespace std;


Lexicon lex;

const string S_LENGTH = "S_LENGTH";
const string T_LENGTH = "T_LENGTH";
const string DIFF = "DIFF";
const string PERC_TRANS_S = "PERC_TRANS_S";
const string PERC_TRANS_T = "PERC_TRANS_T";
const string NUM_WORDS = "NUM_WORDS";
const string PERC_NO_TRANS_S = "PERC_NO_TRANS_S";
const string PERC_NO_TRANS_T = "PERC_NO_TRANS_T";
const string NUM_NO_TRANS = "NUM_NO_TRANS";
const string FERT1 = "FERT1";
const string FERT2 = "FERT2";
const string FERT3 = "FERT3";
const string HIGH_CONTIG_PERC_S = "HIGH_CONTIG_PERC_S";
const string HIGH_CONTIG_PERC_T = "HIGH_CONTIG_PERC_T";
const string LONG_UNCONNECT_STR_S = "LONG_UNCONNECT_STR_S";
const string LONG_UNCONNECT_STR_T = "LONG_UNCONNECT_STR_T";

//--------------------------------------------------------
// Remove all non alpha-numeric chars from a string (including spaces)
//--------------------------------------------------------
string clean_punctuation(string word)
{
	char a;
	string nword = "";
	int start = 0;
	int end = word.length();
	//remvoe white space from the beginning
	for (int i=0; i<nword.length(); i++){
		a = word.at(i);
		if(isalnum(a))
			break;
		if (a == ' ' || a == '\t')
			start++;
	}

	//remove whitespace from the end
	for (int i=word.length()-1; i>0; i--){
		a = word.at(i);
		if(isalnum(a))
			break;
		if (a == ' ' || a == '\t')
			end--;
	}

	//go through the rest, keeping spaces
	for (int i=start; i<end; i++){
		a = word.at(i);
		if (isalnum(a) || a == ' ')
			nword += word[i];
	}

	return nword;
}

//--------------------------------------------------------
// Cut a sentence into words based on spaces
//--------------------------------------------------------
vector<string> get_words(string line){
	vector<string> words;
	string word = "";
	line += " ";
	int index = -1;
	while (line.find(" ") != string::npos){
		index = line.find(" ");
		word = line.substr(0, index);
		line = line.substr(index+1);
		if (word.length() > 0 )
			words.push_back(word);
	}
	return words;
}


//--------------------------------------------------------
// Convert all characters to lowercase
//--------------------------------------------------------
void to_lowercase(string&s)
{
	/*for (string::iterator i = s.begin(); i != s.end(); ++i)
		*i = tolower(*i);*/
}

//--------------------------------------------------------
// Read in a training or test corpus
//--------------------------------------------------------
void read_in_file(string path, vector<string> &source, string s, int start, int end, bool is_find_par){
	  int str_len = 4096;
	  string path_copy = path;

	  //Concat the path and the source
	  string fullpath = path_copy.insert(path_copy.length(),s);
	  ifstream confileS;
	  vector<string> chunks;
	  int i = -1;
	  confileS.open(fullpath.c_str());

	  //increase until the start point
	  char line[str_len];
	  if (end != -1) {
		  while (!confileS.eof() && i <= start){
				i++;
				confileS.getline(line, str_len);
		  }
	  }
	  else i = -2;

	  //then continue to the end point
	  int num_empty = 0;
	  string l = "";
	  while (!confileS.eof() && i <= end)
	  {
	    char line[str_len];
	    confileS.getline(line, str_len);
	    l = line;
	    if (l == "")
	    	num_empty++;
	    //Sometimes files don't have an eof, so keep track of empty lines
	    if (num_empty == 25)
	    	break;
	    to_lowercase(l);
	    source.push_back(l);

	    if (end != -1)
	    	i++;
	  }
	  confileS.close();
}



//IBM Model 1
void ibm_model1(char* problexTable){
	lex=Lexicon(problexTable);
}

double getProb(string sWord, string tWord){
	return lex.getProb(sWord, tWord);
}




//--------------------------------------------------------
// Gather the features for a set of samples
//--------------------------------------------------------
void run_set(string label, const vector<string>& source, const vector<string>& target, vector<ME_Sample> &samples, ofstream &o1, ofstream &o2, bool find){

	vector<string> s_words;
	vector<string> t_words;
	string s, t = "";

	string ws = "";
	string wt = "";
	double prob = 0.0;
	int maxSize = 100;
	for (int i=0; i<source.size(); i++){
			  s = source.at(i);
			  t = target.at(i);

			  s_words = get_words(s);
			  t_words = get_words(t);
	/*	
		cout << s_words.size()<<endl;
		cout << s << endl<<"+++++++++++++++"<<endl;
		
		cout << t_words.size() <<endl;
		cout << t << endl;
		cout << endl;
	*/	 
		
			  if (s_words.size() > maxSize || s_words.size() == 0) continue;
			  if (t_words.size() > maxSize || t_words.size() == 0) continue;
			
			  double ratio = (double) s_words.size()/t_words.size();
			  //don't even bother with ones that are too different, but only when label is unknown
			   if (find && label == "" && (ratio > 1.25 || ratio < 0.8)){
				   continue;
				}
	/*		
		cout << "source" <<s << endl;
		cout << "target" <<t << endl;
*/
			  //---------------------------------------------
			  //First, determine what it's classified as
			  //---------------------------------------------
			  //cout << "....feature ini...." << endl;
			  ME_Sample featureVector(label);
			  //Add a copy of the source and target sentences for easy retrieval later (I added this to the ME_Sample class)
			  featureVector.set_source(s);
			  featureVector.set_target(t);

			  //This loop calculates features 2-6
			 int contig = 0;
			 int high_contig_s = 0;
			 int high_contig_t = 0;
			 int s_num_words = 0;
			 int t_num_words = 0;
			 int fert = 0;
			 int fert1_s = 0;
			 int fert2_s = 0;
			 int fert3_s = 0;
			 int fert1_t = 0;
			 int fert2_t = 0;
			 int fert3_t = 0;
			 int unconnected = 0;
			 int  high_unconnected_s = 0;
			 int  high_unconnected_t = 0;
			 bool found = false;
     		 int s_digit = 0;
			 double prob = 0.0;

			 // cout << "::"<<s << "||" << t <<"::"<< endl;
			  //go through each word in source
			  for (int k=0; k<s_words.size(); k++){
				  ws = s_words.at(k);
				  fert = 0;
				  found = false;
				  for (int l=0; l<t_words.size(); l++){
					  wt = t_words.at(l);
					  prob = getProb(ws, wt);
					  if (prob > 0.2){
						  found = true;
						  fert++;
						}
					 }
				  if (fert >= fert1_s){
					  fert3_s = fert2_s;
					  fert2_s = fert1_s;
					  fert1_s = fert;
				  } else if (fert >= fert2_s){
					  fert3_s = fert2_s;
					  fert2_s = fert;
				  } else if (fert >= fert3_s){
					  fert3_s  = fert;
				  }
				  //check numbers...
				  if (isdigit(ws.c_str()[0]) && t.find(ws) != string::npos){
					  found = true;
				  }else {
					  s_digit ++;
				  }
				  //check morph
				  
				  if (found){
					  unconnected = 0;
					  s_num_words++;
					  contig++;
				  } else {
					  contig = 0;
					  unconnected++;
				  }

				  if (unconnected > high_unconnected_s){
					  high_unconnected_s = unconnected;
				  }

				  //keep track of the highest contiguous span
				  if (contig > high_contig_s){
					  high_contig_s = contig;
				  }
			  }
			
			contig =0 ;
			unconnected = 0;
			int t_digit = 0;
		
			for (int k=0; k<t_words.size(); k++){
				wt = t_words.at(k);
				fert = 0;
				found = false;
				for (int l=0; l<s_words.size(); l++){
					ws = s_words.at(l);
					prob = getProb(ws, wt);
			
					if (prob > 0.2){
						found = true;
						fert++;
					}
				}
				if (fert >= fert1_t){
					fert3_t = fert2_t;
					fert2_t = fert1_t;
					fert1_t = fert;
				} else if (fert >= fert2_t){
					fert3_t = fert2_t;
					fert2_t = fert;
				} else if (fert >= fert3_t){
					fert3_t  = fert;
				}
			//check numbers...
				if (isdigit(wt.c_str()[0]) && s.find(wt) != string::npos){
					found = true;
				} else {
					t_digit ++;
				}
				//check morph
			
				if (found){
					unconnected = 0;
					t_num_words++;
					contig++;
				} else {
						contig = 0;
						unconnected++;
						}
			
			if (unconnected > high_unconnected_t){
				high_unconnected_t = unconnected;
			}
			
			//keep track of the highest contiguous span
			if (contig > high_contig_t){
				high_contig_t = contig;
			}
		}
		
			  //Features are named for the Munteanu and Marcu 2004 paper
			  //cout << "....other features..." << endl;
			  //---------------------------------------------
			  //	FEATURE 0: UNMACTCHED DIGITAL NUMBER
			  //---------------------------------------------
				  featureVector.add_feature("UNMATCH_NUMBER_S",s_digit);
				  featureVector.add_feature("UNMATCH_NUMBER_T",t_digit );
			  //---------------------------------------------
			  //	FEATURE 1: SENTENCE LENGTHS
			  //---------------------------------------------
			  featureVector.add_feature(S_LENGTH, s_words.size());
			  featureVector.add_feature(T_LENGTH, t_words.size());
			  featureVector.add_feature(DIFF, fabs((double) t_words.size() / s_words.size()));
		//cout << "__________________________"<<endl;
		//cout << s << endl;
		//cout << t << endl;
			  //---------------------------------------------
			  //    FEATURE 2: % WORDS WITH TRANSLATIONS
			  //---------------------------------------------
			  featureVector.add_feature("PERC_TRANS_S", (double) s_num_words / s_words.size());
			  featureVector.add_feature("PERC_TRANS_T", (double) t_num_words / t_words.size());
			  featureVector.add_feature("NUM_WORDS_S" , (double) s_num_words);
			  featureVector.add_feature("NUM_WORDS_T" , (double) t_num_words);
		//cout << "NUM_WORDS_S " <<s_num_words <<endl;
		//cout << "NUM_WORDS_T " <<t_num_words <<endl;
			  //-----------------------------------------------
			  //FEATURE 3: % and # OF WORDS WITH NO CONNECTION
			  //-----------------------------------------------
			  double s_inv = s_words.size() - s_num_words;
			  double t_inv = t_words.size() - t_num_words;
			//  featureVector.add_feature("NUM_NO_TRANS_S" , s_inv);
			 // featureVector.add_feature("NUM_NO_TRANS_T" , t_inv);
			  featureVector.add_feature("PERC_NO_TRANS_S", (double)s_inv / s_words.size());
			  featureVector.add_feature("PERC_NO_TRANS_T", (double)t_inv / t_words.size());
		
		//---------------------------------------------
		//FEATURE 4: TOP THREE LARGEST FERTILITIES
		//---------------------------------------------
		featureVector.add_feature("FERT1_S", (double)fert1_s /s_words.size() );
		featureVector.add_feature("FERT2_S", (double)fert2_s /s_words.size());
		featureVector.add_feature("FERT3_S", (double)fert3_s /s_words.size());
		featureVector.add_feature("FERT1_T", (double)fert1_t /t_words.size());
		featureVector.add_feature("FERT2_T", (double)fert2_t /t_words.size());
		featureVector.add_feature("FERT3_T", (double)fert3_t /t_words.size());
		
		//cout << "FERT1_S " << fert1_s <<endl;
		//cout << "FERT2_S " << fert2_s <<endl;
		//cout << "FERT3_S " << fert3_s <<endl;
		//cout << "FERT1_T " << fert1_t <<endl;
		//cout << "FERT2_T " << fert2_t <<endl;
		//cout << "FERT3_T " << fert3_t <<endl;
		//---------------------------------------------
		//FEATURE 5: LENGTH OF LONGEST CONTIGUOUS SPAN
		//---------------------------------------------
		featureVector.add_feature("HIGH_CONTIG_SPAN_S", high_contig_s);
		featureVector.add_feature("HIGH_CONTIG_SPAN_T", high_contig_t);
		featureVector.add_feature("HIGH_CONTIG_PERC_S", (double)high_contig_s/s_words.size());
		featureVector.add_feature("HIGH_CONTIG_PERC_T", (double)high_contig_t/t_words.size());
		//cout << "HIGH_CONTIG_SPAN_S " << high_contig_s<<endl;
		//cout << "HIGH_CONTIG_SPAN_T " << high_contig_t<<endl;
		
		//---------------------------------------------
		//FEATURE 6: LENGTH OF LONGEST UNCONNECTED SUBSTRING
		//---------------------------------------------
		featureVector.add_feature("LONG_UNCONNECT_STR_S", high_unconnected_s);
		featureVector.add_feature("LONG_UNCONNECT_STR_T", high_unconnected_t);
		featureVector.add_feature("LONG_UNCONNECT_PERC_S", (double)high_unconnected_s/s_words.size());
		featureVector.add_feature("LONG_UNCONNECT_PERC_T", (double)high_unconnected_t/t_words.size());
	    
		//cout << "LONG_UNCONNECT_STR_S " << high_unconnected_s<<endl;
		//cout << "LONG_UNCONNECT_STR_T " << high_unconnected_t<<endl;
		//cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" <<endl;
			  //---------------------------------
			  // 		TRAIN THE MODEL
			  //---------------------------------
			  //cout << "...adding to list..."<< endl;
				//cout << "Feature size: "<<featureVector.rvfeatures.size() << endl<< endl;
			  samples.push_back(featureVector);
	}
}



//--------------------------------------------------------
// Get real and fake samples for training
// For fake samples, we add a filler sentence as the first
// item of one of the lists, so everything is shifted by
// 1 item compared to the real set
//--------------------------------------------------------
vector<ME_Sample> get_pos_samples(int start, int end, string dir, string s, string t,  bool is_find_par){
	
	vector<ME_Sample> samples;
	
	vector<string> sourceR;
	vector<string> targetR;
	read_in_file(dir, sourceR, s, start, end, is_find_par);
	read_in_file(dir, targetR, t, start, end, is_find_par);
	ofstream a("");
	ofstream b("");
	run_set("TRANS", sourceR, targetR, samples, a, b, false);
	return samples;
}

vector<ME_Sample> get_neg_samples(int start, int end, string dir, string s, string t,  bool is_find_par){
	
	vector<ME_Sample> samples;
	
	vector<string> sourceR;
	vector<string> targetR;
	read_in_file(dir, sourceR, s, start, end, is_find_par);
	read_in_file(dir, targetR, t, start, end, is_find_par);
	ofstream a("");
	ofstream b("");
	run_set("NOTRANS", sourceR, targetR, samples, a, b, false);
	return samples;
}
//--------------------------------------------------------
// Test the model against some sentences in the corpora
//--------------------------------------------------------
void test_the_model(ME_Model & model, vector<ME_Sample> & samples, int tag){
	
	double found = 0;
	double total = 0;
	string label = "";
	for (int i=0; i < samples.size(); i++){
		ME_Sample s = samples.at(i);
		label = s.label;
		vector<double> vp = model.classify(s);
		/*
		 cout << "SOURCE: "<< s.source_string << endl;
		 cout << "Target: "<< s.target_string << endl; 
		 cout << "value:  "<< vp[0] << endl;
		 cout << "Label: " << s.label << endl<<endl;
		 */
		if (s.label == label)
			found++;
		total++;
		
	}
	string output;
	if (tag == 1 ) output = "Recall: " ;
			else output = "Precision: " ;
	cout << output << found/total << endl;
}

//--------------------------------------------------------
// Add all training samples to the model
//--------------------------------------------------------
void train_the_model_with_pos_samples(ME_Model & model, int start, int end, string path, string s, string t, string dict_path){
	vector<ME_Sample> samples = get_pos_samples(start, end, path, s, t, false);
	for (int i=0; i < samples.size(); i++){
		model.add_training_sample(samples.at(i));
	}
}

void train_the_model_with_neg_samples(ME_Model & model, int start, int end, string path, string s, string t, string dict_path){
	vector<ME_Sample> samples = get_neg_samples(start, end, path, s, t, false);
	for (int i=0; i < samples.size(); i++){
		model.add_training_sample(samples.at(i));
	}
}

int main(int numArgs, char * Args[] )
{
	string source, target, problexTable_str, in_file, out_file, train_switch, TRAIN_SIZE,TEST_SIZE;
	double threshold = 0.8;
	for (int i = 1; i < numArgs ; i ++){
		if (strcmp(Args[i], "--source") == 0)     source = string(Args[++i]);
		else if (strcmp(Args[i], "--target") == 0) target = string(Args[++i]);
		else if (strcmp(Args[i], "--param") == 0) {
			string str;
			for (int j =0 ; j < 4; j++){ 
				str =  string(Args[++i]);
				int pos = (int) str.find("=");
				string leftHand = str.substr(0,pos);
				string rightHand = str.substr(pos+1);
				
				if ( leftHand == "LEX" ) problexTable_str = rightHand;
				else if ( leftHand == "TRAIN" ) train_switch = rightHand;
				else if (leftHand == "TRAIN_SIZE") TRAIN_SIZE = rightHand;
				else if (leftHand == "TEST_SIZE")  TEST_SIZE  = rightHand;
			}
						
		}
		else if (strcmp(Args[i], "--input") == 0)  in_file = string (Args[++i]);
		else if (strcmp(Args[i], "--output") == 0) out_file = string (Args[++i]);
	}
	
	char *problexTable = new char[problexTable_str.length()+1];
	strcpy(problexTable, problexTable_str.c_str());	
	
	cout << "Lexical table Loading..."<<endl;
	ibm_model1(problexTable);
	ME_Model model;
	
	if ( train_switch == "0"){    
		/* EXTRACT */
		cout << "Model Loading..."<<endl;
		model.load_from_file("Script/model");

		ifstream in(in_file.c_str());
		ofstream out(out_file.c_str());

		while (!in.eof()){
			string line;
			getline(in, line);
			int pos1 = (int) line.find("\t");
			int pos2 = (int) line.find("\t", pos1+1);
			
			string infile_E = line.substr(0, pos1);
			string infile_F = line.substr(pos1+1, pos2-pos1-1);
			if (infile_E.length() == 0 ||infile_F.length() == 0) continue;
			
			cout << infile_E << endl;
			cout << infile_F << endl <<endl;
			
			vector<ME_Sample> samples;
			vector<string> set_E,set_F;
			string str_E, str_F;
			
			vector<string> source_set;
			vector<string> target_set;
			read_in_file(infile_E, source_set, "", 0, -1, true);
			read_in_file(infile_F, target_set, "", 0, -1, true);

			for (int i = 0 ; i < source_set.size(); i++)
				for (int j = 0; j < target_set.size(); j++){
					set_E.push_back(source_set[i]);
					set_F.push_back(target_set[j]);
				}
			run_set("", set_E, set_F, samples, out, out, true);
			vector<double> vp;
			for (int i=0; i<samples.size(); i++){
				vp = model.classify(samples[i]);
				if (vp[0] < threshold){
					continue;
				}
				if (samples.at(i).label == "TRANS") {
					string source = samples[i].source_string;
					string target = samples[i].target_string;
					out << source << endl;
					out << target << endl;
					out << vp[0] << endl <<endl;
				}
			}
		}

		in.close();
		out.close();
	} else {
		/*TRAIN*/
		cout << "Training..."<<endl;
		string train_pos = in_file + "/train.pos.";
		string train_neg = in_file + "/train.neg.";
		string test_pos =  in_file + "/test.pos.";
		string test_neg =  in_file + "/test.neg.";
		 
		int train_size = atoi(TRAIN_SIZE.c_str());
		int test_size = atoi(TEST_SIZE.c_str());
		cout << "Loading samples..." <<endl;		
		train_the_model_with_pos_samples(model, 0, train_size, train_pos , source, target, "");
		train_the_model_with_neg_samples(model, 0, train_size, train_neg , source, target, "");
		cout << "Parameterizing..." <<endl;
		model.train();
		vector<ME_Sample> samples = get_neg_samples(0, test_size, test_neg, source, target, false);
		test_the_model(model, samples, -1);
		samples = get_pos_samples(0, test_size, test_pos, source, target, false);
		test_the_model(model, samples, 1);
		//Save for later reference
		model.save_to_file("Script/model");
		
		
	}
}

#endif
