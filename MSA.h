#ifndef MSA_H
#define MSA_H



#include <string>
#include <vector>
#include <iostream>

#include <fstream>

using namespace std;



bool char_in_string(char c, const string& str);
string apply_mask(string seq, const string& mask);


class MSA
{
    public:
    vector<string> text;
    vector<string> name;
    string alphabet="ACGTN-actgn";
    int length;
    int lines;

    void add_sequence(const string& str);
    void concat_sequence(const string& str,int i);
    bool check_sane() const ;
    void printing() const;
    bool perfect_column(int indice)const;
    MSA get_compacted() const;
    char major_nuc(int indice)const;
    bool perfect_column_quasi(int indice,int max_errors) const;
    MSA get_compacted_quasi(int max_errors ) const;
    void parser_fasta(string file);
    string consensus(int threshold) const;
    string consensus_IG(int threshold) const; //consensus IUPAC + gap
    string conversion_IG(string nucleotides_and_gap) const; //IG = code IUPAC + gap
    MSA apply_mask_MSA(const string& mask) const;
    void count_agreements() const;


    MSA(){
		length=0;
		lines=0;
	}

};





#endif
