#ifndef MSA_H
#define MSA_H



#include <string>
#include <vector>
#include <iostream>



using namespace std;



class MSA
{
    public:
    vector<string> text;
    string alphabet="ACTG-actg";
    int length;
    int lines;
    
    void add_sequence(const string& str);
    void concat_sequence(const string& str,int i);
    bool check_sane() const ;
    void printing() const;
    bool perfect_column(int indice)const;
    MSA get_compacted() const;
    bool perfect_column_quasi(int indice,int max_errors)const;
    MSA get_compacted_quasi(int max_errors ) const;

    MSA(){
		length=0;
		lines=0;
	}

};





#endif
