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
    string alphabet="ACTG_actg";
    int length;
    int lines;
    
    void add_sequence(const string& str);
    void concat_sequence(const string& str,int i);
    bool check_sane() const ;
    void printing() const;
    bool perfect_column(int indice)const;
    MSA get_compacted() const;
    MSA(){
		length=0;
		lines=0;
	}

};





#endif
