#ifndef MSA_H
#define MSA_H



#include <string>
#include <vector>
#include <iostream>

#include <fstream>

using namespace std;



bool char_in_string(char c, const string& str);
string apply_mask(string seq, const string& mask);
string get_best_consensus_from_prefix(const string& prefix, vector<vector<uint>> matrix, uint& score);
string get_best_consensus(vector<vector<uint>> matrix, uint prefix_size);
vector<string> get_best_consensus_vector(vector<vector<uint>> matrix, uint prefix_size,uint number_consensus);
void reverse_shuffle_haplotype(vector<uint>& result,string& haplotype);
void reverse_shuffle_haplotypes(vector<uint>& result,vector<string>& haplotypes);



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
    void compare_consensus() const;
    void get_diploid();
    vector<vector<uint>> calculates_distance_matrix();
    vector<string> haplotype_merge(vector<string> haplo_snip);
    vector<uint> shuffle_msa();
    void reverse_shuffle_msa(vector<uint>& result);
    vector<uint> order_colonnes(double pb_error);
    MSA cleaning(double pb_error);
    double pb_colonne(int i, int ploidy, double pb_error);
    void clustering(double pb_error);

    MSA(){
		length=0;
		lines=0;
	}

};





#endif
