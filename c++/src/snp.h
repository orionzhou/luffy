#include <vector>
#include <map>
#include <Sequence/PolyTable.hpp>
using namespace std;

namespace Sequence {
    class SNP: public PolyTable {
      private:
        unsigned nind, npos;
        bool haveOG;
        unsigned outgroup_index;
        vector<string> labels;
      public:
        SNP (void): PolyTable(), nind(0), npos(0), haveOG(false) {};
        ~SNP (void) {};

        istream& read (istream& s) throw (Sequence::badFormat, std::exception);
        istream& import (istream& s, const string& ifmt);
        istream& import_vcftable (istream& s);

        ostream& print (ostream& o) const;
        ostream& write (ostream& o, const string& ofmt); 
        ostream& write_phylip (ostream& o);
        ostream& write_table (ostream& o);
        ostream& write_structure (ostream& o);

        unsigned get_nind(void);
        unsigned get_npos(void);
        string get_label(const unsigned& idx);
        vector<string> get_labels(void);
        bool get_haveOG(void);
        unsigned get_outgroup_index(void);
        void set_label(unsigned i, const string& label);
        void set_labels(const vector<string>& labelsN);
        void set_outgroup(const unsigned& i); 
        void set_outgroup(const string& og); 

        unsigned get_label_idx(const string& label);
        vector<unsigned> get_labels_idx(const vector<string>& labels);

        void ApplyFreqFilter(unsigned mincount, bool haveOutgroup=false, unsigned outgroup=0);
        void RemoveMultiHits(bool skipOutgroup=false, unsigned outgroup=0);
        void RemoveMissing(bool skipOutgroup=false, unsigned outgroup=0);
        void RemoveAmbiguous(bool skipOutgroup=false, unsigned outgroup=0);

        void RemoveMono(bool skipOutgroup=false);
        void FilterMissing(const int& co_mis, bool skipOutgroup=false);
        void subset(vector<uint32_t>& idxs_ind, vector<uint32_t>& idxs_pos);
    };
}

