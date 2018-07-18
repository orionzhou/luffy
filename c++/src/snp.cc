#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <map>
#include <cctype>
#include <cassert>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/assign/std/vector.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>
#include <Sequence/PolyTableFunctions.hpp>
#include <Sequence/stateCounter.hpp>
#include "snp.h" 
using namespace std;
using namespace boost::assign;
using std::vector;
using std::string;
using boost::format;

template< class T >
struct next {
  T seed;
  next( T seed ) : seed(seed) { }
  T operator()() {
    return seed++;
  }
};

namespace Sequence {
  unsigned SNP::get_nind(void) { return nind; }
  unsigned SNP::get_npos(void) { return npos; }
  string SNP::get_label(const unsigned& idx) { return labels[idx]; }
  vector<string> SNP::get_labels(void) { return labels; }
  bool SNP::get_haveOG(void) { return haveOG; }
  unsigned SNP::get_outgroup_index(void) { return outgroup_index; }
  
  void SNP::set_outgroup(const unsigned& i) {
    haveOG = true;
    outgroup_index = i;
  }
  void SNP::set_outgroup(const string& og) {
    haveOG = true;
    outgroup_index = this->get_label_idx(og);
  }
  void SNP::set_label(unsigned i, const string& label) { labels[i-1] = label; }
  void SNP::set_labels(const vector<string>& labelsN) {
    if(labels.size() != labelsN.size()) {
      cerr << "label size error" << endl;
      exit(1);
    } else {
      labels = labelsN;
    }
  }
  unsigned SNP::get_label_idx(const string& label) {
    int idx = -1;
    for(unsigned i = 0; i < nind; ++i) {
      if(labels[i] == label) {
        idx = i;
        break;
      }
    }
    if(idx < 0) {
      cerr << "label not find: " << label << endl;
      exit(1);
    }
    return idx;
  }
  vector<unsigned> SNP::get_labels_idx(const vector<string>& labels) {
    vector<unsigned> idxs(labels.size());
    for(unsigned i = 0; i < labels.size(); ++ i)
      idxs[i] = get_label_idx( labels[i] );
    return idxs;
  }
  
  istream& SNP::read (istream& s) throw (Sequence::badFormat, std::exception) {
    string line;
    vector<string> parts;
    unsigned i, j;
    
    getline(s, line);
    boost::trim_if(line, boost::is_any_of("\t "));
    boost::split(parts, line, boost::is_any_of("\t "), boost::token_compress_on);
    if(parts.size() != 2) {
      throw badFormat("snp.cc: line 1 [nind npos] error");
    } else {
      nind = boost::lexical_cast<unsigned>(parts[0]);
      npos = boost::lexical_cast<unsigned>(parts[1]);
    }
    
    getline(s, line);
    boost::trim_if(line, boost::is_any_of("\t "));
    boost::split(parts, line, boost::is_any_of("\t "), boost::token_compress_on);
    vector<double> _poss(npos);
    if(parts.size() != npos) {
      throw badFormat("snp.cc: line 2 [positions] error");
    } else {
      for (j = 0; j < npos; ++j) 
        _poss[j] = boost::lexical_cast<double>(parts[j]);
    }

    labels.resize(nind);
    vector<string> _data(nind);
    for (i = 0; i < nind; ++i) {
      getline(s, line);
      boost::trim_if(line, boost::is_any_of("\t "));
      boost::split(parts, line, boost::is_any_of("\t "), boost::token_compress_on);
      
      string label;
      unsigned haveLabel = 0;
      if(parts.size() == npos) {
        labels[i] = "ind" + boost::lexical_cast<string>(i+1);
      } else if(parts.size() == npos+1) {
        haveLabel = 1;
        labels[i] = parts[0];
      } else {
        throw (Sequence::badFormat("snp.cc::read() not enough alleles"));
      }
      
      char *temp = new char[npos+1];
      for(j = 0; j < npos; ++j) {
        string str = parts[haveLabel+j];
        if(str.size() != 1)
          throw (Sequence::badFormat("ssp.cc:read() allele length error"));
        char ch = toupper(str[0]);
        switch( ch ) {
            case '?':
            temp[j] = 'N';
            break;
            default:
            temp[j] = ch;
            break;
        }
      }
      temp[j] = '\0';
      _data[i] = temp;
      delete [] temp;
    }
    if (_data.size() != nind)
      throw (Sequence::badFormat("snp::read() nind+haveOG != line number"));
    
    PolyTable::assign(&_poss[0], _poss.size(), &_data[0], _data.size());
    RemoveInvariantColumns(this);
    return s;
  }
  istream& SNP::import_vcftable (istream& s) {
    string line;
    vector<string> parts;
    boost::regex re("([A-Za-z0-9_-]+)=([.01])/([.01])");
    boost::match_results<string::const_iterator> what;
    
    bool first = true;
    vector<string> _data;
    vector<double> _poss;
    char ref, alt;
    while( getline(s, line) ) {
      boost::trim_if(line, boost::is_any_of("\t "));
      boost::split(parts, line, boost::is_any_of("\t "), boost::token_compress_on);
      if(first) {
        if(parts.size() <= 5) {
          cerr << "line 1 < 5 columns" << endl;
          exit(1);
        }
        nind = parts.size() - 4;
        labels.resize(nind);
        _data.resize(nind);
      } else {
        if(parts.size() != 4+nind) {
          cerr << "line columns error" << endl;
          exit(1);
        }
      }
      _poss.push_back(boost::lexical_cast<double>(parts[1]));
      assert(parts[2].size() == 1);
      assert(parts[3].size() == 1);
      ref = parts[2][0];
      alt = parts[3][0];
      for(unsigned i = 0; i < nind; ++ i) {
        boost::regex_match(parts[4+i], what, re);
        if(first) {
          labels[i] = what[1];
        } else if(what[1] != labels[i]) {
          cerr << format("column %d[%s] not %s\n") % (4+i) % what[1] % labels[i]; 
          exit(1);
        }
        string ch1 (what[2]), ch2 (what[3]);
        if(ch1 == "0" && ch2 == "0") {
          _data[i] += ref;
        } else if(ch1 == "1" && ch2 == "1") {
          _data[i] += alt;
        } else {
          _data[i] += 'N';
        }
      }
      if(first) first = false;
    }
    
    PolyTable::assign(&_poss[0], _poss.size(), &_data[0], _data.size());
    RemoveInvariantColumns(this);
    npos = this->numsites();
    return s;
  }
  istream& SNP::import (istream& s, const string& ifmt) {
    if(ifmt == "snp") {
      this->read(s);
    } else if (ifmt == "vcftable") {
      this->import_vcftable(s);
    }
    return s;
  }

  ostream& SNP::print (ostream& o) const {
    vector<double> poss = this->GetPositions();
    vector<string> data = this->GetData();
    o << nind << '\t' << format("%.0f") % npos << endl;
    
    for(unsigned j = 0; j < npos-1; j++)
      o << format("%.0f") % poss[j] << '\t';
    o << format("%.0f") % poss[npos-1] << endl;
    
    for(unsigned i = 0; i < nind; i ++) {
      o << labels[i];
      for(unsigned j = 0; j < npos; j ++)
        o << '\t' << data[i][j];
      o << endl;
    }
    return o;
  }
  ostream& SNP::write_phylip (ostream& o) {
    vector<double> poss = this->GetPositions();
    vector<string> data = this->GetData();
    o << format("%d  %d\n") % nind % npos;
    for (uint32_t r = 0; r < ceil(double(npos) / 250); ++r) {
      for (uint32_t i = 0; i < nind; i++) {
        string data_row_string = "";
        for (size_t j = 250*r; j < min(npos, 250*r+250); j++)
          data_row_string += data[i][j];
        if(r == 0) {
          o << format("%-10s%s\n") % labels[i] % data_row_string;
        } else {
          o << format("%s\n") % data_row_string;
        }
      }
      o << endl;
    }
    return o;
  }
  ostream& SNP::write_table (ostream& o) {
    for (PolyTable::const_site_iterator it=sbegin(); it!=send(); ++it) {
      double pos = it->first;
      string str = it->second;
      string str_new = str.substr(str.length()-1, 1);
      for(uint32_t k = 0; k < str.length() - 1; k ++)
        str_new += "\t" + str.substr(k, 1);
      o << format("%.0f\t%s") % pos % str_new << endl;
    }
    return o;
  }
  ostream& SNP::write_structure (ostream& o) {
    vector<double> poss = this->GetPositions();
    vector<string> data = this->GetData();
    string row1;
    for(unsigned j = 0; j < npos; j++) 
      row1 += boost::lexical_cast<string>(poss[j]) + "\t";
    o << format("%s\n") % row1;
    for(unsigned i = 0; i < nind; i ++) {
      string row = str(format("%s\t%s\t") % labels[i] % (i+1)) ;
      for(unsigned j = 0; j < npos; j ++) {
        switch (char(std::toupper(data[i][j]))) {
          case 'A':
            row += "1\t";
            break;
          case 'T':
            row += "2\t";
            break;
          case 'C':
            row += "3\t";
            break;
          case 'G':
            row += "4\t";
            break;
          case 'N':
            row += "9\t";
            break;
          default:
            row += "9\t";
        }
      }
      o << format("%s\n%s") % row % row << endl;
    }
    return o;
  }
/*  ostream& SNP::write_ped (ostream& o, const string& chr) {
    vector<double> poss = this->GetPositions();
    vector<string> data = this->GetData();
    for(uint32_t j = 0; j < npos; j++)
      o << format("%s\t%d\t0\t%d\n") % chr % (j+1) % poss[j];
    for(unsigned i = 0; i < nind; i++) {
      string id = labels[i];
      string row = id + "\t" + boost::lexical_cast<string>(i+1) + "\t0\t0\t0\t0\t";
      for(unsigned j = 0; j < npos; j ++) {
        string nt = boost::lexical_cast<string>(data[i][j]);
        row += nt + " " + nt + "\t";
      }
      o << row << endl;
    }
    return o;
  }*/
  ostream& SNP::write (ostream& o, const string& ofmt) {
    if(ofmt == "snp") {
      this->print(o);
    } else if(ofmt == "phylip") {
      this->write_phylip(o);
    } else if(ofmt == "table") {
      this->write_table(o);
    } else if(ofmt == "structure") {
      this->write_structure(o);
    } else {
      cerr << format("unsupported format: %s\n") % ofmt;
      exit(1);
    }
    return o;
  }

  void SNP::ApplyFreqFilter(unsigned mincount,bool haveOutgroup, unsigned outgroup) {
    this->ApplyFreqFilter(mincount, haveOutgroup, outgroup);
    npos = this->numsites();
  }
  void SNP::RemoveMultiHits(bool skipOutgroup, unsigned outgroup) {
    this->RemoveMultiHits(skipOutgroup, outgroup);
    npos = this->numsites();
  }
  void SNP::RemoveMissing(bool skipOutgroup, unsigned outgroup) {
    this->RemoveMissing(skipOutgroup, outgroup);
    npos = this->numsites();
  }
  void SNP::RemoveAmbiguous(bool skipOutgroup, unsigned outgroup) {
    this->RemoveAmbiguous(skipOutgroup, outgroup);
    npos = this->numsites();
  }

  void SNP::FilterMissing(const int& co_mis, bool skipOutgroup) {
    vector<double> poss = this->GetPositions();
    vector<string> data = this->GetData();
    vector<double> possN;
    vector<string> dataN(nind);
    for(unsigned j = 0; j < npos; ++ j) {
      stateCounter Counts;
      for(unsigned i = 0; i < nind; ++ i) {
        if(!skipOutgroup || !haveOG || (skipOutgroup && haveOG && i != outgroup_index))
          Counts( data[i][j] );
      }
      if( (int)Counts.n <= co_mis ) {
        possN.push_back(poss[j]);
        for(unsigned i = 0; i < nind; ++ i) 
          dataN[i] += data[i][j];
      }
    }
    assign(&possN[0], possN.size(), &dataN[0], dataN.size());
    npos = this->numsites();
  }
  void SNP::RemoveMono(bool skipOutgroup) {
    vector<double> poss = this->GetPositions();
    vector<string> data = this->GetData();
    vector<double> possN;
    vector<string> dataN(data.size());
    for(unsigned j = 0; j < npos; ++ j) {
      stateCounter Counts;
      for(unsigned i = 0; i < nind; ++ i) {
        if(!skipOutgroup || !haveOG || (skipOutgroup && haveOG && i != outgroup_index))
          Counts( data[i][j] );
      }
      if( Counts.nStates() > 1 ) {
        possN.push_back(poss[j]);
        for(unsigned i = 0; i < nind; ++ i) 
          dataN[i] += data[i][j];
      }
    }
    assign(&possN[0], possN.size(), &dataN[0], dataN.size());
    npos = this->numsites();
  }

  void SNP::subset(vector<uint32_t>& idxs_ind, vector<uint32_t>& idxs_pos) {
    vector<double> poss = this->GetPositions();
    vector<string> data = this->GetData();
    
    if(idxs_ind.empty()) push_back(idxs_ind).repeat_fun(data.size(), next<int>(0));
    if(idxs_pos.empty()) push_back(idxs_pos).repeat_fun(poss.size(), next<int>(0));
    nind = idxs_ind.size();
    npos = idxs_pos.size();
    
    vector<double> possN(npos);
    for(unsigned j = 0; j < npos; ++ j)
      possN[j] = poss[idxs_pos[j]];
    
    vector<string> labelsN(nind);
    vector<string> dataN(nind);
    for(unsigned i = 0; i < nind; ++ i) {
      labelsN[i] = labels[idxs_ind[i]];
      for(unsigned j = 0; j < npos; ++ j) 
        dataN[i] += data[idxs_ind[i]][idxs_pos[j]];
    }

    assign(&possN[0], possN.size(), &dataN[0], dataN.size());
    labels = labelsN;
    nind = this->size();
    npos = this->numsites();
  }

/*  void SNP::write2(const string& dir, const string& seqid) {
    std::vector<double> poss = PolyTable::GetPositions();
    std::vector<std::string> data = PolyTable::GetData();
    uint32_t n_ind = data.size(), n_pos = poss.size();

    for(unsigned i = 0; i < n_ind; i++) {
      string id = label(i+1);
      string fo = dirO + "/" + id + "/" + seqid;
      
      ofstream fho(fo.c_str());
      if (!fho.is_open()) {
        cout << format("cannot write output: %s\n") % fo;
        exit(1);
      }
      for (unsigned j = 0; j < n_pos; j ++) {
        fho << format("%d\t%s\n") % int(poss[j]) % data[i][j];
      }
      fho.close();
    }
  }
*/
}



