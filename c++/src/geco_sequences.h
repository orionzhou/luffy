#ifndef __GECO__SEQUENCES__
#define __GECO__SEQUENCES__

#include <geco/geco_element.h>
#include <vector>

using namespace geco;

class gLocalTxtSequenceRetriever:public gSequenceRetrieverImplementation{
 private:
  gString fileprefix;
  gString filepostfix;
  virtual gArrayRetrieverImplementation<gChar> * clone() const;
  virtual gSequence getSequence_Internal(const gString & reference,gPos start,gPos end) const;  
 protected:
 public:
  gLocalTxtSequenceRetriever();
  gLocalTxtSequenceRetriever(const char*prefix,const char *postfix=".txt");
  gLocalTxtSequenceRetriever(const gLocalTxtSequenceRetriever &retriever);
  virtual ~gLocalTxtSequenceRetriever();

};

class gLocal2bitSequenceRetriever:public gSequenceRetrieverImplementation{
 private:
  gString i_filename;
  std::vector<gString> i_references;
  gArray<gSize> i_dnaSize;
  gArray<gPos> i_offset;
  std::vector< gArray<gPos> > i_nBlockStarts;
  std::vector< gArray<gPos> > i_nBlockEnds;
  std::vector< gArray<gPos> > i_maskBlockStarts;
  std::vector< gArray<gPos> > i_maskBlockEnds;
  virtual gArrayRetrieverImplementation<gChar> * clone() const;
  virtual gSequence getSequence_Internal(const gString & reference,gPos start,gPos end) const;  
 protected:
 public:
  gLocal2bitSequenceRetriever(const gString & filename,bool useMask);
  gLocal2bitSequenceRetriever(const gLocal2bitSequenceRetriever & retriever);
  virtual ~gLocal2bitSequenceRetriever();  
  gSequence getRandomSequence(gSize length) const;
  gSize getReferenceLength(const gString & reference) const;
};
#endif
