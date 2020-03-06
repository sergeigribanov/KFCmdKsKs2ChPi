#ifndef __KsKs2ChPi_TrPh__
#define __KsKs2ChPi_TrPh__

#include <KFCmd/TrPh.hpp>

class TrPh : public KFCmd::TrPh {
public:
  TrPh(TTree *tree = 0);
  virtual ~TrPh();
  virtual Int_t Cut(Long64_t energy) override final;
  virtual void Loop(const std::string&, double mfield = 1.3) override final;
private:
  bool cutTracks();
  static const double _dZ;
  static const double _dRho;
  static const double _mindEdX;
  static const double _maxdEdX;
  static const double _minTPtot;
  static const double _maxTPtot;
  std::vector<std::size_t> _trackIndices;
};

#endif
