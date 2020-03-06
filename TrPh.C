#include <cmath>
#include <iostream>
#include <algorithm>
#include <limits>
#include <set>

#include <TCanvas.h>
#include <TH2.h>
#include <TStyle.h>

#include <KFCmd/Hypo2ChPionsKsKs.hpp>

#include "TrPh.h"

const double TrPh::_dZ = 30;
const double TrPh::_dRho = 30;
const double TrPh::_mindEdX = 0;
const double TrPh::_maxdEdX = 15000;
const double TrPh::_minTPtot = 5;
const double TrPh::_maxTPtot = 1000;

TrPh::TrPh(TTree *tree) :
  KFCmd::TrPh(tree) {
}

TrPh::~TrPh() {
}

bool TrPh::cutTracks() {
  _trackIndices.clear();
  for (int i = 0; i < nt; i++) {
    bool point = (std::fabs(tz[i]) < _dZ) && (std::fabs(trho[i]) < _dRho);
    bool dedx = (tdedx[i] > _mindEdX) && (tdedx[i] < _maxdEdX);
    bool ptot = (tptot[i] > _minTPtot) && (tptot[i] < _maxTPtot);
    if (point && dedx && ptot) _trackIndices.push_back(i);
  }
  if (_trackIndices.size() == 6) {
    int totalCharge = 0;
    for (int i = 0; i < 6; ++i) totalCharge += tcharge[_trackIndices[i]];
    return (totalCharge == 0);
  }
  return false;
}

Int_t TrPh::Cut(Long64_t) {
  if (nt < 6) return -1;
  if (!cutTracks()) return -1;
  // if (!cutPhotons()) return -1;
  std::vector<Int_t> charges(nt);
  std::copy(tcharge, tcharge + nt, charges.begin());
  std::sort(_trackIndices.begin(), _trackIndices.end(),
            [&charges](int i, int j) { return charges[i] < charges[j]; });
  return 1;
}

void TrPh::Loop(const std::string& outpath, double mfield) {
  if (fChain == 0) return;
  fChain->GetEntry(0);
  auto outfl = TFile::Open(outpath.c_str(), "recreate");
  TH1F h_kf_mks1("h_kf_mks1", "", 512, 0, 2000);
  TH1F h_in_mks1("h_in_mks1", "", 512, 0, 2000);
  TH1F h_kf_mks2("h_kf_mks2", "", 512, 0, 2000);
  TH1F h_in_mks2("h_in_mks2", "", 512, 0, 2000);
  TH1F h_kf_chi2("h_kf_chi2", "", 1024, 0, 1024);
  TH1F h_vtx1_dr("h_vtx1_dr", "", 600, 0, 60);
  TH1F h_vtx2_dr("h_vtx2_dr", "", 600, 0, 60);
  std::set<std::string> sKs1 = {"pi+_1", "pi-_1"};
  std::set<std::string> sKs2 = {"pi+_2", "pi-_2"};
  KFCmd::Hypo2ChPionsKsKs hypo(2 * emeas, mfield);
  hypo.setBeamXY(xbeam, ybeam);
  hypo.fixVertexComponent("vtx0", xbeam, KFBase::VERTEX_X);
  hypo.fixVertexComponent("vtx0", ybeam, KFBase::VERTEX_Y);
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry = 0; jentry < nentries; jentry++) {
    if (jentry % 100 == 0)
      std::cout << jentry << std::endl;
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;
    if (Cut(ientry) < 0) continue;
    hypo.setBeamXY(xbeam, ybeam);
    hypo.fixVertexComponent("vtx0", xbeam, KFBase::VERTEX_X);
    hypo.fixVertexComponent("vtx0", ybeam, KFBase::VERTEX_Y);
    std::vector<int> mi_perm = {0, 1, 2};
    double kf_chi2 = std::numeric_limits<double>::infinity();
    double tchi2;
    bool flag = false;
    double v_kf_mks1 = 0;
    double v_kf_mks2 = 0;
    double v_in_mks1 = 0;
    double v_in_mks2 = 0;
    double v_vtx1_dr = 0;
    double v_vtx2_dr = 0;
    do {
      std::vector<int> pl_perm = {3, 4, 5};
      do {
	if (!hypo.fillTrack("pi-_0", _trackIndices[mi_perm[0]], *this)) continue;
	if (!hypo.fillTrack("pi-_1", _trackIndices[mi_perm[1]], *this)) continue;
	if (!hypo.fillTrack("pi-_2", _trackIndices[mi_perm[2]], *this)) continue;
	if (!hypo.fillTrack("pi+_0", _trackIndices[pl_perm[0]], *this)) continue;
	if (!hypo.fillTrack("pi+_1", _trackIndices[pl_perm[1]], *this)) continue;
	if (!hypo.fillTrack("pi+_2", _trackIndices[pl_perm[2]], *this)) continue;
	hypo.optimize();
	if (hypo.getErrorCode() != 0) continue;
	tchi2 = hypo.getChiSquare();
	if (tchi2 < kf_chi2) {
	  flag = true;
	  kf_chi2 = tchi2;
	  v_in_mks1 = hypo.getInitialMomentum(sKs1).M();
	  v_in_mks2 = hypo.getInitialMomentum(sKs2).M();
	  v_kf_mks1 = hypo.getFinalMomentum(sKs1).M();
	  v_kf_mks2 = hypo.getFinalMomentum(sKs2).M();
	  v_vtx1_dr = (hypo.getFinalVertex("vtx1") -  hypo.getFinalVertex("vtx0")).Mag();
	  v_vtx2_dr = (hypo.getFinalVertex("vtx2") -  hypo.getFinalVertex("vtx0")).Mag();
	}
      } while (std::next_permutation(pl_perm.begin(), pl_perm.end()));
    } while (std::next_permutation(mi_perm.begin(), mi_perm.end()));
    if (flag && kf_chi2 < 200) {
      h_kf_mks1.Fill(v_kf_mks1);
      h_kf_mks2.Fill(v_kf_mks2);
      h_in_mks1.Fill(v_in_mks1);
      h_in_mks2.Fill(v_in_mks2);
      h_kf_chi2.Fill(kf_chi2);
      h_vtx1_dr.Fill(v_vtx1_dr);
      h_vtx2_dr.Fill(v_vtx2_dr);
    }
  }
  outfl->cd();
  h_kf_mks1.Write();
  h_kf_mks2.Write();
  h_in_mks1.Write();
  h_in_mks2.Write();
  h_kf_chi2.Write();
  h_vtx1_dr.Write();
  h_vtx2_dr.Write();
  outfl->Close();
}
