#ifndef KFPARTICLESPHENIX_KFPARTICLECENTRALITYINFO_H
#define KFPARTICLESPHENIX_KFPARTICLECENTRALITYINFO_H

class PHCompositeNode;
class TTree;

class KFParticle_centralityInfo
{
 public:
  KFParticle_centralityInfo() = default;
  virtual ~KFParticle_centralityInfo() = default;

 protected:
  //! Returns false if there is no CentralityInfo node. In that case no branches
  //! are made and the caller skips the fill and reset calls, so running on a
  //! sample without centrality costs nothing and prints nothing.
  bool buildCentralityBranches(PHCompositeNode *topNode, TTree *m_tree);

  //! crossingIsZero must be false when any daughter track comes from a bunch
  //! crossing other than the triggered one. The MBD is only read out at the
  //! triggered crossing, so its centrality describes a different collision and
  //! must not be assigned to such a candidate.
  void fillCentralityBranches(PHCompositeNode *topNode, bool crossingIsZero);

  void resetCentralityBranches();

 private:
  int m_centrality_mbd_bin{-999};
  bool m_is_min_bias{false};
};

#endif  // KFPARTICLESPHENIX_KFPARTICLECENTRALITYINFO_H
