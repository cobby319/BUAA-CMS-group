#ifndef TLORENTZVECTORWITHINDEX_H
#define TLORENTZVECTORWITHINDEX_H

#include "TLorentzVector.h"

/**
 * \brief A \c TLorentzVector with an additional index.
 */
class TLorentzVectorWithIndex : public TLorentzVector
{
  private:
    int kIndex;

  public:
    /// \brief Constructor that takes a \c TLorentzVector and the \c index.
    explicit inline TLorentzVectorWithIndex(const TLorentzVector &vec,
                                            int index);

    /// \brief Constructor that takes cartesian coordinates and the \c index.
    explicit inline TLorentzVectorWithIndex(double x, double y,
                                            double z, double t,
                                            int index);

    /**
     * \brief Constructor that takes a three-vector, the `t`-component and the
     *        \c index.
     */
    explicit inline TLorentzVectorWithIndex(const TVector3 &vec3, double t,
                                            int index);

    /// \brief Retrieves the \c index.
    inline int GetIndex() const { return kIndex; }

    /// \brief Retrieves the \c index.
    inline void SetIndex(int index) { kIndex = index; }

    /**
     * \brief Creates a \c TLorentzVectorWithIndex from its \c pt, \c eta,
     *        \c energy and \c index.
     */
    static inline TLorentzVectorWithIndex PtEtaPhiEIndex(
                                                         double pt, double eta, double phi, double energy, int index);

    /**
     * \brief Creates a \c TLorentzVectorWithIndex from its \c pt, \c eta,
     *        \c mass and \c index.
     */
    static inline TLorentzVectorWithIndex PtEtaPhiMIndex(
                                                         double pt, double eta, double phi, double mass, int index);
};

TLorentzVectorWithIndex::TLorentzVectorWithIndex(const TLorentzVector &vec,
                                                 int index) :
  TLorentzVector(vec),
  kIndex(index)
{}

TLorentzVectorWithIndex::TLorentzVectorWithIndex(double x,
                                                 double y,
                                                 double z,
                                                 double t,
                                                 int index) :
  TLorentzVector(x, y, z, t),
  kIndex(index)
{}

TLorentzVectorWithIndex::TLorentzVectorWithIndex(const TVector3 &vec3,
                                                 double t,
                                                 int index) :
  TLorentzVector(vec3, t),
  kIndex(index)
{}

  TLorentzVectorWithIndex
TLorentzVectorWithIndex::PtEtaPhiEIndex(double pt, double eta,
                                        double phi, double energy,
                                        int index)
{
  TLorentzVector tv;
  tv.SetPtEtaPhiE(pt, eta, phi, energy);
  return TLorentzVectorWithIndex(tv, index);
}

  TLorentzVectorWithIndex
TLorentzVectorWithIndex::PtEtaPhiMIndex(double pt, double eta,
                                        double phi, double mass,
                                        int index)
{
  TLorentzVector tv;
  tv.SetPtEtaPhiM(pt, eta, phi, mass);
  return TLorentzVectorWithIndex(tv, index);
}

#endif // TLORENTZVECTORWITHINDEX_H
