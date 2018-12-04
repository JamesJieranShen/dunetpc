#ifndef GENVECTORHELPER_H
#define GENVECTORHELPER_H
/// ROOT GenVector Helper functions
/**
 *  Helper functions to convert old ROOT classes to GenVector 
 *  ones compatible with recob::Track. For example, convert 
 *  the TLorentzVector in simb::MCParticle to a point, 
 *  vector, or XYZTVector.
 *
 *  See https://root.cern/doc/v612/Vector.html
 */

//Framework includes

//LArSoft includes
#include "lardataobj/RecoBase/Track.h"

//ROOT includes
#include "Math/Vector4Dfwd.h"

namespace lsu
{
  /**
   * Convert to a recob::Track::Point_t, a Point3D
   */
  template <class T>
  auto toPoint(T vec)
  {
      return recob::Track::Point_t(vec.X(),vec.Y(),vec.Z());
  }
  
  /**
   * Convert to a recob::Track::Vector_t, a Vector3D
   */
  template <class T>
  auto toVector(T vec)
  {
      return recob::Track::Vector_t(vec.X(),vec.Y(),vec.Z());
  }

  /**
   * Convert to a LorentzVector
   */
  template <class T>
  auto toVector4(T vec)
  {
      return ROOT::Math::XYZTVector(vec.X(),vec.Y(),vec.Z(),vec.T());
  }
} // namespace lsu

#endif
