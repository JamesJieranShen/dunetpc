#ifndef PLANEINTERSECTIONFINDER_H
#define PLANEINTERSECTIONFINDER_H
/*!
 * Title:   Plane Intersection Finder Functions
 * Author:  Justin Hugon (jhugon@fnal.gov)
 *
 * Description: Algorithms for finding where something intersects a plane
 * Input:       
 * Output:
*/

//Framework includes

//LArSoft includes
//#include "nusimdata/SimulationBase/MCTrajectory.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/RecoBase/Track.h"

//ROOT includes
//#include "TLorentzVector.h"
#include "TVector3.h"

//c++ includes
//#include <vector>
//#include <utility>

/// Plane Intersection Finder
/**
 *
 */
namespace lsu
{
    /// Finds the point on a plane intersected by a line
    /**
     * planePoint is any point on the plane (TLorentzVector, TVector3, and similar)
     * planeNormal is the normal vector of the plane (TLorentzVector, TVector3, and similar)
     * linePoint is any point on the line (TLorentzVector, TVector3, and similar)
     * lineDirection is the direction of the line (TLorentzVector, TVector3, and similar)
     *
     * returns a TVector3 with the point. If the point 
     * doesn't go through the plane, returns a vector 
     * with -9999999 for all three components. This 
     * means the vector is parallel to the plane.
     */
    template <class A, class B, class C, class D>
    inline const TVector3 linePlane(
                const A& planePoint,
                const B& planeNormal,
                const C& linePoint,
                const D& lineDirection)
        {
            const TVector3 pp(planePoint.X(),planePoint.Y(),planePoint.Z());
            TVector3 pn(planeNormal.X(),planeNormal.Y(),planeNormal.Z());
            const TVector3 lp(linePoint.X(),linePoint.Y(),linePoint.Z());
            TVector3 ld(lineDirection.X(),lineDirection.Y(),lineDirection.Z());
            pn = pn.Unit();
            ld = ld.Unit();
            const double normalDotDir = pn.Dot(ld);
            if (normalDotDir == 0.) return TVector3(-9999999.,-9999999.,-9999999.);
            const double d = (pp-lp)*pn / (normalDotDir);
            return (d*ld) + lp;
        }

    /// Finds the point on a constant x-plane intersected by a line
    /**
     * planeX is the x-coordinate of the plane
     * linePoint is any point on the line (TLorentzVector, TVector3, and similar)
     * lineDirection is the direction of the line (TLorentzVector, TVector3, and similar)
     *
     * returns a TVector3 with the point. If the point 
     * doesn't go through the plane, returns a vector 
     * with -9999999 for all three components. This 
     * means the vector is parallel to the plane.
     */
    template <class A, class B>
    inline const TVector3 lineXPlane(
                const double planeX,
                const A& linePoint,
                const B& lineDirection)
        {
            const TVector3 planePoint(planeX,0.,0.);
            const TVector3 planeNormal(1.,0.,0.);
            return linePlane(planePoint,planeNormal,linePoint,lineDirection);
        }
    
    /// Finds the point on a constant y-plane intersected by a line
    /**
     * planeY is the y-coordinate of the plane
     * linePoint is any point on the line (TLorentzVector, TVector3, and similar)
     * lineDirection is the direction of the line (TLorentzVector, TVector3, and similar)
     *
     * returns a TVector3 with the point. If the point 
     * doesn't go through the plane, returns a vector 
     * with -9999999 for all three components. This 
     * means the vector is parallel to the plane.
     */
    template <class A, class B>
    inline const TVector3 lineYPlane(
                const double planeY,
                const A& linePoint,
                const B& lineDirection)
        {
            const TVector3 planePoint(0.,planeY,0.);
            const TVector3 planeNormal(0.,1.,0.);
            return linePlane(planePoint,planeNormal,linePoint,lineDirection);
        }
    
    /// Finds the point on a constant z-plane intersected by a line
    /**
     * planeZ is the z-coordinate of the plane
     * linePoint is any point on the line (TLorentzVector, TVector3, and similar)
     * lineDirection is the direction of the line (TLorentzVector, TVector3, and similar)
     *
     * returns a TVector3 with the point. If the point 
     * doesn't go through the plane, returns a vector 
     * with -9999999 for all three components. This 
     * means the vector is parallel to the plane.
     */
    template <class A, class B>
    inline const TVector3 lineZPlane(
                const double planeZ,
                const A& linePoint,
                const B& lineDirection)
        {
            const TVector3 planePoint(0.,0.,planeZ);
            const TVector3 planeNormal(0.,0.,1.);
            return linePlane(planePoint,planeNormal,linePoint,lineDirection);
        }

    /// Finds the point on a constant z-plane intersected by a track
    /**
     * Similar to lineZPlane, but figures out if the plane is > or < the track
     * and extrapolates from the closer end. If the plane is in the middle
     * of the track, will just return a bogus point.
     *
     * planeZ is the z-coordinate of the plane
     * linePoint is any point on the line (TLorentzVector, TVector3, and similar)
     * lineDirection is the direction of the line (TLorentzVector, TVector3, and similar)
     *
     *
     * returns a TVector3 with the point. If the point 
     * doesn't go through the plane, returns a vector 
     * with -9999999 for all three components. This 
     * means the vector is parallel to the plane.
     */
    inline const TVector3 trackZPlane(
                const double planeZ,
                const recob::Track & track)
        {
            if (track.NumberTrajectoryPoints() == 0) 
            {
                return TVector3(-9999999.,-9999999.,-9999999.);
            }
            const TVector3 planePoint(0.,0.,planeZ);
            const TVector3 planeNormal(0.,0.,1.);
            const auto & traj = track.Trajectory();
            const auto startPoint = traj.Start();
            const auto startDir = traj.StartDirection();
            const auto & endPoint = traj.End();
            const auto & endDir = traj.EndDirection();
            const double startZ = startPoint.Z();
            const double endZ = endPoint.Z();
            if (planeZ < startZ && planeZ < endZ)
            {
              if (startZ < endZ)
              {
                return linePlane(planePoint,planeNormal,startPoint,startDir);
              }
              else
              {
                return linePlane(planePoint,planeNormal,endPoint,endDir);
              }
            }
            else if (planeZ > startZ && planeZ > endZ)
            {
              if (startZ > endZ)
              {
                return linePlane(planePoint,planeNormal,startPoint,startDir);
              }
              else
              {
                return linePlane(planePoint,planeNormal,endPoint,endDir);
              }
            }
            else // plane inbetween points, we don't do that
            {
                return TVector3(-9999999.,-9999999.,-9999999.);
            }
        }

    /// Finds the point on a constant z-plane intersected by the MCParticle Start Direction
    /**
     * Similar to lineZPlane, but takes the MCParticle start point and start direction as the line
     *
     * planeZ is the z-coordinate of the plane
     * particle is the MCParticle
     *
     *
     * returns a TVector3 with the point. If the point 
     * doesn't go through the plane, returns a vector 
     * with -9999999 for all three components. This 
     * means the vector is parallel to the plane.
     */
    inline const TVector3 mcPartStartZPlane(
                const double planeZ,
                const simb::MCParticle & particle)
        {
            return lineZPlane(planeZ,particle.Position(),particle.Momentum());
        }

    /// Finds if the point is within a box
    /**
     * point is the point
     * minX is the min x coordinate, similar for minY, minZ
     * maxX is the max x coordinate, similar for maxY, maxZ
     *
     */
    inline const bool pointInXYBox(
                const TVector3& point,
                const double minX,
                const double maxX,
                const double minY,
                const double maxY,
                const double minZ,
                const double maxZ,
                const bool includeBoundaries=false
              )
        {
            if (includeBoundaries)
            {
              return point.X() >= minX && point.X() <= maxX 
                && point.Y() >= minY && point.Y() <= maxY
                && point.Z() >= minZ && point.Z() <= maxZ;
            }
            else
            {
              return point.X() > minX && point.X() < maxX 
                && point.Y() > minY && point.Y() < maxY
                && point.Z() > minZ && point.Z() < maxZ;
            }
        }


    /// Finds if the x-y coordinates of a point are within an x-y box
    /**
     * point is the point
     * minX is the min x coordinate, similar for minY
     * maxX is the max x coordinate, similar for maxY
     *
     */
    inline const bool pointInXYBox(
                const TVector3& point,
                const double minX,
                const double maxX,
                const double minY,
                const double maxY
              )
        {
            return point.X() > minX && point.X() < maxX && point.Y() > minY && point.Y() < maxY;
        }

    /// Finds if the x-z coordinates of a point are within an x-z box
    /**
     * point is the point
     * minX is the min x coordinate, similar for minZ
     * maxX is the max x coordinate, similar for maxZ
     *
     */
    inline const bool pointInXZBox(
                const TVector3& point,
                const double minX,
                const double maxX,
                const double minZ,
                const double maxZ
              )
        {
            return point.X() > minX && point.X() < maxX && point.Z() > minZ && point.Z() < maxZ;
        }

    /// Finds if the y-z coordinates of a point are within an y-z box
    /**
     * point is the point
     * minY is the min y coordinate, similar for minZ
     * maxY is the max y coordinate, similar for maxZ
     *
     */
    inline const bool pointInYZBox(
                const TVector3& point,
                const double minY,
                const double maxY,
                const double minZ,
                const double maxZ
              )
        {
            return point.Y() > minY && point.Y() < maxY && point.Z() > minZ && point.Z() < maxZ;
        }
    
}

#endif
