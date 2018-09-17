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
     * planePoint is any point on the plane
     * planeNormal is the normal vector of the plane
     * linePoint is any point on the line
     * lineDirection is the direction of the line
     *
     * returns a TVector3 with the point. If the point 
     * doesn't go through the plane, returns a vector 
     * with -9999999 for all three components. This 
     * means the vector is parallel to the plane.
     */
    inline const TVector3 linePlane(
                const TVector3& planePoint,
                const TVector3& planeNormal,
                const TVector3& linePoint,
                const TVector3& lineDirection)
        {
            const double normalDotDir = planeNormal.Dot(lineDirection);
            if (normalDotDir == 0.) return TVector3(-9999999.,-9999999.,-9999999.);
            const double d = (planePoint-linePoint)*planeNormal / (normalDotDir);
            return (d*lineDirection) + linePoint;
        };

    /// Finds the point on a constant x-plane intersected by a line
    /**
     * planeX is the x-coordinate of the plane
     * linePoint is any point on the line
     * lineDirection is the direction of the line
     *
     * returns a TVector3 with the point. If the point 
     * doesn't go through the plane, returns a vector 
     * with -9999999 for all three components. This 
     * means the vector is parallel to the plane.
     */
    inline const TVector3 lineXPlane(
                const double planeX,
                const TVector3& linePoint,
                const TVector3& lineDirection)
        {
            const TVector3 planePoint(planeX,0.,0.);
            const TVector3 planeNormal(1.,0.,0.);
            return linePlane(planePoint,planeNormal,linePoint,lineDirection);
        };
    
    /// Finds the point on a constant y-plane intersected by a line
    /**
     * planeY is the y-coordinate of the plane
     * linePoint is any point on the line
     * lineDirection is the direction of the line
     *
     * returns a TVector3 with the point. If the point 
     * doesn't go through the plane, returns a vector 
     * with -9999999 for all three components. This 
     * means the vector is parallel to the plane.
     */
    inline const TVector3 lineYPlane(
                const double planeY,
                const TVector3& linePoint,
                const TVector3& lineDirection)
        {
            const TVector3 planePoint(0.,planeY,0.);
            const TVector3 planeNormal(0.,1.,0.);
            return linePlane(planePoint,planeNormal,linePoint,lineDirection);
        };
    
    /// Finds the point on a constant z-plane intersected by a line
    /**
     * planeZ is the z-coordinate of the plane
     * linePoint is any point on the line
     * lineDirection is the direction of the line
     *
     * returns a TVector3 with the point. If the point 
     * doesn't go through the plane, returns a vector 
     * with -9999999 for all three components. This 
     * means the vector is parallel to the plane.
     */
    inline const TVector3 lineZPlane(
                const double planeZ,
                const TVector3& linePoint,
                const TVector3& lineDirection)
        {
            const TVector3 planePoint(0.,0.,planeZ);
            const TVector3 planeNormal(0.,0.,1.);
            return linePlane(planePoint,planeNormal,linePoint,lineDirection);
        };

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
        };


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
        };

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
        };

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
        };
    
};

#endif
