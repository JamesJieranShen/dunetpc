#ifndef TRAJECTORYINTERPEXTRAPALG_H
#define TRAJECTORYINTERPEXTRAPALG_H
/*!
 * Title:   Trjectory Interpolation and Extrapolation Algorithms Class
 * Author:  Justin Hugon (jhugon@fnal.gov)
 *
 * Description: Algorithm for interpolating and extrapolating trajectories
 * Input:       simb::MCTrajectory trajectory, TVector3 point, bool extraploate or not
 * Output:      TVector3 point of closest approach, distance of closest approach, and interplated momentum
*/

//Framework includes
#include "fhiclcpp/ParameterSet.h"

//LArSoft includes
#include "nusimdata/SimulationBase/MCTrajectory.h"

//ROOT includes
#include "TLorentzVector.h"
#include "TVector3.h"

//c++ includes
#include <vector>
#include <utility>

namespace mctrue {
  class TrajectoryInterpExtrapAlg;
}

/// Algorithm for interpolating and extrapolating trajectories
/**
 * This is a class to make matching space hits to trajectories easier.
 * It allows finding the point of closest approach of the trajectory
 * to the space point, and interpolating the trajectory momentum to 
 * that point.
 */
class mctrue::TrajectoryInterpExtrapAlg
{
  public:

    /// Finds the point of closest approach of the trajectory to the point
    /**
     * distance is set to the distance between the point of closest
     * approach and the input point. It will be < 0 if finding the point 
     * of closest approach fails. If extrapolate is false and
     * distance is large, this may mean that the point is nowhere near
     * the trajectory.
     *
     * If extrapolate is true, will try projecting the trajectory
     * start point backwards and end point forwards in search of
     * a point of closer approach.
     */
    inline const TVector3 pointOfClosestApproach(
                const simb::MCTrajectory& trajectory,
                const TVector3& point,
                double& distance,
                bool extrapolate=false) 
        {
            TLorentzVector mom;
            return pointOfClosestApproach(trajectory,point,distance,mom,extrapolate);
        };

    /// Finds the point of closest approach of the trajectory to the point & interpolates momentum
    /**
     * distance is set to the distance between the point of closest
     * approach and the input point. It will be < 0 if finding the point 
     * of closest approach fails. If extrapolate is false and
     * distance is large, this may mean that the point is nowhere near
     * the trajectory.
     *
     * interpolatedMomentum will be set to the approximate four-momentum 
     * at the point of closest approach. If the point of closest approach
     * is between trajectory points, the momentum is computed by linearly
     * interpolating the momentum between them. Otherwise, it is set to
     * the momentum of the nearest trajectory point.
     *
     * If extrapolate is true, will try projecting the trajectory
     * start point backwards and end point forwards in search of
     * a point of closer approach.
     */
    inline const TVector3 pointOfClosestApproach(
                const simb::MCTrajectory& trajectory,
                const TVector3& point,
                double& distance,
                TLorentzVector& interpolatedMomentum,
                bool extrapolate=false)
        {
            size_t i;
            double d;
            return pointOfClosestApproach(trajectory,point,distance,interpolatedMomentum,i,d,extrapolate);
        };

    /// Finds the point of closest approach of the trajectory to the point & interpolates momentum with debug info
    /**
     * distance is set to the distance between the point of closest
     * approach and the input point. It will be < 0 if finding the point 
     * of closest approach fails. If extrapolate is false and
     * distance is large, this may mean that the point is nowhere near
     * the trajectory.
     *
     * interpolatedMomentum will be set to the approximate four-momentum 
     * at the point of closest approach. If the point of closest approach
     * is between trajectory points, the momentum is computed by linearly
     * interpolating the momentum between them. Otherwise, it is set to
     * the momentum of the nearest trajectory point.
     *
     * iClosestTrajPoint will be set to the index of the closest
     * trajectory point to point
     *
     * distanceToClosestTrajPoint will be set to the distance
     * from point to the closest trajectory point
     *
     * If extrapolate is true, will try projecting the trajectory
     * start point backwards and end point forwards in search of
     * a point of closer approach.
     */
    const TVector3 pointOfClosestApproach(
                const simb::MCTrajectory& trajectory,
                const TVector3& point,
                double& distance,
                TLorentzVector& interpolatedMomentum,
                size_t& iClosestTrajPoint,
                double& distanceToClosestTrajPoint,
                bool extrapolate=false);


  private:

    /// Find point of closest approach when closest traj point in middle of traj points
    /**
     * requires that the closest trajectory point has already been found, iClosest
     */
    const TVector3 interplateInMiddle(
                const simb::MCTrajectory& trajectory,
                const TVector3& point,
                double& distance,
                TLorentzVector& interpolatedMomentum,
                const int iClosest);

    /// Find point of closest approach when closest traj point is the first
    /**
     * Assumes the closest trajectory point is the first
     */
    const TVector3 interpExtrapAtBeginning(
                const simb::MCTrajectory& trajectory,
                const TVector3& point,
                double& distance,
                TLorentzVector& interpolatedMomentum,
                bool extrapolate=false);

    /// Find point of closest approach when closest traj point is the last
    /**
     * Assumes the closest trajectory point is the last
     */
    const TVector3 interpExtrapAtEnd(
                const simb::MCTrajectory& trajectory,
                const TVector3& point,
                double& distance,
                TLorentzVector& interpolatedMomentum,
                bool extrapolate=false);
    
    /// Finds the trajectory point that is closest
    /**
     * Returns -1 if none found, and the index of the closest point otherwise
     */
    int findClosestTrajPoint(
                const simb::MCTrajectory& trajectory,
                const TVector3& point,
                double& closestDistance);

    /// Finds the point of closest approach to the line going through the segment
    /**
     * The resulting vector is extrapolated out along the line
     * going through the points
     *
     * The values of locationMeasure are set such that:
     *    location Measure < 0: segment1 is the closest point on the segment
     *    location Measure > 1: segment2 is the closest point on the segment
     *    otherwise: the returned point is the closest point on the segment
     */
    const TVector3 segmentPointOfClosestApproach(
                const TVector3& segment1,
                const TVector3& segment2,
                const TVector3& point,
                double& locationMeasure);

    /// Find the momentum at the point of closest approach
    /**
     * distance is set to the distance between the point of closest
     * approach and the input point. It will be < 0 if finding the point 
     * of closest approach fails. If extrapolate is false and
     * distance is large, this may mean that the point is nowhere near
     * the trajectory.
     */
    const TLorentzVector interpolateMomentum(
                const TVector3& pointOfClosestApproach,
                const TVector3& point1,
                const TVector3& point2,
                const TLorentzVector& momentum1,
                const TLorentzVector& momentum2);
    
};

#endif
