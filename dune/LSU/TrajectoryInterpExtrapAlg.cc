/*!
 * Title:   Trjectory Interpolation and Extrapolation Algorithms Class
 * Author:  Justin Hugon (jhugon@fnal.gov)
 *
 * Description: Algorithm for interpolating and extrapolating trajectories
 * Input:       simb::MCTrajectory trajectory, TVector3 point, bool extraploate or not
 * Output:      TVector3 point of closest approach, distance of closest approach, and interplated momentum
*/

#include "TrajectoryInterpExtrapAlg.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

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
const TVector3 
mctrue::TrajectoryInterpExtrapAlg::pointOfClosestApproach(
            const simb::MCTrajectory& trajectory,
            const TVector3& point,
            double& distance,
            TLorentzVector& interpolatedMomentum,
            size_t& iClosestTrajPoint,
            double& distanceToClosestTrajPoint,
            bool extrapolate)
{
  const size_t ntrajectory = trajectory.size();
  const size_t iLastTrajectory = ntrajectory-1;
  if (ntrajectory == 0)
  {
    distance = -1.;
    iClosestTrajPoint=0;
    return TVector3();
  }
  else if (ntrajectory == 1)
  {
    distance = (point-trajectory.Position(0).Vect()).Mag();
    interpolatedMomentum = trajectory.Momentum(0);
    iClosestTrajPoint=0;
    return trajectory.Position(0).Vect();
  }
  const int iClosest = findClosestTrajPoint(trajectory,point,distanceToClosestTrajPoint);
  if (iClosest < 0)
  {
    distance = -1.;
    iClosestTrajPoint=0;
    return TVector3();
  }
  iClosestTrajPoint = iClosest;
  if (iClosestTrajPoint == 0)
  {
    return interpExtrapAtBeginning(trajectory,point,distance,interpolatedMomentum,extrapolate);
  }
  else if (iClosestTrajPoint == iLastTrajectory)
  {
    return interpExtrapAtEnd(trajectory,point,distance,interpolatedMomentum,extrapolate);
  }
  else
  {
    return interplateInMiddle(trajectory,point,distance,interpolatedMomentum,iClosestTrajPoint);
  }
}

const TVector3 
mctrue::TrajectoryInterpExtrapAlg::pointClosestToPlane(
                const simb::MCTrajectory& trajectory,
                const double& planeZ,
                TLorentzVector& interpolatedMomentum,
                size_t& iClosestTrajPoint,
                double& distanceToClosestTrajPoint)
{
  const size_t ntrajectory = trajectory.size();
  if (ntrajectory == 0)
  {
    distanceToClosestTrajPoint = -1.;
    iClosestTrajPoint=0;
    return TVector3();
  }
  else if (ntrajectory == 1)
  {
    distanceToClosestTrajPoint = fabs(trajectory.Z(0)-planeZ);
    interpolatedMomentum = trajectory.Momentum(0);
    iClosestTrajPoint=0;
    return trajectory.Position(0).Vect();
  }
  int iClosest = iPointClosestToPlane(trajectory,planeZ);
  if (iClosest < 0)
  {
    distanceToClosestTrajPoint = -1.;
    iClosestTrajPoint=0;
    return TVector3();
  }
  iClosestTrajPoint = iClosest;
  const auto& closestPoint = trajectory.Position(iClosest).Vect();
  interpolatedMomentum = trajectory.Momentum(iClosest);
  distanceToClosestTrajPoint = fabs(closestPoint.Z()-planeZ);
  return closestPoint;
}

/// Find point of closest approach when closest traj point in middle of traj points
/**
 * requires that the closest trajectory point has already been found, iClosest
 */
const TVector3 
mctrue::TrajectoryInterpExtrapAlg::interplateInMiddle(
            const simb::MCTrajectory& trajectory,
            const TVector3& point,
            double& distance,
            TLorentzVector& interpolatedMomentum,
            const int iClosest)
{
  double locationMeasureBefore = 0.;
  double locationMeasureAfter = 0.;
  TVector3 segmentBefore = trajectory.Position(iClosest-1).Vect();
  TVector3 segmentCenter = trajectory.Position(iClosest).Vect();
  TVector3 segmentAfter = trajectory.Position(iClosest+1).Vect();
  TVector3 closestBefore = segmentPointOfClosestApproach(segmentBefore,segmentCenter,point,locationMeasureBefore);
  TVector3 closestAfter = segmentPointOfClosestApproach(segmentAfter,segmentCenter,point,locationMeasureAfter);
  if (locationMeasureBefore < 0)
  {
    throw cet::exception("IvalidValue","locationMeasureBefore < 0 when segmentCenter should be closer");
  }
  if (locationMeasureBefore > 1)
  {
      closestBefore = segmentCenter;
  }
  if (locationMeasureAfter < 0)
  {
    throw cet::exception("IvalidValue","locationMeasureAfter < 0 when segmentCenter should be closer");
  }
  if (locationMeasureAfter > 1)
  {
    closestAfter = segmentCenter;
  }
  double distanceBefore = (point-closestBefore).Mag();
  double distanceAfter = (point-closestAfter).Mag();
  if (distanceBefore < distanceAfter)
  {
    if(locationMeasureBefore > 1)
    {
      interpolatedMomentum = trajectory.Momentum(iClosest);
    }
    else
    {
      interpolatedMomentum = interpolateMomentum(
                              closestBefore,
                              trajectory.Position(iClosest-1).Vect(),
                              trajectory.Position(iClosest).Vect(),
                              trajectory.Momentum(iClosest-1),
                              trajectory.Momentum(iClosest));
    }
    distance = distanceBefore;
    return closestBefore;
  }
  else
  {
    if(locationMeasureAfter > 1)
    {
      interpolatedMomentum = trajectory.Momentum(iClosest);
    }
    else
    {
      interpolatedMomentum = interpolateMomentum(
                              closestAfter,
                              trajectory.Position(iClosest).Vect(),
                              trajectory.Position(iClosest+1).Vect(),
                              trajectory.Momentum(iClosest),
                              trajectory.Momentum(iClosest+1));
    }
    distance = distanceAfter;
    return closestAfter;
  }
}

/// Find point of closest approach when closest traj point is the first
/**
 * Assumes the closest trajectory point is the first
 */
const TVector3 
mctrue::TrajectoryInterpExtrapAlg::interpExtrapAtBeginning(
            const simb::MCTrajectory& trajectory,
            const TVector3& point,
            double& distance,
            TLorentzVector& interpolatedMomentum,
            bool extrapolate)
{
  double locationMeasure = 0.;
  TVector3 segment1 = trajectory.Position(0).Vect();
  TVector3 segment2 = trajectory.Position(1).Vect();
  TVector3 thisPoint = segmentPointOfClosestApproach(segment1,segment2,point,locationMeasure);
  if ((!extrapolate) && locationMeasure < 0.)
  {
    distance = (point-trajectory.Position(0).Vect()).Mag();
    interpolatedMomentum = trajectory.Momentum(0);
    return trajectory.Position(0).Vect();
  }
  else if (locationMeasure > 1.)
  {
    mf::LogError lerr("TrajectoryInterpExtrapAlg vars: \n");
    lerr << " locationMeasure: " << locationMeasure << "\n";
    lerr << " segment1: (";
    lerr << segment1.X() << ",";
    lerr << segment1.Y() << ",";
    lerr << segment1.Z() << ")\n";
    lerr << " segment2: (";
    lerr << segment2.X() << ",";
    lerr << segment2.Y() << ",";
    lerr << segment2.Z() << ")\n";
    lerr << " point: (";
    lerr << point.X() << ",";
    lerr << point.Y() << ",";
    lerr << point.Z() << ")\n";
    lerr << " thisPoint: (";
    lerr << thisPoint.X() << ",";
    lerr << thisPoint.Y() << ",";
    lerr << thisPoint.Z() << ")\n";
    lerr << " point-segment1 distance: " << (point-segment1).Mag() << '\n';
    lerr << " point-segment2 distance: " << (point-segment2).Mag() << '\n';
    lerr << " segment2-segment1 distance: " << (segment2-segment1).Mag() << '\n';
    std::string errMsg;
    errMsg += "TrajectoryInterpExtrapAlg::interpExtrapAtBeginning:";
    errMsg += " locationMeasure > 1 when segment1 should be closer.";
    throw cet::exception("IvalidValue",errMsg);
  }
  else
  {
    distance = (point-thisPoint).Mag();
    if (extrapolate && locationMeasure < 0.)
    {
      interpolatedMomentum = trajectory.Momentum(0);
    }
    else
    {
      interpolatedMomentum = interpolateMomentum(
                                thisPoint,
                                trajectory.Position(0).Vect(),
                                trajectory.Position(1).Vect(),
                                trajectory.Momentum(0),
                                trajectory.Momentum(1));
    }
    return thisPoint;
  }
}

/// Find point of closest approach when closest traj point is the last
/**
 * Assumes the closest trajectory point is the last
 */
const TVector3 
mctrue::TrajectoryInterpExtrapAlg::interpExtrapAtEnd(
            const simb::MCTrajectory& trajectory,
            const TVector3& point,
            double& distance,
            TLorentzVector& interpolatedMomentum,
            bool extrapolate)
{
  double locationMeasure = 0.;
  const size_t ntrajectory = trajectory.size();
  const int iLastTrajectory = ntrajectory-1;
  TVector3 segment1 = trajectory.Position(iLastTrajectory).Vect();
  TVector3 segment2 = trajectory.Position(iLastTrajectory-1).Vect();
  TVector3 thisPoint = segmentPointOfClosestApproach(segment1,segment2,point,locationMeasure);
  if ((!extrapolate) && locationMeasure < 0.)
  {
    distance = (point-trajectory.Position(iLastTrajectory).Vect()).Mag();
    interpolatedMomentum = trajectory.Momentum(iLastTrajectory);
    return trajectory.Position(iLastTrajectory).Vect();
  }
  else if (locationMeasure > 1.)
  {
    std::string errMsg;
    errMsg += "TrajectoryInterpExtrapAlg::interpExtrapAtEnd:";
    errMsg += " locationMeasure > 1 when segment1 should be closer.";
    errMsg += " locationMeasure: " + std::to_string(locationMeasure);
    errMsg += " segment1: (";
    errMsg += std::to_string(segment1.X()) + ",";
    errMsg += std::to_string(segment1.Y()) + ",";
    errMsg += std::to_string(segment1.Z()) + ")";
    errMsg += " segment2: (";
    errMsg += std::to_string(segment2.X()) + ",";
    errMsg += std::to_string(segment2.Y()) + ",";
    errMsg += std::to_string(segment2.Z()) + ")";
    errMsg += " point: (";
    errMsg += std::to_string(point.X()) + ",";
    errMsg += std::to_string(point.Y()) + ",";
    errMsg += std::to_string(point.Z()) + ")";
    errMsg += " thisPoint: (";
    errMsg += std::to_string(thisPoint.X()) + ",";
    errMsg += std::to_string(thisPoint.Y()) + ",";
    errMsg += std::to_string(thisPoint.Z()) + ")";
    throw cet::exception("IvalidValue",errMsg);
  }
  else
  {
    distance = (point-thisPoint).Mag();
    if (extrapolate && locationMeasure < 0.)
    {
      interpolatedMomentum = trajectory.Momentum(iLastTrajectory);
    }
    else
    {
      interpolatedMomentum = interpolateMomentum(
                              thisPoint,
                              trajectory.Position(iLastTrajectory).Vect(),
                              trajectory.Position(iLastTrajectory-1).Vect(),
                              trajectory.Momentum(iLastTrajectory),
                              trajectory.Momentum(iLastTrajectory-1));
    }
    return thisPoint;
  }
}

/// Finds the trajectory point that is closest
/**
 * Returns -1 if none found, and the index of the closest point otherwise
 */
int
mctrue::TrajectoryInterpExtrapAlg::findClosestTrajPoint(
            const simb::MCTrajectory& trajectory,
            const TVector3& point,
            double& closestDistance)
{
  const size_t ntrajectory = trajectory.size();
  int result = -1;
  closestDistance = 1e15;
  for(size_t iPoint=0; iPoint<ntrajectory; iPoint++)
  {
    double thisDistance = (point-trajectory.Position(iPoint).Vect()).Mag();
    if (thisDistance < closestDistance)
    {
      result = iPoint;
      closestDistance = thisDistance;
    }
  }
  return result;
}

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
const TVector3 
mctrue::TrajectoryInterpExtrapAlg::segmentPointOfClosestApproach(
            const TVector3& segment1,
            const TVector3& segment2,
            const TVector3& point,
            double& locationMeasure)
{
   const TVector3 u = point-segment1;
   const TVector3 v = segment2-segment1;
   const double vMagInv = 1./v.Mag();
   locationMeasure = u.Dot(v)*vMagInv*vMagInv;
   const TVector3 vhat = v * vMagInv;
   const TVector3 result = segment1 + u.Dot(vhat)*vhat;
   //mf::LogError ldb("TrajectoryInterpExtrapAlg::segmentPointOfClosestApproach");
   //ldb << " segment1: (";
   //ldb << segment1.X() << ",";
   //ldb << segment1.Y() << ",";
   //ldb << segment1.Z() << ")\n";
   //ldb << " segment2: (";
   //ldb << segment2.X() << ",";
   //ldb << segment2.Y() << ",";
   //ldb << segment2.Z() << ")\n";
   //ldb << " point: (";
   //ldb << point.X() << ",";
   //ldb << point.Y() << ",";
   //ldb << point.Z() << ")\n";
   //ldb << " u: (";
   //ldb << u.X() << ",";
   //ldb << u.Y() << ",";
   //ldb << u.Z() << ")\n";
   //ldb << " v: (";
   //ldb << v.X() << ",";
   //ldb << v.Y() << ",";
   //ldb << v.Z() << ")\n";
   //ldb << "u dot v: "<< u.Dot(v) << "\n";
   //ldb << "|v|: "<< v.Mag() << "\n";
   //ldb << "|u|: "<< u.Mag() << "\n";
   //ldb << "(u dot v ) / |u|: "<< u.Dot(v)/u.Mag() << "\n";
   //ldb << "(u dot v ) / |v|: "<< u.Dot(v)/v.Mag()/v.Mag() << "\n";
   //ldb << "locationMeasure: "<< locationMeasure << "\n";
   //ldb << " result: (";
   //ldb << result.X() << ",";
   //ldb << result.Y() << ",";
   //ldb << result.Z() << ")\n";
   return result;
}

/// Find the momentum at the point of closest approach
/**
 * distance is set to the distance between the point of closest
 * approach and the input point. If extrapolate is false and
 * distance is large, this may mean that the point is nowhere near
 * the trajectory.
 */
const TLorentzVector 
mctrue::TrajectoryInterpExtrapAlg::interpolateMomentum(
            const TVector3& pointOfClosestApproach,
            const TVector3& point1,
            const TVector3& point2,
            const TLorentzVector& momentum1,
            const TLorentzVector& momentum2)
{
    //mf::LogError ldb("interpolateMomentum vars: ");
    //ldb << "pointOfClosestApproach " << pointOfClosestApproach.X() 
    //    << ", " << pointOfClosestApproach.Y() << ", " << pointOfClosestApproach.Z() << " ";
    //ldb << "point1 " << point1.X() 
    //    << ", " << point1.Y() << ", " << point1.Z() << " ";
    //ldb << "point2 " << point2.X() 
    //    << ", " << point2.Y() << ", " << point2.Z() << " ";
    //ldb << "momentum1 " << momentum1.X() 
    //    << ", " << momentum1.Y() << ", " << momentum1.Z() << ", " << momentum1.T() << " ";
    //ldb << "momentum2 " << momentum2.X() 
    //    << ", " << momentum2.Y() << ", " << momentum2.Z() << ", " << momentum2.T() << " ";
    if (pointOfClosestApproach == point1)
    {
      return momentum1;
    }
    else if (pointOfClosestApproach == point2)
    {
      return momentum2;
    }
    const double pointDistance = (point2-point1).Mag();
    const double approachDistance = (pointOfClosestApproach-point1).Mag();
    const double distanceFraction = approachDistance/pointDistance;
    const TLorentzVector momDiff = momentum2-momentum1;
    const TLorentzVector result = distanceFraction*momDiff + momentum1;
    //ldb << '\n';
    //ldb << "pointDistance " << pointDistance << " ";
    //ldb << "approachDistance " << approachDistance << " ";
    //ldb << "distanceFraction " << distanceFraction << " ";
    //ldb << "momDiff " << momDiff.X() 
    //    << ", " << momDiff.Y() << ", " << momDiff.Z() << ", " << momDiff.T() << " ";
    //ldb << "result " << result.X() 
    //    << ", " << result.Y() << ", " << result.Z() << ", " << result.T() << " ";
    return result;
}

const int
mctrue::TrajectoryInterpExtrapAlg::iPointClosestToPlane(
                const simb::MCTrajectory& trajectory,
                const double& planeZ)
{
  int result = -1;
  double minDistance = 1e20;
  for(size_t iPoint=0; iPoint < trajectory.size(); iPoint++)
  {
    const auto& thisDistance = fabs(trajectory.Z(iPoint)-planeZ);
    if (thisDistance < minDistance)
    {
      minDistance = thisDistance;
      result = iPoint;
    }
  }
  return result;
}
