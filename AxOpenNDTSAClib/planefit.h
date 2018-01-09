//##########################################################################
//#                                                                        #
//#                       CLOUDCOMPARE PLUGIN: qOpenNDTSACPlugin           #
//#                                                                        #
//#  This program is free software; you can redistribute it and/or modify  #
//#  it under the terms of the GNU General Public License as published by  #
//#  the Free Software Foundation; version 2 of the License.               #
//#                                                                        #
//#  This program is distributed in the hope that it will be useful,       #
//#  but WITHOUT ANY WARRANTY; without even the implied warranty of        #
//#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
//#  GNU General Public License for more details.                          #
//#                                                                        #
//#                             COPYRIGHT: Fan Yang                        #
//#                                                                        #
//##########################################################################
#ifndef AX_PLANEFIT
#define AX_PLANEFIT
#pragma once
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include "dataUtility.h"
#include "Eigen\Eigenvalues"
#include <vector>
#define SMALL_NUM  0.00000001
using namespace pcl;
class lxLine 
{
public:
	CCVector3 P0, P1;//直线上两个点
	std::vector<CCVector3> *lxPointCloud;
	lxLine():lxPointCloud(NULL)
	{}
	lxLine Compute()
	{

	}
	~lxLine(){}
};
//直线段
class lxLineSegment
{
public:
	CCVector3 P0, P1;
	std::vector<CCVector3> *lxPointCloud;
	lxLineSegment():lxPointCloud(NULL)
	{}
	lxLineSegment Compute()
	{

	}
	double Distance()
	{
		double d=(P1-P0).norm();
		return d;
	}
	CCVector3 Direction()
	{
		return P1-P0;
	}
	~lxLineSegment(){}
};

class AxPlane
{
public:
	double nx;
	double ny;
	double nz;
	double theta;
	double phi;
	double rho;
	Eigen::Vector3d center_;
	Eigen::Vector3d evals;
	Eigen::Matrix3d evecs;
	void setCenter(double x,double y,double z)
	{
		center_[0]=x;
		center_[1]=y;
		center_[2]=z;
	}
	void setNormal(double nx_,double ny_,double nz_)
	{
		nx= nx_;
		ny= ny_;
		nz= nz_;
	}
	Eigen::Vector3d getCenter()
	{
		return center_;
	}
	Eigen::Vector3d getNormal()
	{
		Eigen::Vector3d normal1;
		normal1[0]=nx;
		normal1[1]=ny;
		normal1[2]=nz;
		return normal1;
	}
	std::vector<pcl::PointXYZ> points_;
	 AxPlane():	nx(-1)
		 ,ny(-1)
		 ,nz(-1)
		 ,theta(-1)
		 ,phi(-1)
		 ,rho(-1)
		 ,angleVertialThreshold(2)
	 {

	 }
	 Eigen::Vector3d getLSPlaneX()
	 {
		 Eigen::Vector3d x;
		 x[0]=evecs(0,0);
		 x[1]=evecs(1,0);
		 x[2]=evecs(2,0);
		 return x;
	 }
	 float angleVertialThreshold;
	 float getAngleVertialThreshold() const { return angleVertialThreshold; }
	 void setAngleVertialThreshold(float val) { angleVertialThreshold = val; }
	 bool IsVertialPlanes()
	 {
		 Eigen::Vector3d normal1;
		 normal1[0]=nx;
		 normal1[1]=ny;
		 normal1[2]=nz;
		 Eigen::Vector3d normal2;
		 normal2[0]=0;
		 normal2[1]=0;
		 normal2[2]=1;
		 double a=normal1.dot(normal2)/normal1.norm();
		 double Threshold=cos(90-angleVertialThreshold);
		 if (abs(a)< Threshold)
		 {
			 return true;
		 }
		 return false;
	 }
	 bool IsHorizontalPlanes()
	 {
		 Eigen::Vector3d normal1;
		 normal1[0]=nx;
		 normal1[1]=ny;
		 normal1[2]=nz;
		 Eigen::Vector3d normal2;
		 normal2[0]=0;
		 normal2[1]=0;
		 normal2[2]=1;
		 double a=normal1.dot(normal2)/normal1.norm();
		 float Threshold=cos(angleVertialThreshold);
		 if (abs(a)> Threshold)
		 {
			 return true;
		 }
		 return false;
	 }
};
// intersect3D_SegmentPlane(): find the 3D intersection of a segment and a plane
//    Input:  S = a segment, and Pn = a plane = {Point V0;  Vector n;}
//    Output: *I0 = the intersect point (when it exists)
//    Return: 0 = disjoint (no intersection)
//            1 =  intersection in the unique point *I0
//            2 = the  segment lies in the plane
inline int	intersect3D_SegmentPlane( lxLineSegment S, AxPlane Pn, CCVector3* I )
{
	Eigen::Vector3d normal1=Pn.getNormal();
	CCVector3 Pn_n(normal1[0],normal1[1],normal1[2]);
	Eigen::Vector3d center1=Pn.getCenter();
	CCVector3 PnV0(center1[0],center1[1],center1[2]);
	CCVector3    u = S.P1 - S.P0;
	CCVector3    w = S.P0 - PnV0;

	float     D = Pn_n.dot(u);
	float     N = -Pn_n.dot(w);

	if (fabs(D) < SMALL_NUM) {           // segment is parallel to plane
		if (N == 0)                      // segment lies in plane
			return 2;
		else
			return 0;                    // no intersection
	}
	// they are not parallel
	// compute intersect param
	float sI = N / D;
	if (sI < 0 || sI > 1)
		return 0;                        // no intersection

	*I = S.P0 + sI * u;                  // compute segment intersect point
	return 1;
}
inline int	intersect3D_RayPlane( lxLineSegment S, AxPlane Pn, CCVector3* I )
{
	Eigen::Vector3d normal1=Pn.getNormal();
	CCVector3 Pn_n(normal1[0],normal1[1],normal1[2]);
	Eigen::Vector3d center1=Pn.getCenter();
	CCVector3 PnV0(center1[0],center1[1],center1[2]);
	CCVector3    u = S.P1 - S.P0;
	CCVector3    w = S.P0 - PnV0;

	float     D = Pn_n.dot(u);
	float     N = -Pn_n.dot(w);

	if (fabs(D) < SMALL_NUM) {           // segment is parallel to plane
		if (N == 0)                      // segment lies in plane
			return 2;
		else
			return 0;                    // no intersection
	}
	// they are not parallel
	// compute intersect param
	float sI = N / D;
	if (sI < 0 )
		return 0;                        // no intersection

	*I = S.P0 + sI * u;                  // compute segment intersect point
	return 1;
}
// intersect3D_2Planes(): the 3D intersect of two planes
//    Input:  two planes Pn1 and Pn2
//    Output: *L = the intersection line (when it exists)
//    Return: 0 = disjoint (no intersection)
//            1 = the two planes coincide
//            2 = intersection in the unique line *L
inline int	intersect3D_2Planes( AxPlane Pn1, AxPlane Pn2, lxLine* L )
{
	Eigen::Vector3d normal1=Pn1.getNormal();
	Eigen::Vector3d normal2=Pn2.getNormal();
	Eigen::Vector3d center1=Pn1.getCenter();
	Eigen::Vector3d center2=Pn2.getCenter();
	CCVector3 Pn1_n(normal1[0],normal1[1],normal1[2]);
	CCVector3 Pn2_n(normal2[0],normal2[1],normal2[2]);
	CCVector3 Pn1V0(center1[0],center1[1],center1[2]);
	CCVector3 Pn2V0(center2[0],center2[1],center2[2]);
	CCVector3   u = Pn1_n.cross(Pn2_n);         // cross product
	float    ax = (u[0] >= 0 ? u[0] : -u[0]);
	float    ay = (u[1] >= 0 ? u[1] : -u[1]);
	float    az = (u[2] >= 0 ? u[2] : -u[2]);

	// test if the two planes are parallel
	if ((ax+ay+az) < SMALL_NUM) 
	{   // Pn1 and Pn2 are near parallel
		// test if disjoint or coincide
		CCVector3   v = Pn2V0 - Pn1V0;
		if (Pn1_n.dot(v) == 0)         // Pn2.V0 lies in Pn1
			return 1;                   // Pn1 and Pn2 coincide
		else 
			return 0;                   // Pn1 and Pn2 are disjoint
	}

	// Pn1 and Pn2 intersect in a line
	// first determine max abs coordinate of cross product
	int      maxc;                      // max coordinate
	if (ax > ay) {
		if (ax > az)
			maxc = 1;
		else maxc = 3;
	}
	else {
		if (ay > az)
			maxc = 2;
		else maxc = 3;
	}

	// next, to get a point on the intersect line
	// zero the max coord, and solve for the other two
	CCVector3    iP;               // intersect point
	float    d1, d2;           // the constants in the 2 plane equations
	d1 = -Pn1_n.dot(Pn1V0);  // note: could be pre-stored with plane
	d2 = -Pn2_n.dot(Pn2V0);  // ditto

	switch (maxc) {            // select max coordinate
	case 1:                    // intersect with x=0
		iP[0]= 0;
		iP[1] = (d2 * Pn1_n[2] - d1 * Pn2_n[2]) / u[0];
		iP[2] = (d1 * Pn2_n[1] - d2 * Pn1_n[1]) / u[0];
		break;
	case 2:                    // intersect with y=0
		iP[0] = (d1 * Pn2_n[2] - d2 * Pn1_n[2]) / u[1];
		iP[1] = 0;
		iP[2] = (d2 * Pn1_n[0] - d1 * Pn2_n[0]) / u[1];
		break;
	case 3:                    // intersect with z=0
		iP[0] = (d2 * Pn1_n[1] - d1 * Pn2_n[1]) / u[2];
		iP[1] = (d1 * Pn2_n[0] - d2 * Pn1_n[0]) / u[2];
		iP[2] = 0;
	}
	L->P0 = iP;
	L->P1 = iP + u;
	return 2;
}

class planefit
{
public:

	planefit(void);
	~planefit(void);
	void setInputPointCloud(std::vector<pcl::PointXYZ> &input0_);
	void extract();
	void extractLmed();
	AxPlane getPlane();
public:
	std::vector<pcl::PointXYZ> input_;
	Eigen::Matrix3d cov_;		/// Contains the covatiance of the normal distribution
	Eigen::Matrix3d icov_;  /// Precomputed inverse covariance (updated every time the cell is updated)
	Eigen::Matrix3d evecs_; /// Eigen vectors
	Eigen::Vector3d mean_;  /// Mean of the normal distribution
	Eigen::Vector3d evals_; /// Eigen values
	AxPlane *plane_;

	double threshold2stop;
	double Threshold2stop() const { return threshold2stop; }
	void SetThreshold2stop(double val) { threshold2stop = val; }

	int iterations_;

	int max_iterations_;
	int Max_iterations() const { return max_iterations_; }
	void SetMax_iterations(int val) { max_iterations_ = val; }

	double min_iterations;
	double Min_iterations() const { return min_iterations; }
	void SetMin_iterations(double val) { min_iterations = val; }
};
#endif

