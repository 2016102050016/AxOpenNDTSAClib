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
#ifndef Q_DATAUTILITY_HEADER
#define Q_DATAUTILITY_HEADER
#include <ccPointCloud.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_cloud.h>

inline void CC2PCL_PointCloud(ccPointCloud& theCloud, pcl::PointCloud< pcl::PointXYZ> &pc)
{
	unsigned pointCount = theCloud.size();
	pc.resize(pointCount);

	for (unsigned i = 0; i < pointCount; ++i)
	{
		const CCVector3 *P = theCloud.getPoint(i);
		pc.at(i).x = P->x;
		pc.at(i).y = P->y;
		pc.at(i).z = P->z;
	}
}

inline void CCVec2PCL_PointCloud(std::vector<CCVector3>& theCloud, pcl::PointCloud< pcl::PointXYZ> &pc)
{
	unsigned pointCount = theCloud.size();
	pc.resize(pointCount);

	for (unsigned i = 0; i < pointCount; ++i)
	{
		CCVector3 P = theCloud.at(i);
		pc.at(i).x = P.x;
		pc.at(i).y = P.y;
		pc.at(i).z = P.z;
	}
}

inline void CCVec2PCL_PointXYZVec(std::vector<CCVector3>& theCloud, std::vector< pcl::PointXYZ> &pc)
{
	unsigned pointCount = theCloud.size();
	pc.resize(pointCount);

	for (unsigned i = 0; i < pointCount; ++i)
	{
		CCVector3 P = theCloud.at(i);
		pc.at(i).x = P.x;
		pc.at(i).y = P.y;
		pc.at(i).z = P.z;
	}
}

inline ccPointCloud* PCL2CC_PointCloud(pcl::PointCloud<pcl::PointXYZ> &ptcl, QString group_name)
{
	ccPointCloud* ccCloud = new ccPointCloud();
	int pointCount = ptcl.size();
	if (!ccCloud->reserve(static_cast<unsigned>(pointCount)))
		return NULL;
	for (size_t i = 0; i < pointCount; ++i)
	{
		CCVector3 P(ptcl.at(i).x, ptcl.at(i).y, ptcl.at(i).z);
		ccCloud->addPoint(P);
	}
	ccCloud->setName(group_name);
}

inline bool SortEigenValuesAndVectors(Eigen::Matrix3d& eigenVectors, Eigen::Vector3d& eigenValues)
{
	if (eigenVectors.cols() < 2 || eigenVectors.cols() != eigenValues.rows())
	{
		assert(false);
		return false;
	}

	unsigned n = eigenVectors.cols();
	for (unsigned i = 0; i < n - 1; i++)
	{
		unsigned maxValIndex = i;
		for (unsigned j = i + 1; j<n; j++)
			if (eigenValues[j] > eigenValues[maxValIndex])
				maxValIndex = j;

		if (maxValIndex != i)
		{
			std::swap(eigenValues[i], eigenValues[maxValIndex]);
			for (unsigned j = 0; j < n; ++j)
				std::swap(eigenVectors(j, i), eigenVectors(j, maxValIndex));
		}
	}

	return true;
}
#endif