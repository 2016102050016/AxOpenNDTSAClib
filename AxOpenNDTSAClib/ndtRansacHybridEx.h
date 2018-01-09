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
#pragma once
#include <QWidget>
#include <ccHObject.h>
#include <xfunctional>
#include "ccPointCloud.h"
#include "DgmOctreeNDT.h"
#include "CCGeom.h"
#include "CellDepict.h"
#include "ccMainAppInterface.h"
#include "ccProgressDialog.h"
#include "qelapsedtimer.h"
#include "DgmOctreeNDTReferenceCloud.h"
#include "Jacobi.h"
#include<time.h>
#include "..\qCC\ccConsole.h"

class QMainWindow;

struct ComponentIndexAndSize2
{
	unsigned index;
	unsigned size;

	ComponentIndexAndSize2(unsigned i, unsigned s) : index(i), size(s) {}

	static bool DescendingCompOperator(const ComponentIndexAndSize2& a, const ComponentIndexAndSize2& b)
	{
		return a.size > b.size;
	}
};
struct AxPointXYZNormal
{
	CCVector3 XYZ;
	CCVector3 normal;
};

class ndtRansacHybridEx
{
public:
	ndtRansacHybridEx(void);
	~ndtRansacHybridEx(void);
public:
	QString lyrName;
	QString getLyrName() const { return lyrName; }
	void setLyrName(QString val) { lyrName = val; }
	//Cell的尺寸
	float kernelRadius;
	float getKernelRadius() const { return kernelRadius; }
	void setKernelRadius(float val) { kernelRadius = val; }
	//面点和非面点分类阈值
	double teThreshold;
	double getTeThreshold() const { return teThreshold; }
	void setTeThreshold(double val) { teThreshold = val; }

	int minSupportCell;
	int getMinSupportCell() const { return minSupportCell; }
	void setMinSupportCell(int val) { minSupportCell = val; }

	//
	double probability_;
	double getProbability() const { return probability_; }
	void setProbability(double val) { probability_ = val; }
	int iterations_;

	int max_iterations_;
	int getMax_iterations() const { return max_iterations_; }
	void setMax_iterations(int val) { max_iterations_ = val; }

	double cellDist2Plane;
	double getCellDist2Plane() const { return cellDist2Plane; }
	void setCellDist2Plane(double val) { cellDist2Plane = val; }

	double maxNormDevAngle;
	double getMaxNormDevAngle() const { return maxNormDevAngle; }
	void setMaxNormDevAngle(double val) { maxNormDevAngle = val; }

	double maxDistancetoPlane;
	double getMaxDistancetoPlane() const { return maxDistancetoPlane; }
	void setMaxDistancetoPlane(double val) { maxDistancetoPlane = val; }

	double bitmapEpsilonDouble;
	double scale;
	std::vector<unsigned> *m_indices;
	std::vector<unsigned> samples;

	std::vector<int> model_;
	std::vector<unsigned> inliers_;
	std::vector<double> error_sqr_dists_;

	CCNDTLib::DgmOctreeNDT::CellDescriptor model_coefficients_;
	//这一层的Cell集合
	std::vector<CCNDTLib::DgmOctreeNDT::CellDescriptor> planeCellsAtLevel;

	CellDepict model_coefficients_Ex;
	std::vector<CellDepict> planeCellsAtLevelEx;

	ccPointCloud* theCloud;//输入点云
	ccPointCloud* getInputCloud() const { return theCloud; }
	void setInputCloud(ccPointCloud* val) { theCloud = val; }

	CCNDTLib::DgmOctreeNDT* theOctree;//构建点云八叉树

	CCVector3 bbMin;
	CCVector3 bbMax;
	ccHObject* group;

	//min Points per Component for CC
	double minPointsPerComponent;
	double getMinPointsPerComponent() const { return minPointsPerComponent; }
	void setMinPointsPerComponent(double val) { minPointsPerComponent = val; }
	//min points per cell
	double minPointsPerCell;
	double getMinPointsPerCell() const { return minPointsPerCell; }
	void setMinPointsPerCell(double val) { minPointsPerCell = val; }
public:
	void InitializeComptutation();

	void extract(GenericProgressCallback* progressCb);

	float CandidateFailureProbability(float candidateSize,	float numberOfCells, float drawnCandidates) const
	{
		return std::min(std::pow(1.f - candidateSize/ numberOfCells,drawnCandidates), 1.f);
	}

	int countWithinDistanceEx(CellDepict CellDis, double threshold_,float normalThresh);
	bool computeModelEx(double threshold_,double angleThreshold,size_t *drawnCandidates,unsigned minSupportCells,unsigned maxOverlapCellCount_);
	void selectWithinDistanceEx (CellDepict  model_coefficients, const double threshold,float normalThresh,std::vector<unsigned> &inliers);
};

