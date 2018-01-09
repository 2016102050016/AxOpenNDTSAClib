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
#include "ndtRansacHybridEx.h"
#include "planefit.h"


ndtRansacHybridEx::ndtRansacHybridEx(void):
	theOctree(NULL),theCloud(NULL),group(NULL)
	,m_indices(NULL) 
	,model_()
	,inliers_ ()
	,model_coefficients_ ()
	,probability_ (0.01)
	,iterations_ (0)
	,max_iterations_ (1000)
	,error_sqr_dists_ ()
	, minPointsPerCell(10)
	,minPointsPerComponent(1000)
{
	//输入参数，半径大小
	kernelRadius=1;
	teThreshold=0.01;
	minSupportCell=5;
	maxDistancetoPlane=0.01;
	cellDist2Plane=0.001;
	maxNormDevAngle=25;
}


void ndtRansacHybridEx::InitializeComptutation()
{
	theCloud->getBoundingBox(bbMin,bbMax);
	const CCVector3d& globalShift = theCloud->getGlobalShift();
	double globalScale = theCloud->getGlobalScale();

	CCVector3 diff = bbMax - bbMin;
	scale=std::max(std::max(diff[0], diff[1]), diff[2]);
	bitmapEpsilonDouble= .01f * scale;
}


ndtRansacHybridEx::~ndtRansacHybridEx(void)
{
	theOctree=NULL;
	theCloud=NULL;
	group=NULL;
}

void ndtRansacHybridEx::extract(GenericProgressCallback* progressCb)
{
	InitializeComptutation();
	unsigned count = theCloud->size();
	bool hasNorms = theCloud->hasNormals();

	double angleThreshold = static_cast<float>(cos(maxNormDevAngle * CC_DEG_TO_RAD));

	planeCellsAtLevelEx.clear();

	unsigned numberOfPoints = theCloud->size();
	if (numberOfPoints < 5)
	{
		return;
	}
	//构建八叉树
	if (!theOctree)
	{
		theOctree = new CCNDTLib::DgmOctreeNDT(theCloud);
		if (theOctree->build(progressCb) < 1)
		{
			delete theOctree;
			return;
		}
	}
	int result = 0;
	QString sfName = "Eigen_Value";
	int sfIdx = -1;
	ccPointCloud* pc = 0;
	//注意先添加ScalarField，并设置为当前的
	if (theCloud->isA(CC_TYPES::POINT_CLOUD))
	{
		pc = static_cast<ccPointCloud*>(theCloud);

		sfIdx = pc->getScalarFieldIndexByName(qPrintable(sfName));
		if (sfIdx < 0)
			sfIdx = pc->addScalarField(qPrintable(sfName));
		if (sfIdx >= 0)
			pc->setCurrentInScalarField(sfIdx);
		else
		{
			//ccConsole::Error(QString("Failed to create scalar field on cloud '%1' (not enough memory?)").arg(pc->getName()));
		}
	}
	//开启ScalarField模式
	theCloud->enableScalarField();
	//给定的半径，寻找最佳的Level
	unsigned char level = theOctree->findBestLevelForAGivenNeighbourhoodSizeExtraction(kernelRadius);

	//number of cells for this level
	unsigned cellCount = theOctree->getCellNumber(level);
	unsigned maxCellPopulation = theOctree->getmaxCellPopulation(level);

	//cell descriptor (initialize it with first cell/point)
	CCNDTLib::DgmOctreeNDT::octreeCellNDT cell(theOctree);
	if (!cell.points->reserve(maxCellPopulation)) //not enough memory
		return;
	cell.level = level;
	cell.index = 0;

	CCNDTLib::DgmOctreeNDT::cellIndexesContainer vec;
	try
	{
		vec.resize(cellCount);
	}
	catch (const std::bad_alloc&)
	{
		return;
	}
	QElapsedTimer eTimer;
	eTimer.start();
	unsigned segmentedPlanes = 0;
	std::vector<AxPointXYZNormal> ccNonPlaneCloud;//非面的点集
	//binary shift for cell code truncation
	unsigned char bitDec = GET_NDT_BIT_SHIFT(level);
	CCNDTLib::DgmOctreeNDT::cellsContainer::const_iterator p = theOctree->pointsAndTheirCellCodes().begin();
	CCNDTLib::DgmOctreeNDT::OctreeCellCodeType predCode = (p->theCode >> bitDec) + 1; //pred value must be different than the first element's
	//Cell最小点数
	std::vector<unsigned> cellPtsCount;
	//所有的点对应的Point-Cell，n：1，theOctree->getNumberOfProjectedPoints()已经排序了
	for (unsigned i = 0, j = 0; i < theOctree->getNumberOfProjectedPoints(); ++i, ++p)
	{
		CCNDTLib::DgmOctreeNDT::OctreeCellCodeType currentCode = (p->theCode >> bitDec);

		if (predCode != currentCode)
		{
			vec[j++] = i;//存储索引
			ScalarType curv = NAN_VALUE;
			int n = cell.points->size();
			if (n > 0 && n < minPointsPerCell)
			{
				CCVector3d Psum(0, 0, 0);
				for (unsigned ix = 0; ix < n; ++ix)
				{
					CCVector3 P;
					cell.points->getPoint(ix, P);
					Psum.x += P.x;
					Psum.y += P.y;
					Psum.z += P.z;
				}
				CCVector3 G(static_cast<PointCoordinateType>(Psum.x / n), static_cast<PointCoordinateType>(Psum.y / n),
					static_cast<PointCoordinateType>(Psum.z / n));
				CCNDTLib::DgmOctreeNDT::NearestNeighboursSphericalSearchStruct nNSS;
				nNSS.level = cell.level;
				nNSS.queryPoint = G;
				nNSS.prepare(2 * kernelRadius, cell.parentOctree->getCellSize(nNSS.level));
				cell.parentOctree->getCellPos(cell.truncatedCode, cell.level, nNSS.cellPos, true);
				cell.parentOctree->computeCellCenter(nNSS.cellPos, cell.level, nNSS.cellCenter);
				double radius = 0.01f * scale;
				if (kernelRadius > 0.01f * scale)
				{
					radius = 1.414*kernelRadius;
				}
				unsigned neighborCount = theOctree->findNeighborsInASphereStartingFromCell(nNSS, radius, false);
				if (neighborCount < minPointsPerCell)
				{
					radius = 2 * kernelRadius;
					neighborCount = theOctree->findNeighborsInASphereStartingFromCell(nNSS, radius, false);
				}
				/*unsigned neighborCount2=theOctree->findTheNearestNeighborStartingFromCell(nNSS);
				nNSS.minNumberOfNeighbors=5;
				unsigned neighborCount3=theOctree->findNearestNeighborsStartingFromCell(nNSS,false);*/
				if (neighborCount >= 5)
				{
					CCNDTLib::DgmOctreeNDTReferenceCloud neighboursCloud(&nNSS.pointsInNeighbourhood, neighborCount);
					Neighbourhood Z(&neighboursCloud);
					const CCVector3* Gr = Z.getGravityCenter();
					CCLib::SquareMatrixd C = Z.computeCovarianceMatrix();

					CCLib::SquareMatrixd eigVectors;
					std::vector<double> eigValues;
					if (!Jacobi<double>::ComputeEigenValuesAndVectors(C, eigVectors, eigValues))
					{
						curv = NAN_VALUE;
					}
					Jacobi<double>::SortEigenValuesAndVectors(eigVectors, eigValues);
					//compute eigen value
					double e0 = eigValues[0];
					double e1 = eigValues[1];
					double e2 = eigValues[2];//排序之后，e0最大，e2最小
					CCVector3 eigvalues1(e0, e1, e2);
					cell.evals_ = eigvalues1;
					if (e2 / e1 < teThreshold)//是面特征
					{
						CellDepict cellplaneDis;
						cellplaneDis.index = cell.index;
						cellplaneDis.center = G;
						cellplaneDis.mean_ = G;
						cellplaneDis.cov_ = C;
						cellplaneDis.evals_ = eigvalues1;
						cellplaneDis.evecs_ = eigVectors;
						double uOther[3];
						Jacobi<double>::GetEigenVector(eigVectors, 2, uOther);
						CCVector3 normal(uOther[0], uOther[1], uOther[2]);
						cellplaneDis.normal = normal;//赋值法向量

						cellplaneDis.cellIndex = i;
						for (unsigned ix = 0; ix < n; ++ix)
						{
							//current point index
							unsigned index = cell.points->getPointGlobalIndex(ix);
							CCVector3 CellP;
							cell.points->getPoint(ix, CellP);

							cellplaneDis.points.push_back(CellP);
							curv = 1;
							cell.points->setPointScalarValue(ix, curv);
						}
						planeCellsAtLevelEx.push_back(cellplaneDis);//添加到Cell集合
						cellPtsCount.push_back(n);
					}
					else//非面特征
					{
						for (unsigned ix = 0; ix < n; ++ix)
						{
							//current point index
							unsigned index = cell.points->getPointGlobalIndex(ix);
							CCVector3 CellP;
							cell.points->getPoint(ix, CellP);
							AxPointXYZNormal pt;
							pt.XYZ = CellP;
							double uOther[3];
							Jacobi<double>::GetEigenVector(eigVectors, 2, uOther);
							pt.normal = CCVector3(uOther[0], uOther[1], uOther[2]);
							ccNonPlaneCloud.push_back(pt);
							curv = 0;
							cell.points->setPointScalarValue(ix, curv);
						}
					}
				}

			}
			else if (n >= minPointsPerCell)//Cell中的点数大于等于10
			{
				CCVector3d Psum(0, 0, 0);
				for (unsigned ix = 0; ix < n; ++ix)
				{
					CCVector3 P;
					cell.points->getPoint(ix, P);
					Psum.x += P.x;
					Psum.y += P.y;
					Psum.z += P.z;
				}

				CCVector3 G(static_cast<PointCoordinateType>(Psum.x / n),
					static_cast<PointCoordinateType>(Psum.y / n),
					static_cast<PointCoordinateType>(Psum.z / n));

				double mXX = 0.0;
				double mYY = 0.0;
				double mZZ = 0.0;
				double mXY = 0.0;
				double mXZ = 0.0;
				double mYZ = 0.0;
				//for each point in the cell
				for (unsigned ix = 0; ix < n; ++ix)
				{
					CCVector3 CellP;
					cell.points->getPoint(ix, CellP);
					CCVector3 P = CellP - G;

					mXX += static_cast<double>(P.x)*P.x;
					mYY += static_cast<double>(P.y)*P.y;
					mZZ += static_cast<double>(P.z)*P.z;
					mXY += static_cast<double>(P.x)*P.y;
					mXZ += static_cast<double>(P.x)*P.z;
					mYZ += static_cast<double>(P.y)*P.z;
				}
				//协方差矩阵对称
				CCLib::SquareMatrixd covMat(3);
				covMat.m_values[0][0] = mXX / n;
				covMat.m_values[1][1] = mYY / n;
				covMat.m_values[2][2] = mZZ / n;
				covMat.m_values[1][0] = covMat.m_values[0][1] = mXY / n;
				covMat.m_values[2][0] = covMat.m_values[0][2] = mXZ / n;
				covMat.m_values[2][1] = covMat.m_values[1][2] = mYZ / n;

				CCLib::SquareMatrixd eigVectors;
				std::vector<double> eigValues;
				if (!Jacobi<double>::ComputeEigenValuesAndVectors(covMat, eigVectors, eigValues))
				{
					curv = NAN_VALUE;//failure
				}
				Jacobi<double>::SortEigenValuesAndVectors(eigVectors, eigValues);
				//compute eigen value
				double e0 = eigValues[0];
				double e1 = eigValues[1];
				double e2 = eigValues[2];//排序之后，e0最大，e2最小

				cell.mean_ = G;
				cell.cov_ = covMat;
				CCVector3 eigvalues1(e0, e1, e2);
				cell.evals_ = eigvalues1;

				if (e2 / e1 < teThreshold)//是面特征
				{
					CellDepict cellplaneDis;
					cellplaneDis.index = cell.index;
					cellplaneDis.center = G;
					cellplaneDis.mean_ = G;
					cellplaneDis.cov_ = covMat;
					cellplaneDis.evals_ = eigvalues1;
					cellplaneDis.evecs_ = eigVectors;
					double uOther[3];
					Jacobi<double>::GetEigenVector(eigVectors, 2, uOther);
					CCVector3 normal(uOther[0], uOther[1], uOther[2]);

					cellplaneDis.cellIndex = i;

					double* arr = new double[n];//转存距离参数，用于去中值
					float  sumd = 0;
					for (unsigned ix = 0; ix < n; ++ix)
					{
						//current point index
						unsigned index = cell.points->getPointGlobalIndex(ix);
						CCVector3 CellP;
						cell.points->getPoint(ix, CellP);
						float distance = fabs(normal.dot(CellP - G));//点到平面的距离
						arr[ix] = sqrt(distance);
						sumd += distance;
					}
					std::sort(arr, arr + n);
					double median = arr[n / 2];//中值
					double* arrTheta = new double[n];
					for (unsigned ix = 0; ix < n; ++ix)
					{
						arrTheta[ix] = abs(arr[ix] - median);
					}
					std::sort(arrTheta, arrTheta + n);
					double medianTheta = arrTheta[n / 2];
					double scale0 = medianTheta / .6745;

					double kRob = 2.985;;
					double *wghs = new double[n];
					float meand = sumd / n;
					float S = 0;
					for (unsigned ix = 0; ix < n; ++ix)
					{
						//current point index
						unsigned index = cell.points->getPointGlobalIndex(ix);
						CCVector3 CellP;
						cell.points->getPoint(ix, CellP);
						float distance = fabs(normal.dot(CellP - G));//点到平面的距离
						float distance_n = distance / scale0;//标准化	
						wghs[ix] = exp(-pow(distance_n / kRob, 2));

						S = S + (distance - meand)*(distance - meand);
						cellplaneDis.points.push_back(CellP);
						curv = 1;
						cell.points->setPointScalarValue(ix, curv);
					}
					S = sqrt(S / (n - 1));

					std::vector<PointXYZ> tmpPointCloud;
					planefit pl;
					CCVec2PCL_PointXYZVec(cellplaneDis.points, tmpPointCloud);
					pl.setInputPointCloud(tmpPointCloud);
					pl.extractLmed();
					AxPlane  plane0 = pl.getPlane();
					Eigen::Vector3d vn = plane0.getNormal();
					cellplaneDis.normal = CCVector3(vn[0], vn[1], vn[2]);//赋值法向量

					//m_app->dispToConsole(QString("[Compute Covariance]: %6f m.").arg(S),ccMainAppInterface::STD_CONSOLE_MESSAGE);
					planeCellsAtLevelEx.push_back(cellplaneDis);//添加到Cell集合
					cellPtsCount.push_back(n);
					delete[] arr;
					delete[] arrTheta;
					delete[] wghs;
				}
				else//非面特征
				{
					for (unsigned ix = 0; ix < n; ++ix)
					{
						//current point index
						unsigned index = cell.points->getPointGlobalIndex(ix);
						CCVector3 CellP;
						cell.points->getPoint(ix, CellP);
						AxPointXYZNormal pt;
						pt.XYZ = CellP;
						double uOther[3];
						Jacobi<double>::GetEigenVector(eigVectors, 2, uOther);
						pt.normal[0] = uOther[0];
						pt.normal[1] = uOther[1];
						pt.normal[2] = uOther[2];
						ccNonPlaneCloud.push_back(pt);
						curv = 0;
						cell.points->setPointScalarValue(ix, curv);
					}
				}
			}

			//and we start a new cell开始新的Cell
			cell.index += cell.points->size();
			cell.points->clear(false);
			cell.truncatedCode = currentCode;
		}
		cell.points->addPointIndex(p->theIndex); //can't fail (see above)
		predCode = currentCode;
	}
	qint64 elapsedTime_ms0 = eTimer.elapsed();
	float time0 = static_cast<float>(elapsedTime_ms0) / 1.0e3;

	if (result == 0)
	{
		if (pc && sfIdx >= 0)
		{
			//设置当前显示的ScalarField
			pc->setCurrentDisplayedScalarField(sfIdx);
			pc->showSF(sfIdx >= 0);
			pc->getCurrentInScalarField()->computeMinAndMax();
		}
		theCloud->prepareDisplayForRefresh();
	}
	else
	{
		ccConsole::Warning(QString("Failed to apply processing to cloud '%1'").arg(theCloud->getName()));
		if (pc && sfIdx >= 0)
		{
			pc->deleteScalarField(sfIdx);
			sfIdx = -1;
		}
	}
	size_t drawnCandidates = 0;
	do //执行RANSAC流程，检测面
	{
		if (!m_indices || m_indices->size() == 0)
		{
			m_indices = new std::vector <unsigned>(static_cast<int> (planeCellsAtLevelEx.size()));
			int index = 0;
			for (std::vector <unsigned>::iterator it = m_indices->begin(), it_e = m_indices->end(); it != it_e; it++)
				*it = index++;
		}
		else
		{
			m_indices->clear();
			m_indices->resize(planeCellsAtLevelEx.size());
			int index = 0;
			for (std::vector <unsigned>::iterator it = m_indices->begin(), it_e = m_indices->end(); it != it_e; it++)
				*it = index++;
		}
		std::sort(cellPtsCount.begin(), cellPtsCount.end(), std::greater<unsigned>());//降序排列
		unsigned maxOverlapCount = static_cast<unsigned>(1.0 *cellPtsCount.size());//注意这个数值
		unsigned maxOverlapCellCount = cellPtsCount.at(maxOverlapCount - 1);
		bool res = computeModelEx(cellDist2Plane, angleThreshold, &drawnCandidates, minSupportCell, maxOverlapCellCount);

		if (res)
		{
			if (!group)
				group = new ccHObject(QString("NDT_Ransac(%1)").arg(lyrName));
			//提取已经获取的点云Cell
			std::vector<CCVector3> m_PointCloud;
			unsigned selectCellsCount = 0;
			for (size_t i = 0; i < m_indices->size(); i++)
			{
				if (inliers_[i] == 1)
				{
					selectCellsCount++;
					unsigned indx = (*m_indices)[i];
					std::vector<CCVector3> pts = planeCellsAtLevelEx[indx].points;
					for (std::vector<CCVector3>::const_iterator it = pts.begin(); it < pts.end(); it++)
					{
						CCVector3 pt = (CCVector3)(*it);
						m_PointCloud.push_back(pt);
					}

				}
			}
			/*if (m_PointCloud.size()<maxOverlapCellCount*selectCellsCount)
			{
			continue;
			}*/

			CCVector3 C, N;
			//IRLS拟合平面
			planefit plane;
			std::vector<PointXYZ> tmpPointCloud;
			planefit pl;
			CCVec2PCL_PointXYZVec(m_PointCloud, tmpPointCloud);
			pl.setInputPointCloud(tmpPointCloud);
			pl.extractLmed();
			AxPlane  plane0 = pl.getPlane();
			Eigen::Vector3d vn = plane0.getNormal();
			Eigen::Vector3d center0 = plane0.getCenter();
			//if (pPlane)
			{
				segmentedPlanes++;
				unsigned nonPlanePtCount = ccNonPlaneCloud.size();
				std::vector<unsigned> inliersNonPlanePoint;
				inliersNonPlanePoint.resize(nonPlanePtCount);
				for (size_t idx = 0; idx < nonPlanePtCount; idx++)
				{
					inliersNonPlanePoint[idx] = 0;//初始化标签，标识点云是否已经分类
				}
				CCVector3 G(center0[0], center0[1], center0[2]);
				CCVector3 normaltmp(vn[0], vn[1], vn[2]);

				//判断剩余点是否在平面内	
				for (size_t i = 0; i < nonPlanePtCount; i++)
				{
					AxPointXYZNormal pt = ccNonPlaneCloud[i];
					float distance = fabsf(normaltmp.dot(pt.XYZ - G));
					double angletest = normaltmp.dot(pt.normal);
					//if (distance < maxDistancetoPlane  && angletest > angleThreshold)
					if (distance < maxDistancetoPlane)
					{
						m_PointCloud.push_back(pt.XYZ);
						inliersNonPlanePoint[i] = 1;//标识为面上点
					}
				}

				unsigned lastPoint = 0;
				for (unsigned i = 0; i < nonPlanePtCount; ++i)
				{
					//i持续增长，而lastPoint遇到!=POINT_VISIBLE则跳过，起到迁移的效果
					if (inliersNonPlanePoint[i] == 0)
					{
						if (i != lastPoint)
						{
							std::swap(ccNonPlaneCloud[lastPoint], ccNonPlaneCloud[i]);
						}
						++lastPoint;
					}
				}
				ccNonPlaneCloud.resize(lastPoint);

				size_t numOfPlane = m_PointCloud.size();
				//似乎仅仅是连通分析的需要
				std::vector<size_t> * m_indicesCurrent = new std::vector<size_t>();
				m_indicesCurrent->resize(numOfPlane);
				for (size_t i = 0; i < numOfPlane; i++)
				{
					(*m_indicesCurrent)[i] = i;
				}

				ccPointCloud *ccCloud = new ccPointCloud();
				if (!ccCloud->reserve(static_cast<unsigned>(m_indicesCurrent->size())))
					return;
				for (unsigned j = 0; j < m_indicesCurrent->size(); ++j)
				{
					size_t idx = (*m_indicesCurrent)[j];//注意int型vector数组指针的使用
					CCVector3 vc = m_PointCloud[idx];
					ccCloud->addPoint(vc);
				}
				ccCloud->setName(QString("PlaneNoCC"));
				ccColor::Rgb col = ccColor::Generator::Random();
				ccCloud->setRGBColor(col);
				ccCloud->showColors(true);
				ccCloud->setPointSize(1);
				group->addChild(ccCloud);
				//连通分析
				CCNDTLib::DgmOctreeNDT *theOctreeTmp = new CCNDTLib::DgmOctreeNDT(ccCloud);
				if (theOctreeTmp->build(NULL) < 1)
				{
					return;
				}
				if (theOctreeTmp)
				{
					ccCloud->enableScalarField();
					theOctreeTmp->extractCCs(level, false, NULL);
					std::vector<ReferenceCloud*> components;
					unsigned numberOfPoints = (ccCloud ? ccCloud->size() : 0);
					if (numberOfPoints == 0)
						return;

					for (unsigned i = 0; i < numberOfPoints; ++i)
					{
						ScalarType slabel = ccCloud->getPointScalarValue(i);
						if (slabel >= 1) //labels start from 1! (this test rejects NaN values as well)
						{
							int ccLabel = static_cast<int>(ccCloud->getPointScalarValue(i)) - 1;
							try
							{
								while (static_cast<size_t>(ccLabel) >= components.size())
									components.push_back(new ReferenceCloud(ccCloud));
							}
							catch (const std::bad_alloc&)
							{
								components.clear();
								return;
							}

							//add the point to the current component
							if (!components[ccLabel]->addPointIndex(i))
							{
								//not enough memory
								while (!components.empty())
								{
									delete components.back();
									components.pop_back();
								}
								return;
							}
						}
					}
					if (components.size() > 0)
					{
						std::vector<ComponentIndexAndSize2> sortedIndexes;
						std::vector<ComponentIndexAndSize2>* _sortedIndexes = 0;
						bool sortBysize = true;
						try
						{
							sortedIndexes.reserve(components.size());
						}
						catch (const std::bad_alloc&)
						{
							ccLog::Warning("[CreateComponentsClouds] Not enough memory to sort components by size!");
							sortBysize = false;
						}

						if (sortBysize) //still ok?
						{
							unsigned compCount = static_cast<unsigned>(components.size());
							for (unsigned i = 0; i < compCount; ++i)
							{
								sortedIndexes.push_back(ComponentIndexAndSize2(i, components[i]->size()));
							}

							std::sort(sortedIndexes.begin(), sortedIndexes.end(), ComponentIndexAndSize2::DescendingCompOperator);
							_sortedIndexes = &sortedIndexes;
						}
						for (unsigned ik = 0; ik < components.size(); ++ik)
						{
							CCLib::ReferenceCloud* compIndexes = _sortedIndexes ? components[_sortedIndexes->at(ik).index] : components[ik];
							if (compIndexes->size() >= minPointsPerComponent)
							{
								ccPointCloud* ccPlaneCCCloud = (ccCloud ? ccCloud->partialClone(compIndexes) : ccPointCloud::From(compIndexes));

								ccPlaneCCCloud->setName(QString("Plane"));
								ccColor::Rgb colCC = ccColor::Generator::Random();
								ccPlaneCCCloud->setRGBColor(colCC);
								ccPlaneCCCloud->showColors(true);
								ccPlaneCCCloud->setPointSize(1);
								group->addChild(ccPlaneCCCloud);//将产生的点云增加到分组
							}
							else
							{
								ccPointCloud* ccPlaneCCCloud = (ccCloud ? ccCloud->partialClone(compIndexes) : ccPointCloud::From(compIndexes));
								for (int jk = 0; jk < ccPlaneCCCloud->size(); jk++)
								{
									AxPointXYZNormal pt;
									const CCVector3* vv = ccPlaneCCCloud->getPoint(jk);
									pt.XYZ[0] = vv->x;
									pt.XYZ[1] = vv->y;
									pt.XYZ[2] = vv->z;
									ccNonPlaneCloud.push_back(pt);
								}
							}
						}
					}
				}
				delete theOctreeTmp;
			}

			unsigned count = m_indices->size();
			unsigned lastPoint = 0;
			for (unsigned i = 0; i < count; ++i)
			{
				//i持续增长，而lastPoint遇到!=POINT_VISIBLE则跳过，起到迁移的效果
				if (inliers_[i] == 0)
				{
					if (i != lastPoint)
					{
						std::swap(planeCellsAtLevelEx.at(lastPoint), planeCellsAtLevelEx.at(i));
					}
					++lastPoint;
				}
			}
			planeCellsAtLevelEx.resize(lastPoint);
		}
	} while (CandidateFailureProbability(1, planeCellsAtLevelEx.size(), drawnCandidates) > probability_  && planeCellsAtLevelEx.size() > 5);
	//主要上面应该是while (CandidateFailureProbability(minSupportCell, planeCellsAtLevelEx.size(),drawnCandidates)> probability_）
	qint64 elapsedTime_ms = eTimer.elapsed();
	float timell = static_cast<float>(elapsedTime_ms) / 1.0e3;

	QString log("D:\\ndtsaclog.txt");
	FILE* in = fopen(log.toStdString().c_str(), "a");
	fprintf(in, "NDT_SACH, %s, CELLS_t: %f, SAC_t: %f, PLANES: %d, CELL_S: %f, MINSUPPORT: %d, CELL_D: %f,POINT_D: %f, NORMALTHRES: %f\n", lyrName.toStdString().c_str(), time0, timell, segmentedPlanes,
		kernelRadius, minSupportCell, cellDist2Plane, maxDistancetoPlane, angleThreshold);
	fclose(in);
	//保存剩余的点
	ccPointCloud *ccRemainedCloud = NULL;
	if (planeCellsAtLevelEx.size() > 0)
	{
		int left_num_points = 0;
		for (int i = 0; i < planeCellsAtLevelEx.size(); i++)
		{
			left_num_points += planeCellsAtLevelEx[i].points.size();
		}
		if (ccNonPlaneCloud.size()>0)
		{
			left_num_points += ccNonPlaneCloud.size();
		}
		if (!group)
			group = new ccHObject(QString("NDT_Ransac(%1)").arg(lyrName));
		ccRemainedCloud = new ccPointCloud();
		if (!ccRemainedCloud->reserve(static_cast<unsigned>(left_num_points)))
			return;
		for (int i = 0; i < planeCellsAtLevelEx.size(); i++)
		{
			std::vector<CCVector3> pts = planeCellsAtLevelEx[i].points;
			for (std::vector<CCVector3>::const_iterator it = pts.begin(); it < pts.end(); it++)
			{
				CCVector3 pt = (CCVector3)(*it);
				ccRemainedCloud->addPoint(pt);
			}
		}
		if (ccNonPlaneCloud.size() > 0)
		{
			for (unsigned j = 0; j < ccNonPlaneCloud.size(); ++j)
			{
				AxPointXYZNormal vc = ccNonPlaneCloud[j];
				ccRemainedCloud->addPoint(vc.XYZ);
			}
		}
		ccRemainedCloud->setName(QString("Remained"));
		ccColor::Rgb col = ccColor::Generator::Random();
		ccRemainedCloud->setRGBColor(col);
		ccRemainedCloud->showColors(true);
		ccRemainedCloud->setPointSize(1);
		group->addChild(ccRemainedCloud);
	}

}

int ndtRansacHybridEx::countWithinDistanceEx(CellDepict CellDis, double threshold_,float normalThresh)
{
	int nr_p = 0;

	// Iterate through the 3d points and calculate the distances from them to the plane
	for (size_t i = 0; i < m_indices->size (); ++i)
	{
		// Calculate the distance from the point to the plane normal as the dot product
		// D = (P-A).N/|N|
		CCVector3 pt (planeCellsAtLevelEx[(*m_indices)[i]].center.x,planeCellsAtLevelEx[(*m_indices)[i]].center.y,planeCellsAtLevelEx[(*m_indices)[i]].center.z);
		CCVector3 n=planeCellsAtLevelEx[(*m_indices)[i]].normal;//其它的Cell

		CCVector3 a(CellDis.center.x,CellDis.center.y,CellDis.center.z);
		CCVector3 normal=CellDis.normal;

		double distance=fabs(normal.dot(pt-a));
		double angletest =normal.dot(n)/(normal.norm()*n.norm());
		if (distance < threshold_ && fabs(angletest) > normalThresh)
			nr_p++;
	}
	return (nr_p);
}

bool ndtRansacHybridEx::computeModelEx(double threshold_,double angleThreshold,size_t *drawnCandidates,unsigned minSupportCells,unsigned maxOverlapCellCount_)
{
	iterations_ = 0;
	int n_best_inliers_count = -INT_MAX;
	double k = 1.0;
	std::vector<int> selection;
	unsigned slectionPoints=0;
	model_.clear();

	//double log_probability  = log (1.0 - probability_);//probability_=0.99
	double log_probability  = log (probability_);
	double one_over_indices = 1.0 / static_cast<double> (m_indices->size ());

	int n_inliers_count = 0;

	while (iterations_ < k )
	{
		selection.clear();
		srand(time(NULL));
		unsigned selectedPointIndex =-1;
		size_t iter = 0;
		CellDepict	CellDis;
		do
		{
			selectedPointIndex= static_cast<unsigned>(rand() % m_indices->size ());	//随机选择一个Cell	
			CellDis= planeCellsAtLevelEx.at(selectedPointIndex);
			slectionPoints=CellDis.points.size();
		}
		while(slectionPoints<maxOverlapCellCount_ && iter++ < 40);
		selection.push_back(selectedPointIndex);

		++iterations_;
		n_inliers_count = countWithinDistanceEx (CellDis, threshold_,angleThreshold);	

		// Better match ?
		if (n_inliers_count > n_best_inliers_count)
		{
			n_best_inliers_count = n_inliers_count;

			// Save the current model/inlier/coefficients selection as being the best so far
			model_              = selection;
			model_coefficients_Ex = CellDis;

			// Compute the k parameter (k=log(z)/log(1-w^n))
			double w = static_cast<double> (n_best_inliers_count) * one_over_indices;
			double p_no_outliers = 1.0 - pow (w, static_cast<double> (selection.size ()));
			p_no_outliers = (std::max) (std::numeric_limits<double>::epsilon (), p_no_outliers);       // Avoid division by -Inf
			p_no_outliers = (std::min) (1.0 - std::numeric_limits<double>::epsilon (), p_no_outliers);   // Avoid division by 0.
			k = log_probability / log (p_no_outliers);
		}
		if (iterations_ > max_iterations_)
		{
			break;
		}
	}
	*drawnCandidates += 1;
	//这儿还要考虑一下是否合理
	if (model_.empty ()|| n_best_inliers_count<= minSupportCells)
	{
		inliers_.clear ();
		return (false);
	}

	// Get the set of inliers that correspond to the best model found so far
	selectWithinDistanceEx (model_coefficients_Ex, threshold_,angleThreshold, inliers_);
	return true;
}

void ndtRansacHybridEx::selectWithinDistanceEx (CellDepict  model_coefficients, const double threshold,float normalThresh,std::vector<unsigned> &inliers)
{
	int nr_p = 0;
	inliers.resize (m_indices->size ());
	for (size_t idx=0;idx<m_indices->size();idx++)
	{
		inliers[idx]=0;
	}
	error_sqr_dists_.resize (m_indices->size ());

	// Iterate through the 3d points and calculate the distances from them to the plane
	for (size_t i = 0; i < m_indices->size (); ++i)
	{
		// Calculate the distance from the point to the plane normal as the dot product
		// D = (P-A).N/|N|
		CCVector3 pt (planeCellsAtLevelEx[(*m_indices)[i]].center.x,planeCellsAtLevelEx[(*m_indices)[i]].center.y,planeCellsAtLevelEx[(*m_indices)[i]].center.z);
		CCVector3 n=planeCellsAtLevelEx[(*m_indices)[i]].normal;//其它的Cell

		CCVector3 a(model_coefficients.center.x,model_coefficients.center.y,model_coefficients.center.z);
		CCVector3 normal=model_coefficients.normal;

		float distance=fabsf(normal.dot(pt-a));
		double angletest =normal.dot(n)/(normal.norm()*n.norm());//法向量没有标准化
		if (distance < threshold && fabs(angletest) > normalThresh)
		{
			// Returns the indices of the points whose distances are smaller than the threshold
			inliers[i] = 1;
			error_sqr_dists_[nr_p] = static_cast<double> (distance);
			++nr_p;
		}
	}
	//inliers.resize (nr_p);
	error_sqr_dists_.resize (nr_p);
}
