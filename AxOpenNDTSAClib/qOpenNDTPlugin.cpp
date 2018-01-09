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

#include "qOpenNDTPlugin.h"

qOpenNDTPlugin::qOpenNDTPlugin(QObject* parent/*=0*/)
	: QObject(parent)
	,m_action_test(0)
	, m_action_ndtRansac(0)
	, m_action_ndtRegionGrowing(0)
{
}

void qOpenNDTPlugin::onNewSelection(const ccHObject::Container& selectedEntities)
{
	//if (m_action)
	//	m_action->setEnabled(!selectedEntities.empty());
}

void qOpenNDTPlugin::getActions(QActionGroup& group)
{
	if (!m_action_test)
	{
		m_action_test = new QAction(getName(), this);
		m_action_test->setToolTip(getDescription());
		m_action_test->setIcon(getIcon());
		//connect appropriate signal
		connect(m_action_test, SIGNAL(triggered()), this, SLOT(doTestAction()));
	}
	group.addAction(m_action_test);
	
	if (!m_action_ndtRansac)
	{
		m_action_ndtRansac = new QAction("NDTRansac", this);
		m_action_ndtRansac->setToolTip("NDT Ransac Plane segmentaion");
		m_action_ndtRansac->setIcon(QIcon(":/CC/qOpenNDTPlugin/Resources/qNdtHybridSac.png"));
		connect(m_action_ndtRansac, SIGNAL(triggered()), this, SLOT(doNDTRansac()));
	}
	group.addAction(m_action_ndtRansac);

	if (!m_action_ndtRegionGrowing)
	{
		m_action_ndtRegionGrowing = new QAction("NDTRG", this);
		m_action_ndtRegionGrowing->setToolTip("NDT RegionGrowing Plane segmentaion");
		m_action_ndtRegionGrowing->setIcon(QIcon(":/CC/qMIMSPlugin/Resources/Iron.png"));
		connect(m_action_ndtRegionGrowing, SIGNAL(triggered()), this, SLOT(doNDTRegionGrowing()));
	}
	//group.addAction(m_action_ndtRegionGrowing);
}

void qOpenNDTPlugin::doNDTRansac()
{
	assert(m_app);
	if (!m_app)
		return;

	const ccHObject::Container& selectedEntities = m_app->getSelectedEntities();
	size_t selNum = selectedEntities.size();
	if (selNum != 1)
	{
		m_app->dispToConsole("Select only one cloud!", ccMainAppInterface::ERR_CONSOLE_MESSAGE);
		return;
	}

	ccHObject* ent = selectedEntities[0];
	assert(ent);
	if (!ent || !ent->isA(CC_TYPES::POINT_CLOUD))
	{
		m_app->dispToConsole("Select a real point cloud!", ccMainAppInterface::ERR_CONSOLE_MESSAGE);
		return;
	}

	ccPointCloud* theCloud = static_cast<ccPointCloud*>(ent);
	if (!theCloud)
		return;
	unsigned count = theCloud->size();
	bool hasNorms = theCloud->hasNormals();
	CCVector3 bbMin, bbMax;
	theCloud->getBoundingBox(bbMin, bbMax);
	const CCVector3d& globalShift = theCloud->getGlobalShift();
	double globalScale = theCloud->getGlobalScale();

	CCVector3 diff = bbMax - bbMin;
	double scale = std::max(std::max(diff[0], diff[1]), diff[2]);
	//CloudCompare点云转换城PCL点云l
	pcl::PointCloud<PointXYZ>::Ptr pcl_t_cloud(new pcl::PointCloud<PointXYZ>);
	CC2PCL_PointCloud(*theCloud, *pcl_t_cloud);

	//输入参数，半径大小
	float kernelRadius = 1;
	int minPtsPerCell = 10;
	int minSupportCell = 2;
	double maxDistancetoPlane = 0.01;
	double cellDist2Plane = 0.001;

	ccNDTRansacSDDlg dlg;
	dlg.cellRaidiusSpinBox->setValue(.02f * scale);//注意此处设置计算法向量的参数
	dlg.cellepsilonDoubleSpinBox->setValue(.005f * scale);
	dlg.epsilonDoubleSpinBox->setValue(.001f * scale);		// set distance threshold to 0.5% of bounding box width

	if (!dlg.exec())
		return;
	kernelRadius = dlg.cellRaidiusSpinBox->value();
	minPtsPerCell = dlg.minPtsPerCellSpinBox->value();//每个NDT单元最小点数
	minSupportCell = dlg.supportPointsSpinBox->value();
	double teThreshold = dlg.tuThresholdSpinBox->value();//平面分割阈值
	cellDist2Plane = dlg.cellepsilonDoubleSpinBox->value();//Cell的距离阈值
	maxDistancetoPlane = dlg.epsilonDoubleSpinBox->value();//其他点到平面阈值

	double maxNormDevAngle = dlg.maxNormDevAngleSpinBox->value();
	double	probability_ = dlg.probaDoubleSpinBox->value();
	double angleThreshold = static_cast<float>(cos(maxNormDevAngle * CC_DEG_TO_RAD));
	int minPtsPerComponent = dlg.minPtsPerCCSpinBox->value();//连通分析每个连通区最小点数

	unsigned numberOfPoints = theCloud->size();
	if (numberOfPoints < 5)
	{
		m_app->dispToConsole("too little points!", ccMainAppInterface::ERR_CONSOLE_MESSAGE);
		return;
	}
	//进度条
	ccProgressDialog pDlg(true, (QWidget *)m_app->getMainWindow());
	QElapsedTimer eTimer;
	eTimer.start();
	ndtRansacHybridEx ndtrac;
	ndtrac.setInputCloud(theCloud);
	ndtrac.setLyrName(ent->getName());
	ndtrac.setKernelRadius(kernelRadius);
	ndtrac.setMinPointsPerCell(minPtsPerCell);
	ndtrac.setTeThreshold(teThreshold);
	ndtrac.setMinSupportCell(minSupportCell);
	ndtrac.setCellDist2Plane(cellDist2Plane);
	ndtrac.setMaxNormDevAngle(maxNormDevAngle);
	ndtrac.setMaxDistancetoPlane(maxDistancetoPlane);
	ndtrac.setProbability(probability_);
	ndtrac.setMinPointsPerComponent(minPtsPerComponent);
	ndtrac.extract(&pDlg);

	qint64 elapsedTime_ms = eTimer.elapsed();
	float timell = static_cast<float>(elapsedTime_ms) / 1.0e3;
	ccHObject* group = ndtrac.group;

	m_app->dispToConsole(QString("[Compute Time]: %2.3f s.").arg(timell), ccMainAppInterface::STD_CONSOLE_MESSAGE);
	//隐藏之前的点云
	theCloud->setEnabled(false);
	m_app->dispToConsole("[qNDTRansacSD] Input cloud has been automtically hidden!", ccMainAppInterface::STD_CONSOLE_MESSAGE);
	if (group)
	{
		//we add new group to DB/display
		group->setVisible(true);
		m_app->addToDB(group);
		m_app->refreshAll();
	}
}

void qOpenNDTPlugin::doNDTRegionGrowing()
{

}

void qOpenNDTPlugin::doTestAction()
{
	//m_app should have already been initialized by CC when plugin is loaded!
	//(--> pure internal check)
	assert(m_app);
	if (!m_app)
		return;

	/*** HERE STARTS THE ACTION ***/

	//put your code here
	//--> you may want to start by asking parameters (with a custom dialog, etc.)

	//This is how you can output messages
	m_app->dispToConsole("[qOpenNDTPlugin] Hello world!",ccMainAppInterface::STD_CONSOLE_MESSAGE); //a standard message is displayed in the console
	m_app->dispToConsole("[qOpenNDTPlugin] Warning: dummy plugin shouldn't be used as is!",ccMainAppInterface::WRN_CONSOLE_MESSAGE); //a warning message is displayed in the console
	m_app->dispToConsole("qOpenNDTPlugin plugin shouldn't be used as is!",ccMainAppInterface::ERR_CONSOLE_MESSAGE); //an error message is displayed in the console AND an error box will pop-up!

	/*** HERE ENDS THE ACTION ***/

}

QIcon qOpenNDTPlugin::getIcon() const
{
	return QIcon(":/CC/qOpenNDTPlugin/Resources/icon.png");
}
