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

#ifndef Q_OPENNDTSAC_PLUGIN_HEADER
#define Q_OPENNDTSAC_PLUGIN_HEADER
#include <QtGui>
#include <qmessagebox.h>
#include <qfiledialog.h>
#include <ccPointCloud.h>
#include "ccStdPluginInterface.h"
#include <pcl/io/pcd_io.h>
#include <pcl/point_cloud.h>
#include <pcl/point_traits.h>
#include <pcl/common/transforms.h>
#include "ccNDTRansacSDDlg.h"
#include "ndtRansacHybridEx.h"
#include "dataUtility.h"

using namespace pcl;
using namespace std;
//! MIMS plugin

class qOpenNDTPlugin : public QObject, public ccStdPluginInterface
{
	Q_OBJECT
	Q_INTERFACES(ccStdPluginInterface)
	//replace qDummy by the plugin name (IID should be unique - let's hope your plugin name is unique ;)
	Q_PLUGIN_METADATA(IID "cccorp.cloudcompare.plugin.qOpenNDTPlugin")

public:

	//! Default constructor
	explicit qOpenNDTPlugin(QObject* parent = 0);

	//inherited from ccPluginInterface
	virtual QString getName() const override { return "qOpenNDTPlugin"; }
	virtual QString getDescription() const override { return "qOpenNDTPlugin plugin (add description here)"; }
	virtual QIcon getIcon() const override;

	//inherited from ccStdPluginInterface
	void onNewSelection(const ccHObject::Container& selectedEntities) override;
	virtual void getActions(QActionGroup& group) override;

protected slots:
	/*** ADD YOUR CUSTOM ACTIONS' SLOTS HERE ***/
	void doTestAction();
	void doNDTRansac();
	void doNDTRegionGrowing();
public:
	void display(pcl::PointCloud<pcl::PointXYZ> &ptcl, QString group_name);
	void display(pcl::PointCloud<pcl::PointXYZRGB> &ptcl, QString group_name);
protected:
	QAction* m_action_test;
	QAction* m_action_ndtRansac;
	QAction* m_action_ndtRegionGrowing;
};

#endif
