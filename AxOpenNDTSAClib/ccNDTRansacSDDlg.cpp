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

#include "ccNDTRansacSDDlg.h"

static int    s_minSupport       = 2;		// this is the minimal numer of points required for a primitive
static double s_maxNormalDev_deg = 25.0;	// maximal normal deviation from ideal shape (in degrees)
static double s_probability      = 0.01;	// probability that no better candidate was overlooked during sampling
static double s_cellRadius      = 1;//��Ԫ��뾶
static double s_tuThreshold      = 0.01;//ndt ��ָ���ֵ
static int s_minPointsPerCell = 10;//ndt ��ָ���ֵ
static int s_minPointsPerComponent = 500;//ndt��ͨ����ÿ����ͨ����С����

ccNDTRansacSDDlg::ccNDTRansacSDDlg(QWidget* parent)
	: QDialog(parent, Qt::Tool)
	, Ui::NdtRansacDlg()
{
	setupUi(this);

	connect(buttonBox, SIGNAL(accepted()), this, SLOT(saveSettings()));

	supportPointsSpinBox->setValue(s_minSupport);
	maxNormDevAngleSpinBox->setValue(s_maxNormalDev_deg);
	probaDoubleSpinBox->setValue(s_probability);
	cellRaidiusSpinBox->setValue(s_cellRadius);
	tuThresholdSpinBox->setValue(s_tuThreshold);
	minPtsPerCCSpinBox->setValue(s_minPointsPerComponent);
	minPtsPerCellSpinBox->setValue(s_minPointsPerCell);
}

void ccNDTRansacSDDlg::saveSettings()
{
	s_cellRadius=cellRaidiusSpinBox->value();
	s_minPointsPerCell = minPtsPerCellSpinBox->value();
	s_tuThreshold=tuThresholdSpinBox->value();
	s_minSupport = supportPointsSpinBox->value();
	s_maxNormalDev_deg = maxNormDevAngleSpinBox->value();
	s_probability = probaDoubleSpinBox->value();
	s_minPointsPerComponent = minPtsPerCCSpinBox->value();
}