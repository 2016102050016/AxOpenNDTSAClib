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

#ifndef CC_NDTRANSAC_SD_DLG_HEADER
#define CC_NDTRANSAC_SD_DLG_HEADER

#include "ui_NdtRansacDlg.h"

//! Dialog for qRansacSD plugin
class ccNDTRansacSDDlg : public QDialog, public Ui::NdtRansacDlg
{
	Q_OBJECT

public:

	//! Default constructor
	explicit ccNDTRansacSDDlg(QWidget* parent = 0);

protected slots:

	//! Saves (temporarily) the dialog paramters on acceptation
	void saveSettings();

};

#endif
