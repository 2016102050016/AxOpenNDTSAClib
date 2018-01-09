#ifndef CC_CELLDEPICT
#define CC_CELLDEPICT
//system
#include <algorithm>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <set>
#include <vector>
#include <assert.h>
#include <string.h>
#include "SquareMatrix.h"
#include "DgmOctreeNDT.h"
using namespace CCLib;
struct CellDepict
{
	bool isPlanar;
	//! Cell center
	CCVector3 center;
	CCVector3 normal;
	CCLib::SquareMatrixd cov_;		/// Contains the covatiance of the normal distribution
	CCLib::SquareMatrixd icov_;  /// Precomputed inverse covariance (updated every time the cell is updated)
	CCLib::SquareMatrixd evecs_; /// Eigen vectors

	CCVector3 mean_;  /// Mean of the normal distribution
	CCVector3 evals_; /// Eigen values

	float sigma;//�㵽���ƽ�����ı�׼��

	CCNDTLib::DgmOctreeNDT::OctreeCellCodeType Code;
	//! First point index in associated NeighboursSet
	unsigned index;
	unsigned cellIndex;
	//�����ĵ㼯��
	std::vector<CCVector3> points;
	//���������������
	std::vector<unsigned> pointIncides;
	//! Default empty constructor
	CellDepict() {}

	//! Constructor from a point and an index
	CellDepict(const CCVector3& C, unsigned i)
		: center(C)	, index(i)
	{}
};
#endif