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
#include "planefit.h"

planefit::planefit(void)
	:plane_(NULL)
	,min_iterations (2)
	,max_iterations_ (1000)
	,threshold2stop(1e-5)
{
}


planefit::~planefit(void)
{

}

void planefit::setInputPointCloud(std::vector<pcl::PointXYZ> &input0_)
{
	input_=input0_;
}

void planefit::extract()
{
	if (input_.size()>0)
	{
		Eigen::Vector3d meanSum_;
		Eigen::Matrix3d covSum_;

		mean_<<0,0,0;
		for(unsigned int i=0; i< input_.size(); i++)
		{
			Eigen::Vector3d tmp;
			tmp<<input_[i].x,input_[i].y,input_[i].z;
			mean_ += tmp;
		}
		meanSum_ = mean_;
		mean_ /= (input_.size());
		Eigen::MatrixXd mp;
		mp.resize(input_.size(),3);
		for(unsigned int i=0; i< input_.size(); i++)
		{
			mp(i,0) = input_[i].x - mean_(0);
			mp(i,1) = input_[i].y - mean_(1);
			mp(i,2) = input_[i].z - mean_(2);
		}
		covSum_ = mp.transpose()*mp;
		cov_ = covSum_/(input_.size()-1);
		Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> Sol (cov_);

		evecs_ = Sol.eigenvectors().real();
		evals_ = Sol.eigenvalues().real();
		SortEigenValuesAndVectors(evecs_, evals_);
		//compute eigen value
		double e0 = evals_[0];
		double e1 = evals_[1];
		double e2 = evals_[2];

		Eigen::Vector3d n;
		//n[0]=evecs_(2,0);//其它的Cell
		//n[1]=evecs_(2,1);
		//n[2]=evecs_(2,2);
		n[0]=evecs_(0,2);//其它的Cell
		n[1]=evecs_(1,2);
		n[2]=evecs_(2,2);
		if (this->plane_==NULL)
		{
			this->plane_=new AxPlane();
		}
		this->plane_->nx=n[0];
		this->plane_->ny=n[1];
		this->plane_->nz=n[2];
		this->plane_->center_=mean_;
	}
}


void planefit::extractLmed()
{
	if (input_.size()>0)
	{
		unsigned num_of_points=input_.size();
		Eigen::Vector3d meanSum_;
		Eigen::Matrix3d covSum_;
		
		mean_<<0,0,0;
		for(unsigned int i=0; i< input_.size(); i++)
		{
			Eigen::Vector3d tmp;
			tmp<<input_[i].x,input_[i].y,input_[i].z;
			mean_ += tmp;
		}
		meanSum_ = mean_;
		mean_ /= (input_.size());
		Eigen::MatrixXd mp;
		mp.resize(input_.size(),3);
		for(unsigned int i=0; i< input_.size(); i++)
		{
			mp(i,0) = input_[i].x - mean_(0);
			mp(i,1) = input_[i].y - mean_(1);
			mp(i,2) = input_[i].z - mean_(2);
		}
		covSum_ = mp.transpose()*mp;
		cov_ = covSum_/(input_.size()-1);
		Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> Sol (cov_);

		evecs_ = Sol.eigenvectors().real();
		evals_ = Sol.eigenvalues().real();
		SortEigenValuesAndVectors(evecs_, evals_);
		//compute eigen value
		double e0 = evals_[0];
		double e1 = evals_[1];
		double e2 = evals_[2];

		Eigen::Vector3d normal;
		normal[0]=evecs_(0,2);//其它的Cell
		normal[1]=evecs_(1,2);
		normal[2]=evecs_(2,2);
		if (this->plane_==NULL)
		{
			this->plane_=new AxPlane();
		}
		this->plane_->nx=normal[0];
		this->plane_->ny=normal[1];
		this->plane_->nz=normal[2];
		this->plane_->center_=mean_;

		//进行Robust平面获取
		Eigen::Vector3d oldnormal=normal;
		//IRLS主循环
		for (int iter=0;iter<max_iterations_;iter++)
		{
			oldnormal=normal;
			float  sum_dist=0;
			//循环1，求距离均值
			for (unsigned i=0; i<num_of_points; ++i)
			{
				const pcl::PointXYZ CellP=input_[i];
				Eigen::Vector3d pt(CellP.x,CellP.y,CellP.z);
				float distance=fabs(normal.dot(pt-mean_));//点到平面的距离
				sum_dist+=distance;
			}
			//距离均值
			double mean_dist_=sum_dist/num_of_points;
			double s_dist_mean=0;
			//循环2.求距离方差
			for (unsigned i=0; i<num_of_points; ++i)
			{
				const pcl::PointXYZ CellP=input_[i];
				Eigen::Vector3d pt(CellP.x,CellP.y,CellP.z);
				float distance=fabs(normal.dot(pt-mean_));//点到平面的距离
				s_dist_mean+=(distance-mean_dist_)*(distance-mean_dist_);
			}
			//距离方差
			double sigma=sqrt(s_dist_mean/(num_of_points-1));
			if (sigma<0.000001)
			{
				break;
			}
			//循环3.过滤满足小于2*sigma的点，根据这些点求新的均值
			//求距离均值
			int num_of_ok_points=0;
			Eigen::Vector3d mean_new;
			Eigen::Matrix3d covSum_new;

			mean_new<<0,0,0;
			for (unsigned i=0; i<num_of_points; ++i)
			{
				const pcl::PointXYZ CellP=input_[i];
				Eigen::Vector3d pt(CellP.x,CellP.y,CellP.z);
				float distance=fabs(normal.dot(pt-mean_));//点到平面的距离
				if (distance-mean_dist_<sigma*2)
				{
					//记录该点索引
					num_of_ok_points++;
					mean_new+=pt;
				}
			}
			mean_new/=num_of_ok_points;
			//根据标记的点重新计算均值和方差，产生新的均值点和法向量
			//计算新的方差
			Eigen::MatrixXd mp_1;
			mp_1.resize(num_of_ok_points,3);
			int indx_of_ok_point=0;
			for (unsigned i=0; i<num_of_points; ++i)
			{
				const pcl::PointXYZ CellP=input_[i];
				Eigen::Vector3d pt(CellP.x,CellP.y,CellP.z);
				float distance=fabs(normal.dot(pt-mean_));//点到平面的距离
				if (distance-mean_dist_<sigma*2)
				{
					//记录该点索引
					mp_1(indx_of_ok_point,0) = input_[i].x - mean_new(0);
					mp_1(indx_of_ok_point,1) = input_[i].y - mean_new(1);
					mp_1(indx_of_ok_point,2) = input_[i].z - mean_new(2);
					indx_of_ok_point++;
				}
			}
			covSum_new = mp_1.transpose()*mp_1;
			cov_ = covSum_new/(num_of_ok_points-1);
			Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> Sol_1 (cov_);

			evecs_ = Sol_1.eigenvectors().real();
			evals_ = Sol_1.eigenvalues().real();
			SortEigenValuesAndVectors(evecs_, evals_);
			//compute eigen value
			double e00 = evals_[0];
			double e01 = evals_[1];
			double e02 = evals_[2];

			normal[0]=evecs_(0,2);
			normal[1]=evecs_(1,2);
			normal[2]=evecs_(2,2);
			mean_= mean_new;
			//判断是否终止循环
			Eigen::Vector3d temp=oldnormal-normal;
			double* con=new double[3];
			con[0]=std::abs(temp[0])/std::abs(oldnormal[0]);
			con[1]=std::abs(temp[1])/std::abs(oldnormal[1]);
			con[2]=std::abs(temp[2])/std::abs(oldnormal[2]);
			std::sort(con,con+3);
			double	convg = con[2];
			if (iter>=min_iterations)
			{
				if (convg < threshold2stop)
				{
					break;
				}
			}	
		}
		//结束循环，获取了优化的均值点和法向量
		this->plane_->nx=normal[0];
		this->plane_->ny=normal[1];
		this->plane_->nz=normal[2];
		this->plane_->center_=mean_;
		this->plane_->evals=evals_;
		this->plane_->evecs=evecs_;
	}
}


AxPlane planefit::getPlane()
{
	this->plane_->points_=input_;
	return *this->plane_;
}
