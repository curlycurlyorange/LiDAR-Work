#ifndef __LIDARBASETOOLS__
#define __LIDARBASETOOLS__

#include <iostream>
#include <deque>
#include <vector>
#include <string>
#include <fstream>

#include <iomanip>


namespace LiDARBaseTools{

	const float NODATA = -999.99;
	const float eps = 0.001; 

	struct LasPoint{
		float x;
		float y;
		float z;
		short classification;
		short intensity;
	};

	struct GrideIndex//格网索引
	{
		std::vector<int> ptIndex;
	};
	
	bool readTXTPointData(const char* filename,std::vector<LasPoint>& inputCloud)
	/*
	todo:      read txt point data to a LasPoint vector
	date:      2015.3.8
	author:    JianPing(lijianping@whu.edu.cn)
	see:
	*/	
	{
		std::cout<<"begin reading..."<<std::endl;

		std::fstream infile;
		infile.open(filename);

		if (infile.fail())
		{
			std::cout<<"bad file name"<<std::endl;
			return false;
		}
		
		int nCount(0);
		char temp[256];
		LasPoint tempLasPoint;
		float x(0),y(0),z(0);
		

		while (!infile.eof())
		{
			nCount++;
			infile.getline(temp,256);
			//std::cout<<temp<<std::endl;
			std::sscanf(temp,"%f %f %f",&x,&y,&z);
			tempLasPoint.x = x;
			tempLasPoint.y = y;
			tempLasPoint.z = z;
			inputCloud.push_back(tempLasPoint);
		}

		infile.close();

		std::cout<<"have read "<<nCount<<" points.\nDone."<<std::endl;

		return true;
	}

	bool writeTXTPointData(const char* filename,const std::vector<LasPoint> &outputCloud)
	/*
	todo:      write point data of txt formate from a point cloud vector
	out:       
	date:      2015.3.9
	author:    JianPing(lijianping@whu.edu.cn)
	see:
	*/
	{
		std::cout<<"begin to write..."<<std::endl;

		std::ofstream outfile;
		outfile.open(filename);

		if (outfile.fail())
		{
			std::cout<<"bad filename"<<std::endl;
			return false;
		}

		for (std::vector<LasPoint>::const_iterator it = outputCloud.begin();it!=outputCloud.end();++it)
		{
			//outfile.setprecision(2);
			//std::setprecision(2);
			outfile<<std::setprecision(2)<<std::fixed<<it->x<<" "<<it->y<<" "<<it->z<<"\n";
		}

		outfile.close();

		std::cout<<"have write "<<outputCloud.size()<<" points.\nDone."<<std::endl;

		return true;
	}

	bool writeTXTPointData(const char* filename,const std::vector<LasPoint> &outputCloud , const std::vector<int> &index)
	{
		std::cout<<"begin to write..."<<std::endl;

		std::ofstream outfile;
		outfile.open(filename);

		if (outfile.fail())
		{
			std::cout<<"bad filename"<<std::endl;
			return false;
		}

		const LasPoint *ptempLasPoint;
		int count = 0;

		for (std::vector<int>::const_iterator it = index.begin() ; it != index.end() ; it++ )
		{
			ptempLasPoint = &outputCloud[*it];
			outfile<<std::setprecision(2)<<std::fixed<<ptempLasPoint->x<<" "<<ptempLasPoint->y<<" "<<ptempLasPoint->z<<"\n";
			count++;
		}

		outfile.close();

		std::cout<<"have write "<<count<<" points.\nDone."<<std::endl;

		return true;
	}

	bool createGride(const std::vector<LasPoint> &inputCloud,
		int &ncols,int &nrows,
		float &maxx,float &minx,float &maxy,float &miny,float &maxz,float &minz,
		const float cellsize,
		GrideIndex* &pGrideIndex)
	/*
	todo:      create a gride of a point cloud
	out:       size and range of the gride
	date:      2015.3.8
	author:    JianPing(lijianping@whu.edu.cn)
	see:
	*/
	{
		std::cout<<"begin create gride..."<<std::endl;
		maxx = inputCloud[0].x;
		maxy = inputCloud[0].y;
		maxz = inputCloud[0].z;
		minx = inputCloud[0].x;
		miny = inputCloud[0].y;
		minz = inputCloud[0].z;
		for ( std::vector<LasPoint>::const_iterator it = inputCloud.begin(); it != inputCloud.end(); ++it)
		{
			if (it->x > maxx)
			{
				maxx = it->x;
			}
			if (it->x < minx)
			{
				minx = it->x;
			}
			if (it->y > maxy)
			{
				maxy = it->y;
			}
			if (it->y < miny)
			{
				miny = it->y;
			}
			if (it->z > maxz)
			{
				maxz = it->z;
			}
			if (it->z < minz)
			{
				minz = it->z;
			}
		}

		ncols = (maxx - minx)/cellsize + 1;
		nrows = (maxy - miny)/cellsize + 1;

		pGrideIndex = new GrideIndex[ncols*nrows];

		int nCount = 0;
		for ( std::vector<LasPoint>::const_iterator it = inputCloud.begin(); it != inputCloud.end(); ++it)
		{
			int tempCol,tempRow;
			tempCol = (it->x - minx)/cellsize;
			tempRow = (it->y - miny)/cellsize;
			pGrideIndex[tempCol + tempRow*ncols].ptIndex.push_back(nCount);
			nCount++;
		}

		std::cout<<"have created a "<<nrows<<" * "<<ncols<<" gride.\nDone."<<std::endl;

		return true;
	}

	bool raster2asc(const char* filename,
		const float *rasterData, 
		const int nrows, const int ncols,
		const float xllcorner,const float yllcorner,
		const float cellsize)
	/*
	todo:      write raster data to a asc file
	date:      2015.3.9
	author:    JianPing(lijianping@whu.edu.cn)
	see:
	*/	
	{
		std::cout<<"begin to write asc file..."<<std::endl;

		std::ofstream outfile;
		outfile.open(filename);

		if (outfile.fail())
		{
			std::cout<<"bad filename"<<std::endl;
		}

		outfile<<"ncols "<<ncols<<"\n";
		outfile<<"nrows "<<nrows<<"\n";
		outfile<<"xllcorner "<<xllcorner<<"\n";
		outfile<<"yllcorner "<<yllcorner<<"\n";
		outfile<<"cellsize "<<cellsize<<"\n";
		outfile<<"NODATA_value "<<-9999<<"\n";

		for (int i = nrows-1 ; i >= 0 ; i-- )
		{
			for (int j = 0 ; j < ncols ; j++ )
			{
				outfile<<std::setprecision(2)<<std::fixed<<rasterData[i*ncols+j]<<" ";
			}
			outfile<<"\n";
		}

		std::cout<<"have written asc file.\nDone."<<std::endl;

		return true;

	}

	float getMinZvalueofRasterCell(const LiDARBaseTools::GrideIndex* pGrideIndex,
		const int nrows,const int ncols,
		const std::vector<LiDARBaseTools::LasPoint> &inputCloud,
		const int i ,const int j)
	/*
	todo:      get mini Z value of points in a Raster cell
	date:      2015.3.9
	author:    JianPing(lijianping@whu.edu.cn)
	see:       
	*/	
	{
		float minZ = inputCloud[pGrideIndex[i*ncols+j].ptIndex[0]].z;
		float tempZ;
		for (std::vector<int>::const_iterator it = pGrideIndex[i*ncols+j].ptIndex.begin();
			it!=pGrideIndex[i*ncols+j].ptIndex.end();++it)
		{
			tempZ = inputCloud[*it].z;
			if (minZ>tempZ)
			{
				minZ = tempZ;
			}
		}
		return minZ;
	}


	bool CreateBreakCluster(const float *rasterData, const int nrows, const int ncols,
		std::vector<std::vector<int>> &cluster_result)
	/*
	todo:      cluster nodata areas in a raster
	date:      2015.3.10
	author:    JianPing(lijianping@whu.edu.cn)
	see:
	*/
	{
		bool* breakMask = new bool[nrows*ncols];
		for (int i=0;i<nrows*ncols;i++)
		{
			breakMask[i] = false;
		}

		//找到无数据格网
		for (int i=0;i<nrows;i++)
		{
			for (int j=0;j<ncols;j++)
			{
				if (std::fabs((rasterData[i*ncols+j] - NODATA)) < eps )
				{
					breakMask[i*ncols+j] = true;
				}
			}
		}

		//对无数据网格进行聚类
		bool* breakMaskCopy = new bool[nrows*ncols];
		for (int i=0;i<nrows*ncols;i++ )
		{
			breakMaskCopy[i] = false;
		}

		std::deque<int> seed_x;//聚类种子
		std::deque<int> seed_y;//聚类种子

		//收集无数据格网种子
		for (int i=0;i<nrows;i++)
		{
			for (int j=0;j<ncols;j++)
			{
				if (breakMask[i*ncols+j])
				{
					seed_x.push_back(j);
					seed_y.push_back(i);
				}
			}
		}

		std::vector<int> tempCluster;
		std::deque<int> tempseed_x;
		std::deque<int> tempseed_y;//当前聚类

		while (seed_x.size()>0)
		{
			int x_num,y_num;
			x_num = seed_x.front();

			y_num = seed_y.front();

			//判断当前网格是否已经被聚类
			if (breakMaskCopy[y_num*ncols+x_num])
			{
				seed_x.pop_front();
				seed_y.pop_front();
				continue;
			}
			else
			{
				seed_x.pop_front();
				seed_y.pop_front();
			}

			tempseed_x.push_back(x_num);
			tempseed_y.push_back(y_num);
			tempCluster.push_back(y_num*ncols+x_num);
			breakMaskCopy[y_num*ncols+x_num] = true;//设置当前点已经被聚类

			while (tempseed_x.size()>0)
			{
				x_num = tempseed_x.front();
				tempseed_x.pop_front();

				y_num = tempseed_y.front();
				tempseed_y.pop_front();

				for (int i=y_num-1;i<=y_num+1;i++)//8临域
				{
					if (i<0||i>=nrows)
						continue;
					for (int j=x_num-1;j<=x_num+1;j++)
					{
						if(j<0||j>=ncols)
							continue;
						if (breakMaskCopy[i*ncols+j])
							continue;
						if (breakMask[i*ncols+j])
						{
							tempseed_x.push_back(j);
							tempseed_y.push_back(i);
							tempCluster.push_back(i*ncols+j);
							breakMaskCopy[i*ncols+j] = true;
						}
					}
				}
			}
			cluster_result.push_back(tempCluster);
			tempCluster.clear();
		}
		
		seed_x.clear(); seed_y.clear();
		delete[] breakMaskCopy;
		delete[] breakMask;

		return true;
	}

	bool Ruster(std::vector<std::vector<int>> &cluster_raster,
		std::vector<std::vector<int>> &rust_raster,
		std::vector<std::vector<int>> &boundary_raster,
		const int nrows,const int ncols)
	{
		bool* RustMask = new bool[ncols*nrows];//构建掩膜，记录待腐蚀点
		std::vector<int> temp_cluster;//当前待腐蚀
		for(int i = 0;i<ncols*nrows;i++)
		{
			RustMask[i] = false;
		}

		for(int i= 0;i<cluster_raster.size();i++)//找到待腐蚀点，记录在掩膜中
		{
			temp_cluster = cluster_raster[i];
			for(int j = 0;j<temp_cluster.size();j++)
			{
				RustMask[temp_cluster[j]] = true;
			}
			temp_cluster.clear();
		}


		std::vector<int> temp_rust_raster;//当前腐蚀结果
		std::vector<int> temp_boundary_raster;//当前腐蚀边界
		for(int i= 0;i<cluster_raster.size();i++)//腐蚀运算
		{
			temp_cluster = cluster_raster[i];//当前待腐蚀
			for(int j = 0;j<temp_cluster.size();j++)
			{
				int num_row,num_col;
				num_row = temp_cluster[j]/ncols;
				num_col = temp_cluster[j]%ncols;

				bool flag = true;
				for(int ii = num_row-1;ii<=num_row+1;ii++)//八临域
				{
					if(ii<0||ii>=nrows)
						continue;
					for(int jj = num_col-1;jj<=num_col+1;jj++)
					{
						if(jj<0||jj>=ncols)
							continue;
						if(!RustMask[ii*ncols+jj])
						{
							flag = false;
							break;
						}
					}
					if(!flag)
					{
						break;
					}
				}
				if(flag)
				{
					temp_rust_raster.push_back(num_row*ncols+num_col);
				}
				else
				{
					temp_boundary_raster.push_back(num_row*ncols+num_col);
				}
			}
			if(temp_rust_raster.size())
				rust_raster.push_back(temp_rust_raster);
			if(temp_boundary_raster.size())	
				boundary_raster.push_back(temp_boundary_raster);

			temp_cluster.clear();
			temp_rust_raster.clear();
			temp_boundary_raster.clear();
		}
		delete[] RustMask;
		return true;
	}

	inline float getAverageValueof3x3RasterWindow(const float *rasterData,const int nrows,const int ncols,
		const int row, const int col)
	/*
	todo:      transform grideIndex to a raster
	date:      2015.3.10
	author:    JianPing(lijianping@whu.edu.cn)
	see:
	*/
	{
		float sum = 0;
		int count = 0;
		for (int i = row-1 ; i<=row+1 ; i++)
		{
			if (i<0 || i>=nrows)
			{
				continue;
			}
			for (int j = col-1 ; j<= col+1 ; j++)
			{
				if (j<0 || j>=ncols)
				{
					continue;
				}
				if (std::fabs(rasterData[j+i*ncols]-NODATA)<eps)
				{
					continue;
				}
				sum += rasterData[j+i*ncols];
				count++;
			}
		}
		return sum/count;
	}

	bool grideIndex2raster(const std::vector<LasPoint> &inputCloud,
		const GrideIndex* pGrideIndex,
		const int nrows,const int ncols,
		float *&rasterData)
	/*
	todo:      transform grideIndex to a raster
	date:      2015.3.10
	author:    JianPing(lijianping@whu.edu.cn)
	see:
	*/
	{
		std::cout<<"begin to transform gride to a raster..."<<std::endl;
		rasterData = new float[nrows*ncols];
		for (int i = 0 ; i< ncols*nrows ; ++i)
		{
			rasterData[i] = NODATA;//init
		}

		const std::vector<int> *ptempIndex;
		for (int i = 0 ; i < nrows ; ++i)//init the raster with the lowest value of the cell
		{
			for (int j = 0 ; j < ncols ; ++j)
			{
				ptempIndex = &(pGrideIndex[j+i*ncols].ptIndex);
				if(0 == ptempIndex->size())
				{
					rasterData[j+i*ncols] = NODATA;
				}
				else
				{
					rasterData[j+i*ncols] = getMinZvalueofRasterCell(pGrideIndex,
						nrows,ncols,inputCloud,i,j);
				}
			}
		}

		//cluster nodata value
		std::vector<std::vector<int>> cluster_result;//restore the index set of nodata value in the raw raster
		CreateBreakCluster(rasterData,nrows,ncols,cluster_result);

		//Ruster
		std::vector<std::vector<int>> rust_begin,rust_result,rust_boundary;
		
		rust_result = cluster_result;
		while(rust_result.size())
		{
			rust_begin = rust_result;
			rust_result.clear();
			rust_boundary.clear();
			Ruster(rust_begin,rust_result,rust_boundary,nrows,ncols);

			int row,col;

			for (std::vector<std::vector<int>>::iterator it_i = rust_boundary.begin() ; it_i!=rust_boundary.end();++it_i)
			{
				for (std::vector<int>::iterator it_j = it_i->begin(); it_j!=it_i->end();++it_j)
				{
					row = (*it_j)/ncols;
					col = (*it_j)-row*ncols;
					rasterData[*it_j] = getAverageValueof3x3RasterWindow(rasterData,nrows,ncols,row,col);
				}
			}

			/*for(int i = 0;i<rust_boundary.size();i++)
			{
				for(int j = 0;j<rust_boundary[i].size();j++)
				{
					IDW_Interpolate_break(rust_boundary[i][j],this->radius);
				}
			}*/
		}
		std::cout<<"transform raster success!\nDone."<<std::endl;
		return true;
	}

	bool lidarRaster2demRaster(const float *raster,const int nrows,const int ncols,
		const int *flag,
		float *&demraster)
	/*
	todo:      transform the lidar measurements to dem raster according to the flag(0-nonground points)
	date:      2015.3.12
	author:    JianPing(lijianping@whu.edu.cn)
	see:
	*/
	{
		std::cout<<"begin to create dem..."<<std::endl;
		demraster = new float[nrows*ncols];
		for (int i = 0 ; i < nrows*ncols ; i++)
		{
			if (flag[i]>0)
			{
				demraster[i] = NODATA;
			}
			else
			{
				demraster[i] = raster[i];
			}
		}
		//cluster nodata value
		std::vector<std::vector<int>> cluster_result;//restore the index set of nodata value in the raw raster
		CreateBreakCluster(demraster,nrows,ncols,cluster_result);

		//Ruster
		std::vector<std::vector<int>> rust_begin,rust_result,rust_boundary;

		rust_result = cluster_result;
		while(rust_result.size())
		{
			rust_begin = rust_result;
			rust_result.clear();
			rust_boundary.clear();
			Ruster(rust_begin,rust_result,rust_boundary,nrows,ncols);

			int row,col;

			for (std::vector<std::vector<int>>::iterator it_i = rust_boundary.begin() ; it_i!=rust_boundary.end();++it_i)
			{
				for (std::vector<int>::iterator it_j = it_i->begin(); it_j!=it_i->end();++it_j)
				{
					row = (*it_j)/ncols;
					col = (*it_j)-row*ncols;
					demraster[*it_j] = getAverageValueof3x3RasterWindow(demraster,nrows,ncols,row,col);
				}
			}
		}
		std::cout<<"have created dem.\n"<<std::endl;
		return true;
	}

	

		
}

#endif
