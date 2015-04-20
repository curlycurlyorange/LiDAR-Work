#ifndef __MORPHOLOGICALFILTER__
#define __MORPHOLOGICALFILTER__

#include "BaseTools.h"
#include <math.h>
namespace MorphologicalFilter
{
	bool Erosion(const float* rasterData,const int nrows,const int ncols,
		float* &eroResultData, const int windowSize)
	/*
	todo:      erosion a raster file with a 2-D window 
	date:      2015.3.10
	author:    JianPing(lijianping@whu.edu.cn)
	see:       A Progressive Morphological Filter for Removing
			   Nonground Measurements From Airborne LiDAR Data(2003)
	*/
	{
		int r = windowSize/2;
		eroResultData = new float[nrows*ncols];
		float tempMin;
		float tempValue;
		for (int row = 0 ; row < nrows ; row++)
		{
			for (int col = 0 ; col < ncols ; col++)
			{
				tempMin = rasterData[col+row*ncols];
				for (int i = row-r ; i <= row+r ; i++)
				{
					if (i<0||i>=nrows)
					{
						continue;
					}
					for (int j = col-r ; j <= col+r ; j++)
					{
						if (j<0||j>=ncols)
						{
							continue;
						}
						tempValue = rasterData[j+i*ncols];
						if (tempValue<tempMin)
						{
							tempMin = tempValue;
						}
					}
				}
				eroResultData[col+row*ncols] = tempMin;
			}
		}
		return true;
	}

	bool Dilation(const float* rasterData,const int nrows,const int ncols,
		float* &dilatResultData, const int windowSize)
	/*
	todo:      dilate a raster file with a 2-D window 
	date:      2015.3.10
	author:    JianPing(lijianping@whu.edu.cn)
	see:       A Progressive Morphological Filter for Removing
			   Nonground Measurements From Airborne LiDAR Data(2003)
	*/
	{
		int r = windowSize/2;
		dilatResultData = new float[nrows*ncols];
		float tempMax;
		float tempValue;
		for (int row = 0 ; row<nrows ; row++)
		{
			for (int col = 0 ; col < ncols ; col++)
			{
				tempMax = rasterData[col+row*ncols];
				for (int i = row-r ; i<=row+r ; i++)
				{
					if (i<0||i>=nrows)
					{
						continue;
					}
					for (int j = col-r ; j<=col+r ; j++)
					{
						if (j<0||j>=ncols)
						{
							continue;
						}
						tempValue = rasterData[j+i*ncols];
						if (tempValue>tempMax)
						{
							tempMax = tempValue;
						}
					}
				}
				dilatResultData[col+row*ncols] = tempMax;
			}
		}
		return true;
	}



	/************************************************************************/
	/*=================================main=================================*/
	/************************************************************************/
	bool ProgessiveMorphologicalFilter(
		const float* rasterData,
		const int nrows,const int ncols,
		const float cellSize,
		const float b,
		const float maxWinSize,
		const float terrSlope,
		const float dh0,
		const float dhmax,
		int* &flag)
	/*
	todo:      the main function
	in:        lidar raster data
	           cell size
			   Parameter b (init window size)
			   Maximum window size
			   Terrain slope s
			   Initial elevation difference threshold dh0
			   Maximum elevation difference threshold dhmax
	out:       flag = 0 ground points index 
	           nonground points index
	date:      2015.3.10
	author:    JianPing(lijianping@whu.edu.cn)
	see:       A Progressive Morphological Filter for Removing
	           Nonground Measurements From Airborne LiDAR Data(2003)
	*/
	{
		std::cout<<"begin to MorphologicalFilter..."<<std::endl;
		//1.create a gride
		/*LiDARBaseTools::GrideIndex *pGrideIndex = NULL;
		int ncols,nrows;
		float maxx,minx,maxy,miny,maxz,minz;
		LiDARBaseTools::createGride(inputCloud,ncols,nrows,
			maxx,minx,maxy,miny,maxz,minz,cellSize,pGrideIndex);*/

		/*float* rasterData = NULL;
		LiDARBaseTools::grideIndex2raster(inputCloud,pGrideIndex,nrows,ncols,rasterData);*/

		//2.determine series of window size
		float dhT;//threshold (differences between terrian an building )
		flag = new int[nrows*ncols];
		for (int i = 0 ; i < ncols*nrows ; i++)
		{
			flag[i] = 0;
		}

		float *lasterRaster = new float[nrows*ncols];
		float *tempRaster1 = NULL;
		float *tempRaster2 = NULL;

		for (int i = 0 ; i<nrows*ncols ; i++)
		{
			lasterRaster[i] = rasterData[i];
		}

		int tempWindowSize = b;
		for (int k = 1; k <= 9 ; k++)//each window size
		{
			
			if (tempWindowSize <= 3)
			{
				dhT = dh0;
			}
			else
			{
				dhT = terrSlope*2*b*k*cellSize+dh0;
				if (dhT > dhmax)
				{
					dhT = dhmax;
				}
			}

			//dhT = 2;

			Erosion(lasterRaster,nrows,ncols,tempRaster1,tempWindowSize);
			Dilation(tempRaster1,nrows,ncols,tempRaster2,tempWindowSize);
			delete[] tempRaster1;
			tempRaster1 = NULL;

			for (int i = 0 ; i < nrows ; i++)
			{
				for (int j = 0 ; j < ncols ; j++)
				{
					if (lasterRaster[j+i*ncols]-tempRaster2[j+i*ncols]>dhT)
					{
						if (flag[j+i*ncols] > 0)
						{
							continue;
						}
						flag[j+i*ncols] = tempWindowSize;
					}
				}
			}
			delete[] lasterRaster;
			lasterRaster = tempRaster2;
			tempRaster2 = NULL;
			tempWindowSize = 2*pow(b,k)+1;
			if (tempWindowSize > maxWinSize)
			{
				break;
			}
		}
		
		//delete []rasterData;
		//delete []pGrideIndex;
		std::cout<<"MorphologicalFilter Done."<<std::endl;
		return true;
	}
}

#endif