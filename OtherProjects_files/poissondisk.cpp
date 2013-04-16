
/*
 * pbrt source code Copyright(c) 1998-2007 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

// poissondisk.cpp*
#include "sampling.h"
#include "paramset.h"
#include "film.h"

#include <fstream>

// PoissonDiskSampler Declarations
class PoissonDiskSampler : public Sampler {
public:
	// PoissonDiskSampler Public Methods
	PoissonDiskSampler(int xstart, int xend,
	          int ystart, int yend,
			  int nsamp);
	~PoissonDiskSampler() {
		delete[] imageSamples;
		for (int i = 0; i < n1D; ++i)
			delete[] oneDSamples[i];
		for (int i = 0; i < n2D; ++i)
			delete[] twoDSamples[i];
		delete[] oneDSamples;
		delete[] twoDSamples;

		delete []downSamples;
		delete []right;
		delete []down;
		delete []activeGrid;
		delete []activeList;
	}
	int RoundSize(int size) const {
		return size;
	}
	bool GetNextSample(Sample *sample);
private:
	// PoissonDiskSampler Private Data
	int xPos, yPos, pixelSamples;
	int samplePos;
	float *imageSamples, *lensSamples, *timeSamples;
	float **oneDSamples, **twoDSamples;
	int n1D, n2D;

	float radius2D, radius1D;
	float gridSize1D, gridSize2D;
	int *activeGrid;
	int *right;
	int *down;
	int *downBackup;
	int downIndex;
	int nX, nY;
	int maxAttemp;
	int* activeList;
	float *rightSamples;
	float *downSamples;
	float *downSamplesBackup;

	int genPattern1D(float* samples);
	int genPattern2D(float* samples, bool isImage=false);
};
// PoissonDiskSampler Method Definitions
PoissonDiskSampler::PoissonDiskSampler(int xstart, int xend,
		int ystart, int yend, int ps)
	: Sampler(xstart, xend, ystart, yend, ps) {
	std::cout<<"xstart: "<<xstart<<"xend: "<<xend<<std::endl;
	xPos = xPixelStart-1;
	yPos = yPixelStart;
	pixelSamples = ps;
	samplePos = 0;
	oneDSamples = twoDSamples = NULL;
	n1D = n2D = 0;
	maxAttemp = 100;

	// expected distance between 2 points is 1.5*radius, where radius is the min distance between two points
	radius1D = 1.0/(1.5*(pixelSamples));	
	radius2D = 1.0/sqrt(0.625*(pixelSamples)*M_PI);
	gridSize1D = radius1D;
	gridSize2D = radius2D/sqrt(2.0);
	nX = nY = Ceil2Int(1.0/gridSize2D); 
	activeGrid = new int[nX*nY];
	activeList = new int[nX*nY];
	right = new int[2*nY];
	down = new int[4*nX*(xend - xstart+1)];
	downBackup = down + 2*nX*(xend-xstart+1);
	downSamples = new float[ 4*nY + 8*nX*(xend-xstart+1) ];
	downSamplesBackup = downSamples + 4*nX*(xend-xstart+1);
	rightSamples = downSamples + 8*nX*(xend-xstart+1);
	memset(right, -1, 2*nY*sizeof(int));
	memset(downBackup, -1, 2*nX*(xend-xstart+1)*sizeof(int));
	downIndex = 0;

	// num of grids is the upper bound for pixel samples
	imageSamples = new float[5*nX*nY];
	lensSamples = imageSamples + 2*nX*nY;
	timeSamples = imageSamples + 4*nX*nY;

	pixelSamples = 0;
	std::cout<<"num of samples per pixel: "<<ps<<std::endl
		<<"radius1D: "<<radius1D<<"; radius2D: "<<radius2D<<std::endl
		<<"gridSize1D: "<<gridSize1D<<"; gridSize2D: "<<gridSize2D<<std::endl;
}

bool PoissonDiskSampler::GetNextSample(Sample *sample) {
	if (!oneDSamples) {
		// Allocate space for pixel's poisson disk sample tables
		oneDSamples = new float *[sample->n1D.size()];
		n1D = sample->n1D.size();
		for (u_int i = 0; i < sample->n1D.size(); ++i)
			oneDSamples[i] = new float[sample->n1D[i] *
		                              2* samplesPerPixel];
		twoDSamples = new float *[sample->n2D.size()];
		n2D = sample->n2D.size();
		for (u_int i = 0; i < sample->n2D.size(); ++i)
			twoDSamples[i] = new float[2 * sample->n2D[i] *
		                              2* samplesPerPixel];
	}
	if (samplePos == pixelSamples) {
		// Advance to next pixel for poisson disk sampling
		if (++xPos == xPixelEnd) {
			xPos = xPixelStart;
			++yPos;
		}
		if (yPos == yPixelEnd)
			return false;
		samplePos = 0;
		
		// Generate poisson disk samples for pixel
		///////////
		// initialization
		int nSamples = 0;
		
		//copy the backup down samples over
		if(xPos == xPixelStart){
			memcpy(down, downBackup, 2*nX*(xPixelEnd-xPixelStart+1)*sizeof(int));
			memcpy(downSamples, downSamplesBackup, 4*nX*(xPixelEnd-xPixelStart+1)*sizeof(float));
			downIndex = 0;
		}
		
		memset(activeGrid, -1, nX*nY*sizeof(int));
		pixelSamples = genPattern2D(imageSamples, true);
/*
		for(int i=0; i<pixelSamples; i++)
			std::cout<<"[ "<<imageSamples[2*i]<<", "<<imageSamples[2*i+1]<<" ]"<<std::endl;
*/
		
//		std::cout<<"\nnumber of image samples generated: "<<pixelSamples<<std::endl;
	
		//store information to right and down
		//right
		int rightIndex=0;
		if(xPos != xPixelEnd ){
			for(int i=1; i<=nY; i++){
				right[2*(i-1)] = activeGrid[i*nX-2];
				if(right[2*(i-1)]>=0){
					rightSamples[rightIndex*2] = imageSamples[2*right[2*(i-1)]]-1.0;
					rightSamples[rightIndex*2+1] = imageSamples[2*right[2*(i-1)]+1];
					right[2*(i-1)] = rightIndex++;
				}
				right[2*(i-1)+1] = activeGrid[i*nX-1];
				if(right[2*(i-1)+1]>=0){
					rightSamples[rightIndex*2] = imageSamples[2*right[2*(i-1)+1]]-1.0;
					rightSamples[rightIndex*2+1] = imageSamples[2*right[2*(i-1)+1]+1];
					right[2*(i-1)+1] = rightIndex++;
				}
			}
		}else{
			memset(right, -1, 2*nY*sizeof(int));
		}
/*
		for(int i=0; i<pixelSamples; i++)
			std::cout<<"[ "<<imageSamples[2*i]<<", "<<imageSamples[2*i+1]<<" ]"<<std::endl;

*/	
/*			
		std::cout<<std::endl;	
		for(int i=0; i<nY; i++){
			for(int j=0; j<nX; j++)
				std::cout<<activeGrid[i*nX+j]<<' ';
			std::cout<<std::endl;
		}

		std::cout<<std::endl;
		for(int i=0; i<2*nY; i++)
			std::cout<<right[i]<<' ';
		std::cout<<std::endl;
*/
		//down
		int nCol = nX*(xPixelEnd-xPixelStart+1);
		for(int j=0; j<2; j++){
			for(int i=0; i<nX; i++){
				int index = j*nCol+(xPos-xPixelStart)*nX+i;
				downBackup[index] = activeGrid[(nY-2+j)*nX+i];
				if(downBackup[index]>=0){
					downSamplesBackup[downIndex*2] = imageSamples[2*downBackup[index]];
					downSamplesBackup[downIndex*2+1] = imageSamples[2*downBackup[index]+1]-1.0;
					downBackup[index] = downIndex++;
				}
			}
		}
/*	
		for(int i=0; i<pixelSamples; i++)
			std::cout<<"[ "<<imageSamples[2*i]<<", "<<imageSamples[2*i+1]<<" ]"<<std::endl;
*/
	
		
		Shuffle(imageSamples, pixelSamples, 2);
/*	
		memset(activeGrid, -1, nX*nY*sizeof(int));
		nSamples = genPattern2D(lensSamples);
		Shuffle(lensSamples, nSamples, 2);
//		std::cout<<"number of lensSamples generated: "<<nSamples<<std::endl;
		if(nSamples<pixelSamples) pixelSamples = nSamples;
		
		memset(activeGrid, -1, nX*nY*sizeof(int));
		nSamples = genPattern1D(timeSamples);
		Shuffle(timeSamples, nSamples, 1);
//		std::cout<<"number of timeSamples generated: "<<nSamples<<std::endl;
		/////////
*/		samplesPerPixel = pixelSamples;

		for (u_int i = 0; i < sample->n1D.size(); ++i)
			LDShuffleScrambled1D(sample->n1D[i], pixelSamples, oneDSamples[i]);
		for (u_int i = 0; i < sample->n2D.size(); ++i)
			LDShuffleScrambled2D(sample->n2D[i], pixelSamples, twoDSamples[i]);
	}
	// Copy poisson disk samples from tables
	sample->imageX = xPos + imageSamples[2*samplePos];
	sample->imageY = yPos + imageSamples[2*samplePos+1];
	sample->time = 0; // timeSamples[samplePos];
	sample->lensU = lensSamples[2*samplePos];
	sample->lensV = lensSamples[2*samplePos+1];
	for (u_int i = 0; i < sample->n1D.size(); ++i) {
		int startSamp = sample->n1D[i] * samplePos;
		for (u_int j = 0; j < sample->n1D[i]; ++j)
			sample->oneD[i][j] = oneDSamples[i][startSamp+j];
	}
	for (u_int i = 0; i < sample->n2D.size(); ++i) {
		int startSamp = 2 * sample->n2D[i] * samplePos;
		for (u_int j = 0; j < 2*sample->n2D[i]; ++j)
			sample->twoD[i][j] = twoDSamples[i][startSamp+j];
	}
/*
	std::ofstream xCoord("x.txt", std::ios_base::app);
	xCoord<<sample->imageX<<'\n';
	xCoord.close();

	std::ofstream yCoord("y.txt", std::ios_base::app);
	yCoord<<sample->imageY<<'\n';
	yCoord.close();
*/	
	++samplePos;
	return true;
}

int PoissonDiskSampler::genPattern1D(float* samples)
{
// WARNING: initialize activeGrid before entering this function
// TODO: check the use of ceil and floor

	int numSamples=0;

	// min sqrt distance between two points
	float sqr_radius = radius1D*radius1D;
	
	//declare active list
	int nGrids = Ceil2Int(1.0/gridSize1D);
//	int* activeList	= new int[nGrids];
	int activeSize = 0;
	
	//get first sample and put it into samples, activeList and activeGrid
	float f = RandomFloat(); //[0,1)
	activeList[activeSize++] = numSamples;
	
	int gridIndex = Floor2Int(f/gridSize1D);
	if (gridIndex>=nGrids) gridIndex = nGrids-1;
	if (gridIndex<0) gridIndex = 0;
	activeGrid[gridIndex] = numSamples;
	
	samples[numSamples++] = f;
//	std::cout<<"sample: "<<f<<std::endl;
	
	//while( activeList not empty )
	while( activeSize > 0 ){
	//  get a random sample p
		int t = Floor2Int(RandomFloat()*activeSize);
		float p = samples[activeList[t]];
		float newSample;
		bool foundSample = false;
	//  try maxAttemp to find samples around p
		for(int attemp=0; attemp<maxAttemp; attemp++){
			float l;
			//  new sample is p+l*radius; l within [-2,-1) or (1,2]
			do{
			//  get l within [-2,2]
				l = 4.0*(genrand_real1()-0.5); //TODO: rand need to include 1
			//  reject if l within [-1,1]
			}while(l>=-1 && l<=1);
			newSample = p + l*radius1D;
			//check against boundary
			if(newSample<0 || newSample>1)
				continue;
			//  find neighbors
			int minIndex = Floor2Int((newSample - radius1D)/gridSize1D);
			int maxIndex = Floor2Int((newSample + radius1D)/gridSize1D);
			if(minIndex<0) minIndex = 0;
			if(maxIndex >= nGrids) maxIndex = nGrids-1;
			//  check all the neighbors against this sample
			int i;
			for(i=minIndex; i<=maxIndex; i++){
				if(activeGrid[i]>=0){ // only check nonempty cells
					float neighbor = samples[activeGrid[i]];
					float sqr_dist = (neighbor-newSample)*(neighbor-newSample);
					if( sqr_dist < sqr_radius ) // reject this sample
						break;
				}
			}
			if(i>maxIndex){ // this sample is far from all the neighbors; accept it
				foundSample = true;
				break;
			}
		}
		if(foundSample){
	//  put away the new sample	
			int index = Floor2Int(newSample/gridSize1D);
			if(index<0) index = 0;
			if(index>=nGrids) index = nGrids-1;
			activeGrid[index] = numSamples;
			activeList[activeSize++] = numSamples;
			samples[numSamples++] = newSample;
//			std::cout<<"sample: "<<newSample<<std::endl;
		}else{
	//  remove the current sample if no new valid sample found around it
			activeList[t] = activeList[--activeSize];
		}
	}//end while
//	delete []activeList;
	return numSamples;
}

int PoissonDiskSampler::genPattern2D(float* samples, bool isImage)
{
// WARNING: initialize activeGrid before entering this function
// TODO: check the use of ceil and floor

	int numSamples=0;
	
	// min sqrt distance between two points
	float sqr_radius = radius2D*radius2D;
	
	//declare active list
//	int* activeList	= new int[nX*nY];
	int activeSize = 0;
	
	//get first sample and put it into samples, activeList and activeGrid
	float fx = genrand_real1()*0.5+0.25; //[0.25, 0.75]
	float fy = genrand_real1()*0.5+0.25;
	activeList[activeSize++] = numSamples;
	
	int gridIndexX = Floor2Int(fx/gridSize2D);
	if (gridIndexX>=nX) gridIndexX = nX-1;
	if (gridIndexX<0) gridIndexX = 0;
	int gridIndexY = Floor2Int(fy/gridSize2D);
	if (gridIndexY>=nY) gridIndexY = nY-1;
	if (gridIndexY<0) gridIndexY=0;
	activeGrid[nX*gridIndexY+gridIndexX] = numSamples;
	
	samples[numSamples*2] = fx;
	samples[numSamples*2+1] = fy;
	numSamples++;

//	std::cout<<"\nsample: ["<<fx<<", "<<fy<<"] "<<std::endl;
//	std::cout<<"index: ("<<gridIndexX<<", "<<gridIndexY<<") "<<std::endl;
	//while( activeList not empty )
	while( activeSize > 0 ){
	//  get a random sample p
		int t = Floor2Int(RandomFloat()*activeSize);
		float px = samples[2*activeList[t]];
		float py = samples[2*activeList[t]+1];
		float newSampleX;
		float newSampleY;
		bool foundSample = false;
	//  try maxAttemp to find samples around p
		for(int attemp=0; attemp<maxAttemp; attemp++){
			float lx, ly, sqr_len;
			//  new sample is p+l*radius; len(l) within [-2,-1) or (1,2]
			do{
			//  get lx, ly within [-2,2]
				lx = 4.0*(genrand_real1()-0.5); //TODO: rand need to include 1
				ly = 4.0*(genrand_real1()-0.5);
				sqr_len = lx*lx +ly*ly;
			//  reject if len within [-1,1]
			}while(sqr_len<=1 || sqr_len>4);
			newSampleX = px + lx*radius2D;
			newSampleY = py + ly*radius2D;
		
			//check against boundary
			if(newSampleX<0 || newSampleX>1 || newSampleY<0 || newSampleY>1)
				continue;
		
			//  find neighbors
			int minIndexX= Floor2Int((newSampleX - radius2D)/gridSize2D);
			int maxIndexX = Floor2Int((newSampleX + radius2D)/gridSize2D);
			
			int minIndexY = Floor2Int((newSampleY - radius2D)/gridSize2D);
			int maxIndexY = Floor2Int((newSampleY + radius2D)/gridSize2D);
			
			if(maxIndexY >= nY) maxIndexY = nY-1; 
			if(minIndexY<=0){ // need to check against down for iamge samples
				minIndexY = 0;
				if(isImage){
					int i,j;
					if(xPos == xPixelStart && minIndexX<0)
						minIndexX = 0;
					if(xPos == xPixelEnd && maxIndexX>=nX)
						maxIndexX = nX-1;
					//out of range min and max IndexY help to check diagonal grids
					for(j=minIndexY; j<2; j++){
						for(i=minIndexX; i<=maxIndexX; i++){
							int index = j*nX*(xPixelEnd-xPixelStart+1)+(xPos-xPixelStart)*nX + i;
							if(down[index]>=0){
								float neighborX = downSamples[down[index]*2];
								float neighborY = downSamples[down[index]*2+1];
								float sqr_dist = (neighborX-newSampleX)*(neighborX-newSampleX)+(neighborY-newSampleY)*(neighborY-newSampleY);
								if( sqr_dist < sqr_radius ) // reject this sample
									break;
							}
						}
						if(i<=maxIndexX) break;
					}
					if(j<2) continue;
				}
			}
			if(minIndexX<=0){
				minIndexX = 0;
				if( isImage && xPos != xPixelStart ){ // need to check against right
					int i,j;
					for(j=minIndexX; j<2; j++){
						for(i=minIndexY; i<=maxIndexY; i++){
							if(right[i*2+j]>=0){
								float neighborX = rightSamples[right[i*2+j]*2];
								float neighborY = rightSamples[right[i*2+j]*2+1];
								float sqr_dist = (neighborX-newSampleX)*(neighborX-newSampleX)+(neighborY-newSampleY)*(neighborY-newSampleY);
								if( sqr_dist < sqr_radius ) // reject this sample
									break;
							}
						}
						if(i<=maxIndexY) break;
					}
					if(j<2) continue;
				}
				minIndexX=0;
			}
			if(maxIndexX>=nX) maxIndexX = nX-1;
			
			//  check all the neighbors against this sample
			int i, j;
			for(i=minIndexY; i<=maxIndexY; i++){
				for(j=minIndexX; j<=maxIndexX; j++){
					if(activeGrid[nX*i+j]>=0){ // only check nonempty cells
						int index = activeGrid[nX*i+j];
						float neighborX = samples[2*index];
						float neighborY = samples[2*index+1];
						float sqr_dist = (neighborX-newSampleX)*(neighborX-newSampleX)+(neighborY-newSampleY)*(neighborY-newSampleY);
						if( sqr_dist < sqr_radius ) // reject this sample
							break;
					}
				}
				if(j<=maxIndexX) break;
			}
			if(i>maxIndexY){ // this sample is far from all the neighbors; accept it
				foundSample = true;
				break;
			}
		} // end attemp
		if(foundSample){
	//  put away the new sample	
			int indexX = Floor2Int(newSampleX/gridSize2D);
			int indexY = Floor2Int(newSampleY/gridSize2D);
			if(indexX<0) indexX = 0;
			if(indexX>=nX) indexX = nX-1;
			if(indexY<0) indexY = 0;
			if(indexY>=nY) indexY = nY-1;
			activeGrid[nX*indexY+indexX] = numSamples;
			activeList[activeSize++] = numSamples;
			samples[2*numSamples] = newSampleX;
			samples[2*numSamples+1] = newSampleY;
			numSamples++;
//			std::cout<<"sample: ["<<newSampleX<<", "<<newSampleY<<"] "<<std::endl;
		}else{
	//  remove the current sample if no new valid sample found around it
			activeList[t] = activeList[--activeSize];
		}
	}//end while
//	delete []activeList;
	return numSamples;
}

extern "C" DLLEXPORT Sampler *CreateSampler(const ParamSet &params, const Film *film) {
	// Initialize common sampler parameters
	int xstart, xend, ystart, yend;
	film->GetSampleExtent(&xstart, &xend, &ystart, &yend);
	int nsamp = params.FindOneInt("pixelsamples", 4);
	return new PoissonDiskSampler(xstart, xend, ystart, yend, nsamp);
}
