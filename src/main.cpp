#include "ONICamera.h"
#include "ImageFunc.h"
#include "ICP.h"
#include "Drawing.h"
#include "detect.h"
#include "detectMS.h"

#define GLEW_STATIC
#include <GL/glew.h>
//#include "argh.h"

#define SUBSAMPLE_FACTOR (1)
#define DIST_THRESH (40)
using namespace linalg::aliases;
int3 hashCoeff;
unsigned int collisions = 0;
namespace detail {
	struct VoxelCandidate {
		int3 pos = { 0, 0, 0 };
		float3 sum = { 0, 0, 0 };
		int cnt = 0;
	};
}
//must be a power of two
template <int numVoxels>
void voxelSubsample(const std::vector<float3> & input, std::vector<float3> &output, float voxelSize, int minVoxelNum)
{
	using detail::VoxelCandidate;
	VoxelCandidate cand[numVoxels];
	auto voxelMask = numVoxels - 1;
	auto iVS = 1.0f / voxelSize; //inverse voxel size
	//auto hashCoeff = int3(3863,2029,49139); //good for mm
	//auto hashCoeff = int3(54851,11909,24781); //good for m
	for (const auto & pt : input) {
		//not using floor since floor would map (-1,1) to 0. 
		int3 intPos(static_cast<int>(pt.x*iVS), static_cast<int>(pt.y*iVS), static_cast<int>(pt.z*iVS));
		unsigned int hash = dot(hashCoeff, intPos);
		hash&=voxelMask;
		auto & bucket = cand[hash];
		if (bucket.cnt && bucket.pos != intPos) { //flush on collision
			bucket.cnt = 0;
			bucket.sum = { 0, 0, 0 };
			collisions++;
		}
		if (!bucket.cnt) {
			bucket.pos = intPos;
		}
		bucket.sum += pt;
		bucket.cnt++;
	}
	for (const auto &pt : cand)
	{
		if (pt.cnt >= minVoxelNum)
			output.emplace_back(pt.sum / static_cast<float>(pt.cnt));
	}
}
//must be a power of two
template <int numVoxels>
std::vector<float3> voxelSubsample(const std::vector<float3> & input, float voxelSize, int minVoxelNum)
{
	std::vector<float3> output;
	voxelSubsample<numVoxels>(input, output, voxelSize, minVoxelNum);
	return output;
}
#include <random>
int main(int argc, char *args[]) {
	//std::random_device rd;
	//std::default_random_engine rng(rd());
	//std::uniform_int_distribution<> rp(0, 6492);
	//std::uniform_int_distribution<> depthStart(500, 1500);
	//std::uniform_int_distribution<> depthRang(0, 170);
	//std::vector<int3> oldHashes = { { 7171, 3079, 4231 }, { 12799, 13241, 29947 }, { 7127, 19891, 11159 }, { 44281, 4517, 30851 }, 
	//{ 12373, 42293, 8999 }, { 56599, 11399, 33359 }, { 7723, 13241, 53717 }, { 3863, 49139, 2029 }, { 3863, 2029, 49139 }, { 54851, 11909, 24781 } };
	//int3 bstCoeff;
	//int bstColl = INT_MAX;
	//int hashCnt = 0;
	//while (true) {
	//	if (hashCnt < oldHashes.size())
	//		hashCoeff = oldHashes[hashCnt++];
	//	else
	//		hashCoeff = int3(primes[rp(rng)], primes[rp(rng)], primes[rp(rng)]);
	//	collisions = 0;
	//	for (int iter = 0; iter < 25; iter++){
	//		std::vector<float3> depth, subsample;
	//		auto start = depthStart(rng);
	//		for (int y = 0; y < 480; y++){
	//			for (int x = 0; x < 640; x++) {
	//				auto z = (start + depthRang(rng));
	//				depth.emplace_back((x - 320.0f)*z / 480.0f, (x - 240)*z / 480.0f, z);
	//			}
	//		}
	//		voxelSubsample<1024>(depth, subsample, 10.0f, 1);
	//	}
	//	if (collisions < bstColl) {
	//		bstColl = collisions;
	//		bstCoeff = hashCoeff;
	//		std::cout << "\n" << hashCoeff[0] << ',' << hashCoeff[1] << ',' << hashCoeff[2] << std::endl;
	//		std::cout << collisions << '\n' << std::endl;
	//	}
	//	else {
	//		if (hashCnt < oldHashes.size()) {
	//			std::cout << "\n" << hashCoeff[0] << ',' << hashCoeff[1] << ',' << hashCoeff[2] << std::endl;
	//			std::cout << collisions << '\n' << std::endl;
	//		}
	//		else {
	//			std::cout << ".";
	//		}
	//	}
	//}

	SDL_Init(SDL_INIT_EVERYTHING);

	int posX = 100, posY = 100, width = 320, height = 240;
	uint16_t *hDepth = new uint16_t[width * height];
	width /= SUBSAMPLE_FACTOR;
	height /= SUBSAMPLE_FACTOR;

	Draw2DImage depthWin("Depth", posX, posY, width, height, 2 * SUBSAMPLE_FACTOR);
	Draw2DImage planesWin("PlanesColor", posX, posY, width, height, 2*SUBSAMPLE_FACTOR);
	Draw2DImage normalWin("Normals", posX, posY, width, height, 4 * SUBSAMPLE_FACTOR);
	Draw2DImage visWin("Planes", posX + width * 2 * SUBSAMPLE_FACTOR, posY, 31, 31, 15);
	Draw3DImage glTest("OpenGLWindow", posX, posY + height * SUBSAMPLE_FACTOR, width, height, SUBSAMPLE_FACTOR, Draw3DImage::DRAW_MODE_POINTS);


	std::vector<float> points(3 * width * height, 0);
	std::vector<float> normals(3 * width * height, 0);

	std::vector<float> pointsPrev(3 * width * height, 0);
	std::vector<float> normalsPrev(3 * width * height, 0);

	std::vector<planeCandidate> prevPlanes;
	std::vector<uint8_t> planeDebugImage_Color(width*height * 4,0);

	ONICamera cam;
	auto quit = !cam.Start();
	auto fx = cam.getFx() / SUBSAMPLE_FACTOR;
	auto fy = cam.getFy() / SUBSAMPLE_FACTOR;
	auto px = cam.getPx() / SUBSAMPLE_FACTOR;
	auto py = cam.getPy() / SUBSAMPLE_FACTOR;
	//planeDetectorFast<15> pdC;
	planeDetectorDisjointModes<15> pd;
	//planeDetectorDisjointGMM<15> pd;
	planeDetectorMS pdMS;
	iteratedICP icp(width, height);
	while (!quit) {
		//get camera data
		quit = !cam.syncNext();
		memset(normals.data(), 0, sizeof(float) * width * height * 3);
		memset(points.data(), 0, sizeof(float) * width * height * 3);

		//subsample camera data & generate data
		generateHalfImageDepth<uint16_t, SUBSAMPLE_FACTOR>(cam.getDepth(), hDepth, width, height);
		generatePoints(hDepth, width, height, fx, fy, px,py, points.data());
		generateNormals_FromDepth<1>(hDepth, width, height, fx, fy, px, py, normals.data());
		//generateNormals_fromPoints<1>(points.data(), width, height, normals.data());

		//auto candidatesC = pdC.detectPlanes(512, 3.0f, width, height, points.data(), normals.data());
		auto candidates = pd.detectPlanes(512, 100.0f, width, height, points.data(), normals.data());
		//auto candidates = pdMS.detectPlanes(512, 0.15f, (float)(0.125 / 64.0), width, height, points.data(), normals.data());

		memset(planeDebugImage_Color.data(), 0, width*height * 4);
		for (int i = 0; i < candidates.size(); i++) {
			const auto color = getNewColor(i);
			const auto plane = candidates[i];
			const auto error = std::min<float>(25.0f, plane.stddev * 3);
			for (int j = 0; j < width*height; j++) {
				const auto distFromPlane = plane.n[0] * points[3 * j] + plane.n[1] * points[3 * j + 1] + plane.n[2] * points[3 * j + 2] + plane.d;
				if (fabs(distFromPlane) < error) {
					planeDebugImage_Color[4 * j] = color[0];
					planeDebugImage_Color[4 * j+1] = color[1];
					planeDebugImage_Color[4 * j+2] = color[2];	  
				}							  
			}
		}
		////visualize output
		depthWin.setGrayScale<4>(hDepth);
		normalWin.setNormals(normals.data());
		visWin.setGrayScale<1>(pd.getDebugImg());
		planesWin.setRGBA (planeDebugImage_Color.data());
		glTest.setPoints(points.data(), normals.data(), nullptr,false);
		//figure out which plane came from where
		//auto corrPairs = generateCorrespondencesSVD(candidates, prevPlanes);
		//auto corrPairs2 = generateCorrespondencesMP(candidates, prevPlanes);
		auto transform = icp.runICPIter(10, fx, fy, 0.5, 15, 100, points.data(), normals.data(), pointsPrev.data(), normalsPrev.data());
		//std::cout << transform << std::endl;
		//save frame state
		prevPlanes = candidates;
		normalsPrev = normals;
		pointsPrev = points;

		SDL_Event e;
		while (SDL_PollEvent(&e)){
			if (e.type == SDL_QUIT)
				quit = true;
			if (e.type == SDL_KEYUP)
				if (e.key.keysym.sym == SDLK_ESCAPE)
					quit = true;
			glTest.handleEvent(e);
		}
		depthWin.Draw();
		normalWin.Draw();
		visWin.Draw();
		planesWin.Draw();
		glTest.Draw();
	}
	SDL_Quit();
	return 0;
}