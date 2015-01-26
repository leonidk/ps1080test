

#include "ONICamera.h"
#include "linalg.h"
#include "detect.h"
#include "ICP.h"
#include "Drawing.h"

#define SUBSAMPLE_FACTOR (4)

int main(int argc, char *args[]) {

	int posX = 100, posY = 100, width = 320, height = 240;
	uint16_t *hDepth = new uint16_t[width * height];

	width /= SUBSAMPLE_FACTOR;
	height /= SUBSAMPLE_FACTOR;

	SDLWrapper depthWin("Depth", posX, posY, width, height, 8);
	SDLWrapper normalWin("Normals", posX + width, posY, width, height, 8);
	SDLWrapper visWin("Planes", posX + width * 2, posY, RAD_FULL, RAD_FULL, 15);

	std::vector<float> points(3 * width * height, 0);
	std::vector<float> normals(3 * width * height, 0);

	std::vector<float> pointsPrev(3 * width * height, 0);
	std::vector<float> normalsPrev(3 * width * height, 0);

	std::vector<planeCandidate> prevPlanes;
	uint16_t planeDebugImage[RAD_FULL * RAD_FULL];

	ONICamera cam;
	auto camRunning = cam.Start();
	auto fx = cam.getFx() / SUBSAMPLE_FACTOR;
	auto fy = cam.getFy() / SUBSAMPLE_FACTOR;
	planeDetector<15> pd;
	while (camRunning) {
		//get camera data
		camRunning = cam.syncNext();
		memset(normals.data(), 0, sizeof(float) * width * height * 3);
		memset(points.data(), 0, sizeof(float) * width * height * 3);

		//subsample camera data & generate data
		generateHalfImage<uint16_t, SUBSAMPLE_FACTOR>(cam.getDepth(), hDepth, width, height);
		generatePoints(hDepth, width, height, fx, fy, points.data());
		generateNormals<1>(hDepth, width, height, fx, fy, normals.data());

		auto candidates =  pd.detectPlanes(15, 75, width, height, points.data(), normals.data());
		//auto candidates = generateNormalsAndPlanes<1>(points.data(), width, height, normals.data(), planeDebugImage);

		//visualize output
		depthWin.setGrayScale<4>(hDepth);
		normalWin.setNormals(normals.data());
		visWin.setGrayScale<1>(pd.getDebugImg());

		//figure out which plane came from where
		//auto corrPairs = generateCorrespondencesSVD(candidates, prevPlanes);
		//auto corrPairs2 = generateCorrespondencesMP(candidates, prevPlanes);

		for (int i = 0; i < 4; i++) {
			auto transf = computeLinearApproxICP(width, height, 0.5, 15, 100, points.data(), normals.data(), pointsPrev.data(), normalsPrev.data());
			std::cout << transf << std::endl;
			//transformImage(width, height, fx, fy, points.data(), normals.data(), pointsPrev.data(), normalsPrev.data(), rotation, translation);
		}

		//save frame state
		prevPlanes = candidates;
		normalsPrev = normals;
		pointsPrev = points;

		SDL_Event e;
		if (SDL_PollEvent(&e)) {
			if (e.type == SDL_QUIT) {
				break;
			}
		}
		depthWin.Draw();
		normalWin.Draw();
		visWin.Draw();
	}

	return 0;
}
