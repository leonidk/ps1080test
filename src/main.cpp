#include "ONICamera.h"
#include "linalg.h"
#include "detect.h"
#include "ICP.h"
#include "Drawing.h"

#define SUBSAMPLE_FACTOR (4)
#define DIST_THRESH (40)

int main(int argc, char *args[]) {

	int posX = 100, posY = 100, width = 320, height = 240;
	uint16_t *hDepth = new uint16_t[width * height];

	width /= SUBSAMPLE_FACTOR;
	height /= SUBSAMPLE_FACTOR;

	SDLWrapper depthWin("Depth", posX, posY, width, height, 8);
	SDLWrapper planesWin("PlanesColor", posX, posY, width, height, 8);

	SDLWrapper normalWin("Normals", posX + width*8, posY, width, height, 8);
	SDLWrapper visWin("Planes", posX + width * 2*8, posY, 31, 31, 15);

	std::vector<float> points(3 * width * height, 0);
	std::vector<float> normals(3 * width * height, 0);

	std::vector<float> pointsPrev(3 * width * height, 0);
	std::vector<float> normalsPrev(3 * width * height, 0);

	std::vector<planeCandidate> prevPlanes;
	std::vector<uint8_t> planeDebugImage_Color(width*height * 3,0);

	ONICamera cam;
	auto camRunning = cam.Start();
	auto fx = cam.getFx() / SUBSAMPLE_FACTOR;
	auto fy = cam.getFy() / SUBSAMPLE_FACTOR;
	auto px = cam.getPx() / SUBSAMPLE_FACTOR;
	auto py = cam.getPy() / SUBSAMPLE_FACTOR;
	planeDetector<15> pd;
	iteratedICP icp(width, height);
	while (camRunning) {
		//get camera data
		camRunning = cam.syncNext();
		memset(normals.data(), 0, sizeof(float) * width * height * 3);
		memset(points.data(), 0, sizeof(float) * width * height * 3);

		//subsample camera data & generate data
		generateHalfImage<uint16_t, SUBSAMPLE_FACTOR>(cam.getDepth(), hDepth, width, height);
		generatePoints(hDepth, width, height, fx, fy, px,py, points.data());
		generateNormals<1>(hDepth, width, height, fx, fy, px, py, normals.data());

		auto candidates =  pd.detectPlanes(512, 2.5, width, height, points.data(), normals.data());
		memset(planeDebugImage_Color.data(), 0, width*height * 3);
		for (size_t i = 0; i < candidates.size(); i++) {
			const auto color = visWin.getNewColor(i);
			const auto plane = candidates[i];
			const auto error = std::min<float>(150.0f, plane.stddev * 3);
			for (int j = 0; j < width*height; j++) {
				const auto distFromPlane = plane.n[0] * points[3 * j] + plane.n[1] * points[3 * j + 1] + plane.n[2] * points[3 * j + 2] + plane.d;
				if (fabs(distFromPlane) < error) {
					planeDebugImage_Color[3 * j] = color[0];
					planeDebugImage_Color[3 * j+1] = color[1];
					planeDebugImage_Color[3 * j+2] = color[2];	  
				}							  
			}
		}
		//visualize output
		depthWin.setGrayScale<4>(hDepth);
		normalWin.setNormals(normals.data());
		visWin.setGrayScale<1>(pd.getDebugImg());
		planesWin.setRGB(planeDebugImage_Color.data());
		//figure out which plane came from where
		//auto corrPairs = generateCorrespondencesSVD(candidates, prevPlanes);
		//auto corrPairs2 = generateCorrespondencesMP(candidates, prevPlanes);
		//auto trans = icp.runICPIter(10, fx, fy, 0.5, 15, 100, points.data(), normals.data(), pointsPrev.data(), normalsPrev.data());

		//save frame state
		//prevPlanes = candidates;
		//normalsPrev = normals;
		//pointsPrev = points;

		SDL_Event e;
		if (SDL_PollEvent(&e)) {
			if (e.type == SDL_QUIT) {
				break;
			}
		}
		depthWin.Draw();
		normalWin.Draw();
		visWin.Draw();
		planesWin.Draw();
	}

	return 0;
}