#include "ONICamera.h"
#include "linalg.h"
#include "ICP.h"
#include "Drawing.h"
#include "detect.h"
#include "detectMS.h"
#define GLEW_STATIC
#include <GL/glew.h>
//#include "argh.h"

#define SUBSAMPLE_FACTOR (4)
#define DIST_THRESH (40)

int main(int argc, char *args[]) {

	SDL_Init(SDL_INIT_EVERYTHING);

	int posX = 100, posY = 100, width = 320, height = 240;
	uint16_t *hDepth = new uint16_t[width * height];
	auto a = std::vector<int>({ 1, 2, 3, 4 });
	auto b = imap2(a, [](int a) { return (float)a; });
	width /= SUBSAMPLE_FACTOR;
	height /= SUBSAMPLE_FACTOR;

	SDLWrapperGL depthWin("Depth", posX, posY, width, height, 8);
	SDLWrapperGL planesWin("PlanesColor", posX, posY, width, height, 8);
	SDLWrapperGL normalWin("Normals", posX + width * 8, posY, width, height, 8);
	SDLWrapperGL visWin("Planes", posX + width * 2 * 8, posY, 31, 31, 15); 
	SDLWrapperGL glTest("OpenGLWindow", posX, posY + height * 8, width, height, 8, SDLWrapperGL::DRAW_MODE_POINTS);


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
	planeDetectorFast<15> pd;
	//planeDetectorDisjoint<15> pd;
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


		auto candidates =  pd.detectPlanes(512, 2.5, width, height, points.data(), normals.data());
		//auto candidates = pdMS.detectPlanes(512, 0.15f, (float)(0.125 / 64.0), width, height, points.data(), normals.data());

		memset(planeDebugImage_Color.data(), 0, width*height * 3);
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
		glTest.setPoints(points.data(), normals.data(), nullptr);
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