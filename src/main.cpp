#include <SDL.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include <chrono>
#include "ONICamera.h"
#include "linalg.h"
#include "detect.h"

#include <iostream>

int main(int argc, char *args[]){
	
	SDL_Window *win = NULL;
	SDL_Renderer *renderer = NULL;

	SDL_Window *win2 = NULL;
	SDL_Renderer *renderer2 = NULL;

	SDL_Window *win3 = NULL;
	SDL_Renderer *renderer3 = NULL;

	SDL_Texture *bitmapTex = NULL;
	SDL_Surface *bitmapSurface = NULL;
	int posX = 100, posY = 100, width = 320, height = 240;
	DepthPixel* hDepth = new DepthPixel[width*height];

	width /= 4;
	height /= 4;

	win = SDL_CreateWindow("Hello World", posX, posY, 640, 480, 0);
	win2 = SDL_CreateWindow("Hello World2", posX + width, posY, 640, 480, 0);
	win3 = SDL_CreateWindow("Hello World3", posX + width * 2, posY, 500, 500, 0);

	renderer = SDL_CreateRenderer(win, -1, SDL_RENDERER_ACCELERATED);
	renderer2 = SDL_CreateRenderer(win2, -1, SDL_RENDERER_ACCELERATED);
	renderer3 = SDL_CreateRenderer(win3, -1, SDL_RENDERER_ACCELERATED);

	SDL_Surface *depthFrameSur = SDL_CreateRGBSurface(0, width, height, 32, 0, 0, 0, 0);
	SDL_Surface *normsFrameSur = SDL_CreateRGBSurface(0, width, height, 32, 0, 0, 0, 0);
	SDL_Surface *visBinsSur = SDL_CreateRGBSurface(0, RAD_FULL, RAD_FULL, 32, 0, 0, 0, 0);


	std::vector<float> points(3 * width*height, 0);
	std::vector<float> normals(3 * width*height, 0);


	std::vector<float> pointsPrev(3 * width*height, 0);
	std::vector<float> normalsPrev(3 * width*height, 0);

	
	std::vector<planeCandidate> prevPlanes;
	while (true)
	{
		memset(normals.data(), 0, sizeof(float)*width*height * 3);
		memset(points.data(), 0, sizeof(float)*width*height * 3);



		//LARGE_INTEGER StartingTime, EndingTime, MiddleTime,ElapsedMicroseconds;
		//LARGE_INTEGER Frequency;
		//QueryPerformanceFrequency(&Frequency);

		//QueryPerformanceCounter(&StartingTime);
		generateHalfImage<DepthPixel, 4>(pDepth, hDepth, width, height);
		generatePoints(hDepth, width, height, fx/4, fy/4, points.data());
		//QueryPerformanceCounter(&MiddleTime);

		static uint16_t ret[RAD_FULL * RAD_FULL];
		auto candidates = generateNormals<1>(points.data(), width, height, normals.data(),ret);
		//generateNormals<2>(pDepth, width, height, fx, fy, normals.data());
		//QueryPerformanceCounter(&EndingTime);

		drawGrayScale(depthFrameSur, hDepth);
		//float rot[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
		//float trans[3] = { 0, 0, 150 };
		//for (auto &&d : pointsPrev) {
		//	d = 0;
		//}
		//for (auto &&d : normalsPrev) {
		//	d = 0;
		//}
		//transformImage(width, height, fx/4, fy/4, points.data(), normals.data(), pointsPrev.data(), normalsPrev.data(), rot, trans);
		//drawNormals(depthFrameSur, normalsPrev.data());

		drawGrayScale(visBinsSur, ret);
		drawNormals(normsFrameSur, normals.data());

		auto corrPairs = generateCorrespondencesSVD(candidates, prevPlanes);
		auto corrPairs2 = generateCorrespondencesMP(candidates, prevPlanes);

		computeLinearApproxICP(width, height, 0.5, 15, points.data(), normals.data(), pointsPrev.data(), normalsPrev.data());

		std::cout << candidates.size()<< std::endl;
		//std::cout << std::endl;

		prevPlanes = candidates;
		normalsPrev = normals;
		pointsPrev = points;
		//ElapsedMicroseconds.QuadPart = MiddleTime.QuadPart - StartingTime.QuadPart;
		//ElapsedMicroseconds.QuadPart *= 1000000;
		//ElapsedMicroseconds.QuadPart /= Frequency.QuadPart;

		//printf("%lf ", static_cast<double>(ElapsedMicroseconds.QuadPart));

		//ElapsedMicroseconds.QuadPart = EndingTime.QuadPart - MiddleTime.QuadPart;
		//ElapsedMicroseconds.QuadPart *= 1000000;
		//ElapsedMicroseconds.QuadPart /= Frequency.QuadPart;
		//printf("%lf\n", static_cast<double>(ElapsedMicroseconds.QuadPart));

		SDL_Event e;
		if (SDL_PollEvent(&e)) {
			if (e.type == SDL_QUIT) {
				break;
			}
		}
		SDL_Texture *depthFrameTex = SDL_CreateTextureFromSurface(renderer, depthFrameSur);

		SDL_RenderClear(renderer);
		SDL_RenderCopy(renderer, depthFrameTex, NULL, NULL);
		SDL_RenderPresent(renderer);

		SDL_Texture *normFrameTex = SDL_CreateTextureFromSurface(renderer2, normsFrameSur);

		SDL_RenderClear(renderer2);
		SDL_RenderCopy(renderer2, normFrameTex, NULL, NULL);
		SDL_RenderPresent(renderer2);

		SDL_Texture *visBinsTex = SDL_CreateTextureFromSurface(renderer3, visBinsSur);

		SDL_RenderClear(renderer3);
		SDL_RenderCopy(renderer3, visBinsTex, NULL, NULL);
		SDL_RenderPresent(renderer3);

		SDL_DestroyTexture(visBinsTex);
		SDL_DestroyTexture(depthFrameTex);
		SDL_DestroyTexture(normFrameTex);
	}


	SDL_DestroyTexture(bitmapTex);
	SDL_DestroyRenderer(renderer);
	SDL_DestroyWindow(win);

	return 0;
}
