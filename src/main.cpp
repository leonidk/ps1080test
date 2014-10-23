#include <stdio.h>
#include <OpenNI.h>
#include <SDL.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include <chrono>
using namespace openni;

void drawGrayScale(SDL_Surface* sur, DepthPixel * data) {
	uint8_t* dest = (uint8_t*)sur->pixels;
	for (int i = 0; i < sur->h*sur->w; i++)
		dest[4 * i + 2] = dest[4 * i + 1] = dest[4 * i] = (data[i] >> 1) > 255 ? 255 : (data[i] >> 1);
}

void drawGrayScale(SDL_Surface* sur, uint8_t * data) {
	uint8_t* dest = (uint8_t*)sur->pixels;
	for (int i = 0; i < sur->h*sur->w; i++)
		dest[4 * i + 2] = dest[4 * i + 1] = dest[4 * i] = data[i];
}

void drawNormals(SDL_Surface* sur, float* normals) {
	uint8_t* dest = (uint8_t*)sur->pixels;
	for (int i = 0; i < sur->h*sur->w; i++)
	{
		dest[4 * i] = 128 + normals[3 * i] * 127;
		dest[4 * i + 1] = 128 + normals[3 * i + 1] * 127;
		dest[4 * i + 2] = 128 + normals[3 * i + 2] * 127;

	}
}

template <typename T, int size>
void generateHalfImage(const T * in, T * out, const int outW, const int outH) {
	const auto inH = size*outH;
	const auto inW = size*outW;
	for (int i = 0; i < inH; i += size) {
		for (int j = 0; j < inW; j += size) {
			auto tout = 0;
			for (int m = 0; m < size; m++) {
				for (int n = 0; n < size; n++) {
					tout += in[(i + m)*inW + (j + n)];
				}
			}
			out[(i / size)*outW + (j / size)] = tout / (size*size);
		}
	}
}

void generatePoints(const DepthPixel* depth, const int width, const int height, const float fx, const float fy, float* points)
{
	auto halfX = width / 2;
	auto halfY = height / 2;
	auto cX = 1.0f / fx;
	auto cY = 1.0f / fy;
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			const auto z = depth[(i*width + j)];
			points[3 * (i*width + j)] = (j - halfX)*z*cX;
			points[3 * (i*width + j) + 1] = (i - halfY)*z*cY;
			points[3 * (i*width + j) + 2] = z;
		}
	}
}
// cX*((5*z +j*(z-z2));
inline void cross(const float* a, const float *b, float *c)
{
	c[0] = (a[1] * b[2]) - (a[2] * b[1]);
	c[1] = (a[2] * b[0]) - (a[0] * b[2]);
	c[2] = (a[0] * b[1]) - (a[1] * b[0]);
}
inline float dot(const float* a, const float* b)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
inline void normalize(float* a)
{
	const auto norm = sqrt(dot(a, a));
	a[0] /= norm;
	a[1] /= norm;
	a[2] /= norm;
}
#define RAD_SIZE (15)
#define RAD_FULL (2*RAD_SIZE+1)

#define MAX_CAND (5)
#define MIN_PTS  (50)

struct planeCandidate {
	double d;
	float n[3];
};

template <int size>
std::vector<planeCandidate> generateNormals(const float* points, const int width, const int height, float* normals, uint16_t* image)
{
	std::vector<planeCandidate> returnData;
	static std::vector<planeCandidate> planes[RAD_FULL*RAD_FULL];
	for (auto &p : planes) {
		p.clear();
	}

	memset(image, 0, 2 * RAD_FULL * RAD_FULL);
	for (int i = size; i < height - size; i++) {
		for (int j = size; j < width - size; j++) {
			if (!points[3 * (i*width + j) + 2])		continue;
			const float* pc = &points[3 * (i*width + j)];

			float outNorm[3] = { 0, 0, 0 };
			int count = 0;
			if (points[3 * (i*width + j + size) + 2] && points[3 * ((i + size)*width + j) + 2]){
				const float* px = &points[3 * (i*width + j + size)];
				const float* py = &points[3 * ((i + size)*width + j)];

				float v1[] = { px[0] - pc[0], px[1] - pc[1], px[2] - pc[2] };
				float v2[] = { py[0] - pc[0], py[1] - pc[1], py[2] - pc[2] };

				float v3[3];
				cross(v1, v2, v3);
				normalize(v3);
				outNorm[0] += v3[0];
				outNorm[1] += v3[1];
				outNorm[2] += v3[2];
				count++;
			}
			if (points[3 * (i*width + j - size) + 2] && points[3 * ((i + size)*width + j) + 2]){
				const float* px = &points[3 * (i*width + j - size)];
				const float* py = &points[3 * ((i + size)*width + j)];

				float v1[] = { pc[0] - px[0], pc[1] - px[1], pc[2] - px[2] };
				float v2[] = { py[0] - pc[0], py[1] - pc[1], py[2] - pc[2] };

				float v3[3];
				cross(v1, v2, v3);
				normalize(v3);
				outNorm[0] += v3[0];
				outNorm[1] += v3[1];
				outNorm[2] += v3[2];
				count++;
			}
			if (points[3 * (i*width + j + size) + 2] && points[3 * ((i - size)*width + j) + 2]){
				const float* px = &points[3 * (i*width + j + size)];
				const float* py = &points[3 * ((i - size)*width + j)];

				float v1[] = { px[0] - pc[0], px[1] - pc[1], px[2] - pc[2] };
				float v2[] = { pc[0] - py[0], pc[1] - py[1], pc[2] - py[2] };

				float v3[3];
				cross(v1, v2, v3);
				normalize(v3);
				outNorm[0] += v3[0];
				outNorm[1] += v3[1];
				outNorm[2] += v3[2];
				count++;
			}
			if (points[3 * (i*width + j - size) + 2] && points[3 * ((i - size)*width + j) + 2]){
				const float* px = &points[3 * (i*width + j - size)];
				const float* py = &points[3 * ((i - size)*width + j)];

				float v1[] = { pc[0] - px[0], pc[1] - px[1], pc[2] - px[2] };
				float v2[] = { pc[0] - py[0], pc[1] - py[1], pc[2] - py[2] };

				float v3[3];
				cross(v1, v2, v3);
				normalize(v3);
				outNorm[0] += v3[0];
				outNorm[1] += v3[1];
				outNorm[2] += v3[2];
				count++;
			}
			if (count) {
				float v3[3] = { outNorm[0] / count, outNorm[1] / count, outNorm[2] / count };
				normalize(v3);

				normals[3 * (i*width + j)] = v3[0];
				normals[3 * (i*width + j) + 1] = v3[1];
				normals[3 * (i*width + j) + 2] = v3[2];

				double d = -v3[0] * pc[0] - v3[1] * pc[1] - v3[2] * pc[2];
				int idx1 = RAD_SIZE + (int)(v3[0] * RAD_SIZE);
				int idx2 = RAD_SIZE + (int)(v3[1] * RAD_SIZE);
				assert(idx1 >= 0);
				assert(idx2 >= 0);
				assert(idx1 <RAD_FULL);
				assert(idx2 <RAD_FULL);

				planes[RAD_FULL*idx2 + idx1].push_back({ d, { v3[0], v3[1], v3[2] } });
			}
		}
	}
	std::vector<int> idx(RAD_FULL*RAD_FULL);
	for (int i = 0; i < RAD_FULL*RAD_FULL; i++)
		idx[i] = i;
	auto xySize = [&](int x, int y) { return planes[y*RAD_FULL + x].size();  };
	for (int i = 1; i < RAD_FULL - 1; i++) {
		for (int j = 1; j < RAD_FULL - 1; j++) {
			auto c = xySize(j, i);
			//if (c < MIN_PTS) idx[i*RAD_FULL + j]=0;

			if (c < xySize(j - 1, i - 1)) idx[i*RAD_FULL + j] = 0;
			if (c < xySize(j + 0, i - 1)) idx[i*RAD_FULL + j] = 0;
			if (c < xySize(j + 1, i - 1)) idx[i*RAD_FULL + j] = 0;
			if (c < xySize(j - 1, i + 0)) idx[i*RAD_FULL + j] = 0;
			if (c < xySize(j + 1, i + 0)) idx[i*RAD_FULL + j] = 0;
			if (c < xySize(j - 1, i + 1)) idx[i*RAD_FULL + j] = 0;
			if (c < xySize(j + 0, i + 1)) idx[i*RAD_FULL + j] = 0;
			if (c < xySize(j + 1, i + 1)) idx[i*RAD_FULL + j] = 0;
		}
	}

	std::partial_sort(std::begin(idx), std::begin(idx) + MAX_CAND, std::end(idx), [&](const int i1, const int i2) {
		return planes[i1].size() > planes[i2].size();
	});
	int realCand = 0;
	for (; realCand < MAX_CAND; realCand++) {
		if (planes[idx[realCand]].size() < MIN_PTS)
			break;
	}

	for (auto i = idx.begin(); i < (idx.begin() + realCand); i++) {
		std::vector<planeCandidate> candList;
		int currY = (*i / RAD_FULL);
		int currX = (*i % RAD_FULL);

		candList.insert(candList.end(), planes[currY*RAD_FULL + currX].begin(), planes[currY*RAD_FULL + currX].end());
		if(currY - 1 >= 0)		 candList.insert(candList.end(), planes[(currY-1)*RAD_FULL + currX].begin(), planes[(currY-1)*RAD_FULL + currX].end());
		if(currY + 1 < RAD_FULL) candList.insert(candList.end(), planes[(currY+1)*RAD_FULL + currX].begin(), planes[(currY+1)*RAD_FULL + currX].end());
		if(currX - 1 >= 0)		 candList.insert(candList.end(), planes[(currY)*RAD_FULL + currX-1].begin(), planes[(currY)*RAD_FULL + currX-1].end());
		if(currX + 1 < RAD_FULL) candList.insert(candList.end(), planes[(currY)*RAD_FULL + currX+1].begin(), planes[(currY)*RAD_FULL + currX+1].end());


		//Idea1: lets GMM this?

		//Idea2: sort lists by 'd'
		//Idea2: linear search, runnnig average, when startIdx < Mean -2sigma, stop. Repeat
		//Idea2: Use threshold instead of computing sigma

		//Idea3: be lazy
		planeCandidate avgCand = { 0, { 0, 0, 0 } };
		for (const auto c : candList) {
			avgCand.d += c.d;
			avgCand.n[0] += c.n[0];
			avgCand.n[1] += c.n[1];
			avgCand.n[2] += c.n[2];
		}
		avgCand.d /= candList.size();
		avgCand.n[0] /= candList.size();
		avgCand.n[1] /= candList.size();
		avgCand.n[2] /= candList.size();
		normalize(avgCand.n);
		returnData.push_back(std::move(avgCand));
		image[*i] = planes[*i].size();
	}

	//for (int i = 0; i < RAD_FULL*RAD_FULL; i++) {
	//	image[i] = planes[i].size();
	//}

	return returnData;
}

template <int size>
void generateNormals(const DepthPixel* depth, const int width, const int height, const float fx, const float fy, float* normals)
{
	const auto cX = 1.0f / fx;
	const auto cY = 1.0f / fy;
	const auto xStep = cX*size;
	const auto yStep = cY*size;
	auto halfX = width / 2;
	auto halfY = height / 2;
	const auto xyStep = xStep*yStep;
	for (int i = size; i < height - size; i++) {
		for (int j = size; j < width - size; j++) {
			if (!depth[i*width + j])		continue;
			const auto cDepth = depth[i*width + j];

			float outNorm[3] = { 0, 0, 0 };
			int count = 0;

			if (depth[i*width + j + size] && depth[(i + size)*width + j])
			{
				const auto xDepth = depth[i*width + j + size];
				const auto yDepth = depth[(i + size)*width + j];
				const auto diffXZ = xDepth - cDepth;
				const auto diffYZ = yDepth - cDepth;

				// cX*(5*z +j*(z-z2));
				const float v1[] = { cX*(size*xDepth + (j - halfX)*diffXZ), 0, diffXZ };
				const float v2[] = { 0, cY*(size*yDepth + (halfY - i)*diffYZ), diffYZ };
				float v3[] = { -v1[2] * v2[1], -v1[0] * v2[2], v1[0] * v2[1] };
				normalize(v3);
				outNorm[0] += v3[0];
				outNorm[1] += v3[1];
				outNorm[2] += v3[2];
				count++;
			}
			if (depth[i*width + j - size] && depth[(i + size)*width + j])
			{
				const auto xDepth = depth[i*width + j - size];
				const auto yDepth = depth[(i + size)*width + j];
				const auto diffXZ = cDepth - xDepth;
				const auto diffYZ = yDepth - cDepth;

				// cX*(5*z +j*(z-z2));
				const float v1[] = { cX*(size*xDepth + (j - halfX)*diffXZ), 0, diffXZ };
				const float v2[] = { 0, cY*(size*yDepth + (halfY - i)*diffYZ), diffYZ };
				float v3[] = { -v1[2] * v2[1], -v1[0] * v2[2], v1[0] * v2[1] };
				normalize(v3);
				outNorm[0] += v3[0];
				outNorm[1] += v3[1];
				outNorm[2] += v3[2];
				count++;
			}
			if (depth[i*width + j + size] && depth[(i - size)*width + j])
			{
				const auto xDepth = depth[i*width + j + size];
				const auto yDepth = depth[(i - size)*width + j];
				const auto diffXZ = xDepth - cDepth;
				const auto diffYZ = cDepth - yDepth;

				// cX*(5*z +j*(z-z2));
				const float v1[] = { cX*(size*xDepth + (j - halfX)*diffXZ), 0, diffXZ };
				const float v2[] = { 0, cY*(size*yDepth + (halfY - i)*diffYZ), diffYZ };
				float v3[] = { -v1[2] * v2[1], -v1[0] * v2[2], v1[0] * v2[1] };
				normalize(v3);
				outNorm[0] += v3[0];
				outNorm[1] += v3[1];
				outNorm[2] += v3[2];
				count++;
			}
			if (depth[i*width + j - size] && depth[(i - size)*width + j])
			{
				const auto xDepth = depth[i*width + j - size];
				const auto yDepth = depth[(i - size)*width + j];
				const auto diffXZ = cDepth - xDepth;
				const auto diffYZ = cDepth - yDepth;

				// cX*(5*z +j*(z-z2));
				const float v1[] = { cX*(size*xDepth + (j - halfX)*diffXZ), 0, diffXZ };
				const float v2[] = { 0, cY*(size*yDepth + (halfY - i)*diffYZ), diffYZ };
				float v3[] = { -v1[2] * v2[1], -v1[0] * v2[2], v1[0] * v2[1] };
				normalize(v3);
				outNorm[0] += v3[0];
				outNorm[1] += v3[1];
				outNorm[2] += v3[2];
				count++;
			}
			if (count) {
				float v3[3] = { outNorm[0] / count, outNorm[1] / count, outNorm[2] / count };
				normalize(v3);

				normals[3 * (i*width + j)] = v3[0];
				normals[3 * (i*width + j) + 1] = v3[1];
				normals[3 * (i*width + j) + 2] = v3[2];
			}
		}
	}
}

int main(int argc, char *args[]){
	Status rc = OpenNI::initialize();
	if (rc != STATUS_OK)
	{
		printf("Initialize failed\n%s\n", OpenNI::getExtendedError());
		return 1;
	}

	Device device;
	rc = device.open(ANY_DEVICE);
	if (rc != STATUS_OK)
	{
		printf("Couldn't open device\n%s\n", OpenNI::getExtendedError());
		return 2;
	}

	VideoStream depth;
	auto sensorInfo = device.getSensorInfo(SENSOR_DEPTH);
	auto desiredConfig = sensorInfo->getSupportedVideoModes()[0];

	if (sensorInfo != NULL)
	{

		rc = depth.create(device, SENSOR_DEPTH);
		rc = depth.setVideoMode(desiredConfig);
		if (rc != STATUS_OK)
		{
			printf("Couldn't create depth stream\n%s\n", OpenNI::getExtendedError());
			return 3;
		}
	}

	rc = depth.start();
	if (rc != STATUS_OK)
	{
		printf("Couldn't start the depth stream\n%s\n", OpenNI::getExtendedError());
		return 4;
	}

	VideoFrameRef frame;
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

	auto hfov = depth.getHorizontalFieldOfView();
	auto vfov = depth.getVerticalFieldOfView();

	auto fx = (width / 2) / tan(hfov / 2);
	auto fy = (height / 2) / tan(vfov / 2);

	while (true)
	{
		memset(normals.data(), 0, sizeof(float)*width*height * 3);
		memset(points.data(), 0, sizeof(float)*width*height * 3);

		int changedStreamDummy;
		VideoStream* pStream = &depth;
		rc = OpenNI::waitForAnyStream(&pStream, 1, &changedStreamDummy, 2000);
		if (rc != STATUS_OK)
		{
			printf("Wait failed! (timeout is %d ms)\n%s\n", 2000, OpenNI::getExtendedError());
			continue;
		}

		rc = depth.readFrame(&frame);
		if (rc != STATUS_OK)
		{
			printf("Read failed!\n%s\n", OpenNI::getExtendedError());
			continue;
		}

		if (frame.getVideoMode().getPixelFormat() != PIXEL_FORMAT_DEPTH_1_MM && frame.getVideoMode().getPixelFormat() != PIXEL_FORMAT_DEPTH_100_UM)
		{
			printf("Unexpected frame format\n");
			continue;
		}

		LARGE_INTEGER StartingTime, EndingTime, MiddleTime,ElapsedMicroseconds;
		LARGE_INTEGER Frequency;
		QueryPerformanceFrequency(&Frequency);

		DepthPixel* pDepth = (DepthPixel*)frame.getData();
		QueryPerformanceCounter(&StartingTime);
		generateHalfImage<DepthPixel, 4>(pDepth, hDepth, width, height);
		generatePoints(hDepth, width, height, fx/4, fy/4, points.data());
		QueryPerformanceCounter(&MiddleTime);

		static uint16_t ret[RAD_FULL * RAD_FULL];
		auto candidates = generateNormals<1>(points.data(), width, height, normals.data(),ret);
		//generateNormals<2>(pDepth, width, height, fx, fy, normals.data());
		QueryPerformanceCounter(&EndingTime);

		drawGrayScale(depthFrameSur, hDepth);
		drawGrayScale(visBinsSur, ret);
		drawNormals(normsFrameSur, normals.data());

		ElapsedMicroseconds.QuadPart = MiddleTime.QuadPart - StartingTime.QuadPart;
		ElapsedMicroseconds.QuadPart *= 1000000;
		ElapsedMicroseconds.QuadPart /= Frequency.QuadPart;

		printf("%lf ", static_cast<double>(ElapsedMicroseconds.QuadPart));

		ElapsedMicroseconds.QuadPart = EndingTime.QuadPart - MiddleTime.QuadPart;
		ElapsedMicroseconds.QuadPart *= 1000000;
		ElapsedMicroseconds.QuadPart /= Frequency.QuadPart;
		printf("%lf\n", static_cast<double>(ElapsedMicroseconds.QuadPart));

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

	depth.stop();
	depth.destroy();
	device.close();
	OpenNI::shutdown();

	SDL_DestroyTexture(bitmapTex);
	SDL_DestroyRenderer(renderer);
	SDL_DestroyWindow(win);

	return 0;
}
