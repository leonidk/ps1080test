#include <stdio.h>
#include <OpenNI.h>
#include <SDL.h>
#include <math.h>
#include <vector>
using namespace openni;

template <typename T>
void drawGrayScale(SDL_Surface* sur, T * data) {
	uint8_t* dest = (uint8_t*)sur->pixels;
	for (int i = 0; i < sur->h*sur->w; i++)
		dest[4 * i + 2] = dest[4 * i + 1] = dest[4 * i] = (data[i] >> 4) > 255 ? 255 : (data[i] >> 4);
}

void drawNormals(SDL_Surface* sur, float* normals) {
	uint8_t* dest = (uint8_t*)sur->pixels;
	for (int i = 0; i < sur->h*sur->w; i++)
	{
		dest[4 * i] = 128+normals[3 * i] * 127;
		dest[4 * i + 1] = 128+normals[3 * i + 1] * 127;
		dest[4 * i + 2] = 128+normals[3 * i + 2] * 127;

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
			points[3 * (i*width+j)]		= (j - halfX)*z*cX;
			points[3 * (i*width+j) + 1] = (i-halfY)*z*cY;
			points[3 * (i*width+j) + 2] = z;
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

template <int size>
void generateNormals(const float* points, const int width, const int height, float* normals)
{
	for (int i = size; i < height - size; i++) {
		for (int j = size; j < width - size; j++) {
			if (!points[3 * (i*width + j) + 2])		continue;
			if (!points[3 * (i*width + j + size) + 2])	continue;
			if (!points[3 * ((i + size)*width + j) + 2])	continue;
		
			const float* pc = &points[3 * (i*width + j)];
			const float* px = &points[3 * (i*width + j + size)];
			const float* py = &points[3 * ((i + size)*width + j)];

			float v1[] = { px[0] - pc[0], px[1] - pc[1], px[2] - pc[2] };
			float v2[] = { py[0] - pc[0], py[1] - pc[1], py[2] - pc[2] };
			
			float v3[3];
			cross(v1, v2, v3);
			normalize(v3);
			normals[3 * (i*width + j)] = v3[0];
			normals[3 * (i*width + j) + 1] = v3[1];
			normals[3 * (i*width + j) + 2] = v3[2];
		}
	}
}

template <int size>
void generateNormals(const DepthPixel* depth, const int width, const int height, const float fx, const float fy,float* normals)
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
			if (!depth[i*width + j+size])	continue;
			if (!depth[(i+size)*width + j])	continue;
			const auto cDepth = depth[i*width + j];
			const auto xDepth = depth[i*width + j + size];
			const auto yDepth = depth[(i + size)*width + j];
			const auto diffXZ = xDepth - cDepth;
			const auto diffYZ = yDepth - cDepth;

			//const float v1[] = { xStep, 0, diffXZ };
			//const float v2[] = { 0, yStep, diffYZ };
			//float v3[] = { -diffXZ*yStep, -xStep*diffYZ, xyStep };

			//const float v1[] = { cX*(size*xDepth), 0, diffXZ };
			//const float v2[] = { 0, cY*(size*yDepth), diffYZ };
			//float v3[] = { -diffXZ*yStep*yDepth, -xStep*diffYZ*xDepth, xyStep*xDepth*yDepth };

			// cX*(5*z +j*(z-z2));
			const float v1[] = { cX*(size*xDepth+(j-halfX)*diffXZ), 0, diffXZ };
			const float v2[] = { 0, cY*(size*yDepth + (halfY - i)*diffYZ), diffYZ };
			float v3[] = { -v1[2] * v2[1], -v1[0] * v2[2], v1[0] * v2[1] };
			normalize(v3);
			normals[3 * (i*width + j)]	 = v3[0];
			normals[3 * (i*width + j)+1] = v3[1];
			normals[3 * (i*width + j)+2] = v3[2];
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

	SDL_Texture *bitmapTex = NULL;
	SDL_Surface *bitmapSurface = NULL;
	int posX = 100, posY = 100, width = 320, height = 240;

	win = SDL_CreateWindow("Hello World", posX, posY, width, height, 0);
	win2 = SDL_CreateWindow("Hello World2", posX + width, posY, width, height, 0);

	renderer = SDL_CreateRenderer(win, -1, SDL_RENDERER_ACCELERATED);
	renderer2 = SDL_CreateRenderer(win2, -1, SDL_RENDERER_ACCELERATED);

	SDL_Surface *depthFrameSur = SDL_CreateRGBSurface(0, width, height, 32, 0, 0, 0, 0);
	SDL_Surface *normsFrameSur = SDL_CreateRGBSurface(0, width, height, 32, 0, 0, 0, 0);

	std::vector<float> points (3 * width*height,0);
	std::vector<float> normals(3 * width*height,0);
	
	auto hfov = depth.getHorizontalFieldOfView();
	auto vfov = depth.getVerticalFieldOfView();

	auto fx = (width/2)/tan(hfov / 2);
	auto fy = (height / 2)/tan(vfov / 2);

	while (true)
	{
		memset(normals.data(), 0, sizeof(float)*width*height * 3);
		memset(points.data(), 0, sizeof(float)*width*height * 3);

		int changedStreamDummy;
		VideoStream* pStream = &depth;
		rc = OpenNI::waitForAnyStream (&pStream, 1, &changedStreamDummy, 2000);
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

		DepthPixel* pDepth = (DepthPixel*)frame.getData();
		drawGrayScale(depthFrameSur, pDepth);
		//generatePoints(pDepth, width, height, fx, fy, points.data());
		//generateNormals<5>(points.data(), width, height, normals.data());

		generateNormals<5>(pDepth, width, height, fx, fy, normals.data());
		drawNormals(normsFrameSur, normals.data());
		//memcpy(depthFrameSur->pixels, (void*)pDepth, 320 * 240 * sizeof(DepthPixel));
		int middleIndex = (frame.getHeight() + 1)*frame.getWidth() / 2;

		printf("[%08llu] %8d\n", (long long)frame.getTimestamp(), pDepth[middleIndex]);
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
