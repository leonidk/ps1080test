#pragma once
#include <stdint.h>
#include <math.h>

class Camera {
public:
	Camera();
	~Camera();
	//Camera(size_t w, size_t h, size_t t);
	virtual bool		Start() = 0;
	virtual bool		syncNext() = 0;
	virtual uint16_t*	getDepth() = 0;
	virtual uint16_t*	getRGB() = 0;
	size_t		getXDim() { return width; };
	size_t		getYDim() { return height; };
	size_t		getFx() { return fx; };
	size_t		getFy() { return fy; };
	size_t		getPx() { return px; };
	size_t		getPy() { return py; };
protected:
	size_t width, height;
	float px, py, fx, fy;
};

void generatePoints(const uint16_t* depth, const int width, const int height, const float fx, const float fy, float* points)
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