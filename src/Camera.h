#pragma once
#include <stdint.h>
#include <math.h>
#include "linalg.h"
class Camera {
public:
	virtual bool Start() = 0;
	virtual bool syncNext() = 0;
	virtual uint16_t *getDepth() = 0;
	virtual uint16_t *getRGB() = 0;
	size_t getXDim() {
		return width;
	};
	size_t getYDim() {
		return height;
	};
	float getFx() {
		return fx;
	};
	float getFy() {
		return fy;
	};
	float getPx() {
		return px;
	};
	float getPy() {
		return py;
	};

protected:
	size_t width, height;
	float px, py, fx, fy;
};

void generatePoints(const uint16_t *depth, const int width, const int height, const float fx, const float fy, const float px, const float py, float *points) {
	auto halfX = px;
	auto halfY = py;
	auto cX = 1.0f / fx;
	auto cY = 1.0f / fy;
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			const auto z = depth[(i * width + j)];
			points[3 * (i * width + j)] = (j - halfX) * z * cX;
			points[3 * (i * width + j) + 1] = (i - halfY) * z * cY;
			points[3 * (i * width + j) + 2] = z;
		}
	}
}

template <int size>
void generateNormals(const uint16_t *depth, const int width, const int height, const float fx, const float fy, const float px, const float py, float *normals) {
	const auto cX = 1.0f / fx;
	const auto cY = 1.0f / fy;
	const auto xStep = cX * size;
	const auto yStep = cY * size;
	auto halfX = px;
	auto halfY = py;
	const auto xyStep = xStep * yStep;
	for (int i = size; i < height - size; i++) {
		for (int j = size; j < width - size; j++) {
			if (!depth[i * width + j])
				continue;
			const auto cDepth = depth[i * width + j];

			float outNorm[3] = { 0, 0, 0 };
			int count = 0;

			if (depth[i * width + j + size] && depth[(i + size) * width + j]) {
				const auto xDepth = depth[i * width + j + size];
				const auto yDepth = depth[(i + size) * width + j];
				const auto diffXZ = xDepth - cDepth;
				const auto diffYZ = yDepth - cDepth;

				// cX*(5*z +j*(z-z2));
				const float v1[] = { cX * (size * xDepth + (j - halfX) * diffXZ), 0, (float)diffXZ };
				const float v2[] = { 0, cY * (size * yDepth + (halfY - i) * diffYZ), (float)diffYZ };
				float v3[] = { -v1[2] * v2[1], -v1[0] * v2[2], v1[0] * v2[1] };
				normalize(v3);
				outNorm[0] += v3[0];
				outNorm[1] += v3[1];
				outNorm[2] += v3[2];
				count++;
			}
			if (depth[i * width + j - size] && depth[(i + size) * width + j]) {
				const auto xDepth = depth[i * width + j - size];
				const auto yDepth = depth[(i + size) * width + j];
				const auto diffXZ = cDepth - xDepth;
				const auto diffYZ = yDepth - cDepth;

				// cX*(5*z +j*(z-z2));
				const float v1[] = { cX * (size * xDepth + (j - halfX) * diffXZ), 0, (float)diffXZ };
				const float v2[] = { 0, cY * (size * yDepth + (halfY - i) * diffYZ), (float)diffYZ };
				float v3[] = { -v1[2] * v2[1], -v1[0] * v2[2], v1[0] * v2[1] };
				normalize(v3);
				outNorm[0] += v3[0];
				outNorm[1] += v3[1];
				outNorm[2] += v3[2];
				count++;
			}
			if (depth[i * width + j + size] && depth[(i - size) * width + j]) {
				const auto xDepth = depth[i * width + j + size];
				const auto yDepth = depth[(i - size) * width + j];
				const auto diffXZ = xDepth - cDepth;
				const auto diffYZ = cDepth - yDepth;

				// cX*(5*z +j*(z-z2));
				const float v1[] = { cX * (size * xDepth + (j - halfX) * diffXZ), 0, (float)diffXZ };
				const float v2[] = { 0, cY * (size * yDepth + (halfY - i) * diffYZ), (float)diffYZ };
				float v3[] = { -v1[2] * v2[1], -v1[0] * v2[2], v1[0] * v2[1] };
				normalize(v3);
				outNorm[0] += v3[0];
				outNorm[1] += v3[1];
				outNorm[2] += v3[2];
				count++;
			}
			if (depth[i * width + j - size] && depth[(i - size) * width + j]) {
				const auto xDepth = depth[i * width + j - size];
				const auto yDepth = depth[(i - size) * width + j];
				const auto diffXZ = cDepth - xDepth;
				const auto diffYZ = cDepth - yDepth;

				// cX*(5*z +j*(z-z2));
				const float v1[] = { cX * (size * xDepth + (j - halfX) * diffXZ), 0, (float)diffXZ };
				const float v2[] = { 0, cY * (size * yDepth + (halfY - i) * diffYZ), (float)diffYZ };
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

				normals[3 * (i * width + j)] = v3[0];
				normals[3 * (i * width + j) + 1] = v3[1];
				normals[3 * (i * width + j) + 2] = v3[2];
			}
		}
	}
}

template <int size>
void generateNormals_FromDepth(const uint16_t *depth, const int width, const int height, const float fx, const float fy, const float ppx, const float ppy, float *normals) {
	const auto cX = 1.0f / fx;
	const auto cY = 1.0f / fy;
	const auto xStep = cX * size;
	const auto yStep = cY * size;
	auto halfX = ppx;
	auto halfY = ppy;
	const auto xyStep = xStep * yStep;
	for (int i = size; i < height - size; i++) {
		for (int j = size; j < width - size; j++) {
			if (!depth[i * width + j])
				continue;
			const auto cDepth = depth[i * width + j];
			float pc[] = { cX*(j - halfX)*cDepth, cY*(i - halfY)*cDepth,cDepth };
			float outNorm[3] = { 0, 0, 0 };
			int count = 0;

			if (depth[i * width + j + size] && depth[(i + size) * width + j]) {
				const auto xDepth = depth[i * width + j + size];
				const auto yDepth = depth[(i + size) * width + j];
				const auto diffXZ = xDepth - cDepth;
				const auto diffYZ = yDepth - cDepth;
				float px[] = { cX*(j - halfX+size)*xDepth, cY*(i - halfY)*xDepth, xDepth };
				float py[] = { cX*(j - halfX)*yDepth, cY*(i - halfY + size)*yDepth, yDepth };

				// cX*(5*z +j*(z-z2));
				float v1[] = { px[0] - pc[0], px[1] - pc[1], px[2] - pc[2] };
				float v2[] = { py[0] - pc[0], py[1] - pc[1], py[2] - pc[2] };

				float v3[3];
				cross(v1, v2, v3);

				//possible optimization to unroll all operations
				//const float v11[] = { cX*(diffXZ*(j - halfX) + size*xDepth), cY*diffXZ*(i - halfY), (float)diffXZ };
				//const float v22[] = { cX*diffYZ*(j - halfX), cY*(diffYZ*(i - halfY) + size*yDepth), (float)diffYZ };
				//float v33[3];
				//cross(v11, v22, v33);

				normalize(v3);
				outNorm[0] += v3[0];
				outNorm[1] += v3[1];
				outNorm[2] += v3[2];
				count++;
			}
			if (depth[i * width + j - size] && depth[(i + size) * width + j]) {
				const auto xDepth = depth[i * width + j - size];
				const auto yDepth = depth[(i + size) * width + j];
				const auto diffXZ = xDepth - cDepth;
				const auto diffYZ = yDepth - cDepth;
				float px[] = { cX*(j - halfX - size)*xDepth, cY*(i - halfY)*xDepth, xDepth };
				float py[] = { cX*(j - halfX)*yDepth, cY*(i - halfY + size)*yDepth, yDepth };

				// cX*(5*z +j*(z-z2));
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
			if (depth[i * width + j + size] && depth[(i - size) * width + j]) {
				const auto xDepth = depth[i * width + j + size];
				const auto yDepth = depth[(i - size) * width + j];
				const auto diffXZ = xDepth - cDepth;
				const auto diffYZ = yDepth - cDepth;
				float px[] = { cX*(j - halfX + size)*xDepth, cY*(i - halfY)*xDepth, xDepth };
				float py[] = { cX*(j - halfX)*yDepth, cY*(i - halfY - size)*yDepth, yDepth };

				// cX*(5*z +j*(z-z2));
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
			if (depth[i * width + j - size] && depth[(i - size) * width + j]) {
				const auto xDepth = depth[i * width + j - size];
				const auto yDepth = depth[(i - size) * width + j];
				const auto diffXZ = xDepth - cDepth;
				const auto diffYZ = yDepth - cDepth;
				float px[] = { cX*(j - halfX - size)*xDepth, cY*(i - halfY)*xDepth, xDepth };
				float py[] = { cX*(j - halfX)*yDepth, cY*(i - halfY - size)*yDepth, yDepth };

				// cX*(5*z +j*(z-z2));
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

				normals[3 * (i * width + j)] = v3[0];
				normals[3 * (i * width + j) + 1] = v3[1];
				normals[3 * (i * width + j) + 2] = v3[2];
			}
		}
	}
}


template <int size>
void generateNormals_FromDepthInv(const uint16_t *depth, const int width, const int height, const float fx, const float fy, const float ppx, const float ppy, float *normals) {
	const auto cX = 1.0f / fx;
	const auto cY = 1.0f / fy;
	const auto xStep = cX * size;
	const auto yStep = cY * size;
	auto halfX = ppx;
	auto halfY = ppy;
	const auto xyStep = xStep * yStep;
	for (int i = size; i < height - size; i++) {
		for (int j = size; j < width - size; j++) {
			if (!depth[i * width + j])
				continue;
			const auto cDepth = static_cast<float>(UINT16_MAX)/depth[i * width + j];
			float pc[] = { cX*(j - halfX)*cDepth, cY*(i - halfY)*cDepth, cDepth };
			float outNorm[3] = { 0, 0, 0 };
			int count = 0;

			if (depth[i * width + j + size] && depth[(i + size) * width + j]) {
				const auto xDepth = static_cast<float>(UINT16_MAX) / depth[i * width + j + size];
				const auto yDepth = static_cast<float>(UINT16_MAX) / depth[(i + size) * width + j];
				const auto diffXZ = xDepth - cDepth;
				const auto diffYZ = yDepth - cDepth;
				float px[] = { cX*(j - halfX + size)*xDepth, cY*(i - halfY)*xDepth, xDepth };
				float py[] = { cX*(j - halfX)*yDepth, cY*(i - halfY + size)*yDepth, yDepth };

				// cX*(5*z +j*(z-z2));
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
			if (depth[i * width + j - size] && depth[(i + size) * width + j]) {
				const auto xDepth = static_cast<float>(UINT16_MAX) / depth[i * width + j - size];
				const auto yDepth = static_cast<float>(UINT16_MAX) / depth[(i + size) * width + j];
				const auto diffXZ = xDepth - cDepth;
				const auto diffYZ = yDepth - cDepth;
				float px[] = { cX*(j - halfX - size)*xDepth, cY*(i - halfY)*xDepth, xDepth };
				float py[] = { cX*(j - halfX)*yDepth, cY*(i - halfY + size)*yDepth, yDepth };

				// cX*(5*z +j*(z-z2));
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
			if (depth[i * width + j + size] && depth[(i - size) * width + j]) {
				const auto xDepth = static_cast<float>(UINT16_MAX) / depth[i * width + j + size];
				const auto yDepth = static_cast<float>(UINT16_MAX) / depth[(i - size) * width + j];
				const auto diffXZ = xDepth - cDepth;
				const auto diffYZ = yDepth - cDepth;
				float px[] = { cX*(j - halfX + size)*xDepth, cY*(i - halfY)*xDepth, xDepth };
				float py[] = { cX*(j - halfX)*yDepth, cY*(i - halfY - size)*yDepth, yDepth };

				// cX*(5*z +j*(z-z2));
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
			if (depth[i * width + j - size] && depth[(i - size) * width + j]) {
				const auto xDepth = static_cast<float>(UINT16_MAX) / depth[i * width + j - size];
				const auto yDepth = static_cast<float>(UINT16_MAX) / depth[(i - size) * width + j];
				const auto diffXZ = xDepth - cDepth;
				const auto diffYZ = yDepth - cDepth;
				float px[] = { cX*(j - halfX - size)*xDepth, cY*(i - halfY)*xDepth, xDepth };
				float py[] = { cX*(j - halfX)*yDepth, cY*(i - halfY - size)*yDepth, yDepth };

				// cX*(5*z +j*(z-z2));
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

				normals[3 * (i * width + j)] = -v3[0];
				normals[3 * (i * width + j) + 1] = -v3[1];
				normals[3 * (i * width + j) + 2] = v3[2];
			}
		}
	}
}

template <int size>
void generateNormals_fromPoints(const float *points, const int width, const int height, float *normals) {
	for (int i = size; i < height - size; i++) {
		for (int j = size; j < width - size; j++) {
			if (!points[3 * (i * width + j) + 2])
				continue;
			const float *pc = &points[3 * (i * width + j)];

			float outNorm[3] = { 0, 0, 0 };
			int count = 0;
			if (points[3 * (i * width + j + size) + 2] && points[3 * ((i + size) * width + j) + 2]) {
				const float *px = &points[3 * (i * width + j + size)];
				const float *py = &points[3 * ((i + size) * width + j)];

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
			if (points[3 * (i * width + j - size) + 2] && points[3 * ((i + size) * width + j) + 2]) {
				const float *px = &points[3 * (i * width + j - size)];
				const float *py = &points[3 * ((i + size) * width + j)];

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
			if (points[3 * (i * width + j + size) + 2] && points[3 * ((i - size) * width + j) + 2]) {
				const float *px = &points[3 * (i * width + j + size)];
				const float *py = &points[3 * ((i - size) * width + j)];

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
			if (points[3 * (i * width + j - size) + 2] && points[3 * ((i - size) * width + j) + 2]) {
				const float *px = &points[3 * (i * width + j - size)];
				const float *py = &points[3 * ((i - size) * width + j)];

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

				normals[3 * (i * width + j)] = v3[0];
				normals[3 * (i * width + j) + 1] = v3[1];
				normals[3 * (i * width + j) + 2] = v3[2];

				//double d = -v3[0] * pc[0] - v3[1] * pc[1] - v3[2] * pc[2];
			}
		}
	}
}
