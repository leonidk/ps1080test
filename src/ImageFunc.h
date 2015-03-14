#pragma once
#include <stdio.h>
#include <type_traits>
#include <memory>
#include <array>


template <typename T, int size>
inline void generateHalfImage(const T *in, T *out, const int outW, const int outH) {
    const auto inH = size * outH;
    const auto inW = size * outW;
    for (int i = 0; i < inH; i += size) {
        for (int j = 0; j < inW; j += size) {
            auto tout = 0;
            for (int m = 0; m < size; m++) {
                for (int n = 0; n < size; n++) {
                    tout += in[(i + m) * inW + (j + n)];
                }
            }
            out[(i / size) * outW + (j / size)] = tout / (size * size);
        }
    }
}


template <typename T, int size>
inline void generateHalfImageDepth(const T *in, T *out, const int outW, const int outH) {
	const auto inH = size * outH;
	const auto inW = size * outW;
	for (int i = 0; i < inH; i += size) {
		for (int j = 0; j < inW; j += size) {
			auto tout = 0;
			auto cnt = 0;
			for (int m = 0; m < size; m++) {
				for (int n = 0; n < size; n++) {
					const auto val = in[(i + m) * inW + (j + n)];
					tout += val ? val : 0;
					cnt += val ? 1 : 0;
				}
			}
			out[(i / size) * outW + (j / size)] = cnt ? (int)std::round(tout / ((float)cnt)) : 0 ;
		}
	}
}

//in units of z
template <typename T, int size>
inline void generateHalfImageDepthVarMed(const T *in, T *out, const int outW, const int outH) {
	const auto inH = size * outH;
	const auto inW = size * outW;
	for (int i = 0; i < inH; i += size) {
		for (int j = 0; j < inW; j += size) {
			auto cnt = 0;
			float mean = 0;
			float M2 = 0;
			//loop over candidates points, computing mean and variance
			T tmp[size*size] = { 0 };
			for (int m = 0; m < size; m++) {
				for (int n = 0; n < size; n++) {
					const auto val = in[(i + m) * inW + (j + n)];
					if (val) {
						tmp[cnt++] = val;
						float delta = val - mean;
						mean += delta / cnt;
						M2 += delta*(val - mean);
					}
				}
			}
			//insertion sort
			for (int c = 1; c < (cnt ); c++) {
				auto _x = tmp[c];
				auto _j = c;
				while (_j > 0 && tmp[_j - 1] > _x) {
					tmp[_j] = tmp[_j - 1];
					_j = _j - 1;
				}
				tmp[_j] = _x;
			}
			int ret = 0;
			if (cnt) {
				auto var = M2 / (cnt - 1);
				//int val_int = *(int*)&var;
				//val_int = (1 << 29) + (val_int >> 1) - (1 << 22) + -0x4C000;
				//float stdDev = *(float*)&val_int;
				auto stdDev = std::sqrt(var)+0.5;
				ret = abs(mean - tmp[cnt / 2]) < stdDev/2 ? (int)(mean+0.5) : tmp[cnt / 2];
			}

			out[(i / size) * outW + (j / size)] = ret;
		}
	}
}
// An idea on subsampling depth maps: sort+split in two. Compute Mean+Var of both halves. Compare MeanA+DevA to MeanB-DevB with "Reasonable differences" 
// (perhaps a var of uniform or parabolic distribution) 
template <typename T, int size>
inline void generateHalfImageDepthMed(const T *in, T *out, const int outW, const int outH) {
	const auto inH = size * outH;
	const auto inW = size * outW;
	for (int i = 0; i < inH; i += size) {
		for (int j = 0; j < inW; j += size) {
			auto cnt = 0;
			T tmp[size*size] = { 0 };
			for (int m = 0; m < size; m++) {
				for (int n = 0; n < size; n++) {
					const auto val = in[(i + m) * inW + (j + n)];
					if (val) {
						tmp[cnt++] = val;
					}
				}
			}
			//insertion sort
			for (int c = 1; c < (cnt ); c++) {
				auto _x = tmp[c];
				auto _j = c;
				while (_j > 0 && tmp[_j - 1] > _x) {
					tmp[_j] = tmp[_j - 1];
					_j = _j - 1;
				}
				tmp[_j] = _x;
			}
			out[(i / size) * outW + (j / size)] = tmp[cnt/2];
		}
	}
}

template <typename T, int size>
inline void generateHalfImageDepthTMean(const T *in, T *out, const int outW, const int outH) {
	const auto inH = size * outH;
	const auto inW = size * outW;
	for (int i = 0; i < inH; i += size) {
		for (int j = 0; j < inW; j += size) {
			auto cnt = 0;
			T tmp[size*size] = { 0 };
			for (int m = 0; m < size; m++) {
				for (int n = 0; n < size; n++) {
					const auto val = in[(i + m) * inW + (j + n)];
					if (val) {
						tmp[cnt++] = val;
					}
				}
			}
			//insertion sort
			for (int c = 1; c < (cnt - 1); c++) {
				auto _x = tmp[c];
				auto _j = c;
				while (_j > 0 && tmp[_j - 1] > _x) {
					tmp[_j] = tmp[_j - 1];
					_j = _j - 1;
				}
				tmp[_j] = _x;
			}
			auto sum = 0;
			auto lBnd = cnt/8;
			auto rBnd = (cnt*7)/8;
			auto div = (rBnd - lBnd);
			for (int c = lBnd; c < rBnd; c++)
				sum += tmp[c];
			out[(i / size) * outW + (j / size)] = (sum + div/2) / (div ? div : 1);
		}
	}
}

template <typename T, int size>
inline void generateHalfImageDepthInv(const T *in, T *out, const int outW, const int outH) {
	const auto inH = size * outH;
	const auto inW = size * outW;
	for (int i = 0; i < inH; i += size) {
		for (int j = 0; j < inW; j += size) {
			float tout = 0;
			auto cnt = 0;
			for (int m = 0; m < size; m++) {
				for (int n = 0; n < size; n++) {
					const auto val = in[(i + m) * inW + (j + n)];
					tout += val ? 1.0f / val : 0;
					cnt += val ? 1 : 0;
				}
			}
			out[(i / size) * outW + (j / size)] = cnt ? (int)std::round(((float)cnt)/tout) : 0;
		}
	}
}


inline void generatePoints(const uint16_t *depth, const int width, const int height, const float fx, const float fy, const float px, const float py, float *points) {
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
inline void generateNormals(const uint16_t *depth, const int width, const int height, const float fx, const float fy, const float px, const float py, float *normals) {
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
inline void generateNormals_FromDepth(const uint16_t *depth, const int width, const int height, const float fx, const float fy, const float ppx, const float ppy, float *normals) {
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
			float pc[] = { cX*(j - halfX)*cDepth, cY*(i - halfY)*cDepth, cDepth };
			float outNorm[3] = { 0, 0, 0 };
			int count = 0;

			if (depth[i * width + j + size] && depth[(i + size) * width + j]) {
				const auto xDepth = depth[i * width + j + size];
				const auto yDepth = depth[(i + size) * width + j];
				const auto diffXZ = xDepth - cDepth;
				const auto diffYZ = yDepth - cDepth;
				float px[] = { cX*(j - halfX + size)*xDepth, cY*(i - halfY)*xDepth, xDepth };
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
inline  void generateNormals_FromDepthInv(const uint16_t *depth, const int width, const int height, const float fx, const float fy, const float ppx, const float ppy, float *normals) {
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
			const auto cDepth = static_cast<float>(UINT16_MAX) / depth[i * width + j];
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
inline void generateNormals_fromPoints(const float *points, const int width, const int height, float *normals) {
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

template <
	template <typename, typename> class Container,
	typename Value,
	typename Allocator = std::allocator<Value>, typename Func >
	inline auto imap(const Container<Value, Allocator> & input, const Func &f)
	-> Container<decltype(f(std::declval<Value>())), std::allocator<decltype(f(std::declval<Value>()))>> {
	Container<decltype(f(std::declval<Value>())), std::allocator<decltype(f(std::declval<Value>()))>> ret;
	std::transform(std::begin(input), std::end(input), std::back_inserter(ret), f);
	//for (const auto & v : input) {
	//	ret.emplace_back(f(v));
	//}
	return ret;
}

template <
	template <typename, typename> class Container,
	typename Value,
	typename Allocator = std::allocator<Value>,
	typename Func,
	typename OutType = std::result_of<Func(Value)>::type >
	inline auto imap2(const Container<Value, Allocator> & input, const Func &f)
	->Container<OutType, std::allocator<OutType>> {
	Container<OutType, std::allocator<OutType>> ret;
	std::transform(std::begin(input), std::end(input), std::back_inserter(ret), f);
	//for (const auto & v : input) {
	//	ret.emplace_back(f(v));
	//}
	return ret;
}

static std::array<uint8_t, 3> HsvToRgb(std::array<uint8_t, 3> hsv)
{
	std::array<uint8_t, 3> rgb;
	unsigned char region, remainder, p, q, t;

	if (hsv[1] == 0)
	{
		rgb[0] = hsv[2];
		rgb[1] = hsv[2];
		rgb[2] = hsv[2];
		return rgb;
	}

	region = hsv[0] / 43;
	remainder = (hsv[0] - (region * 43)) * 6;

	p = (hsv[2] * (255 - hsv[1])) >> 8;
	q = (hsv[2] * (255 - ((hsv[1] * remainder) >> 8))) >> 8;
	t = (hsv[2] * (255 - ((hsv[1] * (255 - remainder)) >> 8))) >> 8;

	switch (region)
	{
	case 0:
		rgb[0] = hsv[2]; rgb[1] = t; rgb[2] = p;
		break;
	case 1:
		rgb[0] = q; rgb[1] = hsv[2]; rgb[2] = p;
		break;
	case 2:
		rgb[0] = p; rgb[1] = hsv[2]; rgb[2] = t;
		break;
	case 3:
		rgb[0] = p; rgb[1] = q; rgb[2] = hsv[2];
		break;
	case 4:
		rgb[0] = t; rgb[1] = p; rgb[2] = hsv[2];
		break;
	default:
		rgb[0] = hsv[2]; rgb[1] = p; rgb[2] = q;
		break;
	}

	return rgb;
}

static std::array<uint8_t, 3> getNewColor(int value, int total = 10) {
	return HsvToRgb({ (uint8_t)(255.0f*(static_cast<float>(value) / total)), 128, 255 });
}