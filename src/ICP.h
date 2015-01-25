#pragma once
#include "linalg.h"
#include <vector>

void transformImage(
    const int width, const int height,
    const float fx, const float fy,
    const float *depthSrc, const float *normalsSrc,
    float *depthDst, float *normalsDst,
    float rotation[9], float translation[3]) {
    using namespace Eigen;
    Matrix3f rM(rotation);
    Vector3f tV(translation);
    for (int i = 0; i < width * height; i++) {
        const auto p = Vector3f(depthSrc + 3 * i);
        const auto n = Vector3f(normalsSrc + 3 * i);

        const auto p_new = rM * p + tV;
        const auto n_new = rM * n;

        const auto x = static_cast<int>(p_new[0] * fx / p_new[2] + width / 2 + 0.5f);
        const auto y = static_cast<int>(p_new[1] * fy / p_new[2] + height / 2 + 0.5f);
        if (x < 0 || x >= width)
            continue;
        if (y < 0 || y >= height)
            continue;
        const auto existingZ = depthDst[3 * (y * width + x) + 2];
        if (existingZ != 0 && existingZ < p_new[2])
            continue;
        depthDst[3 * (y * width + x)] = p_new[0];
        depthDst[3 * (y * width + x) + 1] = p_new[1];
        depthDst[3 * (y * width + x) + 2] = p_new[2];
        normalsDst[3 * (y * width + x)] = n_new[0];
        normalsDst[3 * (y * width + x) + 1] = n_new[1];
        normalsDst[3 * (y * width + x) + 2] = n_new[2];
    }
}

void computeLinearApproxICP(
    const int width, const int height,
    const float normThresh, const float distThresh,
    const float *depthSrc, const float *normalsSrc,
    const float *depthDst, const float *normalsDst) {
    using namespace Eigen;
    using namespace std;
    vector<int> goodIndex;
    for (int i = 0; i < width * height; i++) {
        if (depthSrc[3 * i + 2] == 0)
            continue;
        if (depthDst[3 * i + 2] == 0)
            continue;
        if (normalsSrc[3 * i + 2] == 0)
            continue;
        if (normalsDst[3 * i + 2] == 0)
            continue;
        if (square((int)depthSrc[3 * i + 2] - (int)depthSrc[3 * i + 2]) > distThresh)
            continue;
        if (sqrNorm(normalsSrc + 3 * i, normalsDst + 3 * i) > normThresh)
            continue;

        goodIndex.push_back(i);
    }

    if (goodIndex.size() > 100) {

        //LARGE_INTEGER StartingTime, EndingTime, MiddleTime,ElapsedMicroseconds;
        //LARGE_INTEGER Frequency;
        //QueryPerformanceFrequency(&Frequency);
        //QueryPerformanceCounter(&StartingTime);

        //{
        //	//Matrix<float, Dynamic, Dynamic, RowMajor>
        //	MatrixXf A(goodIndex.size(), 6);
        //	VectorXf b(goodIndex.size());

        //	for (int i = 0; i < goodIndex.size(); i++) {
        //		const auto idx = goodIndex[i];
        //		A(i, 0) = normalsDst[3 * idx + 2] * depthSrc[3 * idx + 1] - normalsDst[3 * idx + 1] * depthSrc[3 * idx + 2];
        //		A(i, 1) = normalsDst[3 * idx + 0] * depthSrc[3 * idx + 2] - normalsDst[3 * idx + 2] * depthSrc[3 * idx + 0];
        //		A(i, 2) = normalsDst[3 * idx + 1] * depthSrc[3 * idx + 0] - normalsDst[3 * idx + 0] * depthSrc[3 * idx + 1];
        //		A(i, 3) = normalsDst[3 * idx + 0];
        //		A(i, 4) = normalsDst[3 * idx + 1];
        //		A(i, 5) = normalsDst[3 * idx + 2];

        //		b(i) =	  normalsDst[3 * idx + 0] * (depthDst[3 * idx + 0] - depthSrc[3 * idx + 0])
        //				+ normalsDst[3 * idx + 1] * (depthDst[3 * idx + 1] - depthSrc[3 * idx + 1])
        //				+ normalsDst[3 * idx + 2] * (depthDst[3 * idx + 2] - depthSrc[3 * idx + 2]);

        //	}
        //	JacobiSVD<MatrixXf> svd(A, ComputeThinU | ComputeThinV);
        //	VectorXf x = svd.solve(b);

        //	Matrix3f rotation;
        //	rotation << 1, x[0] * x[1] - x[2], x[0] * x[2] + x[1],
        //		x[2], x[0] * x[1] * x[2] + 1, x[1] * x[2] - x[0],
        //		-x[1], x[0], 1;
        //	Vector3f translation(x[3], x[4], x[5]);
        //	//std::cout << rotation << std::endl << translation << std::endl;

        //}
        //QueryPerformanceCounter(&MiddleTime);

        {
            float A[6 * 6] = { 0 };
            float b[6] = { 0 };
            float covar[6] = { 0 };

            float cPt[3] = { 0 };

            for (int i = 0; i < goodIndex.size(); i++) {
                const auto idx = goodIndex[i];

                cPt[0] = normalsDst[3 * idx + 2] * depthSrc[3 * idx + 1] - normalsDst[3 * idx + 1] * depthSrc[3 * idx + 2];
                cPt[1] = normalsDst[3 * idx + 0] * depthSrc[3 * idx + 2] - normalsDst[3 * idx + 2] * depthSrc[3 * idx + 0];
                cPt[2] = normalsDst[3 * idx + 1] * depthSrc[3 * idx + 0] - normalsDst[3 * idx + 0] * depthSrc[3 * idx + 1];

                float diff = normalsDst[3 * idx + 0] * (depthSrc[3 * idx + 0] - depthDst[3 * idx + 0]) + normalsDst[3 * idx + 1] * (depthSrc[3 * idx + 1] - depthDst[3 * idx + 1]) + normalsDst[3 * idx + 2] * (depthSrc[3 * idx + 2] - depthDst[3 * idx + 2]);

                covar[0] = cPt[0];
                covar[1] = cPt[1];
                covar[2] = cPt[2];
                covar[3] = normalsDst[3 * idx];
                covar[4] = normalsDst[3 * idx + 1];
                covar[5] = normalsDst[3 * idx + 2];

                for (int i = 0; i < 6; i++) {
                    for (int j = 0; j < 6; j++) {
                        A[i * 6 + j] += covar[i] * covar[j];
                    }
                }

                b[0] += cPt[0] * diff;
                b[1] += cPt[1] * diff;
                b[2] += cPt[2] * diff;
                b[3] += normalsDst[3 * idx] * diff;
                b[4] += normalsDst[3 * idx + 1] * diff;
                b[5] += normalsDst[3 * idx + 2] * diff;
            }
            MatrixXf AMat = Map<MatrixXf, 0, InnerStride<0> >(A, 6, 6, InnerStride<0>());
            VectorXf bMat = Map<VectorXf, 0, InnerStride<0> >(b, 6);
            //std::cout << AMat << std::endl << bMat << std::endl;
            LDLT<MatrixXf> ch(AMat);

            VectorXf x = ch.solve(-bMat);

            Matrix3f rotation;
            rotation << 1, x[0] * x[1] - x[2], x[0] * x[2] + x[1],
                x[2], x[0] * x[1] * x[2] + 1, x[1] * x[2] - x[0],
                -x[1], x[0], 1;
            Vector3f translation(x[3], x[4], x[5]);
            //std::cout << x << std::endl;

            std::cout << rotation << std::endl << translation << std::endl;
        }
        //QueryPerformanceCounter(&EndingTime);
        //ElapsedMicroseconds.QuadPart = MiddleTime.QuadPart - StartingTime.QuadPart;
        //ElapsedMicroseconds.QuadPart *= 1000000;
        //ElapsedMicroseconds.QuadPart /= Frequency.QuadPart;

        //printf("%lf ", static_cast<double>(ElapsedMicroseconds.QuadPart));

        //ElapsedMicroseconds.QuadPart = EndingTime.QuadPart - MiddleTime.QuadPart;
        //ElapsedMicroseconds.QuadPart *= 1000000;
        //ElapsedMicroseconds.QuadPart /= Frequency.QuadPart;
        //printf("%lf\n", static_cast<double>(ElapsedMicroseconds.QuadPart));
    }
}
