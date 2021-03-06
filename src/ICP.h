#pragma once
#include <vector>
#include <memory>
#include <iostream>
#include "linalg.h"
#include <Eigen/Geometry>
#include <Eigen/Cholesky>


void transformImage(
	const int width, const int height,
	const float fx, const float fy,
	const float *depthSrc, const float *normalsSrc,
	float *depthDst, float *normalsDst,
	const Eigen::Matrix4f & transformation) {
	using namespace Eigen;

	Matrix3f rM = transformation.block(0, 0, 3, 3);
	Vector3f tV = transformation.block(0, 3, 3, 1);
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

Eigen::Matrix4f computeLinearApproxICP(
	const int width, const int height,
	const float normThresh, const float distThresh, const size_t ptThresh,
	const float *depthSrc, const float *normalsSrc,
	const float *depthDst, const float *normalsDst) {
	using namespace Eigen;
	using namespace std;
	vector<int> goodIndex;
	for (int i = 0; i < width * height; i++)
	{
		if (depthSrc[3 * i + 2] == 0)
			continue;
		if (depthDst[3 * i + 2] == 0)
			continue;
		if (normalsSrc[3 * i + 2] == 0)
			continue;
		if (normalsDst[3 * i + 2] == 0)
			continue;
		if (square(depthSrc[3 * i + 2] - depthSrc[3 * i + 2]) > distThresh)
			continue;
		if (sqrNorm(normalsSrc + 3 * i, normalsDst + 3 * i) > normThresh)
			continue;

		goodIndex.push_back(i);
	}

	if (goodIndex.size() > ptThresh) {

		// solving the linear system with an SVD
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

		// solving the linear system with normal equations
		{
			float A[6 * 6] = { 0 };
			float b[6] = { 0 };
			float covar[6] = { 0 };

			float cPt[3] = { 0 };

			for (size_t i = 0; i < goodIndex.size(); i++) {
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
#define I_DONT_LIKE_EIGEN
#ifdef I_DONT_LIKE_EIGEN
			//equiv to LLT from Eigen. Adapted from dlib
			float L[6 * 6] = { 0 };
			float LT[6 * 6] = { 0 };

			float z[6] = { 0 };

			float xA[6] = { 0 };

			const float eps = 1e-3;
			bool isspd = true;

			// Main loop.
			for (int j = 0; j < 6; j++)	{
				float d(0.0);
				for (int k = 0; k < j; k++)	{
					float s(0.0);
					for (int i = 0; i < k; i++)	{
						s += L[6*k+ i]*L[6*j+i];
					}

					// if L_(k,k) != 0
					if (std::abs(L[6*k+k]) > eps)
					{
						s = (A[j * 6 + k] - s) / L[k * 6 + k];
					} else	{
						s = (A[j * 6 + k] - s);
						isspd = false;
					}
					L[j*6+k] = s;
					d = d + s*s;
					// this is approximately doing: isspd = isspd && ( A(k,j) == A(j,k))
					isspd = isspd && (std::abs(A[k*6+j] - A[j*6+k]) < eps);
				}
				d = A[j*6+j] - d;
				isspd = isspd && (d > eps);
				L[j*6+j] = sqrt(d > 0.0 ? d : 0.0);
				for (int k = j + 1; k < 6; k++)
				{
					L[j*6+ k] = 0.0;
				}
			}
			//forw sub
			for (int i = 0; i < 6; i++)	{
				z[i] = b[i];
				for (int j = 0; j < i; j++)
				{
					z[i] -= L[i * 6 + j] * z[j];
				}
				z[i] /= L[i*6+i];
			}

			//transpose
			for (int i = 0; i < 6; i++) {
				for (int j = 0; j < 6; j++) {
					LT[i * 6 + j] = L[j * 6 + i];
				}
			}

			//back sub
			for (int i = 6 - 1; i >= 0; i--) {
				xA[i] = z[i];
				for (int j = i + 1; j < 6; j++)	{
					xA[i] -= LT[i * 6 + j] * xA[j];
				}
				xA[i] /= LT[i * 6 + i];
			}
#endif
			MatrixXf AMat = Map<MatrixXf, 0, InnerStride<0> >(A, 6, 6, InnerStride<0>());
			VectorXf bMat = Map<VectorXf, 0, InnerStride<0> >(b, 6);
			LDLT<MatrixXf> ch(AMat); 
			MatrixXf L2 = ch.matrixU();
			VectorXf x = ch.solve(bMat);

			//Matrix3f rotation;
			//rotation << 1, x[0] * x[1] - x[2], x[0] * x[2] + x[1],
			//    x[2], x[0] * x[1] * x[2] + 1, x[1] * x[2] - x[0],
			//    -x[1], x[0], 1;
			//Vector3f translation(x[3], x[4], x[5]);
			// std::cout << rotation << std::endl << translation << std::endl;


			Matrix4f transform;
			transform << 1, x[0] * x[1] - x[2], x[0] * x[2] + x[1], x[3],
				x[2], x[0] * x[1] * x[2] + 1, x[1] * x[2] - x[0], x[4],
				-x[1], x[0], 1, x[5],
				0, 0, 0, 1;
			return transform;
		}

	}
	Affine3f transform(Translation3f(0, 0, 0));
	Matrix4f matrix = transform.matrix();
	return matrix;
}

class iteratedICP
{
public:
	iteratedICP(const int width, const int height)
		: width(width), height(height), dataPtr(std::unique_ptr<float[]>(new float[3 * 2 * width*height]))
	{
		points = static_cast<float*>(dataPtr.get());
		normals = static_cast<float*>(dataPtr.get()) + 3 * width*height;
	}
	Eigen::Matrix4f runICPIter(const int iter, const float fx, const float fy,
		const float normThresh, const float distThresh, const size_t ptThresh,
		const float *depthSrc, const float *normalsSrc,
		const float *depthDst, const float *normalsDst)
	{
		using namespace Eigen;
		if (iter < 1) {
			Affine3f transform(Translation3f(0, 0, 0));
			Matrix4f matrix = transform.matrix();
			return matrix;
		}

		Eigen::Transform<float, 3, Affine> trans(computeLinearApproxICP(width, height, normThresh, distThresh, ptThresh, depthSrc, normalsSrc, depthDst, normalsDst));
		for (int i = 1; i < iter; ++i) {
			memset(dataPtr.get(), 0, 6 * width*height*sizeof(float));
			transformImage(width, height, fx, fy, depthSrc, normalsSrc, points, normals, trans.matrix());
			trans = Eigen::Transform<float, 3, Affine>(computeLinearApproxICP(width, height, normThresh, distThresh, ptThresh, points, normals, depthDst, normalsDst)) * trans;
		}
		return trans.matrix();
	}
private:
	int width, height;
	float *normals;
	float *points;
	std::unique_ptr<float[]> dataPtr;
};

//class iteratedICP_MS
//{
//public:
//	iteratedICP_MS(const int width, const int height, const std::vector<int> iters)
//		: width(width), height(height), iters(iters)
//	{
//		for (int scale = 0; scale < iters.size(); scale++) {
//			auto scaleF = 2 << scale;
//			scales.emplace_back(width/scaleF, height / scaleF);
//		}
//	}
//	Eigen::Matrix4f runICPIter(const float fx, const float fy,
//		const float normThresh, const float distThresh, const size_t ptThresh,
//		const float *depthSrc, const float *normalsSrc,
//		const float *depthDst, const float *normalsDst)
//	{
//		for (int scale = 0; scale < iters.size(); scale++) {
//			auto scaleF = 2 << scale;
//			auto iter = iters[scale];
//			scales[iter].runICPIter(iter, fx / scaleF, fy / scaleF, normThresh, distThresh, ptThresh / (scaleF*scaleF), depthSrc, normalsSrc, depthDst, normalsDst);
//		}
//	}
//private:
//	int width, height;
//	std::vector<iteratedICP> scales;
//	std::vector<int> iters;
//};
