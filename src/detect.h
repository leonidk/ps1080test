#pragma once
#include "linalg.h"
#include <vector>
#include <memory>


struct planeCandidate {
    double d;
    float n[3];
};

template <int numAngles>
class planeDetector {
public:
	planeDetector() : numAngleSym(2 * numAngles + 1), idx((2 * numAngles + 1)*(2 * numAngles + 1))
	{
	}
	std::vector<planeCandidate> detectPlanes(const int maxPlanes, const int minPts, const int width, const int height, const float *points, const float *normals)
	{
		std::vector<planeCandidate> returnData;
		for (auto &p : planes) {
			p.clear();
		}
		memset(debugImg, 0, 2 * numAngleSym * numAngleSym);
		auto numPixels = height*width;
		for (int i = 0; i < numPixels; i++) {
			if (normals[3 * i + 2] != 0) {
				int idx1 = numAngles + (int)(normals[3 * i + 0] * numAngles);
				int idx2 = numAngles + (int)(normals[3 * i + 1] * numAngles);
				double d = -normals[3 * i] * points[3*i] - normals[3 * i + 1] * points[3*i+1] - normals[3 * i + +2] * points[3*i+2];
				planes[numAngleSym * idx2 + idx1].push_back({ d, { normals[3*i], normals[3*i+1], normals[3*i+2] } });
			}
		}
		
		for (int i = 0; i < numAngleSym * numAngleSym; i++)
			idx[i] = i;
		auto xySize = [&](int x, int y) { return planes[y*numAngleSym + x].size(); };
		for (int i = 1; i < numAngleSym - 1; i++) {
			for (int j = 1; j < numAngleSym - 1; j++) {
				auto c = xySize(j, i);
				//if (c < minPts) idx[i*numAngleSym + j]=0;

				if (c < xySize(j - 1, i - 1))
					idx[i * numAngleSym + j] = 0;
				if (c < xySize(j + 0, i - 1))
					idx[i * numAngleSym + j] = 0;
				if (c < xySize(j + 1, i - 1))
					idx[i * numAngleSym + j] = 0;
				if (c < xySize(j - 1, i + 0))
					idx[i * numAngleSym + j] = 0;
				if (c < xySize(j + 1, i + 0))
					idx[i * numAngleSym + j] = 0;
				if (c < xySize(j - 1, i + 1))
					idx[i * numAngleSym + j] = 0;
				if (c < xySize(j + 0, i + 1))
					idx[i * numAngleSym + j] = 0;
				if (c < xySize(j + 1, i + 1))
					idx[i * numAngleSym + j] = 0;
			}
		}

		std::partial_sort(std::begin(idx), std::begin(idx) + maxPlanes, std::end(idx), [&](const int i1, const int i2) {
			return planes[i1].size() > planes[i2].size();
		});
		int realCand = 0;
		for (; realCand < maxPlanes; realCand++) {
			if (planes[idx[realCand]].size() < minPts)
				break;
		}

		for (auto i = idx.begin(); i < (idx.begin() + maxPlanes); i++) {
			std::vector<planeCandidate> candList;
			int currY = (*i / numAngleSym);
			int currX = (*i % numAngleSym);

			candList.insert(candList.end(), planes[currY * numAngleSym + currX].begin(), planes[currY * numAngleSym + currX].end());
			if (currY - 1 >= 0)
				candList.insert(candList.end(), planes[(currY - 1) * numAngleSym + currX].begin(), planes[(currY - 1) * numAngleSym + currX].end());
			if (currY + 1 < numAngleSym)
				candList.insert(candList.end(), planes[(currY + 1) * numAngleSym + currX].begin(), planes[(currY + 1) * numAngleSym + currX].end());
			if (currX - 1 >= 0)
				candList.insert(candList.end(), planes[(currY)* numAngleSym + currX - 1].begin(), planes[(currY)* numAngleSym + currX - 1].end());
			if (currX + 1 < numAngleSym)
				candList.insert(candList.end(), planes[(currY)* numAngleSym + currX + 1].begin(), planes[(currY)* numAngleSym + currX + 1].end());

			//Idea1: lets GMM this?

			//Idea2: sort lists by 'd'
			//Idea2: linear search, runnnig average, when startIdx < Mean -2sigma, stop. Repeat
			//Idea2: Use threshold instead of computing sigma

			//Idea3: be lazy
			if (candList.size() < minPts)
				continue;
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
			debugImg[*i] = planes[*i].size();
		}

		return returnData;
	}
	uint16_t* getDebugImg() 
	{
		return debugImg;
	}
private:
	size_t numAngleSym;
	uint16_t debugImg[(2 * numAngles + 1)*(2 * numAngles + 1)];
	std::vector<planeCandidate> planes[(2 * numAngles + 1)*(2 * numAngles + 1)];
	std::vector<int> idx;
};


#define RAD_SIZE (15)
#define RAD_FULL (2 * RAD_SIZE + 1)

#define MAX_CAND (10)
#define MIN_PTS (75)

//old, original, debugging function. 
//template <int size>
//std::vector<planeCandidate> generateNormalsAndPlanes(const float *points, const int width, const int height, float *normals, uint16_t *image) {
//    std::vector<planeCandidate> returnData;
//    static std::vector<planeCandidate> planes[RAD_FULL * RAD_FULL];
//    for (auto &p : planes) {
//        p.clear();
//    }
//
//    memset(image, 0, 2 * RAD_FULL * RAD_FULL);
//    for (int i = size; i < height - size; i++) {
//        for (int j = size; j < width - size; j++) {
//            if (!points[3 * (i * width + j) + 2])
//                continue;
//            const float *pc = &points[3 * (i * width + j)];
//
//            float outNorm[3] = { 0, 0, 0 };
//            int count = 0;
//            if (points[3 * (i * width + j + size) + 2] && points[3 * ((i + size) * width + j) + 2]) {
//                const float *px = &points[3 * (i * width + j + size)];
//                const float *py = &points[3 * ((i + size) * width + j)];
//
//                float v1[] = { px[0] - pc[0], px[1] - pc[1], px[2] - pc[2] };
//                float v2[] = { py[0] - pc[0], py[1] - pc[1], py[2] - pc[2] };
//
//                float v3[3];
//                cross(v1, v2, v3);
//                normalize(v3);
//                outNorm[0] += v3[0];
//                outNorm[1] += v3[1];
//                outNorm[2] += v3[2];
//                count++;
//            }
//            if (points[3 * (i * width + j - size) + 2] && points[3 * ((i + size) * width + j) + 2]) {
//                const float *px = &points[3 * (i * width + j - size)];
//                const float *py = &points[3 * ((i + size) * width + j)];
//
//                float v1[] = { pc[0] - px[0], pc[1] - px[1], pc[2] - px[2] };
//                float v2[] = { py[0] - pc[0], py[1] - pc[1], py[2] - pc[2] };
//
//                float v3[3];
//                cross(v1, v2, v3);
//                normalize(v3);
//                outNorm[0] += v3[0];
//                outNorm[1] += v3[1];
//                outNorm[2] += v3[2];
//                count++;
//            }
//            if (points[3 * (i * width + j + size) + 2] && points[3 * ((i - size) * width + j) + 2]) {
//                const float *px = &points[3 * (i * width + j + size)];
//                const float *py = &points[3 * ((i - size) * width + j)];
//
//                float v1[] = { px[0] - pc[0], px[1] - pc[1], px[2] - pc[2] };
//                float v2[] = { pc[0] - py[0], pc[1] - py[1], pc[2] - py[2] };
//
//                float v3[3];
//                cross(v1, v2, v3);
//                normalize(v3);
//                outNorm[0] += v3[0];
//                outNorm[1] += v3[1];
//                outNorm[2] += v3[2];
//                count++;
//            }
//            if (points[3 * (i * width + j - size) + 2] && points[3 * ((i - size) * width + j) + 2]) {
//                const float *px = &points[3 * (i * width + j - size)];
//                const float *py = &points[3 * ((i - size) * width + j)];
//
//                float v1[] = { pc[0] - px[0], pc[1] - px[1], pc[2] - px[2] };
//                float v2[] = { pc[0] - py[0], pc[1] - py[1], pc[2] - py[2] };
//
//                float v3[3];
//                cross(v1, v2, v3);
//                normalize(v3);
//                outNorm[0] += v3[0];
//                outNorm[1] += v3[1];
//                outNorm[2] += v3[2];
//                count++;
//            }
//            if (count) {
//                float v3[3] = { outNorm[0] / count, outNorm[1] / count, outNorm[2] / count };
//                normalize(v3);
//
//                normals[3 * (i * width + j)] = v3[0];
//                normals[3 * (i * width + j) + 1] = v3[1];
//                normals[3 * (i * width + j) + 2] = v3[2];
//
//                double d = -v3[0] * pc[0] - v3[1] * pc[1] - v3[2] * pc[2];
//                int idx1 = RAD_SIZE + (int)(v3[0] * RAD_SIZE);
//                int idx2 = RAD_SIZE + (int)(v3[1] * RAD_SIZE);
//                assert(idx1 >= 0);
//                assert(idx2 >= 0);
//                assert(idx1 < RAD_FULL);
//                assert(idx2 < RAD_FULL);
//
//                planes[RAD_FULL * idx2 + idx1].push_back({ d, { v3[0], v3[1], v3[2] } });
//            }
//        }
//    }
//    std::vector<int> idx(RAD_FULL * RAD_FULL);
//    for (int i = 0; i < RAD_FULL * RAD_FULL; i++)
//        idx[i] = i;
//    auto xySize = [&](int x, int y) { return planes[y*RAD_FULL + x].size(); };
//    for (int i = 1; i < RAD_FULL - 1; i++) {
//        for (int j = 1; j < RAD_FULL - 1; j++) {
//            auto c = xySize(j, i);
//            //if (c < MIN_PTS) idx[i*RAD_FULL + j]=0;
//
//            if (c < xySize(j - 1, i - 1))
//                idx[i * RAD_FULL + j] = 0;
//            if (c < xySize(j + 0, i - 1))
//                idx[i * RAD_FULL + j] = 0;
//            if (c < xySize(j + 1, i - 1))
//                idx[i * RAD_FULL + j] = 0;
//            if (c < xySize(j - 1, i + 0))
//                idx[i * RAD_FULL + j] = 0;
//            if (c < xySize(j + 1, i + 0))
//                idx[i * RAD_FULL + j] = 0;
//            if (c < xySize(j - 1, i + 1))
//                idx[i * RAD_FULL + j] = 0;
//            if (c < xySize(j + 0, i + 1))
//                idx[i * RAD_FULL + j] = 0;
//            if (c < xySize(j + 1, i + 1))
//                idx[i * RAD_FULL + j] = 0;
//        }
//    }
//
//    std::partial_sort(std::begin(idx), std::begin(idx) + MAX_CAND, std::end(idx), [&](const int i1, const int i2) {
//		return planes[i1].size() > planes[i2].size();
//    });
//    int realCand = 0;
//    for (; realCand < MAX_CAND; realCand++) {
//        if (planes[idx[realCand]].size() < MIN_PTS)
//            break;
//    }
//
//    for (auto i = idx.begin(); i < (idx.begin() + MAX_CAND); i++) {
//        std::vector<planeCandidate> candList;
//        int currY = (*i / RAD_FULL);
//        int currX = (*i % RAD_FULL);
//
//        candList.insert(candList.end(), planes[currY * RAD_FULL + currX].begin(), planes[currY * RAD_FULL + currX].end());
//        if (currY - 1 >= 0)
//            candList.insert(candList.end(), planes[(currY - 1) * RAD_FULL + currX].begin(), planes[(currY - 1) * RAD_FULL + currX].end());
//        if (currY + 1 < RAD_FULL)
//            candList.insert(candList.end(), planes[(currY + 1) * RAD_FULL + currX].begin(), planes[(currY + 1) * RAD_FULL + currX].end());
//        if (currX - 1 >= 0)
//            candList.insert(candList.end(), planes[(currY) * RAD_FULL + currX - 1].begin(), planes[(currY) * RAD_FULL + currX - 1].end());
//        if (currX + 1 < RAD_FULL)
//            candList.insert(candList.end(), planes[(currY) * RAD_FULL + currX + 1].begin(), planes[(currY) * RAD_FULL + currX + 1].end());
//
//        //Idea1: lets GMM this?
//
//        //Idea2: sort lists by 'd'
//        //Idea2: linear search, runnnig average, when startIdx < Mean -2sigma, stop. Repeat
//        //Idea2: Use threshold instead of computing sigma
//
//        //Idea3: be lazy
//        if (candList.size() < MIN_PTS)
//            continue;
//        planeCandidate avgCand = { 0, { 0, 0, 0 } };
//        for (const auto c : candList) {
//            avgCand.d += c.d;
//            avgCand.n[0] += c.n[0];
//            avgCand.n[1] += c.n[1];
//            avgCand.n[2] += c.n[2];
//        }
//        avgCand.d /= candList.size();
//        avgCand.n[0] /= candList.size();
//        avgCand.n[1] /= candList.size();
//        avgCand.n[2] /= candList.size();
//        normalize(avgCand.n);
//        returnData.push_back(std::move(avgCand));
//        image[*i] = planes[*i].size();
//    }
//
//    return returnData;
//}
//
//std::vector<std::pair<int, int> > generateCorrespondencesSVD(const std::vector<planeCandidate> &src, const std::vector<planeCandidate> &dst) {
//    using namespace Eigen;
//    using namespace std;
//    if (src.size() * dst.size() == 0)
//        return std::move(std::vector<std::pair<int, int> >({}));
//    MatrixXf m(src.size(), dst.size());
//    auto sqr = [](const float a) { return a*a; };
//    const auto planeDiff = [sqr](const planeCandidate &a, const planeCandidate &b) {
//		return exp(-sqr(a.n[0] - b.n[0]) / 0.25 - sqr(a.n[1] - b.n[1]) / 0.25 - sqr(a.d - b.d) / sqr(50));
//    };
//    const auto planeDiffOrig = [sqr](const planeCandidate &a, const planeCandidate &b) {
//		return abs(a.n[0] - b.n[0]) + abs(a.n[1] - b.n[1]) + abs(a.d - b.d) / 50;
//    };
//    for (int i = 0; i < src.size(); i++) {
//        for (int j = 0; j < dst.size(); j++) {
//            m(i, j) = planeDiff(src[i], dst[j]);
//        }
//    }
//    //cout << "Here is the matrix m:" << endl << m << endl;
//    JacobiSVD<MatrixXf> svd(m, ComputeThinU | ComputeThinV);
//    MatrixXf u = svd.matrixU();
//    MatrixXf s = svd.singularValues().asDiagonal();
//    MatrixXf v = svd.matrixV();
//
//    MatrixXf res = u * s * v.adjoint();
//
//    res = u * MatrixXf::Identity(s.rows(), s.cols()) * v.adjoint();
//    //cout << "It's P matrix is" << endl << res << endl;
//
//    VectorXf colMax = res.colwise().maxCoeff();
//    VectorXf rowMax = res.rowwise().maxCoeff();
//    //cout << colMax << endl << endl <<  rowMax << endl;
//    std::vector<std::pair<int, int> > corrPairs;
//
//    for (int i = 0; i < res.rows(); i++) {
//        for (int j = 0; j < res.cols(); j++) {
//            auto v = res(i, j);
//            if (v == colMax(j) && v == rowMax(i))
//                corrPairs.push_back(std::pair<int, int>(i, j));
//        }
//    }
//
//    return std::move(corrPairs);
//}
//
//std::vector<std::pair<int, int> > generateCorrespondencesMP(const std::vector<planeCandidate> &src, const std::vector<planeCandidate> &dst) {
//    if (src.size() * dst.size() == 0)
//        return std::move(std::vector<std::pair<int, int> >({}));
//    auto sqr = [](const float a) { return a*a; };
//    const auto planeDiff = [sqr](const planeCandidate &a, const planeCandidate &b) {
//		return exp(-sqr(a.n[0] - b.n[0]) / 0.25 - sqr(a.n[1] - b.n[1]) / 0.25 - sqr(a.d - b.d) / sqr(50));
//    };
//    std::vector<std::pair<int, int> > corrPairs;
//
//    std::vector<int> occMask(dst.size(), 0);
//
//    for (int i = 0; i < src.size(); i++) {
//        int maxIdx = -1;
//        double maxVal = -1;
//        for (int j = 0; j < dst.size(); j++) {
//            if (occMask[j] == 0) {
//                const auto p = planeDiff(src[i], dst[j]);
//                if (p > maxVal) {
//                    maxVal = p;
//                    maxIdx = j;
//                }
//            }
//        }
//        if (maxIdx >= 0) {
//            corrPairs.push_back({ i, maxIdx });
//            occMask[maxIdx] = 1;
//        }
//    }
//
//    return std::move(corrPairs);
//}
