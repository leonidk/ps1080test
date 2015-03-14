#pragma once
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include "plane.h"

class DisjointSets {
public:
	DisjointSets() {}
	~DisjointSets() {}
	void makeParent(int child, int parent) {
		if (map.find(parent) == map.end() ) {
			map[parent] = parent;
		}
		map[child] = find(parent);
	}
	int find(int node){
		if (map[node] != node) {
			map[node] = find(map[node]);
		}
		return map[node];
	}
	void Union(int node1, int node2){
		int root1 = find(node1);
		int root2 = find(node2);
		map[root1] = root2;
	}
	const std::unordered_map<int, int> & getMap() { return map; }
private:
	std::unordered_map<int, int> map;
};

template <int numAngles>
class planeDetectorDisjointModes {
public:
	planeDetectorDisjointModes() : numAngleSym(2 * numAngles + 1)
	{
	}
	std::vector<planeCandidate> detectPlanes(const int minPts, const float stdDev, const int width, const int height, const float *points, const float *normals)
	{
		std::vector<planeCandidate> returnData;
		DisjointSets ds;
		auto xySize = [&](int x, int y) { return (int)planes[y*numAngleSym + x].size(); };
		auto planeCompare = [](const planeCandidate & a, const planeCandidate & b)	{	return a.d < b.d;	};

		for (auto &p : planes) {
			p.clear();
		}
		memset(debugImg, 0, 2 * numAngleSym * numAngleSym);
		auto numPixels = height*width;
		for (int i = 0; i < numPixels; i++) {
			if (normals[3 * i + 2] != 0) {
				int idx1 = numAngles + (int)(normals[3 * i + 0] * numAngles);
				int idx2 = numAngles + (int)(normals[3 * i + 1] * numAngles);
				double d = -normals[3 * i] * points[3 * i] - normals[3 * i + 1] * points[3 * i + 1] - normals[3 * i + +2] * points[3 * i + 2];
				planes[numAngleSym * idx2 + idx1].push_back({ d, { normals[3 * i], normals[3 * i + 1], normals[3 * i + 2] } });
			}
		}

		for (int i = 0; i < numAngleSym * numAngleSym; i++) {
			if (planes[i].size())
				ds.makeParent(i, i);
		}

		auto getSizeIdxPair = [&](int x, int y) { return std::pair<int, int>({ y*numAngleSym + x, xySize(x, y) }); };
		for (int i = 1; i < numAngleSym - 1; i++) {
			for (int j = 1; j < numAngleSym - 1; j++) {
				auto c = xySize(j, i);
				auto idx = i*numAngleSym + j;

				if (c == 0) continue;
				const std::array < std::pair<int, int>, 8> neighbors = {
					getSizeIdxPair(j - 1, i - 1),
					getSizeIdxPair(j + 0, i - 1),
					getSizeIdxPair(j + 1, i - 1),
					getSizeIdxPair(j - 1, i + 0),
					getSizeIdxPair(j + 1, i + 0),
					getSizeIdxPair(j - 1, i + 1),
					getSizeIdxPair(j + 0, i + 1),
					getSizeIdxPair(j + 1, i + 1) };
				auto maxE = std::max_element(neighbors.begin(), neighbors.end(), [](const std::pair<int, int> a, const std::pair<int, int> b) {
					return a.second < b.second;
				});
				if (maxE->second > c)
					ds.makeParent(idx, maxE->first);
			}
		}
		std::unordered_map<int, std::vector<planeCandidate>> planeRoots;

		auto map = ds.getMap();
		//insert into root list
		for (const auto &node : ds.getMap()) {
			auto& list = planeRoots[ds.find(node.first)];
			list.insert(list.end(), planes[node.first].begin(), planes[node.first].end());
			debugImg[node.first] = (uint16_t)(planes[node.first].size());

		}

		for (auto& root : planeRoots) {
			debugImg[root.first] = (uint16_t)(root.second.size());

			auto bounds = std::minmax_element(root.second.begin(), root.second.end(), planeCompare);
			auto min = bounds.first->d;
			auto max = bounds.second->d;
			auto maxMin = max - min;
			std::vector<planeCandidate> buckets(std::ceil(maxMin/stdDev + 0.001));
			for (const auto & p : root.second) {
				auto diff = (p.d - min);
				auto idx = (int)(diff / (stdDev*1.00001)); //to avoid single element last bucket
				auto &  bucket = buckets[idx];
				bucket.n[0] += p.n[0];
				bucket.n[1] += p.n[1];
				bucket.n[2] += p.n[2];
				bucket.d += p.d;
				bucket.cnt++;
			}
			auto bNum = buckets.size();
			if (bNum == 1) {
				auto b = buckets[0];
				if (b.cnt >= minPts) {
					b.n[0] /= b.cnt;
					b.n[1] /= b.cnt;
					b.n[2] /= b.cnt;
					b.d /= b.cnt;
					normalize(b.n);
					b.stddev = stdDev;
					returnData.push_back(b);
				}
			}
			else if (bNum == 2){
				auto b1 = buckets[0];
				auto b2 = buckets[1];
				auto tPts = b1.cnt + b2.cnt;
				if (tPts >= minPts) {
					planeCandidate avgCand = { 0, { 0, 0, 0 } };
					avgCand.n[0] = (b1.n[0] + b2.n[0]) / tPts;
					avgCand.n[1] = (b1.n[1] + b2.n[1]) / tPts;
					avgCand.n[2] = (b1.n[2] + b2.n[2]) / tPts;
					avgCand.d = (b1.d + b2.d) / tPts;
					normalize(avgCand.n);
					avgCand.stddev = stdDev;
					returnData.push_back(avgCand);
				}
			}
			else if (bNum >= 3) {
				for (int i = 1; i < bNum - 1; i++) {
					auto & b1 = buckets[i - 1];
					auto & b2 = buckets[i];
					auto & b3 = buckets[i + 1];
					auto tPts = b1.cnt + b2.cnt + b3.cnt;
					if (tPts >= minPts) {
						b1.cnt = 0;
						b2.cnt = 0;
						b3.cnt = 0;
						planeCandidate avgCand = { 0, { 0, 0, 0 } };
						avgCand.n[0] = (b1.n[0] + b2.n[0] + b3.n[0]) / tPts;
						avgCand.n[1] = (b1.n[1] + b2.n[1] + b3.n[1]) / tPts;
						avgCand.n[2] = (b1.n[2] + b2.n[2] + b3.n[2]) / tPts;
						avgCand.d = (b1.d + b2.d + b3.d) / tPts;
						normalize(avgCand.n);
						avgCand.stddev = stdDev;
						returnData.push_back(avgCand);
					}
				}
			}
			
		}

		return returnData;
	}
	uint16_t* getDebugImg()
	{
		return debugImg;
	}
private:
	int numAngleSym;
	uint16_t debugImg[(2 * numAngles + 1)*(2 * numAngles + 1)];
	std::vector<planeCandidate> planes[(2 * numAngles + 1)*(2 * numAngles + 1)];
};

template <int numAngles>
class planeDetectorDisjointGMM {
public:
	planeDetectorDisjointGMM() : numAngleSym(2 * numAngles + 1)
	{
	}
	std::vector<planeCandidate> detectPlanes(const int minPts, const float stdDev, const int width, const int height, const float *points, const float *normals)
	{
		std::vector<planeCandidate> returnData;
		DisjointSets ds;
		auto xySize = [&](int x, int y) { return (int)planes[y*numAngleSym + x].size(); };
		auto planeCompare = [](const planeCandidate & a, const planeCandidate & b)	{	return a.d < b.d;	}	;

		for (auto &p : planes) {
			p.clear();
		}
		memset(debugImg, 0, 2 * numAngleSym * numAngleSym);
		auto numPixels = height*width;
		for (int i = 0; i < numPixels; i++) {
			if (normals[3 * i + 2] != 0) {
				int idx1 = numAngles + (int)(normals[3 * i + 0] * numAngles);
				int idx2 = numAngles + (int)(normals[3 * i + 1] * numAngles);
				double d = -normals[3 * i] * points[3 * i] - normals[3 * i + 1] * points[3 * i + 1] - normals[3 * i + +2] * points[3 * i + 2];
				planes[numAngleSym * idx2 + idx1].push_back({ d, { normals[3 * i], normals[3 * i + 1], normals[3 * i + 2] } });
			}
		}

		for (int i = 0; i < numAngleSym * numAngleSym; i++) {
			if(planes[i].size()) 
				ds.makeParent(i, i);
		}

		auto getSizeIdxPair = [&](int x, int y) { return std::pair<int, int>({ y*numAngleSym + x, xySize(x, y) }); };
		for (int i = 1; i < numAngleSym - 1; i++) {
			for (int j = 1; j < numAngleSym - 1; j++) {
				auto c = xySize(j, i);
				auto idx = i*numAngleSym + j;

				if (c == 0) continue;
				const std::array < std::pair<int, int>, 8> neighbors = {
					getSizeIdxPair(j - 1, i - 1),
					getSizeIdxPair(j + 0, i - 1),
					getSizeIdxPair(j + 1, i - 1),
					getSizeIdxPair(j - 1, i + 0),
					getSizeIdxPair(j + 1, i + 0),
					getSizeIdxPair(j - 1, i + 1),
					getSizeIdxPair(j + 0, i + 1),
					getSizeIdxPair(j + 1, i + 1) };
				auto maxE = std::max_element(neighbors.begin(), neighbors.end(), [](const std::pair<int, int> a, const std::pair<int, int> b) {
					return a.second < b.second; 
				});
				if (maxE->second > c)
					ds.makeParent(idx, maxE->first);
			}
		}
		std::unordered_map<int,std::vector<planeCandidate>> planeRoots;
		//USE A DISJOINT SET
		//COULD ALSO USE A POSITIVE GRAD SEGMENTATION
		//This is an approximation that could also be solved with a 3D meanshift (nx/ny/d). Hopefully this is faster?
		auto map = ds.getMap();
		//insert into root list
		for (const auto &node : ds.getMap()) {
			auto& list = planeRoots[ds.find(node.first)];
			list.insert(list.end(), planes[node.first].begin(), planes[node.first].end());
			debugImg[node.first] = (uint16_t)(planes[node.first].size());

		}

		int cntP = 0;
		//std::cout << std::endl;
		for (auto& root : planeRoots) {
			//make a local for in-place editing
			debugImg[root.first] = (uint16_t)(root.second.size());
			std::vector<planeCandidate> rootCopy = std::move(root.second);
			cntP++;
			
			while ((int)rootCopy.size() > minPts){
				planeCandidate avgCand = { 0, { 0, 0, 0 } };
				int n = 0;
				double mean = 0, M2 = 0;
				auto bounds = std::minmax_element(rootCopy.begin(), rootCopy.end(), planeCompare);
				auto min = bounds.first->d;
				auto max = bounds.second->d;
				auto center = min + (max - min) / 2;
				auto lim = (max - min);
				int nPrev = (int)(rootCopy.size());
				n = (int)(rootCopy.size() + 1);
				//loop over until we get a stable subset
				do {
					nPrev = n - 1;
					n = 0;
					mean = 0;
					M2 = 0;
					avgCand = { 0, { 0, 0, 0 } };
					//loop over candidates points, computing mean and variance
					for (const auto c : rootCopy) {
						if (fabs(c.d - center) < lim) {
							n++;
							double delta = c.d - mean;
							mean += delta / n;
							M2 += delta*(c.d - mean);
							avgCand.d += c.d;
							avgCand.n[0] += c.n[0];
							avgCand.n[1] += c.n[1];
							avgCand.n[2] += c.n[2];
						}
					}
					auto var = M2 / (n - 1);
					center = mean;
					lim = stdDev * sqrt(var);
				} while (n >= 2 && nPrev != (n - 1));
				//remove this subset from global list
				rootCopy.erase(std::remove_if(rootCopy.begin(), rootCopy.end(),
					[&](const planeCandidate& a) {return fabs(a.d - center) < lim ? true : false;  }),
					rootCopy.end());
				avgCand.d /= n - 1;
				normalize(avgCand.n);
				avgCand.stddev = (float)lim / stdDev;
				//std::cout << cntP << '\t' << n << '\t' << rootCopy.size() << '\t' << avgCand.stddev << std::endl;
				avgCand.cnt = n;
				if (n > minPts)
					returnData.push_back(avgCand);
			}
		}
		std::sort(returnData.begin(), returnData.end(), [](const  planeCandidate &a, const  planeCandidate&b) 
		{
			return a.cnt < b.cnt;
		});
		return returnData;
	}
	uint16_t* getDebugImg()
	{
		return debugImg;
	}
private:
	int numAngleSym;
	uint16_t debugImg[(2 * numAngles + 1)*(2 * numAngles + 1)];
	std::vector<planeCandidate> planes[(2 * numAngles + 1)*(2 * numAngles + 1)];
};

template <int numAngles>
class planeDetectorFast {
public:
	planeDetectorFast() : numAngleSym(2 * numAngles + 1), idx((2 * numAngles + 1)*(2 * numAngles + 1))
	{
	}
	std::vector<planeCandidate> detectPlanes(const int minPts, const float stdDev, const int width, const int height, const float *points, const float *normals)
	{
		auto maxPlanes = 10;
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
				double d = -normals[3 * i] * points[3 * i] - normals[3 * i + 1] * points[3 * i + 1] - normals[3 * i + +2] * points[3 * i + 2];
				planes[numAngleSym * idx2 + idx1].push_back({ d, { normals[3 * i], normals[3 * i + 1], normals[3 * i + 2] } });
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
			int currY = (int)(*i / numAngleSym);
			int currX = (int)(*i % numAngleSym);
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
			avgCand.stddev = 25.0f;
			avgCand.cnt = candList.size();

			returnData.push_back(std::move(avgCand));
			debugImg[*i] = (uint16_t)planes[*i].size();
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
// THNKING OF GENERAL IDEAS
//IDEA: QUICK_CLUSTER. partion-based clustering.
// IDEA: Iterate top-down and bottom up. Top-down to get variances. Bottom up for form subclusters. 
// THINKING OF PLANES:
//IDEA: Grid-based clustering in general(x,y,d)? Or at least for 2D normals. Partion-based for d values
//IDEA: What to do with 

//Uncommented for Eigen dependency
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

//find matching pairs
std::vector<std::pair<int, int> > generateCorrespondencesMP(const std::vector<planeCandidate> &src, const std::vector<planeCandidate> &dst) {
    if (src.size() * dst.size() == 0)
        return std::move(std::vector<std::pair<int, int> >({}));
    auto sqr = [](const float a) { return a*a; };
    const auto planeDiff = [sqr](const planeCandidate &a, const planeCandidate &b) {
		return exp(-sqr(a.n[0] - b.n[0]) / 0.25 - sqr(a.n[1] - b.n[1]) / 0.25 - sqr(a.d - b.d) / sqr(50));
    };
    std::vector<std::pair<int, int> > corrPairs;

    std::vector<int> occMask(dst.size(), 0);

    for (int i = 0; i < src.size(); i++) {
        int maxIdx = -1;
        double maxVal = -1;
        for (int j = 0; j < dst.size(); j++) {
            if (occMask[j] == 0) {
                const auto p = planeDiff(src[i], dst[j]);
                if (p > maxVal) {
                    maxVal = p;
                    maxIdx = j;
                }
            }
        }
        if (maxIdx >= 0) {
            corrPairs.push_back({ i, maxIdx });
            occMask[maxIdx] = 1;
        }
    }

    return std::move(corrPairs);
}