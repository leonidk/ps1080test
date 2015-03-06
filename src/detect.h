#pragma once
#include <vector>
#include <memory>
#include <unordered_map>
#include <map>
#include <unordered_set>
#include <cmath>
#include "linalg.h"
#include "nanoflann.hpp"

struct planeCandidate {
	double d;
	float n[3];
	float stddev;
};

struct PointCloud
{
	struct Point
	{
		float x, y, d;
	};
	std::vector<Point> pts;
	inline size_t kdtree_get_point_count() const { return pts.size(); }
	inline float kdtree_distance(const float *p1, const size_t idx_p2, size_t /*size*/) const
	{
		const float d0 = p1[0] - pts[idx_p2].x;
		const float d1 = p1[1] - pts[idx_p2].y;
		const float d2 = p1[2] - pts[idx_p2].d;
		return d0*d0 + d1*d1 + d2*d2;
	}
	inline float kdtree_get_pt(const size_t idx, int dim) const
	{
		if (dim == 0) return pts[idx].x;
		else if (dim == 1) return pts[idx].y;
		else return pts[idx].d;
	}
	template <class BBOX>
	bool kdtree_get_bbox(BBOX& /*bb*/) const { return false; }
};

class planeDetectorMS {
public:
	planeDetectorMS() {}
	std::vector<planeCandidate> detectPlanes(const size_t minPts, const float radius, const float distToNormal, const int width, const int height, const float *points, const float *normals)
	{
		std::vector<planeCandidate> returnData;
		PointCloud planeCloud; 
		std::vector<std::pair<PointCloud::Point,size_t>> centers;
		const auto tolerance = square(radius*0.0625);
		auto numPixels = height*width;
		for (int i = 0; i < numPixels; i++) {
			if (std::isnormal(normals[3 * i + 2])) {
				float d = -normals[3 * i] * points[3 * i] - normals[3 * i + 1] * points[3 * i + 1] - normals[3 * i + +2] * points[3 * i + 2];
				planeCloud.pts.push_back({ normals[3 * i], normals[3 * i + 1], distToNormal*d });
			}
		} 
		typedef nanoflann::KDTreeSingleIndexAdaptor<
			nanoflann::L2_Simple_Adaptor<float, PointCloud >,
			PointCloud, 3 > my_kd_tree_t;
		my_kd_tree_t index(3, planeCloud, nanoflann::KDTreeSingleIndexAdaptorParams(20));
		index.buildIndex();

		for (const auto p : planeCloud.pts) {
			const float search_radius = static_cast<float>(radius);
			nanoflann::SearchParams params;
			params.sorted = false;
			float diff = 0;
			PointCloud::Point pointCenter = p;
			int iters = 0;
			size_t nPts = 0;
			do {
				iters++;
				const float query_pt[3] = { pointCenter.x, pointCenter.y, pointCenter.d };
				std::vector<std::pair<size_t, float> > ret_matches;
				nPts = index.radiusSearch(&query_pt[0], search_radius, ret_matches, params);
				for (const auto matchIdx : ret_matches) {
					const auto match = planeCloud.pts[matchIdx.first];
					pointCenter.x += match.x;
					pointCenter.y += match.y;
					pointCenter.d += match.d;
				}
				//printf("%d\n", ret_matches.size());
				pointCenter.x /= (float)ret_matches.size();
				pointCenter.y /= (float)ret_matches.size();
				pointCenter.d /= (float)ret_matches.size();
				diff = square(pointCenter.x - query_pt[0]) + square(pointCenter.y - query_pt[1]) + square(pointCenter.d - query_pt[2]);
			} while (diff > tolerance);
			//printf("%d\n", iters);
				auto candP = std::find_if(centers.begin(), centers.end(), [&]
					(std::pair<PointCloud::Point, size_t> a){
					return (square(a.first.x - pointCenter.x) + square(a.first.y - pointCenter.y) + square(a.first.d - pointCenter.d) < radius*radius)
						? true : false;
				});
				if (candP != centers.end()) {
					auto weightOld = candP->second / (candP->second + 1.0);
					auto weightNew = 1.0 / (candP->second + 1.0);
					candP->first.x = (weightOld*candP->first.x + weightNew*pointCenter.x);
					candP->first.y = (weightOld*candP->first.y + weightNew*pointCenter.y);
					candP->first.d = (weightOld*candP->first.d + weightNew*pointCenter.d);
					candP->second += 1;
				}
				else {
					centers.push_back({ pointCenter, 1 });
				}
			//centers[pointCenter] = centers[pointCenter] + 1;
		}
		for (const auto kv : centers) {
			if (kv.second > minPts) {
				planeCandidate pc = { 0 };
				pc.d = kv.first.d / distToNormal;
				pc.n[0] = kv.first.x;
				pc.n[1] = kv.first.y;
				pc.n[2] = sqrtf(1.0f - square(pc.n[0]) - square(pc.n[1]));
				pc.stddev = 25;
				returnData.push_back(pc);
			}
		}
		
		return returnData;
	}

private:
};

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
class planeDetector {
public:
	planeDetector() : numAngleSym(2 * numAngles + 1)
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
			debugImg[node.first] = planes[node.first].size();

		}

		int cntP = 0;
		//std::cout << std::endl;
		for (auto& root : planeRoots) {
			//make a local for in-place editing
			debugImg[root.first] = root.second.size();
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
				auto lim = max - min;
				int nPrev = rootCopy.size();
				n = rootCopy.size() + 1;
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

				if (n > minPts)
					returnData.push_back(avgCand);
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
