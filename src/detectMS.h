#pragma once
#include <vector>
#include "plane.h"
#include "nanoflann.hpp"

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
		std::vector<std::pair<PointCloud::Point, size_t>> centers;
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
				candP->first.x = (float)(weightOld*candP->first.x + weightNew*pointCenter.x);
				candP->first.y = (float)(weightOld*candP->first.y + weightNew*pointCenter.y);
				candP->first.d = (float)(weightOld*candP->first.d + weightNew*pointCenter.d);
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
				pc.cnt = kv.second;
				returnData.push_back(pc);
			}
		}

		return returnData;
	}

private:
	template <typename T>
	inline T square(const T a) {
		return a * a;
	}

};