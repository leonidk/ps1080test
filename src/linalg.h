#pragma once
#include <stdio.h>
#include <type_traits>
#include <memory>
struct planeCandidate {
	double d;
	float n[3];
	float stddev;
};

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

template <typename T>
inline T square(const T a) {
    return a * a;
}

inline float sqrNorm(const float *a, const float *b) {
    return square(a[0] - b[0]) + square(a[1] - b[1]) + square(a[2] - b[2]);
}
// cX*((5*z +j*(z-z2));
inline void cross(const float *a, const float *b, float *c) {
    c[0] = (a[1] * b[2]) - (a[2] * b[1]);
    c[1] = (a[2] * b[0]) - (a[0] * b[2]);
    c[2] = (a[0] * b[1]) - (a[1] * b[0]);
}
inline float dot(const float *a, const float *b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
inline void normalize(float *a) {
    const auto norm = sqrt(dot(a, a));
    a[0] /= norm;
    a[1] /= norm;
    a[2] /= norm;
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