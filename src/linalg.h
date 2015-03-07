#pragma once
#include <stdio.h>

struct planeCandidate {
	double d;
	float n[3];
	float stddev;
};

template <typename T, int size>
void generateHalfImage(const T *in, T *out, const int outW, const int outH) {
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

template <typename T>
T square(const T a) {
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
