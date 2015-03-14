// Copyright (c) 2015 Sterling Orsten
// 
// This software is provided 'as-is', without any express or implied
// warranty. In no event will the author be held liable for any damages
// arising from the use of this software. You are granted a perpetual, 
// irrevocable, world-wide license to copy, modify, and redistribute
// this software for any purpose, including commercial applications.
#pragma once
#ifndef LINALG_H_INCLUDE_GUARD
#define LINALG_H_INCLUDE_GUARD

#include <cstdint>
#include <cmath>


template <typename T>
inline T square(const T a) {
	return a * a;
}

inline float sqrNorm(const float *a, const float *b) {
	return square(a[0] - b[0]) + square(a[1] - b[1]) + square(a[2] - b[2]);
}

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

namespace linalg
{
	/////////////////////
	// Data structures //
	/////////////////////

	// A vector with exactly M elements, each of which is an instance of type T. Can also be thought of as an M x 1 matrix.
	template<class T, int M> struct vec;
	template<class T> struct vec<T, 2>
	{
		T           x, y;

		vec() : x(), y() {}
		vec(T x, T y) : x(x), y(y) {}
		explicit    vec(T s) : x(s), y(s) {}

		const T &   operator [] (int i) const           { return (&x)[i]; }
		T &         operator [] (int i)                 { return (&x)[i]; }
	};
	template<class T> struct vec<T, 3>
	{
		T           x, y, z;

		vec() : x(), y(), z() {}
		vec(T x, T y, T z) : x(x), y(y), z(z) {}
		explicit    vec(T s) : x(s), y(s), z(s) {}

		const T &   operator [] (int i) const           { return (&x)[i]; }
		T &         operator [] (int i)                 { return (&x)[i]; }
	};
	template<class T> struct vec<T, 4>
	{
		T           x, y, z, w;

		vec() : x(), y(), z(), w() {}
		vec(T x, T y, T z, T w) : x(x), y(y), z(z), w(w) {}
		explicit    vec(T s) : x(s), y(s), z(s), w(s) {}

		const T &   operator [] (int i) const           { return (&x)[i]; }
		T &         operator [] (int i)                 { return (&x)[i]; }
	};

	// A matrix with exactly M rows and N columns, where each element is an instance of type T
	template<class T, int M, int N> struct mat;
	template<class T, int M> struct mat<T, M, 2>
	{
		typedef     vec<T, M> C;
		C           x, y;

		mat()                               {}
		mat(C x, C y) : x(x), y(y) {}
		explicit    mat(T s) : x(s), y(s) {}

		const C &   operator [] (int j) const           { return (&x)[j]; }
		const T &   operator () (int i, int j) const    { return (*this)[j][i]; }
		C &         operator [] (int j)                 { return (&x)[j]; }
		T &         operator () (int i, int j)          { return (*this)[j][i]; }
	};

	// Specialization for an M x 3 matrix, with exactly three columns, each of which is an M-element vector
	template<class T, int M> struct mat<T, M, 3>
	{
		typedef     vec<T, M> C;
		C           x, y, z;

		mat()                               {}
		mat(C x, C y, C z) : x(x), y(y), z(z) {}
		explicit    mat(T s) : x(s), y(s), z(s) {}

		const C &   operator [] (int j) const           { return (&x)[j]; }
		const T &   operator () (int i, int j) const    { return (*this)[j][i]; }
		C &         operator [] (int j)                 { return (&x)[j]; }
		T &         operator () (int i, int j)          { return (*this)[j][i]; }
	};

	// Specialization for an M x 4 matrix, with exactly four columns, each of which is an M-element vector
	template<class T, int M> struct mat<T, M, 4>
	{
		typedef     vec<T, M> C;
		C           x, y, z, w;

		mat()                               {}
		mat(C x, C y, C z, C w) : x(x), y(y), z(z), w(w) {}
		explicit    mat(T s) : x(s), y(s), z(s), w(s) {}

		const C &   operator [] (int j) const           { return (&x)[j]; }
		const T &   operator () (int i, int j) const    { return (*this)[j][i]; }
		C &         operator [] (int j)                 { return (&x)[j]; }
		T &         operator () (int i, int j)          { return (*this)[j][i]; }
	};

	/////////////////
	// Type traits //
	/////////////////

	// elem_t<T> is a type trait for the type of the underlying elements of a tensor (scalar, vector, or matrix)
	template<class T> struct elem_type { typedef T type; };
	template<class T, int M> struct elem_type<vec<T, M  >> { typedef T type; };
	template<class T, int M, int N> struct elem_type<mat<T, M, N>> { typedef T type; };
	template<class T> using elem_t = typename elem_type<T>::type;

	// binop_result_type<T,L,R>::type specifies the type of applying an elementwise operation to L and R, if said operation returns values of type T
	template<class T, class L, class R> struct binop_result_type;
	template<class T, class U, int M> struct binop_result_type<T, vec<U, M>, vec<U, M>>             { typedef vec<T, M> type; };     // vector x vector -> vector
	template<class T, class U, int M> struct binop_result_type<T, vec<U, M>, U>                    { typedef vec<T, M> type; };     // vector x scalar -> vector
	template<class T, class U, int M> struct binop_result_type<T, U, vec<U, M>>                    { typedef vec<T, M> type; };     // scalar x vector -> vector
	template<class T, class U, int M, int N> struct binop_result_type<T, mat<U, M, N>, mat<U, M, N>>  { typedef mat<T, M, N> type; };   // matrix x matrix -> matrix
	template<class T, class U, int M, int N> struct binop_result_type<T, mat<U, M, N>, U>           { typedef mat<T, M, N> type; };   // matrix x scalar -> matrix
	template<class T, class U, int M, int N> struct binop_result_type<T, U, mat<U, M, N>>           { typedef mat<T, M, N> type; };   // scalar x matrix -> matrix
	template<class L, class R> using comp_t = typename binop_result_type<bool, L, R>::type;
	template<class L, class R> using arith_t = typename binop_result_type<elem_t<L>, L, R>::type; // Note: We use the lefthand element type so that a x= b is well formed

	////////////////////////
	// Operator overloads //
	////////////////////////

	// Unary operators return a vector or matrix of the original type, by applying the operator to each element
	template<class T> T operator + (const T & v) { return map(v, [](elem_t<T> a) { return +a; }); }
	template<class T> T operator - (const T & v) { return map(v, [](elem_t<T> a) { return -a; }); }
	template<class T> T operator ! (const T & v) { return map(v, [](elem_t<T> a) { return !a; }); }
	template<class T> T operator ~ (const T & v) { return map(v, [](elem_t<T> a) { return ~a; }); }

	// Equality operators return a bool indicating exact equality or inequality of two values
	template<class T> bool operator == (const vec<T, 2> & l, const vec<T, 2> & r) { return l.x == r.x && l.y == r.y; }
	template<class T> bool operator == (const vec<T, 3> & l, const vec<T, 3> & r) { return l.x == r.x && l.y == r.y && l.z == r.z; }
	template<class T> bool operator == (const vec<T, 4> & l, const vec<T, 4> & r) { return l.x == r.x && l.y == r.y && l.z == r.z && l.w == r.w; }
	template<class T, int M> bool operator == (const mat<T, M, 2> & l, const mat<T, M, 2> & r) { return l.x == r.x && l.y == r.y; }
	template<class T, int M> bool operator == (const mat<T, M, 3> & l, const mat<T, M, 3> & r) { return l.x == r.x && l.y == r.y && l.z == r.z; }
	template<class T, int M> bool operator == (const mat<T, M, 4> & l, const mat<T, M, 4> & r) { return l.x == r.x && l.y == r.y && l.z == r.z && l.w == r.w; }
	template<class T> bool operator != (const T & l, const T & r) { return !(l == r); }

	// Assignment operators are defined trivially
	template<class L, class R> L & operator += (L & l, const R & r) { return l = l + r; }
	template<class L, class R> L & operator -= (L & l, const R & r) { return l = l - r; }
	template<class L, class R> L & operator *= (L & l, const R & r) { return l = l*r; }
	template<class L, class R> L & operator /= (L & l, const R & r) { return l = l / r; }
	template<class L, class R> L & operator %= (L & l, const R & r) { return l = l%r; }
	template<class L, class R> L & operator |= (L & l, const R & r) { return l = l | r; }
	template<class L, class R> L & operator &= (L & l, const R & r) { return l = l&r; }
	template<class L, class R> L & operator ^= (L & l, const R & r) { return l = l^r; }

	/////////////////////////////////
	// General aggregate functions //
	/////////////////////////////////

	// Return true if any element of vector or matrix v would evaluate to true in a boolean context
	template<class T> bool any(const T & v) { return reduce(v, [](bool a, bool b) { return a || b; }); }

	// Return true if all elements of vector or matrix v would evaluate to true in a boolean context
	template<class T> bool all(const T & v) { return reduce(v, [](bool a, bool b) { return a && b; }); }

	// Compute the sum of all elements in vector or matrix v
	template<class T> elem_t<T> sum(const T & v) { return reduce(v, [](elem_t<T> a, elem_t<T> b) { return a + b; }); }

	// Form a vector or matrix by taking the absolute value of each component of vector or matrix v
	template<class T> T abs(const T & v) { return map(v, [](elem_t<T> a) { return std::abs(a); }); }

	// Form a vector or matrix by taking the ceiling of each component of vector or matrix v
	template<class T> T ceil(const T & v) { return map(v, [](elem_t<T> a) { return std::ceil(a); }); }

	// Form a vector or matrix by taking the floor of each component of vector or matrix v
	template<class T> T floor(const T & v) { return map(v, [](elem_t<T> a) { return std::floor(a); }); }

	// Form a vector or matrix by rounding each component of vector or matrix v
	template<class T> T round(const T & v) { return map(v, [](elem_t<T> a) { return std::round(a); }); }

	// Form a vector or matrix by taking the componentwise maximum of two vectors or matrices
	//template<class L, class R> arith_t<L, R> max(const L & l, const R & r) { return zip<arith_t<L, R>>(l, r, [](elem_t<L> a, elem_t<R> b) { return a > b ? a : b; }); }

	// Form a vector or matrix by taking the componentwise minimum of two vectors or matrices
	//template<class L, class R> arith_t<L, R> min(const L & l, const R & r) { return zip<arith_t<L, R>>(l, r, [](elem_t<L> a, elem_t<R> b) { return a < b ? a : b; }); }

	// Return a when t = 0, b when t = 1, and linearly interpolate between a and b when t is in the interval (0,1)
	template<class T> T lerp(const T & a, const T & b, elem_t<T> t) { return a*(1 - t) + b*t; }

	//////////////////////////////
	// Vector algebra functions //
	//////////////////////////////

	// Compute the dot product of two vectors
	template<class T, int N> T dot(const vec<T, N> & l, const vec<T, N> & r) { return sum(l*r); }

	// Compute the cross product of two vectors
	template<class T> vec<T, 3> cross(const vec<T, 3> & l, const vec<T, 3> & r) { return{ l.y*r.z - l.z*r.y, l.z*r.x - l.x*r.z, l.x*r.y - l.y*r.x }; }

	// Compute the magnitude of a vector
	template<class T, int N> T mag(const vec<T, N> & v) { return std::sqrt(mag2(v)); }

	// Compute the squared magnitude of a vector
	template<class T, int N> T mag2(const vec<T, N> & v) { return dot(v, v); }

	// Compute the distance between two points, expressed as vectors
	template<class T, int N> T dist(const vec<T, N> & l, const vec<T, N> & r) { return mag(r - l); }

	// Compute the squared distance between two points, expressed as vectors
	template<class T, int N> T dist2(const vec<T, N> & l, const vec<T, N> & r) { return mag2(r - l); }

	// Compute a normalized vector with the same direction as an original vector
	template<class T, int N> vec<T, N> norm(const vec<T, N> & v) { return v / mag(v); }
	
	//////////////////////////////
	// Matrix algebra functions //
	//////////////////////////////

	// Compute the adjugate of a square matrix (equivalent to the transpose of the cofactor matrix)
	template<class T> mat<T, 2, 2> adj(const mat<T, 2, 2> & a) { return{ { a.y.y, -a.x.y }, { -a.y.x, a.x.x } }; }
	template<class T> mat<T, 3, 3> adj(const mat<T, 3, 3> & a)
	{
		return{ { a.y.y*a.z.z - a.z.y*a.y.z, a.z.y*a.x.z - a.x.y*a.z.z, a.x.y*a.y.z - a.y.y*a.x.z },
		{ a.y.z*a.z.x - a.z.z*a.y.x, a.z.z*a.x.x - a.x.z*a.z.x, a.x.z*a.y.x - a.y.z*a.x.x },
		{ a.y.x*a.z.y - a.z.x*a.y.y, a.z.x*a.x.y - a.x.x*a.z.y, a.x.x*a.y.y - a.y.x*a.x.y } };
	}
	template<class T> mat<T, 4, 4> adj(const mat<T, 4, 4> & a)
	{
		return{ { a.y.y*a.z.z*a.w.w + a.w.y*a.y.z*a.z.w + a.z.y*a.w.z*a.y.w - a.y.y*a.w.z*a.z.w - a.z.y*a.y.z*a.w.w - a.w.y*a.z.z*a.y.w,
			a.x.y*a.w.z*a.z.w + a.z.y*a.x.z*a.w.w + a.w.y*a.z.z*a.x.w - a.w.y*a.x.z*a.z.w - a.z.y*a.w.z*a.x.w - a.x.y*a.z.z*a.w.w,
			a.x.y*a.y.z*a.w.w + a.w.y*a.x.z*a.y.w + a.y.y*a.w.z*a.x.w - a.x.y*a.w.z*a.y.w - a.y.y*a.x.z*a.w.w - a.w.y*a.y.z*a.x.w,
			a.x.y*a.z.z*a.y.w + a.y.y*a.x.z*a.z.w + a.z.y*a.y.z*a.x.w - a.x.y*a.y.z*a.z.w - a.z.y*a.x.z*a.y.w - a.y.y*a.z.z*a.x.w },
			{ a.y.z*a.w.w*a.z.x + a.z.z*a.y.w*a.w.x + a.w.z*a.z.w*a.y.x - a.y.z*a.z.w*a.w.x - a.w.z*a.y.w*a.z.x - a.z.z*a.w.w*a.y.x,
			a.x.z*a.z.w*a.w.x + a.w.z*a.x.w*a.z.x + a.z.z*a.w.w*a.x.x - a.x.z*a.w.w*a.z.x - a.z.z*a.x.w*a.w.x - a.w.z*a.z.w*a.x.x,
			a.x.z*a.w.w*a.y.x + a.y.z*a.x.w*a.w.x + a.w.z*a.y.w*a.x.x - a.x.z*a.y.w*a.w.x - a.w.z*a.x.w*a.y.x - a.y.z*a.w.w*a.x.x,
			a.x.z*a.y.w*a.z.x + a.z.z*a.x.w*a.y.x + a.y.z*a.z.w*a.x.x - a.x.z*a.z.w*a.y.x - a.y.z*a.x.w*a.z.x - a.z.z*a.y.w*a.x.x },
			{ a.y.w*a.z.x*a.w.y + a.w.w*a.y.x*a.z.y + a.z.w*a.w.x*a.y.y - a.y.w*a.w.x*a.z.y - a.z.w*a.y.x*a.w.y - a.w.w*a.z.x*a.y.y,
			a.x.w*a.w.x*a.z.y + a.z.w*a.x.x*a.w.y + a.w.w*a.z.x*a.x.y - a.x.w*a.z.x*a.w.y - a.w.w*a.x.x*a.z.y - a.z.w*a.w.x*a.x.y,
			a.x.w*a.y.x*a.w.y + a.w.w*a.x.x*a.y.y + a.y.w*a.w.x*a.x.y - a.x.w*a.w.x*a.y.y - a.y.w*a.x.x*a.w.y - a.w.w*a.y.x*a.x.y,
			a.x.w*a.z.x*a.y.y + a.y.w*a.x.x*a.z.y + a.z.w*a.y.x*a.x.y - a.x.w*a.y.x*a.z.y - a.z.w*a.x.x*a.y.y - a.y.w*a.z.x*a.x.y },
			{ a.y.x*a.w.y*a.z.z + a.z.x*a.y.y*a.w.z + a.w.x*a.z.y*a.y.z - a.y.x*a.z.y*a.w.z - a.w.x*a.y.y*a.z.z - a.z.x*a.w.y*a.y.z,
			a.x.x*a.z.y*a.w.z + a.w.x*a.x.y*a.z.z + a.z.x*a.w.y*a.x.z - a.x.x*a.w.y*a.z.z - a.z.x*a.x.y*a.w.z - a.w.x*a.z.y*a.x.z,
			a.x.x*a.w.y*a.y.z + a.y.x*a.x.y*a.w.z + a.w.x*a.y.y*a.x.z - a.x.x*a.y.y*a.w.z - a.w.x*a.x.y*a.y.z - a.y.x*a.w.y*a.x.z,
			a.x.x*a.y.y*a.z.z + a.z.x*a.x.y*a.y.z + a.y.x*a.z.y*a.x.z - a.x.x*a.z.y*a.y.z - a.y.x*a.x.y*a.z.z - a.z.x*a.y.y*a.x.z } };
	}

	// Compute the determinant of a square matrix
	template<class T> T det(const mat<T, 2, 2> & a) { return a.x.x*a.y.y - a.x.y*a.y.x; }
	template<class T> T det(const mat<T, 3, 3> & a) { return a.x.x*(a.y.y*a.z.z - a.z.y*a.y.z) + a.x.y*(a.y.z*a.z.x - a.z.z*a.y.x) + a.x.z*(a.y.x*a.z.y - a.z.x*a.y.y); }
	template<class T> T det(const mat<T, 4, 4> & a)
	{
		return a.x.x*(a.y.y*a.z.z*a.w.w + a.w.y*a.y.z*a.z.w + a.z.y*a.w.z*a.y.w - a.y.y*a.w.z*a.z.w - a.z.y*a.y.z*a.w.w - a.w.y*a.z.z*a.y.w)
			+ a.x.y*(a.y.z*a.w.w*a.z.x + a.z.z*a.y.w*a.w.x + a.w.z*a.z.w*a.y.x - a.y.z*a.z.w*a.w.x - a.w.z*a.y.w*a.z.x - a.z.z*a.w.w*a.y.x)
			+ a.x.z*(a.y.w*a.z.x*a.w.y + a.w.w*a.y.x*a.z.y + a.z.w*a.w.x*a.y.y - a.y.w*a.w.x*a.z.y - a.z.w*a.y.x*a.w.y - a.w.w*a.z.x*a.y.y)
			+ a.x.w*(a.y.x*a.w.y*a.z.z + a.z.x*a.y.y*a.w.z + a.w.x*a.z.y*a.y.z - a.y.x*a.z.y*a.w.z - a.w.x*a.y.y*a.z.z - a.z.x*a.w.y*a.y.z);
	}

	// Compute the inverse of a square matrix
	template<class T, int N> mat<T, N, N> inv(const mat<T, N, N> & a) { return adj(a) / det(a); }

	// Compute the product of a matrix postmultiplied by a column vector
	template<class T, int M> vec<T, M> mul(const mat<T, M, 2> & a, const vec<T, 2> & b) { return a.x*b.x + a.y*b.y; }
	template<class T, int M> vec<T, M> mul(const mat<T, M, 3> & a, const vec<T, 3> & b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
	template<class T, int M> vec<T, M> mul(const mat<T, M, 4> & a, const vec<T, 4> & b) { return a.x*b.x + a.y*b.y + a.z*b.z + a.w*b.w; }

	// Compute the product of two matrices
	template<class T, int M, int N> mat<T, M, 2> mul(const mat<T, M, N> & a, const mat<T, N, 2> & b) { return{ mul(a, b.x), mul(a, b.y) }; }
	template<class T, int M, int N> mat<T, M, 3> mul(const mat<T, M, N> & a, const mat<T, N, 3> & b) { return{ mul(a, b.x), mul(a, b.y), mul(a, b.z) }; }
	template<class T, int M, int N> mat<T, M, 4> mul(const mat<T, M, N> & a, const mat<T, N, 4> & b) { return{ mul(a, b.x), mul(a, b.y), mul(a, b.z), mul(a, b.w) }; }

	// Retrieve a row (by index) from a matrix
	template<class T, int M> vec<T, 2> row(const mat<T, M, 2> & m, int i) { return{ m.x[i], m.y[i] }; }
	template<class T, int M> vec<T, 3> row(const mat<T, M, 3> & m, int i) { return{ m.x[i], m.y[i], m.z[i] }; }
	template<class T, int M> vec<T, 4> row(const mat<T, M, 4> & m, int i) { return{ m.x[i], m.y[i], m.z[i], m.w[i] }; }

	// Compute the tranpose of a matrix
	template<class T, int M> mat<T, M, 2> transpose(const mat<T, 2, M> & m) { return{ row(m, 0), row(m, 1) }; }
	template<class T, int M> mat<T, M, 3> transpose(const mat<T, 3, M> & m) { return{ row(m, 0), row(m, 1), row(m, 2) }; }
	template<class T, int M> mat<T, M, 4> transpose(const mat<T, 4, M> & m) { return{ row(m, 0), row(m, 1), row(m, 2), row(m, 3) }; }

	//////////////////////////////////
	// Quaternion algebra functions //
	//////////////////////////////////

	// Compute the conjugate of a quaternion represented by a four-element vector
	template<class T> vec<T, 4> qconj(const vec<T, 4> & q) { return{ -q.x, -q.y, -q.z, q.w }; }

	// Compute the inverse of a quaternion represented by a four-element vector
	template<class T> vec<T, 4> qinv(const vec<T, 4> & q) { return qconj(q) / mag2(q); }

	// Compute the product of two quaternions, each represented by a four-element vector
	template<class T> vec<T, 4> qmul(const vec<T, 4> & a, const vec<T, 4> & b) { return{ a.x*b.w + a.w*b.x + a.y*b.z - a.z*b.y, a.y*b.w + a.w*b.y + a.z*b.x - a.x*b.z, a.z*b.w + a.w*b.z + a.x*b.y - a.y*b.x, a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z }; }

	// Compute the vector part of the quaternion product q [1,0,0,0] q*
	template<class T> vec<T, 3> qxdir(const vec<T, 4> & q) { return{ q.w*q.w + q.x*q.x - q.y*q.y - q.z*q.z, (q.x*q.y + q.z*q.w) * 2, (q.z*q.x - q.y*q.w) * 2 }; }

	// Compute the vector part of the quaternion product q [0,1,0,0] q*
	template<class T> vec<T, 3> qydir(const vec<T, 4> & q) { return{ (q.x*q.y - q.z*q.w) * 2, q.w*q.w - q.x*q.x + q.y*q.y - q.z*q.z, (q.y*q.z + q.x*q.w) * 2 }; }

	// Compute the vector part of the quaternion product q [0,0,1,0] q*
	template<class T> vec<T, 3> qzdir(const vec<T, 4> & q) { return{ (q.z*q.x + q.y*q.w) * 2, (q.y*q.z - q.x*q.w) * 2, q.w*q.w - q.x*q.x - q.y*q.y + q.z*q.z }; }

	// Compute the vector part of the quaternion product q [v.x,v.y,v.z,0] q*
	template<class T> vec<T, 3> qtransform(const vec<T, 4> & q, const vec<T, 3> & v) { return qxdir(q)*v.x + qydir(q)*v.y + qzdir(q)*v.z; }

	////////////////////////
	// Functional helpers //
	////////////////////////


	// Form a scalar by applying function f to adjacent components of vector or matrix v
	template<class T, class F> T reduce(const vec<T, 2> & v, F f) { return f(v.x, v.y); }
	template<class T, class F> T reduce(const vec<T, 3> & v, F f) { return f(f(v.x, v.y), v.z); }
	template<class T, class F> T reduce(const vec<T, 4> & v, F f) { return f(f(f(v.x, v.y), v.z), v.w); }
	template<class T, int M, class F> T reduce(const mat<T, M, 2> & m, F f) { return f(reduce(m.x, f), reduce(m.y, f)); }
	template<class T, int M, class F> T reduce(const mat<T, M, 3> & m, F f) { return f(f(reduce(m.x, f), reduce(m.y, f)), reduce(m.z, f)); }
	template<class T, int M, class F> T reduce(const mat<T, M, 4> & m, F f) { return f(f(f(reduce(m.x, f), reduce(m.y, f)), reduce(m.z, f)), reduce(m.w, f)); }

	// Form a vector or matrix by applying function f to corresponding pairs of components from l and r
	template<class U, class T, class F> U zip(const vec<T, 2> & l, const vec<T, 2> & r, F f) { return{ f(l.x, r.x), f(l.y, r.y) }; }
	template<class U, class T, class F> U zip(const vec<T, 3> & l, const vec<T, 3> & r, F f) { return{ f(l.x, r.x), f(l.y, r.y), f(l.z, r.z) }; }
	template<class U, class T, class F> U zip(const vec<T, 4> & l, const vec<T, 4> & r, F f) { return{ f(l.x, r.x), f(l.y, r.y), f(l.z, r.z), f(l.w, r.w) }; }
	template<class U, class T, class F> U zip(const vec<T, 2> & l, const T & r, F f) { return{ f(l.x, r), f(l.y, r) }; }
	template<class U, class T, class F> U zip(const vec<T, 3> & l, const T & r, F f) { return{ f(l.x, r), f(l.y, r), f(l.z, r) }; }
	template<class U, class T, class F> U zip(const vec<T, 4> & l, const T & r, F f) { return{ f(l.x, r), f(l.y, r), f(l.z, r), f(l.w, r) }; }
	template<class U, class T, class F> U zip(const T & l, const vec<T, 2> & r, F f) { return{ f(l, r.x), f(l, r.y) }; }
	template<class U, class T, class F> U zip(const T & l, const vec<T, 3> & r, F f) { return{ f(l, r.x), f(l, r.y), f(l, r.z) }; }
	template<class U, class T, class F> U zip(const T & l, const vec<T, 4> & r, F f) { return{ f(l, r.x), f(l, r.y), f(l, r.z), f(l, r.w) }; }
	template<class U, class T, int M, class F> U zip(const mat<T, M, 2> & l, const mat<T, M, 2> & r, F f) { typedef typename U::C C; return{ zip<C>(l.x, r.x, f), zip<C>(l.y, r.y, f) }; }
	template<class U, class T, int M, class F> U zip(const mat<T, M, 3> & l, const mat<T, M, 3> & r, F f) { typedef typename U::C C; return{ zip<C>(l.x, r.x, f), zip<C>(l.y, r.y, f), zip<C>(l.z, r.z, f) }; }
	template<class U, class T, int M, class F> U zip(const mat<T, M, 4> & l, const mat<T, M, 4> & r, F f) { typedef typename U::C C; return{ zip<C>(l.x, r.x, f), zip<C>(l.y, r.y, f), zip<C>(l.z, r.z, f), zip<C>(l.w, r.w, f) }; }
	template<class U, class T, int M, class F> U zip(const mat<T, M, 2> & l, const T & r, F f) { typedef typename U::C C; return{ zip<C>(l.x, r, f), zip<C>(l.y, r, f) }; }
	template<class U, class T, int M, class F> U zip(const mat<T, M, 3> & l, const T & r, F f) { typedef typename U::C C; return{ zip<C>(l.x, r, f), zip<C>(l.y, r, f), zip<C>(l.z, r, f) }; }
	template<class U, class T, int M, class F> U zip(const mat<T, M, 4> & l, const T & r, F f) { typedef typename U::C C; return{ zip<C>(l.x, r, f), zip<C>(l.y, r, f), zip<C>(l.z, r, f), zip<C>(l.w, r, f) }; }
	template<class U, class T, int M, class F> U zip(const T & l, const mat<T, M, 2> & r, F f) { typedef typename U::C C; return{ zip<C>(l, r.x, f), zip<C>(l, r.y, f) }; }
	template<class U, class T, int M, class F> U zip(const T & l, const mat<T, M, 3> & r, F f) { typedef typename U::C C; return{ zip<C>(l, r.x, f), zip<C>(l, r.y, f), zip<C>(l, r.z, f) }; }
	template<class U, class T, int M, class F> U zip(const T & l, const mat<T, M, 4> & r, F f) { typedef typename U::C C; return{ zip<C>(l, r.x, f), zip<C>(l, r.y, f), zip<C>(l, r.z, f), zip<C>(l, r.w, f) }; }

	// Form a vector or matrix by applying function f to each element of vector or matrix v
	template<class T, class F> T map(const T & v, F f) { return zip<T>(v, elem_t<T>(), [f](elem_t<T> a, elem_t<T>) { return f(a); }); }

	// Ordering operators return a vector or matrix of bools, by comparing componentwise pairs
	template<class L, class R> comp_t<L, R> operator <  (const L & l, const R & r) { return zip<comp_t<L, R>>(l, r, [](elem_t<L> a, elem_t<R> b) { return a <  b; }); }
	template<class L, class R> comp_t<L, R> operator <= (const L & l, const R & r) { return zip<comp_t<L, R>>(l, r, [](elem_t<L> a, elem_t<R> b) { return a <= b; }); }
	template<class L, class R> comp_t<L, R> operator >(const L & l, const R & r) { return zip<comp_t<L, R>>(l, r, [](elem_t<L> a, elem_t<R> b) { return a >  b; }); }
	template<class L, class R> comp_t<L, R> operator >= (const L & l, const R & r) { return zip<comp_t<L, R>>(l, r, [](elem_t<L> a, elem_t<R> b) { return a >= b; }); }

	// Arithmetic and bitwise binary operators return a vector or matrix of the original type, by applying the operator to componentwise pairs
	template<class L, class R> arith_t<L, R> operator + (const L & l, const R & r) { return zip<arith_t<L, R>>(l, r, [](elem_t<L> a, elem_t<R> b) { return a + b; }); }
	template<class L, class R> arith_t<L, R> operator - (const L & l, const R & r) { return zip<arith_t<L, R>>(l, r, [](elem_t<L> a, elem_t<R> b) { return a - b; }); }
	template<class L, class R> arith_t<L, R> operator * (const L & l, const R & r) { return zip<arith_t<L, R>>(l, r, [](elem_t<L> a, elem_t<R> b) { return a*b; }); }
	template<class L, class R> arith_t<L, R> operator / (const L & l, const R & r) { return zip<arith_t<L, R>>(l, r, [](elem_t<L> a, elem_t<R> b) { return a / b; }); }
	template<class L, class R> arith_t<L, R> operator % (const L & l, const R & r) { return zip<arith_t<L, R>>(l, r, [](elem_t<L> a, elem_t<R> b) { return a%b; }); }
	template<class L, class R> arith_t<L, R> operator | (const L & l, const R & r) { return zip<arith_t<L, R>>(l, r, [](elem_t<L> a, elem_t<R> b) { return a | b; }); }
	template<class L, class R> arith_t<L, R> operator & (const L & l, const R & r) { return zip<arith_t<L, R>>(l, r, [](elem_t<L> a, elem_t<R> b) { return a&b; }); }
	template<class L, class R> arith_t<L, R> operator ^ (const L & l, const R & r) { return zip<arith_t<L, R>>(l, r, [](elem_t<L> a, elem_t<R> b) { return a^b; }); }
	
	/////////////////////////
	// Convenience aliases //
	/////////////////////////

	// These aliases resemble those used in HLSL, Cg, OpenCL, etc.
	namespace aliases
	{
		typedef vec<int32_t, 2> int2; typedef vec<uint32_t, 2> uint2; typedef vec<bool, 2> bool2; typedef vec<float, 2> float2; typedef vec<double, 2> double2;
		typedef vec<int32_t, 3> int3; typedef vec<uint32_t, 3> uint3; typedef vec<bool, 3> bool3; typedef vec<float, 3> float3; typedef vec<double, 3> double3;
		typedef vec<int32_t, 4> int4; typedef vec<uint32_t, 4> uint4; typedef vec<bool, 4> bool4; typedef vec<float, 4> float4; typedef vec<double, 4> double4;
		typedef mat<float, 2, 2> float2x2; typedef mat<float, 2, 3> float2x3; typedef mat<float, 2, 4> float2x4;
		typedef mat<float, 3, 2> float3x2; typedef mat<float, 3, 3> float3x3; typedef mat<float, 3, 4> float3x4;
		typedef mat<float, 4, 2> float4x2; typedef mat<float, 4, 3> float4x3; typedef mat<float, 4, 4> float4x4;
		typedef mat<double, 2, 2> double2x2; typedef mat<double, 2, 3> double2x3; typedef mat<double, 2, 4> double2x4;
		typedef mat<double, 3, 2> double3x2; typedef mat<double, 3, 3> double3x3; typedef mat<double, 3, 4> double3x4;
		typedef mat<double, 4, 2> double4x2; typedef mat<double, 4, 3> double4x3; typedef mat<double, 4, 4> double4x4;
	}
}

#endif
