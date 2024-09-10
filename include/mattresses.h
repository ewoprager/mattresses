#pragma once

#include <iostream>
#include <concepts>
#include <utility>

#include <macros.h>

// ----- Internal tools -----

namespace mattresses {

template <typename T, typename return_t, typename... param_ts>
concept function_c = requires (const T &function, param_ts&&... params) {
	{ function(std::forward<param_ts>(params)...) } -> std::convertible_to<return_t>;
};

} // namespace mattresses


// ----- Vectors -----
// -----
#define MATTRESSES_VECTOR_MUTATION_OPERATOR_OVERLOADS(_operator) \
vec &operator _operator (const vec &other){ \
	[&]<size_t... indices>(std::index_sequence<indices...>){ \
		(void([&](){ \
			(*this)[indices] _operator other[indices]; \
		}()), ...); \
	}(std::make_index_sequence<dimensionality>{}); \
	return *this; \
} \
vec &operator _operator (T f){ \
	[&]<size_t... indices>(std::index_sequence<indices...>){ \
		(void([&](){ \
			(*this)[indices] _operator f; \
		}()), ...); \
	}(std::make_index_sequence<dimensionality>{}); \
	return *this; \
}

template <size_t N, typename T=float>
struct vec {
	static constexpr size_t dimensionality = N;
	
	T data[dimensionality];
	
#include "vector_shared.inline"
};


// General vector functions
// -----
template <size_t N, typename T> inline std::ostream &operator<<(std::ostream &stream, const vec<N, T> &vector){
	stream << "(";
	for(int n=0; n<N - 1; n++) stream << vector[n] << ", ";
	stream << vector[N - 1] << ")";
	return stream;
}

template <size_t N, typename T> inline T Dot(const vec<N, T> &lhs, const vec<N, T> &rhs){
	return [&]<size_t... indices>(std::index_sequence<indices...>){
		return ((lhs[indices] * rhs[indices]) + ...);
	}(std::make_index_sequence<N>{});
}

template <size_t N1, size_t N2, typename T> inline vec<N1 + N2, T> operator|(const vec<N1, T> &lhs, const vec<N2, T> &rhs){
	vec<N1 + N2, T> ret;
	memcpy(&ret, &lhs, N1*sizeof(T));
	memcpy((T *)(&ret) + N1, &rhs, N2*sizeof(T));
	return ret;
}
template <size_t N, typename T> inline vec<N + 1, T> operator|(const vec<N, T> &lhs, T rhs){
	vec<N + 1, T> ret;
	memcpy(&ret, &lhs, N*sizeof(T));
	ret[N] = rhs;
	return ret;
}
template <size_t N, typename T> inline vec<N + 1, T> operator|(T lhs, const vec<N, T> &rhs){
	vec<N + 1, T> ret;
	ret[0] = lhs;
	memcpy((T *)(&ret) + 1, &rhs, N*sizeof(T));
	return ret;
}

#define MATTRESSES_VECTOR_BINARY_OPERATOR_OVERLOADS(_operator) \
template <size_t N, typename T> vec<N, T> operator _operator (const vec<N, T> &lhs, const vec<N, T> &rhs){ \
	return [&]<size_t... indices>(std::index_sequence<indices...>){ \
		return vec<N, T>{(void(indices), lhs[indices] _operator rhs[indices]) ...}; \
	}(std::make_index_sequence<N>{}); \
} \
template <size_t N, typename T> vec<N, T> operator _operator (const vec<N, T> &lhs, T rhs){ \
	return [&]<size_t... indices>(std::index_sequence<indices...>){ \
		return vec<N, T>{(void(indices), lhs[indices] _operator rhs) ...}; \
	}(std::make_index_sequence<N>{}); \
} \
template <size_t N, typename T> vec<N, T> operator _operator (T lhs, const vec<N, T> &rhs){ \
	return [&]<size_t... indices>(std::index_sequence<indices...>){ \
		return vec<N, T>{(void(indices), lhs _operator rhs[indices]) ...}; \
	}(std::make_index_sequence<N>{}); \
}
												
MATTRESSES_VECTOR_BINARY_OPERATOR_OVERLOADS(*)
MATTRESSES_VECTOR_BINARY_OPERATOR_OVERLOADS(+)
MATTRESSES_VECTOR_BINARY_OPERATOR_OVERLOADS(-)
MATTRESSES_VECTOR_BINARY_OPERATOR_OVERLOADS(/)

template <size_t N, typename T> inline bool operator==(const vec<N, T> &lhs, const vec<N, T> &rhs){
	return [&]<size_t... indices>(std::index_sequence<indices...>){
		return ((lhs[indices] == rhs[indices]) && ...);
	}(std::make_index_sequence<N>{});
}

#define MATTRESSES_VECTOR_SUBSET_N_REFERENCE_PARENT(N) \
struct subset##N##_reference_parent { \
	subset##N##_reference_parent(vec &_owner) : owner(_owner) {} \
	 \
	virtual void operator=(const vec<N, T> &other) = 0; \
	virtual operator vec<N, T> () const = 0; \
	 \
protected: \
	vec &owner; \
};

#define MATTRESSES_DEFINE_VECTOR_SUBSET_2(element1, element2) \
vec<2, T> element1##element2() const { return {element1, element2}; } \
struct element1##element2##_reference : public subset2_reference_parent { \
	element1##element2##_reference(vec &_owner) : subset2_reference_parent(_owner) {} \
	void operator=(const vec<2, T> &other) override { subset2_reference_parent::owner.element1 = other.x; subset2_reference_parent::owner.element2 = other.y; } \
	operator vec<2, T> () const override { return subset2_reference_parent::owner.element1##element2(); } \
}; \
element1##element2##_reference element1##element2##_r(){ return element1##element2##_reference(*this); }

#define MATTRESSES_DEFINE_VECTOR_SUBSET_3(element1, element2, element3) \
vec<3, T> element1##element2##element3() const { return {element1, element2, element3}; } \
struct element1##element2##element3##_reference : public subset3_reference_parent { \
	element1##element2##element3##_reference(vec &_owner) : subset3_reference_parent(_owner) {} \
	void operator=(const vec<3, T> &other) override { subset3_reference_parent::owner.element1 = other.x; subset3_reference_parent::owner.element2 = other.y; subset3_reference_parent::owner.element3 = other.z; } \
	operator vec<3, T> () const override { return subset3_reference_parent::owner.element1##element2##element3(); } \
}; \
element1##element2##element3##_reference element1##element2##element3##_r(){ return element1##element2##element3##_reference(*this); }


// 2-vector
// -----
template <typename T>
struct vec<2, T> {
	static constexpr size_t dimensionality = 2;
	
	T x, y;
	
	vec Crossed() const {
		return {-y, x};
	}
	double Angle() const {
		return atan2(y, x);
	}
	static vec UnitInDirection(double angle) {
		return {(T)cos(angle), (T)sin(angle)};
	}
	static vec RandomUnit(size_t precision=1000) {
		return UnitInDirection(2.0 * M_PI * double(rand() % precision) / double(precision));
	}
	
#include "vector_shared.inline"
};


// 3-vector
// -----
template <typename T>
struct vec<3, T> {
	static constexpr size_t dimensionality = 3;
	
	T x, y, z;
	
	MATTRESSES_VECTOR_SUBSET_N_REFERENCE_PARENT(2)

	FOR_EACH_PAIR(MATTRESSES_DEFINE_VECTOR_SUBSET_2,
				  x, y,
				  y, x,
				  x, z,
				  z, x,
				  y, z,
				  z, y
	)
	
#include "vector_shared.inline"
};

template<typename T> vec<3, T> Cross(const vec<3, T> &lhs, const vec<3, T> &rhs){
	return {lhs.y * rhs.z - lhs.z * rhs.y, lhs.z * rhs.x - lhs.x * rhs.z, lhs.x * rhs.y - lhs.y * rhs.x};
}


// 4-vector
// -----
template <typename T>
struct vec<4, T>{
	static constexpr size_t dimensionality = 4;
	
	T x, y, z, w;
	
	
	MATTRESSES_VECTOR_SUBSET_N_REFERENCE_PARENT(2)

	FOR_EACH_PAIR(MATTRESSES_DEFINE_VECTOR_SUBSET_2,
				  x, y,
				  y, x,
				  
				  x, z,
				  z, x,
				  
				  x, w,
				  w, x,
				  
				  y, z,
				  z, y,
				  
				  y, w,
				  w, y,
				  
				  z, w,
				  w, z
				  )

	MATTRESSES_VECTOR_SUBSET_N_REFERENCE_PARENT(3)

	FOR_EACH_TRIPLE(MATTRESSES_DEFINE_VECTOR_SUBSET_3,
				    x, y, z,
				    z, x, y,
				    y, z, x,
				    x, z, y,
				    y, x, z,
				    z, y, x,
				    
				    x, y, w,
				    w, x, y,
				    y, w, x,
				    x, w, y,
				    y, x, w,
				    w, y, x,
				    
				    x, w, z,
				    z, x, w,
				    w, z, x,
				    x, z, w,
				    w, x, z,
				    z, w, x,
				    
				    w, y, z,
				    z, w, y,
				    y, z, w,
				    w, z, y,
				    y, w, z,
				    z, y, w
					)
	
#include "vector_shared.inline"
};


// ----- Matrices -----
// -----

// N x M, meaning N rows & M columns
#define MATTRESSES_MATRIX_UNARY_OPERATOR_OVERLOADS(_operator) \
mat &operator _operator (const mat &other){ \
	for(int m=0; m<columnCount; m++) columns[m] _operator other[m]; \
	return *this; \
} \
mat &operator _operator (T f){ \
	for(int m=0; m<columnCount; m++) columns[m] _operator f; \
	return *this; \
}

#define MATTRESSES_MATRIX_BINARY_OPERATOR_OVERLOADS(_operator) \
friend mat operator _operator (const mat &lhs, const mat &rhs){ \
	mat ret; \
	for(int m=0; m<columnCount; m++) ret[m] = lhs[m] _operator rhs[m]; \
	return ret; \
} \
friend mat operator _operator (T lhs, const mat &rhs){ \
	mat ret; \
	for(int m=0; m<columnCount; m++) ret[m] = lhs _operator rhs[m]; \
	return ret; \
} \
friend mat operator _operator (const mat &lhs, T rhs){ \
	mat ret; \
	for(int m=0; m<columnCount; m++) ret[m] = lhs[m] _operator rhs; \
	return ret; \
}

template <size_t N, size_t M, typename T=float>
struct mat {
	static constexpr size_t rowCount = N;
	static constexpr size_t columnCount = M;
	
#include "matrix_shared.inline"
};


// General matrix functions
// -----
template <size_t N, size_t M, typename T> inline std::ostream &operator<<(std::ostream &stream, const mat<N, M, T> &matrix){
	stream << "\n";
	stream << "/ ";
	for(int m=0; m<M - 1; m++) stream << matrix[m][0] << "  ";
	stream << matrix[M - 1][0] << " \\\n";
	for(int n=1; n<N-1; n++){
		stream << "| ";
		for(int m=0; m<M - 1; m++) stream << matrix[m][n] << "  ";
		stream << matrix[M - 1][n] << " |\n";
	}
	stream << "\\ ";
	for(int m=0; m<M - 1; m++) stream << matrix[m][N - 1] << "  ";
	stream << matrix[M - 1][N - 1] << " /\n";
	return stream;
}

template <size_t M, typename T> inline std::ostream &operator<<(std::ostream &stream, const mat<1, M, T> &matrix){
	stream << "\n";
	stream << "[ ";
	for(int m=0; m<M - 1; m++) stream << matrix[m][0] << "  ";
	stream << matrix[M - 1][0] << " ]\n";
	return stream;
}

template <size_t N, size_t M, typename T> inline T Dot(const typename mat<N, M, T>::const_row &lhs, const vec<M, T> &rhs){
	T ret = 0.0f;
	for(int i=0; i<M; i++) ret += lhs[i] * rhs[i];
	return ret;
}

template <size_t N1, size_t MN, size_t M2, typename T> inline mat<N1, M2, T> operator &(const mat<N1, MN, T> &lhs, const mat<MN, M2, T> &rhs){
	mat<N1, M2, T> ret;
	for(int m=0; m<M2; m++) for(int n=0; n<N1; n++) ret[m][n] = Dot<N1, MN, T>(lhs(n), rhs[m]);
	return ret;
}

template <size_t N, size_t M, typename T> inline vec<N, T> operator &(const mat<N, M, T> &matrix, const vec<M, T> &vector){
	vec<N, T> ret;
	for(int n=0; n<N; n++) ret[n] = Dot<N, M, T>(matrix(n), vector);
	return ret;
}

template <size_t N, size_t M1, size_t M2, typename T> inline mat<N, M1 + M2, T> operator |(const mat<N, M1, T> &lhs, const mat<N, M2, T> &rhs){
    mat<N, M1 + M2, T> ret;
    memcpy(&ret[0][0], &lhs[0][0], sizeof(lhs));
    memcpy(&ret[M1][0], &rhs[0][0], sizeof(rhs));
    return ret;
}
template <size_t N, size_t M, typename T> inline mat<N, M + 1, T> operator |(const mat<N, M, T> &lhs, const vec<N, T> &rhs){
    mat<N, M + 1, T> ret;
    memcpy(&ret[0][0], &lhs[0][0], sizeof(lhs));
    memcpy(&ret[M][0], &rhs[0], sizeof(rhs));
    return ret;
}
template <size_t N, size_t M, typename T> inline mat<N, M + 1, T> operator |(const vec<N, T> &lhs, const mat<N, M, T> &rhs){
    mat<N, M + 1, T> ret;
    memcpy(&ret[0][0], &lhs[0], sizeof(lhs));
    memcpy(&ret[1][0], &rhs[0][0], sizeof(rhs));
    return ret;
}


// Square matrices
// -----
template <size_t N, typename T>
struct mat<N, N, T> {
	static constexpr size_t rowCount = N;
	static constexpr size_t columnCount = N;
	
#include "square_matrix_shared.inline"
	
#include "matrix_shared.inline"
};


// Single columns matrices
// -----
template <size_t N, typename T>
struct mat<N, 1, T> {
	static constexpr size_t rowCount = N;
	static constexpr size_t columnCount = 1;
	
	operator vec<N, T> () const {
		return columns[0];
	}
	
#include "matrix_shared.inline"
};


// Single row matrices
// -----
template <size_t M, typename T>
struct mat<1, M, T> {
	static constexpr size_t rowCount = 1;
	static constexpr size_t columnCount = M;
	
	operator vec<M, T> () const {
		return [this]<size_t... indices>(std::index_sequence<indices...>){
			return vec<M, T>{ columns[indices][0] ... };
		}(std::make_index_sequence<columnCount>{});
	}
	
#include "matrix_shared.inline"
};


// 4x4 matrices
// -----
template <typename T>
struct mat<4, 4, T> {
	static constexpr size_t rowCount = 4;
	static constexpr size_t columnCount = 4;
	
	mat ScaledExcludingTranslation(const vec<3, T> &scaling){
		const vec<3, T> savedTranslation = columns[3].xyz();
		mat ret = Scaling(scaling) & *this;
		ret[3].xyz_r() = savedTranslation;
		return ret;
	}
	mat Inverted() const {
		const T m00 = columns[0][0];
		const T m01 = columns[0][1];
		const T m02 = columns[0][2];
		const T m03 = columns[0][3];
		const T m10 = columns[1][0];
		const T m11 = columns[1][1];
		const T m12 = columns[1][2];
		const T m13 = columns[1][3];
		const T m20 = columns[2][0];
		const T m21 = columns[2][1];
		const T m22 = columns[2][2];
		const T m23 = columns[2][3];
		const T m30 = columns[3][0];
		const T m31 = columns[3][1];
		const T m32 = columns[3][2];
		const T m33 = columns[3][3];
		const T tmp_0  = m22 * m33;
		const T tmp_1  = m32 * m23;
		const T tmp_2  = m12 * m33;
		const T tmp_3  = m32 * m13;
		const T tmp_4  = m12 * m23;
		const T tmp_5  = m22 * m13;
		const T tmp_6  = m02 * m33;
		const T tmp_7  = m32 * m03;
		const T tmp_8  = m02 * m23;
		const T tmp_9  = m22 * m03;
		const T tmp_10 = m02 * m13;
		const T tmp_11 = m12 * m03;
		const T tmp_12 = m20 * m31;
		const T tmp_13 = m30 * m21;
		const T tmp_14 = m10 * m31;
		const T tmp_15 = m30 * m11;
		const T tmp_16 = m10 * m21;
		const T tmp_17 = m20 * m11;
		const T tmp_18 = m00 * m31;
		const T tmp_19 = m30 * m01;
		const T tmp_20 = m00 * m21;
		const T tmp_21 = m20 * m01;
		const T tmp_22 = m00 * m11;
		const T tmp_23 = m10 * m01;

		const T t0 = (tmp_0 * m11 + tmp_3 * m21 + tmp_4 * m31) - (tmp_1 * m11 + tmp_2 * m21 + tmp_5 * m31);
		const T t1 = (tmp_1 * m01 + tmp_6 * m21 + tmp_9 * m31) - (tmp_0 * m01 + tmp_7 * m21 + tmp_8 * m31);
		const T t2 = (tmp_2 * m01 + tmp_7 * m11 + tmp_10 * m31) - (tmp_3 * m01 + tmp_6 * m11 + tmp_11 * m31);
		const T t3 = (tmp_5 * m01 + tmp_8 * m11 + tmp_11 * m21) - (tmp_4 * m01 + tmp_9 * m11 + tmp_10 * m21);

		const T d = T(1.0) / (m00 * t0 + m10 * t1 + m20 * t2 + m30 * t3);

		return d * (mat){{
			{t0, t1, t2, t3},
			{
				((tmp_1 * m10 + tmp_2 * m20 + tmp_5 * m30) - (tmp_0 * m10 + tmp_3 * m20 + tmp_4 * m30)),
				((tmp_0 * m00 + tmp_7 * m20 + tmp_8 * m30) - (tmp_1 * m00 + tmp_6 * m20 + tmp_9 * m30)),
				((tmp_3 * m00 + tmp_6 * m10 + tmp_11 * m30) - (tmp_2 * m00 + tmp_7 * m10 + tmp_10 * m30)),
				((tmp_4 * m00 + tmp_9 * m10 + tmp_10 * m20) - (tmp_5 * m00 + tmp_8 * m10 + tmp_11 * m20))
			},
			{
				((tmp_12 * m13 + tmp_15 * m23 + tmp_16 * m33) - (tmp_13 * m13 + tmp_14 * m23 + tmp_17 * m33)),
				((tmp_13 * m03 + tmp_18 * m23 + tmp_21 * m33) - (tmp_12 * m03 + tmp_19 * m23 + tmp_20 * m33)),
				((tmp_14 * m03 + tmp_19 * m13 + tmp_22 * m33) - (tmp_15 * m03 + tmp_18 * m13 + tmp_23 * m33)),
				((tmp_17 * m03 + tmp_20 * m13 + tmp_23 * m23) - (tmp_16 * m03 + tmp_21 * m13 + tmp_22 * m23))
			},
			{
				((tmp_14 * m22 + tmp_17 * m32 + tmp_13 * m12) - (tmp_16 * m32 + tmp_12 * m12 + tmp_15 * m22)),
				((tmp_20 * m32 + tmp_12 * m02 + tmp_19 * m22) - (tmp_18 * m22 + tmp_21 * m32 + tmp_13 * m02)),
				((tmp_18 * m12 + tmp_23 * m32 + tmp_15 * m02) - (tmp_22 * m32 + tmp_14 * m02 + tmp_19 * m12)),
				((tmp_22 * m22 + tmp_16 * m02 + tmp_21 * m12) - (tmp_20 * m12 + tmp_23 * m22 + tmp_17 * m02))
			}
		}};
	}

	// now static methods
	static mat XRotation(double angle){
		const T c = cos(angle);
		const T s = sin(angle);
		return {{
			{1.0, 0.0, 0.0, 0.0},
			{0.0,   c,   s, 0.0},
			{0.0,  -s,   c, 0.0},
			{0.0, 0.0, 0.0, 1.0}
		}};
	}
	static mat YRotation(double angle){
		const T c = cos(angle);
		const T s = sin(angle);
		return {{
			{  c, 0.0,  -s, 0.0},
			{0.0, 1.0, 0.0, 0.0},
			{  s, 0.0,   c, 0.0},
			{0.0, 0.0, 0.0, 1.0}
		}};
	}
	static mat ZRotation(double angle){
		const T c = cos(angle);
		const T s = sin(angle);
		return {{
			{  c,   s, 0.0, 0.0},
			{ -s,   c, 0.0, 0.0},
			{0.0, 0.0, 1.0, 0.0},
			{0.0, 0.0, 0.0, 1.0}
		}};
	}
	static mat PerspectiveProjection(const double angleOfView, const double aspectRatio, const double frustumNear, const double frustumFar){
		const T f = tan(0.5 * (M_PI - angleOfView));
		const T frustumDepthInverse = 1.0 / (frustumNear - frustumFar);
		return {{
			{T(f / aspectRatio), 0.0,                                                  0.0,  0.0},
			{            0.0,   f,                                                  0.0,  0.0},
			{            0.0, 0.0,     T((frustumNear + frustumFar) * frustumDepthInverse), -1.0},
			{            0.0, 0.0, T(frustumNear * frustumFar * frustumDepthInverse * 2.0),  0.0}
		}};
	}
	static mat OrthographicProjection(const float left, const float right, const float bottom, const float top, const float zNear, const float zFar){
		const T dirInv = 1.0 / (right - left);
		const T heightInv = 1.0 / (top - bottom);
		const T zDeltaInv = 1.0 / (zFar - zNear);
		return {{
			{T(2.0 * dirInv), 			   0.0,        0.0,    T(-(right + left)*dirInv)},
			{	  	  0.0, T(2.0 * heightInv),        0.0, T(-(top + bottom)*heightInv)},
			{		  0.0,             0.0, T(-zDeltaInv),           T(zNear*zDeltaInv)},
			{		  0.0,		       0.0,        0.0,		                  1.0}
		}};
	}
	static mat Translation(const vec<3, T> &vector){
		return {{
			vec<4, T>::PositiveCartesianUnit(0),
			vec<4, T>::PositiveCartesianUnit(1),
			vec<4, T>::PositiveCartesianUnit(2),
			vector | T(1.0)
		}};
	}
	static mat LookAt(const vec<3, T> &position, const vec<3, T> &target, const vec<3, T> &unitUp){
		const vec<3, T> zAxis = (position - target).Normalised();
		const vec<3, T> xAxis = Cross(unitUp, zAxis).Normalised();
		const vec<3, T> yAxis = Cross(zAxis, xAxis).Normalised();
		return {{
			xAxis | T(0.0),
			yAxis | T(0.0),
			zAxis | T(0.0),
			position | T(1.0)
		}};
	}
	static mat Scaling(const vec<3, T> &scaling){
		return {{
			{scaling.x, 0.0, 0.0, 0.0},
			{0.0, scaling.y, 0.0, 0.0},
			{0.0, 0.0, scaling.z, 0.0},
			{0.0, 0.0, 0.0, 1.0}
		}};
	}
	
#include "square_matrix_shared.inline"
	
#include "matrix_shared.inline"
};
