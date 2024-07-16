#pragma once

#include <iostream>
#include <array>
#include <utility>
#include <concepts>
#include <type_traits>
#include <functional>

#include <macros.h>

// ----- Internal tools -----
namespace mattresses {

struct empty_t {};

template <size_t d> consteval bool DimensionalityHasX(){ return d >= 1 && d <= 4; }
template <size_t d> consteval bool DimensionalityHasY(){ return d >= 2 && d <= 4; }
template <size_t d> consteval bool DimensionalityHasZ(){ return d >= 3 && d <= 4; }
template <size_t d> consteval bool DimensionalityHasW(){ return d == 4; }

template <typename T, typename return_t, typename... param_ts>
concept function_c = requires (const T &function, param_ts&&... params) {
	{ function(std::forward<param_ts>(params)...) } -> std::convertible_to<return_t>;
};

} // namespace mattresses

#define MATTRESSES_VECTOR_MUTATION_OPERATOR_OVERLOADS(_operator) \
vec &operator _operator (const vec &other){ \
	[&]<size_t... indices>(std::index_sequence<indices...>){ \
		(void([&](){ \
			(*this)[indices] _operator other[indices]; \
		}()), ...); \
	}(std::make_index_sequence<N>{}); \
	return *this; \
} \
vec &operator _operator (T f){ \
	[&]<size_t... indices>(std::index_sequence<indices...>){ \
		(void([&](){ \
			(*this)[indices] _operator f; \
		}()), ...); \
	}(std::make_index_sequence<N>{}); \
	return *this; \
}

// ----- Vectors -----

template <size_t N, typename T=float> struct vec : std::array<T, N> {
	
	static constexpr bool have_x = mattresses::DimensionalityHasX<N>();
	
	std::conditional_t<have_x, T &, mattresses::empty_t> x = [this]() -> std::conditional_t<have_x, T &, mattresses::empty_t> {
		if constexpr (have_x) {
			return std::array<T, N>::at(0);
		} else {
			return mattresses::empty_t{};
		}
	}();
	
	static constexpr bool have_y = mattresses::DimensionalityHasY<N>();
	
	std::conditional_t<have_y, T &, mattresses::empty_t> y = [this]() -> std::conditional_t<have_y, T &, mattresses::empty_t> {
		if constexpr (have_y) {
			return std::array<T, N>::at(1);
		} else {
			return mattresses::empty_t{};
		}
	}();
	
	static constexpr bool have_z = mattresses::DimensionalityHasZ<N>();
	
	std::conditional_t<have_z, T &, mattresses::empty_t> z = [this]() -> std::conditional_t<have_z, T &, mattresses::empty_t> {
		if constexpr (have_z) {
			return std::array<T, N>::at(2);
		} else {
			return mattresses::empty_t{};
		}
	}();
	
	static constexpr bool have_w = mattresses::DimensionalityHasW<N>();
	
	std::conditional_t<have_w, T &, mattresses::empty_t> w = [this]() -> std::conditional_t<have_w, T &, mattresses::empty_t> {
		if constexpr (have_w) {
			return std::array<T, N>::at(3);
		} else {
			return mattresses::empty_t{};
		}
	}();
	
	T SqMag() const {
		return [&]<size_t... indices>(std::index_sequence<indices...>){
			return (((*this)[indices] * (*this)[indices]) + ...);
		}(std::make_index_sequence<N>{});
	}
	
	vec Normalised() const { return *this / sqrt(SqMag()); }
	
	
	// 2-dimensional vectors
	// -----
	std::conditional_t<N == 2, vec, mattresses::empty_t> Crossed() const {
		if constexpr (N == 2) {
			return vec{-y, x};
		} else {
			return mattresses::empty_t{};
		}
	}
	std::conditional_t<N == 2, double, mattresses::empty_t> Angle() const {
		if constexpr (N == 2) {
			return atan2(y, x);
		} else {
			return mattresses::empty_t{};
		}
	}
	static std::conditional_t<N == 2, vec, mattresses::empty_t> UnitInDirection(double angle) {
		if constexpr (N == 2) {
			return {(T)cos(angle), (T)sin(angle)};
		} else {
			return mattresses::empty_t{};
		}
	}
	static std::conditional_t<N == 2, vec, mattresses::empty_t> RandomUnit(size_t precision=1000) {
		if constexpr (N == 2) {
			return UnitInDirection(2.0 * M_PI * double(rand() % precision) / double(precision));
		} else {
			return mattresses::empty_t{};
		}
	}
	
	MATTRESSES_VECTOR_MUTATION_OPERATOR_OVERLOADS(*=)
	MATTRESSES_VECTOR_MUTATION_OPERATOR_OVERLOADS(+=)
	MATTRESSES_VECTOR_MUTATION_OPERATOR_OVERLOADS(-=)
	MATTRESSES_VECTOR_MUTATION_OPERATOR_OVERLOADS(/=)
	
	template <typename other_t> operator vec<N, other_t>() const {
		return [&]<size_t... indices>(std::index_sequence<indices...>){
			return vec<N, other_t>{ static_cast<other_t>((*this)[indices]) ... };
		}(std::make_index_sequence<N>{});
	}
	
	template <typename function_t>
	requires mattresses::function_c<function_t, T, T>
	vec &ForEach(const function_t &function){
		[&]<size_t... indices>(std::index_sequence<indices...>){
			(void([&](){
				(*this)[indices] = function((*this)[indices]);
			}()), ...);
		}(std::make_index_sequence<N>{});
		return *this;
	}
	
	static vec Zero(){
		return [&]<size_t... indices>(std::index_sequence<indices...>){
			return vec{ (void(indices), T(0)) ... };
		}(std::make_index_sequence<N>{});
	}
	
	static vec Ones(){
		return [&]<size_t... indices>(std::index_sequence<indices...>){
			return vec{ (void(indices), T(1)) ... };
		}(std::make_index_sequence<N>{});
	}
	
	static vec PositiveCartesianUnit(size_t n){ vec ret = vec::Zero(); ret[n] = T(1); return ret; }
	
	template <typename function_t, typename... param_ts>
	requires mattresses::function_c<function_t, T, param_ts...>
	static vec GenerateEach(const function_t &function, param_ts&&... params){
		return [&]<size_t... indices>(std::index_sequence<indices...>){
			return vec{ (void(indices), function(std::forward<param_ts>(params)...)) ... };
		}(std::make_index_sequence<N>{});
	}
};

template <size_t N, typename T> std::ostream &operator<<(std::ostream &stream, const vec<N, T> &vector){
	stream << "(";
	for(int n=0; n<N - 1; ++n) stream << vector[n] << ", ";
	stream << vector[N - 1] << ")";
	return stream;
}

template <size_t N, typename T> T Dot(const vec<N, T> &lhs, const vec<N, T> &rhs){
	return [&]<size_t... indices>(std::index_sequence<indices...>){
		return ((lhs[indices] * rhs[indices]) + ...);
	}(std::make_index_sequence<N>{});
}

template <size_t N1, size_t N2, typename T> vec<N1 + N2, T> operator|(const vec<N1, T> &lhs, const vec<N2, T> &rhs){
	vec<N1 + N2, T> ret;
	memcpy(ret.data(), lhs.data(), N1 * sizeof(T));
	memcpy((T *)(ret.data()) + N1, rhs.data(), N2 * sizeof(T));
	return ret;
}
template <size_t N, typename T> vec<N + 1, T> operator|(const vec<N, T> &lhs, T rhs){
	vec<N + 1, T> ret;
	memcpy(ret.data(), lhs.data(), N * sizeof(T));
	ret[N] = rhs;
	return ret;
}
template <size_t N, typename T> vec<N + 1, T> operator|(T lhs, const vec<N, T> &rhs){
	vec<N + 1, T> ret;
	ret[0] = lhs;
	memcpy((T *)(ret.data()) + 1, rhs.data(), N * sizeof(T));
	return ret;
}

#define DEFINE_VECTOR_BINARY_OPERATOR_OVERLOADS(_operator) \
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
												
DEFINE_VECTOR_BINARY_OPERATOR_OVERLOADS(*)
DEFINE_VECTOR_BINARY_OPERATOR_OVERLOADS(+)
DEFINE_VECTOR_BINARY_OPERATOR_OVERLOADS(-)
DEFINE_VECTOR_BINARY_OPERATOR_OVERLOADS(/)

template <size_t N, typename T> bool operator==(const vec<N, T> &lhs, const vec<N, T> &rhs){
	return [&]<size_t... indices>(std::index_sequence<indices...>){
		return ((lhs[indices] == rhs[indices]) && ...);
	}(std::make_index_sequence<N>{});
}

template <typename T> vec<3, T> Cross(const vec<3, T> &lhs, const vec<3, T> &rhs){
	return {lhs.y * rhs.z - lhs.z * rhs.y, lhs.z * rhs.x - lhs.x * rhs.z, lhs.x * rhs.y - lhs.y * rhs.x};
}


// ----- Quaternions -----

template <typename T> struct quat : public vec<4, T> {

};

//#define _N_VECTOR_ELEMENT(element) T element;
//#define _N_VECTOR_OSTREAM(element) vector.element
//#define _N_VECTOR_DOT(element) lhs.element * rhs.element
//#define _N_VECTOR_EQUAL(element) (lhs.element == rhs.element)
//#define _N_VECTOR_ZERO(element) T(0)
//#define _N_VECTOR_ONE(element) T(1)
//#define _N_VECTOR_SQ(element) element * element
//#define _N_VECTOR_FUNCTION(element) function()
//#define _N_VECTOR_CAST(element) static_cast<other_t>(element)
//
//#define _N_VECTOR_MUTATION_OPERATION(_operator, element) element _operator other.element;
//#define _N_VECTOR_FLOAT_MUTATION_OPERATION(_operator, element) element _operator f;
//#define DEFINE_N_VECTOR_MUTATION_OPERATOR_OVERLOADS(_operator, vector_elements...) \
//vec &operator _operator (const vec &other){ \
//	FOR_EACH_EX(_N_VECTOR_MUTATION_OPERATION, _operator, vector_elements) \
//	return *this; \
//} \
//vec &operator _operator (T f){ \
//	FOR_EACH_EX(_N_VECTOR_FLOAT_MUTATION_OPERATION, _operator, vector_elements) \
//	return *this; \
//}
//
//#define START_DEFINE_N_VECTOR_STRUCT_SPECIALISATION(size_t_size, vector_elements...) \
//template <typename T> struct vec<size_t_size, T> { \
//	const T &operator[](size_t n) const { return ((T *)this)[n]; } \
//	T &operator[](size_t n){ return ((T *)this)[n]; } \
//\
//	FOR_EACH(_N_VECTOR_ELEMENT, vector_elements) \
//\
//	T SqMag() const { return CUSTOM_SEPARATED_FOR_EACH(_N_VECTOR_SQ, +, vector_elements); } \
//    vec Normalised() const { vec ret = *this; ret /= sqrtf(ret.SqMag()); return ret; } \
//\
//	DEFINE_N_VECTOR_MUTATION_OPERATOR_OVERLOADS(*=, vector_elements) \
//	DEFINE_N_VECTOR_MUTATION_OPERATOR_OVERLOADS(+=, vector_elements) \
//	DEFINE_N_VECTOR_MUTATION_OPERATOR_OVERLOADS(-=, vector_elements) \
//	DEFINE_N_VECTOR_MUTATION_OPERATOR_OVERLOADS(/=, vector_elements) \
//\
//	template <typename other_t> operator vec<size_t_size, other_t>() const { return {COMMA_SEPARATED_FOR_EACH(_N_VECTOR_CAST, vector_elements)}; } \
//\
//	static vec Zero(){ return {COMMA_SEPARATED_FOR_EACH(_N_VECTOR_ZERO, vector_elements)}; } \
//	static vec One(){ return {COMMA_SEPARATED_FOR_EACH(_N_VECTOR_ONE, vector_elements)}; } \
//	static vec PositiveCartesianUnit(uint8_t n){ vec ret = vec::Zero(); ret[n] = T(1); return ret; } \
//	template <typename T_function_t> static vec FromFunction(T_function_t function){ return {COMMA_SEPARATED_FOR_EACH(_N_VECTOR_FUNCTION, vector_elements)}; }
//
//#define FINISH_DEFINE_N_VECTOR_STRUCT_SPECIALISATION(uint8_t_size, vector_elements...) \
//};
//
//#define _N_VECTOR_BINARY_OPERATION(_operator, element) lhs.element _operator rhs.element
//#define _N_FLOAT_VECTOR_BINARY_OPERATION(_operator, element) lhs _operator rhs.element
//#define _N_VECTOR_FLOAT_BINARY_OPERATION(_operator, element) lhs.element _operator rhs
//#define _N_VECTOR_BINARY_OPERATOR_OVERLOADS(_operator, uint8_t_size, vector_elements...) \
//template <typename T> inline vec<uint8_t_size, T> operator _operator (const vec<uint8_t_size, T> &lhs, const vec<uint8_t_size, T> &rhs){ \
//	return { COMMA_SEPARATED_FOR_EACH_EX(_N_VECTOR_BINARY_OPERATION, _operator, vector_elements) }; \
//} \
//template <typename T> inline vec<uint8_t_size, T> operator _operator (T lhs, const vec<uint8_t_size, T> &rhs){ \
//	return { COMMA_SEPARATED_FOR_EACH_EX(_N_FLOAT_VECTOR_BINARY_OPERATION, _operator, vector_elements) }; \
//} \
//template <typename T> inline vec<uint8_t_size, T> operator _operator (const vec<uint8_t_size, T> &lhs, T rhs){ \
//	return { COMMA_SEPARATED_FOR_EACH_EX(_N_VECTOR_FLOAT_BINARY_OPERATION, _operator, vector_elements) }; \
//}
//
//#define DEFINE_N_VECTOR_GLOBAL_SPECIALISATIONS(uint8_t_size, vector_elements...) \
//template <typename T> inline T Dot(const vec<uint8_t_size, T> &lhs, const vec<uint8_t_size, T> &rhs){ \
//   return CUSTOM_SEPARATED_FOR_EACH(_N_VECTOR_DOT, +, vector_elements); \
//} \
//template <typename T> inline std::ostream &operator<<(std::ostream &stream, const vec<uint8_t_size, T> &vector){ \
//	stream << "(" << CUSTOM_SEPARATED_FOR_EACH(_N_VECTOR_OSTREAM, << ", " <<, vector_elements) << ")"; \
//	return stream; \
//} \
//template <typename T> inline bool operator==(const vec<uint8_t_size, T> &lhs, const vec<uint8_t_size, T> &rhs){ \
//	return CUSTOM_SEPARATED_FOR_EACH(_N_VECTOR_EQUAL, &&, vector_elements); \
//} \
//_N_VECTOR_BINARY_OPERATOR_OVERLOADS(*, uint8_t_size, vector_elements) \
//_N_VECTOR_BINARY_OPERATOR_OVERLOADS(+, uint8_t_size, vector_elements) \
//_N_VECTOR_BINARY_OPERATOR_OVERLOADS(-, uint8_t_size, vector_elements) \
//_N_VECTOR_BINARY_OPERATOR_OVERLOADS(/, uint8_t_size, vector_elements)
//
//
//
//#define START_DEFINE_N_VECTOR_SPECIALISATION(uint8_t_size, vector_elements...) \
//START_DEFINE_N_VECTOR_STRUCT_SPECIALISATION(uint8_t_size, vector_elements)
//
//#define FINISH_DEFINE_N_VECTOR_SPECIALISATION(uint8_t_size, vector_elements...) \
//FINISH_DEFINE_N_VECTOR_STRUCT_SPECIALISATION(uint8_t_size, vector_elements) \
//DEFINE_N_VECTOR_GLOBAL_SPECIALISATIONS(uint8_t_size, vector_elements)
//
//
//#define DEFINE_VECTOR_SUBSET_N_REFERENCE_PARENT(N) \
//struct subset##N##_reference_parent { \
//	subset##N##_reference_parent(vec &_owner) : owner(_owner) {} \
//	 \
//	virtual void operator=(const vec<N, T> &other) = 0; \
//	virtual operator vec<N, T> () const = 0; \
//	 \
//protected: \
//	vec &owner; \
//};
//
//#define _VECTOR_DEFINE_SUBSET_2(element1, element2) \
//vec<2, T> element1##element2() const { return {element1, element2}; } \
//struct element1##element2##_reference : public subset2_reference_parent { \
//	element1##element2##_reference(vec &_owner) : subset2_reference_parent(_owner) {} \
//	void operator=(const vec<2, T> &other) override { subset2_reference_parent::owner.element1 = other.x; subset2_reference_parent::owner.element2 = other.y; } \
//	operator vec<2, T> () const override { return subset2_reference_parent::owner.element1##element2(); } \
//}; \
//element1##element2##_reference element1##element2##_r(){ return element1##element2##_reference(*this); }
//
//#define _VECTOR_DEFINE_SUBSET_3(element1, element2, element3) \
//vec<3, T> element1##element2##element3() const { return {element1, element2, element3}; } \
//struct element1##element2##element3##_reference : public subset3_reference_parent { \
//	element1##element2##element3##_reference(vec &_owner) : subset3_reference_parent(_owner) {} \
//	void operator=(const vec<3, T> &other) override { subset3_reference_parent::owner.element1 = other.x; subset3_reference_parent::owner.element2 = other.y; subset3_reference_parent::owner.element3 = other.z; } \
//	operator vec<3, T> () const override { return subset3_reference_parent::owner.element1##element2##element3(); } \
//}; \
//element1##element2##element3##_reference element1##element2##element3##_r(){ return element1##element2##element3##_reference(*this); }
//
//

//
//START_DEFINE_N_VECTOR_SPECIALISATION(3, x, y, z)
//
//DEFINE_VECTOR_SUBSET_N_REFERENCE_PARENT(2)
//
//FOR_EACH_PAIR(_VECTOR_DEFINE_SUBSET_2,
//			  x, y,
//			  y, x,
//			  x, z,
//			  z, x,
//			  y, z,
//			  z, y
//)
//
//
//START_DEFINE_N_VECTOR_SPECIALISATION(4, x, y, z, w)
//
//DEFINE_VECTOR_SUBSET_N_REFERENCE_PARENT(2)
//
//FOR_EACH_PAIR(_VECTOR_DEFINE_SUBSET_2,
//			  x, y,
//			  y, x,
//			  
//			  x, z,
//			  z, x,
//			  
//			  x, w,
//			  w, x,
//			  
//			  y, z,
//			  z, y,
//			  
//			  y, w,
//			  w, y,
//			  
//			  z, w,
//			  w, z
//)
//
//DEFINE_VECTOR_SUBSET_N_REFERENCE_PARENT(3)
//
//FOR_EACH_TRIPLE(_VECTOR_DEFINE_SUBSET_3,
//				x, y, z,
//				z, x, y,
//				y, z, x,
//				x, z, y,
//				y, x, z,
//				z, y, x,
//				
//				x, y, w,
//				w, x, y,
//				y, w, x,
//				x, w, y,
//				y, x, w,
//				w, y, x,
//				
//				x, w, z,
//				z, x, w,
//				w, z, x,
//				x, z, w,
//				w, x, z,
//				z, w, x,
//				
//				w, y, z,
//				z, w, y,
//				y, z, w,
//				w, z, y,
//				y, w, z,
//				z, y, w
//)
//
//FINISH_DEFINE_N_VECTOR_SPECIALISATION(4, x, y, z, w)
//




// ----- Matrices -----

// N x M, meaning N rows & M columns
//#define _MATRIX_TEMPLATE_PARAMETER(param_name) uint8_t param_name
//#define START_MATRIX_TEMPLATE(N, M, ...) \
//template < __VA_OPT__( COMMA_SEPARATED_FOR_EACH(_MATRIX_TEMPLATE_PARAMETER, __VA_ARGS__) ) __VA_OPT__(,) typename T=float> struct mat { \
//_CONTINUE_MATRIX_TEMPLATE(N, M)

//#define START_MATRIX_TEMPLATE_SPECIALISATION(N, M, ...) \
// \
//_CONTINUE_MATRIX_TEMPLATE(N, M)

#define _MATRIX_MUTATION_OPERATOR_OVERLOADS(_operator) \
mat &operator _operator (const mat &other){ \
	[&]<size_t... indices>(std::index_sequence<indices...>){ \
		(void([&](){ \
			(*this)[indices] _operator other[indices]; \
		}()), ...); \
	}(std::make_index_sequence<M>{}); \
	return *this; \
} \
mat &operator _operator (T f){ \
	[&]<size_t... indices>(std::index_sequence<indices...>){ \
		(void([&](){ \
			(*this)[indices] _operator f; \
		}()), ...); \
	}(std::make_index_sequence<M>{}); \
	return *this; \
}

#define _MATRIX_BINARY_OPERATOR_OVERLOADS(_operator) \
friend mat operator _operator (const mat &lhs, const mat &rhs){ \
	return [&]<size_t... indices>(std::index_sequence<indices...>){ \
		return vec{(void(indices), lhs[indices] _operator rhs[indices]) ...}; \
	}(std::make_index_sequence<M>{}); \
} \
friend mat operator _operator (T lhs, const mat &rhs){ \
	return [&]<size_t... indices>(std::index_sequence<indices...>){ \
		return vec{(void(indices), lhs _operator rhs[indices]) ...}; \
	}(std::make_index_sequence<M>{}); \
} \
friend mat operator _operator (const mat &lhs, T rhs){ \
	return [&]<size_t... indices>(std::index_sequence<indices...>){ \
		return vec{(void(indices), lhs[indices] _operator rhs) ...}; \
	}(std::make_index_sequence<M>{}); \
}

template <size_t N, size_t M, typename T=float> struct mat : public std::array<vec<N, T>, M> {

	struct const_row {
		const_row(const mat &_matrix, uint8_t _n) : matrix(_matrix), n(_n) {}

		const T &operator[](uint8_t m) const { return matrix[m][n]; }

        operator vec<M, T> () const {
            vec<M, T> ret;
            for(int m=0; m<M; m++) ret[m] = matrix[m][n];
            return ret;
        }

		friend std::ostream &operator<<(std::ostream &stream, const const_row &matrix_row){
			stream << "[";
			for(int m=0; m<M - 1; m++) stream << matrix_row[m] << ", ";
			stream << matrix_row[M - 1] << "]";
			return stream;
		}

	private:
		const mat &matrix;
		const uint8_t n;
	};
	struct row {
		row(mat &_matrix, uint8_t _n) : matrix(_matrix), n(_n) {}

		row &operator=(const vec<M, T> &col){
			for(int i=0; i<M; i++) matrix[i][n] = col[i];
			return *this;
		}

		T &operator[](uint8_t m) const { return matrix[m][n]; }

		operator const_row() const { return const_row(matrix, n); }
        operator vec<M, T> () const {
            vec<M, T> ret;
            for(int m=0; m<M; m++) ret[m] = matrix[m][n];
            return ret;
        }

		friend std::ostream &operator<<(std::ostream &stream, const row &matrix_row){
			stream << "[";
			for(int m=0; m<M - 1; m++) stream << matrix_row[m] << ", ";
			stream << matrix_row[M - 1] << "]";
			return stream;
		}

	private:
		mat &matrix;
		const uint8_t n;
	};

	const_row operator()(uint8_t n) const { return const_row(*this, n); }
	row operator()(uint8_t n){ return row(*this, n); }

	_MATRIX_MUTATION_OPERATOR_OVERLOADS(*=)
	_MATRIX_MUTATION_OPERATOR_OVERLOADS(+=)
	_MATRIX_MUTATION_OPERATOR_OVERLOADS(-=)
	_MATRIX_MUTATION_OPERATOR_OVERLOADS(/=)

	mat<M, N, T> Transposed() const {
		mat<M, N, T> ret;
		for(int i=0; i<M; i++) ret(i) = (*this)[i];
		return ret;
	}

	_MATRIX_BINARY_OPERATOR_OVERLOADS(+)
	_MATRIX_BINARY_OPERATOR_OVERLOADS(-)
	_MATRIX_BINARY_OPERATOR_OVERLOADS(*)
	_MATRIX_BINARY_OPERATOR_OVERLOADS(/)

	static mat Zeros(){
		return [&]<size_t... indices>(std::index_sequence<indices...>){
			return mat{ (void(indices), vec<N, T>::Zero()) ... };
		}(std::make_index_sequence<M>{});
	}
	
	// Single columns matrices
	// -----
	operator vec<N, T> () const requires(M == 1) {
		return (*this)[0];
	}
	
	// Single row matrices
	// -----
	operator vec<M, T> () const requires(N == 1) {
		return [this]<size_t... indices>(std::index_sequence<indices...>){
			return vec<M, T>{ (*this)[indices][0] ... };
		}(std::make_index_sequence<M>{});
	}
	
	// Square matrices
	// -----
	mat &TransposeInPlace() requires (M == N) {
		for(int m=0; m<N; m++) for(int n=0; n<N; n++){
			if(m == n) continue;
			const T temp = (*this)[m][n];
			(*this)[m][n] = (*this)[n][m];
			(*this)[n][m] = temp;
		}
		return *this;
	}
	static mat Identity() requires (M == N) {
		mat ret;
		for(int i=0; i<N; i++) ret[i] = vec<N, T>::PositiveCartesianUnit(i);
		return ret;
	}
	
	
	// 4x4 matrices
	// -----
//	mat ScaledExcludingTranslation(const vec<3, T> &scaling) const requires (M == 4 && N == 4) {
//		const vec<3, T> savedTranslation = columns[3].xyz();
//		mat ret = Scaling(scaling) & *this;
//		ret[3].xyz_r() = savedTranslation;
//		return ret;
//	}
	mat Inverted() const requires (M == 4 && N == 4) {
		const T m00 = (*this)[0][0];
		const T m01 = (*this)[0][1];
		const T m02 = (*this)[0][2];
		const T m03 = (*this)[0][3];
		const T m10 = (*this)[1][0];
		const T m11 = (*this)[1][1];
		const T m12 = (*this)[1][2];
		const T m13 = (*this)[1][3];
		const T m20 = (*this)[2][0];
		const T m21 = (*this)[2][1];
		const T m22 = (*this)[2][2];
		const T m23 = (*this)[2][3];
		const T m30 = (*this)[3][0];
		const T m31 = (*this)[3][1];
		const T m32 = (*this)[3][2];
		const T m33 = (*this)[3][3];
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
		
		const T d = T(1) / (m00 * t0 + m10 * t1 + m20 * t2 + m30 * t3);
		
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
	
	static mat XRotation(double angle) requires (M == 4 && N == 4) {
		const T c = cos(angle);
		const T s = sin(angle);
		return {{
			{1.0, 0.0, 0.0, 0.0},
			{0.0,   c,   s, 0.0},
			{0.0,  -s,   c, 0.0},
			{0.0, 0.0, 0.0, 1.0}
		}};
	}
	static mat YRotation(double angle) requires (M == 4 && N == 4) {
		const T c = cos(angle);
		const T s = sin(angle);
		return {{
			{  c, 0.0,  -s, 0.0},
			{0.0, 1.0, 0.0, 0.0},
			{  s, 0.0,   c, 0.0},
			{0.0, 0.0, 0.0, 1.0}
		}};
	}
	static mat ZRotation(double angle) requires (M == 4 && N == 4) {
		const T c = cos(angle);
		const T s = sin(angle);
		return {{
			{  c,   s, 0.0, 0.0},
			{ -s,   c, 0.0, 0.0},
			{0.0, 0.0, 1.0, 0.0},
			{0.0, 0.0, 0.0, 1.0}
		}};
	}
	static mat PerspectiveProjection(const double angleOfView, const double aspectRatio, const double frustumNear, const double frustumFar) requires (M == 4 && N == 4) {
		const T f = tan(0.5 * (M_PI - angleOfView));
		const T frustumDepthInverse = 1.0 / (frustumNear - frustumFar);
		return {{
			{T(f / aspectRatio), 0.0,                                                  0.0,  0.0},
			{            0.0,   f,                                                  0.0,  0.0},
			{            0.0, 0.0,     T((frustumNear + frustumFar) * frustumDepthInverse), -1.0},
			{            0.0, 0.0, T(frustumNear * frustumFar * frustumDepthInverse * 2.0),  0.0}
		}};
	}
	static mat OrthographicProjection(const float left, const float right, const float bottom, const float top, const float zNear, const float zFar) requires (M == 4 && N == 4) {
		const T dirInv = 1.0 / (right - left);
		const T heightInv = 1.0 / (top - bottom);
		const T zDeltaInv = 1.0 / (zFar - zNear);
		return {{
			{T(2.0 * dirInv), 			   	 0.0,       	0.0, T(-(right + left)*dirInv)   },
			{	  	 	 0.0, T(2.0 * heightInv),       	0.0, T(-(top + bottom)*heightInv)},
			{		 	 0.0,             	 0.0, T(-zDeltaInv), T(zNear*zDeltaInv)          },
			{		 	 0.0,		         0.0,       	0.0,						  1.0}
		}};
	}
	static mat Translation(const vec<3, T> &vector) requires (M == 4 && N == 4) {
		return {{
			vec<4, T>::PositiveCartesianUnit(0),
			vec<4, T>::PositiveCartesianUnit(1),
			vec<4, T>::PositiveCartesianUnit(2),
			vector | T(1.0)
		}};
	}
	static mat LookAt(const vec<3, T> &position, const vec<3, T> &target, const vec<3, T> &unitUp) requires (M == 4 && N == 4) {
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
	static mat Scaling(const vec<3, T> &scaling) requires (M == 4 && N == 4) {
		return {{
			{scaling.x, 0.0, 0.0, 0.0},
			{0.0, scaling.y, 0.0, 0.0},
			{0.0, 0.0, scaling.z, 0.0},
			{0.0, 0.0, 0.0, 1.0}
		}};
	}

	static mat AxisRotation(const vec<3, T> &unitAxis, const double angle) requires (M == 4 && N == 4) {
		const vec<3, T> &u = unitAxis;
		const T c = cos(angle);
		const T s = sin(angle);
		const T t = T(1) - c;
		return {{
			{t * u.x * u.x + c, t * u.x * u.y + s * u.z, t * u.x * u.z - s * u.y, 0.0},
			{t * u.x * u.y - s * u.z, t * u.y * u.y + c, t * u.y * u.z + s * u.x, 0.0},
			{t * u.x * u.z + s * u.y, t * u.y * u.z - s * u.x, t * u.z * u.z + c, 0.0},
			{0.0, 0.0, 0.0, 1.0}
		}};
	}
	static mat SO3(const vec<3, T> &omega) requires (M == 4 && N == 4) {
		if(omega.SqMag() < 1e-10)
			return Identity();
		const T angle = sqrt(omega.SqMag());
		return AxisRotation(omega / angle, angle);
	}
	static mat SE3(const vec<3, T> &angleAxis, const vec<3, T> &translation) requires (M == 4 && N == 4) {
		if(angleAxis.SqMag() < 1e-10){
			 return Translation(translation);
		}
		const T angle = sqrt(angleAxis.SqMag());
		const vec<3, T> u = angleAxis / angle;
		const T c = cos(angle);
		const T s = sin(angle);
		const T t = T(1) - c;
		mat ret = {{
			{t * u.x * u.x + c, t * u.x * u.y + s * u.z, t * u.x * u.z - s * u.y, 0.0},
			{t * u.x * u.y - s * u.z, t * u.y * u.y + c, t * u.y * u.z + s * u.x, 0.0},
			{t * u.x * u.z + s * u.y, t * u.y * u.z - s * u.x, t * u.z * u.z + c, 0.0},
			{0.0, 0.0, 0.0, 1.0}
		}};
		mat Omega = Identity() + (1.0 - c) / (angle * angle) * ret + (angle - s) / (angle * angle * angle) * ret & ret;
		Omega[3].w = 1.0;
		ret[3] = Omega & (translation | 1.0);
		return ret;
	}
};

template <size_t N, size_t M, typename T> std::ostream &operator<<(std::ostream &stream, const mat<N, M, T> &matrix){
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
	T ret = 0;
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
