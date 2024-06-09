#ifndef mattresses_h
#define mattresses_h

#include <iostream>

#include <macros.h>

// macros for counting number of arguments
#define ELEVENTH_ARGUMENT(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, ...) a11
#define COUNT_ARGUMENTS(...) ELEVENTH_ARGUMENT(dummy, ## __VA_ARGS__, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0)



#define DEFINE_VECTOR_UNARY_OPERATOR_OVERLOADS(_operator) \
vec &operator _operator (const vec &other){ \
	for(int i=0; i<N; i++) elements[i] _operator other[i]; \
	return *this; \
} \
vec &operator _operator (T f){ \
	for(int i=0; i<N; i++) elements[i] _operator f; \
	return *this; \
}

template <uint8_t N, typename T=float> struct vec {
	const T &operator[](uint8_t n) const { return elements[n]; }
	T &operator[](uint8_t n){ return elements[n]; }
	
	T elements[N];
	
	T SqMag() const { T ret = 0.0f; for(int i=0; i<N; i++) ret += elements[i] * elements[i]; return ret; }
	vec Normalised() const { return *this / sqrtf(SqMag()); }
	
	DEFINE_VECTOR_UNARY_OPERATOR_OVERLOADS(*=)
	DEFINE_VECTOR_UNARY_OPERATOR_OVERLOADS(+=)
	DEFINE_VECTOR_UNARY_OPERATOR_OVERLOADS(-=)
	DEFINE_VECTOR_UNARY_OPERATOR_OVERLOADS(/=)
	
	template <typename other_t> operator vec<N, other_t>() const { vec<N, other_t> ret; for(int i=0; i<N; ++i) ret[i] = static_cast<other_t>(elements[i]); return ret; }
	
	static vec Zero(){ vec ret; memset(ret.elements, 0, N * sizeof(T)); return ret; }
	static vec One(){ vec ret; for(int i=0; i<N; i++) ret[i] = 1.0; return ret; }
	static vec PositiveCartesianUnit(uint8_t n){ vec ret = vec::Zero(); ret[n] = 1.0; return ret; }
	template <typename T_function_t> static vec FromFunction(T_function_t function){ vec ret; for(int i=0; i<N; ++i) ret[i] = function(); return ret; }
};

template <uint8_t N, typename T> inline std::ostream &operator<<(std::ostream &stream, const vec<N, T> &vector){
	stream << "(";
	for(int n=0; n<N - 1; n++) stream << vector[n] << ", ";
	stream << vector[N - 1] << ")";
	return stream;
}

template <uint8_t N, typename T> inline T Dot(const vec<N, T> &lhs, const vec<N, T> &rhs){
	T ret = 0.0f;
	for(int n=0; n<N; n++) ret += lhs.elements[n] * rhs.elements[n];
	return ret;
}

template <uint8_t N1, uint8_t N2, typename T> inline vec<N1 + N2, T> operator|(const vec<N1, T> &lhs, const vec<N2, T> &rhs){
	vec<N1 + N2, T> ret;
	memcpy(&ret, &lhs, N1*sizeof(T));
	memcpy((T *)(&ret) + N1, &rhs, N2*sizeof(T));
	return ret;
}
template <uint8_t N, typename T> inline vec<N + 1, T> operator|(const vec<N, T> &lhs, T rhs){
	vec<N + 1, T> ret;
	memcpy(&ret, &lhs, N*sizeof(T));
	ret[N] = rhs;
	return ret;
}
template <uint8_t N, typename T> inline vec<N + 1, T> operator|(T lhs, const vec<N, T> &rhs){
	vec<N + 1, T> ret;
	ret[0] = lhs;
	memcpy((T *)(&ret) + 1, &rhs, N*sizeof(T));
	return ret;
}

#define DEFINE_VECTOR_BINARY_OPERATOR_OVERLOADS(_operator) \
template <uint8_t N, typename T> inline vec<N, T> operator _operator (const vec<N, T> &lhs, const vec<N, T> &rhs){ \
	vec<N, T> ret; \
	for(int i=0; i<N; i++) ret[i] = lhs[i] _operator rhs[i]; \
	return ret; \
} \
template <uint8_t N, typename T> inline vec<N, T> operator _operator (const vec<N, T> &lhs, T rhs){ \
	vec<N, T> ret; \
	for(int i=0; i<N; i++) ret[i] = lhs[i] _operator rhs; \
	return ret; \
} \
template <uint8_t N, typename T> inline vec<N, T> operator _operator (T lhs, const vec<N, T> &rhs){ \
	vec<N, T> ret; \
	for(int i=0; i<N; i++) ret[i] = lhs _operator rhs[i]; \
	return ret; \
}
												
DEFINE_VECTOR_BINARY_OPERATOR_OVERLOADS(*)
DEFINE_VECTOR_BINARY_OPERATOR_OVERLOADS(+)
DEFINE_VECTOR_BINARY_OPERATOR_OVERLOADS(-)
DEFINE_VECTOR_BINARY_OPERATOR_OVERLOADS(/)

template <uint8_t N, typename T> inline bool operator==(const vec<N, T> &lhs, const vec<N, T> &rhs){
	bool ret = true;
	for(int n=0; n<N; ++n) ret &= (lhs.elements[n] == rhs.elements[n]);
	return ret;
}

#define _N_VECTOR_ELEMENT(element) T element;
#define _N_VECTOR_OSTREAM(element) vector.element
#define _N_VECTOR_DOT(element) lhs.element * rhs.element
#define _N_VECTOR_EQUAL(element) (lhs.element == rhs.element)
#define _N_VECTOR_ZERO(element) 0.0
#define _N_VECTOR_ONE(element) 1.0
#define _N_VECTOR_SQ(element) element * element
#define _N_VECTOR_FUNCTION(element) function()
#define _N_VECTOR_CAST(element) static_cast<other_t>(element)

#define _N_VECTOR_UNARY_OPERATION(_operator, element) element _operator other.element;
#define _N_VECTOR_FLOAT_UNARY_OPERATION(_operator, element) element _operator f;
#define DEFINE_N_VECTOR_UNARY_OPERATOR_OVERLOADS(_operator, vector_elements...) \
vec &operator _operator (const vec &other){ \
	FOR_EACH_EX(_N_VECTOR_UNARY_OPERATION, _operator, vector_elements) \
	return *this; \
} \
vec &operator _operator (T f){ \
	FOR_EACH_EX(_N_VECTOR_FLOAT_UNARY_OPERATION, _operator, vector_elements) \
	return *this; \
}

#define START_DEFINE_N_VECTOR_STRUCT_SPECIALISATION(uint8_t_size, vector_elements...) \
template <typename T> struct vec<uint8_t_size, T> { \
	const T &operator[](uint8_t n) const { return ((T *)this)[n]; } \
	T &operator[](uint8_t n){ return ((T *)this)[n]; } \
\
	FOR_EACH(_N_VECTOR_ELEMENT, vector_elements) \
\
	T SqMag() const { return CUSTOM_SEPARATED_FOR_EACH(_N_VECTOR_SQ, +, vector_elements); } \
    vec Normalised() const { vec ret = *this; ret /= sqrtf(ret.SqMag()); return ret; } \
\
	DEFINE_N_VECTOR_UNARY_OPERATOR_OVERLOADS(*=, vector_elements) \
	DEFINE_N_VECTOR_UNARY_OPERATOR_OVERLOADS(+=, vector_elements) \
	DEFINE_N_VECTOR_UNARY_OPERATOR_OVERLOADS(-=, vector_elements) \
	DEFINE_N_VECTOR_UNARY_OPERATOR_OVERLOADS(/=, vector_elements) \
\
	template <typename other_t> operator vec< COUNT_ARGUMENTS(vector_elements) , other_t>() const { return {COMMA_SEPARATED_FOR_EACH(_N_VECTOR_CAST, vector_elements)}; } \
\
	static vec Zero(){ return {COMMA_SEPARATED_FOR_EACH(_N_VECTOR_ZERO, vector_elements)}; } \
	static vec One(){ return {COMMA_SEPARATED_FOR_EACH(_N_VECTOR_ONE, vector_elements)}; } \
	static vec PositiveCartesianUnit(uint8_t n){ vec ret = vec::Zero(); ret[n] = 1.0; return ret; } \
	template <typename T_function_t> static vec FromFunction(T_function_t function){ return {COMMA_SEPARATED_FOR_EACH(_N_VECTOR_FUNCTION, vector_elements)}; }

#define FINISH_DEFINE_N_VECTOR_STRUCT_SPECIALISATION(uint8_t_size, vector_elements...) \
};

#define _N_VECTOR_BINARY_OPERATION(_operator, element) lhs.element _operator rhs.element
#define _N_FLOAT_VECTOR_BINARY_OPERATION(_operator, element) lhs _operator rhs.element
#define _N_VECTOR_FLOAT_BINARY_OPERATION(_operator, element) lhs.element _operator rhs
#define _N_VECTOR_BINARY_OPERATOR_OVERLOADS(_operator, uint8_t_size, vector_elements...) \
template <typename T> inline vec<uint8_t_size, T> operator _operator (const vec<uint8_t_size, T> &lhs, const vec<uint8_t_size, T> &rhs){ \
	return { COMMA_SEPARATED_FOR_EACH_EX(_N_VECTOR_BINARY_OPERATION, _operator, vector_elements) }; \
} \
template <typename T> inline vec<uint8_t_size, T> operator _operator (T lhs, const vec<uint8_t_size, T> &rhs){ \
	return { COMMA_SEPARATED_FOR_EACH_EX(_N_FLOAT_VECTOR_BINARY_OPERATION, _operator, vector_elements) }; \
} \
template <typename T> inline vec<uint8_t_size, T> operator _operator (const vec<uint8_t_size, T> &lhs, T rhs){ \
	return { COMMA_SEPARATED_FOR_EACH_EX(_N_VECTOR_FLOAT_BINARY_OPERATION, _operator, vector_elements) }; \
}

#define DEFINE_N_VECTOR_GLOBAL_SPECIALISATIONS(uint8_t_size, vector_elements...) \
template <typename T> inline T Dot(const vec<uint8_t_size, T> &lhs, const vec<uint8_t_size, T> &rhs){ \
   return CUSTOM_SEPARATED_FOR_EACH(_N_VECTOR_DOT, +, vector_elements); \
} \
template <typename T> inline std::ostream &operator<<(std::ostream &stream, const vec<uint8_t_size, T> &vector){ \
	stream << "(" << CUSTOM_SEPARATED_FOR_EACH(_N_VECTOR_OSTREAM, << ", " <<, vector_elements) << ")"; \
	return stream; \
} \
template <typename T> inline bool operator==(const vec<uint8_t_size, T> &lhs, const vec<uint8_t_size, T> &rhs){ \
	return CUSTOM_SEPARATED_FOR_EACH(_N_VECTOR_EQUAL, &&, vector_elements); \
} \
_N_VECTOR_BINARY_OPERATOR_OVERLOADS(*, uint8_t_size, vector_elements) \
_N_VECTOR_BINARY_OPERATOR_OVERLOADS(+, uint8_t_size, vector_elements) \
_N_VECTOR_BINARY_OPERATOR_OVERLOADS(-, uint8_t_size, vector_elements) \
_N_VECTOR_BINARY_OPERATOR_OVERLOADS(/, uint8_t_size, vector_elements)



#define START_DEFINE_N_VECTOR_SPECIALISATION(uint8_t_size, vector_elements...) \
START_DEFINE_N_VECTOR_STRUCT_SPECIALISATION(uint8_t_size, vector_elements)

#define FINISH_DEFINE_N_VECTOR_SPECIALISATION(uint8_t_size, vector_elements...) \
FINISH_DEFINE_N_VECTOR_STRUCT_SPECIALISATION(uint8_t_size, vector_elements) \
DEFINE_N_VECTOR_GLOBAL_SPECIALISATIONS(uint8_t_size, vector_elements)


#define DEFINE_VECTOR_SUBSET_N_REFERENCE_PARENT(N) \
struct subset##N##_reference_parent { \
	subset##N##_reference_parent(vec &_owner) : owner(_owner) {} \
	 \
	virtual void operator=(const vec<N, T> &other) = 0; \
	virtual operator vec<N, T> () const = 0; \
	 \
protected: \
	vec &owner; \
};

#define _VECTOR_DEFINE_SUBSET_2(element1, element2) \
vec<2, T> element1##element2() const { return {element1, element2}; } \
struct element1##element2##_reference : public subset2_reference_parent { \
	element1##element2##_reference(vec &_owner) : subset2_reference_parent(_owner) {} \
	void operator=(const vec<2, T> &other) override { subset2_reference_parent::owner.element1 = other.x; subset2_reference_parent::owner.element2 = other.y; } \
	operator vec<2, T> () const override { return subset2_reference_parent::owner.element1##element2(); } \
}; \
element1##element2##_reference element1##element2##_r(){ return element1##element2##_reference(*this); }

#define _VECTOR_DEFINE_SUBSET_3(element1, element2, element3) \
vec<3, T> element1##element2##element3() const { return {element1, element2, element3}; } \
struct element1##element2##element3##_reference : public subset3_reference_parent { \
	element1##element2##element3##_reference(vec &_owner) : subset3_reference_parent(_owner) {} \
	void operator=(const vec<3, T> &other) override { subset3_reference_parent::owner.element1 = other.x; subset3_reference_parent::owner.element2 = other.y; subset3_reference_parent::owner.element3 = other.z; } \
	operator vec<3, T> () const override { return subset3_reference_parent::owner.element1##element2##element3(); } \
}; \
element1##element2##element3##_reference element1##element2##element3##_r(){ return element1##element2##element3##_reference(*this); }


START_DEFINE_N_VECTOR_SPECIALISATION(1, value)
FINISH_DEFINE_N_VECTOR_SPECIALISATION(1, value)

START_DEFINE_N_VECTOR_SPECIALISATION(2, x, y)
vec Crossed() const {
	return {-y, x};
}
double Angle() const {
	return atan2(y, x);
}

static vec UnitInDirection(double angle){
	return {(T)cos(angle), (T)sin(angle)};
}
static vec RandomUnit(int precision=1000){
	return UnitInDirection(2.0 * M_PI * double(rand() % precision) / double(precision));
}
FINISH_DEFINE_N_VECTOR_SPECIALISATION(2, x, y)

START_DEFINE_N_VECTOR_SPECIALISATION(3, x, y, z)

DEFINE_VECTOR_SUBSET_N_REFERENCE_PARENT(2)

FOR_EACH_PAIR(_VECTOR_DEFINE_SUBSET_2,
			  x, y,
			  y, x,
			  x, z,
			  z, x,
			  y, z,
			  z, y
)

FINISH_DEFINE_N_VECTOR_SPECIALISATION(3, x, y, z)
template<typename T> vec<3, T> Cross(const vec<3, T> &lhs, const vec<3, T> &rhs){
	return {lhs.y * rhs.z - lhs.z * rhs.y, lhs.z * rhs.x - lhs.x * rhs.z, lhs.x * rhs.y - lhs.y * rhs.x};
}

START_DEFINE_N_VECTOR_SPECIALISATION(4, x, y, z, w)

DEFINE_VECTOR_SUBSET_N_REFERENCE_PARENT(2)

FOR_EACH_PAIR(_VECTOR_DEFINE_SUBSET_2,
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

DEFINE_VECTOR_SUBSET_N_REFERENCE_PARENT(3)

FOR_EACH_TRIPLE(_VECTOR_DEFINE_SUBSET_3,
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

FINISH_DEFINE_N_VECTOR_SPECIALISATION(4, x, y, z, w)



// N x M, meaning N rows & M columns
#define _MATRIX_TEMPLATE_PARAMETER(param_name) uint8_t param_name
#define START_MATRIX_TEMPLATE(N, M, ...) \
template < __VA_OPT__( COMMA_SEPARATED_FOR_EACH(_MATRIX_TEMPLATE_PARAMETER, __VA_ARGS__) ) __VA_OPT__(,) typename T=float> struct mat { \
_CONTINUE_MATRIX_TEMPLATE(N, M)

#define START_MATRIX_TEMPLATE_SPECIALISATION(N, M, ...) \
template < __VA_OPT__( COMMA_SEPARATED_FOR_EACH(_MATRIX_TEMPLATE_PARAMETER, __VA_ARGS__) ) __VA_OPT__(,) typename T> struct mat<N, M, T> { \
_CONTINUE_MATRIX_TEMPLATE(N, M)

#define _MATRIX_UNARY_OPERATOR_OVERLOADS(M, _operator) \
mat &operator _operator (const mat &other){ \
	for(int m=0; m<M; m++) columns[m] _operator other[m]; \
	return *this; \
} \
mat &operator _operator (T f){ \
	for(int m=0; m<M; m++) columns[m] _operator f; \
	return *this; \
}

#define _MATRIX_BINARY_OPERATOR_OVERLOADS(M, _operator) \
friend mat operator _operator (const mat &lhs, const mat &rhs){ \
	mat ret; \
	for(int m=0; m<M; m++) ret[m] = lhs[m] _operator rhs[m]; \
	return ret; \
} \
friend mat operator _operator (T lhs, const mat &rhs){ \
	mat ret; \
	for(int m=0; m<M; m++) ret[m] = lhs _operator rhs[m]; \
	return ret; \
} \
friend mat operator _operator (const mat &lhs, T rhs){ \
	mat ret; \
	for(int m=0; m<M; m++) ret[m] = lhs[m] _operator rhs; \
	return ret; \
}

#define _CONTINUE_MATRIX_TEMPLATE(N, M) \
	const vec<N, T> &operator[](uint8_t m) const { return columns[m]; } \
	vec<N, T> &operator[](uint8_t m){ return columns[m]; } \
\
	struct const_row { \
		const_row(const mat &_matrix, uint8_t _n) : matrix(_matrix), n(_n) {} \
\
		const T &operator[](uint8_t m) const { return matrix[m][n]; } \
\
        operator vec<M, T> () const { \
            vec<M, T> ret; \
            for(int m=0; m<M; m++) ret[m] = matrix[m][n];\
            return ret; \
        } \
\
		friend std::ostream &operator<<(std::ostream &stream, const const_row &matrix_row){ \
			stream << "["; \
			for(int m=0; m<M - 1; m++) stream << matrix_row[m] << ", "; \
			stream << matrix_row[M - 1] << "]"; \
			return stream; \
		} \
\
	private: \
		const mat &matrix; \
		const uint8_t n; \
	}; \
	struct row { \
		row(mat &_matrix, uint8_t _n) : matrix(_matrix), n(_n) {} \
\
		row &operator=(const vec<M, T> &col){ \
			for(int i=0; i<M; i++) matrix[i][n] = col[i]; \
			return *this; \
		} \
\
		T &operator[](uint8_t m) const { return matrix[m][n]; } \
\
		operator const_row() const { return const_row(matrix, n); } \
        operator vec<M, T> () const { \
            vec<M, T> ret; \
            for(int m=0; m<M; m++) ret[m] = matrix[m][n];\
            return ret; \
        } \
\
		friend std::ostream &operator<<(std::ostream &stream, const row &matrix_row){ \
			stream << "["; \
			for(int m=0; m<M - 1; m++) stream << matrix_row[m] << ", "; \
			stream << matrix_row[M - 1] << "]"; \
			return stream; \
		} \
\
	private: \
		mat &matrix; \
		const uint8_t n; \
	}; \
\
	const_row operator()(uint8_t n) const { return const_row(*this, n); } \
	row operator()(uint8_t n){ return row(*this, n); } \
\
	_MATRIX_UNARY_OPERATOR_OVERLOADS(M, *=) \
	_MATRIX_UNARY_OPERATOR_OVERLOADS(M, +=) \
	_MATRIX_UNARY_OPERATOR_OVERLOADS(M, -=) \
	_MATRIX_UNARY_OPERATOR_OVERLOADS(M, /=) \
\
	vec<N, T> columns[M]; \
\
	mat<M, N, T> Transposed() const { \
		mat<M, N, T> ret; \
		for(int i=0; i<M; i++) ret(i) = columns[i]; \
		return ret; \
	} \
\
	_MATRIX_BINARY_OPERATOR_OVERLOADS(M, +) \
	_MATRIX_BINARY_OPERATOR_OVERLOADS(M, -) \
	_MATRIX_BINARY_OPERATOR_OVERLOADS(M, *) \
	_MATRIX_BINARY_OPERATOR_OVERLOADS(M, /) \
\
	static mat Zeros(){ mat ret; for(int i=0; i<M; i++) ret[i] = vec<N, T>::Zero(); return ret; }

#define FINISH_MATRIX_TEMPLATE };

// the main template definition
START_MATRIX_TEMPLATE(N, M, N, M)
FINISH_MATRIX_TEMPLATE

// general matrix template functions
template <uint8_t N, uint8_t M, typename T> inline std::ostream &operator<<(std::ostream &stream, const mat<N, M, T> &matrix){
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

template <uint8_t M, typename T> inline std::ostream &operator<<(std::ostream &stream, const mat<1, M, T> &matrix){
	stream << "\n";
	stream << "[ ";
	for(int m=0; m<M - 1; m++) stream << matrix[m][0] << "  ";
	stream << matrix[M - 1][0] << " ]\n";
	return stream;
}

template <uint8_t N, uint8_t M, typename T> inline T Dot(const typename mat<N, M, T>::const_row &lhs, const vec<M, T> &rhs){
	T ret = 0.0f;
	for(int i=0; i<M; i++) ret += lhs[i] * rhs[i];
	return ret;
}

template <uint8_t N1, uint8_t MN, uint8_t M2, typename T> inline mat<N1, M2, T> operator &(const mat<N1, MN, T> &lhs, const mat<MN, M2, T> &rhs){
	mat<N1, M2, T> ret;
	for(int m=0; m<M2; m++) for(int n=0; n<N1; n++) ret[m][n] = Dot<N1, MN, T>(lhs(n), rhs[m]);
	return ret;
}

template <uint8_t N, uint8_t M, typename T> inline vec<N, T> operator &(const mat<N, M, T> &matrix, const vec<M, T> &vector){
	vec<N, T> ret;
	for(int n=0; n<N; n++) ret[n] = Dot<N, M, T>(matrix(n), vector);
	return ret;
}

template <uint8_t N, uint8_t M1, uint8_t M2, typename T> inline mat<N, M1 + M2, T> operator |(const mat<N, M1, T> &lhs, const mat<N, M2, T> &rhs){
    mat<N, M1 + M2, T> ret;
    memcpy(&ret[0][0], &lhs[0][0], sizeof(lhs));
    memcpy(&ret[M1][0], &rhs[0][0], sizeof(rhs));
    return ret;
}
template <uint8_t N, uint8_t M, typename T> inline mat<N, M + 1, T> operator |(const mat<N, M, T> &lhs, const vec<N, T> &rhs){
    mat<N, M + 1, T> ret;
    memcpy(&ret[0][0], &lhs[0][0], sizeof(lhs));
    memcpy(&ret[M][0], &rhs[0], sizeof(rhs));
    return ret;
}
template <uint8_t N, uint8_t M, typename T> inline mat<N, M + 1, T> operator |(const vec<N, T> &lhs, const mat<N, M, T> &rhs){
    mat<N, M + 1, T> ret;
    memcpy(&ret[0][0], &lhs[0], sizeof(lhs));
    memcpy(&ret[1][0], &rhs[0][0], sizeof(rhs));
    return ret;
}

// for square matrices
#define DEFINE_SQUARE_MATRIX_SPECIALISATIONS(N) \
mat &Transpose(){ \
	for(int m=0; m<N; m++) for(int n=0; n<N; n++){ \
		if(m == n) continue; \
		const T temp = columns[m][n]; \
		columns[m][n] = columns[n][m]; \
		columns[n][m] = temp; \
	} \
	return *this; \
} \
static mat Identity(){ \
	mat ret; \
	for(int i=0; i<N; i++) ret[i] = vec<N, T>::PositiveCartesianUnit(i); \
	return ret; \
}

START_MATRIX_TEMPLATE_SPECIALISATION(N, N, N)
DEFINE_SQUARE_MATRIX_SPECIALISATIONS(N)
FINISH_MATRIX_TEMPLATE

// for 4x4 matrices
START_MATRIX_TEMPLATE_SPECIALISATION(4, 4)
DEFINE_SQUARE_MATRIX_SPECIALISATIONS(4)
// member functions
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
FINISH_MATRIX_TEMPLATE

// for single columns matrices
START_MATRIX_TEMPLATE_SPECIALISATION(N, 1, N)
operator vec<N, T> () const {
    return columns[0];
}
FINISH_MATRIX_TEMPLATE

// for single row matrices
START_MATRIX_TEMPLATE_SPECIALISATION(1, M, M)
operator vec<M, T> () const {
    return (*this)(0);
}
FINISH_MATRIX_TEMPLATE

#endif /* mattresses_h */
