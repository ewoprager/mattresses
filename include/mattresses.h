#ifndef mattresses_h
#define mattresses_h

#include <iostream>

#include <macros.h>

#define DEFINE_VECTOR_UNARY_OPERATOR_OVERLOADS(_operator) \
vec &operator _operator (const vec &other){ \
	for(int i=0; i<N; i++) elements[i] _operator other[i]; \
	return *this; \
} \
vec &operator _operator (float f){ \
	for(int i=0; i<N; i++) elements[i] _operator f; \
	return *this; \
}

template <uint8_t N> struct vec {
	const float &operator[](uint8_t n) const { return elements[n]; }
	float &operator[](uint8_t n){ return elements[n]; }
	
	float elements[N];
	
	float SqMag() const { float ret = 0.0f; for(int i=0; i<N; i++) ret += elements[i] * elements[i]; return ret; }
	vec Normalised() const { return *this / sqrtf(SqMag()); }
	
	DEFINE_VECTOR_UNARY_OPERATOR_OVERLOADS(*=)
	DEFINE_VECTOR_UNARY_OPERATOR_OVERLOADS(+=)
	DEFINE_VECTOR_UNARY_OPERATOR_OVERLOADS(-=)
	DEFINE_VECTOR_UNARY_OPERATOR_OVERLOADS(/=)
	
	static vec Zero(){ vec ret; memset(ret.elements, 0, N*sizeof(float)); return ret; }
	static vec One(){ vec ret; for(int i=0; i<N; i++) ret[i] = 1.0f; return ret; }
	static vec PositiveCartesianUnit(uint8_t n){ vec ret = vec::Zero(); ret[n] = 1.0; return ret; }
};

template <uint8_t N> inline std::ostream &operator<<(std::ostream &stream, const vec<N> &vector){
	stream << "(";
	for(int n=0; n<N - 1; n++) stream << vector[n] << ", ";
	stream << vector[N - 1] << ")";
	return stream;
}

template <uint8_t N> inline float Dot(const vec<N> &lhs, const vec<N> &rhs){
	float ret = 0.0f;
	for(int n=0; n<N; n++) ret += lhs.elements[n] * rhs.elements[n];
	return ret;
}

template <uint8_t N1, uint8_t N2> inline vec<N1 + N2> operator|(const vec<N1> &lhs, const vec<N2> &rhs){
	vec<N1 + N2> ret;
	memcpy(&ret, &lhs, N1*sizeof(float));
	memcpy((float *)(&ret) + N1, &rhs, N2*sizeof(float));
	return ret;
}
template <uint8_t N> inline vec<N + 1> operator|(const vec<N> &lhs, float rhs){
	vec<N + 1> ret;
	memcpy(&ret, &lhs, N*sizeof(float));
	ret[N] = rhs;
	return ret;
}
template <uint8_t N> inline vec<N + 1> operator|(float lhs, const vec<N> &rhs){
	vec<N + 1> ret;
	ret[0] = lhs;
	memcpy((float *)(&ret) + 1, &rhs, N*sizeof(float));
	return ret;
}

#define DEFINE_VECTOR_BINARY_OPERATOR_OVERLOADS(_operator) \
template <uint8_t N> inline vec<N> operator _operator (const vec<N> &lhs, const vec<N> &rhs){ \
	vec<N> ret; \
	for(int i=0; i<N; i++) ret[i] = lhs[i] _operator rhs[i]; \
	return ret; \
} \
template <uint8_t N> inline vec<N> operator _operator (const vec<N> &lhs, float rhs){ \
	vec<N> ret; \
	for(int i=0; i<N; i++) ret[i] = lhs[i] _operator rhs; \
	return ret; \
} \
template <uint8_t N> inline vec<N> operator _operator (float lhs, const vec<N> &rhs){ \
	vec<N> ret; \
	for(int i=0; i<N; i++) ret[i] = lhs _operator rhs[i]; \
	return ret; \
}
												
DEFINE_VECTOR_BINARY_OPERATOR_OVERLOADS(*)
DEFINE_VECTOR_BINARY_OPERATOR_OVERLOADS(+)
DEFINE_VECTOR_BINARY_OPERATOR_OVERLOADS(-)
DEFINE_VECTOR_BINARY_OPERATOR_OVERLOADS(/)


#define _N_VECTOR_ELEMENT(element) float element;
#define _N_VECTOR_OSTREAM(element) vector.element
#define _N_VECTOR_DOT(element) lhs.element * rhs.element
#define _N_VECTOR_ZERO(element) 0.0f
#define _N_VECTOR_ONE(element) 1.0f
#define _N_VECTOR_SQ(element) element * element

#define _N_VECTOR_UNARY_OPERATION(_operator, element) element _operator other.element;
#define _N_VECTOR_FLOAT_UNARY_OPERATION(_operator, element) element _operator f;
#define DEFINE_N_VECTOR_UNARY_OPERATOR_OVERLOADS(_operator, uint8_t_size, vector_elements...) \
vec<uint8_t_size> &operator _operator (const vec<uint8_t_size> &other){ \
	FOR_EACH_EX(_N_VECTOR_UNARY_OPERATION, _operator, vector_elements) \
	return *this; \
} \
vec<uint8_t_size> &operator _operator (float f){ \
	FOR_EACH_EX(_N_VECTOR_FLOAT_UNARY_OPERATION, _operator, vector_elements) \
	return *this; \
}

#define START_DEFINE_N_VECTOR_STRUCT_SPECIALISATION(uint8_t_size, vector_elements...) \
template <> struct vec<uint8_t_size> { \
	const float &operator[](uint8_t n) const { return ((float *)this)[n]; } \
	float &operator[](uint8_t n){ return ((float *)this)[n]; } \
\
	FOR_EACH(_N_VECTOR_ELEMENT, vector_elements) \
\
	float SqMag() const { return CUSTOM_SEPARATED_FOR_EACH(_N_VECTOR_SQ, +, vector_elements); } \
    vec<uint8_t_size> Normalised() const { vec<uint8_t_size> ret = *this; ret /= sqrtf(ret.SqMag()); return ret; } \
\
	DEFINE_N_VECTOR_UNARY_OPERATOR_OVERLOADS(*=, uint8_t_size, vector_elements) \
	DEFINE_N_VECTOR_UNARY_OPERATOR_OVERLOADS(+=, uint8_t_size, vector_elements) \
	DEFINE_N_VECTOR_UNARY_OPERATOR_OVERLOADS(-=, uint8_t_size, vector_elements) \
	DEFINE_N_VECTOR_UNARY_OPERATOR_OVERLOADS(/=, uint8_t_size, vector_elements) \
\
	static vec Zero(){ return {COMMA_SEPARATED_FOR_EACH(_N_VECTOR_ZERO, vector_elements)}; } \
	static vec One(){ return {COMMA_SEPARATED_FOR_EACH(_N_VECTOR_ONE, vector_elements)}; } \
	static vec PositiveCartesianUnit(uint8_t n){ vec ret = vec::Zero(); ret[n] = 1.0f; return ret; }

#define FINISH_DEFINE_N_VECTOR_STRUCT_SPECIALISATION(uint8_t_size, vector_elements...) \
};

#define _N_VECTOR_BINARY_OPERATION(_operator, element) lhs.element _operator rhs.element
#define _N_FLOAT_VECTOR_BINARY_OPERATION(_operator, element) lhs _operator rhs.element
#define _N_VECTOR_FLOAT_BINARY_OPERATION(_operator, element) lhs.element _operator rhs
#define _N_VECTOR_BINARY_OPERATOR_OVERLOADS(_operator, uint8_t_size, vector_elements...) \
template <> inline vec<uint8_t_size> operator _operator (const vec<uint8_t_size> &lhs, const vec<uint8_t_size> &rhs){ \
	return { COMMA_SEPARATED_FOR_EACH_EX(_N_VECTOR_BINARY_OPERATION, _operator, vector_elements) }; \
} \
template <> inline vec<uint8_t_size> operator _operator (float lhs, const vec<uint8_t_size> &rhs){ \
	return { COMMA_SEPARATED_FOR_EACH_EX(_N_FLOAT_VECTOR_BINARY_OPERATION, _operator, vector_elements) }; \
} \
template <> inline vec<uint8_t_size> operator _operator (const vec<uint8_t_size> &lhs, float rhs){ \
	return { COMMA_SEPARATED_FOR_EACH_EX(_N_VECTOR_FLOAT_BINARY_OPERATION, _operator, vector_elements) }; \
}

#define DEFINE_N_VECTOR_GLOBAL_SPECIALISATIONS(uint8_t_size, vector_elements...) \
template <> inline float Dot(const vec<uint8_t_size> &lhs, const vec<uint8_t_size> &rhs){ \
   return CUSTOM_SEPARATED_FOR_EACH(_N_VECTOR_DOT, +, vector_elements); \
} \
template <> inline std::ostream &operator<<(std::ostream &stream, const vec<uint8_t_size> &vector){ \
	stream << "(" << CUSTOM_SEPARATED_FOR_EACH(_N_VECTOR_OSTREAM, << ", " <<, vector_elements) << ")"; \
	return stream; \
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
	virtual void operator=(const vec<N> &other) = 0; \
	virtual operator vec<N> () const = 0; \
	 \
protected: \
	vec &owner; \
};

#define _VECTOR_DEFINE_SUBSET_2(element1, element2) \
vec<2> element1##element2() const { return {element1, element2}; } \
struct element1##element2##_reference : public subset2_reference_parent { \
	element1##element2##_reference(vec &owner) : subset2_reference_parent(owner) {} \
	void operator=(const vec<2> &other) override { owner.element1 = other.x; owner.element2 = other.y; } \
operator vec<2> () const override { return owner.element1##element2(); } \
}; \
element1##element2##_reference element1##element2##_r(){ return element1##element2##_reference(*this); }

#define _VECTOR_DEFINE_SUBSET_3(element1, element2, element3) \
vec<3> element1##element2##element3() const { return {element1, element2, element3}; } \
struct element1##element2##element3##_reference : public subset3_reference_parent { \
	element1##element2##element3##_reference(vec &owner) : subset3_reference_parent(owner) {} \
	void operator=(const vec<3> &other) override { owner.element1 = other.x; owner.element2 = other.y; owner.element3 = other.z; } \
	operator vec<3> () const override { return owner.element1##element2##element3(); } \
}; \
element1##element2##element3##_reference element1##element2##element3##_r(){ return element1##element2##element3##_reference(*this); }


START_DEFINE_N_VECTOR_SPECIALISATION(1, value)
FINISH_DEFINE_N_VECTOR_SPECIALISATION(1, value)

START_DEFINE_N_VECTOR_SPECIALISATION(2, x, y)
vec Crossed() const {
	return {-y, x};
}
float Angle() const {
	return atan2f(y, x);
}

static vec UnitInDirection(float angle){
	return {cosf(angle), sinf(angle)};
}
static vec RandomUnit(int precision=1000){
	return UnitInDirection(2.0f * (float)M_PI * (rand() % precision) / (float)precision);
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
inline vec<3> Cross(const vec<3> &lhs, const vec<3> &rhs){
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
template < __VA_OPT__( COMMA_SEPARATED_FOR_EACH(_MATRIX_TEMPLATE_PARAMETER, __VA_ARGS__) ) > struct mat { \
_CONTINUE_MATRIX_TEMPLATE(N, M)

#define START_MATRIX_TEMPLATE_SPECIALISATION(N, M, ...) \
template < __VA_OPT__( COMMA_SEPARATED_FOR_EACH(_MATRIX_TEMPLATE_PARAMETER, __VA_ARGS__) ) > struct mat<N, M> { \
_CONTINUE_MATRIX_TEMPLATE(N, M)

#define _MATRIX_UNARY_OPERATOR_OVERLOADS(M, _operator) \
mat &operator _operator (const mat &other){ \
	for(int m=0; m<M; m++) columns[m] _operator other[m]; \
	return *this; \
} \
mat &operator _operator (float f){ \
	for(int m=0; m<M; m++) columns[m] _operator f; \
	return *this; \
}

#define _MATRIX_BINARY_OPERATOR_OVERLOADS(M, _operator) \
friend mat operator _operator (const mat &lhs, const mat &rhs){ \
	mat ret; \
	for(int m=0; m<M; m++) ret[m] = lhs[m] _operator rhs[m]; \
	return ret; \
} \
friend mat operator _operator (float lhs, const mat &rhs){ \
	mat ret; \
	for(int m=0; m<M; m++) ret[m] = lhs _operator rhs[m]; \
	return ret; \
} \
friend mat operator _operator (const mat &lhs, float rhs){ \
	mat ret; \
	for(int m=0; m<M; m++) ret[m] = lhs[m] _operator rhs; \
	return ret; \
}

#define _CONTINUE_MATRIX_TEMPLATE(N, M) \
	const vec<N> &operator[](uint8_t m) const { return columns[m]; } \
	vec<N> &operator[](uint8_t m){ return columns[m]; } \
\
	struct const_row { \
		const_row(const mat &_matrix, uint8_t _n) : matrix(_matrix), n(_n) {} \
\
		const float &operator[](uint8_t m) const { return matrix[m][n]; } \
\
        operator vec<M> () const { \
            vec<M> ret; \
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
		row &operator=(const vec<M> &col){ \
			for(int i=0; i<M; i++) matrix[i][n] = col[i]; \
			return *this; \
		} \
\
		float &operator[](uint8_t m) const { return matrix[m][n]; } \
\
		operator const_row() const { return const_row(matrix, n); } \
        operator vec<M> () const { \
            vec<M> ret; \
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
	vec<N> columns[M]; \
\
	mat<M, N> Transposed() const { \
		mat<M, N> ret; \
		for(int i=0; i<M; i++) ret(i) = columns[i]; \
		return ret; \
	} \
\
	_MATRIX_BINARY_OPERATOR_OVERLOADS(M, +) \
	_MATRIX_BINARY_OPERATOR_OVERLOADS(M, -) \
	_MATRIX_BINARY_OPERATOR_OVERLOADS(M, *) \
	_MATRIX_BINARY_OPERATOR_OVERLOADS(M, /) \
\
	static mat Zeros(){ mat ret; for(int i=0; i<M; i++) ret[i] = vec<N>::Zero(); return ret; }

#define FINISH_MATRIX_TEMPLATE };

// the main template definition
START_MATRIX_TEMPLATE(N, M, N, M)
FINISH_MATRIX_TEMPLATE

// general matrix template functions
template <uint8_t N, uint8_t M> inline std::ostream &operator<<(std::ostream &stream, const mat<N, M> &matrix){
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

template <uint8_t M> inline std::ostream &operator<<(std::ostream &stream, const mat<1, M> &matrix){
	stream << "\n";
	stream << "[ ";
	for(int m=0; m<M - 1; m++) stream << matrix[m][0] << "  ";
	stream << matrix[M - 1][0] << " ]\n";
	return stream;
}

template <uint8_t N, uint8_t M> inline float Dot(const typename mat<N, M>::const_row &lhs, const vec<M> &rhs){
	float ret = 0.0f;
	for(int i=0; i<M; i++) ret += lhs[i] * rhs[i];
	return ret;
}

template <uint8_t N1, uint8_t MN, uint8_t M2> inline mat<N1, M2> Dot(const mat<N1, MN> &lhs, const mat<MN, M2> &rhs){
	mat<N1, M2> ret;
	for(int m=0; m<M2; m++) for(int n=0; n<N1; n++) ret[m][n] = Dot<N1, MN>(lhs(n), rhs[m]);
	return ret;
}

template <uint8_t N, uint8_t M> inline vec<N> Dot(const mat<N, M> &matrix, const vec<M> &vector){
	vec<N> ret;
	for(int n=0; n<N; n++) ret[n] = Dot<N, M>(matrix(n), vector);
	return ret;
}

template <uint8_t N, uint8_t M1, uint8_t M2> inline mat<N, M1 + M2> operator |(const mat<N, M1> &lhs, const mat<N, M2> &rhs){
    mat<N, M1 + M2> ret;
    memcpy(&ret[0][0], &lhs[0][0], sizeof(lhs));
    memcpy(&ret[M1][0], &rhs[0][0], sizeof(rhs));
    return ret;
}
template <uint8_t N, uint8_t M> inline mat<N, M + 1> operator |(const mat<N, M> &lhs, const vec<N> &rhs){
    mat<N, M + 1> ret;
    memcpy(&ret[0][0], &lhs[0][0], sizeof(lhs));
    memcpy(&ret[M][0], &rhs[0], sizeof(rhs));
    return ret;
}
template <uint8_t N, uint8_t M> inline mat<N, M + 1> operator |(const vec<N> &lhs, const mat<N, M> &rhs){
    mat<N, M + 1> ret;
    memcpy(&ret[0][0], &lhs[0], sizeof(lhs));
    memcpy(&ret[1][0], &rhs[0][0], sizeof(rhs));
    return ret;
}

// for square matrices
#define DEFINE_SQUARE_MATRIX_SPECIALISATIONS(N) \
mat &Transpose(){ \
	for(int m=0; m<N; m++) for(int n=0; n<N; n++){ \
		if(m == n) continue; \
		const float temp = columns[m][n]; \
		columns[m][n] = columns[n][m]; \
		columns[n][m] = temp; \
	} \
	return *this; \
} \
static mat Identity(){ \
	mat ret; \
	for(int i=0; i<N; i++) ret[i] = vec<N>::PositiveCartesianUnit(i); \
	return ret; \
}

START_MATRIX_TEMPLATE_SPECIALISATION(N, N, N)
DEFINE_SQUARE_MATRIX_SPECIALISATIONS(N)
FINISH_MATRIX_TEMPLATE

// for 4x4 matrices
START_MATRIX_TEMPLATE_SPECIALISATION(4, 4)
DEFINE_SQUARE_MATRIX_SPECIALISATIONS(4)
// member functions
mat ScaledExcludingTranslation(const vec<3> &scaling){
	const vec<3> savedTranslation = columns[3].xyz();
	mat ret = Dot(Scaling(scaling), *this);
	ret[3].xyz_r() = savedTranslation;
	return ret;
}
mat Inverted(){
	const float m00 = columns[0][0];
	const float m01 = columns[0][1];
	const float m02 = columns[0][2];
	const float m03 = columns[0][3];
	const float m10 = columns[1][0];
	const float m11 = columns[1][1];
	const float m12 = columns[1][2];
	const float m13 = columns[1][3];
	const float m20 = columns[2][0];
	const float m21 = columns[2][1];
	const float m22 = columns[2][2];
	const float m23 = columns[2][3];
	const float m30 = columns[3][0];
	const float m31 = columns[3][1];
	const float m32 = columns[3][2];
	const float m33 = columns[3][3];
	const float tmp_0  = m22 * m33;
	const float tmp_1  = m32 * m23;
	const float tmp_2  = m12 * m33;
	const float tmp_3  = m32 * m13;
	const float tmp_4  = m12 * m23;
	const float tmp_5  = m22 * m13;
	const float tmp_6  = m02 * m33;
	const float tmp_7  = m32 * m03;
	const float tmp_8  = m02 * m23;
	const float tmp_9  = m22 * m03;
	const float tmp_10 = m02 * m13;
	const float tmp_11 = m12 * m03;
	const float tmp_12 = m20 * m31;
	const float tmp_13 = m30 * m21;
	const float tmp_14 = m10 * m31;
	const float tmp_15 = m30 * m11;
	const float tmp_16 = m10 * m21;
	const float tmp_17 = m20 * m11;
	const float tmp_18 = m00 * m31;
	const float tmp_19 = m30 * m01;
	const float tmp_20 = m00 * m21;
	const float tmp_21 = m20 * m01;
	const float tmp_22 = m00 * m11;
	const float tmp_23 = m10 * m01;

	const float t0 = (tmp_0 * m11 + tmp_3 * m21 + tmp_4 * m31) - (tmp_1 * m11 + tmp_2 * m21 + tmp_5 * m31);
	const float t1 = (tmp_1 * m01 + tmp_6 * m21 + tmp_9 * m31) - (tmp_0 * m01 + tmp_7 * m21 + tmp_8 * m31);
	const float t2 = (tmp_2 * m01 + tmp_7 * m11 + tmp_10 * m31) - (tmp_3 * m01 + tmp_6 * m11 + tmp_11 * m31);
	const float t3 = (tmp_5 * m01 + tmp_8 * m11 + tmp_11 * m21) - (tmp_4 * m01 + tmp_9 * m11 + tmp_10 * m21);

	const float d = 1.0f / (m00 * t0 + m10 * t1 + m20 * t2 + m30 * t3);

	return d * (mat<4, 4>){{
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
static mat XRotation(float angle){
	const float c = cosf(angle);
	const float s = sinf(angle);
	return {{
		{1.0f, 0.0f, 0.0f, 0.0f},
		{0.0f,    c,    s, 0.0f},
		{0.0f,   -s,    c, 0.0f},
		{0.0f, 0.0f, 0.0f, 1.0f}
	}};
}
static mat YRotation(float angle){
	const float c = cosf(angle);
	const float s = sinf(angle);
	return {{
		{   c, 0.0f,   -s, 0.0f},
		{0.0f, 1.0f, 0.0f, 0.0f},
		{   s, 0.0f,    c, 0.0f},
		{0.0f, 0.0f, 0.0f, 1.0f}
	}};
}
static mat ZRotation(float angle){
	const float c = cosf(angle);
	const float s = sinf(angle);
	return {{
		{   c,    s, 0.0f, 0.0f},
		{  -s,    c, 0.0f, 0.0f},
		{0.0f, 0.0f, 1.0f, 0.0f},
		{0.0f, 0.0f, 0.0f, 1.0f}
	}};
}
static mat PerspectiveProjection(const float angleOfView, const float aspectRatio, const float frustumNear, const float frustumFar){
	const float f = tanf(0.5f * ((float)M_PI - angleOfView));
	const float frustumDepthInverse = 1.0f / (frustumNear - frustumFar);
	return {{
		{f / aspectRatio, 0.0f,                                                  0.0f,  0.0f},
		{           0.0f,    f,                                                  0.0f,  0.0f},
		{           0.0f, 0.0f,      (frustumNear + frustumFar) * frustumDepthInverse, -1.0f},
		{           0.0f, 0.0f, frustumNear * frustumFar * frustumDepthInverse * 2.0f,  0.0f}
	}};
}
static mat Translation(const vec<3> &vector){
	return {{
		vec<4>::PositiveCartesianUnit(0),
		vec<4>::PositiveCartesianUnit(1),
		vec<4>::PositiveCartesianUnit(2),
		vector | 1.0f
	}};
}
static mat LookAt(const vec<3> &position, const vec<3> &target, const vec<3> &unitUp){
	const vec<3> zAxis = (position - target).Normalised();
	const vec<3> xAxis = Cross(unitUp, zAxis).Normalised();
	const vec<3> yAxis = Cross(zAxis, xAxis).Normalised();
	return {{
		xAxis | 0.0f,
		yAxis | 0.0f,
		zAxis | 0.0f,
		position | 1.0f
	}};
}
static mat Scaling(const vec<3> &scaling){
	return {{
		{scaling.x, 0.0f, 0.0f, 0.0f},
		{0.0f, scaling.y, 0.0f, 0.0f},
		{0.0f, 0.0f, scaling.z, 0.0f},
		{0.0f, 0.0f, 0.0f, 1.0f}
	}};
}
FINISH_MATRIX_TEMPLATE

// for single columns matrices
START_MATRIX_TEMPLATE_SPECIALISATION(N, 1, N)
operator vec<N> () const {
    return columns[0];
}
FINISH_MATRIX_TEMPLATE

// for single row matrices
START_MATRIX_TEMPLATE_SPECIALISATION(1, M, M)
operator vec<M> () const {
    return (*this)(0);
}
FINISH_MATRIX_TEMPLATE

#endif /* mattresses_h */
