#include "../include/mattresses.h"

int main(int argc, char *argv[]){
	
	std::cout << "Hello, World\n";
	
	vec<4> v = vec<4>::Zero();
	std::cout << v << "\n";
	
	v.x += 1.0f;
	std::cout << v << "\n";
	
	v *= 2.0f;
	std::cout << v << "\n";
	
	v -= 1.0f;
	std::cout << v << "\n";
	
	std::cout << v.SqMag() << ", " << v.Normalised() << "\n";
	
	v.ForEach([](float f){ return f * 2.0f; });
	std::cout << v << "\n";
	
//	v.ForEach([](float f) -> float { return 2.0f * f; });
//	std::cout << v << "\n";
	
	vec<5> v5 = vec<5>::GenerateEach([](int i) -> float { return rand() % i; }, 10);
	std::cout << v5 << "\n";
	
	vec<4, uint8_t> iv = v;
	std::cout << static_cast<vec<4, int>>(iv) << "\n";
	
	vec<2> v2 = vec<2>::UnitInDirection(1.5);
	std::cout << v2.Crossed() << "\n";
	std::cout << v2.Angle() << "\n";
	
	vec<3> v3a = {1.0f, 1.1f, 1.2f};
	vec<3> v3b = {-2.0f, -1.1f, 1.2f};
	std::cout << Cross(v3a, v3b) << "\n";
	std::cout << Dot(v3a, v3b) << "\n";
	std::cout << (v3a | v3b | 0.0f) << "\n";
	std::cout << (1.0f | v3a | v) << "\n";
	std::cout << (v3a + v3b) << "\n";
	std::cout << (v3a * v3b) << "\n";
	
	mat<3, 3> m = {vec<3>{1.0f, 2.0f, 3.0f}, vec<3>{4.0f, 5.0f, 6.0f}, vec<3>{7.0f, 8.0f, 9.0f}};
	std::cout << m << "\n";
	std::cout << m.Transposed() << "\n";
	
	m.TransposeInPlace();
	std::cout << m << "\n";
	
	mat<3, 1> mv = {{1.0f, 2.0f, 3.0f}};
	std::cout << vec<3>(mv) << "\n";
	
	mat<1, 3> mvT = mv.Transposed();
	std::cout << vec<3>(mvT) << "\n";
	
	return 0;
}
