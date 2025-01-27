//#define ACCELERATE_NEW_LAPACK

#include <Accelerate/Accelerate.h>

#ifndef CBLAS_INDEX
#	define CBLAS_INT __LAPACK_int
#	define CBLAS_INDEX size_t
#	define CBLAS_S __LAPACK_float_complex
#	define CBLAS_Z __LAPACK_double_complex
#else // old Accelerate version
#	define CBLAS_INT int
#	define CBLAS_S void
#	define CBLAS_Z void
#endif
extern "C" {
	void cblas_scopy(const CBLAS_INT N, const float *X, const CBLAS_INT incX, float *Y, const CBLAS_INT incY);
	void cblas_dcopy(const CBLAS_INT N, const double *X, const CBLAS_INT incX, double *Y, const CBLAS_INT incY);
	void cblas_ccopy(const CBLAS_INT N, const CBLAS_S *X, const CBLAS_INT incX, CBLAS_S *Y, const CBLAS_INT incY);
	void cblas_zcopy(const CBLAS_INT N, const CBLAS_Z *X, const CBLAS_INT incX, CBLAS_Z *Y, const CBLAS_INT incY);
};

namespace signalsmith { namespace linear {

template<class T>
std::string typeName() {
	return typeid(T).name();
}
template<>
std::string typeName<float>() {
	return "float";
}
template<>
std::string typeName<double>() {
	return "double";
}
template<>
std::string typeName<std::complex<float>>() {
	return "complex<float>";
}
template<>
std::string typeName<std::complex<double>>() {
	return "complex<double>";
}
static size_t _basicFillWarningCounter = 1;
static void basicFillWarningReset() {
	// All warnings will print again
	++_basicFillWarningCounter;
}
template<class Expr>
static void basicFillWarning() {
	static size_t counter = 0;
	if (counter != _basicFillWarningCounter) {
		counter = _basicFillWarningCounter;
		std::cerr << "used basic .fill() for " << Expr::name() << " -> " << typeName<decltype(std::declval<Expr>().get(0))>() << "\n";
	}
}


template<>
struct LinearImpl<true> : public LinearImplBase<true> {
	using Base = LinearImplBase<true>;

	LinearImpl() : Base(this), cached(*this) {
		basicFillWarningReset();
	}

	template<class V>
	void reserve(size_t) {}
	
	template<>
	void reserve<float>(size_t size) {
		cached.reserveFloats(size*4);
	}
	template<>
	void reserve<double>(size_t size) {
		cached.reserveDoubles(size*4);
	}

	template<class Pointer, class Expr>
	void fill(Pointer pointer, Expr expr, size_t size) {
		fillBasic(pointer, expr, size);
	}

	template<class Pointer, class Expr>
	void fill(Pointer pointer, Expression<Expr> expr, size_t size) {
		return self().fill(pointer, (Expr &)expr, size);
	};

private:
	template<class Pointer, class Expr>
	void fillBasic(Pointer pointer, Expr expr, size_t size) {
		basicFillWarning<Expr>();
		for (size_t i = 0; i < size; ++i) {
			pointer[i] = expr.get(i);
		}
	}

// Forwards .fill() to .fillName(), but doesn't define that
#define SIGNALSMITH_AUDIO_LINEAR_OP2_R(Name) \
public: \
	template<class A, class B> \
	void fill(RealPointer<float> pointer, expression::Name<A, B> expr, size_t size) { \
		fill##Name(pointer, expr, size); \
	} \
	template<class A, class B> \
	void fill(RealPointer<double> pointer, expression::Name<A, B> expr, size_t size) { \
		fill##Name(pointer, expr, size); \
	} \
// Real -> real operators where the vDSP function arguments are the other way around for some reason
#define SIGNALSMITH_AUDIO_LINEAR_TREE2FLIP_RR(Name, vDSP_func) \
	template<class A, class B> \
	void fill##Name(RealPointer<float> pointer, expression::Name<A, B> expr, size_t size) { \
		auto floats = cached.floatScope(); \
		auto *a = floats.real(expr.a, size); \
		auto *b = floats.real(expr.b, size); \
		vDSP_func(b, 1, a, 1, pointer, 1, size); \
	} \
	template<class A, class B> \
	void fill##Name(RealPointer<double> pointer, expression::Name<A, B> expr, size_t size) { \
		auto doubles = cached.doubleScope(); \
		auto *a = doubles.real(expr.a, size); \
		auto *b = doubles.real(expr.b, size); \
		vDSP_func##D(b, 1, a, 1, pointer, 1, size); \
	}
	SIGNALSMITH_AUDIO_LINEAR_OP2_R(Add)
	SIGNALSMITH_AUDIO_LINEAR_TREE2FLIP_RR(Add, vDSP_vadd);
	SIGNALSMITH_AUDIO_LINEAR_OP2_R(Sub)
	SIGNALSMITH_AUDIO_LINEAR_TREE2FLIP_RR(Sub, vDSP_vsub);
	SIGNALSMITH_AUDIO_LINEAR_OP2_R(Mul)
	SIGNALSMITH_AUDIO_LINEAR_TREE2FLIP_RR(Mul, vDSP_vmul);
	SIGNALSMITH_AUDIO_LINEAR_OP2_R(Div)
	SIGNALSMITH_AUDIO_LINEAR_TREE2FLIP_RR(Div, vDSP_vdiv);
#undef SIGNALSMITH_AUDIO_LINEAR_TREE2FLIP_RR
#undef SIGNALSMITH_AUDIO_LINEAR_OP2_R

// Forwards .fillName() to .fillNameName2(), but doesn't define that
#define SIGNALSMITH_AUDIO_LINEAR_Op3L_R(Name, NameL) \
	template<class A, class B, class C> \
	void fill##Name(RealPointer<float> pointer, expression::Name<expression::NameL<A, B>, C> expr, size_t size) { \
		fill##Name##NameL(pointer, expr, size); \
	} \
	template<class A, class B, class C> \
	void fill##Name(RealPointer<double> pointer, expression::Name<expression::NameL<A, B>, C> expr, size_t size) { \
		fill##Name##NameL(pointer, expr, size); \
	}
#define SIGNALSMITH_AUDIO_LINEAR_TREE3L_RRR(Name, NameL, vDSP_func) \
	template<class A, class B, class C> \
	void fill##Name##NameL(RealPointer<float> pointer, expression::Name<expression::NameL<A, B>, C> expr, size_t size) { \
		auto floats = cached.floatScope(); \
		auto *a = floats.real(expr.a.a, size); \
		auto *b = floats.real(expr.a.b, size); \
		auto *c = floats.real(expr.b, size); \
		vDSP_func(a, 1, b, 1, c, 1, pointer, 1, size); \
	} \
	template<class A, class B, class C> \
	void fill##Name##NameL(RealPointer<double> pointer, expression::Name<expression::NameL<A, B>, C> expr, size_t size) { \
		auto doubles = cached.doubleScope(); \
		auto *a = doubles.real(expr.a.a, size); \
		auto *b = doubles.real(expr.a.b, size); \
		auto *c = doubles.real(expr.b, size); \
		vDSP_func##D(a, 1, b, 1, c, 1, pointer, 1, size); \
	}
	SIGNALSMITH_AUDIO_LINEAR_Op3L_R(Mul, Add)
	SIGNALSMITH_AUDIO_LINEAR_TREE3L_RRR(Mul, Add, vDSP_vam)

#undef SIGNALSMITH_AUDIO_LINEAR_TREE3L_RRR
#undef SIGNALSMITH_AUDIO_LINEAR_Op3L_R

//	template<class A, class B, class C> \
//	void fill#Name(RealPointer<float> pointer, expression::Name<expression::NameInner<A, B>, C> expr, size_t size, const priority<4*priorityN> &) { \
//		auto floats = cached.floatScope(); \
//		auto *a = floats.real(expr.a.a, size); \
//		auto *b = floats.real(expr.a.b, size); \
//		auto *c = floats.real(expr.b, size); \
//		vDSP_func(a, 1, b, 1, c, 1, pointer, 1, size); \
//	} \
//	template<class A, class B, class C> \
//	void fastFill(RealPointer<double> pointer, expression::Name<expression::NameInner<A, B>, C> expr, size_t size, const priority<4*priorityN> &) { \
//		auto doubles = cached.doubleScope(); \
//		auto *a = doubles.real(expr.a.a, size); \
//		auto *b = doubles.real(expr.a.b, size); \
//		auto *c = doubles.real(expr.b, size); \
//		vDSP_func##D(a, 1, b, 1, c, 1, pointer, 1, size); \
//	}

//#define SIGNALSMITH_AUDIO_LINEAR_TREE3COMMUTATIVE_RRR(NameInner, Name, vDSP_func, priorityN) \
//	SIGNALSMITH_AUDIO_LINEAR_TREE3L_RRR(NameInner, Name, vDSP_func, priorityN) \
//	template<class A, class B, class C> \
//	void fastFill(RealPointer<float> pointer, expression::Name<C, expression::NameInner<A, B>> expr, size_t size, const priority<4*priorityN + 1> &) { \
//		auto floats = cached.floatScope(); \
//		auto *a = floats.real(expr.b.a, size); \
//		auto *b = floats.real(expr.b.b, size); \
//		auto *c = floats.real(expr.a, size); \
//		vDSP_func(a, 1, b, 1, c, 1, pointer, 1, size); \
//	} \
//	template<class A, class B, class C> \
//	void fastFill(RealPointer<double> pointer, expression::Name<C, expression::NameInner<A, B>> expr, size_t size, const priority<4*priorityN + 1> &) { \
//		auto doubles = cached.doubleScope(); \
//		auto *a = doubles.real(expr.b.a, size); \
//		auto *b = doubles.real(expr.b.b, size); \
//		auto *c = doubles.real(expr.a, size); \
//		vDSP_func##D(a, 1, b, 1, c, 1, pointer, 1, size); \
//	}
//	SIGNALSMITH_AUDIO_LINEAR_TREE3COMMUTATIVE_RRR(Add, Mul, vDSP_vam, 0)
//	SIGNALSMITH_AUDIO_LINEAR_TREE3COMMUTATIVE_RRR(Mul, Add, vDSP_vma, 1)
//	SIGNALSMITH_AUDIO_LINEAR_TREE3COMMUTATIVE_RRR(Sub, Mul, vDSP_vsbm, 2)
//	SIGNALSMITH_AUDIO_LINEAR_TREE3L_RRR(Mul, Sub, vDSP_vmsb, 3)
//#undef SIGNALSMITH_AUDIO_LINEAR_TREE3L_RRR

protected:
	CachedResults<LinearImpl, 32> cached;
/*
	using Base::copy;
	void copy(const int N, const float *x, const int xStride, float *y, const int yStride) {
		cblas_scopy(CBLAS_INT(N), x, CBLAS_INT(xStride), y, CBLAS_INT(yStride));
	}
	void copy(const int N, const std::complex<float> *x, const int xStride, std::complex<float> *y, const int yStride) {
		cblas_ccopy(CBLAS_INT(N), x, CBLAS_INT(xStride), y, CBLAS_INT(yStride));
	}

	using Base::norm2;
	float norm2(const int N, const float *x, const int xStride) {
		float sum2 = 0;
		vDSP_svesq(x, std::abs(xStride), &sum2, N);
		return sum2;
	}
	float norm2(const int N, const std::complex<float> *x, const int xStride) {
		if (std::abs(xStride) == 1) {
			float sum2 = 0;
			// Aliasing crimes, but this specific one is as kosher as it gets
			vDSP_svesq((const float *)x, std::abs(xStride), &sum2, N*2);
			return sum2;
		} else {
			// criiiimes
			float real2 = 0, imag2 = 0;
			vDSP_svesq((const float *)x, std::abs(xStride)*2, &real2, N);
			vDSP_svesq((const float *)x + 1, std::abs(xStride)*2, &imag2, N);
			return real2 + imag2;
		}
	}

	float norm2(const int N, ConstSplitComplex<float> x, const int xStride) {
		float real2 = 0, imag2 = 0;
		vDSP_svesq(x.real, std::abs(xStride), &real2, N);
		vDSP_svesq(x.imag, std::abs(xStride), &imag2, N);
		return real2 + imag2;
	}

	void norm2(const int N, const std::complex<float> *x, const int xStride, float *y, const int yStride) {
		if (xStride < 0) x += (1 - N)*xStride;
		if (yStride < 0) y += (1 - N)*yStride;
		
		auto splitX = dspSplit(x);
		vDSP_zvmags(&splitX, 2*xStride, y, yStride, N);
	}

	void norm2(const int N, ConstSplitComplex<float> x, const int xStride, float *y, const int yStride) {
		if (xStride < 0) x += (1 - N)*xStride;
		if (yStride < 0) y += (1 - N)*yStride;

		auto splitX = dspSplit(x);
		vDSP_zvmags(&splitX, xStride, y, yStride, N);
	}

	using Base::mul;
	void mul(const int N, const float *a, const int aStride, const float *b, const int bStride, float *c, const int cStride) {
		if (aStride < 0) a += (1 - N)*aStride;
		if (bStride < 0) b += (1 - N)*bStride;
		if (cStride < 0) c += (1 - N)*cStride;

		vDSP_vmul(a, aStride, b, bStride, c, cStride, N);
	}
	void mul(const int N, const float *a, const int aStride, float *b, const int bStride) {
		if (aStride < 0) a += (1 - N)*aStride;
		if (bStride < 0) b += (1 - N)*bStride;

		vDSP_vmul(a, aStride, b, bStride, b, bStride, N);
	}

	void mul(const int N, const std::complex<float> *a, const int aStride, const std::complex<float> *b, const int bStride, std::complex<float> *c, const int cStride) {
		if (aStride < 0) a += (1 - N)*aStride;
		if (bStride < 0) b += (1 - N)*bStride;
		if (cStride < 0) c += (1 - N)*cStride;

		auto splitA = dspSplit(a), splitB = dspSplit(b), splitC = dspSplit(c);
		vDSP_zvmul(&splitA, aStride*2, &splitB, bStride*2, &splitC, cStride*2, N, 1);
	}

	void mul(const int N, ConstSplitComplex<float> a, const int aStride, ConstSplitComplex<float> b, const int bStride, SplitComplex<float> c, const int cStride) {
		if (aStride < 0) a += (1 - N)*aStride;
		if (bStride < 0) b += (1 - N)*bStride;
		if (cStride < 0) c += (1 - N)*cStride;

		auto splitA = dspSplit(a), splitB = dspSplit(b), splitC = dspSplit(c);
		vDSP_zvmul(&splitA, aStride, &splitB, bStride, &splitC, cStride, N, 1);
	}

	template<class V>
	void reserve(size_t) {}
	
	template<>
	void reserve<float>(size_t) {}
	template<>
	void reserve<double>(size_t) {}
*/
	static DSPSplitComplex dspSplit(ConstSplitPointer<float> x) {
		DSPSplitComplex dsp;
		dsp.realp = (float *)x.real;
		dsp.imagp = (float *)x.imag;
		return dsp;
	}
	static DSPSplitComplex dspSplit(ConstComplexPointer<float> x) {
		DSPSplitComplex dsp;
		dsp.realp = (float *)x;
		dsp.imagp = (float *)x + 1;
		return dsp;
	}
	static DSPDoubleSplitComplex dspSplit(ConstSplitPointer<double> x) {
		DSPDoubleSplitComplex dsp;
		dsp.realp = (double *)x.real;
		dsp.imagp = (double *)x.imag;
		return dsp;
	}
	static DSPDoubleSplitComplex dspSplit(const std::complex<double> *x) {
		DSPDoubleSplitComplex dsp;
		dsp.realp = (double *)x;
		dsp.imagp = (double *)x + 1;
		return dsp;
	}
};

/*
template<>
struct Linear<double> : public LinearImplBase<double> {
	using Base = LinearImplBase<double>;

	Linear() : Base(this) {}

	using Base::copy;
	void copy(const int N, const double *x, const int xStride, double *y, const int yStride) {
		cblas_dcopy(CBLAS_INT(N), x, CBLAS_INT(xStride), y, CBLAS_INT(yStride));
	}
	void copy(const int N, const std::complex<double> *x, const int xStride, std::complex<double> *y, const int yStride) {
		cblas_zcopy(CBLAS_INT(N), x, CBLAS_INT(xStride), y, CBLAS_INT(yStride));
	}

	using Base::norm2;
	double norm2(const int N, const double *x, const int xStride) {
		double sum2 = 0;
		vDSP_svesqD(x, std::abs(xStride), &sum2, N);
		return sum2;
	}
	double norm2(const int N, const std::complex<double> *x, const int xStride) {
		if (std::abs(xStride) == 1) {
			double sum2 = 0;
			// Aliasing crimes, but this specific one is as kosher as it gets
			vDSP_svesqD((const double *)x, std::abs(xStride), &sum2, N*2);
			return sum2;
		} else {
			// criiiimes
			double real2 = 0, imag2 = 0;
			vDSP_svesqD((const double *)x, std::abs(xStride)*2, &real2, N);
			vDSP_svesqD((const double *)x + 1, std::abs(xStride)*2, &imag2, N);
			return real2 + imag2;
		}
	}
	double norm2(const int N, ConstSplitComplex<double> x, const int xStride) {
		double real2 = 0, imag2 = 0;
		vDSP_svesqD(x.real, std::abs(xStride), &real2, N);
		vDSP_svesqD(x.imag, std::abs(xStride), &imag2, N);
		return real2 + imag2;
	}

	void norm2(const int N, const std::complex<double> *x, const int xStride, double *y, const int yStride) {
		if (xStride < 0) x += (1 - N)*xStride;
		if (yStride < 0) y += (1 - N)*yStride;
		
		auto splitX = dspSplit(x);
		vDSP_zvmagsD(&splitX, 2*xStride, y, yStride, N);
	}

	void norm2(const int N, ConstSplitComplex<double> x, const int xStride, double *y, const int yStride) {
		if (xStride < 0) x += (1 - N)*xStride;
		if (yStride < 0) y += (1 - N)*yStride;

		auto splitX = dspSplit(x);
		vDSP_zvmagsD(&splitX, xStride, y, yStride, N);
	}
	
	void reserve(size_t) {}
private:
};
*/
}}; // namespace

#undef CBLAS_INT
#undef CBLAS_S
#undef CBLAS_Z
