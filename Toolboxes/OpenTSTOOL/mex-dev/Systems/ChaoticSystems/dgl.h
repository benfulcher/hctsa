#ifndef DGL_H
#define DGL_H

#include <cmath>

using namespace std;

class DGL {
	protected:
		const long number_equations;	// we don't expect the number of equations to change
	public:
		double* v;		// vector of function values
		DGL(const long n) : number_equations(n), v(0) {
			v = new double[n];	
		};
		~DGL() {
			delete[] v;
		}
		long get_number_equations() const { return number_equations; };
};

class LorenzDGL : public DGL {
        private:
			const double a, b, c;	// parameters            
        public:
			LorenzDGL(const double A = -10, const double B = 28, const double C = -2.6666666) : 
        		DGL(3), a(A), b(B), c(C) { 
			};
			~LorenzDGL() {};                 
			void operator()(const double* const x, double* y, double* f) const {
        		f[0] = a*(y[0]-y[1]); 
        		f[1] = b*y[0]-y[1]-y[0]*y[2];
        		f[2] = y[0]*y[1]+c*y[2];                
			};
};

// From Johan A.K. Suykens, Joos Vanderwalle : Nonlinear Modelling
// The K.U. Leuven Time Series Prediction Competition
// The test data set was constructed by integrating the generalized 5-scroll (see derived class FiveScroll)
// Chua system and applying the follwing transformation R^3 -> R
//
// W = [-.0124 0.3267 1.2288];
// V = [-0.1004 -0.1102 -0.2784; 0.0009 0.5792 0.6892; 0.1063 -0.0042 0.0943];
// 
// If x is a state vector :  y = W * tanh(V * x');
//
class GeneralizedChuaDGL  : public DGL {		// this is an abstract base class			
	private:
		const double alpha, beta;
		
		const long q;		// degree of 
		
		double* m;  	  // double vector of size 2*q 
		double* c;  	  // double vector of size 2*q - 1
		
	public:
		GeneralizedChuaDGL(const double a, const double b, const long Q, double* M, double* C) : 
			alpha(a), beta(b), DGL(3), q(Q), m(M), c(C) { 
			
			v[0] = 0.1; v[1] = -0.2; v[2] = -0.3;	// give some default initial values		
		};	
		~GeneralizedChuaDGL() {};
				
		inline double h(const double x) const {
			double y = m[2*q-1] * x;
			
			for (long i = 1; i <= 2*q-1; i++) {
				y += 0.5 * ((m[i-1]-m[i]) * (fabs(x + c[i-1]) - fabs(x - c[i-1]))); 
			}
			
			return y; 	
		}
		
		void operator()(double* x, double* y, double* f) const {
			f[0] = alpha * ( y[1] - h(y[0]) );
			f[1] = y[0] - y[1]  + y[2];
			f[2] = - beta * y[1];
		}
};

class DoubleScroll : public GeneralizedChuaDGL {
	private:	
		double M[2];
		double C[1];
		
	public:
		DoubleScroll(const double a = 9.0, const double b = 14.286) :
			GeneralizedChuaDGL(a, b, 1, M, C) {
			
			M[0] = -1.0/7; M[1] = +2.0/7; 
			C[0] = 1; 
		}	
		~DoubleScroll() {};
};

class FiveScroll : public GeneralizedChuaDGL {
	private:	
		double M[6];
		double C[5];
	
	public:
		FiveScroll(const double a = 9.0, const double b = 14.286) :
			GeneralizedChuaDGL(a, b, 3, M, C) {
			
			M[0] = 0.9/7; M[1] = -3.0/7; M[2] = 3.5/7; M[3] = -2.7/7; M[4] = 4.0/7; M[5] = -2.4/7;
			C[0] = 1; C[1] = 2.15; C[2] = 3.6; C[3] = 6.2; C[4] = 9.0;	
		}	
		~FiveScroll() {};
};

class DuffingDGL : public DGL {
	private:
		const double A , B, C;

	public:	
		DuffingDGL(const double a = 40, const double b = 0.2, const double c = 1) : 
			A(a), B(b), C(c), DGL(3) {};
		
		~DuffingDGL() {};

		void operator()(double* x, double* y, double* f) const {
			f[0] = y[1];
			f[1] = - y[0] - y[0]*y[0]*y[0] - B*y[1] + A*cos(y[2]);
			f[2] = C;
		}
};

class RoesslerDGL : public DGL {
	private:
		const double A, B, C;
				
	public:	
		RoesslerDGL(const double a = 0.45, const double b = 2.0, const double c = 4) : 
			A(a), B(b), C(c), DGL(3) {};
		~RoesslerDGL() {};
	
		void operator()(double* x, double* y, double* f) const {
			f[0] = -y[1] - y[2];
			f[1] = y[0] + A * y[1];
			f[2] = B + y[2] * ( y[0] - C);
		}
};


// Normalized Colpitts, see "Bifurcation Phenomena in the Colpitts Oscillator: A robustness analysis", Oscar De Feo, Gian Maria Maggio

 class NormalizedColpittsDGL : public DGL {
	private:
		const double G, Q, K;
				
	public:	
		NormalizedColpittsDGL(const double g, const double q, const double k) : 
			G(g), Q(q), K(k), DGL(3) {};
		~NormalizedColpittsDGL() {};
	
		void operator()(double* x, double* y, double* f) const {
		    f[0] = G/(Q*(1-K))*(-exp(-y[1])+1+y[2]);
		    f[1] = G/(Q*K)*y[2];
		    f[2] = -Q*K*(1-K)/G*(y[0]+y[1])-1/Q*y[2];
		}
};

// Colpitts, see "Bifurcation Phenomena in the Colpitts Oscillator: A robustness analysis", Oscar De Feo, Gian Maria Maggio

class ColpittsDGL : public DGL {
	private:
		const double C1, C2, L, R, R1, Vc, Ve;
				
	public:	
		ColpittsDGL(const double c1, const double c2, const double l, const double r, const double r1, const double vc,
			    const double ve) : 
			C1(c1), C2(c2), L(l), R(r), R1(r1), Vc(vc), Ve(ve), DGL(3) {};
		~ColpittsDGL() {};

		void operator()(double* x, double* y, double* f) const {
		    f[0] = (y[2]-14e-15/200*exp(y[1]/0.026))/C1;
		    f[1] = (-(Ve+y[1])/R-y[2]-14e-15/200*exp(y[1]/0.026))/C2;
		    f[2] = (-Vc-y[0]+y[1]-R1*y[2])/L;
		}
};

class TodaDGL : public DGL {
	private:
		const double R, A, W;

	public:
		TodaDGL(const double r, const double a, const double w) :
			R(r), A(a), W(w), DGL(2) {};
		~TodaDGL() {};
		
		void operator()(double* x, double* y, double* f) const {
			f[0] = 1 + A * sin(W * (*x)) - R * y[1] - exp(y[0]);
			f[1] = y[0];			
		}
};


class VanDerPolDGL : public DGL {
	private:
		const double E, W0, W, A;

	public:
		VanDerPolDGL(const double e, const double wo, const double w, const double a) :
			E(e), W0(wo), W(w), A(a), DGL(2) {};
		~VanDerPolDGL() {};
		
		void operator()(double* x, double* y, double* f) const {
			f[0] = A * sin(W * (*x)) - E * (y[1] * y[1] - 1) - W0 * W0 * y[0];
			f[1] = y[0];			
		}
};


class PendelumDGL : public DGL {
	private:
		const double ALPHA, BETA, A, W;

	public:
		PendelumDGL(const double alpha, const double beta, const double a, const double w) :
			ALPHA(alpha) , BETA(beta), A(a), W(w), DGL(2) {};
		~PendelumDGL() {};
		
		void operator()(double* x, double* y, double* f) const {
			f[0] = A * sin(W * (*x)) - ALPHA * y[1] - BETA * sin(y[0]);
			f[1] = y[0];			
		}
};

class HarmonicDGL : public DGL {
	private:
		double w;
	
	public:	
		HarmonicDGL(const double W) : w(W), DGL(2) {};
		
		~HarmonicDGL() {};
		
		void increase_w(const double increment) { w += increment; }
				
		void operator()(double* x, double* y, double* f) const {
			f[0] = y[1];
			f[1] = - (w * w * y[0]);
		}
};

class BaierSahleDGL : public DGL {
	private:
		const double A, B, D, EPSILON;
				
	public:	
		BaierSahleDGL(const long n = 3, const double a = 0.28, const double b = 4.0, const double d = 2.0, 
			const double epsilon = 0.1) : 
			A(a), B(b), D(d), EPSILON(epsilon), DGL(n) {};
		~BaierSahleDGL() {};
	
		void operator()(double* x, double* y, double* f) const {
			f[0] = -y[1] + A * y[0];
			for (long i=1; i < number_equations-1; i++)
				f[i] = y[i-1] - y[i+1]; 
				
			f[number_equations-1] = EPSILON + B * y[number_equations-1] * (y[number_equations-2]  - D);
		}
};




#endif
