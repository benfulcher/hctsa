function s = genbyode(system, len, samplerate, initcond, params, yunit)

% signal/genbyode
%
% Generate signal by integrating a set of ordinary differential equations,
% taken from a predefined set of systems
%
% si = genbyode(system, len, samplerate, initcond, params, yunit)
%
% len - length of signal, in seconds
% samplerate  - in Hertz
% initcond - vector of inital conditions
% params - vector of parameter values
%
% yunit is an optional argument
%
% system may be one of the following :
%      'Lorenz'
%      'Chua'
%      'Chua5Scroll'
%      'Duffing'
%      'Roessler'
%      'Toda'
%      'VanDerPol'
%      'Pendelum'
%      'BaierSahle'
%      'Colpitts'
%
% Integration is done by a ode solver using a adams pece scheme with local extrapolation 
% (see "The initial value problem"  by L. F. Shampine and M. K. Gordon.
%
% Here are the DGL used for the systems :
% 
% 
% class LorenzDGL : public DGL {
%         private:
% 			const double a, b, c;	// parameters            
%         public:
% 			LorenzDGL(const double A = -10, const double B = 28, const double C = -2.6666666) : 
%         		DGL(3), a(A), b(B), c(C) { 
% 			};
% 			~LorenzDGL() {};                 
% 			void operator()(const double* const x, double* y, double* f) const {
%         		f[0] = a*(y[0]-y[1]); 
%         		f[1] = b*y[0]-y[1]-y[0]*y[2];
%         		f[2] = y[0]*y[1]+c*y[2];                
% 			};
% };
% 
% // From Johan A.K. Suykens, Joos Vanderwalle : Nonlinear Modelling
% // The K.U. Leuven Time Series Prediction Competition
% // The test data set was constructed by integrating the generalized 5-scroll (see derived class FiveScroll)
% // Chua system and applying the follwing transformation R^3 -> R
% //
% // W = [-.0124 0.3267 1.2288];
% // V = [-0.1004 -0.1102 -0.2784; 0.0009 0.5792 0.6892; 0.1063 -0.0042 0.0943];
% // 
% // If x is a state vector :  y = W * tanh(V * x');
% //
% class GeneralizedChuaDGL  : public DGL {		// this is an abstract base class			
% 	private:
% 		const double alpha, beta;
% 		
% 		const long q;		// degree of 
% 		
% 		double* m;  	  // double vector of size 2*q 
% 		double* c;  	  // double vector of size 2*q - 1
% 		
% 	public:
% 		GeneralizedChuaDGL(const double a, const double b, const long Q, double* M, double* C) : 
% 			alpha(a), beta(b), DGL(3), q(Q), m(M), c(C) { 
% 			
% 			v[0] = 0.1; v[1] = -0.2; v[2] = -0.3;	// give some default initial values		
% 		};	
% 		~GeneralizedChuaDGL() {};
% 				
% 		inline double h(const double x) const {
% 			double y = m[2*q-1] * x;
% 			
% 			for (long i = 1; i <= 2*q-1; i++) {
% 				y += 0.5 * ((m[i-1]-m[i]) * (fabs(x + c[i-1]) - fabs(x - c[i-1]))); 
% 			}
% 			
% 			return y; 	
% 		}
% 		
% 		void operator()(double* x, double* y, double* f) const {
% 			f[0] = alpha * ( y[1] - h(y[0]) );
% 			f[1] = y[0] - y[1]  + y[2];
% 			f[2] = - beta * y[1];
% 		}
% };
% 
% class DoubleScroll : public GeneralizedChuaDGL {
% 	private:	
% 		double M[2];
% 		double C[1];
% 		
% 	public:
% 		DoubleScroll(const double a = 9.0, const double b = 14.286) :
% 			GeneralizedChuaDGL(a, b, 1, M, C) {
% 			
% 			M[0] = -1.0/7; M[1] = +2.0/7; 
% 			C[0] = 1; 
% 		}	
% 		~DoubleScroll() {};
% };
% 
% class FiveScroll : public GeneralizedChuaDGL {
% 	private:	
% 		double M[6];
% 		double C[5];
% 	
% 	public:
% 		FiveScroll(const double a = 9.0, const double b = 14.286) :
% 			GeneralizedChuaDGL(a, b, 3, M, C) {
% 			
% 			M[0] = 0.9/7; M[1] = -3.0/7; M[2] = 3.5/7; M[3] = -2.7/7; M[4] = 4.0/7; M[5] = -2.4/7;
% 			C[0] = 1; C[1] = 2.15; C[2] = 3.6; C[3] = 6.2; C[4] = 9.0;	
% 		}	
% 		~FiveScroll() {};
% };
% 
% class DuffingDGL : public DGL {
% 	private:
% 		const double A , B, C;
% 
% 	public:	
% 		DuffingDGL(const double a = 40, const double b = 0.2, const double c = 1) : 
% 			A(a), B(b), C(c), DGL(3) {};
% 		
% 		~DuffingDGL() {};
% 
% 		void operator()(double* x, double* y, double* f) const {
% 			f[0] = y[1];
% 			f[1] = - y[0] - y[0]*y[0]*y[0] - B*y[1] + A*cos(y[2]);
% 			f[2] = C;
% 		}
% };
% 
% class RoesslerDGL : public DGL {
% 	private:
% 		const double A, B, C;
% 				
% 	public:	
% 		RoesslerDGL(const double a = 0.45, const double b = 2.0, const double c = 4) : 
% 			A(a), B(b), C(c), DGL(3) {};
% 		~RoesslerDGL() {};
% 	
% 		void operator()(double* x, double* y, double* f) const {
% 			f[0] = -y[1]- y[2];
% 			f[1] = y[0] + A * y[1];
% 			f[2] = B + y[2] * ( y[0] - C);
% 		}
% };
% 
% class TodaDGL : public DGL {
% 	private:
% 		const double R, A, W;
% 
% 	public:
% 		TodaDGL(const double r, const double a, const double w) :
% 			R(r), A(a), W(w), DGL(2) {};
% 		~TodaDGL() {};
% 		
% 		void operator()(double* x, double* y, double* f) const {
% 			f[0] = 1 + A * sin(W * (*x)) - R * y[1] - exp(y[0]);
% 			f[1] = y[0];			
% 		}
% };
% 
% 
% class VanDerPolDGL : public DGL {
% 	private:
% 		const double E, W0, W, A;
% 
% 	public:
% 		VanDerPolDGL(const double e, const double wo, const double w, const double a) :
% 			E(e), W0(wo), W(w), A(a), DGL(2) {};
% 		~VanDerPolDGL() {};
% 		
% 		void operator()(double* x, double* y, double* f) const {
% 			f[0] = A * sin(W * (*x)) - E * (y[1] * y[1] - 1) - W0 * W0 * y[0];
% 			f[1] = y[0];			
% 		}
% };
% 
% 
% class PendelumDGL : public DGL {
% 	private:
% 		const double ALPHA, BETA, A, W;
% 
% 	public:
% 		PendelumDGL(const double alpha, const double beta, const double a, const double w) :
% 			ALPHA(alpha) , BETA(beta), A(a), W(w), DGL(2) {};
% 		~PendelumDGL() {};
% 		
% 		void operator()(double* x, double* y, double* f) const {
% 			f[0] = A * sin(W * (*x)) - ALPHA * y[1] - BETA * sin(y[0]);
% 			f[1] = y[0];			
% 		}
% };
% 
% class HarmonicDGL : public DGL {
% 	private:
% 		double w;
% 	
% 	public:	
% 		HarmonicDGL(const double W) : w(W), DGL(2) {};
% 		
% 		~HarmonicDGL() {};
% 		
% 		void increase_w(const double increment) { w += increment; }
% 				
% 		void operator()(double* x, double* y, double* f) const {
% 			f[0] = y[1];
% 			f[1] = - (w * w * y[0]);
% 		}
% };
% 
% 
%
% cmerk 1999

if nargin < 1
	system = 'Lorenz';		
end
if nargin < 2
	len = 1;			% one second 
end
if nargin < 3
	samplerate = 100;	% Hertz
end
if nargin < 4
	initcond = [0.1 -0.1 0];
end
if nargin < 5
	params = [-10 28 -2.666666];
end
if nargin < 6
	yunit = unit;		% no dimension
end


len = ceil(len*samplerate);
optinaltext = '';
disp(system);
switch system
	case	'Lorenz'
		tmp = chaosys(len, 1/samplerate, initcond, 0, params);
		optinaltext = [''];
	case	'Chua'
		tmp = chaosys(len, 1/samplerate, initcond, 1, params);
		optinaltext = [''];
	case	'Chua5Scroll'
		tmp = chaosys(len, 1/samplerate, initcond, 2, params);
		optinaltext = [''];
	case	'Duffing'
		tmp = chaosys(len, 1/samplerate, initcond, 3, params);
		optinaltext = [''];
	case	'Roessler'
		tmp = chaosys(len, 1/samplerate, initcond, 4, params);
		optinaltext = [''];
	case	'Toda'
		tmp = chaosys(len, 1/samplerate, initcond, 5, params);
		optinaltext = [''];
	case	'VanDerPol'
		tmp = chaosys(len, 1/samplerate, initcond, 6, params);
		optinaltext = [''];												
	case	'Pendelum'
		tmp = chaosys(len, 1/samplerate, initcond, 7, params);
		optinaltext = [''];			
	case	'BaierSahle'
		tmp = chaosys(len, 1/samplerate, initcond, 8, params);
		optinaltext = [''];	
	case	'Colpitts'
		tmp = chaosys(len, 1/samplerate, initcond, 9, params);
		optinaltext = [''];
	otherwise	
		error(['system ' system ' not supported']);
end	
	
s = signal(tmp, ['Generated ' system ' signal' optinaltext]);
s = setaxis(s, 1,achse(unit('s'), 0, 1/samplerate));	
s = setyunit(s, yunit);
if (size(tmp, 2) == 3) 
  s = setplothint(s, '3dcurve');
end
if (size(tmp, 2) == 2)
  s=setplothint(s,'xyplot');
end
      

