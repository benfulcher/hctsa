#define SMALL 0.0001

struct list{
   long n;
   double *tcm;
   double **tse;
   struct list *pre;
};

long neigh(double *, long, long, long, long, long, double, long *);
long train(long, double *, long, long, long, double **, long, long, long, 
   double, double *);
void clean(long, double *, double *, long, long, long, long, long, long, 
   double, double, double *, long, long *);
