#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "helper_functions.h"

double min_(const double a[], const int size)
{
    double m = a[0];
    for (int i = 1; i < size; i++) {
        if (a[i] < m) {
            m = a[i];
        }
    }
    return m;
}

double max_(const double a[], const int size)
{
    double m = a[0];
    for (int i = 1; i < size; i++) {
        if (a[i] > m) {
            m = a[i];
        }
    }
    return m;
}

double mean(const double a[], const int size)
{
    double m = 0.0;
    for (int i = 0; i < size; i++) {
        m += a[i];
    }
    m /= size;
    return m;
}

double sum(const double a[], const int size)
{
    double m = 0.0;
    for (int i = 0; i < size; i++) {
        m += a[i];
    }
    return m;
}

void cumsum(const double a[], const int size, double b[])
{
    b[0] = a[0];
    for (int i = 1; i < size; i++) {
        b[i] = a[i] + b[i-1];
        //printf("b[%i]%1.3f = a[%i]%1.3f + b[%i-1]%1.3f\n", i, b[i], i, a[i], i, a[i-1]);
    }
    
}

void icumsum(const int a[], const int size, int b[])
{
    b[0] = a[0];
    for (int i = 1; i < size; i++) {
        b[i] = a[i] + b[i-1];
        //printf("b[%i]%1.3f = a[%i]%1.3f + b[%i-1]%1.3f\n", i, b[i], i, a[i], i, a[i-1]);
    }
    
}

double isum(const int a[], const int size)
{
    double m = 0.0;
    for (int i = 0; i < size; i++) {
        m += a[i];
    }
    return m;
}

double median(const double a[], const int size)
{
    double m;
    double * b = malloc(size * sizeof *b);
    memcpy(b, a, size * sizeof *b);
    sort(b, size);
    if (size % 2 == 1) {
        m = b[size / 2];
    } else {
        int m1 = size / 2;
        int m2 = m1 - 1;
        m = (b[m1] + b[m2]) / (double)2.0;
    }
    free(b);
    return m;
}

double stddev(const double a[], const int size)
{
    double m = mean(a, size);
    double sd = 0.0;
    for (int i = 0; i < size; i++) {
        sd += pow(a[i] - m, 2);
    }
    sd = sqrt(sd / (size - 1));
    return sd;
}

double cov(const double x[], const double y[], const int size){
    
    double covariance = 0;
    
    double meanX = mean(x, size);
    double meanY = mean(y, size);
    
    for(int i = 0; i < size; i++){
        // double xi =x[i];
        // double yi =y[i];
        covariance += (x[i] - meanX) * (y[i] - meanY);
        
    }
    
    return covariance/(size-1);
    
}

double cov_mean(const double x[], const double y[], const int size){
    
    double covariance = 0;
    
    for(int i = 0; i < size; i++){
        // double xi =x[i];
        // double yi =y[i];
        covariance += x[i] * y[i];
        
    }
    
    return covariance/size;
    
}

double corr(const double x[], const double y[], const int size){
    
    double nom = 0;
    double denomX = 0;
    double denomY = 0;
    
    double meanX = mean(x, size);
    double meanY = mean(y, size);
    
    for(int i = 0; i < size; i++){
        nom += (x[i] - meanX) * (y[i] - meanY);
        denomX += (x[i] - meanX) * (x[i] - meanX);
        denomY += (y[i] - meanY) * (y[i] - meanY);
        
        //printf("x[%i]=%1.3f, y[%i]=%1.3f, nom[%i]=%1.3f, denomX[%i]=%1.3f, denomY[%i]=%1.3f\n", i, x[i], i, y[i], i, nom, i, denomX, i, denomY);
    }
    
    return nom/sqrt(denomX * denomY);
    
}

double autocorr_lag(const double x[], const int size, const int lag){
    
    return corr(x, &(x[lag]), size-lag);
    
}

double autocov_lag(const double x[], const int size, const int lag){
    
    return cov_mean(x, &(x[lag]), size-lag);
    
}

void zscore_norm(double a[], int size)
{
    double m = mean(a, size);
    double sd = stddev(a, size);
    for (int i = 0; i < size; i++) {
        a[i] = (a[i] - m) / sd;
    }
    return;
}

void zscore_norm2(const double a[], const int size, double b[])
{
    double m = mean(a, size);
    double sd = stddev(a, size);
    for (int i = 0; i < size; i++) {
        b[i] = (a[i] - m) / sd;
    }
    return;
}

double moment(const double a[], const int size, const int start, const int end, const int r)
{
    int win_size = end - start + 1;
    a += start;
    double m = mean(a, win_size);
    double mr = 0.0;
    for (int i = 0; i < win_size; i++) {
        mr += pow(a[i] - m, r);
    }
    mr /= win_size;
    mr /= stddev(a, win_size); //normalize
    return mr;
}

void diff(const double a[], const int size, double b[])
{
    for (int i = 1; i < size; i++) {
        b[i - 1] = a[i] - a[i - 1];
    }
}

int linreg(const int n, const double x[], const double y[], double* m, double* b) //, double* r)
{
    double   sumx = 0.0;                      /* sum of x     */
    double   sumx2 = 0.0;                     /* sum of x**2  */
    double   sumxy = 0.0;                     /* sum of x * y */
    double   sumy = 0.0;                      /* sum of y     */
    double   sumy2 = 0.0;                     /* sum of y**2  */
    
    /*
    for (int i = 0; i < n; i++)
    {
        fprintf(stdout, "x[%i] = %f, y[%i] = %f\n", i, x[i], i, y[i]);
    }
    */
    
    for (int i=0;i<n;i++){
        sumx  += x[i];
        sumx2 += x[i] * x[i];
        sumxy += x[i] * y[i];
        sumy  += y[i];
        sumy2 += y[i] * y[i];
    }
    
    double denom = (n * sumx2 - sumx * sumx);
    if (denom == 0) {
        // singular matrix. can't solve the problem.
        *m = 0;
        *b = 0;
        //if (r) *r = 0;
        return 1;
    }
    
    *m = (n * sumxy  -  sumx * sumy) / denom;
    *b = (sumy * sumx2  -  sumx * sumxy) / denom;
    
    /*if (r!=NULL) {
        *r = (sumxy - sumx * sumy / n) /    // compute correlation coeff
        sqrt((sumx2 - sumx * sumx/n) *
             (sumy2 - sumy * sumy/n));
    }
    */
    
    return 0;
}

double norm_(const double a[], const int size)
{
    
    double out = 0.0;
    
    for (int i = 0; i < size; i++)
    {
        out += a[i]*a[i];
    }
    
    out = sqrt(out);
    
    return out;
}
