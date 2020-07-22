#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "SB_CoarseGrain.h"
#include "helper_functions.h"

double SB_MotifThree_quantile_hh(const double y[], const int size)
{
    // NaN check
    for(int i = 0; i < size; i++)
    {
        if(isnan(y[i]))
        {
            return NAN;
        }
    }
    
    int tmp_idx, r_idx;
    int dynamic_idx;
    int alphabet_size = 3;
    int array_size;
    int * yt = malloc(size * sizeof(yt)); // alphabetized array
    double hh; // output
    double * out = malloc(124 * sizeof(out)); // output array
    
    // transfer to alphabet
    sb_coarsegrain(y, size, "quantile", 3, yt);
    
    // words of length 1
    array_size = alphabet_size;
    int ** r1 = malloc(array_size * sizeof(*r1));
    int * sizes_r1 = malloc(array_size * sizeof(sizes_r1));
    double * out1 = malloc(array_size * sizeof(out1));
    for (int i = 0; i < alphabet_size; i++) {
        r1[i] = malloc(size * sizeof(r1[i])); // probably can be rewritten
        // using selfresizing array for memory efficiency. Time complexity
        // should be comparable due to ammotization.
        r_idx = 0;
        sizes_r1[i] = 0;
        for (int j = 0; j < size; j++) {
            if (yt[j] == i + 1) {
                r1[i][r_idx++] = j;
                sizes_r1[i]++;
            }
        }
    }
    
    // words of length 2
    array_size *= alphabet_size;
    // removing last item if it is == max possible idx since later we are taking idx + 1
    // from yt
    for (int i = 0; i < alphabet_size; i++) {
        if (sizes_r1[i] != 0 && r1[i][sizes_r1[i] - 1] == size - 1) {
            //int * tmp_ar = malloc((sizes_r1[i] - 1) * sizeof(tmp_ar));
            int* tmp_ar = malloc(sizes_r1[i] * sizeof(tmp_ar));
            subset(r1[i], tmp_ar, 0, sizes_r1[i]);
            memcpy(r1[i], tmp_ar, (sizes_r1[i] - 1) * sizeof(tmp_ar));
            sizes_r1[i]--;
            free(tmp_ar);
        }
    }
    
    /*
    int *** r2 = malloc(array_size * sizeof(**r2));
    int ** sizes_r2 = malloc(array_size * sizeof(*sizes_r2));
    double ** out2 = malloc(array_size * sizeof(*out2));
    */
    int*** r2 = malloc(alphabet_size * sizeof(**r2));
    int** sizes_r2 = malloc(alphabet_size * sizeof(*sizes_r2));
    double** out2 = malloc(alphabet_size * sizeof(*out2));
    

    // allocate separately
    for (int i = 0; i < alphabet_size; i++) {
        r2[i] = malloc(alphabet_size * sizeof(*r2[i]));
        sizes_r2[i] = malloc(alphabet_size * sizeof(*sizes_r2[i]));
        //out2[i] = malloc(alphabet_size * sizeof(out2[i]));
        out2[i] = malloc(alphabet_size * sizeof(**out2));
        for (int j = 0; j < alphabet_size; j++) {
            r2[i][j] = malloc(size * sizeof(*r2[i][j]));
        }
    }

    // fill separately
    for (int i = 0; i < alphabet_size; i++) {
    // for (int i = 0; i < array_size; i++) {
        //r2[i] = malloc(alphabet_size * sizeof(r2[i]));
        //sizes_r2[i] = malloc(alphabet_size * sizeof(sizes_r2[i]));
        //out2[i] = malloc(alphabet_size * sizeof(out2[i]));
        for (int j = 0; j < alphabet_size; j++) {
            //r2[i][j] = malloc(size * sizeof(r2[i][j]));
            sizes_r2[i][j] = 0;
            dynamic_idx = 0; //workaround as you can't just add elements to array
            // like in python (list.append()) for example, so since for some k there will be no adding,
            // you need to keep track of the idx at which elements will be inserted
            for (int k = 0; k < sizes_r1[i]; k++) {
                tmp_idx = yt[r1[i][k] + 1];
                if (tmp_idx == (j + 1)) {
                    r2[i][j][dynamic_idx++] = r1[i][k];
                    sizes_r2[i][j]++;
                    // printf("dynamic_idx=%i, size = %i\n", dynamic_idx, size);
                }
            }
            double tmp = (double)sizes_r2[i][j] / ((double)(size) - (double)(1.0));
            out2[i][j] =  tmp;
        }
    }

    hh = 0.0;
    for (int i = 0; i < alphabet_size; i++) {
        hh += f_entropy(out2[i], alphabet_size);
    }

    free(yt);
    free(out);

    free(sizes_r1);

    // free nested array
    for (int i = 0; i < alphabet_size; i++) {
        free(r1[i]);
    }
    free(r1);
    // free(sizes_r1);
    
    for (int i = 0; i < alphabet_size; i++) {
    //for (int i = alphabet_size - 1; i >= 0; i--) {

        free(sizes_r2[i]);
        free(out2[i]);
    }

    //for (int i = alphabet_size-1; i >= 0 ; i--) {
    for(int i = 0; i < alphabet_size; i++) {
        for (int j = 0; j < alphabet_size; j++) {
            free(r2[i][j]);
        }
        free(r2[i]);
    }
    
    free(r2);
    free(sizes_r2);
    free(out2);
    
    
    return hh;
    
}

double * sb_motifthree(const double y[], int size, const char how[])
{
    int tmp_idx, r_idx, i, j, k, l, m, array_size;
    int dynamic_idx;
    int * tmp_ar;
    int alphabet_size = 3;
    int out_idx = 0;
    int * yt = malloc(size * sizeof(yt));
    double tmp;
    double * out = malloc(124 * sizeof(out)); // output array
    if (strcmp(how, "quantile") == 0) {
        sb_coarsegrain(y, size, how, alphabet_size, yt);
    } else if (strcmp(how, "diffquant") == 0) {
        double * diff_y = malloc((size - 1) * sizeof(diff_y));
        diff(y, size, diff_y);
        sb_coarsegrain(diff_y, size, how, alphabet_size, yt);
        size--;
    } else {
        fprintf(stdout, "ERROR in sb_motifthree: Unknown how method");
        exit(1);
    }

    // words of length 1
    array_size = alphabet_size;
    int ** r1 = malloc(array_size * sizeof(*r1));
    int * sizes_r1 = malloc(array_size * sizeof(sizes_r1));
    double * out1 = malloc(array_size * sizeof(out1));
    for (i = 0; i < array_size; i++) {
        r1[i] = malloc(size * sizeof(r1[i])); // probably can be rewritten
        // using selfresizing array for memory efficiency. Time complexity
        // should be comparable due to ammotization.
        r_idx = 0;
        sizes_r1[i] = 0;
        for (j = 0; j < size; j++) {
            if (yt[j] == i + 1) {
                r1[i][r_idx++] = j;
                sizes_r1[i]++;
            }
        }
        tmp = (double)sizes_r1[i] / size;

        out1[i] = tmp;
        out[out_idx++] = tmp;
    }
    out[out_idx++] = f_entropy(out1, array_size);

    // words of length 2
    array_size *= alphabet_size;
    // removing last item if it is == max possible idx since later we are taking idx + 1
    // from yt
    for (i = 0; i < alphabet_size; i++) {
        if (sizes_r1[i] != 0 && r1[i][sizes_r1[i] - 1] == size - 1) {
            tmp_ar = malloc((sizes_r1[i] - 1) * sizeof(tmp_ar));
            subset(r1[i], tmp_ar, 0, sizes_r1[i]);
            memcpy(r1[i], tmp_ar, (sizes_r1[i] - 1) * sizeof(tmp_ar));
            sizes_r1[i]--;
        }
    }

    int *** r2 = malloc(array_size * sizeof(**r2));
    int ** sizes_r2 = malloc(array_size * sizeof(*sizes_r2));
    double ** out2 = malloc(array_size * sizeof(*out2));
    for (i = 0; i < alphabet_size; i++) {
        r2[i] = malloc(alphabet_size * sizeof(r2[i]));
        sizes_r2[i] = malloc(alphabet_size * sizeof(sizes_r2[i]));
        out2[i] = malloc(alphabet_size * sizeof(out2[i]));
        for (j = 0; j < alphabet_size; j++) {
            r2[i][j] = malloc(size * sizeof(r2[i][j]));
            sizes_r2[i][j] = 0;
            dynamic_idx = 0; //workaround as you can't just add elements to array
            // like in python (list.append()) for example, so since for some k there will be no adding,
            // you need to keep track of the idx at which elements will be inserted
            for (k = 0; k < sizes_r1[i]; k++) {
                tmp_idx = yt[r1[i][k] + 1];
                if (tmp_idx == (j + 1)) {
                    r2[i][j][dynamic_idx++] = r1[i][k];
                    sizes_r2[i][j]++;
                }
            }
            tmp = (double)sizes_r2[i][j] / (size - 1);
            out2[i][j] = tmp;
            out[out_idx++] = tmp;
        }
    }
    tmp = 0.0;
    for (i = 0; i < alphabet_size; i++) {
        tmp += f_entropy(out2[i], alphabet_size);
    }
    out[out_idx++] = tmp;

    // words of length 3
    array_size *= alphabet_size;
    for (i = 0; i < alphabet_size; i++) {
        for (j = 0; j < alphabet_size; j++) {
            if (sizes_r2[i][j] != 0 && r2[i][j][sizes_r2[i][j] - 1] == size - 2) {
                subset(r2[i][j], tmp_ar, 0, sizes_r2[i][j]);
                memcpy(r2[i][j], tmp_ar, (sizes_r2[i][j] - 1) * sizeof(tmp_ar));
                sizes_r2[i][j]--;
            }
        }
    }

    int **** r3 = malloc(array_size * sizeof(***r3));
    int *** sizes_r3 = malloc(array_size * sizeof(**sizes_r3));
    double *** out3 = malloc(array_size * sizeof(**out3));
    for (i = 0; i < alphabet_size; i++) {
        r3[i] = malloc(alphabet_size * sizeof(r3[i]));
        sizes_r3[i] = malloc(alphabet_size * sizeof(sizes_r3[i]));
        out3[i] = malloc(alphabet_size * sizeof(out3[i]));
        for (j = 0; j < alphabet_size; j++) {
            r3[i][j] = malloc(alphabet_size * sizeof(r3[i][j]));
            sizes_r3[i][j] = malloc(alphabet_size * sizeof(sizes_r3[i][j]));
            out3[i][j] = malloc(alphabet_size * sizeof(out3[i][j]));
            for (k = 0; k < alphabet_size; k++) {
                r3[i][j][k] = malloc(size * sizeof(r3[i][j][k]));
                sizes_r3[i][j][k] = 0;
                dynamic_idx = 0;
                for (l = 0; l < sizes_r2[i][j]; l++) {
                    tmp_idx = yt[r2[i][j][l] + 2];
                    if (tmp_idx == (k + 1)) {
                        r3[i][j][k][dynamic_idx++] = r2[i][j][l];
                        sizes_r3[i][j][k]++;
                    }
                }
                tmp = (double)sizes_r3[i][j][k] / (size - 2);
                out3[i][j][k] = tmp;
                out[out_idx++] = tmp;
            }
        }
    }
    tmp = 0.0;
    for (i = 0; i < alphabet_size; i++) {
        for (j = 0; j < alphabet_size; j++) {
            tmp += f_entropy(out3[i][j], alphabet_size);
        }
    }
    out[out_idx++] = tmp;

    // words of length 4
    array_size *= alphabet_size;
    for (i = 0; i < alphabet_size; i++) {
        for (j = 0; j < alphabet_size; j++) {
            for (k = 0; k < alphabet_size; k++) {
                if (sizes_r3[i][j][k] != 0 && r3[i][j][k][sizes_r3[i][j][k] - 1] == size - 3) {
                    subset(r3[i][j][k], tmp_ar, 0, sizes_r3[i][j][k]);
                    memcpy(r3[i][j][k], tmp_ar, (sizes_r3[i][j][k] - 1) * sizeof(tmp_ar));
                    sizes_r3[i][j][k]--;
                }
            }
        }
    }

    int ***** r4 = malloc(array_size * sizeof(****r4));
    // just an array of pointers of array of pointers of array of pointers
    // of array of pointers of array of ints... We need to go deeper (c)
    int **** sizes_r4 = malloc(array_size * sizeof(***sizes_r3));
    double **** out4 = malloc(array_size * sizeof(***out4));
    for (i = 0; i < alphabet_size; i++) {
        r4[i] = malloc(alphabet_size * sizeof(r4[i]));
        sizes_r4[i] = malloc(alphabet_size * sizeof(sizes_r4[i]));
        out4[i] = malloc(alphabet_size * sizeof(out4[i]));
        for (j = 0; j < alphabet_size; j++) {
            r4[i][j] = malloc(alphabet_size * sizeof(r4[i][j]));
            sizes_r4[i][j] = malloc(alphabet_size * sizeof(sizes_r4[i][j]));
            out4[i][j] = malloc(alphabet_size * sizeof(out4[i][j]));
            for (k = 0; k < alphabet_size; k++) {
                r4[i][j][k] = malloc(alphabet_size * sizeof(r4[i][j][k]));
                sizes_r4[i][j][k] = malloc(alphabet_size * sizeof(sizes_r4[i][j][k]));
                out4[i][j][k] = malloc(alphabet_size * sizeof(out4[i][j][k]));
                for (l = 0; l < alphabet_size; l++) {
                    r4[i][j][k][l] = malloc(size * sizeof(r4[i][j][k][l]));
                    sizes_r4[i][j][k][l] = 0;
                    dynamic_idx = 0;
                    for (m = 0; m < sizes_r3[i][j][k]; m++) {
                        tmp_idx = yt[r3[i][j][k][m] + 3];
                        if (tmp_idx == l + 1) {
                            r4[i][j][k][l][dynamic_idx++] = r3[i][j][k][m];
                            sizes_r4[i][j][k][l]++;
                        }
                    }
                    tmp = (double)sizes_r4[i][j][k][l] / (size - 3);
                    out4[i][j][k][l] = tmp;
                    out[out_idx++] = tmp;
                }
            }
        }
    }
    tmp = 0.0;
    for (i = 0; i < alphabet_size; i++) {
        for (j = 0; j < alphabet_size; j++) {
            for (k = 0; k < alphabet_size; k++) {
                tmp += f_entropy(out4[i][j][k], alphabet_size);
            }
        }
    }
    out[out_idx++] = tmp;

    return out;
}

























