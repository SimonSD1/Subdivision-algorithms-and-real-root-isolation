#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include <stdlib.h>


void x_plus_1_powerNewton(fmpz_poly_t result, slong PowerOf2){
    // void fmpz_bin_uiui(fmpz_t f, ulong n, ulong k)
    // Sets f to the binomial coefficient (n,k)

    // void _fmpz_poly_shift_left(fmpz *res, const fmpz *poly, slong len, slong n)
    // Sets (res, len + n) to (poly, len) shifted left by coefficients.

    fmpz_poly_init(result);

    fmpz_t coef_binom;
    fmpz_init(coef_binom);

    for(slong i=0; i<PowerOf2+1; i++){
        fmpz_bin_uiui(coef_binom,PowerOf2,i);
        fmpz_poly_set_coeff_fmpz(result, i, coef_binom);
    }

    fmpz_clear(coef_binom);
}



void write_x_plus_1_power_Newton(slong maxSize) {
    mkdir("DATA", 0777);  // Create DATA if not exists
    mkdir("DATA/Precomputation_DivConq_Newton", 0777);

    FILE* filePrecomp;
    char pathFile[80];

    fmpz_poly_t temp;
    fmpz_poly_init(temp);

    for(slong i=0; i<=maxSize; i++) {
        slong exponent = 1 << i;
        sprintf(pathFile, "DATA/Precomputation_DivConq_Newton/precomputation%d.txt", (int) i);
        filePrecomp = fopen(pathFile, "w");
        
        if (filePrecomp == NULL) {
            printf("The file is not opened. The program will "
                "now exit.\n");
            exit(0);
        }

        x_plus_1_powerNewton(temp, exponent);
        fmpz_poly_fprint(filePrecomp, temp);
        fclose(filePrecomp);
    }

    fmpz_poly_clear(temp);
}





void write_x_plus_1_power(slong maxSize) {
    fmpz_poly_t *precomputed = flint_malloc((maxSize+1) * sizeof(fmpz_poly_t));
    if(!precomputed){
        printf("error allocating\n");
    }
    for (slong i = 0; i <= maxSize; i++)
    {
        fmpz_poly_init(precomputed[i]);
    }

    fmpz_poly_set_coeff_ui(precomputed[0], 1, 1);
    fmpz_poly_set_coeff_ui(precomputed[0], 0, 1);

    // square (x+a) m times by using the previous one
    for (slong i = 1; i <= maxSize; i++)
    {
        // square the previous one in the i-th
        fmpz_poly_sqr(precomputed[i], precomputed[i - 1]);
    }


    mkdir("DATA", 0777);  // Create DATA if not exists
    mkdir("DATA/Precomputation_DivConq", 0777);

    FILE* filePrecomp;
    char pathFile[80];

    for(slong i=0; i<=maxSize; i++) {
        sprintf(pathFile, "DATA/Precomputation_DivConq/precomputation%d.txt", (int) i);
        filePrecomp = fopen(pathFile, "w");

        if (filePrecomp == NULL) {
            printf("The file is not opened. The program will "
                "now exit.\n");
            exit(0);
        }

        fmpz_poly_fprint(filePrecomp, precomputed[i]);
        fclose(filePrecomp);
    }


    for (slong i = 0; i <= maxSize; i++)
        fmpz_poly_clear(precomputed[i]);
    flint_free(precomputed);
}




int main() {
    struct timespec start_ts, end_ts;
    clock_gettime(CLOCK_MONOTONIC, &start_ts);

    write_x_plus_1_power(14);       //because the sizes of our test polynomial database goes up to 2**15
    
    clock_gettime(CLOCK_MONOTONIC, &end_ts);
    double elapsed_time = (end_ts.tv_sec - start_ts.tv_sec) + (double)(end_ts.tv_nsec - start_ts.tv_nsec) / 1000000000.0;
    printf("Time taken for naive : %f seconds\n", elapsed_time);
    /*start = clock();
    write_x_plus_1_power_Newton(14);
    end = clock();
    elapsed_time = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Time taken for Newton: %f seconds\n", elapsed_time);*/

    //Time taken for naive: 4.769573 seconds
    //Time taken for Newton: 6.618666 seconds

    return 0;
}