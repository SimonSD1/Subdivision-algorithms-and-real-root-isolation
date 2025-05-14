#include <stdio.h>
#include <flint/fmpz.h>
#include "../HeaderFiles/coeff_truncation.h"
#include "../HeaderFiles/taylorShift_implem.h"

void print_bits(const fmpz_t x)
{
    char *bin = fmpz_get_str(NULL, 2, x);
    printf("binaire : %s\n", bin);
    flint_free(bin);
}

int main()
{
    fmpz_poly_t poly, poly_trunc;

    fmpz_poly_init(poly);
    fmpz_poly_init(poly_trunc);
    fmpz_t precision_limit, one;
    fmpz_init(precision_limit);

    fmpz_init_set_si(one, 1);

    double tps[101];
    double tps_trunc[101];

    clock_t start, end;

    for (int i = 0; i < 101; i++)
    {

        int nb_unreliable = 0;
        int unreliable_mismatch = 0;
        int nb_mismatch = 0;

        readPolyDATA(poly, 0, i);
        readPolyDATA(poly_trunc, 0, i);

        // on tronque le poly_trunc de d coeff

        slong max_bits = fmpz_poly_max_bits(poly);

        if (max_bits < 0)
        {
            max_bits = -max_bits;
        }

        slong degree = poly_trunc->length - 1;

        slong remainding_bit_target = max_bits / 2;

        slong to_truncate = max_bits - remainding_bit_target;

        if (to_truncate < 0)
        {
            to_truncate = 0;
        }

        fmpz_set_si(precision_limit, 1);
        fmpz_mul_2exp(precision_limit, precision_limit, degree + 1);

        truncate_coefficients(poly_trunc, poly_trunc, to_truncate);

        // on fait les taylor shift

        start = clock();
        fmpz_poly_taylor_shift(poly, poly, one);
        end = clock();
        tps[i] = (double)(end - start);

        start = clock();
        fmpz_poly_taylor_shift(poly_trunc, poly_trunc, one);
        end = clock();
        tps_trunc[i] = (double)(end - start);

        // on compare les signes obtenus en prenant en compte la perte de precision
        // on ne prend en compte que les coeffs superieur a 2^(d+1)

        for (slong j = 0; j < poly->length; j++)
        {

            fmpz *coeff_poly = &(poly->coeffs[j]);
            fmpz *coeff_poly_trunc = &(poly_trunc->coeffs[j]);

            int sign_poly = fmpz_cmp_ui(coeff_poly, 0);

            // printf("\n");
            // fmpz_print(coeff_poly);
            // printf("\n");
            int sign_truncated = fmpz_cmp_ui(coeff_poly_trunc, 0);

            if (fmpz_cmpabs(coeff_poly_trunc, precision_limit) < 0)
            {
                nb_unreliable++;
            }

            if (sign_poly != sign_truncated)
            {
                nb_mismatch++;
                if (fmpz_cmpabs(coeff_poly_trunc, precision_limit) < 0)
                {
                    unreliable_mismatch++;
                }
                else
                {
                    printf("erreur !!!!!");
                }
            }
        }

        printf("\nnb_mismatch=%d, nb_unreliable=%d, nb_unreliable_mismatch=%d\n, degree=%ld", nb_mismatch, nb_unreliable, unreliable_mismatch, degree);
    }

    FILE *fileResults;
    fileResults = fopen("src/EfficiencyTests/Results/Truncation/Truncation_ChangingDegree.txt", "w");
    if (fileResults == NULL)
    {
        printf("The file is not opened. The program will "
               "now exit.\n");
        exit(0);
    }

    fprintf(fileResults, "Coeff truncation\n"); // Title of the plot
    fprintf(fileResults, "untruncated\n");      // labels of the plot
    fprintTab(tps, 101, fileResults);
    fprintf(fileResults, "truncated\n"); // labels of the plot
    fprintTab(tps_trunc, 101, fileResults);
    fclose(fileResults);

    system("python3 src/EfficiencyTests/plotGenerator.py src/EfficiencyTests/Results/Truncation/Truncation_ChangingDegree.txt 'time' 0");
}
