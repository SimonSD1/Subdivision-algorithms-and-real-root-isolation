#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include "../HeaderFiles/taylorShift_implem.h"
#include "../HeaderFiles/functionsForTests.h"
#include "../HeaderFiles/coeff_truncation.h"

int same_signs(const fmpz_poly_t poly1, const fmpz_poly_t poly2) {
    slong len_1 = fmpz_poly_length(poly1);
    slong len_2 = fmpz_poly_length(poly2);

    // Vérification rapide : mêmes longueurs ?
    if (len_1 != len_2) {
        return 0;
    }

    fmpz_t coeff_1, coeff_2;
    fmpz_init(coeff_1);
    fmpz_init(coeff_2);

    for (slong i = 0; i < len_1; i++) {
        fmpz_poly_get_coeff_fmpz(coeff_1, poly1, i);
        fmpz_poly_get_coeff_fmpz(coeff_2, poly2, i);

        int sign_a = fmpz_sgn(coeff_1);
        int sign_b = fmpz_sgn(coeff_2);

        if (sign_a != sign_b) {
            fmpz_clear(coeff_1);
            fmpz_clear(coeff_2);
            return 0;
        }
    }

    fmpz_clear(coeff_1);
    fmpz_clear(coeff_2);
    return 1;
}

void fmpz_trunc(fmpz_t rop, const fmpz_t op, slong keep_bits)
{
    slong total_bits = fmpz_bits(op);

    // Ne rien faire si l'entier est déjà plus petit que la précision cible
    if (total_bits <= keep_bits)
        return;

    slong drop_bits = total_bits - keep_bits;

    // Effectuer une troncature en supprimant les bits de poids faible
    //fmpz_fdiv_q_2exp(rop, op, drop_bits);  // shift right
    fmpz_mul_2exp(rop, op, drop_bits);    // shift back left
}

void truncate_coefficients(fmpz_poly_t result, const fmpz_poly_t poly) {
    slong len = fmpz_poly_length(poly);
    fmpz_poly_set(result, poly);
    fmpz_t coeff;
    fmpz_init(coeff);
    slong bit_coeff;
    fmpz_t deg;
    fmpz_init_set_si(deg, len-1);
    slong truncate_bits = (slong) (fmpz_bits(deg));
    //printf("truncating %ld bits\n", truncate_bits);

    for (slong i = 0; i < len; i++) {
        fmpz_poly_get_coeff_fmpz(coeff, result, i);
        bit_coeff = fmpz_bits(coeff);

        //Tronquer les bits de poids fort : on garde (bit_coeff-cutoff) bits de poids faible
        fmpz_tdiv_r_2exp(coeff, coeff, bit_coeff-truncate_bits);
        
        // Réassigner le coefficient tronqué dans le polynôme
        fmpz_poly_set_coeff_fmpz(result, i, coeff);
    }
    
    fmpz_clear(coeff);
    fmpz_clear(deg);
}