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
    fmpz_fdiv_q_2exp(rop, op, drop_bits);  // shift right
    fmpz_mul_2exp(rop, rop, drop_bits);    // shift back left
}

void truncate_coefficients(fmpz_poly_t result, const fmpz_poly_t poly, slong precision_bits) {
    slong len = fmpz_poly_length(poly);
    fmpz_poly_set(result, poly);
    fmpz_t coeff;
    fmpz_init(coeff);

    for (slong i = 0; i < len; i++) {
        fmpz_poly_get_coeff_fmpz(coeff, result, i);
        
        // Tronquer les coefficients à la précision donnée
        // La troncature ici consiste à éliminer les bits au-delà de la précision spécifiée
        fmpz_trunc(coeff, coeff, precision_bits);  // Reduit à la précision spécifiée
        
        // Réassigner le coefficient tronqué dans le polynôme
        fmpz_poly_set_coeff_fmpz(result, i, coeff);
    }
    
    fmpz_clear(coeff);
}