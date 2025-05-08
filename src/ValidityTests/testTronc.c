#include <stdio.h>
#include <flint/fmpz.h>

void print_bits(const fmpz_t x) {
    char *bin = fmpz_get_str(NULL, 2, x);
    printf("binaire : %s\n", bin);
    flint_free(bin);
}

int main() {
    fmpz_t coeff1, coeff2;
    fmpz_init(coeff1);
    fmpz_init(coeff2);

    // Exemple avec un entier positif
    fmpz_set_si(coeff1, 123456789);  // essaie aussi avec un négatif

    flint_bitcnt_t bit_coeff = fmpz_bits(coeff1);
    slong truncate_bits = 5;

    printf("=== Avant troncature ===\n");
    gmp_printf("coeff1 = %ld\n", coeff1);
    printf("bit length = %lu\n", (unsigned long)bit_coeff);
    print_bits(coeff1);

    // Appliquer la troncature (tentative de ne garder QUE les bits de poids faible)
    fmpz_tdiv_r_2exp(coeff2, coeff1, bit_coeff - truncate_bits);

    printf("\n=== Après troncature (tdiv_r_2exp) ===\n");
    gmp_printf("coeff2 = %ld\n", coeff2);
    print_bits(coeff2);

    fmpz_clear(coeff1);
    fmpz_clear(coeff2);
    return 0;
}
