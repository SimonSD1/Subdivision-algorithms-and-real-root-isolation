#include "../HeaderFiles/isolation.h"

int main()
{
    fmpz_poly_t pol;
    fmpz_poly_init(pol);
    
    // Définition du polynôme comme dans SageMath (*840 pour ne pas avoir de div)
    fmpz_poly_set_str(pol, "6  840 -26616 26546 -9253 1354 -71");
    
    printf("Polynôme : ");
    fmpz_poly_print_pretty(pol, "x");
    printf("\n");
    // Transformation pol(x*20)
    fmpz_poly_t transformed_pol;
    fmpz_poly_init(transformed_pol);
    
    fmpz_poly_t x_scaled;
    fmpz_poly_init(x_scaled);
    fmpz_poly_set_coeff_si(x_scaled, 1, 20); // 20 = borne donner a la main selon desmos
    
    fmpz_poly_compose(transformed_pol, pol, x_scaled);
    
    printf("Polynôme transformé : ");
    fmpz_poly_print_pretty(transformed_pol, "x");
    printf("\n");
    
    solution* isol = malloc(10 * sizeof(solution)); 
    int nb_sol = 0;
    
    isolation(transformed_pol, 0, 0, isol, &nb_sol);
    


    // Affichage des intervalles isolés
    printf("Intervalles isolés :\n");
    for (int i = 0; i < nb_sol; i++) {
        double denom = (1 << isol[i].k); 
        printf("[ %lf , %lf ]\n", isol[i].c / denom, (isol[i].c + isol[i].h) / denom);
    }
    
    // Libération de la mémoire
    free(isol);
    fmpz_poly_clear(pol);
    fmpz_poly_clear(transformed_pol);
    fmpz_poly_clear(x_scaled);
    
    return 0;
}
