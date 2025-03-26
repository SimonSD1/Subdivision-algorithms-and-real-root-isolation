#include "../HeaderFiles/isolation.h"


int main()
{
    fmpz_poly_t pol;

    fmpz_poly_init(pol);

    fmpz_poly_set_str(pol,"4  2 -3 0 5");

    fmpz_poly_print_pretty(pol,"x");

    div_by_x(pol);

    printf("\n");

    fmpz_poly_print_pretty(pol,"x");

    solution* solutions = malloc(10 *sizeof(solution));
    
    return 0;
}
