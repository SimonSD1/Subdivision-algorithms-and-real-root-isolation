# Compilateur et options
CC = gcc
CFLAGS = -Wall -Wextra -I src/HeaderFiles $(shell pkg-config --cflags flint)
LIBS = $(shell pkg-config --libs flint) -lmpfr -lgmp

# Répertoires
SRC_DIR = src
BIN_DIR = $(SRC_DIR)/bin				# Pour fichiers temporaires crées par les fonctions test
IMPL_DIR = $(SRC_DIR)/Implementations
VALIDITY_DIR = $(SRC_DIR)/ValidityTests
EFFICIENCY_DIR = $(SRC_DIR)/EfficiencyTests
OBJ_DIR = $(SRC_DIR)/obj
DATA_DIR = DATA

# Répertoires de database à clean
DATA_DIRS = $(DATA_DIR)/Poly_ChangingCoeffSize \
            $(DATA_DIR)/Poly_ChangingDegree \
            $(DATA_DIR)/Precomputation_DivConq


# Fichiers d'implémentation
IMPL_FILES = $(wildcard $(IMPL_DIR)/*.c)
IMPL_OBJECTS = $(IMPL_FILES:$(IMPL_DIR)/%.c=$(OBJ_DIR)/%.o)

# Fichiers de test
VALIDITY_TESTS = $(wildcard $(VALIDITY_DIR)/*.c)
EFFICIENCY_TESTS = $(wildcard $(EFFICIENCY_DIR)/*.c)

# Exécutables de test (produits dans le répertoire courant)
VALIDITY_EXECS = $(VALIDITY_TESTS:$(VALIDITY_DIR)/%.c=%)
EFFICIENCY_EXECS = $(EFFICIENCY_TESTS:$(EFFICIENCY_DIR)/%.c=%)

# Liste complète des exécutables pour le nettoyage
CLEAN_EXECS = $(VALIDITY_EXECS) $(EFFICIENCY_EXECS)


# Cibles principales
all: validity efficiency

# Cibles pour compiler tous les tests de validité et d'efficacité
validity: $(BIN_DIR) $(VALIDITY_EXECS)
efficiency: $(BIN_DIR) $(EFFICIENCY_EXECS)


# Création du dossier src/bin
$(BIN_DIR):
	mkdir -p $@

# Création du dossier obj/
$(OBJ_DIR):
	mkdir -p $@



# Compilation des fichiers d'implémentation .o
$(OBJ_DIR)/%.o: $(IMPL_DIR)/%.c | $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

%: $(VALIDITY_DIR)/%.c $(IMPL_OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS) $< $(IMPL_OBJECTS) $(LIBS) -o $@

%: $(EFFICIENCY_DIR)/%.c $(IMPL_OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS) $< $(IMPL_OBJECTS) $(LIBS) -o $@



clean:
	rm -rf $(OBJ_DIR)
	rm -rf $(BIN_DIR)
	rm -f $(CLEAN_EXECS)

clean_database:
	@for dir in $(DATA_DIRS); do \
		if [ -d "$$dir" ]; then \
			echo "Nettoyage du répertoire $$dir..."; \
			rm -f "$$dir"/*; \
			echo "Contenu supprimé de $$dir"; \
		else \
			echo "Le répertoire $$dir n'existe pas, ignoré."; \
		fi \
	done


# Rend les cibles non-fichiers
.PHONY: all validity efficiency clean clean_database