# Compilateur et options
CC = gcc
CFLAGS = -Wall -Wextra -I src/HeaderFiles $(shell pkg-config --cflags flint)
LIBS = $(shell pkg-config --libs flint) -lmpfr -lgmp

# Répertoires
SRC_DIR = src
IMPL_DIR = $(SRC_DIR)/Implementations
VALIDITY_DIR = $(SRC_DIR)/ValidityTests
EFFICIENCY_DIR = $(SRC_DIR)/EfficiencyTests
OBJ_DIR = $(SRC_DIR)/obj
BIN_DIR = $(SRC_DIR)/bin
DATA_DIR = DATA

# Répertoires de données à nettoyer
DATA_DIRS = $(DATA_DIR)/Poly_ChangingCoeffSize \
            $(DATA_DIR)/Poly_ChangingDegree \
            $(DATA_DIR)/Precomputation_DivConq


# Fichiers d'implémentation
IMPL_FILES = $(wildcard $(IMPL_DIR)/*.c)
IMPL_OBJECTS = $(IMPL_FILES:$(IMPL_DIR)/%.c=$(OBJ_DIR)/%.o)

# Fichiers de test
VALIDITY_TESTS = $(wildcard $(VALIDITY_DIR)/*.c)
EFFICIENCY_TESTS = $(wildcard $(EFFICIENCY_DIR)/*.c)

# Exécutables de test
VALIDITY_EXECS = $(VALIDITY_TESTS:$(VALIDITY_DIR)/%.c=$(BIN_DIR)/%)
EFFICIENCY_EXECS = $(EFFICIENCY_TESTS:$(EFFICIENCY_DIR)/%.c=$(BIN_DIR)/%)

# Cibles principales
all: validity efficiency

validity: $(VALIDITY_EXECS)

efficiency: $(EFFICIENCY_EXECS)

# Création des dossiers nécessaires
$(OBJ_DIR) $(BIN_DIR):
	mkdir -p $@

# Compilation des fichiers d'implémentation en .o
$(OBJ_DIR)/%.o: $(IMPL_DIR)/%.c | $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Compilation des fichiers de test en exécutables dans bin/
$(BIN_DIR)/%: $(VALIDITY_DIR)/%.c $(IMPL_OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS) $< $(IMPL_OBJECTS) $(LIBS) -o $@

$(BIN_DIR)/%: $(EFFICIENCY_DIR)/%.c $(IMPL_OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS) $< $(IMPL_OBJECTS) $(LIBS) -o $@

# Nettoyage des exécutables
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

# Nettoyage des bases de données
# Utilise 'rm -f' après un cd pour supprimer de manière forcée le contenu, ce qui est plus robuste.
clean_database:
	@echo "Nettoyage des fichiers de données..."
	@for dir in $(DATA_DIRS); do \
		if [ -d "$$dir" ]; then \
			echo "Nettoyage du répertoire $$dir..."; \
			rm -f "$$dir"/*; \
			echo "Contenu supprimé de $$dir"; \
		else \
			echo "Le répertoire $$dir n'existe pas, ignoré."; \
		fi \
	done