# Compilateur et options
CC = gcc
CFLAGS = -Wall -Wextra -IHeaderFiles
LIBS = -lflint -lgmp

# Répertoires
IMPL_DIR = Implementations
VALIDITY_DIR = ValidityTests
EFFICIENCY_DIR = EfficiencyTests
OBJ_DIR = obj

# Création des fichiers .o pour les implémentations
IMPL_FILES = $(wildcard $(IMPL_DIR)/*.c)
IMPL_OBJECTS = $(IMPL_FILES:$(IMPL_DIR)/%.c=$(OBJ_DIR)/%.o)

# Fichiers de test
VALIDITY_TESTS = $(wildcard $(VALIDITY_DIR)/*.c)
EFFICIENCY_TESTS = $(wildcard $(EFFICIENCY_DIR)/*.c)

# Exécutables (placés dans src/)
VALIDITY_EXECS = $(VALIDITY_TESTS:$(VALIDITY_DIR)/%.c=%)
EFFICIENCY_EXECS = $(EFFICIENCY_TESTS:$(EFFICIENCY_DIR)/%.c=%)

# Cible principale
all: $(VALIDITY_EXECS) $(EFFICIENCY_EXECS)

# Création du dossier obj si nécessaire
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# Compilation des fichiers d'implémentation en .o
$(OBJ_DIR)/%.o: $(IMPL_DIR)/%.c | $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Compilation des fichiers de test en exécutables dans src/
%: $(VALIDITY_DIR)/%.c $(IMPL_OBJECTS) | $(OBJ_DIR)
	$(CC) $(CFLAGS) $< $(IMPL_OBJECTS) $(LIBS) -o $@

%: $(EFFICIENCY_DIR)/%.c $(IMPL_OBJECTS) | $(OBJ_DIR)
	$(CC) $(CFLAGS) $< $(IMPL_OBJECTS) $(LIBS) -o $@

# Nettoyage
clean:
	rm -rf $(OBJ_DIR) $(VALIDITY_EXECS) $(EFFICIENCY_EXECS)
