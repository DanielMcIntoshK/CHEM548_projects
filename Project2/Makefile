SRC_DIR := src
OBJ_DIR := obj
BIN_DIR := bin
INC_DIR := include
LIN_INC := /usr/include
EXE := Project2
SRC := $(wildcard $(SRC_DIR)/*.cxx)
OBJ := $(SRC:$(SRC_DIR)/%.cxx=$(OBJ_DIR)/%.o)
CPPFLAGS := -I$(INC_DIR) 

all: $(EXE)

CC := g++
CC_FLGS := -O2 -std=c++23 -pipe

.PHONY: all

$(EXE): $(OBJ) | $(BIN_DIR)
	$(CC) $^ -o $(BIN_DIR)/$@ 
$(BIN_DIR):
	mkdir -p $@
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cxx | $(OBJ_DIR) 
	$(CC) $(CPPFLAGS) $(CC_FLGS) -w -c $< -o $@
$(OBJ_DIR):
	mkdir -p $@

clean:
	@$(RM) -rv $(BIN_DIR) $(OBJ_DIR)
