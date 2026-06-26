# Makefile for mini-sar-imaging
# SAR Imaging Module

CC := gcc
CFLAGS := -Wall -Wextra -std=c11 -O2 -I./include
LDFLAGS := -lm

SRC := src/sar_core.c src/sar_geometry.c src/sar_algorithm.c src/sar_interferometry.c src/sar_advanced.c
OBJ := $(SRC:.c=.o)

TEST_SRC := tests/test_sar.c
TEST_BIN := test_sar

EX1 := example_point_target
EX2 := example_stripmap
EX3 := example_insar

.PHONY: all test examples clean

all: $(OBJ) $(TEST_BIN)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

$(TEST_BIN): $(OBJ) $(TEST_SRC)
	$(CC) $(CFLAGS) $(TEST_SRC) $(OBJ) -o $(TEST_BIN) $(LDFLAGS)

examples: $(EX1) $(EX2) $(EX3)

$(EX1): $(OBJ) examples/example_point_target.c
	$(CC) $(CFLAGS) examples/example_point_target.c $(OBJ) -o $(EX1) $(LDFLAGS)

$(EX2): $(OBJ) examples/example_stripmap.c
	$(CC) $(CFLAGS) examples/example_stripmap.c $(OBJ) -o $(EX2) $(LDFLAGS)

$(EX3): $(OBJ) examples/example_insar.c
	$(CC) $(CFLAGS) examples/example_insar.c $(OBJ) -o $(EX3) $(LDFLAGS)

test: $(TEST_BIN)
	./$(TEST_BIN)

clean:
	rm -f $(OBJ) $(TEST_BIN) $(EX1) $(EX2) $(EX3) src/*.o
