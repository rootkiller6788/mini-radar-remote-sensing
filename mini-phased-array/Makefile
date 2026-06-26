# Makefile for mini-phased-array
# Phased Array Beamforming and Radar Library

CC      := gcc
CFLAGS  := -Wall -Wextra -Wno-implicit-fallthrough -Wno-misleading-indentation -Wno-unused-function -std=c11 -O2 -I./include
LDFLAGS := -lm

CORE_SRC       := src/phased_array_core.c
STEERING_SRC   := src/phased_array_steering.c
BEAM_SRC       := src/phased_array_beam.c
SYNTHESIS_SRC  := src/phased_array_synthesis.c
ADAPTIVE_SRC   := src/phased_array_adaptive.c
AESA_SRC       := src/phased_array_aesa.c

ALL_SRC := $(CORE_SRC) $(STEERING_SRC) $(BEAM_SRC) \
           $(SYNTHESIS_SRC) $(ADAPTIVE_SRC) $(AESA_SRC)
ALL_OBJ := $(ALL_SRC:.c=.o)

TEST_SRC := tests/test_phased_array.c
TEST_BIN := test_phased_array

EX1_SRC := examples/example_linear_array.c
EX2_SRC := examples/example_planar_array.c
EX3_SRC := examples/example_adaptive_beamforming.c
EX1_BIN := example_linear_array
EX2_BIN := example_planar_array
EX3_BIN := example_adaptive_beamforming

.PHONY: all test examples clean

all: $(ALL_OBJ) $(TEST_BIN)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

$(TEST_BIN): $(ALL_OBJ) $(TEST_SRC)
	$(CC) $(CFLAGS) $(TEST_SRC) $(ALL_OBJ) -o $(TEST_BIN) $(LDFLAGS)

examples: $(EX1_BIN) $(EX2_BIN) $(EX3_BIN)

$(EX1_BIN): $(ALL_OBJ) $(EX1_SRC)
	$(CC) $(CFLAGS) $(EX1_SRC) $(ALL_OBJ) -o $(EX1_BIN) $(LDFLAGS)

$(EX2_BIN): $(ALL_OBJ) $(EX2_SRC)
	$(CC) $(CFLAGS) $(EX2_SRC) $(ALL_OBJ) -o $(EX2_BIN) $(LDFLAGS)

$(EX3_BIN): $(ALL_OBJ) $(EX3_SRC)
	$(CC) $(CFLAGS) $(EX3_SRC) $(ALL_OBJ) -o $(EX3_BIN) $(LDFLAGS)

test: $(TEST_BIN)
	./$(TEST_BIN)

clean:
	rm -f $(ALL_OBJ) $(TEST_BIN) $(EX1_BIN) $(EX2_BIN) $(EX3_BIN)
	rm -f src/*.o tests/*.o examples/*.o *.o
