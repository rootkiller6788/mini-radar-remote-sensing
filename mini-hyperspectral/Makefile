# Makefile — mini-hyperspectral
# Builds library, tests, and examples for hyperspectral remote sensing.

CC       = gcc
CFLAGS   = -std=c11 -Wall -Wextra -O2 -g -Iinclude
LDFLAGS  = -lm

# Library source files
LIB_SRCS = src/hyperspectral_core.c \
           src/hyperspectral_spectroscopy.c \
           src/hyperspectral_classification.c \
           src/hyperspectral_unmixing.c \
           src/hyperspectral_dimensionality.c \
           src/hyperspectral_radiometry.c

LIB_OBJS = $(LIB_SRCS:.c=.o)

# Test
TEST_SRC = tests/test_hyperspectral.c
TEST_BIN = tests/test_hyperspectral

# Examples
EX1_SRC = examples/example_mineral_detection.c
EX2_SRC = examples/example_vegetation_analysis.c
EX3_SRC = examples/example_spectral_unmixing.c
EX1_BIN = examples/example_mineral_detection
EX2_BIN = examples/example_vegetation_analysis
EX3_BIN = examples/example_spectral_unmixing

.PHONY: all clean test examples

all: test examples

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

$(TEST_BIN): $(TEST_SRC) $(LIB_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

test: $(TEST_BIN)
	$(TEST_BIN).exe

$(EX1_BIN): $(EX1_SRC) $(LIB_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(EX2_BIN): $(EX2_SRC) $(LIB_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(EX3_BIN): $(EX3_SRC) $(LIB_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

examples: $(EX1_BIN) $(EX2_BIN) $(EX3_BIN)

clean:
	rm -f src/*.o tests/test_hyperspectral tests/test_hyperspectral.exe
	rm -f examples/example_mineral_detection examples/example_mineral_detection.exe
	rm -f examples/example_vegetation_analysis examples/example_vegetation_analysis.exe
	rm -f examples/example_spectral_unmixing examples/example_spectral_unmixing.exe
