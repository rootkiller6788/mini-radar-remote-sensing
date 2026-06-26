# Makefile for mini-infrared-thermal
# IR Thermal Imaging and Remote Sensing Module
#
# Targets:
#   make all      - Build library (libirthermal.a)
#   make test     - Build and run tests
#   make examples - Build examples
#   make clean    - Remove build artifacts

CC       = gcc
CFLAGS   = -std=c99 -Wall -Wextra -O2 -Iinclude
AR       = ar
ARFLAGS  = rcs

# Source files
SRC = src/ir_core.c       src/ir_detector.c       src/ir_radiometry.c       src/ir_calibration.c       src/ir_image.c       src/ir_detection.c       src/ir_atmosphere.c       src/ir_enhancement.c

OBJ = $(SRC:.c=.o)
LIB = libirthermal.a

# Test sources
TEST_SRC = tests/test_ir_core.c            tests/test_ir_detector.c            tests/test_ir_radiometry.c            tests/test_ir_calibration.c

TEST_BIN = $(TEST_SRC:.c=.exe)

# Example sources
EXAMPLE_SRC = examples/example_blackbody.c               examples/example_fever_screening.c               examples/example_building_inspection.c

EXAMPLE_BIN = $(EXAMPLE_SRC:.c=.exe)

.PHONY: all test examples clean

all: $(LIB)

$(LIB): $(OBJ)
	$(AR) $(ARFLAGS) $@ $^

src/%.o: src/%.c include/*.h
	$(CC) $(CFLAGS) -c $< -o $@

test: $(TEST_BIN)
	@echo "Running tests..."
	@for t in $(TEST_BIN); do ./$$t && echo "  PASS: $$t" || echo "  FAIL: $$t"; done

tests/%.exe: tests/%.c $(LIB)
	$(CC) $(CFLAGS) $< -L. -lirthermal -lm -o $@

examples: $(EXAMPLE_BIN)

examples/%.exe: examples/%.c $(LIB)
	$(CC) $(CFLAGS) $< -L. -lirthermal -lm -o $@

clean:
	rm -f $(OBJ) $(LIB) $(TEST_BIN) $(EXAMPLE_BIN)
	rm -f src/*.o tests/*.o examples/*.o

# Help target
help:
	@echo "mini-infrared-thermal Makefile"
	@echo "  make all      - Build static library"
	@echo "  make test     - Build and run tests"
	@echo "  make examples - Build example programs"
	@echo "  make clean    - Remove build artifacts"
