CC = gcc
CFLAGS = -Wall -Wextra -std=c11 -O2 -I./include
LDFLAGS = -lm
LEAN = lake

SRC = $(wildcard src/*.c)
TEST = $(wildcard tests/*.c)
EXAMPLES = $(wildcard examples/*.c)
LEAN_SRC = $(wildcard src/*.lean)

TEST_BINS = $(TEST:.c=.exe)
EXAMPLE_BINS = $(EXAMPLES:.c=.exe)

.PHONY: all clean test examples check-lean lines

all: test examples

test: $(TEST_BINS)
	@echo "=== Running Tests ==="
	@passed=0; total=0; \
	for t in $(TEST_BINS); do \
		./$$t && passed=$$((passed+1)); \
		total=$$((total+1)); \
	done; \
	echo "=== $$passed/$$total test suites passed ==="

examples: $(EXAMPLE_BINS)
	@echo "=== Running Examples ==="
	@for e in $(EXAMPLE_BINS); do ./$$e; done

$(TEST_BINS): %.exe: %.c $(SRC)
	$(CC) $(CFLAGS) -o $@ $< $(filter-out src/pulse_doppler.lean,$(SRC)) $(LDFLAGS)

$(EXAMPLE_BINS): %.exe: %.c $(SRC)
	$(CC) $(CFLAGS) -o $@ $< $(filter-out src/pulse_doppler.lean,$(SRC)) $(LDFLAGS)

check-lean:
	@echo "=== Lean 4 Check ==="
	@if [ -f src/pulse_doppler.lean ]; then \
		echo "src/pulse_doppler.lean exists"; \
	else \
		echo "Lean file missing"; \
	fi

lines:
	@echo "include/ lines: $$(find include -name '*.h' -exec cat {} + | wc -l)"
	@echo "src/ lines:    $$(find src -name '*.c' -exec cat {} + | wc -l)"
	@echo "total include+src: $$(find include src -name '*.h' -o -name '*.c' | xargs cat | wc -l)"

clean:
	rm -f tests/*.exe examples/*.exe
	rm -f *.o src/*.o tests/*.o examples/*.o
