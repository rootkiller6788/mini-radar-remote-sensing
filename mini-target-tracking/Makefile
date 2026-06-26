# Makefile for mini-target-tracking
# Target Tracking Library — Radar Multi-Target Tracking Algorithms
#
# Reference: Richards, Scheer, Holm "Principles of Modern Radar" (2010)
#            Bar-Shalom "Tracking and Data Fusion" (2011)
#            Blackman & Popoli "Modern Tracking Systems" (1999)
#
# Course coverage: MIT 6.450, Stanford EE359, Michigan EECS 411,
#                  Georgia Tech ECE 6350, TU Munich HF Eng, ETH 227-0455

CC = gcc
CFLAGS = -Wall -Wextra -std=c99 -O2 -I./include
LDFLAGS = -lm

SRCDIR = src
INCDIR = include
TESTDIR = tests
EXAMPLEDIR = examples
BUILDDIR = build

SRCS = $(wildcard $(SRCDIR)/*.c)
OBJS = $(SRCS:$(SRCDIR)/%.c=$(BUILDDIR)/%.o)

TEST_SRCS = $(wildcard $(TESTDIR)/*.c)
TEST_BINS = $(TEST_SRCS:$(TESTDIR)/%.c=$(BUILDDIR)/%)

EXAMPLE_SRCS = $(wildcard $(EXAMPLEDIR)/*.c)
EXAMPLE_BINS = $(EXAMPLE_SRCS:$(EXAMPLEDIR)/%.c=$(BUILDDIR)/%)

.PHONY: all clean test examples

all: $(OBJS)
	@echo "=== mini-target-tracking library built ==="
	@echo "Objects: $(words $(OBJS)) files in build/"

$(BUILDDIR)/%.o: $(SRCDIR)/%.c
	@mkdir -p $(BUILDDIR)
	$(CC) $(CFLAGS) -c $< -o $@

test: $(TEST_BINS)
	@echo "=== Running Tests ==="
	@for t in $(TEST_BINS); do \
		echo "--- $$t ---"; \
		./$$t && echo "PASS: $$t" || echo "FAIL: $$t"; \
	done
	@echo "=== Tests Complete ==="

$(BUILDDIR)/%: $(TESTDIR)/%.c $(OBJS)
	@mkdir -p $(BUILDDIR)
	$(CC) $(CFLAGS) $< $(OBJS) -o $@ $(LDFLAGS)

examples: $(EXAMPLE_BINS)
	@echo "=== Running Examples ==="
	@for e in $(EXAMPLE_BINS); do \
		echo "--- $$e ---"; \
		./$$e || echo "Example $$e returned non-zero"; \
	done
	@echo "=== Examples Complete ==="

$(BUILDDIR)/%: $(EXAMPLEDIR)/%.c $(OBJS)
	@mkdir -p $(BUILDDIR)
	$(CC) $(CFLAGS) $< $(OBJS) -o $@ $(LDFLAGS)

clean:
	rm -rf $(BUILDDIR)
