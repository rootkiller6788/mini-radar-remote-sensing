# Makefile — mini-radar-basics
# Radar Basics: radar range equation, waveform design, detection, Doppler, signal model
#
# Targets:
#   make          — build library and tests
#   make test     — run test suite
#   make examples — build all examples
#   make all      — build everything (lib, test, examples, demo, bench)
#   make clean    — remove build artifacts

CC       = gcc
CFLAGS   = -std=c11 -Wall -Wextra -O2 -g
INCLUDES = -Iinclude
LDFLAGS  = -lm

SRC_DIR  = src
TST_DIR  = tests
EX_DIR   = examples
DM_DIR   = demos
BN_DIR   = benches
BLD_DIR  = build

SRCS = $(SRC_DIR)/radar_core.c \
       $(SRC_DIR)/radar_waveform.c \
       $(SRC_DIR)/radar_detection.c \
       $(SRC_DIR)/radar_doppler.c \
       $(SRC_DIR)/radar_signal_model.c \
       $(SRC_DIR)/radar_array.c

OBJS = $(SRCS:$(SRC_DIR)/%.c=$(BLD_DIR)/%.o)
LIB  = $(BLD_DIR)/libradar.a

TEST_SRC = $(TST_DIR)/test_radar.c
TEST_BIN = $(BLD_DIR)/test_radar

EX1_SRC = $(EX_DIR)/ex1_pulse_radar.c
EX2_SRC = $(EX_DIR)/ex2_fmcw_radar.c
EX3_SRC = $(EX_DIR)/ex3_cfar_detection.c
EX1_BIN = $(BLD_DIR)/ex1_pulse_radar
EX2_BIN = $(BLD_DIR)/ex2_fmcw_radar
EX3_BIN = $(BLD_DIR)/ex3_cfar_detection

DEMO_SRC = $(DM_DIR)/demo_radar.c
DEMO_BIN = $(BLD_DIR)/demo_radar
BENCH_SRC = $(BN_DIR)/bench_radar.c
BENCH_BIN = $(BLD_DIR)/bench_radar

.PHONY: all test examples demo bench clean

all: $(LIB) $(TEST_BIN) $(EX1_BIN) $(EX2_BIN) $(EX3_BIN) $(DEMO_BIN) $(BENCH_BIN)

$(BLD_DIR):
	mkdir -p $(BLD_DIR)

$(LIB): $(BLD_DIR) $(OBJS)
	ar rcs $@ $(OBJS)

$(BLD_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

$(TEST_BIN): $(LIB) $(TEST_SRC)
	$(CC) $(CFLAGS) $(INCLUDES) $(TEST_SRC) -L$(BLD_DIR) -lradar $(LDFLAGS) -o $@

$(EX1_BIN): $(LIB) $(EX1_SRC)
	$(CC) $(CFLAGS) $(INCLUDES) $(EX1_SRC) -L$(BLD_DIR) -lradar $(LDFLAGS) -o $@

$(EX2_BIN): $(LIB) $(EX2_SRC)
	$(CC) $(CFLAGS) $(INCLUDES) $(EX2_SRC) -L$(BLD_DIR) -lradar $(LDFLAGS) -o $@

$(EX3_BIN): $(LIB) $(EX3_SRC)
	$(CC) $(CFLAGS) $(INCLUDES) $(EX3_SRC) -L$(BLD_DIR) -lradar $(LDFLAGS) -o $@

$(DEMO_BIN): $(LIB) $(DEMO_SRC)
	$(CC) $(CFLAGS) $(INCLUDES) $(DEMO_SRC) -L$(BLD_DIR) -lradar $(LDFLAGS) -o $@

$(BENCH_BIN): $(LIB) $(BENCH_SRC)
	$(CC) $(CFLAGS) $(INCLUDES) $(BENCH_SRC) -L$(BLD_DIR) -lradar $(LDFLAGS) -o $@

test: $(TEST_BIN)
	./$(TEST_BIN)

examples: $(EX1_BIN) $(EX2_BIN) $(EX3_BIN)

demo: $(DEMO_BIN)
	./$(DEMO_BIN)

bench: $(BENCH_BIN)
	./$(BENCH_BIN)

clean:
	rm -rf $(BLD_DIR)
