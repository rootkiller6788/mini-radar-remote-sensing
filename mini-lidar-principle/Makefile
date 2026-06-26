# Makefile — mini-lidar-principle
CC       = gcc
CFLAGS   = -std=c11 -Wall -Wextra -O2 -g -Iinclude -DM_PI=3.14159265358979323846
LDFLAGS  = -lm

LIB_SRCS = src/lidar_core.c src/lidar_waveform.c src/lidar_geometry.c src/lidar_scanning.c src/lidar_detection.c src/lidar_registration.c src/lidar_applications.c
LIB_OBJS = src/lidar_core.o src/lidar_waveform.o src/lidar_geometry.o src/lidar_scanning.o src/lidar_detection.o src/lidar_registration.o src/lidar_applications.o

TEST_SRC = tests/test_lidar.c
TEST_BIN = tests/test_lidar.exe

.PHONY: all clean test

all: test

src/%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@

$(TEST_BIN): $(LIB_OBJS) $(TEST_SRC)
	$(CC) $(CFLAGS) -o $@ $(TEST_SRC) $(LIB_OBJS) $(LDFLAGS)

test: $(TEST_BIN)
	./$(TEST_BIN)

clean:
	rm -f src/*.o $(TEST_BIN)
