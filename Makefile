CC = gcc
CFLAGS = -Ofast -pipe
.PHONY: solver
TARGET = solver
SRCS = main.c model.c
SRCS += fpm.c field.c ss_curve.c
SRCS += internal_force.c
SRCS += matrix.c
OBJS = $(SRCS:.c=.o)

$(TARGET): $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

$(OBJS): $(SRCS)
	$(CC) -c $(SRCS)

clean:
	-rm -f $(OBJS) $(TARGET)

all: clean $(OBJS) $(TARGET)