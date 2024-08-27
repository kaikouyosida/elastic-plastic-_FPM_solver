CC = gcc
CFLAGS = -Ofast -pipe -Wall -g
.PHONY: solver
TARGET = solver
SRCS = main.c model.c
SRCS += fpm.c field.c ss_curve.c
SRCS += internal_force.c d_matrix.c b_matrix.c stress.c
SRCS += external_force.c
SRCS += coefficient_matrix.c s_matrix.c
SRCS += matrix.c vector.c tensor.c scalar.c GetGaussPoints.c ImposeDirichretCondition.c LU_decomposition.c
SRCS += Output.c Output_data.c
OBJS = $(SRCS:.c=.o)

run:
	make -C Data_Files_Input
	make all
	make clean
	
all: $(TARGET)
	.\solver.exe

$(TARGET): $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

$(OBJS): $(SRCS)
	$(CC) -c $(SRCS)

clean:
	-rm -f $(OBJS) $(TARGET)