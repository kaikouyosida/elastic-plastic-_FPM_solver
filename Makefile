CC = icx
CFLAGS = -Ofast -Wall  -g  -I$(INCLUDE_DIR) -mkl
CFLAGS += -Wno-unused-variable -Wno-unused-value -Wno-deprecated-declarations -Wno-uninitialized -Wno-format-extra-args -Wno-misleading-indentation
# インクルードディレクトリ
INCLUDE_DIR = "C:/Program Files (x86)/Intel/oneAPI/mkl/latest/include"

# ライブラリディレクトリとファイル
LIB_DIR_MKL = C:/Program Files (x86)/Intel/oneAPI/mkl/2024.2/lib
LIB_DIR_COMPILER = C:/Program Files (x86)/Intel/oneAPI/compiler/latest/lib
LIB_DIR_DECIMAL = C:/Program Files (x86)/Intel/oneAPI/2024.2/lib

LIBS = "$(LIB_DIR_MKL)/mkl_rt.lib" \
       "$(LIB_DIR_COMPILER)/libircmt.lib" \
       "$(LIB_DIR_COMPILER)/svml_dispmt.lib" \
       "$(LIB_DIR_COMPILER)/libmmt.lib" \
	   "$(LIB_DIR_COMPILER)/libirc.lib" \
       "$(LIB_DIR_DECIMAL)/libdecimal.lib"

.PHONY: solver
TARGET = solver
SRCS = main.c model.c
SRCS += fpm.c field.c ss_curve.c
SRCS += internal_force.c d_matrix.c b_matrix.c stress.c
SRCS += external_force.c
SRCS += coefficient_matrix.c s_matrix.c
SRCS += matrix.c vector.c tensor.c scalar.c GetGaussPoints.c ImposeDirichretCondition.c LU_decomposition.c MKL_solver.c
SRCS += Output.c Output_data.c
OBJS = $(SRCS:.c=.o)

run:
	make -C Data_Files_Input
	make all
	make clean
	
all: $(TARGET)
	.\solver.exe

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS) $(LINK_FLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean1:
	-rm -f $(OBJS) $(TARGET)
	
clean2:
	make clean -C Data_Files_Input
	-rm -f $(OBJS) $(TARGET)