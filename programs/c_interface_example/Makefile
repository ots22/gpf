# should contain paths to NLOPT and GPF
INC = ../../include   /home/raid/ots22/include
LIB = ../../lib   /home/raid/ots22/lib
CFLAGS = -Wall -Wextra -O3 $(addprefix -I,$(INC))
LDFLAGS = -lgpf -lnlopt -llapack -lblas -lgfortran -lm $(addprefix -L,$(LIB))

c_interface_example: c_interface_example.c
	$(CC) $(CFLAGS) -o c_interface_example c_interface_example.c $(LDFLAGS)
