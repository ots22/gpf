#include <stdio.h>
#include <gpf.h>

int main(void)
{
	char fname[] = "square.gp";

	gpf_handle handle;
	gp_projected_process_read(&handle, sizeof(fname), fname);

	double input[1], output;

	int Nout = 50;
	int i;
	for (i=0; i<=Nout; i++) {
		input[0] = 5.0 * (double)i / Nout;
		output = gp_projected_process_predict(handle, 1, input, 0);
		printf("%f %f\n", input[0], output);
	}

	gp_projected_process_destroy(handle);

	return 0;
}
