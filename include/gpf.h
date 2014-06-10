#ifdef __cplusplus
extern "C" {
#endif

typedef void *gpf_handle;
void gp_projected_process_read(gpf_handle *handle, int fname_len, char *fname);
void gp_projected_process_predict(double *output, gpf_handle *handle, int input_len, double *input, int input_type);

#ifdef __cplusplus
}
#endif
