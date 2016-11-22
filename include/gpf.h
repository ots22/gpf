#ifdef __cplusplus
extern "C" {
#endif

typedef int gpf_handle;
void gp_projected_process_read(gpf_handle *handle, int fname_len, char *fname);
double gp_projected_process_predict(gpf_handle handle, int input_len, double *input, int input_type);
void gp_projected_process_destroy(gpf_handle handle);

#ifdef __cplusplus
}
#endif
