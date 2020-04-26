/* util file
 *
 * Some functions used for coverting data between matlab and c++
 *
 */
#ifndef __UTILS_H__
#define __UTILS_H__

template<typename T, typename Index>
void csr2csc(
  Index n, Index m, Index nz, 
  T *a, Index *col_idx, Index *row_start,
  T *csc_a, Index *row_idx, Index *col_start)
{
  Index i, j, k, l;
  Index *ptr;

  for (i=0; i<=m; i++) col_start[i] = 0;

  /* determine column lengths */
  for (i=0; i<nz; i++) col_start[col_idx[i]+1]++;
  for (i=0; i<m; i++) col_start[i+1] += col_start[i];

  /* go through the structure once more. Fill in output matrix. */
  for (i=0, ptr=row_start; i<n; i++, ptr++)
    for (j=*ptr; j<*(ptr+1); j++){
      k = col_idx[j];
      l = col_start[k]++;
      row_idx[l] = i;
      if (a) csc_a[l] = a[j];
    }

  /* shift back col_start */
  for (i=m; i>0; i--) col_start[i] = col_start[i-1];

  col_start[0] = 0;
}

// Do CHECK and throw a Mex error if check fails
inline void mxCHECK(bool expr, const char* msg) {
  if (!expr) {
    mexErrMsgTxt(msg);
  }
}

static double init_key = static_cast<double>(rand());

// Create a handle struct vector, without setting up each handle in it
template <typename T>
static mxArray* create_handle_vec(int ptr_num) {
  const int handle_field_num = 2;
  const char* handle_fields[handle_field_num] = { "ptr", "init_key" };
  return mxCreateStructMatrix(ptr_num, 1, handle_field_num, handle_fields);
}

// Set up a handle in a handle struct vector by its index
template <typename T>
static void setup_handle(const T* ptr, int index, mxArray* mx_handle_vec) {
  mxArray* mx_ptr = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
  *reinterpret_cast<uint64_t*>(mxGetData(mx_ptr)) =
      reinterpret_cast<uint64_t>(ptr);
  mxSetField(mx_handle_vec, index, "ptr", mx_ptr);
  mxSetField(mx_handle_vec, index, "init_key", mxCreateDoubleScalar(init_key));
}

// Convert a pointer in C++ to a handle in matlab
template <typename T>
static mxArray* ptr_to_handle(const T* ptr) {
  mxArray* mx_handle = create_handle_vec<T>(1);
  setup_handle(ptr, 0, mx_handle);
  return mx_handle;
}

template <typename T>
static T* handle_to_ptr(const mxArray* mx_handle) {
  mxArray* mx_ptr = mxGetField(mx_handle, 0, "ptr");
  mxArray* mx_init_key = mxGetField(mx_handle, 0, "init_key");
  mxCHECK(mxIsUint64(mx_ptr), "pointer type must be uint64");
  mxCHECK(mxGetScalar(mx_init_key) == init_key,
      "Could not convert handle to pointer due to invalid init_key. "
      "The object might have been cleared.");
  return reinterpret_cast<T*>(*reinterpret_cast<uint64_t*>(mxGetData(mx_ptr)));
}

#endif
