
inline int get_index_x (int nelems, int index ) { 
  if (index == -1)  {
    index = get_global_id(0);
  } else {
    index += get_global_size(0);
  }

  if (index >= nelems) index = -1;

  return index;
}

inline int get_index_y (int nelems, int index ) { 
  if (index == -1)  {
    index = get_global_id(1);
  } else {
    index += get_global_size(1);
  }

  if (index >= nelems) index = -1;

  return index;
}

inline double Gaussian (double x, double y, double sigma) {
  return (1/sqrt(2*M_PI*sigma)) * exp((-pow(x - y,2.0)) / (2*pow(sigma,2.0)));
}

__kernel void MED_AC(__global double *out, __global const  double *x, __global const double *y, double sigma, int ncols, int nrows) {
  int idx = get_index_x (ncols, -1);
  int idy = get_index_y (nrows, -1);
  
  int i;
  double sum = 0;

  while (idy >= 0) {
    while(idx >= 0) {
      sum = 0;

      for (i=0; i<ncols; i++) {
        sum += Gaussian(x[idx + idy*ncols], y[abs(idx-i) + idy*ncols], sigma);
      }


      out[idx + idy*ncols] = (1/(double)ncols) * sum;
      idx = get_index_x (ncols, idx);
    }
    idy = get_index_y (nrows, idy);
  }
}

__kernel void AC(__global double *out, __global const double *x, __global const double *y,  __global const int *marray, double sigma, int msize, int ncols, int nrows) {
  double sum = 0;
  int i = 0;
  int idm = get_index_x(msize, -1);
  int idy = get_index_y(nrows, -1);
  int m = marray[idm];

  while (idy >= 0) {
    while(idm >= 0) {
      sum = 0;

      for (i=m; i < ncols; i++)
        sum += Gaussian (x[i + idy*ncols], y[abs(i-m) + idy*ncols], sigma);


      out[idm + idy*msize] = (1/((double)ncols-m)) * sum;

      idm = get_index_x(msize, idm);
      m = marray[idm];
    }
    idy = get_index_y (nrows, idy);
  }
}

__kernel void CAC(__global double *out, __global const double *x, __global const double *med, int ncols, int nrows) {
  int idx = get_index_x(ncols, -1);
  int idy = get_index_y(nrows, -1);
  
  while(idy >= 0) {
    while(idx >= 0) {
      out[idx + idy*ncols] = x[idx + idy*ncols] - med[idy];


      idx = get_index_x(ncols, idx);
    }
    idy = get_index_y(nrows, idy);
  }
}
