

inline int get_index(int nelems, int index ) { 
  if (index == -1)  {
    index = get_global_id(0);
  } else {
    index += get_global_size(0);
  }

  if (index >= nelems) index = -1;

  return index;
}

inline double Gaussian (double x, double y, double sigma) {
  return (1/sqrt(2*M_PI*sigma)) * exp((-pow(x - y,2)) / (2*pow(sigma,2)));
}

__kernel void GAUS (__global double *out, __global const double *x, __global const double *y, double sigma, int N) {

  int id = get_index(N, -1);
  while(id >= 0) {
    out[id] = (1/sqrt(2*M_PI*sigma)) * exp((-1*pow(x[id] - y[id],2)) / (2*pow(sigma,2)));
    //out[id] = exp(1);
    id = get_index(N, id);
  }

}

__kernel void MED_AC(__global double *out, __global  double *x, __global  double *y, double sigma, int N) {
  int id = get_index(N, -1);
  int i;
  double sum = 0;

  while(id >= 0) {
    sum = 0;

    for (i=0; i<N; i++) {
      sum += Gaussian(x[id], y[abs(id-i)], sigma);
    }

    out[id] = (1/(double)N) * sum;
    id = get_index(N, id);
  }
}

__kernel void AC(__global double *out, __global const double *x, __global const double *y,  __global const int *marray, double sigma, int M, int N) {
  double sum = 0;
  int i = 0;
  int idm = get_index(M, -1);
  int m = marray[idm];

  while(idm >= 0) {
    sum = 0;

    for (i=m; i < N; i++)
      sum += Gaussian (x[i], y[abs(i-m)], sigma);

    out[idm] = (1/((double)N-m)) * sum;

    idm = get_index(M, idm);
    m = marray[idm];
  }
}

__kernel void CAC(__global double *out, __global const double *x, double sigma, double med, int N) {
  int id = get_index(N, -1);
  
  while(id >= 0) {
    out[id] = x[id] - med;

    id = get_index(N, id);
  }
}

__kernel void IP_CC(__global double *out, __global  double *x, __global  double *y, double sigma, int N) {
  int id = get_index(N, -1);
  int i;
  double sum = 0;

  while(id >= 0) {
    sum = 0;

    for (i=0; i<N; i++) {
      sum += Gaussian(x[id], y[i], sigma);
    }

    out[id] = (1/(double)N) * sum;
    id = get_index(N, id);
  }
}

__kernel void IPi_CC(__global double *out, __global const double *x, __global  double *y, double sigma, int N) {
  int id = get_index(N, -1);

  while(id >= 0) {
    out[id] = Gaussian (x[id], y[id], sigma);

    id = get_index(N, id);
  }
}
