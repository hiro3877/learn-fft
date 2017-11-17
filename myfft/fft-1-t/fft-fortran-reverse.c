#include <stdio.h>
#define N 8

int main(){
  double W;
  double F[N];
  int K,L,N2,I;

  N2=N/2;

  for (K=0;K<N;K++){
  F[K]=K;
  printf("%lf\n",F[K]);
}



  //reverse bit

  L=N2;

  for(K=1;K<=N-2;K++){
    if(K<L){
      W=F[L+1-1];
      F[L+1-1]=F[K+1-1];
      F[K+1-1]=W;
    }
    I=N2;
    while(I<L+1){
      L=L-I;
      I=I/2;
    }
    L=L+I;
  }

for (K=0;K<N;K++){
printf("%lf\n",F[K]);
}

}
