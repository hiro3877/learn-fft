#include <stdio.h>
#define n 8

int main(){
  double temp;
  double xr[n];
  int i,j,n2,k;

  n2=n/2;

  for (i=0;i<n;i++){
  xr[i]=i;
  printf("%lf\n",xr[i]);
}



  //reverse bit

  j=n2;

  for(i=1;i<=n-3;i++){
    if(i<j){
      temp =xr[i];
      xr[i]=xr[j];
      xr[j]=temp;
    }
    k=n2;
    while(k<=j){
      j -= k;
      k /=2;
    }
    j+=k;
  }

for (i=0;i<n;i++){
printf("%lf\n",xr[i]);
}

}
