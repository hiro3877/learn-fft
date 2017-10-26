#include<stdio.h>
#include<stdlib.h>

#define WID 16
#define HEI 16


int main(void)
{
        FILE *fp;
        fp=fopen("data256.txt","w");
        if(fp==NULL){
          printf("file open error\n");
        }

        int i,buf;

        for (i=0; i<WID*HEI; i++) {
                buf=rand()%2;
                //printf("%d\n",buf);
                fprintf(fp,"%d\n",buf);
        }
        return 0;
}
