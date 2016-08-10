
#include <cppcolleff/Bunch.h>

void Bunch_t::InsertionSort()
{
	for (long i=1; i < xx.size(); i++){
		double atualxx=xx[i]; double atualxl=xl[i]; double atualde=de[i]; double atualss=ss[i];
        long j = i - 1;
		while ((j>=0) && (atualss < ss[j])){
			ss[j+1]=ss[j];   xx[j+1]=xx[j];    xl[j+1]=xl[j];    de[j+1]=de[j];
            --j;
		}
		ss[j+1] = atualss;    xx[j+1] = atualxx;    xl[j+1] = atualxl;    de[j+1] = atualde;
	}
}
