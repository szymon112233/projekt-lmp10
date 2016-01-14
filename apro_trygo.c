#include "makespl.h"

#include <math.h>

void
make_spl(points_t * pts, spline_t * spl)
{

    double         *x = pts->x;
    double         *y = pts->y;
    int		i, j, k, m n=pts->n;
    double a0=0.0,a[n],b[n]

    m = (n-1)/2;

    for(i=0; i<n ; i++)
    {
        a0+=y[i];
    }
    a0/=n;

    for(i=0; i<m ; i++)
    {
        a[i]=0.0;
        for(j=0; j<n ; j++)
        {
            a[i]+=y[j]*cos(2*M_PI*i*j/n);
        }
        a[i]*=2;
        a[i]/=n;

        b[i]=0.0;
        for(j=0; j<n ; j++)
        {
            b[i]+=y[j]*sin(2*M_PI*i*j/n);
        }
        b[i]*=2;
        b[i]/=n;
    }




    if (alloc_spl(spl, nb) == 0)
    {
        for(i = 0 ; i < n ; i++)
        {
            double xx = spl->x[i] = p + i*(k-p)/(m-1);
            spl->f[i] = 0;
            spl->f1[i] = 0;
            spl->f2[i] = 0;
            spl->f3[i] = 0;
            for (j = 0 ; j < m ; j++)
            {
                dt = 2.0*M_PI*j/n;
                spl->f[i] += a[j]*cos(dt*xx) + b[j]*sin(dt*xx);
                spl->f1[i] += (-1)*dt*a[j]*sin(dt*xx) + dt*b[j]*cos(dt*xx);
                spl->f2[i] += (-1)*dt*dt*a[j]*cos(dt*xx) + (-1)*dt*dt*b[j]*sin(dt*xx);
                spl->f3[i] += pow(dt,3)*a[j]*sin(dt*xx) + (-1)*pow(dt,3)*b[j]*cos(dt*xx);
            }
        }
    }
}

