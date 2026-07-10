//
//  Interactions for IBI
//  
//
//  Created by Haisu Kang 10/15/23
//
//
//Spring force
for(i = 0; i<nchain1; ++i)
{
	for (j = 0; j<Nisc1-1; j++)
	{
		dx = rxa[i*Nisc1+j] - rxa[i*Nisc1+j+1];
		dy = rya[i*Nisc1+j] - rya[i*Nisc1+j+1];
		dz = rza[i*Nisc1+j] - rza[i*Nisc1+j+1];
		if(dx>Lx/2) dx -= Lx; if(dx<-Lx/2) dx += Lx;
		if(dy>Ly/2) dy -= Ly; if(dy<-Ly/2) dy += Ly;
		if(dz>Lz/2) dz -= Lz; if(dz<-Lz/2) dz += Lz;
		rr = dx*dx+dy*dy+dz*dz;
		r = sqrt(rr);
		Fs = -kappa*(r-rcut1);
		fxa[i*Nisc1+j] += Fs*dx/r;
		fya[i*Nisc1+j] += Fs*dy/r;
		fza[i*Nisc1+j] += Fs*dz/r;
		fxa[i*Nisc1+j+1] += -Fs*dx/r;
		fya[i*Nisc1+j+1] += -Fs*dy/r;
		fza[i*Nisc1+j+1] += -Fs*dz/r;

        U += 0.5*kappa*(r-rcut1)*(r-rcut1);
        F += Fs;  
        bondlen += (double)r/((Nisc1-1)*nchain1);
		//printf("bond %d %d\n",i*Nisc+j, i*Nisc+j+1);
	}
}
for(i=0; i<nchain2; i++)
{
    for (j = 0; j<Nisc2-1; j++)
    {
        dx = rxb[i*Nisc2+j] - rxb[i*Nisc2+j+1];
        dy = ryb[i*Nisc2+j] - ryb[i*Nisc2+j+1];
        dz = rzb[i*Nisc2+j] - rzb[i*Nisc2+j+1];
        if(dx>Lx/2) dx -= Lx; if(dx<-Lx/2) dx += Lx;
        if(dy>Ly/2) dy -= Ly; if(dy<-Ly/2) dy += Ly;
        if(dz>Lz/2) dz -= Lz; if(dz<-Lz/2) dz += Lz;
        rr = dx*dx+dy*dy+dz*dz;
        r = sqrt(rr);
        Fs = -kappa*(r-rcut2);
        fxb[i*Nisc2+j] += Fs*dx/r;
        fyb[i*Nisc2+j] += Fs*dy/r;
        fzb[i*Nisc2+j] += Fs*dz/r;
        fxb[i*Nisc2+j+1] += -Fs*dx/r;
        fyb[i*Nisc2+j+1] += -Fs*dy/r;
        fzb[i*Nisc2+j+1] += -Fs*dz/r;

        U += 0.5*kappa*(r-rcut2)*(r-rcut2);
        F += Fs;  
        bondlen += (double)r/((Nisc2-1)*nchain2);
        //printf("bond %d %d\n",i*Nisc+j, i*Nisc+j+1);
    }
}
//Bending force
for(i=0; i<nchain1; i++)
{
	for(j=0; j<(Nisc1-2); j++)
	{
		dx1 = rxa[i*Nisc1+j] - rxa[i*Nisc1+j+1];
        dy1 = rya[i*Nisc1+j] - rya[i*Nisc1+j+1];
        dz1 = rza[i*Nisc1+j] - rza[i*Nisc1+j+1];
        dx2 = rxa[i*Nisc1+j+1] - rxa[i*Nisc1+j+2];
        dy2 = rya[i*Nisc1+j+1] - rya[i*Nisc1+j+2];
        dz2 = rza[i*Nisc1+j+1] - rza[i*Nisc1+j+2];
        if(dx1>Lx/2) dx1 -= Lx; if(dx1<-Lx/2) dx1 += Lx;
        if(dy1>Ly/2) dy1 -= Ly; if(dy1<-Ly/2) dy1 += Ly;
        if(dz1>Lz/2) dz1 -= Lz; if(dz1<-Lz/2) dz1 += Lz;
        if(dx2>Lx/2) dx2 -= Lx; if(dx2<-Lx/2) dx2 += Lx;
        if(dy2>Ly/2) dy2 -= Ly; if(dy2<-Ly/2) dy2 += Ly;
        if(dz2>Lz/2) dz2 -= Lz; if(dz2<-Lz/2) dz2 += Lz;
		rr1 = dx1*dx1+dy1*dy1+dz1*dz1;
		rr2 = dx2*dx2+dy2*dy2+dz2*dz2;
		r1 = sqrt(rr1);
		r2 = sqrt(rr2);
		cos_theta = (dx1*dx2+dy1*dy2+dz1*dz2)/(r1*r2);
		Fb = kappa_bend1*(1-cos_theta);
        fxa[i*Nisc1+j] += Fb*(dx2/r2-cos_theta*dx1/r1);
        fya[i*Nisc1+j] += Fb*(dy2/r2-cos_theta*dy1/r1);
        fza[i*Nisc1+j] += Fb*(dz2/r2-cos_theta*dz1/r1);
        fxa[i*Nisc1+j+1] += Fb*(dx1/r1-cos_theta*dx2/r2);
        fya[i*Nisc1+j+1] += Fb*(dy1/r1-cos_theta*dy2/r2);
        fza[i*Nisc1+j+1] += Fb*(dz1/r1-cos_theta*dz2/r2);
        fxa[i*Nisc1+j+2] += Fb*(cos_theta*dx1/r1+cos_theta*dx2/r2);
        fya[i*Nisc1+j+2] += Fb*(cos_theta*dy1/r1+cos_theta*dy2/r2);
        fza[i*Nisc1+j+2] += Fb*(cos_theta*dz1/r1+cos_theta*dz2/r2);
        U += -0.5*kappa_bend2*(1-cos_theta)*(1-cos_theta);
        F += Fb;
	}
}
for(i=0; i<nchain2; i++)
{
    for(j=0; j<(Nisc2-2); j++)
    {
        dx1 = rxb[i*Nisc2+j] - rxb[i*Nisc2+j+1];
        dy1 = ryb[i*Nisc2+j] - ryb[i*Nisc2+j+1];
        dz1 = rzb[i*Nisc2+j] - rzb[i*Nisc2+j+1];
        dx2 = rxb[i*Nisc2+j+1] - rxb[i*Nisc2+j+2];
        dy2 = ryb[i*Nisc2+j+1] - ryb[i*Nisc2+j+2];
        dz2 = rzb[i*Nisc2+j+1] - rzb[i*Nisc2+j+2];
        if(dx1>Lx/2) dx1 -= Lx; if(dx1<-Lx/2) dx1 += Lx;
        if(dy1>Ly/2) dy1 -= Ly; if(dy1<-Ly/2) dy1 += Ly;
        if(dz1>Lz/2) dz1 -= Lz; if(dz1<-Lz/2) dz1 += Lz;
        if(dx2>Lx/2) dx2 -= Lx; if(dx2<-Lx/2) dx2 += Lx;
        if(dy2>Ly/2) dy2 -= Ly; if(dy2<-Ly/2) dy2 += Ly;
        if(dz2>Lz/2) dz2 -= Lz; if(dz2<-Lz/2) dz2 += Lz;
        rr1 = dx1*dx1+dy1*dy1+dz1*dz1;
        rr2 = dx2*dx2+dy2*dy2+dz2*dz2;
        r1 = sqrt(rr1);
        r2 = sqrt(rr2);
        cos_theta = (dx1*dx2+dy1*dy2+dz1*dz2)/(r1*r2);
        Fb = kappa_bend2*(1-cos_theta);
        fxb[i*Nisc2+j] += Fb*(dx2/r2-cos_theta*dx1/r1);
        fyb[i*Nisc2+j] += Fb*(dy2/r2-cos_theta*dy1/r1);
        fzb[i*Nisc2+j] += Fb*(dz2/r2-cos_theta*dz1/r1);
        fxb[i*Nisc2+j+1] += Fb*(dx1/r1-cos_theta*dx2/r2);
        fyb[i*Nisc2+j+1] += Fb*(dy1/r1-cos_theta*dy2/r2);
        fzb[i*Nisc2+j+1] += Fb*(dz1/r1-cos_theta*dz2/r2);
        fxb[i*Nisc2+j+2] += Fb*(cos_theta*dx1/r1+cos_theta*dx2/r2);
        fyb[i*Nisc2+j+2] += Fb*(cos_theta*dy1/r1+cos_theta*dy2/r2);
        fzb[i*Nisc2+j+2] += Fb*(cos_theta*dz1/r1+cos_theta*dz2/r2);
        U += -0.5*kappa_bend2*(1-cos_theta)*(1-cos_theta);
        F += Fb;
    }
}
//PMF force from table
for(i = 0; i<nchain1*Nisc1-1; i++)
{
	for(j = i+1; j<nchain1*Nisc1; j++)
	{
		// AA pair
        dx = rxa[i] - rxa[j];
        dy = rya[i] - rya[j];
        dz = rza[i] - rza[j];
        if(dx>Lx/2) dx -= Lx; if(dx<-Lx/2) dx += Lx;
        if(dy>Ly/2) dy -= Ly; if(dy<-Ly/2) dy += Ly;
        if(dz>Lz/2) dz -= Lz; if(dz<-Lz/2) dz += Lz;
        rr = dx*dx+dy*dy+dz*dz;
        r = sqrt(rr);
        if(r<rp[dmax-1])// When r is within the table
        {
            for(k=0; k<dmax; k++)
            {
                if(r - rp[k] <= 0.01)// When the pair distance can be found in the table within the error 0.01
                {
                    low_index=k;
                    break;
                }
                else low_index = -1;// when the pair distance cannot be found in the table
            }

            if(low_index != -1)//Use table value
            {
                coeff = mfo_aa[low_index];
                U += pmfo_aa[low_index];
                Unb += pmfo_aa[low_index];
            }
            else if(r > rp[0] && r < Lx)//if r is not in table, linear interpolation
            {
                low_index = (int)r/binsize;
                high_index = low_index +1;
                if(r >= rp[low_index] && r <= rp[high_index])
                {
                    coeff = -(mfo_aa[high_index]-mfo_aa[low_index])/(high_index-low_index);					
                }
                U += pmfo_aa[low_index] + (r-rp[low_index])*(pmfo_aa[high_index]-pmfo_aa[low_index])/(double)(rp[high_index]-rp[low_index]);
                Unb += pmfo_aa[low_index] + (r-rp[low_index])*(pmfo_aa[high_index]-pmfo_aa[low_index])/(double)(rp[high_index]-rp[low_index]);
            }
            else if(r > Lx) // if r is larger than maximum cutoff, 0 
            {
                coeff = 0; 
                U+= 0;
            }	
            else if(r <= rp[0])// if ri is smaller than minimum cutoff, use minimum value
            {
                coeff = mfo_aa[0];
                U += pmfo_aa[0];
                Unb += pmfo_aa[0];
            }
            fxa[i] += coeff*dx/r;
            fya[i] += coeff*dy/r;
            fza[i] += coeff*dz/r;
            fxa[j] += -coeff*dx/r;
            fya[j] += -coeff*dy/r;
            fza[j] += -coeff*dz/r;

            F += coeff;
			
		}
    }
}

for(i = 0; i<nchain2*Nisc2-1; i++)
{
    for(j = i+1; j<nchain2*Nisc2; j++)
    {
       //BB pair
       dx = rxb[i] - rxb[j];
       dy = ryb[i] - ryb[j];
       dz = rzb[i] - rzb[j];
       if(dx>Lx/2) dx -= Lx; if(dx<-Lx/2) dx += Lx;
       if(dy>Ly/2) dy -= Ly; if(dy<-Ly/2) dy += Ly;
       if(dz>Lz/2) dz -= Lz; if(dz<-Lz/2) dz += Lz;
       rr = dx*dx+dy*dy+dz*dz;
       r = sqrt(rr);
       if(r<rp[dmax-1])// When r is within the table
       {
           for(k=0; k<dmax; k++)
           {
               if(r - rp[k] <= 0.01)// When the pair distance can be found in the table within the error 0.01
               {
                   low_index=k;
                   break;
               }
               else low_index = -1;// when the pair distance cannot be found in the table
           }

           if(low_index != -1)//Use table value
           {
               coeff = mfo_bb[low_index];
               U += pmfo_bb[low_index];
               Unb += pmfo_bb[low_index];
           }
           else if(r > rp[0] && r < Lx)//if r is not in table, linear interpolation
           {
               low_index = (int)r/binsize;
               high_index = low_index +1;
               if(r >= rp[low_index] && r <= rp[high_index])
               {
                   coeff = -(mfo_bb[high_index]-mfo_bb[low_index])/(high_index-low_index);					
               }
               U += pmfo_bb[low_index] + (r-rp[low_index])*(pmfo_bb[high_index]-pmfo_bb[low_index])/(double)(rp[high_index]-rp[low_index]);
               Unb += pmfo_bb[low_index] + (r-rp[low_index])*(pmfo_bb[high_index]-pmfo_bb[low_index])/(double)(rp[high_index]-rp[low_index]);
           }
           else if(r > Lx) // if r is larger than maximum cutoff, 0 
           {
               coeff = 0; 
               U+= 0;
           }	
           else if(r <= rp[0])// if ri is smaller than minimum cutoff, use minimum value
           {
               coeff = mfo_bb[0];
               U += pmfo_bb[0];
               Unb += pmfo_bb[0];
           }
           fxb[i] += coeff*dx/r;
           fyb[i] += coeff*dy/r;
           fzb[i] += coeff*dz/r;
           fxb[j] += -coeff*dx/r;
           fyb[j] += -coeff*dy/r;
           fzb[j] += -coeff*dz/r;

           F += coeff;
        }
    }
}  

for(i = 0; i<nchain1*Nisc1; i++)
{
    for(j = 0; j<nchain2*Nisc2; j++)
    {
        dx = rxa[i] - rxb[j];
        dy = rya[i] - ryb[j];
        dz = rza[i] - rzb[j];
        if(dx>Lx/2) dx -= Lx; if(dx<-Lx/2) dx += Lx;
        if(dy>Ly/2) dy -= Ly; if(dy<-Ly/2) dy += Ly;
        if(dz>Lz/2) dz -= Lz; if(dz<-Lz/2) dz += Lz;
        rr = dx*dx+dy*dy+dz*dz;
        r = sqrt(rr);
        if(r<rp[dmax-1])// When r is within the table
        {
            for(k=0; k<dmax; k++)
            {
                if(r - rp[k] <= 0.01)// When the pair distance can be found in the table within the error 0.01
                {
                    low_index=k;
                    break;
                }
                else low_index = -1;// when the pair distance cannot be found in the table
            }

            if(low_index != -1)//Use table value
            {
                coeff = mfi[low_index];
                U += pmfi[low_index];
                Unb += pmfi[low_index];
            }
            else if(r > rp[0] && r < Lx)//if r is not in table, linear interpolation
            {
                low_index = (int)r/binsize;
                high_index = low_index +1;
                if(r >= rp[low_index] && r <= rp[high_index])
                {
                    coeff = -(mfi[high_index]-mfi[low_index])/(high_index-low_index);					
                }
                U += pmfi[low_index] + (r-rp[low_index])*(pmfi[high_index]-pmfi[low_index])/(double)(rp[high_index]-rp[low_index]);
                Unb += pmfi[low_index] + (r-rp[low_index])*(pmfi[high_index]-pmfi[low_index])/(double)(rp[high_index]-rp[low_index]);
            }
            else if(r > Lx) // if r is larger than maximum cutoff, 0 
            {
                coeff = 0; 
                U+= 0;
            }	
            else if(r <= rp[0])// if ri is smaller than minimum cutoff, use minimum value
            {
                coeff = mfo_ab[0];
                U += pmfo_ab[0];
                Unb += pmfo_ab[0];
            }
            fxa[i] += coeff*dx/r;
            fya[i] += coeff*dy/r;
            fza[i] += coeff*dz/r;
            fxb[j] += -coeff*dx/r;
            fyb[j] += -coeff*dy/r;
            fzb[j] += -coeff*dz/r;

            F += coeff;
        }
    }
}
