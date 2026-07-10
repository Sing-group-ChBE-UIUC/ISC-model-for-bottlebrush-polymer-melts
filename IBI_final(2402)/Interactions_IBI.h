//
//  Interactions for IBI
//  
//
//  Created by Haisu Kang 10/15/23
//
//
//Spring force
for(i = 0; i<nchain; ++i)
{
	for (j = 0; j<Nisc-1; j++)
	{
		dx = rx[i*Nisc+j] - rx[i*Nisc+j+1];
		dy = ry[i*Nisc+j] - ry[i*Nisc+j+1];
		dz = rz[i*Nisc+j] - rz[i*Nisc+j+1];
		if(dx>Lx/2) dx -= Lx; if(dx<-Lx/2) dx += Lx;
		if(dy>Ly/2) dy -= Ly; if(dy<-Ly/2) dy += Ly;
		if(dz>Lz/2) dz -= Lz; if(dz<-Lz/2) dz += Lz;
		rr = dx*dx+dy*dy+dz*dz;
		r = sqrt(rr);
		Fs = -kappa*(r-rcut);
		fx[i*Nisc+j] += Fs*dx/r;
		fy[i*Nisc+j] += Fs*dy/r;
		fz[i*Nisc+j] += Fs*dz/r;
		fx[i*Nisc+j+1] += -Fs*dx/r;
		fy[i*Nisc+j+1] += -Fs*dy/r;
		fz[i*Nisc+j+1] += -Fs*dz/r;

        U += 0.5*kappa*(r-rcut)*(r-rcut);
        F += Fs;  
        bondlen += (double)r/((Nisc-1)*nchain);
		//printf("bond %d %d\n",i*Nisc+j, i*Nisc+j+1);
	}
}
//Bending force
for(i=0; i<nchain; i++)
{
	for(j=0; j<(Nisc-2); j++)
	{
		dx1 = rx[i*Nisc+j] - rx[i*Nisc+j+1];
		dy1 = ry[i*Nisc+j] - ry[i*Nisc+j+1];
		dz1 = rz[i*Nisc+j] - rz[i*Nisc+j+1];
		dx2 = rx[i*Nisc+j+1] - rx[i*Nisc+j+2];
		dy2 = ry[i*Nisc+j+1] - ry[i*Nisc+j+2];
		dz2 = rz[i*Nisc+j+1] - rz[i*Nisc+j+2];
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
		Fb = kappa_bend*(1-cos_theta);
		fx[i*Nisc+j] += Fb*(dx2/r2-cos_theta*dx1/r1);
		fy[i*Nisc+j] += Fb*(dy2/r2-cos_theta*dy1/r1);
		fz[i*Nisc+j] += Fb*(dz2/r2-cos_theta*dz1/r1);

		fx[i*Nisc+j+1] += Fb*(dx1/r1-cos_theta*dx2/r2);
		fy[i*Nisc+j+1] += Fb*(dy1/r1-cos_theta*dy2/r2);
		fz[i*Nisc+j+1] += Fb*(dz1/r1-cos_theta*dz2/r2);

		fx[i*Nisc+j+2] += Fb*(cos_theta*dx1/r1+cos_theta*dx2/r2);
		fy[i*Nisc+j+2] += Fb*(cos_theta*dy1/r1+cos_theta*dy2/r2);
		fz[i*Nisc+j+2] += Fb*(cos_theta*dz1/r1+cos_theta*dz2/r2);

		U += -0.5*kappa_bend*(1-cos_theta)*(1-cos_theta);
		F += Fb;
	}
}
//PMF force from table
for(i = 0; i<Ntot-1; i++)
{
	for(j = i+1; j<Ntot; j++)
	{
		dx = rx[i] - rx[j];
		dy = ry[i] - ry[j];
		dz = rz[i] - rz[j];
		if(dx>Lx/2) dx -= Lx; if(dx<-Lx/2) dx += Lx;
		if(dy>Ly/2) dy -= Ly; if(dy<-Ly/2) dy += Ly;
		if(dz>Lz/2) dz -= Lz; if(dz<-Lz/2) dz += Lz;
		rr = dx*dx+dy*dy+dz*dz;
		r = sqrt(rr);
		if(r<rp[dmax-1])
		{
			for(k=0; k<dmax; k++)
			{
				if(r - rp[k] <= 0.01) 
				{
					low_index=k;
					break;
				}
				else low_index = -1;
			}

			if(low_index != -1)//if r is in ro, rp table 10-2
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
				coeff = mfo[0];
				U += pmfo[0];
				Unb += pmfi[0];
			}	
			fx[i] += coeff*dx/r;
			fy[i] += coeff*dy/r;
			fz[i] += coeff*dz/r;
			fx[j] += -coeff*dx/r;
			fy[j] += -coeff*dy/r;
			fz[j] += -coeff*dz/r;

			F += coeff;
		}
	}
}
