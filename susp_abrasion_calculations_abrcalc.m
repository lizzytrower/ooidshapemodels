%compute flow velocity
z0 = 3*D/30; %[m] roughness coefficient set by grainsize
dz = (H-z0)/1000; %[m]
z = z0:dz:H; %[m]
Uf = sum((ustar/0.41)*log(z/z0))*dz/H; %[m/s] depth-averaged flow velocity

%compute bed load height and velocity
hb = D*1.44*(tstage-1)^0.5; %[m] height of the bed load layer 
Us = (R*g*D)^0.5*1.56*(tstage-1)^0.56; %[m/s] bed load velocity

if Us>Uf
Us=Uf;
end

%compute suspended load
if hb < H 
    hb(hb<D)=D; %sets minimum ht of bed load layer
    b = hb; %bottom of the suspended load layer - same as height of bedload
    betta = 2; %Based on Scheingross et al. (2014) best fit
    P = ws/(0.41*ustar*betta); %Rouse number
    res = 1000;
    di5 = (log(H)-log(b))/res;
    i5 = log(b):di5:log(H);
    z = exp(i5);
    z(length(z))=H;
    dz = diff(z);
    dz = [dz(1) dz];
    a1=sum((((1-(z(z>z0)./H))/(1-(b/H))).*(b./z(z>z0))).^P.*log(z(z>z0)./z0).*dz(z>z0))...
        /(Uf*H)*(ustar/0.41);
    cb = 1/(Us*hb +Uf*H*a1);

    %find concentration profile
    c=0;
    c(1) = cb;
    c(2:length(z)+1) = cb.*(((1-(z./H))./(1-(b./H))).*(b./z)).^P ;
    z=[0 z];
    c(z==H)=0;

%calculate the fall distance
    gradc(1:length(c))=0;
    gradc(2:length(c)) = -diff(c);    
    Hfall = (1./cb).*sum(z.*gradc);
else
    hb = H;
    cb = 1/(Us.*hb);
    Hfall = hb;   
    a1=0;
end

if cb == 0
    Hfall = 0;
end

sig = ustar;
dx = sig/100; %the number of bins to subdivide - can play with this for resolution
X = -6*sig:dx:6*sig; %spread distribution for six sigma = w'/ws
f = normpdf(X,0,sig); %centered at zero normal gausian distribution
X = X./ws;  % to normalize as w'/ws same as psi

Scos = 1;  %cosine of the angle of the bed for impacts.

wfall = Scos*((2*(2/3)*D*g/cdrag*R)*(1-exp(-cdrag*rho_f/rho_s*...
    (Hfall/Scos)/(2/3*D))))^0.5;

psifall = wfall./ws;
settlematrix = psifall + X;
settlematrix1=settlematrix;
settlematrix(settlematrix<0) = 0; %no negative impacts
psifall_turb = sum((settlematrix).*f).*dx;
psi_fall3 = sum((settlematrix.^3).*f).*dx;
E1 = psi_fall3;  %erosion with turbulence

%Stokes number correction
wi_st = settlematrix;
wi_st((D.*wi_st.*ws.*rho_s./(18.*nu.*rho_f))<Stc) = 0;
psi_fall3_st = sum((wi_st.^3).*f).*dx;
E1_st = psi_fall3_st;  %erosion with turbulence and stokes correction

ti = D./Hfall; %[dimensionless]
ti(Hfall<=0.5*D)=0;

En_suspt_st = 1./kv.*(ws./(g.*D).^0.5).^3.*E1_st.*ti./6; 
En_suspt_st(En_suspt_st<0)=0;

Efactor = 10^6.*60*60*rho_s.*young.*(g.*D).^(3/2)./(strength.^2);

if tstage <=1
    En_suspt_st = 0;
end

if H < D
    En_suspt_st = 0;
end

Rabrasion = En_suspt_st*Efactor*A1; %[um/hr]
