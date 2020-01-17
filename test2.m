
%ionic liquid dielectric constant
K_IL = 15;
%nanocapacitor gap size
d_nc = 5e-9; %cm
%width of channel
W = 10e-6; %cm
%length of the channel
L = 100e-6; %cm


%elementary charge
e =  1.602e-19;
eV = abs(e);
epsilon_0 = 8.854e-12; %using cm
hbar = 1.05457182e-34; %using cm
m_e = 9.11e-31;
%dielectric constant of STO
K_STO = 303;

Vg = 3;
Vsd = 0.5;
Ef = 0.006 * eV;
ml = 1.2 .* m_e;
gl = 4;
mh = 4.8 .* m_e;
gh = 2;

Constants = [K_IL d_nc W L e epsilon_0 hbar m_e K_STO ml gl mh gh];


%sheet carrier density in m-2 position dependent 
n_2d =  epsilon_0 .* K_IL .* Vg ./ (abs(e) .* d_nc);
% n_2d = 3e13;

z = linspace(0,5e-8,100);

A = 4.097e-5;
B = 4.907e-10;
F_av = (A / B) .* (exp( B .* e .* n_2d ./ (2 .* epsilon_0)) - 1);

E_i_l =@(i) (hbar.^2 / (2.*ml)).^(1/3) .* (3 .* pi .* e .* F_av .* ( i - 0.25) / 2) .^ (2 / 3);
E_i_h =@(i) (hbar.^2 / (2.*mh)).^(1/3) .* (3 .* pi .* e .* F_av .* ( i - 0.25) / 2) .^ (2 / 3);

psi_i_l =@(i) airy(0, z .* (2 .* ml .* e .* F_av ./ hbar .^ 2) .^ (1/3) - (2 * ml * E_i_l(i) ./ hbar .^2) .* ...
    (hbar .^2 / (2 .* ml .* e .* F_av)) .^ (2 / 3));
psi_i_h =@(i) airy(0, z .* (2 .* mh .* e .* F_av ./ hbar .^ 2) .^ (1/3) - (2 * mh * E_i_h(i) ./ hbar .^2) .* ...
    (hbar .^2 / (2 .* mh .* e .* F_av)) .^ (2 / 3));

limit_l = floor((2./(3.*pi.*e.*F_av))*((2.*ml/(hbar.^2)).^(1/3) * Ef).^(2/3) + 0.25);
limit_h = floor((2./(3.*pi.*e.*F_av))*((2.*mh/(hbar.^2)).^(1/3) * Ef).^(2/3) + 0.25);
floor((2/(3*pi)).*(Ef./(e .* F_av .* z_o_h)).^(3/2) + 1/4 ); 

suml = 0;
for i = 1:limit_l
    suml = suml + vpa((Ef - E_i_l(i)).*abs(psi_i_l(i)).^2);
end
sumh = 0;
for i = 1:limit_h
     sumh = sumh + vpa((Ef - E_i_h(i)).*abs(psi_i_h(i)).^2);
end

n_3d =(gl.*ml/(2 .* pi .* hbar.^2)) .* suml + (gh .* mh ./ (2.* pi .* hbar.^2)) .* sumh;

plot(z, n_3d)
