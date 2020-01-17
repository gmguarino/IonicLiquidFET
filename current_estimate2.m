
%ionic liquid dielectric constant
K_IL = 15;
%nanocapacitor gap size
d_nc = 5e-9;
%width of channel
W = 15e-6;
%length of the channel
L= 200e-6;
% Ueno dimensions


%elementary charge
e = - 1.602e-19;
epsilon_0 = 8.854e-12;
hbar = 1.05457182e-34;
m_e = 9.11e-31;
%dielectric constant of STO
K_STO = 300;
%Capacitance of EDL
Ci = epsilon_0 * K_IL * W * L / d_nc;
%sheet carrier density in m-2
n_2d =@(V)  epsilon_0 .* K_IL .* V / (abs(e) .* d_nc);
%Threshold voltage of EDLT according to Ueno
V_th = 1.5;
%average electric field induced in the STO layer
F_av =@(V) e .* n_2d(V) /(2 .* epsilon_0 .* K_STO);

z_o =@(m,V) (hbar.^2 ./ (2 .* m  .* e .* F_av(V))).^(1/3);

E_i =@(m,i, V) e .* F_av(V) .* z_o(m,V) .* (3 .* pi .* (i - 0.25) ./ 2).^(2/3);

syms z k V E

psi_i =@(z,m,i, V) airy(z/z_o(m,V) - E_i(i,m, V)/(e .* F_av(V) .* z_o(m,V)));

ml = 1.2 .* m_e;
gl = 4;
mh = 4.8 .* m_e;
gh = 2;

Constants = [K_IL d_nc W L e epsilon_0 hbar m_e K_STO ml gl mh gh];
%width of channel
W = 10e-6;
%length of the channel
L= 100e-6;


%elementary charge
e = - 1.602e-19;
epsilon_0 = 8.854e-12;
hbar = 1.05457182e-34;
m_e = 9.11e-31;
%dielectric constant of STO
K_STO = 300;

mu =20e-4;

%Defining the limit of the summation
limit = @(Ef, m, V) floor((2/(3*pi)).*(Ef/(e * F_av(V) * z_o(m,V))).^(3/2) + 1/4 ); 

n_3d =@(z, Ef,V) (gl.*ml/(2 .* pi .* hbar.^2)) .* sum((Ef - E_i(ml,1:limit(Ef,ml,V),V)).* abs(psi_i(z,ml,1:limit(Ef,ml,V),V)).^2) + ...
    (gh .* mh / (2.* pi .* hbar.^2)) .* sum((Ef - E_i(mh,1:limit(Ef,ml,V),V)).* abs(psi_i(z,mh,1:limit(Ef,ml,V),V)).^2);

%Mod squared integral of the wavefunctions
int_psi_sq = @(k,E, m, V) E_i(k,m, V)/(e .* F_av(V) .* z_o(m,V)).* airy(0,-E_i(k,m, V)/(e .* F_av(V) .* z_o(m,V))).^2 + airy(1,-E_i(k,m, V)/(e .* F_av(V) .* z_o(m,V))).^2;

Ids =@(Vg, Vds) piecewise(Vds <= Vg - V_th, (W .* mu .* Ci / L) .* (Vg - V_th - Vds./2) .* Vds, Vds > Vg - V_th, (W .* mu .* Ci / L) * (Vg-V_th).^2./2);
figure (1)
hold on
for Vsd=0:0.2:0.6
fplot(Ids(V,Vsd),[1.5, 5])
end



function [K_IL, d_nc, W, L, e, epsilon_0, hbar, m_e, K_STO, ml, gl, mh, gh] = unpack(c)
constant = num2cell(c); 
[K_IL, d_nc, W, L, e, epsilon_0, hbar, m_e, K_STO, ml, gl, mh, gh] = constant{:};

end

