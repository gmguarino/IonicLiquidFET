
%ionic liquid dielectric constant
K_IL = 15;
%nanocapacitor gap size
d_nc = 1e-8;
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



ml = 1.2 .* m_e;
gl = 4;
mh = 4.8 .* m_e;
gh = 2;

Constants = [K_IL d_nc W L e epsilon_0 hbar m_e K_STO ml gl mh gh]

syms z k V 

Vg = 3;
n_int(5e-4*abs(e), 3, Constants)
figure(1)
fplot(n_3d(z,0.005 * abs(e), Vg, Constants), [-1e-9, 5e-7])



%E_f = abs(vpasolve(n_int(E,Vg, Constants) == n_2d(Vg, Constants), E, 5e-3 * abs(e)))

function [K_IL, d_nc, W, L, e, epsilon_0, hbar, m_e, K_STO, ml, gl, mh, gh] = unpack(c)
constant = num2cell(c);
[K_IL, d_nc, W, L, e, epsilon_0, hbar, m_e, K_STO, ml, gl, mh, gh] = constant{:};

end

function [n]=n_int(Ef, V, C)
[K_IL, d_nc, W, L, e, epsilon_0, hbar, m_e, K_STO, ml, gl, mh, gh] = unpack(C);
pl = gl.*ml/(2 .* pi .* hbar.^2);
ph = gh .* mh / (2.* pi .* hbar.^2);

Eli = 0;
suml = 0;
for k = 1:limit(Ef, ml, V, C)
    Eli = E_i(ml, k , V,C);
    suml = suml +  Ef - Eli.* Eli/(e .* F_av(V,C) .* z_o(ml,V,C)).* ...
        (airy(-Eli/(e .* F_av(V,C) .* z_o(ml,V,C)))^2 + ...
        airy(1,-Eli/(e .* F_av(V,C) .* z_o(ml,V,C))).^2);
   
end

Ehi = 0;
sumh = 0;
for k = 1:limit(Ef, mh, V, C)
    Ehi = E_i(mh, k , V,C);
    sumh = sumh + Ef - Ehi.* Ehi/(e .* F_av(V,C) .* z_o(mh,V,C)).* ...
        (airy(-Ehi/(e .* F_av(V,C) .* z_o(mh,V,C)))^2 + ...
        airy(1,-Ehi/(e .* F_av(V,C) .* z_o(mh,V,C))).^2);
    
end

n = pl * suml + ph * sumh;
end

function [n3d]=n_3d(z, Ef,V, C)
%sheet carrier density in cm-3, z -dependent
[K_IL, d_nc, W, L, e, epsilon_0, hbar, m_e, K_STO, ml, gl, mh, gh] = unpack(C);

pl = gl.*ml/(2 .* pi .* hbar.^2);
ph = gh .* mh / (2.* pi .* hbar.^2);

suml = 0;
Eli = 0;

for k = 1:limit(Ef, ml, V, C)
    Eli = E_i(ml, k , V,C);
    suml = suml +  Ef - Eli.* Eli/(e .* F_av(V,C) .* z_o(ml,V,C)).* ...
        (airy(-Eli/(e .* F_av(V,C) .* z_o(ml,V,C)))^2 + ...
        airy(1,-Eli/(e .* F_av(V,C) .* z_o(ml,V,C))).^2);
   
end

Ehi = 0;
sumh = 0;

for k = 1:limit(Ef, mh, V, C)
    Ehi = E_i(mh, k , V,C);
    sumh = sumh + Ef - Ehi.* Ehi/(e .* F_av(V,C) .* z_o(mh,V,C)).* ...
        (airy(-Ehi/(e .* F_av(V,C) .* z_o(mh,V,C)))^2 + ...
        airy(1,-Ehi/(e .* F_av(V,C) .* z_o(mh,V,C))).^2);
    
end

n3d = pl * suml + ph * sumh;
end

function [n] = n_2d(V,C)
%sheet carrier density in cm-2
[K_IL, d_nc, W, L, e, epsilon_0, hbar, m_e, K_STO, ml, gl, mh, gh] = unpack(C);
n = 1e-4 .* epsilon_0 .* K_IL .* V / (abs(e) .* d_nc);
end

function [z0] = z_o(m,V,C)
[K_IL, d_nc, W, L, e, epsilon_0, hbar, m_e, K_STO, ml, gl, mh, gh] = unpack(C);
z0 = (hbar.^2 / (2 .* m  .* e .* F_av(V,C))).^(1/3);
end

function [F] = F_av(V,C)
[K_IL, d_nc, W, L, e, epsilon_0, hbar, m_e, K_STO, ml, gl, mh, gh] = unpack(C);
F = e .* n_2d(V,C) /(2 .* epsilon_0 .* K_STO);
end

function [psi] = psi_i(z,m,i, V,C)
[K_IL, d_nc, W, L, e, epsilon_0, hbar, m_e, K_STO, ml, gl, mh, gh] = unpack(C);
psi = airy(z/z_o(m,V) - E_i(i,m, V,C)/(e .* F_av(V,C) .* z_o(m,V,C)));
end

function [E] = E_i(m,i, V, C)
[K_IL, d_nc, W, L, e, epsilon_0, hbar, m_e, K_STO, ml, gl, mh, gh] = unpack(C);
E = e .* F_av(V,C) .* z_o(m,V,C) .* (3 .* pi .* (i - 0.25) / 2).^(2/3);
end

function[l] = limit(Ef, m, V, C)
[K_IL, d_nc, W, L, e, epsilon_0, hbar, m_e, K_STO, ml, gl, mh, gh] = unpack(C);
l = floor((2/(3*pi)).*(Ef/(e * F_av(V,C) * z_o(m,V,C))).^(3/2) +1/4);
end
