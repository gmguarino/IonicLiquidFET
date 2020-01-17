%ionic liquid dielectric constant
K_IL = 15;
%nanocapacitor gap size
d_nc = 10e-9;
%width of channel
W = 10e-6;
%length of the channel
L= 100e-6;


%elementary charge
e =  -1.602e-19;
eV = abs(e);
epsilon_0 = 8.85e-12;
hbar = 1.05457182e-34;
m_e = 9.11e-31;
%dielectric constant of STO
K_STO = 303;

Vg = 3;
Vsd = 0.5;


ml = 1.2 .* m_e;
gl = 4;
mh = 4.8 .* m_e;
gh = 2;

z = linspace(0,1e-6,1000);


%sheet carrier density in m-2 position dependent 
n_2d =  epsilon_0 .* K_IL .* Vg ./ (abs(e) .* d_nc);
%n_2d = 7.8e17 ./ (nm).^2;

Ef = 0;
n = 0;
while (n<n_2d)
    

z = linspace(0,1e-6,1000);
Ef = Ef + 0.01 .* eV;

A = 4.097e-5;
B = 4.907e-10;
%average electric field induced in the STO layer
F_av = e .* n_2d ./(2 .* epsilon_0 .* K_STO);
%F_av = (A / B) .* (exp( B .* e .* n_2d ./ (2 .* epsilon_0)) - 1);


z_o_l = (hbar.^2 ./ (2 .* ml  .* e .* F_av)).^(1/3);
z_o_h = (hbar.^2 ./ (2 .* mh  .* e .* F_av)).^(1/3);

E_i_l =@(i) e .* F_av .* z_o_l .* (3 .* pi .* (i - 0.25) ./ 2).^(2/3);
E_i_h =@(i) e .* F_av .* z_o_h .* (3 .* pi .* (i - 0.25) ./ 2).^(2/3);


psi_i_l =@(i) airy(0,z./z_o_l - (3 .* pi .* (i - 0.25) ./ 2).^(2/3));
psi_i_h =@(i) airy(0,z./z_o_h - (3 .* pi .* (i - 0.25) ./ 2).^(2/3));

% syms z V E x

Constants = [K_IL d_nc W L e epsilon_0 hbar m_e K_STO ml gl mh gh];
Variables = [n_2d F_av z_o_l z_o_h] ;

figure(2)
hold on
% sumpsi = 0;
% for i = 1:10
%     sumpsi = sumpsi + abs(psi_i_l(i)).^2;
% end
plot(z, n_3d(z, Ef, Constants,Variables));
n = integral(@(z)n_3d(z, Ef, Constants,Variables),0,1e-7,'ArrayValued', true)
end

% set(gca,'YScale', 'log')

% trapz(z, n_3d(z, 10*eV, Constants,Variables))
% 
% n = 0;
% E = 0;
% while n < double(n_2d)
%     E = E + 10 * eV;
%     n = trapz(z, n_3d(z, E, Constants,Variables));
%     n = double(n);
% end
% E
% plot(z, n_3d)



function [K_IL, d_nc, W, L, e, epsilon_0, hbar, m_e, K_STO, ml, gl, mh, gh] = unpack_C(c)
constant = num2cell(c); 
[K_IL, d_nc, W, L, e, epsilon_0, hbar, m_e, K_STO, ml, gl, mh, gh] = constant{:};

end
function [n_2d, F_av, z_o_l, z_o_h] = unpack_V(v)
variables = num2cell(v); 
[n_2d, F_av, z_o_l, z_o_h] = variables{:};

end
function [E] = E_i(z_o, i, C, V)
[K_IL, d_nc, W, L, e, epsilon_0, hbar, m_e, K_STO, ml, gl, mh, gh] = unpack_C(C);
[n_2d, F_av, z_o_l, z_o_h] = unpack_V(V);

E = e .* F_av .* z_o .* (3 .* pi .* (i - 0.25) ./ 2).^(2/3);
end

function [psi] = psi_i(z, z_o, i, C, V)

[K_IL, d_nc, W, L, e, epsilon_0, hbar, m_e, K_STO, ml, gl, mh, gh] = unpack_C(C);
[n_2d, F_av, z_o_l, z_o_h] = unpack_V(V);
psi = airy(0,z./z_o - (3 .* pi .* (i - 0.25) ./ 2).^(2/3));

end
function [n3d] = n_3d(z, Ef, C,V)

[K_IL, d_nc, W, L, e, epsilon_0, hbar, m_e, K_STO, ml, gl, mh, gh] = unpack_C(C);
[n_2d, F_av, z_o_l, z_o_h] = unpack_V(V);

limit_l = floor((2/(3*pi)).*(Ef./(e .* F_av .* z_o_l)).^(3/2) + 1/4 ); 
limit_h = floor((2/(3*pi)).*(Ef./(e .* F_av .* z_o_h)).^(3/2) + 1/4 ); 
suml = 0;
for i = 1:limit_l
    suml = suml + (Ef - E_i(z_o_l, i, C,V)).*abs(psi_i(z,z_o_l,i, C,V)).^2;
end
sumh = 0;
for i = 1:limit_h
     sumh = sumh + (Ef - E_i(z_o_h, i, C,V)).*abs(psi_i(z,z_o_h,i, C,V)).^2;
end

n3d =(gl.*ml/(2 .* pi .* hbar.^2)) .* double(suml) + (gh .* mh ./ (2.* pi .* hbar.^2)) .* double(sumh);
%n3d = double(suml);

end
function [psip] = psipl(z, Ef, C,V)

[K_IL, d_nc, W, L, e, epsilon_0, hbar, m_e, K_STO, ml, gl, mh, gh] = unpack_C(C);
[n_2d, F_av, z_o_l, z_o_h] = unpack_V(V);

limit_l = floor((2/(2.*pi.*e.*F_av)).*(2.*ml/hbar).^(1/2).*Ef.^(3/2) + 1/4); 

suml = 0;
for i = 1:limit_l
    suml = suml + psi_i(z,z_o_l,i, C,V);
end

psip = double(suml);

end