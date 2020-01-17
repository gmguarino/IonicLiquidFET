
%ionic liquid dielectric constant
K_IL = 15;
%nanocapacitor gap size
d_nc = 5e-9;
%width of channel
W = 10e-6;
%length of the channel
L= 100e-6;


%elementary charge
e = - 1.602e-19;
eV = abs(e);
epsilon_0 = 8.854e-12;
hbar = 1.05457182e-34;
m_e = 9.11e-31;
%dielectric constant of STO
K_STO = 265;

Vg = 3;
Vsd = 0.5;
Ef = 0.1 * eV;
ml = 1.2 .* m_e;
gl = 4;
mh = 4.8 .* m_e;
gh = 2;

zz = linspace(0,1e-8,1000);
xx = linspace(0,L,1000);

[x,z] = meshgrid(xx,zz);

V_x = Vg + (Vsd./L) .* x;

%sheet carrier density in m-2 position dependent 
n_2d =  epsilon_0 .* K_IL .* V_x ./ (abs(e) .* d_nc);


%average electric field induced in the STO layer
F_av = e .* n_2d ./(2 .* epsilon_0 .* K_STO);

z_o_l = (hbar.^2 ./ (2 .* ml  .* e .* F_av)).^(1/3);
z_o_h = (hbar.^2 ./ (2 .* mh  .* e .* F_av)).^(1/3);

E_i_l =@(i) e .* F_av .* z_o_l .* (3 .* pi .* (i - 0.25) ./ 2).^(2/3);
E_i_h =@(i) e .* F_av .* z_o_h .* (3 .* pi .* (i - 0.25) ./ 2).^(2/3);

% syms z V E x

psi_i_l =@(i) airy(z./z_o_l - E_i_l(i)./(e .* F_av .* z_o_l));
psi_i_h =@(i) airy(z./z_o_h - E_i_h(i)./(e .* F_av .* z_o_h));



Constants = [K_IL d_nc W L e epsilon_0 hbar m_e K_STO ml gl mh gh];


%Defining the limit of the summation
limit_l = floor((2/(3*pi)).*(Ef./(e .* F_av .* z_o_l)).^(3/2) + 1/4 ); 
limit_h = floor((2/(3*pi)).*(Ef./(e .* F_av .* z_o_h)).^(3/2) + 1/4 ); 
suml = 0;
for i = 1:limit_l
    suml = suml + vpa((Ef - E_i_l(i)).*abs(psi_i_l(i)).^2);
end
sumh = 0;
for i = 1:limit_h
     sumh = sumh + vpa((Ef - E_i_h(i)).*abs(psi_i_h(i)).^2);
end

n_3d =(gl.*ml/(2 .* pi .* hbar.^2)) .* suml + (gh .* mh ./ (2.* pi .* hbar.^2)) .* sumh;



figure(1)

% x = linspace(-2*pi,2*pi);
% z = linspace(0,4*pi);
% [Z,X] = meshgrid(linspace(0,1e-7,100),linspace(0,L,100));
% for i=1:
% N = arrayfun(real((n_3d(Z,Ef,Vds,Vg,X))), Z,X, 'un', 0);
% 
[C,h] = contourf(x, z, n_3d,1000);
h.LineStyle = 'none';
% %Mod squared integral of the wavefunctions
% int_psi_sq = @(k,E, m, V) E_i(k,m, V)/(e .* F_av(V) .* z_o(m,V)).* airy(0,-E_i(k,m, V)/(e .* F_av(V) .* z_o(m,V))).^2 + airy(1,-E_i(k,m, V)/(e .* F_av(V) .* z_o(m,V))).^2;
% 
% n_int = @(Ef,V) (gl.*ml.* z_o(ml,V)/(2 .* pi .* hbar.^2)) .* ...
%     sum((Ef - E_i(ml,1:limit(Ef,ml,V),V)).*int_psi_sq(1:limit(Ef,ml,V),Ef, ml,V)) + ...
%     (gh .* mh .* z_o(mh,V)./ (2.* pi .* hbar.^2)) .* sum((Ef - E_i(ml,1:limit(Ef,ml,V),V)).*int_psi_sq(1:limit(Ef,ml,V),Ef, mh,V));
% 
% % 
% fplot(n_3d(1e-9, 0.01 * eV ,0.05, Vg, x),[0 L])

% fcountour(n_3d(z, 0.01 * eV ,0.05, Vg, x),[0 1e-7 0 L])
% %end
% 
% n=0;
% Ef = 0;
% while n <= n_2d(Vg)
%     Ef = Ef + 0.001*abs(e);
%     n =  vpa(n_int(0.1*abs(e), Vg));
% end
% n
% Ef
% % xlabel("depth (m)")
% ylabel("n_{3d}")
% legend("Vg = 1V", "Vg = 2V","Vg = 3V","Vg = 4V", "Location", "southwest")

% 

function [K_IL, d_nc, W, L, e, epsilon_0, hbar, m_e, K_STO, ml, gl, mh, gh] = unpack(c)
constant = num2cell(c); 
[K_IL, d_nc, W, L, e, epsilon_0, hbar, m_e, K_STO, ml, gl, mh, gh] = constant{:};

end

