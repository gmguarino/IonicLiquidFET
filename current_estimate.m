
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
eV = abs(e)
epsilon_0 = 8.854e-12;
hbar = 1.05457182e-34;
m_e = 9.11e-31;
%dielectric constant of STO
K_STO = 265;

%sheet carrier density in m-2
n_2d =@(V)  epsilon_0 .* K_IL .* V / (abs(e) .* d_nc);

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

%Defining the limit of the summation
limit = @(Ef, m, V) floor((2/(3*pi)).*(Ef/(e * F_av(V) * z_o(m,V))).^(3/2) + 1/4 ); 

n_3d =@(z, Ef,V) (gl.*ml/(2 .* pi .* hbar.^2)) .* sum((Ef - E_i(ml,1:limit(Ef,ml,V),V)).* abs(psi_i(z,ml,1:limit(Ef,ml,V),V)).^2) + ...
    (gh .* mh / (2.* pi .* hbar.^2)) .* sum((Ef - E_i(mh,1:limit(Ef,ml,V),V)).* abs(psi_i(z,mh,1:limit(Ef,ml,V),V)).^2);

% E_f =@(V) solve(integral(@(z) arrayfun(@(z) n_3d(z,E,V),z,'un',0),0,0.001) == n_2d, E);

% Ef =@(V) (gh .* hbar.^2/(2 .* mh )) .* (3 .* pi.^2 .*  n_2d(V)).^(2/3) + ...
%     (gl .* hbar.^2/(2 .* ml )) .* (3 .* pi.^2 .*  n_2d(V)).^(2/3);

%n_3d_av=@(E_f) integral(@(z) arrayfun(@(z) n_3d(z,E_f).*2,z, 'un', 0), 0, inf)/integral(@(z) arrayfun(@(z) n_3d(z,E_f),z, 'un', 0), 0, inf);

%d_2deg_av=@(E_f) integral(@(z) arrayfun(@(z) z * n_3d(z,E_f),z, 'un', 0), 0, inf)/integral(@(z) arrayfun(@(z) n_3d(z,E_f),z, 'un', 0), 0, inf);

%n_3d_av(0.001)

%d_2deg_av(0.001)



figure(1)



Vg = 3;
%Mod squared integral of the wavefunctions
int_psi_sq = @(k,E, m, V) E_i(k,m, V)/(e .* F_av(V) .* z_o(m,V)).* airy(0,-E_i(k,m, V)/(e .* F_av(V) .* z_o(m,V))).^2 + airy(1,-E_i(k,m, V)/(e .* F_av(V) .* z_o(m,V))).^2;

n_int = @(Ef,V) (gl.*ml.* z_o(ml,V)/(2 .* pi .* hbar.^2)) .* ...
    sum((Ef - E_i(ml,1:limit(Ef,ml,V),V)).*int_psi_sq(1:limit(Ef,ml,V),Ef, ml,V)) + ...
    (gh .* mh .* z_o(mh,V)./ (2.* pi .* hbar.^2)) .* sum((Ef - E_i(ml,1:limit(Ef,ml,V),V)).*int_psi_sq(1:limit(Ef,ml,V),Ef, mh,V));


% E_f = abs(vpasolve(feval(n_int, E,Vg) == feval(n_2d,Vg), E))


%for Vg = 1:4

% vpa(n_int(0.1*abs(e), Vg))

% n_2d(Vg)
% psi_i(1e-10,mh,2,Vg)
% % n_3d(1e-10,10e-3 * abs(e), Vg)
%fplot(n_int(E, Vg), [0, 0.01*abs(e)])
fplot(limit(E,ml,Vg),[0, 1e-19])
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

