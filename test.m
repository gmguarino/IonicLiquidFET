
%ionic liquid dielectric constant
K_IL = 15;
%nanocapacitor gap size
d_nc = 5e-9;
%width of channel
W = 10e-6;
%length of the channel
L= 100e-6;
%surface area
A =  W .* L;

%elementary charge
e =  1.602e-19;
eV = abs(e);
epsilon_0 = 8.85e-12 ;
hbar = 1.05457182e-34 ;
m_e = 9.11e-31;
%gas constant
R = 8.3144598;
%Faradaic constant
F = 96485.33289;
% reciprocal molar volume
c = (3.0725e-4)^(-1);
%temperature
T = 298;
%dielectric constant of STO
K_STO = 303;
epsIL = epsilon_0 .* K_IL;
%Debye length
debye = sqrt(R .* T .* epsIL ./ (2 .* F.^2 .* c));

%effective mass and valley degeneracy for light and heavy electrons
ml = 1.2 .* m_e;
gl = 4;
mh = 4.8 .* m_e;
gh = 2;

format short

Vsd = 0.5;
n_list = [1.1; 3; 6.5; 7.8] .* 1e17;
Ef_guess = 5e-3 .* eV;
delta_Ef = 0.00005 .* eV;
tolerance = 0.01;

% q = linspace(1e17 .* e,8e17  .* e, 1000);
syms q

D_V = @(q) d_nc .* q ./ epsIL +  R .* T ./ (F .* abs(q)) .* ...
    acosh(exp( q.^2 ./ (4 .* R .* T .* epsIL .* c)));

n_2d = solve(D_V(q) == 2.5);
n_2d = double(n_2d) ./ e;

z = linspace(0,50e-9,10000);

A = 4.097e-5;
B = 4.097e-10;
%average electric field induced in the STO layer
%F_av = e .* n_2d ./(2 .* epsilon_0 .* K_STO);
F_av = (A / B) .* (exp( B .* e .* n_2d ./ (2 .* epsilon_0)) - 1);


z_o_l = (hbar.^2 ./ (2 .* ml  .* e .* F_av)).^(1/3);
z_o_h = (hbar.^2 ./ (2 .* mh  .* e .* F_av)).^(1/3);

Constants = [K_IL d_nc W L e epsilon_0 hbar m_e K_STO ml gl mh gh];
Variables = [n_2d F_av z_o_l z_o_h] ;
stop = false;
Ef = Ef_guess;
while stop == false
    if E_i(z_o_h, 1, Constants, Variables) <= Ef
        
        int_n3d = trapz(z, n_3d(z, Ef, Constants,Variables));
        if abs(int_n3d - n_2d)  >= tolerance * n_2d
            int_n3d_p = trapz(z, n_3d(z, Ef + delta_Ef, Constants,Variables))
            if abs(int_n3d_p - n_2d)  < abs(int_n3d - n_2d)
                Ef = Ef + delta_Ef
            else
                Ef = Ef - delta_Ef
            end
        else
            stop = true;
        end
    elseif E_i(z_o_h, 1, Constants, Variables) > Ef
        Ef = Ef + delta_Ef
    end
end
n = n_3d(z, Ef, Constants,Variables);
figure(1)
hold on
plot((z .* 1e9),(n ./ 1e6))
set(gca,'YScale', 'log')
ylim([1e18 max(n ./ 1e6)*10])
xlim([0 50])




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
psi = airy(z./z_o - E_i(z_o, i, C,V)./(e .* F_av .* z_o));

end
function [n3d] = n_3d(z, Ef, C,V)

[K_IL, d_nc, W, L, e, epsilon_0, hbar, m_e, K_STO, ml, gl, mh, gh] = unpack_C(C);
[n_2d, F_av, z_o_l, z_o_h] = unpack_V(V);

limit_l = floor((2/(3*pi)).*(Ef./(e .* F_av .* z_o_l)).^(3/2) + 1/4 );
limit_h = floor((2/(3*pi)).*(Ef./(e .* F_av .* z_o_h)).^(3/2) + 1/4 );

light = gl .* ml ./ (2 .* pi .* hbar.^2);
heavy = gh .* mh ./ (2 .* pi .* hbar.^2);


suml = 0;
for i = 1:(limit_l)
    sum_templ = abs(psi_i(z,z_o_l,i, C,V)).^2;
    int_suml = trapz(z, sum_templ);
    delta_E_l = Ef -  E_i(z_o_l, i, C, V);
    suml = suml + light .* delta_E_l .* (sum_templ ./ int_suml);
end
sumh = 0;
for i = 1:(limit_h)
    sum_temph = abs(psi_i(z,z_o_h,i, C,V)).^2;
    int_sumh = trapz(z, sum_temph);
    delta_E_h = Ef -  E_i(z_o_h, i, C, V);
    sumh = sumh + heavy .* delta_E_h .* (sum_temph ./ int_sumh);
end

n3d = suml + sumh;
end
