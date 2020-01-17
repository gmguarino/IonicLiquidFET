
%ionic liquid dielectric constant
K_IL = 15;
%width of channel
W = 10e-6;
%length of the channel
L= 100e-6;
%surface area
Area =  W .* L;

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
%Molar Volume
Vm = 0.4264/1405;
% reciprocal molar volume
c = (Vm)^(-1);
%Avogadro's number
NA = 6.02214086e23;
%temperature
T_list = [4.2;10;15;20;30;40;50;60;65;70;77;100;120;200;280;300];
F_list = zeros(length(T_list),1);
A_list = [4.097;4.782;5.446;6.175;8.430;12.64;19.58;31.84;39.37;44.84;53.19;78.13;109.9;192.3;280.1;303] .* 1e-5;
B_list = [4.907;4.887;4.848;4.817;4.438;3.777;3.156;0.9852;0;0;0;0;0;0;0;0] .* 1e-10;
%dielectric constant of STO
K_STO = 303;
epsIL = epsilon_0 .* K_IL;
%Helmholtz layer size
d_nc = 5e-9;%(3 .* Vm ./(4 .* pi .* NA)) .^ (1/3);

%effective mass and valley degeneracy for light and heavy electrons
ml = 1.2 .* m_e;
gl = 4;
mh = 4.8 .* m_e;
gh = 2;

format short

Vsd = 0.5;
Vg_list = linspace(2,3.5,10);
Ef_guess = 1e-3 .* eV;
delta_Ef = 0.00005 .* eV;
tolerance = 0.01;

qtest = linspace(-5e17 .* e,5e17  .* e, 100000);
syms q

D_V = @(q) d_nc .* q ./ epsIL +  R .* T ./ (F .* abs(q)) .* ...
    acosh(exp( q.^2 ./ (4 .* R .* T .* epsIL .* c)));
figure(2)
plot(qtest .* 1e3 , D_V(qtest),'-b', 'LineWidth',2)
xlabel("q (mC m^{-2})")
ylabel("\Delta\phi (V)")

% xlim([(-1e21 .* e) (1e21  .* e)])

d_av = zeros(10,1);
n_av = zeros(10,1);
Ef_list = zeros(10,1);



for index = 1:length(T_list)

Vg = 2.5;

T = T_list(index);
A = A_list(index);
B = B_list(index);
    


n_2d = vpasolve(D_V(q) == Vg, q, [0 Inf]);
n_2d = double(n_2d) ./ e;

z = linspace(0,100e-9,10000);


%average electric field induced in the STO layer
if T<65
    F_av = (A ./ B) .* (exp( B .* e .* n_2d ./ (2 .* epsilon_0)) - 1);
else
    F_av = e .* n_2d ./(2 .* epsilon_0 .* (1/A));
end
F_list(index) = F_av;

z_o_l = (hbar.^2 ./ (2 .* ml  .* e .* F_av)).^(1/3);
z_o_h = (hbar.^2 ./ (2 .* mh  .* e .* F_av)).^(1/3);

Constants = [K_IL d_nc W L e epsilon_0 hbar m_e K_STO ml gl mh gh];
Variables = [n_2d F_av z_o_l z_o_h] ;
stop = false;
Ef = Ef_guess;

while stop == false
    if E_i(z_o_h, 1, Constants, Variables) < Ef
        E_i(z_o_h, 1, Constants, Variables)
        Ef
        max(n_3d(z, Ef, Constants,Variables))
        int_n3d = trapz(z, n_3d(z, Ef, Constants,Variables));
        if abs(int_n3d - n_2d)  >= tolerance * n_2d
            int_n3d_p = trapz(z, n_3d(z, Ef + delta_Ef, Constants,Variables))
            if abs(int_n3d_p - n_2d)  < abs(int_n3d - n_2d)
                Ef = Ef + delta_Ef;
            else
                Ef = Ef - delta_Ef;
            end
        else
            stop = true;
        end
    elseif E_i(z_o_h, 1, Constants, Variables) >= Ef
        Ef = Ef + delta_Ef;
    end
end

Ef_list(index) = Ef;
d_av(index) = simps(z, (z.*n_3d(z, Ef, Constants,Variables)))/simps(z, (n_3d(z, Ef, Constants,Variables)));
n_av(index) = simps(z, (n_3d(z, Ef, Constants,Variables)).^2)/simps(z, n_3d(z, Ef, Constants,Variables));
end

figure(6)
hold on
plot(T_list,(Ef_list .* 1e3 ./ e), '-r*')
xlabel("T (K)")
ylabel("E_f (meV)")

figure(7)
hold on
plot(T_list,F_list, '-r.','MarkerSize', 20)
xlabel("T (K)")
ylabel("E_{AV} (V/m)")

figure(5)
hold on
plot(T_list,(d_av .* 1e9), '-r*')
plot(T_list,(n_av ./ 1e25), '-bo')
hold off
xlabel("T (K)")
leg = legend(["<d> (nm)";"<n_{3d}> (\times 10^{19} cm^{-3})"]);
legend('boxoff')

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
