
%ionic liquid dielectric constant
K_IL = 15;
%width of channel
W = 10e-6;
%length of the channel
L= 100e-6;

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

%temperature
T = 4;
%dielectric constant of STO
K_STO = 303;
epsIL = epsilon_0 .* K_IL;
%Helmholtz layer size
d_H = 5e-9;

%effective mass and valley degeneracy for light and heavy electrons
ml = 1.2 .* m_e;
gl = 4;
mh = 4.8 .* m_e;
gh = 2;

format short

%Defining the gate voltages
Vg_list = linspace(2,3.5,10);

%Initial guess for the Fermi energy is 1meV
Ef_guess = 1e-3 .* eV;
%Steps in Ef
delta_Ef = 0.000001 .* eV;
%Ef correct to 1%
tolerance = 0.01;

syms q
qtest = linspace(-5e17 .* e,5e17  .* e, 100000);

%initialising result arrays
d_av = zeros(10,1);
n_av = zeros(10,1);

%Defining the change in potential according to GCS theory
D_V = @(q) d_H .* q ./ epsIL +  R .* T ./ (F) .* ...
    acosh(exp( q.^2 ./ (4 .* R .* T .* epsIL .* c)));

%plot of GCS model behaviour
figure(2)
plot(qtest .* 1e3 , D_V(qtest),'-b', 'LineWidth',2)
xlabel("q (mC m^{-2})")
ylabel("\Delta\phi (V)")

%looping over voltages
for index = 1:length(Vg_list)
    
    %defining gate voltage
    Vg = Vg_list(index);
    
    %solve for n_2d as q/e
    n_2d = solve(D_V(q) == Vg);
    n_2d = double(n_2d) ./ e;
    
    %defining depth
    z = linspace(0,50e-9,10000);
    
    %A and B taken from literature
    A = 4.097e-5;
    B = 4.097e-10;
    %average electric field induced in the STO layer
    F_av = (A / B) .* (exp( B .* e .* n_2d ./ (2 .* epsilon_0)) - 1);
    
    %define z_o for both electron types
    z_o_l = (hbar.^2 ./ (2 .* ml  .* e .* F_av)).^(1/3);
    z_o_h = (hbar.^2 ./ (2 .* mh  .* e .* F_av)).^(1/3);
    
    %constants and variables to pass to functions
    Constants = [K_IL d_H W L e epsilon_0 hbar m_e K_STO ml gl mh gh];
    Variables = [n_2d F_av z_o_l z_o_h] ;
    %stop clause
    stop = false;
    %setting Ef to the initial guess
    Ef = Ef_guess;
    
    while stop == false
        if E_i(z_o_h, 1, Constants, Variables) < Ef
            %integrating n_3d
            int_n3d = trapz(z, n_3d(z, Ef, Constants,Variables));
            %checking if within tolerance level
            if abs(int_n3d - n_2d)  >= tolerance * n_2d
                %evaluating integral for Ef+dEf
                int_n3d_p = trapz(z, n_3d(z, Ef + delta_Ef, Constants,Variables));
                if abs(int_n3d_p - n_2d)  < abs(int_n3d - n_2d)
                    %updating Ef
                    Ef = Ef + delta_Ef;
                else
                    %updating Ef
                    Ef = Ef - delta_Ef;
                end
            else
                stop = true;
            end
            %checking if not below lowest eigenstate
        elseif E_i(z_o_h, 1, Constants, Variables) >= Ef
            %updating Ef
            Ef = Ef + delta_Ef;
        end
    end
    
    %calculate n_av and d_av
    d_av(index) = simps(z, (z.*n_3d(z, Ef, Constants,Variables)))/simps(z, (n_3d(z, Ef, Constants,Variables)));
    n_av(index) = simps(z, (n_3d(z, Ef, Constants,Variables)).^2)/simps(z, n_3d(z, Ef, Constants,Variables));
    
    figure (2)
    E_test = linspace(0, Ef.*1.05,1000);
    dos = zeros(length(E_test),1);
    %calculate DOS
    for idx = 1:length(E_test)
        dos(idx) =Dos(E_test(idx), Constants,Variables);
    end
    %plot DOS
    plot(dos ./(1e18) .* (eV),E_test ./ (1e-3 * eV),dos ./(1e18) .* (eV), ones(length(E_test),1) .* Ef ./ (1e-3 * eV))
    ylabel("\epsilon (meV)")
    xlabel("DOS (nm^{-2} eV^{-1})")
    
    %calculate depth profile and plot it for each Vg
    n = n_3d(z, Ef, Constants,Variables);
    figure(4)
    hold on
    plot((z .* 1e9),(n ./ 1e6))
    set(gca,'YScale', 'log')
    ylim([1e17 max(n ./ 1e6)*10])
    xlim([0 50])
end
xlabel("d (nm)")
ylabel("n_{3d} (cm^{-3})")
leg1 = legend(num2str(V_list));
legend('boxoff')
title(leg1, "V_g (V)")

%plot of n_av and d_av vs Vg
figure(3)
hold on
plot(Vg_list,(d_av .* 1e9), '-r*')
plot(Vg_list,(n_av ./ 1e25), '-bo')
hold off
xlabel("V_g (V)")
leg2 = legend(["<d> (nm)";"<n_{3d}> (\times 10^{19} cm^{-2})"]);
legend('boxoff')


function [K_IL, d_nc, W, L, e, epsilon_0, hbar, m_e, K_STO, ml, gl, mh, gh] = unpack_C(c)
%function to unpack constants
constant = num2cell(c);
[K_IL, d_nc, W, L, e, epsilon_0, hbar, m_e, K_STO, ml, gl, mh, gh] = constant{:};

end
function [n_2d, F_av, z_o_l, z_o_h] = unpack_V(v)
%function to unpack variables
variables = num2cell(v);
[n_2d, F_av, z_o_l, z_o_h] = variables{:};

end
function [E] = E_i(z_o, i, C, V)
%calculating eigenenergies
[K_IL, d_nc, W, L, e, epsilon_0, hbar, m_e, K_STO, ml, gl, mh, gh] = unpack_C(C);
[n_2d, F_av, z_o_l, z_o_h] = unpack_V(V);

E = e .* F_av .* z_o .* (3 .* pi .* (i - 0.25) ./ 2).^(2/3);
end

function [psi] = psi_i(z, z_o, i, C, V)
%calculating eigenstates
[K_IL, d_nc, W, L, e, epsilon_0, hbar, m_e, K_STO, ml, gl, mh, gh] = unpack_C(C);
[n_2d, F_av, z_o_l, z_o_h] = unpack_V(V);
psi = airy(z./z_o - E_i(z_o, i, C,V)./(e .* F_av .* z_o));

end
function [n3d] = n_3d(z, Ef, C,V)
%calculating the depth profile
[K_IL, d_nc, W, L, e, epsilon_0, hbar, m_e, K_STO, ml, gl, mh, gh] = unpack_C(C);
[n_2d, F_av, z_o_l, z_o_h] = unpack_V(V);

%defining summation limits
limit_l = floor((2/(3*pi)).*(Ef./(e .* F_av .* z_o_l)).^(3/2) + 1/4 );
limit_h = floor((2/(3*pi)).*(Ef./(e .* F_av .* z_o_h)).^(3/2) + 1/4 );

%defining prefactors
light = gl .* ml ./ (2 .* pi .* hbar.^2);
heavy = gh .* mh ./ (2 .* pi .* hbar.^2);

%summing over light electron contributions
suml = 0;
for i = 1:(limit_l)
    sum_templ = abs(psi_i(z,z_o_l,i, C,V)).^2;
    int_suml = trapz(z, sum_templ);
    delta_E_l = Ef -  E_i(z_o_l, i, C, V);
    suml = suml + light .* delta_E_l .* (sum_templ ./ int_suml);
end

%summing over heavy electron contributions
sumh = 0;
for i = 1:(limit_h)
    sum_temph = abs(psi_i(z,z_o_h,i, C,V)).^2;
    int_sumh = trapz(z, sum_temph);
    delta_E_h = Ef -  E_i(z_o_h, i, C, V);
    sumh = sumh + heavy .* delta_E_h .* (sum_temph ./ int_sumh);
end
%summing up all contributions
n3d = suml + sumh;
end
function [dos] = Dos(E, C,V)

[K_IL, d_nc, W, L, e, epsilon_0, hbar, m_e, K_STO, ml, gl, mh, gh] = unpack_C(C);
[n_2d, F_av, z_o_l, z_o_h] = unpack_V(V);

%defining prefactors
light = gl .* ml ./ (2 .* pi .* hbar.^2);
heavy = gh .* mh ./ (2 .* pi .* hbar.^2);

%setting the energy of the test state to 0
E_state = 0;
i=0;
%looping over light electron eigenstates to find the number of states that
%are allowed for an input energy
while E_state <= E
    i = i + 1;
    E_state = E_i(z_o_l, i, C, V);
end
if i==0
    suml = 0;
else
    suml = light .* (i-1);
end
%setting the energy of the test state to 0
E_state = 0;
%looping over heavy electron eigenstates to find the number of states that
%are allowed for an input energy
i=0;
while E_state <= E
    i = i + 1;
    E_state = E_i(z_o_h, i, C, V);
end
if i==0
    sumh =0;
else
    sumh = heavy .* (i-1);
end
%summing up all contributions
dos = suml + sumh;
end