%ionic liquid dielectric constant
K_IL = 15;
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
%Molar Volume
Vm = 0.4264/1405;
% reciprocal molar volume
c = (Vm)^(-1);
%Avogadro's number
NA = 6.02214086e23;
%temperature
T = 4;
%dielectric constant of STO
K_STO = 303;
epsIL = epsilon_0 .* K_IL;
%Helmholtz layer size
d_H = 5e-9;%(3 .* Vm ./(4 .* pi .* NA)) .^ (1/3);

%effective mass and valley degeneracy for light and heavy electrons
ml = 1.2 .* m_e;
gl = 4;
mh = 4.8 .* m_e;
gh = 2;

format short

Vsd = 0.5;
Vg_list = linspace(2,3.5,10);
V_list = [2.5; 2.7; 3.0; 3.2];
n_list = [1.1; 3.0; 6.5; 7.8] .* 1e17;

Ef_guess = 1e-3 .* eV;
delta_Ef = 0.000001 .* eV;
tolerance = 0.01;

qtest = linspace(1e17 .* e,8e17  .* e, 1000);
syms q


n_gcs = zeros(4,1);

D_V = @(q) d_H .* q ./ epsIL +  R .* T ./ (F) .* ...
    acosh(exp( q.^2 ./ (4 .* R .* T .* epsIL .* c)));

for index = 1:4

Vg = V_list(index);
    
n_2d = solve(D_V(q) == Vg);
n_2d = double(n_2d) ./ e;

n_gcs(index) = n_2d;
end

figure(1)
hold on
plot(V_list, n_list ./ 1e17,'-r.', 'MarkerSize', 20)
plot(V_list, n_gcs ./ 1e17, '-b.', 'MarkerSize', 20)
xlabel("V_g (V)")
ylabel("n_{2d} (cm^{-2} \times 10^{13})")
legend(["Ueno et al.";"GCS Model"]);
legend('boxoff')
