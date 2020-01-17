nm = 1;
%ionic liquid dielectric constant
K_IL = 15;
%nanocapacitor gap size
d_nc = 5e-9 .* nm;
%width of channel
W = 10e-6 .* nm;
%length of the channel
L= 100e-6 .* nm;


%elementary charge
e =  1.602e-19;
eV = abs(e);
epsilon_0 = 8.85e-12 ;
hbar = 1.05457182e-34 ;
m_e = 9.11e-31;
%dielectric constant of STO
K_STO = 303;
ml = 1.2 .* m_e;
gl = 4;
mh = 4.8 .* m_e;
gh = 2;

format short

Vg = 2.5;
Vsd = 0.5;
n_list = [1.1; 3.0; 6.5; 7.8] .* 1e17;
Ef_guess_list = [1e-3 5e-3 20e-3 25e-3] .* eV;

V_list = [2.5; 2.7; 3; 3.2];
delta_Ef = 0.000005 .* eV;
tolerance = 0.005;

d_av = zeros(4,1);
n_av = zeros(4,1);
Ef_list = zeros(4,1);

figure(1)
colormap jet



for index = 1:4
    clear int_n3d int_n3d_p
    z = linspace(0,40e-9,10000);
    
    %sheet carrier density in m-2 position dependent
    % n_2d =  epsilon_0 .* K_IL .* Vg ./ (abs(e) .* d_nc);
    n_2d = n_list(index);
    
    Ef_guess = Ef_guess_list(index);
    
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
    
    Ef_list(index) = Ef;
    d_av(index) = simps(z, (z.*n_3d(z, Ef, Constants,Variables)))/simps(z, (n_3d(z, Ef, Constants,Variables)));
    n_av(index) = simps(z, (n_3d(z, Ef, Constants,Variables)).^2)/simps(z, n_3d(z, Ef, Constants,Variables));
    
    figure (4)
    E_test = linspace(0, Ef.*1.05,1000);
    dos = zeros(length(E_test),1);
    for idx = 1:length(E_test)
        dos(idx) =Dos(E_test(idx), Constants,Variables);
    end
    plot(dos ./(1e18) .* (eV),E_test ./ (1e-3 * eV),dos ./(1e18) .* (eV), ones(length(E_test),1) .* Ef ./ (1e-3 * eV))
    ylabel("\epsilon (meV)")
    xlabel("DOS (nm^{-2} eV^{-1})")
%     set(gca,'xtick',[])
    n = n_3d(z, Ef, Constants,Variables);
    figure(1)
    hold on
    plot((z .* 1e9),(n ./ 1e6))
    set(gca,'YScale', 'log')
    ylim([1e18 max(n ./ 1e6)*10])
    xlim([0 40])
    
end

xlabel("d (nm)")
ylabel("n_{3d} (cm^{-3})")
leg1 = legend(num2str(V_list));
legend('boxoff')
title(leg1, "V_g (V)")

figure(2)
hold on
plot(V_list,(d_av .* 1e9), '-r*')
plot(V_list,(n_av ./ 1e25), '-bo')
hold off
xlabel("V_g (V)")
leg = legend(["<d> (nm)";"<n_{3d}> (\times 10^{19} cm^{-2})"]);
legend('boxoff')

figure(3)
plot((V_list),(Ef_list ./ e), '-ro')
xlabel("V_g (V)")
ylabel("E_f (eV)")



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
    int_suml = simps(z, sum_templ);
    delta_E_l = Ef -  E_i(z_o_l, i, C, V);
    suml = suml + light .* delta_E_l .* (sum_templ ./ int_suml);
end
sumh = 0;
for i = 1:(limit_h)
    sum_temph = abs(psi_i(z,z_o_h,i, C,V)).^2;
    int_sumh = simps(z, sum_temph);
    delta_E_h = Ef -  E_i(z_o_h, i, C, V);
    sumh = sumh + heavy .* delta_E_h .* (sum_temph ./ int_sumh);
end

n3d = suml + sumh;
end
function [dos] = Dos(E, C,V)

[K_IL, d_nc, W, L, e, epsilon_0, hbar, m_e, K_STO, ml, gl, mh, gh] = unpack_C(C);
[n_2d, F_av, z_o_l, z_o_h] = unpack_V(V);

light = gl .* ml ./ (2 .* pi .* hbar.^2);
heavy = gh .* mh ./ (2 .* pi .* hbar.^2);

E_state = 0;
i=0;

while E_state <= E
    i = i + 1;
    E_state = E_i(z_o_l, i, C, V);
end
if i==0
    suml = 0;
else
    suml = light .* (i-1);
end
E_state = 0;
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
dos = suml + sumh;
end
function edos = E_Dos(dos,C,V)
end