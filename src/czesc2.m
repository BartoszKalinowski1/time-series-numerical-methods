function czesc2
clear; close all;
R1 = 0.1;
R2 = 10;
C = 0.5;
L1 = 3;
L2 = 5;

u_data = [20  50  100  150  200  250  280  300];
M_data = [0.46 0.64 0.78 0.68 0.44 0.23 0.18 0.18];
 
% Dane równoodległe (co 30 V) używane przez splajn, zrobiłem ich 10 nie 8, bo sa
% fajne wartosci całkowite (ale moze przez to dawac troche dokładniejsze wyniki niz 8 wezlow)
% Rysunek skad wziałem pkt. równoodległe załączam w raporcie
u_rown = 30:30:300;                          
M_rown = [0.54 0.70 0.78 0.76 0.68 ...
          0.51 0.37 0.26 0.19 0.18];
 
coeff_spline = splajn3_przygotuj(u_rown, M_rown);
coeff_mnk3   = oblicz_wspolczynniki_mnk3(u_data, M_data);
coeff_mnk5   = oblicz_wspolczynniki_mnk5(u_data, M_data);
 
metodyM = { ...
    @(u1) interpolacja_lagrange(u1, u_data, M_data), ...
    @(u1) interpolacja_splajnami(u1, coeff_spline), ...
    @(u1) aproksymacja_mnk3(u1, coeff_mnk3), ...
    @(u1) aproksymacja_mnk5(u1, coeff_mnk5) };
 
nazwyM  = {'Lagrange', 'Splajn 3°', 'MNK 3°', 'MNK 5°'};
kolory  = {'b', 'r', 'g', 'm'};
styl    = {'-', '--', '-.', ':'};
 
e2 = @(x) 240 * sin(x);
e4 = @(x) 120 * sin(2 * pi * 50 * x);
 
wym   = {e2, e4};
nazwyW = {'240\cdotsin(t)', '120\cdotsin(2\pi\cdot50\cdott)'};

h_vec = [0.5 0.009];
tk_vec = [30, 0.1];
 
for w = 1:2
    e_fun = wym{w};

    h2 = h_vec(w);
    t0 = 0;
    tk = tk_vec(w);
    t2 = t0:h2:tk;
    N2 = length(t2);
 
    i1_all  = zeros(4, N2);
    i2_all  = zeros(4, N2);
    uc_all  = zeros(4, N2);
    uR2_all = zeros(4, N2);
    e_plot  = arrayfun(e_fun, t2);   
 
    for m = 1:4
        metoda_M = metodyM{m};
 
        i1 = 0;  i2 = 0;  uc = 0;
 
        i1_all(m,1)  = i1;
        i2_all(m,1)  = i2;
        uc_all(m,1)  = uc;
        uR2_all(m,1) = R2 * i2;
 
        for j = 1:N2-1
            ej = e_fun(t2(j));
 
            uL1 = abs(ej - R1*i1 - uc);
            uL1 = max(20, min(300, uL1));
 
            M_n = metoda_M(uL1);
            M_n = max(0.10, min(1.00, M_n));
 
            [di1, di2, duC] = pochodne(i1, i2, uc, ej, R1, R2, C, L1, L2, M_n);

            i1h  = i1 + (h2/2) * di1;
            i2h  = i2 + (h2/2) * di2;
            uch  = uc + (h2/2) * duC;
            ejh  = e_fun(t2(j) + h2/2);
 
            uL1h = abs(ejh - R1*i1h - uch);
            uL1h = max(20, min(300, uL1h));
            M_nh = metoda_M(uL1h);
            M_nh = max(0.10, min(1.00, M_nh));
 
            [di1h, di2h, duCh] = pochodne(i1h, i2h, uch, ejh, R1, R2, C, L1, L2, M_nh);

            i1 = i1 + h2 * di1h;
            i2 = i2 + h2 * di2h;
            uc = uc + h2 * duCh;
 
            i1_all(m, j+1)  = i1;
            i2_all(m, j+1)  = i2;
            uc_all(m, j+1)  = uc;
            uR2_all(m, j+1) = R2 * i2;
        end
    end

    % i1(t) 
    figure('Name', ['i1(t) – ' nazwyW{w}], 'NumberTitle', 'off', ...
           'Position', [50, 50, 1000, 550]);
    hold on;
    for m = 1:4
        plot(t2, i1_all(m,:), styl{m}, 'Color', kolory{m}, 'LineWidth', 1.3, ...
             'DisplayName', nazwyM{m});
    end
    hold off;
    title(['Prąd i_1(t) – e(t) = ' nazwyW{w}], 'Interpreter', 'tex');
    xlabel('t [s]');  ylabel('i_1 [A]');
    legend('Location', 'best');
    grid on;  xlim([0 tk]);
 
    % i2(t) 
    figure('Name', ['i2(t) – ' nazwyW{w}], 'NumberTitle', 'off', ...
           'Position', [100, 80, 1000, 550]);
    hold on;
    for m = 1:4
        plot(t2, i2_all(m,:), styl{m}, 'Color', kolory{m}, 'LineWidth', 1.3, ...
             'DisplayName', nazwyM{m});
    end
    hold off;
    title(['Prąd i_2(t) – e(t) = ' nazwyW{w}], 'Interpreter', 'tex');
    xlabel('t [s]');  ylabel('i_2 [A]');
    legend('Location', 'best');
    grid on;  xlim([0 tk]);
 
    % uR2(t) 
    figure('Name', ['uR2(t) – ' nazwyW{w}], 'NumberTitle', 'off', ...
           'Position', [150, 110, 1000, 550]);
    hold on;
    for m = 1:4
        plot(t2, uR2_all(m,:), styl{m}, 'Color', kolory{m}, 'LineWidth', 1.3, ...
             'DisplayName', nazwyM{m});
    end
    hold off;
    title(['Napięcie u_{R2}(t) – e(t) = ' nazwyW{w}], 'Interpreter', 'tex');
    xlabel('t [s]');  ylabel('u_{R2} [V]');
    legend('Location', 'best');
    grid on;  xlim([0 tk]);
 
    % uC(t) 
    figure('Name', ['uC(t) – ' nazwyW{w}], 'NumberTitle', 'off', ...
           'Position', [200, 140, 1000, 550]);
    hold on;
    for m = 1:4
        plot(t2, uc_all(m,:), styl{m}, 'Color', kolory{m}, 'LineWidth', 1.3, ...
             'DisplayName', nazwyM{m});
    end
    hold off;
    title(['Napięcie u_C(t) – e(t) = ' nazwyW{w}], 'Interpreter', 'tex');
    xlabel('t [s]');  ylabel('u_C [V]');
    legend('Location', 'best');
    grid on;  xlim([0 tk]);
 
   
end

end

function [di1, di2, duC] = pochodne(i1, i2, uc, e, R1, R2, C, L1, L2, M)
    D = L1*L2 - M^2;

    di1 = ( L2*(e - R1*i1 - uc) + M*R2*i2 ) / D;
    di2 = ( -M*(e - R1*i1 - uc) - L1*R2*i2 ) / D;
    duC = i1/ C;
end
 
% INTERPOLACJA WIELOMIANOWA (JA UŻYŁEM INTERPOLACJI WIELOMIANEM LAGRANGE'A) 
function M_val = interpolacja_lagrange(u1, u_nodes, M_nodes)
    n = length(u_nodes);
    M_val = 0;
    for i = 1:n
        L = 1;
        for j = 1:n
            if i ~= j
                L = L * (u1 - u_nodes(j)) / (u_nodes(i) - u_nodes(j));
            end
        end
        M_val = M_val + L * M_nodes(i);
    end
end
 
% INTERPOLACJA FUNKCJAMI SKLEJANYMI 3. STOPNIA 
function coeff = splajn3_przygotuj(x_nodes, y_nodes)
    n = length(x_nodes) - 1;
    h = x_nodes(2) - x_nodes(1);   % węzły równoodległe
 
    % Warunki brzegowe: zerowe pochodne pierwszego rzędu, założyłem to, bo
    % w projekcie nie było podanych tych wartości
    alpha = 0;
    beta  = 0;
 
    A = zeros(n+1);
    b = zeros(n+1, 1);
 
    A(1,1) = 4;  A(1,2) = 2;
    b(1)   = y_nodes(1) + h/3 * alpha;
 
    for i = 2:n
        A(i, i-1) = 1;
        A(i, i)   = 4;
        A(i, i+1) = 1;
        b(i)      = y_nodes(i);
    end
 
    A(n+1, n+1) = 4;  A(n+1, n) = 2;
    b(n+1)      = y_nodes(n+1) - h/3 * beta;
 
    c = zeros(n+3, 1);
    c(2:n+2) = A \ b;
    c(1)     = c(3)   - h/3 * alpha;
    c(n+3)   = c(n+1) + h/3 * beta;
 
    coeff.c  = c;
    coeff.h  = h;
    coeff.x0 = x_nodes(1);
    coeff.n  = n;
end
 
function y = phi_bspline(xi, h, x)
    % funkcja bazowa
    if x < xi - 2*h || x > xi + 2*h
        y = 0;
    elseif x < xi - h
        y = (x - (xi - 2*h))^3;
    elseif x < xi
        y = (x - (xi - 2*h))^3 - 4*(x - (xi - h))^3;
    elseif x < xi + h
        y = ((xi + 2*h) - x)^3 - 4*((xi + h) - x)^3;
    else
        y = ((xi + 2*h) - x)^3;
    end
    y = y / h^3;
end
 
function M_val = interpolacja_splajnami(u1, coeff)
    c  = coeff.c;
    h  = coeff.h;
    x0 = coeff.x0;
    n  = coeff.n;
    M_val = 0;
    for i = -1:n+1
        xi    = x0 + i * h;
        M_val = M_val + c(i+2) * phi_bspline(xi, h, u1);
    end
end
 
% APROKSYMACJA WIELOMIANEM 3. STOPNIA 
function coeff = oblicz_wspolczynniki_mnk3(u_nodes, M_nodes)
    V     = [ones(numel(u_nodes),1), u_nodes(:), u_nodes(:).^2, u_nodes(:).^3];
    coeff = (V' * V) \ (V' * M_nodes(:));
end
 
function y = aproksymacja_mnk3(u1, coeff)
    y = coeff(1) + coeff(2)*u1 + coeff(3)*u1^2 + coeff(4)*u1^3;
end
 
% APROKSYMACJA WIELOMIANEM 5. STOPNIA 
function coeff = oblicz_wspolczynniki_mnk5(u_nodes, M_nodes)
    u_max  = max(u_nodes);
    us     = u_nodes(:) / u_max;
    V      = [ones(numel(us),1), us, us.^2, us.^3, us.^4, us.^5];
    coeff  = (V' * V) \ (V' * M_nodes(:));
    coeff  = [coeff; u_max];
end
 
function y = aproksymacja_mnk5(u1, coeff)
    u_max  = coeff(end);
    c      = coeff(1:end-1);
    us     = u1 / u_max;
    y = c(1) + c(2)*us + c(3)*us^2 + c(4)*us^3 + c(5)*us^4 + c(6)*us^5;
end