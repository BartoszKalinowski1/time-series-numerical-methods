function czesc3
clear; close all;
R1 = 0.1;
R2 = 10;
C = 0.5;
L1 = 3;
L2 = 5;
M = 0.8;
t0 = 0;
tk = 30;
dt1 = 0.001;   % krok krótki 
dt2 = 0.4;     % krok długi 
kroki_c3 = [dt1, dt2];

u_rown = 30:30:300; 
M_rown = [0.54 0.70 0.78 0.76 0.68 ...
          0.51 0.37 0.26 0.19 0.18];

e_c3 = { @(x) 1, ...
          @(x) (mod(x, 3) < 1.5) * 120, ...
          @(x) 240 * sin(x), ...
          @(x) 210 * sin(2*pi*5*x), ...
          @(x) 120 * sin(2*pi*50*x) };
 
nazwy_c3 = {'1 V (stałe)', ...
            'Prostokąt T=3s, 120V', ...
            '240 sin(t)', ...
            '210 sin(10\pi t)', ...
            '120 sin(100\pi t)'};
 
E_prost_lin = zeros(5, 2);   
E_parab_lin = zeros(5, 2);
E_prost_nlin = zeros(5, 2);  
E_parab_nlin = zeros(5, 2);
 

for k = 1:5
    e_fun = e_c3{k};
    
    figure('Name', ['Część 3 – p(t) – ' nazwy_c3{k}], ...
           'NumberTitle', 'off', 'Position', [50+k*20, 50+k*20, 1100, 460]);
    
    for ki = 1:2
        h  = kroki_c3(ki);
        t  = t0 : h : tk;
        N  = length(t);
 
        i1_lin = zeros(1, N);  i2_lin = zeros(1, N);  uc_lin = zeros(1, N);
        i1_nlin = zeros(1, N); i2_nlin = zeros(1, N); uc_nlin = zeros(1, N);
        
        p_mid_lin = zeros(1, N-1);
        p_mid_nlin = zeros(1, N-1);
        coeff_spline = splajn3_przygotuj(u_rown, M_rown);
 
        for j = 1:N-1
            ej = e_fun(t(j));
            ejh = e_fun(t(j) + h/2);

            i1_L = i1_lin(j); i2_L = i2_lin(j); uc_L = uc_lin(j);
            
            [di1_L, di2_L, duC_L] = pochodne(i1_L, i2_L, uc_L, ej, R1, R2, C, L1, L2, M);
            i1h_L = i1_L + (h/2) * di1_L; i2h_L = i2_L + (h/2) * di2_L; uch_L = uc_L + (h/2) * duC_L;
            
            p_mid_lin(j) = R1 * i1h_L^2 + R2 * i2h_L^2; 
            
            [di1h_L, di2h_L, duCh_L] = pochodne(i1h_L, i2h_L, uch_L, ejh, R1, R2, C, L1, L2, M);
            i1_lin(j+1) = i1_L + h * di1h_L; i2_lin(j+1) = i2_L + h * di2h_L; uc_lin(j+1) = uc_L + h * duCh_L;


            i1_NL = i1_nlin(j); i2_NL = i2_nlin(j); uc_NL = uc_nlin(j);
            
            
            uL1 = max(20, min(300, abs(ej - R1*i1_NL - uc_NL)));
            M_n = max(0.10, min(1.00, interpolacja_splajnami(uL1, coeff_spline)));
            
            [di1_NL, di2_NL, duC_NL] = pochodne(i1_NL, i2_NL, uc_NL, ej, R1, R2, C, L1, L2, M_n);
            i1h_NL = i1_NL + (h/2) * di1_NL; i2h_NL = i2_NL + (h/2) * di2_NL; uch_NL = uc_NL + (h/2) * duC_NL;
            
            
            uL1h = max(20, min(300, abs(ejh - R1*i1h_NL - uch_NL)));
            M_nh = max(0.10, min(1.00, interpolacja_splajnami(uL1h, coeff_spline)));
            
            p_mid_nlin(j) = R1 * i1h_NL^2 + R2 * i2h_NL^2; 
            
            [di1h_NL, di2h_NL, duCh_NL] = pochodne(i1h_NL, i2h_NL, uch_NL, ejh, R1, R2, C, L1, L2, M_nh);
            i1_nlin(j+1) = i1_NL + h * di1h_NL; i2_nlin(j+1) = i2_NL + h * di2h_NL; uc_nlin(j+1) = uc_NL + h * duCh_NL;
        end
 
        p_full_lin = R1 * i1_lin.^2 + R2 * i2_lin.^2;
        p_full_nlin = R1 * i1_nlin.^2 + R2 * i2_nlin.^2;
 
        E_prost_lin(k, ki)  = calka_prostokaty_c3(p_mid_lin,  h);
        E_parab_lin(k, ki)  = calka_parabol_c3(p_full_lin, h);
        
        E_prost_nlin(k, ki) = calka_prostokaty_c3(p_mid_nlin, h);
        E_parab_nlin(k, ki) = calka_parabol_c3(p_full_nlin, h);
 
        subplot(1, 2, ki);
        plot(t, p_full_lin, 'b', 'LineWidth', 1.2); hold on;
        plot(t, p_full_nlin, 'r--', 'LineWidth', 1.2);
        title(sprintf('p(t),  h = %.4f s', h), 'Interpreter', 'none');
        xlabel('t [s]');  ylabel('p [W]');  grid on;  xlim([0 tk]);
        legend('p_{Liniowa}', 'p_{Nieliniowa}', 'Location', 'best');
        hold off;
 
    end
end


sep = repmat('=', 1, 82);
lin_str = repmat('-', 1, 82);


fprintf('WYNIKI DLA M LINIOWEGO\n');
fprintf('%-24s | Prostokąty  | Prostokąty  | Parabol    | Parabol\n',    '');
fprintf('%-24s |  h₁=%.3f s |  h₂=%.3f s | h₁=%.3f s | h₂=%.3f s\n', ...
        'Wymuszenie', dt1, dt2, dt1, dt2);
fprintf('%s\n', lin_str);
for k = 1:5
    fprintf('%-24s | %11.4f | %11.4f | %10.4f | %10.4f  [J]\n', ...
        nazwy_c3{k}, E_prost_lin(k,1), E_prost_lin(k,2), E_parab_lin(k,1), E_parab_lin(k,2));
    fprintf('%-24s | %11.4f | %11.4f | %10.4f | %10.4f  [W]\n', ...
        '', E_prost_lin(k,1)/tk, E_prost_lin(k,2)/tk, E_parab_lin(k,1)/tk, E_parab_lin(k,2)/tk);
    if k<5, fprintf('%s\n', lin_str); end
end

fprintf('\nWYNIKI DLA M NIELINIOWEGO\n');
fprintf('%-24s | Prostokąty  | Prostokąty  | Parabol    | Parabol\n',    '');
fprintf('%-24s |  h₁=%.3f s |  h₂=%.3f s | h₁=%.3f s | h₂=%.3f s\n', ...
        'Wymuszenie', dt1, dt2, dt1, dt2);
fprintf('%s\n', lin_str);
for k = 1:5
    fprintf('%-24s | %11.4f | %11.4f | %10.4f | %10.4f  [J]\n', ...
        nazwy_c3{k}, E_prost_nlin(k,1), E_prost_nlin(k,2), E_parab_nlin(k,1), E_parab_nlin(k,2));
    fprintf('%-24s | %11.4f | %11.4f | %10.4f | %10.4f  [W]\n', ...
        '', E_prost_nlin(k,1)/tk, E_prost_nlin(k,2)/tk, E_parab_nlin(k,1)/tk, E_parab_nlin(k,2)/tk);
    if k<5, fprintf('%s\n', lin_str); end
end
fprintf('%s\n\n', sep);
end

function [di1, di2, duC] = pochodne(i1, i2, uc, e, R1, R2, C, L1, L2, M)
    D = L1*L2 - M^2;

    di1 = ( L2*(e - R1*i1 - uc) + M*R2*i2 ) / D;
    di2 = ( -M*(e - R1*i1 - uc) - L1*R2*i2 ) / D;
    duC = i1 / C;
end

function coeff = splajn3_przygotuj(x_nodes, y_nodes)
    n = length(x_nodes) - 1;
    h = x_nodes(2) - x_nodes(1);   % węzły równoodległe
 
    % Warunki brzegowe: zerowe pochodne pierwszego rzędu
    alpha = 0;
    beta  = 0;
 
    A = zeros(n+1);
    b = zeros(n+1, 1);
 
    A(1,1) = 4;  A(1,2) = 2;
    b(1) = y_nodes(1) + h/3 * alpha;
 
    for i = 2:n
        A(i, i-1) = 1;
        A(i, i)   = 4;
        A(i, i+1) = 1;
        b(i) = y_nodes(i);
    end
 
    A(n+1, n+1) = 4;  A(n+1, n) = 2;
    b(n+1) = y_nodes(n+1) - h/3 * beta;
 
    c = zeros(n+3, 1);
    c(2:n+2) = A \ b;
    c(1) = c(3)   - h/3 * alpha;
    c(n+3) = c(n+1) + h/3 * beta;
 
    coeff.c = c;
    coeff.h = h;
    coeff.x0 = x_nodes(1);
    coeff.n = n;
end
 
function y = phi_bspline(xi, h, x)
    % Bazowa funkcja B-spline 3. stopnia
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
    c = coeff.c;
    h = coeff.h;
    x0 = coeff.x0;
    n = coeff.n;
    M_val = 0;
    for i = -1:n+1
        xi = x0 + i * h;
        M_val = M_val + c(i+2) * phi_bspline(xi, h, u1);
    end
end

function I = calka_prostokaty_c3(p_mid_vals, h)
    m = length(p_mid_vals);   % liczba podprzedziałów
    S = 0;
    for i = 1:m
        S = S + p_mid_vals(i);
    end
    I = S * h;
end
 

function I = calka_parabol_c3(p_vals, h)
    N = length(p_vals);
    m = N - 1;
    if mod(m, 2) ~= 0           
        p_vals = p_vals(1:end-1);
        N = N - 1;
    end
    S = p_vals(1) + p_vals(end);
    for i = 2 : 2 : N-1        
        S = S + 4 * p_vals(i);
    end
    for i = 3 : 2 : N-2       
        S = S + 2 * p_vals(i);
    end
    I = h/3 * S;
end