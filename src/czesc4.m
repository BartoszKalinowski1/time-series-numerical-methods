function czesc4
clear; close all;
R1 = 0.1;
R2 = 10;
C = 0.5;
L1 = 3;
L2 = 5;
M = 0.8;
t0 = 0;
tk = 30;

d_tol = 1e-4;   
max_iter = 50;


F_celu = @(f) oblicz_f(abs(f), R1, R2, C, L1, L2, M, t0, tk, 0.001) - 406;

f_skan = 0.1:0.2:20;
a_start = f_skan(1); 
b_start = f_skan(end);
znalaziono = false;

licznik_wywolan = 1; 
F_poprzednie = F_celu(f_skan(1));
for i = 2:length(f_skan)
    F_aktualne = F_celu(f_skan(i));
    licznik_wywolan = licznik_wywolan + 1;
    if F_poprzednie * F_aktualne <= 0
        a_start = f_skan(i-1);
        b_start = f_skan(i);
        znalaziono = true;
        break;
    end
    F_poprzednie = F_aktualne;
end

if ~znalaziono
    warning('Nie znaleziono zmiany znaku w przedziale 0.1-20 Hz.');
else
    fprintf('Znaleziono zmianę znaku w przedziale [%d Hz, %d Hz].', a_start, b_start);
end

wyniki_f = zeros(1, 3);
wyniki_iter = zeros(1, 3);
wyniki_wywolan = zeros(1, 3);

% bisekcja
licznik_wywolan = 2; 
a = a_start; 
b = b_start;
iter_bisekcja = 0;
f_wynik_bis = 0;

Fa = F_celu(a);
Fb = F_celu(b);

if Fa * Fb > 0
    warning('Brak zmiany znaku w przedziale [a, b] dla bisekcji');
else
    for i = 1:max_iter
        iter_bisekcja = iter_bisekcja + 1;
        f_mid = (a + b) / 2;
        F_mid = F_celu(f_mid);
        licznik_wywolan = licznik_wywolan + 1;
        if (Fa * F_mid < 0)
            b = f_mid;
            Fb = F_mid;
        else
            a = f_mid;
            Fa = F_mid;
        end
        
        if abs(F_mid) < d_tol || (b - a)/2 < d_tol
            f_wynik_bis = f_mid;
            break;
        end
    end
end
wyniki_f(1) = f_wynik_bis;
wyniki_iter(1) = iter_bisekcja;
wyniki_wywolan(1) = licznik_wywolan;


% sieczne
licznik_wywolan = 2;
x = [a_start, b_start];
iter_sieczne = 1;
n = 2;

Fx0 = F_celu(x(n-1));
Fx1 = F_celu(x(n));

while abs(Fx1) > d_tol && iter_sieczne < max_iter
    x_new = x(n) - Fx1 * (x(n) - x(n-1)) / (Fx1 - Fx0);
    x = [x, x_new];
    n = n + 1;
    
    Fx0 = Fx1;
    Fx1 = F_celu(x(n));
    licznik_wywolan = licznik_wywolan + 1;
    iter_sieczne = iter_sieczne + 1;
end
wyniki_f(2) = x(end);
wyniki_iter(2) = iter_sieczne;
wyniki_wywolan(2) = licznik_wywolan;

%{
% Sprawdzanie warunku dokładności przybliżenia pochodnej
f_test = (a_start + b_start) / 2; 
F_test = F_celu(f_test);

% różne wielkości kroku delta f
df_wektor = [1, 0.5, 0.1, 0.01, 0.001, 1e-4, 1e-5];

for i = 1:length(df_wektor)
    df = df_wektor(i);
    df_polowa = df / 2;
    
    F_df = F_celu(f_test + df);
    pochodna_1 = (F_df - F_test) / df;
    
    F_df2 = F_celu(f_test + df_polowa);
    pochodna_2 = (F_df2 - F_test) / df_polowa;
    
    blad_procentowy = abs((pochodna_1 - pochodna_2) / pochodna_1) * 100;
    
    fprintf('%10.1e | %15.4f | %15.4f | %9.4f %%\n', ...
            df, pochodna_1, pochodna_2, blad_procentowy);
end
fprintf('%s\n\n', repmat('-', 1, 65));
%}
% quasi - newton
licznik_wywolan = 0;
f_start = (a_start + b_start) / 2; 
delta_f = 1e-4; 
iter_newton = 0;

x_n = f_start;

for i = 1:max_iter
    iter_newton = iter_newton + 1;
    F_i = F_celu(x_n);
    licznik_wywolan = licznik_wywolan + 1;
    if abs(F_i) < d_tol
        break;
    end
    
    F_i_delta = F_celu(x_n + delta_f);
    licznik_wywolan = licznik_wywolan + 1;
    pochodna = (F_i_delta - F_i) / delta_f;
    
    x_n_new = x_n - F_i / pochodna;
    
    if abs(x_n_new - x_n) < d_tol
        x_n = x_n_new;
        break;
    end
    x_n = x_n_new;
end
wyniki_f(3) = x_n;
wyniki_iter(3) = iter_newton;
wyniki_wywolan(3) = licznik_wywolan;

sep = repmat('=', 1, 80);
fprintf('\n%s\n', sep);
fprintf('%-20s | %-15s | %-15s | %-20s\n', 'Metoda', 'Częstotliwość f', 'Liczba iteracji', 'Wywołań funkcji P(f)');
fprintf('%s\n', repmat('-', 1, 80));
nazwy_metod = {'Bisekcji', 'Siecznych', 'Quasi-Newtona'};
for m = 1:3
    fprintf('%-20s | %12.4f Hz | %15d | %20d\n', ...
        nazwy_metod{m}, wyniki_f(m), wyniki_iter(m), wyniki_wywolan(m));
end
fprintf('%s\n\n', sep);
end

function [di1, di2, duC] = pochodne(i1, i2, uc, e, R1, R2, C, L1, L2, M)
    D = L1*L2 - M^2;

    di1 = ( L2*(e - R1*i1 - uc) + M*R2*i2 ) / D;
    di2 = ( -M*(e - R1*i1 - uc) - L1*R2*i2 ) / D;
    duC = i1 / C;
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

function F_val = oblicz_f(f, R1, R2, C, L1, L2, M, t0, tk, h)
    P_srednia = oblicz_moc_srednia(f, R1, R2, C, L1, L2, M, t0, tk, h);
    F_val = P_srednia;
end

function P_avg = oblicz_moc_srednia(f, R1, R2, C, L1, L2, M, t0, tk, h)
    t = t0:h:tk;
    N = length(t);
    
    i1 = 0; i2 = 0; uc = 0;
    p_full = zeros(1, N);
    
    e_fun = @(x) 210 * sin(2 * pi * abs(f) * x);
    
    for j = 1:N-1
        ej = e_fun(t(j));
        ejh = e_fun(t(j) + h/2);
        
        [di1, di2, duC] = pochodne(i1, i2, uc, ej, R1, R2, C, L1, L2, M);
        
        i1h = i1 + (h/2) * di1;
        i2h = i2 + (h/2) * di2;
        uch = uc + (h/2) * duC;
        
        [di1h, di2h, duCh] = pochodne(i1h, i2h, uch, ejh, R1, R2, C, L1, L2, M);
        
        i1 = i1 + h * di1h;
        i2 = i2 + h * di2h;
        uc = uc + h * duCh;
        
        p_full(j+1) = R1 * i1^2 + R2 * i2^2;
    end
    
    E = calka_parabol_c3(p_full, h);
    P_avg = E / (tk - t0);
end