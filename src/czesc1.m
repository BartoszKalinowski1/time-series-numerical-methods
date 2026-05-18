function czesc1
clear; close all;
R1 = 0.1;
R2 = 10;
C = 0.5;
L1 = 3;
L2 = 5;
M = 0.8;

t0 = 0;
tk = 30;
h_wartosci = [0.1, 0.5, 0.1, 0.009];

e1 = @(x) (mod(x, 3) < 1.5) * 120;
e2 = @(x) (240 * sin(x));
e3 = @(x) (210 * sin(2 * pi * 5 * x));
e4 = @(x) (120 * sin(2 * pi * 50 * x));

wymuszenia = {e1, e2, e3, e4};
nazwy = {'Prostokątne (T=3s, 120V)', '240sin(t)', '210sin(2\pi\cdot5t)', '120sin(2\pi\cdot50t)'};

D1 = L1/M - M/L2;
D2 = M/L1 - L2/M;

for k = 1:4
    h = h_wartosci(k);
    t = t0:h:tk;
    N = length(t);

    e = wymuszenia{k};

    i1_euler = zeros(1, N);
    i2_euler = zeros(1, N);
    uc_euler = zeros(1, N);

    i1_euler_ulepszony = zeros(1, N);
    i2_euler_ulepszony = zeros(1, N);
    uc_euler_ulepszony = zeros(1, N);

    for j = 1:N-1
        i1 = i1_euler(j);
        i2 = i2_euler(j);
        uc = uc_euler(j);
        ej = e(t(j));

        di1_dt = (1/D1) * (-R1/M * i1 + R2/L2 * i2 - 1/M * uc + 1/M * ej);
        di2_dt = (1/D2) * (-R1/L1 * i1 + R2/M * i2 - 1/L1 * uc + 1/L1 * ej);
        duC_dt = (1/C) * i1;

        i1_euler(j+1) = i1 + h * di1_dt;
        i2_euler(j+1) = i2 + h * di2_dt;
        uc_euler(j+1) = uc + h * duC_dt;
    end
    for j = 1:N-1
        i1 = i1_euler_ulepszony(j);
        i2 = i2_euler_ulepszony(j);
        uc = uc_euler_ulepszony(j);
        ej = e(t(j));

        di1_dt = (1/D1) * (-R1/M * i1 + R2/L2 * i2 - 1/M * uc + 1/M * ej);
        di2_dt = (1/D2) * (-R1/L1 * i1 + R2/M * i2 - 1/L1 * uc + 1/L1 * ej);
        duC_dt = (1/C) * i1;

        i1_half = i1 + (h/2) * di1_dt;
        i2_half = i2 + (h/2) * di2_dt;
        uc_half = uc + (h/2) * duC_dt;
        ej_half = e(t(j) + h/2);

        di1_dt_half = (1/D1) * (-R1/M * i1_half + R2/L2 * i2_half - 1/M * uc_half + 1/M * ej_half);
        di2_dt_half = (1/D2) * (-R1/L1 * i1_half + R2/M * i2_half - 1/L1 * uc_half + 1/L1 * ej_half);
        duC_dt_half = (1/C) * i1_half;

        i1_euler_ulepszony(j+1) = i1 + h * di1_dt_half;
        i2_euler_ulepszony(j+1) = i2 + h * di2_dt_half;
        uc_euler_ulepszony(j+1) = uc + h * duC_dt_half;
    end
         
    e_plot = arrayfun(e, t);

    figure('Name', ['Euler – ' nazwy{k}], 'NumberTitle', 'off');

    % i1 i i2 razem
    subplot(3,1,1);
    plot(t, i1_euler, 'b', 'LineWidth', 1.3); hold on;
    plot(t, i2_euler, 'r', 'LineWidth', 1.3);
    hold off;
    title(['Prądy i_1(t), i_2(t) – Euler – ' nazwy{k}]);
    xlabel('t [s]'); ylabel('i [A]');
    legend('i_1','i_2','Location','best');
    grid on;

    % uC
    subplot(3,1,2);
    plot(t, uc_euler, 'm', 'LineWidth', 1.3);
    title(['Napięcie u_C(t) – Euler – ' nazwy{k}]);
    xlabel('t [s]'); ylabel('u_C [V]');
    grid on;

    % e(t)
    subplot(3,1,3);
    plot(t, e_plot, 'k', 'LineWidth', 1.3);
    title(['Sygnał wymuszający e(t) – ' nazwy{k}]);
    xlabel('t [s]'); ylabel('e [V]');
    grid on;

    figure('Name', ['Ulepszony Euler – ' nazwy{k}], 'NumberTitle', 'off');

    % i1 i i2 razem
    subplot(3,1,1);
    plot(t, i1_euler_ulepszony, '--b', 'LineWidth', 1.3); hold on;
    plot(t, i2_euler_ulepszony, '--r', 'LineWidth', 1.3);
    hold off;
    title(['Prądy i_1(t), i_2(t) – Ulepszony Euler – ' nazwy{k}]);
    xlabel('t [s]'); ylabel('i [A]');
    legend('i_1','i_2','Location','best');
    grid on;

    % uC
    subplot(3,1,2);
    plot(t, uc_euler_ulepszony, '--m', 'LineWidth', 1.3);
    title(['Napięcie u_C(t) – Ulepszony Euler – ' nazwy{k}]);
    xlabel('t [s]'); ylabel('u_C [V]');
    grid on;

    % e(t)
    subplot(3,1,3);
    plot(t, e_plot, 'k', 'LineWidth', 1.3);
    title(['Sygnał wymuszający e(t) – ' nazwy{k}]);
    xlabel('t [s]'); ylabel('e [V]');
    grid on;

end 
end