function fpa_optimization_lora()
    %---------------- Parâmetros ----------------
    n = 100; max_iter = 200; switch_prob = 0.8;
    S_min = 7; S_max = 12;
    Nc_initial = 500; Nc_final = 4500; step_Nc = 1000;
    lambda_val = 6; b = 48;
    Toa  = [0.1048, 0.1802, 0.3211, 0.5636, 1.0485, 1.9398];
    Trx1 = [1.1048, 1.1802, 1.3211, 1.5636, 2.0485, 2.9398];
    Trx2 = [2.1048, 2.1802, 2.3211, 2.5636, 3.0485, 3.9398];
    V = 3.3; I_tx = 44; I_rx = 10.5; I_st = 1.4; I_id = 0.0015;
    RD1 = 1; RD2 = 2; T = 720;
    pesos = [1, 0; 0.75, 0.25; 0.5, 0.5; 0.25, 0.75; 0.1, 0.9];
    Nc_values = Nc_initial:step_Nc:Nc_final;
    dim = S_max - S_min + 1; % número de SFs
    num_pesos = size(pesos,1);

    % Matrizes para salvar em .mat
    vazao_fpa       = zeros(num_pesos, numel(Nc_values));
    energia_fpa     = zeros(num_pesos, numel(Nc_values));
    utility_fpa     = zeros(num_pesos, numel(Nc_values));
    EFF_fpa         = zeros(num_pesos, numel(Nc_values));
    convergencia_fpa= zeros(num_pesos, max_iter);

    %---------------- Gráficos ----------------
    figure; sgtitle('Otimização com FPA em Redes LoRa');
    subplot(4,1,1); hold on; title('Número de Nós vs Vazão Total');    xlabel('Número de Nós'); ylabel('Vazão Total (bps)'); grid on;
    subplot(4,1,2); hold on; title('Número de Nós vs Energia Total');   xlabel('Número de Nós'); ylabel('Energia Total (J)');   grid on;
    subplot(4,1,3); hold on; title('Convergência média do FPA');        xlabel('Iterações');     ylabel('Fitness');            grid on;
    subplot(4,1,4); hold on; title('Número de Nós vs Utilidade');       xlabel('Número de Nós'); ylabel('Utilidade');          grid on;

    %---------------- Loop principal ----------------
    for peso_idx = 1:num_pesos
        pesoR = pesos(peso_idx,1); pesoE = pesos(peso_idx,2);
        vazao_total = zeros(1, numel(Nc_values));
        energia_total = zeros(1, numel(Nc_values));
        utilidade_total = zeros(1, numel(Nc_values));
        EFF_total = zeros(1, numel(Nc_values));
        convergencia_media = zeros(max_iter,1);

        for Nc_idx = 1:numel(Nc_values)
            Nc = Nc_values(Nc_idx);

            % Rmax: maximiza utilidade -> p_Rmax
            [p_Rmax, conv_R, Rmax] = flower_pollination_algorithm( ...
                @(pv) sum(utilidade_de_rede(lambda_val, pv, Nc, b, Toa)), ...
                dim, n, max_iter, switch_prob, 'maximize');

            % Energia no p_Rmax (referência)
            E_ref = sum(modelo_de_energia(V, I_id, I_st, I_tx, I_rx, Trx1, Trx2, RD1, RD2, p_Rmax, Nc, T, Toa));

            % Emin: minimiza energia -> p_Emin
            [p_Emin, conv_E, Emin] = flower_pollination_algorithm( ...
                @(pv) sum(modelo_de_energia(V, I_id, I_st, I_tx, I_rx, Trx1, Trx2, RD1, RD2, pv, Nc, T, Toa)), ...
                dim, n, max_iter, switch_prob, 'minimize');

            % Utilidade no p_Emin
            Rmin = sum(utilidade_de_rede(lambda_val, p_Emin, Nc, b, Toa));

            % Intervalos (positivos)
            alfa = max(Rmax - Rmin, eps);
            beta = max(E_ref - Emin, eps);

            % Maximiza eficiência ponderada -> best_solution
            [best_solution, conv_eff, ~] = flower_pollination_algorithm( ...
                @(pv) sum(eficiencia(lambda_val, pv, Nc, b, Toa, V, I_id, I_st, I_tx, I_rx, Trx1, Trx2, RD1, RD2, T, alfa, beta, pesoR, pesoE)), ...
                dim, n, max_iter, switch_prob, 'maximize');

            % Medidas no best_solution
            vazao_total(Nc_idx)     = sum(vazao(lambda_val, best_solution, Nc, b, Toa));
            energia_total(Nc_idx)   = sum(modelo_de_energia(V, I_id, I_st, I_tx, I_rx, Trx1, Trx2, RD1, RD2, best_solution, Nc, T, Toa));
            utilidade_total(Nc_idx) = sum(utilidade_de_rede(lambda_val, best_solution, Nc, b, Toa));
            EFF_total(Nc_idx)       = sum(eficiencia(lambda_val, best_solution, Nc, b, Toa, V, I_id, I_st, I_tx, I_rx, Trx1, Trx2, RD1, RD2, T, alfa, beta, pesoR, pesoE));

            % Acumula convergência
            fprintf('Peso (%.2f, %.2f) - Convergência:\n', pesoR, pesoE);
            disp(convergencia_media');

            convergencia_media = convergencia_media + reshape(conv_eff, [], 1);
        end

        convergencia_media = convergencia_media / numel(Nc_values);

        % Armazena nas matrizes para salvar
        vazao_fpa(peso_idx,:)        = vazao_total;
        energia_fpa(peso_idx,:)      = energia_total;
        utility_fpa(peso_idx,:)      = utilidade_total;
        EFF_fpa(peso_idx,:)          = EFF_total;
        convergencia_fpa(peso_idx,:) = convergencia_media(:)';

        %---------------- Plots ----------------
        subplot(4,1,1); plot(Nc_values, vazao_total, '-o', 'DisplayName', sprintf('Peso (%.2f, %.2f)', pesoR, pesoE));
        subplot(4,1,2); plot(Nc_values, energia_total, '-o', 'DisplayName', sprintf('Peso (%.2f, %.2f)', pesoR, pesoE));
        subplot(4,1,3); plot(1:max_iter, convergencia_media, '-',  'DisplayName', sprintf('Peso (%.2f, %.2f)', pesoR, pesoE));
        subplot(4,1,4); plot(Nc_values, utilidade_total, '-o', 'DisplayName', sprintf('Peso (%.2f, %.2f)', pesoR, pesoE));
    end

    % Legendas
    subplot(4,1,1); legend('Location','best');
    subplot(4,1,2); legend('Location','best');
    subplot(4,1,3); legend('Location','best');
    subplot(4,1,4); legend('Location','best');

    %---------------- Salvamento em .mat ----------------
    pesos_fpa = pesos;
    save('fpa_vazao.mat','Nc_values','vazao_fpa','pesos_fpa');
    save('fpa_energia.mat','Nc_values','energia_fpa','pesos_fpa');
    save('fpa_utility.mat','Nc_values','utility_fpa','pesos_fpa');
    save('fpa_EFF.mat','Nc_values','EFF_fpa','pesos_fpa');
    save('fpa_convergencia.mat','convergencia_fpa','Nc_values','pesos_fpa');
end

%---------------- FPA ----------------
function [best_solution, convergencia, best_fitness] = flower_pollination_algorithm(fitness_func, num_dimensions, pop_size, max_iter, switch_prob, mode)
    % Pop inicial no simplex (p_i >= 0, somatório = 1)
    population = rand(pop_size, num_dimensions);
    population = population ./ sum(population, 2);

    fitness = arrayfun(@(i) fitness_func(population(i, :)), 1:pop_size);

    if strcmpi(mode, 'maximize')
        [best_fitness, idx] = max(fitness);
        better = @(a,b) a > b;
    else
        [best_fitness, idx] = min(fitness);
        better = @(a,b) a < b;
    end
    best_solution = population(idx, :);

    convergencia = zeros(1, max_iter);

    for iter = 1:max_iter
        for i = 1:pop_size
            if rand > switch_prob
                % Polinização global (Levy flight)
                step = levy_flight(num_dimensions);
                new_solution = population(i, :) + step .* (population(i, :) - best_solution);
            else
                % Polinização local
                epsilon = rand;
                jk = randperm(pop_size, 2);
                new_solution = population(i, :) + epsilon .* (population(jk(1), :) - population(jk(2), :));
            end

            % Projeção em [0,1] e normalização no simplex
            new_solution = max(0, min(1, new_solution));
            s = sum(new_solution);
            if s <= 0
                new_solution = ones(1, num_dimensions) / num_dimensions;
            else
                new_solution = new_solution / s;
            end

            % Avaliação
            new_fitness = fitness_func(new_solution);

            % Aceitação
            if better(new_fitness, fitness(i))
                population(i, :) = new_solution;
                fitness(i) = new_fitness;

                if better(new_fitness, best_fitness)
                    best_solution = new_solution;
                    best_fitness = new_fitness;
                end
            end
        end
        convergencia(iter) = best_fitness;
    end
end

%---------------- Métricas ----------------
function Va = vazao(lambda_val, p, Nc, b, Toa)
    Va = (lambda_val .* p .* Nc .* b) .* exp(-2 .* trafego_de_carga(lambda_val, p, Nc, Toa));
end

function En = modelo_de_energia(V, I_id, I_st, I_tx, I_rx, Trx1, Trx2, RD1, RD2, p, Nc, T, Toa)
    En = (0.5 .* p .* Nc .* V .* (Toa .* I_tx + RD1 .* I_st + Trx1 .* I_rx) + ...
          0.5 .* p .* Nc .* V .* (Toa .* I_tx + (RD2 - Trx1) .* I_st + (Trx1 + Trx2) .* I_rx)) + ...
         (0.5 .* p .* Nc .* V .* (T - (Toa + Trx1 + RD1)) .* I_id + ...
          0.5 .* p .* Nc .* V .* (T - (Toa + RD2 + Trx2)) .* I_id);
end

function Re = utilidade_de_rede(lambda_val, p, Nc, b, Toa)
    Re = log(lambda_val .* p .* Nc .* b + eps) - 2 .* trafego_de_carga(lambda_val, p, Nc, Toa);
end

function EFF = eficiencia(lambda_val, p, Nc, b, Toa, V, I_id, I_st, I_tx, I_rx, Trx1, Trx2, RD1, RD2, T, alfa, beta, pesoR, pesoE)
    Re = utilidade_de_rede(lambda_val, p, Nc, b, Toa);
    En = modelo_de_energia(V, I_id, I_st, I_tx, I_rx, Trx1, Trx2, RD1, RD2, p, Nc, T, Toa);
    EFF = (pesoR .* Re ./ alfa) - (pesoE .* En ./ beta);
end

function Gs = trafego_de_carga(lambda_val, p, Nc, Toa)
    Gs = lambda_val .* p .* Nc .* Toa ./ 10000;
end

%---------------- Levy Flight ----------------
function step = levy_flight(num_dimensions)
    beta = 1.5;
    sigma = (gamma(1 + beta) * sin(pi * beta / 2) / ...
            (gamma((1 + beta) / 2) * beta * 2^((beta - 1) / 2)))^(1 / beta);
    u = randn(1, num_dimensions) * sigma;
    v = randn(1, num_dimensions);
    step = u ./ (abs(v).^(1 / beta) + eps);
end
