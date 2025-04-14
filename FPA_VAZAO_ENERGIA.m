function [X_all, Y_all] = fpa_optimization_lora(pesoR, pesoE, Nc_values)
    % Parâmetros principais
    n = 100; % Tamanho da população
    max_iter = 200; % Iterações máximas
    p = 0.8; % Probabilidade de polinização global
    S_min = 7;
    S_max = 12;
    Nc_initial = 100;
    Nc_final = 1000;
    step_Nc = 100; % Incremento no número de nós
    lambda_val = 6;
    b = 48;
    Toa = [0.1048, 0.1802, 0.3211, 0.5636, 1.0485, 1.9398];
    Trx1 = [1.1048, 1.1802, 1.3211, 1.5636, 2.0485, 2.9398];
    Trx2 = [2.1048, 2.1802, 2.3211, 2.5636, 3.0485, 3.9398];
    V = 3.3; I_tx = 44; I_rx = 10.5; I_st = 1.4; I_id = 0.0015;
    RD1 = 1; RD2 = 2; T = 720;
  

    % Inicializa os dados para cada combinação de peso
    X_all = zeros(1, length(Nc_values)); % Vazão (Throughput)
    Y_all = zeros(1, length(Nc_values)); % Energia

    % Executa a otimização para os diferentes valores de Nc
    for Nc_idx = 1:length(Nc_values)
        Nc = Nc_values(Nc_idx);

        % Otimização de vazão máxima e energia mínima
        [Rmax, ~] = flower_pollination_algorithm( ...
            @(p) sum(utilidade_de_rede(lambda_val, p, Nc, b, Toa)), ...
            S_max - S_min + 1, n, max_iter, p, 'maximize');
        Emin = modelo_de_energia(V, I_id, I_st, I_tx, I_rx, Trx1, Trx2, RD1, RD2, p, Nc, T, Toa);

        [Emax, ~] = flower_pollination_algorithm( ...
            @(p) sum(modelo_de_energia(V, I_id, I_st, I_tx, I_rx, Trx1, Trx2, RD1, RD2, p, Nc, T, Toa)), ...
            S_max - S_min + 1, n, max_iter, p, 'maximize');
        Rmin = utilidade_de_rede(lambda_val, p, Nc, b, Toa);

        alfa = sum(Rmax - Rmin); % Normalizador para vazão
        beta = sum(Emin - Emax); % Normalizador para energia

        % Otimização conjunta
        [best_solution, ~] = flower_pollination_algorithm( ...
            @(p) sum(eficiencia(lambda_val, p, Nc, b, Toa, V, I_id, I_st, I_tx, I_rx, Trx1, Trx2, RD1, RD2, T, alfa, beta, pesoR, pesoE)), ...
            S_max - S_min + 1, n, max_iter, p, 'maximize');

        % Calculando a vazão e energia para a solução otimizada
        vazao_atual = sum(vazao(lambda_val, best_solution, Nc, b, Toa));
        energia_atual = sum(modelo_de_energia(V, I_id, I_st, I_tx, I_rx, Trx1, Trx2, RD1, RD2, best_solution, Nc, T, Toa));

        % Armazenando os resultados
        X_all(Nc_idx) = vazao_atual;
        Y_all(Nc_idx) = energia_atual;
    end
end

% Funções Auxiliares 
function [best_fitness, convergencia] = flower_pollination_algorithm(fitness_func, num_dimensions, pop_size, max_iter, p, mode)
    % Inicialização
    population = rand(pop_size, num_dimensions);
    population = population ./ sum(population, 2);
    fitness = arrayfun(@(i) fitness_func(population(i, :)), 1:pop_size);

    if strcmp(mode, 'maximize')
        [best_fitness, idx] = max(fitness);
    else
        [best_fitness, idx] = min(fitness);
    end
    best_solution = population(idx, :);

    convergencia = zeros(1, max_iter);

    for iter = 1:max_iter
        for i = 1:pop_size
            if rand > p
                step = levy_flight(num_dimensions);
                new_solution = population(i, :) + step .* (population(i, :) - best_solution);
            else
                epsilon = rand;
                jk = randperm(pop_size, 2);
                new_solution = population(i, :) + epsilon .* (population(jk(1), :) - population(jk(2), :));
            end

            new_solution = max(0, min(1, new_solution));
            new_solution = new_solution / sum(new_solution);

            new_fitness = fitness_func(new_solution);

            if (strcmp(mode, 'maximize') && new_fitness > fitness(i)) || ...
               (strcmp(mode, 'minimize') && new_fitness < fitness(i))
                population(i, :) = new_solution;
                fitness(i) = new_fitness;

                if (strcmp(mode, 'maximize') && new_fitness > best_fitness) || ...
                   (strcmp(mode, 'minimize') && new_fitness < best_fitness)
                    best_solution = new_solution;
                    best_fitness = new_fitness;
                end
            end
        end
        convergencia(iter) = best_fitness;
    end
end

function Va = vazao(lambda_val, p, Nc, b, Toa)
    Va = (lambda_val * p * Nc * b .* exp(-2 * trafego_de_carga(lambda_val, p, Nc, Toa)));
end

function En = modelo_de_energia(V, I_id, I_st, I_tx, I_rx, Trx1, Trx2, RD1, RD2, p, Nc, T, Toa)
    En = (0.5 .* p .* Nc .* V .* (Toa .* I_tx + RD1 .* I_st + Trx1 .* I_rx) + ...
           0.5 .* p .* Nc .* V .* (Toa .* I_tx + (RD2 - Trx1) .* I_st + (Trx1 + Trx2) .* I_rx)) + ...
          (0.5 .* p .* Nc .* V .* (T - (Toa + Trx1 + RD1)) .* I_id + ...
           0.5 .* p .* Nc .* V .* (T - (Toa + RD2 + Trx2)) .* I_id);
end

function Re = utilidade_de_rede(lambda_val, p, Nc, b, Toa)
    Re = log(lambda_val .* p .* Nc .* b) - 2 .* trafego_de_carga(lambda_val, p, Nc, Toa);
end

function EFF = eficiencia(lambda_val, p, Nc, b, Toa, V, I_id, I_st, I_tx, I_rx, Trx1, Trx2, RD1, RD2, T, alfa, beta, pesoR, pesoE)
    EFF = (pesoR ./ alfa .* utilidade_de_rede(lambda_val, p, Nc, b, Toa)) - ...
          (pesoE ./ beta .* modelo_de_energia(V, I_id, I_st, I_tx, I_rx, Trx1, Trx2, RD1, RD2, p, Nc, T, Toa));
end

function Gs = trafego_de_carga(lambda_val, p, Nc, Toa)
    Gs = lambda_val .* p .* Nc .* Toa ./ 10000;
end
function step = levy_flight(num_dimensions)
    beta = 1.5;
    sigma = (gamma(1 + beta) * sin(pi * beta / 2) / ...
            (gamma((1 + beta) / 2) * beta * 2^((beta - 1) / 2)))^(1 / beta);
    u = randn(1, num_dimensions) * sigma;
    v = randn(1, num_dimensions);
    step = u ./ abs(v).^(1 / beta);
end

