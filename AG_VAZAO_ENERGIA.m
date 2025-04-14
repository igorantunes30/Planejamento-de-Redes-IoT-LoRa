function ga_optimization_lora()
    % Parâmetros principais
    n = 100; max_iter = 200;
    S_min = 7; S_max = 12;
    Nc_initial = 500; Nc_final = 4000; step_Nc = 100;
    lambda_val = 6; b = 48;
    Toa = [0.1048, 0.1802, 0.3211, 0.5636, 1.0485, 1.9398];
    Trx1 = [1.1048, 1.1802, 1.3211, 1.5636, 2.0485, 2.9398];
    Trx2 = [2.1048, 2.1802, 2.3211, 2.5636, 3.0485, 3.9398];
    V = 3.3; I_tx = 44; I_rx = 10.5; I_st = 1.4; I_id = 0.0015;
    RD1 = 1; RD2 = 2; T = 720;
    pesos = [1, 0; 0.75, 0.25; 0.5, 0.5; 0.25, 0.75; 0.1, 0.9];
    Nc_values = Nc_initial:step_Nc:Nc_final;

    figure;
    sgtitle('Otimização com GA em Redes LoRa');

    subplot(3,1,1); hold on; title('Vazão Total'); xlabel('Número de Nós'); ylabel('Bps'); grid on;
    subplot(3,1,2); hold on; title('Energia Total'); xlabel('Número de Nós'); ylabel('Joules'); grid on;
    subplot(3,1,3); hold on; title('Convergência do GA'); xlabel('Iterações'); ylabel('Fitness'); grid on;

    for peso_idx = 1:size(pesos, 1)
        pesoR = pesos(peso_idx, 1); pesoE = pesos(peso_idx, 2);
        vazao_total = []; energia_total = []; convergencia_media = zeros(max_iter, 1);

        for Nc = Nc_values
            [p_Rmax, ~] = genetic_algorithm(@(p) sum(utilidade(lambda_val, p, Nc, b, Toa)), S_max-S_min+1, n, max_iter, 'maximize');
            Rmax = utilidade(lambda_val, p_Rmax, Nc, b, Toa);
            Emin = energia(V, I_id, I_st, I_tx, I_rx, Trx1, Trx2, RD1, RD2, p_Rmax, Nc, T, Toa);
            [p_Emax, ~] = genetic_algorithm(@(p) sum(energia(V, I_id, I_st, I_tx, I_rx, Trx1, Trx2, RD1, RD2, p, Nc, T, Toa)), S_max-S_min+1, n, max_iter, 'maximize');
            Emax = energia(V, I_id, I_st, I_tx, I_rx, Trx1, Trx2, RD1, RD2, p_Emax, Nc, T, Toa);
            Rmin = utilidade(lambda_val, p_Emax, Nc, b, Toa);
            alfa = sum(Rmax - Rmin); 
            beta = sum(Emin - Emax);

            [best_solution, convergencia] = genetic_algorithm(@(p) sum(eficiencia(lambda_val, p, Nc, b, Toa, V, I_id, I_st, I_tx, I_rx, Trx1, Trx2, RD1, RD2, T, alfa, beta, pesoR, pesoE)), S_max-S_min+1, n, max_iter, 'maximize');
            vazao_total = [vazao_total, sum(vazao(lambda_val, best_solution, Nc, b, Toa))];
            energia_total = [energia_total, sum(energia(V, I_id, I_st, I_tx, I_rx, Trx1, Trx2, RD1, RD2, best_solution, Nc, T, Toa))];
            convergencia_media = convergencia_media + convergencia(:);
        end

        convergencia_media = convergencia_media / length(Nc_values);
        subplot(3,1,1); plot(Nc_values, vazao_total, '-o', 'DisplayName', sprintf('Peso (%.2f, %.2f)', pesoR, pesoE));
        subplot(3,1,2); plot(Nc_values, energia_total, '-o', 'DisplayName', sprintf('Peso (%.2f, %.2f)', pesoR, pesoE));
        subplot(3,1,3); plot(1:max_iter, convergencia_media, 'DisplayName', sprintf('Peso (%.2f, %.2f)', pesoR, pesoE));
    end

    subplot(3,1,1); legend;
    subplot(3,1,2); legend;
    subplot(3,1,3); legend;

end

function [best_solution, convergencia] = genetic_algorithm(fitness_func, num_dimensions, pop_size, max_iter, mode)
    population = rand(pop_size, num_dimensions);
    population = population ./ sum(population, 2);
    fitness = arrayfun(@(i) fitness_func(population(i, :)), 1:pop_size);
    if strcmp(mode, 'maximize'); [best_fitness, idx] = max(fitness); else; [best_fitness, idx] = min(fitness); end
    best_solution = population(idx, :);
    convergencia = zeros(1, max_iter);

    for iter = 1:max_iter
        [~, sorted_idx] = sort(fitness, 'descend');
        parents = population(sorted_idx(1:round(pop_size/2)), :);
        offspring = zeros(pop_size - size(parents, 1), num_dimensions);
        for k = 1:size(offspring,1)
            p1 = parents(randi(size(parents,1)), :);
            p2 = parents(randi(size(parents,1)), :);
            alpha = rand;
            offspring(k,:) = alpha*p1 + (1-alpha)*p2;
        end
        mutation_rate = 0.1;
        for k = 1:size(offspring,1)
            for d = 1:num_dimensions
                if rand < mutation_rate
                    offspring(k,d) = rand;
                end
            end
            offspring(k,:) = offspring(k,:) / sum(offspring(k,:));
        end
        population = [parents; offspring];
        fitness = arrayfun(@(i) fitness_func(population(i, :)), 1:pop_size);
        if strcmp(mode, 'maximize'); [current_best_fitness, idx] = max(fitness); else; [current_best_fitness, idx] = min(fitness); end
        if (strcmp(mode, 'maximize') && current_best_fitness > best_fitness) || (strcmp(mode, 'minimize') && current_best_fitness < best_fitness)
            best_fitness = current_best_fitness;
            best_solution = population(idx, :);
        end
        convergencia(iter) = best_fitness;
    end
end

function Va = vazao(lambda_val, p, Nc, b, Toa)
    Va = (lambda_val * p * Nc * b .* exp(-2 * trafego(lambda_val, p, Nc, Toa)));
end

function En = energia(V, I_id, I_st, I_tx, I_rx, Trx1, Trx2, RD1, RD2, p, Nc, T, Toa)
    En = (0.5 .* p .* Nc .* V .* (Toa .* I_tx + RD1 .* I_st + Trx1 .* I_rx) + ...
           0.5 .* p .* Nc .* V .* (Toa .* I_tx + (RD2 - Trx1) .* I_st + (Trx1 + Trx2) .* I_rx)) + ...
          (0.5 .* p .* Nc .* V .* (T - (Toa + Trx1 + RD1)) .* I_id + ...
           0.5 .* p .* Nc .* V .* (T - (Toa + RD2 + Trx2)) .* I_id);
end

function Re = utilidade(lambda_val, p, Nc, b, Toa)
    Re = log(lambda_val .* p .* Nc .* b) - 2 .* trafego(lambda_val, p, Nc, Toa);
end

function EFF = eficiencia(lambda_val, p, Nc, b, Toa, V, I_id, I_st, I_tx, I_rx, Trx1, Trx2, RD1, RD2, T, alfa, beta, pesoR, pesoE)
    EFF = (pesoR ./ alfa .* utilidade(lambda_val, p, Nc, b, Toa)) - ...
          (pesoE ./ beta .* energia(V, I_id, I_st, I_tx, I_rx, Trx1, Trx2, RD1, RD2, p, Nc, T, Toa));
end

function Gs = trafego(lambda_val, p, Nc, Toa)
    Gs = lambda_val .* p .* Nc .* Toa ./ 10000;
end




