function [best, fmin] = vafpa_ag(n, d, iter, p)
    clc;
    clear all;
    close all;
    tic; % Início da medição do tempo de execução
    rng('shuffle');
    if nargin < 1
        n = 100; % Número de soluções
        d = 400; % Número de nós
        iter = 1000; % Interações
        p = 0.8; % Taxa de Mutação
    end

    % Parâmetros
    s_min = 7;
    s_max = 12;
    
    bw = 125000;
    cr = 1;
    h = 0;
    npream = 8;
    pl = [230, 230, 123, 59, 59, 59];
    b = [50, 75, 100, 125, 150, 175];
    toa = [0.1048, 0.1802, 0.3211, 0.5636, 1.0485, 1.9398];
    lambda_val = [343.32, 199.75, 112.1, 63.87, 34.33, 18.55];
    num_nodes_range = 100:100:4000; % Faixa de número de nós
    vazao_total = zeros(length(num_nodes_range), 1);
    proporcoes_otimas = zeros(length(num_nodes_range), s_max - s_min + 1);
    
    % Variáveis para os gráficos
    convergence = zeros(iter, 1);
    diversity = zeros(iter, 1);
    mutation_rate = zeros(iter, 1);
    crossover_rate = zeros(iter, 1);

    for idx = 1:length(num_nodes_range)
        Nc = num_nodes_range(idx);
        lb = zeros(1, s_max - s_min + 1);
        ub = ones(1, s_max - s_min + 1);
        sol = rand(n, s_max - s_min + 1);
        sol = sol ./ sum(sol, 2); % Normalizar para somar 1

        fitness = arrayfun(@(x) network_utility(sol(x, :), lambda_val, Nc, bw, cr, h, npream, pl, b, toa, s_min), 1:n);

        [fmin, best_idx] = max(fitness);
        best = sol(best_idx, :);
        s = sol;

        % Algoritmo Genético (AG)
        for t = 1:iter
            new_pop = zeros(size(s));
            mutation_count = 0;
            crossover_count = 0;
            
            for i = 1:2:n % Seleção e Cruzamento
                % Seleção por torneio
                [parent1, parent2] = tournament_selection(s, fitness, 2);
                
                % Cruzamento (crossover)
                [child1, child2, crossover_done] = crossover(parent1, parent2, cr);
                if crossover_done
                    crossover_count = crossover_count + 1;
                end
                
                % Mutação
                child1 = mutation(child1, p, lb, ub);
                child2 = mutation(child2, p, lb, ub);
                if any(child1 ~= parent1) || any(child2 ~= parent2)
                    mutation_count = mutation_count + 1;
                end
                 
                % Garantir que as proporções somem para 1
                child1 = simple_bounds(child1, lb, ub);
                child1 = child1 / sum(child1); % Normaliza para que a soma seja 1
                child2 = simple_bounds(child2, lb, ub);
                child2 = child2 / sum(child2); % Normaliza para que a soma seja 1
                % Atualização da nova população
                new_pop(i, :) = child1;
                new_pop(i+1, :) = child2;
            end
            
            % Avaliar a população nova
            new_fitness = arrayfun(@(x) network_utility(new_pop(x, :), lambda_val, Nc, bw, cr, h, npream, pl, b, toa, s_min), 1:n);

            % Substituição das piores soluções
            for i = 1:n
                if new_fitness(i) > fitness(i)
                    s(i, :) = new_pop(i, :);
                    fitness(i) = new_fitness(i);
                end
            end
            
            % Melhor solução
            [fmin, best_idx] = max(fitness);
            best = s(best_idx, :);

            % Convergência (fitness)
            convergence(t) = fmin;

            % Diversidade (variabilidade das soluções)
            diversity(t) = std(fitness);

            % Taxas de mutação e cruzamento
            mutation_rate(t) = mutation_count / n;
            crossover_rate(t) = crossover_count / n;
        end

        final_node_proportion = best;
        proporcoes_otimas(idx, :) = final_node_proportion;

        G = lambda_val * Nc .* final_node_proportion .* toa;
        T = lambda_val .* final_node_proportion .* Nc .* b .* exp(-2 .* G/10000);
        vazao_total(idx) = sum(T);

        fprintf('Número de Nós: %d\n', Nc);
        for sf = s_min:s_max
            fprintf('  SF %d: P(s) = %.4f\n', sf, final_node_proportion(sf - s_min + 1));
        end
    end

    % Plot da relação entre o número de nós e a vazão total
    figure;
    plot(num_nodes_range, vazao_total, '-o', 'LineWidth', 2);
    title('Relação entre Número de Nós e Vazão Total para o AG');
    xlabel('Número de Nós');
    ylabel('Vazão Total');
    grid on;

    % Gráfico de Convergência (Fitness ao longo das iterações)
    figure;
    plot(1:iter, convergence, '-o', 'LineWidth', 2);
    title('Convergência do Algoritmo Genético');
    xlabel('Iterações');
    ylabel('Fitness');
    grid on;

    % Gráfico de Diversidade ao longo das iterações
    figure;
    plot(1:iter, diversity, '-o', 'LineWidth', 2);
    title('Controle de Diversidade');
    xlabel('Iterações');
    ylabel('Desvio Padrão da Fitness');
    grid on;

    % Gráfico da Taxa de Mutação
    %figure;
    %plot(1:iter, mutation_rate, '-o', 'LineWidth', 2);
    %title('Taxa de Mutação');
    %xlabel('Iterações');
    %ylabel('Taxa de Mutação');
    %grid on;

    % Gráfico da Taxa de Cruzamento
    %figure;
    %plot(1:iter, crossover_rate, '-o', 'LineWidth', 2);
    %title('Taxa de Cruzamento');
    %xlabel('Iterações');
    %ylabel('Taxa de Cruzamento');
    %grid on;

    toc; % Fim da medição do tempo de execução

    vazao_totalag = vazao_total;
    num_nodes_rangeag = num_nodes_range;
    save('vazaoag3.mat', 'num_nodes_rangeag', 'vazao_totalag');

    load("vazaocvx3.mat");

    figure;
    plot(num_nodes_rangeag, vazao_totalag, 'b-o', 'LineWidth', 2);
    hold on;
    plot(num_nodes_range, vazao_total, 'k-o', 'LineWidth', 2);
    title('Relação entre Número de Nós e Vazão Total para o AG e CVX');
    xlabel('Número de Nós Cobertos');
    ylabel('Vazão Total');
    legend('Vazão AG', 'Vazão CVX');
    grid on;
end

% Funções auxiliares do AG:

% Seleção por Torneio
function [parent1, parent2] = tournament_selection(pop, fitness, tournament_size)
    n = size(pop, 1);
    idx1 = randi([1, n], [tournament_size, 1]);
    idx2 = randi([1, n], [tournament_size, 1]);
    [~, best_idx1] = max(fitness(idx1));
    [~, best_idx2] = max(fitness(idx2));
    parent1 = pop(idx1(best_idx1), :);
    parent2 = pop(idx2(best_idx2), :);
end

% Cruzamento (Crossover)
function [child1, child2, crossover_done] = crossover(parent1, parent2, cr)
    crossover_done = false;
    if rand < cr
        crossover_point = randi(length(parent1));
        child1 = [parent1(1:crossover_point), parent2(crossover_point+1:end)];
        child2 = [parent2(1:crossover_point), parent1(crossover_point+1:end)];
        crossover_done = true;
    else
        child1 = parent1;
        child2 = parent2;
    end
end

% Mutação
function child = mutation(child, p, lb, ub)
    if rand < p
        mutation_point = randi(length(child));
        child(mutation_point) = rand * (ub(mutation_point) - lb(mutation_point)) + lb(mutation_point);
    end
    child = simple_bounds(child, lb, ub); % Garantir que a solução esteja dentro dos limites
end

% Função de Limite Simples
function s = simple_bounds(s, lb, ub)
    s = max(min(s, ub), lb);
end

% Função de utilidade de rede
function z = network_utility(proportions, lambda_val, nc, bw, cr, h, npream, pl, b, toa, s_min)
    gs = lambda_val .* nc .* toa .* proportions;
    z = sum(log(lambda_val .* b .* nc .* proportions) - 2 * gs);
end
