function ag_optimization_lora()
    %============================== PARÂMETROS GERAIS ==============================================
    rng(1); % reprodutibilidade
    pop_size       = 80;
    generations    = 200; %<--- interações
    crossover_rate = 0.9;
    mutation_rate  = 0.2;
    sigma_mut      = 0.10; % desvio da mutação gaussiana no simplex
    tournament_k   = 3;

    %============================== PARÂMETROS DO PROBLEMA =========================================
    S_min = 7; S_max = 12;
    Nc_initial = 500; Nc_final = 4500; step_Nc = 1000;   % manter comparável ao CVX
    lambda_val = 6; b = 48;
    Toa  = [0.1048, 0.1802, 0.3211, 0.5636, 1.0485, 1.9398];   % vetores linha
    Trx1 = [1.1048, 1.1802, 1.3211, 1.5636, 2.0485, 2.9398];
    Trx2 = [2.1048, 2.1802, 2.3211, 2.5636, 3.0485, 3.9398];
    V = 3.3; I_tx = 44; I_rx = 10.5; I_st = 1.4; I_id = 0.0015;
    RD1 = 1; RD2 = 2; T = 720;
    dim = S_max - S_min + 1;

    pesos  = [1, 0; 0.75, 0.25; 0.5, 0.5; 0.25, 0.75; 0.1, 0.9];
    labels = {'a=1, b=0','a=0.75, b=0.25','a=0.5, b=0.5','a=0.25, b=0.75','a=0.1, b=0.9',};
    colors = {'r','b','g','m','k'};

    Nc_values = Nc_initial:step_Nc:Nc_final;

    %============================== ACUMULADORES ===============================================
    X_all   = zeros(length(pesos), numel(Nc_values)); % Vazão total
    Y_all   = zeros(length(pesos), numel(Nc_values)); % Energia total
    R_all   = zeros(length(pesos), numel(Nc_values)); % Utilidade total
    EFF_all = zeros(length(pesos), numel(Nc_values)); % Eficiência ponderada
    conv_avg_all = zeros(length(pesos), generations); % convergência média por peso

    %============================== LOOP PRINCIPAL =============================================
    for pw = 1:size(pesos,1)
        pesoR = pesos(pw,1); pesoE = pesos(pw,2);
        conv_sum = zeros(1, generations);

        for k = 1:numel(Nc_values)
            Nc = Nc_values(k);

            % --------- Rmax e p_Rmax (maximiza utilidade agregada) ----------
            [p_Rmax, conv_R, Rmax] = genetic_algorithm_simplex( ...
                @(p) sum(utilidade_de_rede(lambda_val, p, Nc, b, Toa)), ...
                dim, pop_size, generations, crossover_rate, mutation_rate, sigma_mut, tournament_k, 'maximize');

            % --------- Energia no p_Rmax (Emax para intervalo) --------------
            Emax = sum(modelo_de_energia(V, I_id, I_st, I_tx, I_rx, Trx1, Trx2, RD1, RD2, p_Rmax, Nc, T, Toa));

            % --------- Emin e p_Emin (minimiza energia agregada) ------------
            [p_Emin, conv_E, Emin] = genetic_algorithm_simplex( ...
                @(p) sum(modelo_de_energia(V, I_id, I_st, I_tx, I_rx, Trx1, Trx2, RD1, RD2, p, Nc, T, Toa)), ...
                dim, pop_size, generations, crossover_rate, mutation_rate, sigma_mut, tournament_k, 'minimize');

            % --------- Utilidade no p_Emin (Rmin) ---------------------------
            Rmin = sum(utilidade_de_rede(lambda_val, p_Emin, Nc, b, Toa));

            % --------- Intervalos positivos ---------------------------------
            alfa = max(Rmax - Rmin, eps);
            beta = max(Emax - Emin, eps);

            % --------- Maximiza eficiência ponderada ------------------------
            [p_eff, conv_eff, ~] = genetic_algorithm_simplex( ...
                @(p) sum(eficiencia(lambda_val, p, Nc, b, Toa, V, I_id, I_st, I_tx, I_rx, Trx1, Trx2, RD1, RD2, T, alfa, beta, pesoR, pesoE)), ...
                dim, pop_size, generations, crossover_rate, mutation_rate, sigma_mut, tournament_k, 'maximize');

            % --------- Métricas no p_eff -----------------------------------
            X_all(pw, k)   = sum(vazao(lambda_val, p_eff, Nc, b, Toa));
            Y_all(pw, k)   = sum(modelo_de_energia(V, I_id, I_st, I_tx, I_rx, Trx1, Trx2, RD1, RD2, p_eff, Nc, T, Toa));
            R_all(pw, k)   = sum(utilidade_de_rede(lambda_val, p_eff, Nc, b, Toa));
            EFF_all(pw, k) = sum(eficiencia(lambda_val, p_eff, Nc, b, Toa, V, I_id, I_st, I_tx, I_rx, Trx1, Trx2, RD1, RD2, T, alfa, beta, pesoR, pesoE));

            conv_sum = conv_sum + conv_eff;
        end

        conv_avg_all(pw, :) = conv_sum ./ numel(Nc_values);
    end

    %============================== PLOTAGEM (estilo FPA) ======================================
    figure; sgtitle('Otimização com AG em Redes LoRa');

    % Número de Nós vs Vazão Total
    subplot(4,1,1); hold on; grid on;
    title('Número de Nós vs Vazão Total'); xlabel('Número de Nós'); ylabel('Vazão Total (bps)');
    for pw = 1:size(pesos,1)
        plot(Nc_values, X_all(pw, :), '-o', 'Color', colors{pw}, 'LineWidth', 1.5, 'DisplayName', labels{pw});
    end
    legend('Location','best'); hold off;

    % Número de Nós vs Energia Total
    subplot(4,1,2); hold on; grid on;
    title('Número de Nós vs Energia Total'); xlabel('Número de Nós'); ylabel('Energia Total (J)');
    for pw = 1:size(pesos,1)
        plot(Nc_values, Y_all(pw, :), '-o', 'Color', colors{pw}, 'LineWidth', 1.5, 'DisplayName', labels{pw});
    end
    legend('Location','best'); hold off;

    % Convergência média do GA
    subplot(4,1,3); hold on; grid on;
    title('Convergência média do GA'); xlabel('Iterações'); ylabel('Fitness');
    for pw = 1:size(pesos,1)
        plot(1:generations, conv_avg_all(pw, :), '-', 'Color', colors{pw}, 'LineWidth', 1.5, 'DisplayName', labels{pw});
    end
    legend('Location','best'); hold off;

    % Número de Nós vs Utilidade
    subplot(4,1,4); hold on; grid on;
    title('Número de Nós vs Utilidade'); xlabel('Número de Nós'); ylabel('Utilidade (soma)');
    for pw = 1:size(pesos,1)
        plot(Nc_values, R_all(pw, :), '-o', 'Color', colors{pw}, 'LineWidth', 1.5, 'DisplayName', labels{pw});
    end
    legend('Location','best'); hold off;

    %============================== SALVAMENTO EM .MAT =========================================
    vazao_ga   = X_all;
    energia_ga = Y_all;
    utility_ga = R_all;
    EFF_ga     = EFF_all;
    conv_ga    = conv_avg_all;
    pesos_ga   = pesos;

    save('ga_vazao.mat','Nc_values','vazao_ga','pesos_ga');
    save('ga_energia.mat','Nc_values','energia_ga','pesos_ga');
    save('ga_utility.mat','Nc_values','utility_ga','pesos_ga');
    save('ga_EFF.mat','Nc_values','EFF_ga','pesos_ga');
    save('ga_convergencia.mat','Nc_values','conv_ga','pesos_ga');
end

%================================== GENETIC ALGORITHM NO SIMPLEX ===========================================
function [best_sol, history_best, best_fit] = genetic_algorithm_simplex(fitness_func, dim, pop_size, generations, crossover_rate, mutation_rate, sigma_mut, tournament_k, mode)
    % Inicialização no simplex (p_i >= 0, sum(p)=1)
    P = rand(pop_size, dim);
    P = P ./ sum(P, 2);

    fit = arrayfun(@(i) fitness_func(P(i,:)), 1:pop_size);
    if strcmpi(mode,'maximize')
        [best_fit, idx] = max(fit); better = @(a,b) a > b; worse = @(a,b) a < b;
    else
        [best_fit, idx] = min(fit); better = @(a,b) a < b; worse = @(a,b) a > b;
    end
    best_sol = P(idx,:);
    history_best = zeros(1, generations);

    for g = 1:generations
        newP = zeros(size(P));
        newFit = zeros(size(fit));

        % elitismo
        newP(1,:) = best_sol;
        newFit(1) = best_fit;
        fill = 2;

        while fill <= pop_size
            % seleção por torneio
            p1 = tournament_select(P, fit, tournament_k, better);
            p2 = tournament_select(P, fit, tournament_k, better);

            % crossover aritmético
            if rand < crossover_rate
                w = rand;
                c1 = w .* p1 + (1-w) .* p2;
                c2 = w .* p2 + (1-w) .* p1;
            else
                c1 = p1; c2 = p2;
            end

            % mutação gaussiana + projeção no simplex
            if rand < mutation_rate
                c1 = c1 + sigma_mut .* randn(1, dim);
            end
            if rand < mutation_rate
                c2 = c2 + sigma_mut .* randn(1, dim);
            end
            c1 = project_to_simplex(c1);
            c2 = project_to_simplex(c2);

            % avaliação
            f1 = fitness_func(c1);
            f2 = fitness_func(c2);

            % inserção (respeitando objetivo)
            if fill <= pop_size
                newP(fill,:) = c1; newFit(fill) = f1; fill = fill + 1;
            end
            if fill <= pop_size
                newP(fill,:) = c2; newFit(fill) = f2; fill = fill + 1;
            end
        end

        % substituição geracional
        P = newP; fit = newFit;

        % substituição geracional
        P = newP; fit = newFit;
        
        % atualiza melhor (sem deal)
        if strcmpi(mode,'maximize')
            [curr_best_fit, idxb] = max(fit);
        else
            [curr_best_fit, idxb] = min(fit);
        end

        if better(curr_best_fit, best_fit)
            best_fit = curr_best_fit;
            best_sol = P(idxb,:);
        end

        history_best(g) = best_fit;
    end
end

function sel = tournament_select(P, fit, k, better)
    n = size(P,1);
    idx = randi(n, 1, k);
    best_i = idx(1);
    for j = 2:k
        if better(fit(idx(j)), fit(best_i)), best_i = idx(j); end
    end
    sel = P(best_i,:);
end

function x = project_to_simplex(y)
    % Projeta vetor (linha) y no simplex {x | x>=0, sum(x)=1}
    y = y(:);
    n = numel(y);
    u = sort(y, 'descend');
    cssv = cumsum(u);
    rho = find(u + (1 - cssv) ./ (1:n)' > 0, 1, 'last');
    theta = (cssv(rho) - 1) / rho;
    x = max(y - theta, 0);
    x = x' / sum(x); % retorna linha normalizada
end

%======================================== MÉTRICAS / MODELOS ===============================================
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
