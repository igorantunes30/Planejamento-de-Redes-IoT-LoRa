%=============================================PARAMETROS====================================================
S_min = 7;                                                % valor mínimo do SF
S_max = 12;                                               % valor máximo do SF
Nc_initial = 500;                                         % número inicial de nós cobertos
Nc_final = 4500;                                          % número final de nós cobertos
lambda = 6;                                               % taxa de geração de pacotes por unidade de tempo
b = 48;                                                   % a relação entre SF e número de bits por pacote
Toa = [0.1048, 0.1802, 0.3211, 0.5636, 1.0485, 1.9398]';  % Time on air (seconds)
Trx1 = [1.1048, 1.1802, 1.3211, 1.5636, 2.0485, 2.9398]'; % janela de recepção 1 (em segundos)
Trx2 = [2.1048, 2.1802, 2.3211, 2.5636, 3.0485, 3.9398]'; % janela de recepção 2 (em segundos)
V = 3.3;                                                  % tensão (em Volts)
I_tx = 44;                                                % mA corrente no modo de transmissão
I_rx = 10.5;                                              % mA corrente no modo de recepção
I_st = 1.4;                                               % mA corrente no modo de espera
I_id = 0.0015;                                            % mA corrente no modo inativo
RD1 = 1;                                                  % segundos
RD2 = 2;                                                  % segundos 
T = 720;                                                  % segundos
n_sf = S_max - S_min + 1;

% Pesos para as combinações lineares
pesos = [1, 0; 0.75, 0.25; 0.5, 0.5; 0.25, 0.75; 0.1, 0.9];
labels = {'a=1, b=0', 'a=0.75, b=0.25', 'a=0.5, b=0.5', 'a=0.25, b=0.75', 'a=0.1, b=0.9'};
colors = {'r', 'b', 'g', 'm', 'k'};  % Diferenciar por cor

Nc_values = Nc_initial:1000:Nc_final; % Corrigido para gerar os valores de Nc a cada 100

% Inicializa os dados para cada combinação de peso
X_all = zeros(length(pesos), length(Nc_values));
Y_all = zeros(length(pesos), length(Nc_values));
EFF_all = zeros(length(pesos), length(Nc_values));      % <-- adicionado
R_all = zeros(length(pesos), length(Nc_values));        % <-- adicionado

% Executa os cálculos para todas as combinações
for p_idx = 1:length(pesos)
    % Define os pesos atuais para R(a) e E(b)
    pesoR = pesos(p_idx, 1);
    pesoE = pesos(p_idx, 2);
    X = zeros(size(Nc_values)); % Vazão (Throughput)
    Y = zeros(size(Nc_values)); % Energia
    EFF_curve = zeros(size(Nc_values)); % <-- adicionado
    R_curve = zeros(size(Nc_values));   % <-- adicionado

    for idx = 1:length(Nc_values)
        Nc = Nc_values(idx);

        %===============Rmax e Emin====================================================================    
        cvx_begin 
            variable p(n_sf)   
            expression Rmax(n_sf)
            Rmax = UTILIDADE_DE_REDE(lambda, p, Nc, b, Toa);
            Emin = MODELO_DE_ENERGIA(V, I_id, I_st,I_tx,I_rx,Trx1,Trx2,RD1,RD2, p, Nc,T, Toa);
            maximize(sum(Rmax))        
            subject to
                sum(p) == 1; % A soma das proporções = 1
                p >= 0;
                p <= 1;
        cvx_end   

        % =============Emax e Rmin====================================================================   
        cvx_begin 
            variable p(n_sf)
            expression Emax(n_sf)
            Emax = MODELO_DE_ENERGIA(V, I_id, I_st,I_tx,I_rx,Trx1,Trx2,RD1,RD2, p, Nc,T, Toa);
            Rmin = UTILIDADE_DE_REDE(lambda, p, Nc, b, Toa);
            minimize(sum(Emax))
            subject to
                sum(p) == 1;
                p >= 0;
                p <= 1;
        cvx_end   

        %=============Cálculo de alfa e beta==========================================================
        alfa = Rmax - Rmin;
        beta = Emax - Emin; 

        %============Maximizar a eficiência EFF====================================================== 
        cvx_begin 
            variable p(n_sf)
            expression EFF(n_sf)
            EFF = (sum(pesoR/alfa) * UTILIDADE_DE_REDE(lambda, p, Nc, b, Toa)) - (sum(pesoE/beta) * MODELO_DE_ENERGIA(V,I_id, I_st,I_tx,I_rx,Trx1,Trx2,RD1,RD2, p, Nc , T, Toa));
            maximize(sum(EFF)) 
            subject to
                sum(p) == 1;
                p >= 0;
                p <= 1;
        cvx_end    

        %===========Vazão, Energia, EFF e Utilidade (para plot)======================================    
        X(idx) = sum(VAZAO(lambda, p, Nc, b, Toa));
        Y(idx) = sum(MODELO_DE_ENERGIA(V, I_id, I_st,I_tx,I_rx,Trx1,Trx2,RD1,RD2, p, Nc,T, Toa));
        % valores agregados para os novos subplots
        R_curve(idx) = sum(UTILIDADE_DE_REDE(lambda, p, Nc, b, Toa));                                   % <-- adicionado
        EFF_curve(idx) = sum( (pesoR./alfa) .* UTILIDADE_DE_REDE(lambda, p, Nc, b, Toa) ...             % <-- adicionado
                              - (pesoE./beta) .* MODELO_DE_ENERGIA(V, I_id, I_st,I_tx,I_rx,Trx1,Trx2,RD1,RD2, p, Nc , T, Toa) );
    end
    
    % Armazena os valores para todas as combinações
    X_all(p_idx, :) = X;
    Y_all(p_idx, :) = Y;
    EFF_all(p_idx, :) = EFF_curve;  % <-- adicionado
    R_all(p_idx, :) = R_curve;      % <-- adicionado
end

% Plotagem dos gráficos para todas as combinações
figure;

% Vazão vs Número de Nós
subplot(2, 2, 1); % <-- alterado de (2,1,1) para (2,2,1)
hold on;
for p_idx = 1:length(pesos)
    plot(Nc_values, X_all(p_idx, :), 'Color', colors{p_idx}, 'LineWidth', 1.5); % Espessura reduzida
end
xlabel('Número de Nós');
ylabel('Vazão (Throughput)');
title('Relação Vazão vs. Número de Nós');
legend(labels, 'Location', 'best');
hold off;

% Energia vs Número de Nós
subplot(2, 2, 2); % <-- alterado de (2,1,2) para (2,2,2)
hold on;
for p_idx = 1:length(pesos)
    plot(Nc_values, Y_all(p_idx, :), 'Color', colors{p_idx}, 'LineWidth', 1.5); % Espessura reduzida
end
xlabel('Número de Nós');
ylabel('Energia Consumida (J)');
title('Relação Energia vs. Número de Nós');
legend(labels, 'Location', 'best');
hold off;

% EFF vs Número de Nós  (NOVO)
subplot(2, 2, 3);
hold on;
for p_idx = 1:length(pesos)
    plot(Nc_values, EFF_all(p_idx, :), 'Color', colors{p_idx}, 'LineWidth', 1.5);
end
xlabel('Número de Nós');
ylabel('EFF (adimensional)');
title('Relação EFF vs. Número de Nós');
legend(labels, 'Location', 'best');
hold off;

% UTILIDADE_DE_REDE vs Número de Nós  (NOVO)
subplot(2, 2, 4);
hold on;
for p_idx = 1:length(pesos)
    plot(Nc_values, R_all(p_idx, :), 'Color', colors{p_idx}, 'LineWidth', 1.5);
end
xlabel('Número de Nós');
ylabel('UTILIDADE\_DE\_REDE (soma)');
title('Relação UTILIDADE\_DE\_REDE vs. Número de Nós');
legend(labels, 'Location', 'best');
hold off;

%============================ SALVAMENTO EM .MAT ============================
vazao_cvx   = X_all;
energia_cvx = Y_all;
utility_cvx = R_all;
EFF_cvx     = EFF_all;
pesos_cvx   = pesos;

save('cvx_vazao.mat','Nc_values','vazao_cvx','pesos_cvx');
save('cvx_energia.mat','Nc_values','energia_cvx','pesos_cvx');
save('cvx_utility.mat','Nc_values','utility_cvx','pesos_cvx');
save('cvx_EFF.mat','Nc_values','EFF_cvx','pesos_cvx');

%=============================== FUNÇÕES ====================================
function G = TRAFEGO_DE_CARGA(lambda, p, Nc, Toa)
    G = lambda .* p .* Nc .* Toa ./10000;
end

function X = VAZAO(lambda, p, Nc, b, Toa) 
    X = lambda .* p .* Nc .* b .* exp(-2 .* TRAFEGO_DE_CARGA(lambda, p, Nc, Toa));
end

function E = MODELO_DE_ENERGIA(V, I_id, I_st,I_tx,I_rx,Trx1,Trx2,RD1,RD2, p, Nc,T, Toa)
  E = (0.5 .* p .* Nc .* V .* (Toa .* I_tx + RD1 .* I_st + Trx1 .* I_rx) +...
            0.5 .* p .* Nc .* V .* (Toa .* I_tx + (RD2 - Trx1) .* I_st +...
            (Trx1 + Trx2) .* I_rx)) + (0.5 .* p .* Nc .* V .* (T - (Toa + Trx1 + ...
            RD1)).*I_id +  0.5 .* p .* Nc .* V .* (T - (Toa + RD2 + Trx2)).*I_id);
end

function R = UTILIDADE_DE_REDE(lambda, p, Nc, b, Toa)    
        R = (log(lambda .* p .* Nc .* b)) - 2 * (TRAFEGO_DE_CARGA(lambda, p, Nc, Toa));  
end
