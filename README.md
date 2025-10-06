# README — Otimização Conjunta de Vazão e Energia em Redes LoRaWAN (CVX, GA e FPA)

> Projeto Matlab para **planejamento e otimização multiobjetivo** em redes **LoRaWAN multi-gateway**, comparando **CVX** (programação convexa), **Algoritmo Genético (AG)** e **Flower Pollination Algorithm (FPA)**. O objetivo é determinar a **proporção ótima de nós por SF** (7–12) que **maximiza a utilidade/vazão** e **minimiza o consumo de energia**, com normalização por valores de **Utopia/Nadir** (α, β) e varredura de **pesos (a,b)**. 

---

## 1) Escopo e visão geral

* **Domínio**: LPWAN/LoRaWAN multi-gateway, com foco em **vazão**, **energia** e **eficiência energética**. 
* **Problema**: decidir a distribuição (p(s)) de nós por **spreading factor** (s\in{7,\dots,12}) sob **ALOHA** e **janelas RX** para **equilibrar** desempenho e consumo. 
* **Abordagens**:

  * **CVX** (baseline clássico, solução ótima para formulações convexas);
  * **AG** (metaheurística populacional, contínuo no **simplex**);
  * **FPA** (metaheurística com **polinização global/local** via **Levy flight**). 

---

## 2) Estrutura (arquivos principais)

* `cvx_otimizacao.m`
  Implementa o pipeline em **CVX**: calcula (R_{\max}), (E_{\min}), obtém (\alpha,\beta) e **maximiza a eficiência** ( \text{EFF} = \frac{a}{\alpha}R(p) - \frac{b}{\beta}E(p) ). Gera gráficos e arquivos `.mat` de **vazão**, **energia**, **utilidade** e **EFF**.
* `ag_optimization_lora.m`
  **Algoritmo Genético** no **simplex** (proporções (p)), com **torneio**, **crossover aritmético**, **mutação gaussiana** e **projeção no simplex**. Repete o mesmo pipeline (Rmax/Emin/α/β/EFF) e salva resultados em `.mat`.
* `fpa_optimization_lora.m`
  **FPA** com alternância **global/local** (probabilidade `switch_prob`), **Levy flight**, **normalização/inviabilidade** tratadas com **projeção e renormalização** para o simplex. Mesmo pipeline e `.mat`.
* `TCC_Ivan_Neves2025.pdf`
  Artigo que fundamenta a modelagem, cenário, métricas e referência bibliográfica completa. 

---

## 3) Pré-requisitos

* **MATLAB** R2020b+ (testado); toolboxes padrão.
* **CVX** (se usar o script CVX):

  * Instale **CVX** e um solver compatível (ex.: SDPT3/SeDuMi) conforme a doc oficial do CVX.
* Sistema operacional: Windows/Linux/macOS.

---

## 4) Modelos e formulações (detalhado)

### 4.1 Vazão por SF (ALOHA)

Para SF (s):
[
T(s) ;=; \lambda , p(s), N_c, b ;\exp!\big(-2,G(s)\big),\quad
G(s) ;=; \lambda , p(s), N_c, \text{ToA}(s).
]
**Vazão total**:
[
S = \sum_{s=S_{\min}}^{S_{\max}} T(s).
]
A janela de vulnerabilidade (2\cdot\text{ToA}) penaliza colisões em ALOHA; **ToA** cresce com **SF** (para **BW** fixa), degradando (\exp(-2G)). 

### 4.2 Energia (transmissão/recepção/espera/idle)

Energia **ativa** (por SF):
[
E_{\text{active}}(s) = \tfrac{1}{2} p(s) N_c V\left[
\text{ToA}(s),I_{tx} + RD_1 I_{rx} + \text{Trx1}(s) I_{rx} + (RD_2-\text{Trx1}(s)) I_{st} + (\text{Trx1}(s)+\text{Trx2}(s)) I_{rx}
\right].
]
Energia **idle** (por SF):
[
E_{\text{idle}}(s) = \tfrac{1}{2} p(s) N_c V\left[
T - \text{ToA}(s) + \text{Trx1}(s) + RD_1 + RD_2
\right] I_{id}.
]
Total por SF: ( E(s) = E_{\text{active}}(s) + E_{\text{idle}}(s) ).
Total na rede: ( E = \sum_s E(s)). 

### 4.3 Utilidade proporcionalmente justa

[
R(p) = \sum_{s} \log!\big(T(s)\big)
= \sum_s \log!\big(\lambda p(s) N_c b\big) ;-; 2\sum_s G(s).
]
Promove **justiça proporcional**, evita concentração excessiva em um SF e **mitiga colisões**. 

### 4.4 Multiobjetivo normalizado

[
\max_{p};; \underbrace{\frac{a}{\alpha}R(p)}*{\text{desempenho}} ;-; \underbrace{\frac{b}{\beta}E(p)}*{\text{energia}}
\quad\text{s.a.}\quad \sum_s p(s)=1,; 0\le p(s)\le 1.
]
Com:
[
\alpha = R_{\max}-R_{\min}, \qquad \beta = E_{\max}-E_{\min}, \qquad a,b\in[0,1],; a+b=1.
]
O **pipeline** resolve (R_{\max}), (E_{\min}), computa (\alpha,\beta) e então **maximiza** a eficiência ponderada (**EFF**). 

---

## 5) Parâmetros do cenário (default nos scripts)

* **SF**: 7–12; **Nc**: 500:1000:4500; **λ**: 6; **b**: 48 bits/pacote;
* **ToA**: `[0.1048 0.1802 0.3211 0.5636 1.0485 1.9398]` (s);
* **Trx1/Trx2**: conforme vetores no código; **RD1=1 s**, **RD2=2 s**;
* **V**: 3.3 V; **I_tx=44 mA**, **I_rx=10.5 mA**, **I_st=1.4 mA**, **I_id=0.0015 mA**;
* **T**: 720 s;
* **Pesos (a,b)**: `[(1,0),(0.75,0.25),(0.5,0.5),(0.25,0.75),(0.1,0.9)]`.
  Esses parâmetros refletem o **caso de estudo** com 4 gateways, BW=125 kHz, CR=4/5, operação a 868 MHz etc., conforme o artigo. 

---

## 6) Como executar

> **Importante**: CVX é opcional; GA e FPA não exigem CVX.

1. **Abrir MATLAB** no diretório do projeto.
2. (Opcional – CVX) **Instalar/ativar** o CVX (uma vez no MATLAB):

   ```matlab
   run('caminho/para/cvx/cvx_setup.m');  % ajuste o caminho
   ```
3. **Executar um método** por vez:

   * **CVX**:

     ```matlab
     cvx_otimizacao   % roda o pipeline CVX
     ```
   * **GA**:

     ```matlab
     ag_optimization_lora
     ```
   * **FPA**:

     ```matlab
     fpa_optimization_lora
     ```
4. **Saídas**: figuras + arquivos `.mat` na pasta do projeto.

---

## 7) Saídas (artefatos gerados)

Cada método salva, para todos os **Nc** e **pesos (a,b)**:

* **Vazão** (`*_vazao.mat`):

  * `Nc_values`: vetor de nós;
  * `vazao_*`: matriz `[numPesos × numNc]` com **S** total;
  * `pesos_*`: matriz com as combinações (a,b).
* **Energia** (`*_energia.mat`):

  * `energia_*`: matriz de **E** total.
* **Utilidade** (`*_utility.mat`):

  * `utility_*`: matriz de **R(p)** total.
* **Eficiência (EFF)** (`*_EFF.mat`):

  * `EFF_*`: matriz do **fitness** composto (\tfrac{a}{\alpha}R - \tfrac{b}{\beta}E).
* **Convergência**:

  * GA: `ga_convergencia.mat` (fitness médio por geração);
  * FPA: `fpa_convergencia.mat` (fitness médio por iteração).

Gráficos:

* **CVX**: 4 subplots (Vazão×Nós, Energia×Nós, **EFF**×Nós, **Utilidade**×Nós) para todos os pesos.
* **GA/FPA**: 4 subplots (Vazão×Nós, Energia×Nós, **Convergência**, **Utilidade**×Nós).

---

## 8) Personalização de parâmetros

* **Varredura de rede**: edite `Nc_initial`, `Nc_final`, `step_Nc`.
* **Pesos multiobjetivo**: edite a matriz `pesos` (linhas somam a 1).
* **Métricas físicas**: ajuste vetores `Toa`, `Trx1`, `Trx2`, e correntes/tensão conforme rádio/board usados. 
* **AG**: `pop_size`, `generations`, `crossover_rate`, `mutation_rate`, `sigma_mut`, `tournament_k`.
* **FPA**: `n` (população), `max_iter`, `switch_prob`.
* **Reprodutibilidade**: `rng(1)` já definido no AG.

---

## 9) Detalhes de implementação (pontos críticos)

* **Simplex**: (p\ge 0), (\sum p=1).

  * **AG**: mutação gaussiana + **projeção ao simplex** (método de Michelot/condensado por ordenação/threshold).
  * **FPA**: perturbação global (**Levy**) e local; **clamp** em ([0,1]) e renormalização/ajuste ao simplex.
* **Normalização (α,β)**: calculadas **por cenário (Nc)** para cada par de pesos, usando **Rmax** (max (R)) e **Emin** (min (E)); o código também obtém **Emax**/**Rmin** para estabilidade.
* **Estabilidade numérica**:

  * `log(...)` com **`+eps`** no AG/FPA;
  * proteção de divisão por zero em (\alpha,\beta \leftarrow \max(\cdot,\varepsilon)).
* **ALOHA**: (e^{-2G(s)}) com (G(s)=\lambda p(s)N_c \text{ToA}(s)).
* **Energia**: inclui **TX**, **RX em janelas**, **standby** e **idle** no horizonte **T** (vide equações). 

---

## 10) Interpretação dos resultados (o que observar)

* **Trade-off** esperado: ao **subir (a)** (peso em desempenho), a vazão cresce e a energia também; ao **subir (b)**, a energia cai com perda de vazão — **EFF** reflete o melhor compromisso. 
* **CVX**: tende a entregar **vazões mais altas** em cenários favoráveis, mas com **maior custo energético**; é **mais lento** computacionalmente. 
* **AG**: **bom equilíbrio** entre qualidade e tempo, com **eficiência energética alta** em vários cenários. 
* **FPA**: **tempo menor** e **energia baixa**, com vazão geralmente inferior ao CVX; útil quando **latência de otimização** importa. 

---

## 11) Validação e extensões

* **Valide** os gráficos gerados confrontando **Vazão×Nós** e **Energia×Nós** entre métodos.
* **Convergência**: curvas (GA/FPA) devem estabilizar < 200 iterações/gerações com ordenamento coerente por pesos. 
* **Extensões**:

  * Inclusão de **BW/CR/TP** na decisão (hoje fixos nos scripts);
  * **ADR reativo** por nó (nível de pacote), ou **DRL** para adaptação online;
  * **Posicionamento de GW** e **alocação de canal/frequência**. 

---

## 12) FAQ rápido

* **Preciso do CVX?** Somente para o script CVX. AG e FPA rodam **sem** CVX.
* **As somas de (p) explodem?** Não. Todos os métodos **projetam/normalizam** no simplex.
* **Unidades**: **ToA** em segundos; correntes em **mA**; **V** em volts; **T** em segundos; **b** em **bits/pacote**.
* **Ortogonalidade de SF**: assumida para compor vazão **aditiva** entre SFs (hipótese clássica em LoRa). 

---

## 13) Citações rápidas ao artigo

* **Modelos de Vazão/Energia/Utilidade e Multiobjetivo**: fundamentos e equações no PDF. 
* **Cenário, parâmetros físicos e correntes (SX1272)**: conforme tabelas/texto do PDF. 
* **Discussão CVX vs GA vs FPA (desempenho/tempo)**: seção de resultados do PDF. 

---

## 14) Referências (as do artigo)

> A lista a seguir replica as referências usadas no PDF associado. Para detalhes completos, veja o arquivo `TCC_Ivan_Neves2025.pdf`. 

* Abdelhedi, M. A., Trabelsi, H., Derbel, F. (2025). **Performance evaluation of LoRaWAN physical transmission parameters**. *IEEE SSD 2025*, pp. 1180–1185.
* Al-zamili, J. J., Al-Zubaidi, H. A. (2024). **Optimizing IoT WSNs: PSO vs GA**. *Fusion: Practice and Applications*, 15(02):278–287.
* Banti, K., et al. (2022). **LoRaWAN communication protocols: energy efficiency survey**. *Telecom*, 3(2):322–357.
* Casals, L., et al. (2017). **Modeling the energy performance of LoRaWAN**. *Sensors*, 17(10).
* Faisal, M., et al. (2020). **Review of Flower Pollination Algorithm**. *IEEE*.
* Goldberg, D. E. (1989). **Genetic Algorithms in Search, Optimization, and Machine Learning**. Addison-Wesley.
* Grant, M., Boyd, S. (2014). **CVX: Matlab Software for Disciplined Convex Programming**.
* Holland, J. H. (1975). **Adaptation in Natural and Artificial Systems**. Univ. of Michigan Press.
* Javaid, M., et al. (2022). **Agriculture 4.0 technologies**. *International Journal of Intelligent Networks*, 3:150–164.
* Jiang, S., et al. (2022). **Evolutionary dynamic multi-objective optimisation: survey**. *ACM Computing Surveys*, 55(4):76:1–76:47.
* Loubany, A., Lahoud, S., El Chall, R. (2020). **Adaptive SF selection in LoRaWAN with multiple gateways**. *Computer Networks*.
* Loubany, A., et al. (2023). **Joint throughput-energy optimization in multi-gateway LoRaWAN**. *Telecommunication Systems*, 84(2):271–283.
* Maurya, P., Sørensen, T. B., Sharma, H. (2024). **LoRaWAN Relay mode impact on end node energy**. *IEEE MetroAgriFor 2024*, 683–688.
* Mekki, K., et al. (2019). **Comparative study of LPWAN technologies**. *ICT Express*, 5:1–7.
* Mitchell, M. (1998). **An Introduction to Genetic Algorithms**. MIT Press.
* Neves, I. I. A., et al. (2024). **Determinação ótima de SF em redes IoT-LoRa**. *ENCOM 2024, Brasil*.
* Sahu, R., Tripathi, P. (2024). **GA-based optimization for large-scale LoRaWAN**. *ICACRS*, 553–561.
* Semtech Corporation (2013). **SX1272/3/6/7/8 LoRa Modem Design Guide**, AN1200.13 Rev 1.
* Semtech Corporation (2019). **SX1272/3 datasheet, Rev. 4**.
* Wang, H., et al. (2024). **Improved ADR for mobile agricultural nodes**. *Computers and Electronics in Agriculture*, 219:108773.
* Zheng, J., et al. (2024). **Dynamic parameter tuning for MOEAs**. *Applied Sciences*, 14(8):3481.

---

## 15) Licença e crédito

* Créditos do artigo: **Ivan Neves, Caio Cardoso, Fabricio Barros, Jasmine Araujo — UFPA** (Maio/2025). 
* Verifique a licença dos scripts/artefatos antes de redistribuir.

---

## 16) Contato

Dúvidas técnicas sobre execução, parâmetros ou extensões: abra uma *issue* no repositório associado citado no PDF ou ajuste diretamente os scripts conforme a sua topologia e requisitos. 


Medical References:
1. None — DOI: file_00000000c80c622f9e3670a5706f3dff
