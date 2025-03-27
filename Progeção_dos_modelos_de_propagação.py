import numpy as np
import matplotlib.pyplot as plt

class PathLossCalculator:
    def __init__(self, ht_m, hr_m):
        self.ht_m = ht_m  # altura da antena transmissora em metros
        self.hr_m = hr_m  # altura da antena receptora em metros

    def pathloss(self, distance, frequency_mhz, model_pathloss, fading_type=None):
        """Calcula a perda de percurso usando diferentes modelos de propagação.

        Parâmetros:
        - distance: Distância entre transmissor e receptor (km).
        - frequency_mhz: Frequência de operação (MHz).
        - model_pathloss: Modelo de propagação a ser usado:
            * "okumura_hata"
            * "log_distance"
            * "fspl" (Free-Space Path Loss)
            * "cost_hata"
            * "fading" (Path Loss + Fading Avançado)
        - fading_type: Tipo de fading a ser aplicado (se "fading" for escolhido)
            * "rayleigh"
            * "rician"
            * "nakagami"

        Retorna:
        - path_loss: Perda de percurso em dB.
        """

        if distance <= 0:
            return float('inf')  # Evita erro se a distância for zero ou negativa

        if model_pathloss == "okumura_hata":

            correction_factor = (1.1 * np.log10(frequency_mhz) - 0.7) * self.hr_m - (1.56 * np.log10(frequency_mhz) - 0.8)
            path_loss = (
                69.55 + 26.16 * np.log10(frequency_mhz)
                - 13.82 * np.log10(self.ht_m)
                - correction_factor
                + (44.9 - 6.55 * np.log10(self.ht_m)) * np.log10(distance)
            )

        elif model_pathloss == "log_distance":

            path_loss_exponent = 3.5  # Coeficiente de perda para ambientes urbanos
            path_loss = 20 * np.log10(frequency_mhz) + 10 * path_loss_exponent * np.log10(distance) - 28

        elif model_pathloss == "fspl":

            path_loss = 20 * np.log10(distance) + 20 * np.log10(frequency_mhz) + 32.45

        elif model_pathloss == "cost_hata":

            correction_factor = (1.1 * np.log10(frequency_mhz) - 0.7) * self.hr_m - (1.56 * np.log10(frequency_mhz) - 0.8)
            path_loss = (
                46.3 + 33.9 * np.log10(frequency_mhz)
                - 13.82 * np.log10(self.ht_m)
                - correction_factor
                + (44.9 - 6.55 * np.log10(self.ht_m)) * np.log10(distance)
            )

        elif model_pathloss == "fading":

            base_loss = self.pathloss(distance, frequency_mhz, model_pathloss)

            if fading_type == "rayleigh":
                fading_factor = np.random.rayleigh(scale=2)  # Desvanecimento Rayleigh
            elif fading_type == "rician":
                fading_factor = np.random.rician(3, 1)  # Desvanecimento Rician (LOS)
            elif fading_type == "nakagami":
                fading_factor = np.random.gamma(shape=3, scale=2)  # Nakagami-m
            else:
                fading_factor = 0  # Sem fading

            path_loss = base_loss + fading_factor  # Aplica o efeito do fading

        return path_loss

    def plot_pathloss_vs_distance(self, frequency_mhz, model_pathloss, max_distance=100, fading_type=None):
        """Gera gráficos das perdas de percurso em relação à distância."""
        
        distances = np.linspace(0.1, max_distance, 100)  # Distâncias de 0.1 km até max_distance km
        path_losses = [self.pathloss(d, frequency_mhz, model_pathloss, fading_type) for d in distances]

        # Plotando os gráficos
        plt.figure(figsize=(10, 6))
        plt.plot(distances, path_losses, label=model_pathloss)
        plt.title(f"Perda de Percurso em função da Distância para {model_pathloss} (Frequência: {frequency_mhz} MHz)")
        plt.xlabel("Distância (km)")
        plt.ylabel("Perda de Percurso (dB)")
        plt.grid(True)
        plt.legend()
        plt.show()

# Exemplo de uso
ht_m = 1.5  # altura da antena transmissora (metros)
hr_m = 30  # altura da antena receptora (metros)
frequency_mhz = 915  # Frequência em MHz
model_pathloss = "okumura_hata"  # Pode ser "okumura_hata", "log_distance", "fspl", "cost_hata" ou "fading"
fading_type = None  # Se "fading" for escolhido, coloque o tipo de fading como "rayleigh", "rician" ou "nakagami"

# Criando o objeto e gerando o gráfico
pathloss_calculator = PathLossCalculator(ht_m, hr_m)
pathloss_calculator.plot_pathloss_vs_distance(frequency_mhz, model_pathloss, max_distance=100, fading_type=fading_type)
