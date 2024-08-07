import numpy as np


def get_routh_hurwitz_table(coeffs: np.ndarray) -> np.ndarray:
    """
    Calcula a tabela de Routh-Hurwitz para os coeficientes de um polinômio.

    Args:
    coeffs: np.ndarray, Coeficientes do polinômio.

    Returns:
    table: np.ndarray, A tabela de Routh-Hurwitz.
    """
    if not isinstance(coeffs, np.ndarray):
        coeffs = np.array(coeffs)

    n = len(coeffs)
    table = np.zeros((n, (n + 1) // 2))

    # Preenche as duas primeiras linhas da tabela
    even, odd = coeffs[::2], coeffs[1::2]
    table[0, : len(even)] = even
    table[1, : len(odd)] = odd

    # Preenche o restante da tabela
    for i in range(2, len(coeffs)):
        new_row = (
            table[i - 1, 0] * table[i - 2, 1:] - table[i - 2, 0] * table[i - 1, 1:]
        ) / table[i - 1, 0]

        table[i, : len(new_row)] = new_row

    return table


def count_unstable_poles(table: np.ndarray) -> int:
    """
    Conta o número de polos instáveis de um sistema a partir dos coeficientes de um polinômio.

    Args:
    table: np.ndarray, A tabela de Routh-Hurwitz.

    Returns:
    unstable_poles: int, Número de polos instáveis. -1 em caso de caso especial.
    """
    # Caso especial
    if np.any(table[:, 0] == 0):
        return -1

    unstable_poles = 0
    # Conta o número de mudanças de sinal na primeira coluna
    for i in range(1, len(coeffs)):
        if table[i, 0] * table[i - 1, 0] < 0:
            unstable_poles += 1

    return unstable_poles


def main(coeffs: np.ndarray, print_table: bool = False) -> None:
    """
    Função principal que calcula o número de polos instáveis de um sistema a partir dos coeficientes de um polinômio.
    """
    table = get_routh_hurwitz_table(coeffs)
    unstable_poles = count_unstable_poles(table)

    if print_table:
        print("Tabela de Routh-Hurwitz:")
        print(table)

    if unstable_poles == -1:
        print("Caso especial.")
    else:
        print("O polinômio possuí", unstable_poles, "polos instáveis.")


if __name__ == "__main__":
    # Exemplo
    coeficientes = [1, 3, 2, 7]

    main(coeficientes)
