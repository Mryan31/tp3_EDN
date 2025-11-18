import numpy as np

def f_x(x):
    return np.piecewise(x, 
                        [np.logical_and(x > 0, x <= 0.5), 
                         np.logical_and(x > 0.5, x < 1)], 
                        [lambda x: 2*x, 
                         lambda x: 2 - 2*x]
                       )

def get_bn(n):
    """
    Calcula o coeficiente b_n da série de Fourier para a f(x) dada.
    
    Após integração, o coeficiente é: b_n = (8 * sin(n*pi/2)) / (n*pi)^2
    """
    if n % 2 == 0:
        return 0.0
    
    n_pi = n * np.pi
    return (8 * np.sin(n_pi / 2)) / (n_pi**2)

def find_N_termos(epsilon=0.01):
    N = 0
    erro_max_teorico = 0.0
    
    # A série de b_n converge (Soma(1/n^2)). 
    # Podemos somar os termos em ordem decrescente (mais estável)
    # ou apenas somar os primeiros N e estimar o resto.
    
    # Por simplicidade, vamos calcular N.
    # A soma da cauda (n>N) é <= Sum(8/(n*pi)^2) ~ 8/(pi^2) * (1/N)
    # 8/(pi^2 * N) < 0.01 => N > 8 / (0.01 * pi^2) ~ 800 / 9.87 ~ 81
    # Vamos usar N=85 por segurança.
    N_calculado = int(8 / (epsilon * np.pi**2)) + 5
    print(f"Número de termos N calculado para erro {epsilon}: {N_calculado}")
    return N_calculado

def solucao_fourier(x_vec, t, N):
    if isinstance(t, np.ndarray):
        x_grid, t_grid = np.meshgrid(x_vec, t)
    else:
        x_grid = x_vec
        t_grid = t

    U = np.zeros_like(x_grid, dtype=float)
    
    for n in range(1, N + 1):
        bn = get_bn(n)
        if bn != 0.0:
            termo = bn * np.sin(n * np.pi * x_grid) * np.exp(-(n * np.pi)**2 * t_grid)
            U += termo
            
    return U