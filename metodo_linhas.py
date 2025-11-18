import numpy as np
from scipy.integrate import solve_ivp
from .fourier import f_x

def sistema_MOL(t, U, M, dx_quadrado_inv):
    dUdt = np.zeros_like(U)
    
    for j in range(1, M - 2): # de U_2 a U_{M-2}
        dUdt[j] = (U[j+1] - 2*U[j] + U[j-1]) * dx_quadrado_inv
        
    # Casos especiais perto das fronteiras (usando U_0=0 e U_M=0)
    # Para j=0 (ponto U_1):
    # dU_1/dt = (U_2 - 2*U_1 + U_0) / dx^2 = (U_2 - 2*U_1 + 0) / dx^2
    dUdt[0] = (U[1] - 2*U[0] + 0.0) * dx_quadrado_inv
    
    # Para j=M-2 (ponto U_{M-1}):
    # dU_{M-1}/dt = (U_M - 2*U_{M-1} + U_{M-2}) / dx^2 = (0 - 2*U_{M-1} + U_{M-2}) / dx^2
    dUdt[-1] = (0.0 - 2*U[-1] + U[-2]) * dx_quadrado_inv
    
    return dUdt

def solucao_MOL(delta_x, t_final, t_pontos):
    """
    Resolve a equação do calor usando o Método das Linhas.
    """
    L = 1.0
    M = int(L / delta_x) + 1  # Número de pontos na malha (ex: dx=0.1 -> M=11)
    
    # Vetor de pontos x (incluindo fronteiras 0 e 1)
    x_vec = np.linspace(0, L, M)
    
    # Pegamos apenas os pontos internos para o solver de EDO
    x_internos = x_vec[1:-1]
    
    # Condição inicial u(x,0) = f(x) [cite: 7]
    U0 = f_x(x_internos)
    
    dx_quadrado_inv = 1.0 / (delta_x**2)
    M_internos = len(x_internos) # M-2 pontos
    
    print(f"Executando MOL com Delta_x = {delta_x} (M={M} pontos)")
    
    # Resolve o sistema de EDOs
    sol = solve_ivp(
        fun=sistema_MOL,           
        t_span=[0, t_final],       
        y0=U0,                     
        t_eval=t_pontos,           
        args=(M_internos, dx_quadrado_inv),
        method='RK45'              
    )
    
    # Remonta a solução completa (adicionando as fronteiras 0)
    U_completo = np.zeros((M, len(t_pontos)))
    U_completo[1:-1, :] = sol.y
    
    return x_vec, sol.t, U_completo