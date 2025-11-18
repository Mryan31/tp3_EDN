import numpy as np
from scipy.linalg import solve_banded

def setup_grid(M, N, mu):
    L = 1.0
    
    # Malha espacial
    dx = L / (M - 1)
    x_vec = np.linspace(0, L, M)
    
    # Malha temporal
    # r = dt / dx^2 => dt = r * dx^2
    r = mu # Assumindo que mu é o 'r' padrão
    dt = r * (dx**2)
    
    # Recalcula N baseado no T_final (a ser definido no main)
    # Por enquanto, apenas criamos N passos
    t_vec = np.linspace(0, N * dt, N + 1)
    
    U = np.zeros((M, N + 1))
    return x_vec, t_vec, U, M, N, r, dt, dx

def solver_ftcs(U, M, N, r):
    """
    U_j^{n+1} = U_j^n + r * (U_{j+1}^n - 2*U_j^n + U_{j-1}^n)
    """
    for n in range(N):
        for j in range(1, M - 1):
            U[j, n+1] = U[j, n] + r * (U[j+1, n] - 2*U[j, n] + U[j-1, n])
    return U

def solver_btcs(U, M, N, r):
    """
    Resolve usando BTCS (Regressivo, Implícito).
    -r*U_{j-1}^{n+1} + (1+2r)*U_j^{n+1} - r*U_{j+1}^{n+1} = U_j^n
    
    """
    # Diagonais para solve_banded: (superior, principal, inferior)
    # O sistema é (M-2)x(M-2).
    # A diag. principal tem M-2 elementos.
    # As sub/super-diagonais têm M-3 elementos.
    diag_sup = np.full(M - 3, -r)
    diag_main = np.full(M - 2, 1 + 2*r)
    diag_inf = np.full(M - 3, -r)
    
    # Formato 'ab' para solve_banded:
    # (n_sup, n_inf) -> (1, 1)
    # Linha 0: 0, -r, -r, ... (superior)
    # Linha 1: 1+2r, 1+2r, ... (principal)
    # Linha 2: -r, -r, ..., 0 (inferior)
    A_banded = np.zeros((3, M - 2))
    A_banded[0, 1:] = diag_sup
    A_banded[1, :] = diag_main
    A_banded[2, :-1] = diag_inf
    
    for n in range(N):
        # O vetor 'b' é a solução no tempo anterior (U_j^n)
        b = U[1:-1, n]
        
        # Resolve o sistema A * U^{n+1} = b
        U[1:-1, n+1] = solve_banded((1, 1), A_banded, b)
        
    return U

def solver_crank_nicholson(U, M, N, r):
    """
    Resolve usando Crank-Nicholson.
    
    Lado Esquerdo (tempo n+1):
    -r/2 * U_{j-1} + (1+r) * U_j - r/2 * U_{j+1}
    
    Lado Direito (tempo n):
    r/2 * U_{j-1} + (1-r) * U_j + r/2 * U_{j+1}
    """
    # --- Monta Matriz A (Lado Esquerdo) ---
    r2 = r / 2.0
    diag_sup_A = np.full(M - 3, -r2)
    diag_main_A = np.full(M - 2, 1 + r)
    diag_inf_A = np.full(M - 3, -r2)
    
    A_banded = np.zeros((3, M - 2))
    A_banded[0, 1:] = diag_sup_A
    A_banded[1, :] = diag_main_A
    A_banded[2, :-1] = diag_inf_A

    for n in range(N): # Itera no tempo
        # --- Monta Vetor B (Lado Direito) ---
        U_n = U[1:-1, n] # Pontos internos no tempo n
        b = np.zeros(M - 2)
        
        for j in range(1, M - 3): # Pontos internos de b
            b[j] = r2 * U_n[j-1] + (1 - r) * U_n[j] + r2 * U_n[j+1]
        
        # Fronteiras de b (j=0 e j=M-3)
        # j=0 (ponto U_1)
        b[0] = r2 * U[0, n] + (1 - r) * U_n[0] + r2 * U_n[1] # U[0,n] = 0
        # j=M-3 (ponto U_{M-2})
        b[-1] = r2 * U_n[-2] + (1 - r) * U_n[-1] + r2 * U[M-1, n] # U[M-1,n] = 0

        # Resolve o sistema A * U^{n+1} = b
        U[1:-1, n+1] = solve_banded((1, 1), A_banded, b)

    return U