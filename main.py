import numpy as np

from . import fourier
from . import metodo_linhas
from . import dif_finitas_calor as fdm
from . import graficos

def calcular_erro_l2(U_num, U_exata):
    erro = np.sqrt(np.sum((U_num - U_exata)**2))
    norma_exata = np.sqrt(np.sum(U_exata**2))
    
    if norma_exata == 0:
        return 0.0
    
    return erro / norma_exata

def run_tarefa_1_fourier(x_malha, t_malha):
    print("--- Executando Tarefa 1: Solução de Fourier ---")
    N = fourier.find_N_termos(epsilon=0.01)
    
    U_exata = fourier.solucao_fourier(x_malha, t_malha, N)
    
    x_completo = np.linspace(0, 1, len(x_malha))
    
    graficos.plotar_solucao_3d(x_completo, t_malha, U_exata.T, 
                                   "Tarefa 1 - Solucao Fourier")
    
    t_indices_plot = [0, int(len(t_malha)*0.1), int(len(t_malha)*0.5), len(t_malha)-1]
    graficos.plotar_snapshots(x_completo, U_exata.T, t_malha, t_indices_plot,
                                  "Tarefa 1 - Fourier Snapshots")
    
    print("--- Tarefa 1 Concluída ---")
    return U_exata, x_completo

def run_tarefa_2_mol(U_exata_grid, x_exato, t_malha):
    print("\n--- Executando Tarefa 2: Método das Linhas ---")
    
    dx_lista = [0.1, 0.05, 0.025, 0.0125]
    erros_mol = []
    
    T_FINAL = t_malha[-1]
    
    for dx in dx_lista:
        x_mol, t_mol, U_mol = metodo_linhas.solucao_MOL(
            delta_x=dx, 
            t_final=T_FINAL, 
            t_pontos=t_malha
        )
        
        U_exata_interp = fourier.solucao_fourier(x_mol, t_malha, 85) # (t, x)
        
        erro = calcular_erro_l2(U_mol, U_exata_interp.T)
        erros_mol.append(erro)
        print(f"MOL (dx={dx}): Erro L2 = {erro:.2e}")
        
    graficos.plotar_convergencia_erro(dx_lista, erros_mol, 
                                          "Tarefa 2 - Erro MOL vs Delta_x")
    
    graficos.plotar_solucao_3d(x_mol, t_mol, U_mol, "Tarefa 2 - Solucao MOL")
    print("--- Tarefa 2 Concluída ---")

def run_tarefa_3_fdm(U_exata_grid, x_exato, t_malha):
    print("\n--- Executando Tarefa 3: Diferenças Finitas ---")
    
    M = 21
    T_FINAL = t_malha[-1]
    
    dx = 1.0 / (M - 1)
    
    MU_LISTA = [0.45, 0.55]
    
    for mu in MU_LISTA:
        print(f"\nTestando com mu (r) = {mu:.4f}")
        
        dt = mu * (dx**2)
        N = int(T_FINAL / dt)
        
        x_fdm, t_fdm, U_base, M, N, r, dt, dx = fdm.setup_grid(M, N, mu)
        
        U_base[:, 0] = fourier.f_x(x_fdm)
        
        U_ftcs = fdm.solver_ftcs(U_base.copy(), M, N, r)
        U_btcs = fdm.solver_btcs(U_base.copy(), M, N, r)
        U_cn = fdm.solver_crank_nicholson(U_base.copy(), M, N, r)

        U_exata_fdm = fourier.solucao_fourier(x_fdm, t_fdm, 85) # (t, x)
        
        erro_ftcs = calcular_erro_l2(U_ftcs, U_exata_fdm.T)
        erro_btcs = calcular_erro_l2(U_btcs, U_exata_fdm.T)
        erro_cn = calcular_erro_l2(U_cn, U_exata_fdm.T)
        
        print(f"  Erro FTCS (Progressivo): {erro_ftcs:.2e}")
        print(f"  Erro BTCS (Regressivo): {erro_btcs:.2e}")
        print(f"  Erro CN (Crank-Nich.): {erro_cn:.2e}")
        
        # --- Plots ---
        if not np.isnan(erro_ftcs) and np.isfinite(erro_ftcs):
             graficos.plotar_solucao_3d(x_fdm, t_fdm, U_ftcs, 
                                            f"Tarefa 3 - FTCS (mu={mu:.2f})")
        else:
            print(f"FTCS com mu={mu:.2f} instável. Não foi possível plotar.")

        graficos.plotar_solucao_3d(x_fdm, t_fdm, U_btcs, 
                                        f"Tarefa 3 - BTCS (mu={mu:.2f})")
        graficos.plotar_solucao_3d(x_fdm, t_fdm, U_cn, 
                                        f"Tarefa 3 - CN (mu={mu:.2f})")

    print("--- Tarefa 3 Concluída ---")


if __name__ == "__main__":
    
    L = 1.0
    T_FINAL = 0.5 
    
    DX_EXATO = 0.01
    DT_EXATO = 0.005
    
    x_malha_exata = np.arange(0, L + DX_EXATO, DX_EXATO)
    t_malha_exata = np.arange(0, T_FINAL + DT_EXATO, DT_EXATO)

    U_exata_grid, x_completo = run_tarefa_1_fourier(x_malha_exata, t_malha_exata)
    
    run_tarefa_2_mol(U_exata_grid, x_completo, t_malha_exata)
    
    run_tarefa_3_fdm(U_exata_grid, x_completo, t_malha_exata)

    print("\nSimulação do TP3 concluída. Verifique os arquivos .png gerados.")