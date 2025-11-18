import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import os

SAVE_DIR = os.path.join('tp3_equacao_calor', 'graficos')

os.makedirs(SAVE_DIR, exist_ok=True)


def plotar_solucao_3d(x, t, U, titulo):
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')
    
    T, X = np.meshgrid(t, x)
    
    surf = ax.plot_surface(X, T, U, cmap=cm.viridis,
                           linewidth=0, antialiased=False)
    
    ax.set_title(titulo)
    ax.set_xlabel('x (espaço)')
    ax.set_ylabel('t (tempo)')
    ax.set_zlabel('u(x,t) (temperatura)')
    fig.colorbar(surf, shrink=0.5, aspect=5)
    
    save_name = f"{titulo.replace(' ', '_').lower()}.png"
    save_path = os.path.join(SAVE_DIR, save_name)
    
    plt.savefig(save_path)
    print(f"Gráfico 3D salvo em {save_path}")
    plt.close(fig)

def plotar_snapshots(x, U, t_vec, t_indices, titulo):
    """
    Plota u(x) em diferentes instantes de tempo t.
    U é (M, N)
    t_indices é uma lista de índices de tempo para plotar (ex: [0, 10, 50])
    """
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111)
    
    for i in t_indices:
        if i < len(t_vec):
            ax.plot(x, U[:, i], label=f't = {t_vec[i]:.2f}')
            
    ax.set_title(titulo)
    ax.set_xlabel('x (espaço)')
    ax.set_ylabel('u(x,t)')
    ax.legend()
    ax.grid(True)
    
    save_name = f"{titulo.replace(' ', '_').lower()}_snapshots.png"
    save_path = os.path.join(SAVE_DIR, save_name)

    plt.savefig(save_path)
    print(f"Gráfico de snapshots salvo em {save_path}")
    plt.close(fig)

def plotar_convergencia_erro(dx_lista, erros, titulo):
    fig = plt.figure(figsize=(8, 5))
    ax = fig.add_subplot(111)
    
    ax.plot(dx_lista, erros, 'o-', label='Erro L2')
    ax.set_title(titulo)
    ax.set_xlabel(r'$\Delta x$')
    ax.set_ylabel('Erro (log)')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid(True, which="both", ls="--")
    ax.legend()
    
    save_name = f"{titulo.replace(' ', '_').lower()}_erro.png"
    save_path = os.path.join(SAVE_DIR, save_name)
    
    plt.savefig(save_path)
    print(f"Gráfico de erro salvo em {save_path}")
    plt.close(fig)