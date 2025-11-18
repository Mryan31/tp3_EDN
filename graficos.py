import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import os # 1. Importar a biblioteca 'os'

# 2. Definir o caminho da pasta onde os gráficos serão salvos
# O caminho é relativo ao local de onde você roda o script (C:\EDN)
SAVE_DIR = os.path.join('tp3_equacao_calor', 'graficos')

# 3. Criar a pasta, se ela não existir
os.makedirs(SAVE_DIR, exist_ok=True)


def plotar_solucao_3d(x, t, U, titulo):
    """
    Gera um gráfico 3D da solução u(x,t).
    U deve ser (M, N) -> (pontos em x, pontos em t)
    """
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
    
    # 4. Modificar o 'savefig'
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
    
    # 4. Modificar o 'savefig'
    save_name = f"{titulo.replace(' ', '_').lower()}_snapshots.png"
    save_path = os.path.join(SAVE_DIR, save_name)

    plt.savefig(save_path)
    print(f"Gráfico de snapshots salvo em {save_path}")
    plt.close(fig)

def plotar_convergencia_erro(dx_lista, erros, titulo):
    """
    Plota o erro em função de Delta_x (ou outro parâmetro).
    """
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
    
    # 4. Modificar o 'savefig'
    save_name = f"{titulo.replace(' ', '_').lower()}_erro.png"
    save_path = os.path.join(SAVE_DIR, save_name)
    
    plt.savefig(save_path)
    print(f"Gráfico de erro salvo em {save_path}")
    plt.close(fig)