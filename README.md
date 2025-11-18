# Trabalho Prático 3 – Equação do Calor (EDN)

Este projeto implementa diferentes métodos numéricos para resolver a equação do calor no intervalo (0,1) com condições de fronteira homogêneas. O trabalho segue o enunciado da disciplina de **Equações Diferenciais Numéricas (EDN)** da UFMG.

## 📌 Conteúdo do Projeto

O projeto contém implementações das três abordagens solicitadas:

### 1. Séries de Fourier
Utilizada como solução de referência (solução "exata" numérica).
- Implementação dos coeficientes bn
- Truncamento otimizado para erro ≤ 0.01
- Comparação com métodos numéricos

### 2. Método das Linhas (MOL)
Discretização espacial + integração temporal via **solve_ivp (RK45)**.
- Convergência de ordem 2 observada
- Testes com vários Δx

### 3. Diferenças Finitas
Implementação dos métodos:
- **FTCS (explícito)** – condicionalmente estável
- **BTCS (implícito)** – sempre estável
- **Crank–Nicholson** – estável e de segunda ordem

Foram testados dois valores de:
\(\mu = \Delta t/(\Delta x)^2\)
- μ = 5/11 (estável)
- μ = 5/9 (instável para FTCS)

## 📂 Estrutura do Repositório

```
tp3_EDN/
│
├── codigo/
│   ├── fourier.py
│   ├── mol.py
│   ├── ftcs.py
│   ├── btcs.py
│   ├── crank_nicholson.py
│   └── utils.py
│
├── imagens/
│   ├── tarefa_1_-_fourier_snapshots.png
│   ├── tarefa_2_-_erro_mol_vs_delta_x.png
│   ├── tarefa_3_-_ftcs_(mu=0.45).png
│   ├── tarefa_3_-_ftcs_(mu=0.55).png
│   ├── tarefa_3_-_btcs_(mu=0.45).png
│   ├── tarefa_3_-_btcs_(mu=0.55).png
│   ├── tarefa_3_-_cn_(mu=0.45).png
│   └── tarefa_3_-_cn_(mu=0.55).png
│
├── relatorio/
│   └── trabalho.tex
│
└── README.txt   (este arquivo)
```

## ▶️ Como Executar

1. Instale dependências:
```
pip install numpy scipy matplotlib
```

2. Execute qualquer método, por exemplo:
```
python fourier.py
python mol.py
python ftcs.py
```

3. As imagens são geradas automaticamente na pasta `/imagens`.

## 📊 Resultados

O relatório discute:
- precisão de cada método,
- estabilidade,
- comparação com a solução de Fourier,
- influência do parâmetro μ.

## 👨‍💻 Autor
Mateus Ryan de Castro Lima  
Universidade Federal de Minas Gerais – UFMG  
GitHub: https://github.com/Mryan31/tp3_EDN

