# Otimização do custo de um tanque cilíndrico

Trabalho computacional de Modelagem Computacional (UERJ, Instituto Politécnico, IPRJ01-10760/IPRJ01-11873), prof. Gustavo Barbosa Libotte.

Autores: Ana Clara da Silva Vaz e Danillo Silva Linhares.

## Problema

Minimizar o custo de fabricação de um tanque cilíndrico (base R$ 150/m², parede lateral R$ 100/m², teto R$ 200/m²):

    f(r, h) = 350·π·r² + 200·π·r·h

sujeito ao volume mínimo de 500 m³:

    g(r, h) = 500 − π·r²·h ≤ 0

## Método

- Penalização exterior quadrática: `P(r, h) = f(r, h) + μ·[max(0, g(r, h))]²`
- Derivadas parciais por diferenças finitas adiantadas, com Δ = 1e-5
- Método do gradiente: `x(k+1) = x(k) − α·∇P(x(k))`, com μ e α constantes
- Critério de parada: norma do gradiente < 1,0 (máximo de 3.000.000 iterações)

## Como executar

Requer Python 3.10 ou superior.

    pip install -r requirements.txt
    python trab_mc.py

O programa pergunta o número de execuções com chutes iniciais aleatórios (ex.: 10), imprime os resultados e mostra os gráficos de convergência para diferentes valores de μ e α.

## Parâmetros principais

| Parâmetro | Valor |
|-----------|-------|
| μ (penalidade) | 100 |
| α (passo) | 8e-8 |
| Δ (diferença finita) | 1e-5 |
| Chute inicial | aleatório, r ∈ [0,1; 10] e h ∈ [0,1; 20] |

## Resultados

Solução analítica de referência: r ≈ 3,5693 m, h ≈ 12,4926 m, custo ≈ R$ 42.024,95.

Solução numérica (μ = 100, α = 8e-8): r ≈ 3,5595 m, h ≈ 12,5545 m, V ≈ 499,72 m³, custo ≈ R$ 42.009,47.

O volume fica ligeiramente abaixo de 500 m³ porque a penalização exterior admite uma pequena violação da restrição, que diminui com o aumento de μ.

## Arquivos

- `trab_mc.py`: código-fonte
- `requirements.txt`: dependências
- `Relatorio_Tanque_Cilindrico.pdf`: relatório completo
