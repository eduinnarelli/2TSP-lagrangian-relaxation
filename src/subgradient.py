'''
Nesse módulo consta a implementação do método do subgradiente para a resolução
do problema dual lagrangiano relativo à relaxação lagrangiana do K-TSP.
'''

import gurobipy as gp
from lagrangian_heuristic import lagrangian_heuristic
from subtour_elimination import subtour_elimination
from utils import build_tours_in_sol, shortest_cycle

def subgradient(model, sgvars, dist):
    '''
    Método do subgradiente que visa encontrar os multiplicadores de lagrange
    que otimizam o limitante inferior retornado pela relaxação lagrangiana do 
    K-TSP.

    Args:
        model: modelo do K-TSP.
        sgvars: variável associada ao subgradiente.
        dist: dicionário de custo das arestas (i,j), i >= j.

    Returns:
        Dicionário da solução com melhores limitantes inferior e superior 
        ('best_lb' e 'best_ub') obtidos pelo método, assim como o tempo total
        de execução 'runtime'.
    '''

    pi = 2.0
    max_iterations = 100
    runtime = 0.0
    best_ub = {'cost': float('inf')}

    # Recuperar alguns atributos salvos no modelo
    n = model._n
    K = model._K
    xvars = model._xvars
  
    # Inicializar multiplicadores com 0 (um para cada aresta)
    u = {(i,j): 0 for i in range(n) for j in range(i)}

    # Função objetivo original
    original_obj = gp.quicksum(
        dist[i,j] * xvars[i,j,k] for i in range(n) 
                                 for j in range(i) 
                                 for k in range(K)
    )

    for i in range(max_iterations):
        # Penalidades correspondentes às restrições dualizadas
        obj_penalty = gp.quicksum(
            u[i,j] * sgvars[i,j] for i in range(n) for j in range(i)
        )

        # Penalizar função objetivo e re-otimizar
        model.setObjective(original_obj + obj_penalty)
        model.optimize(subtour_elimination)
        runtime += model.Runtime

        # Recuperar solução
        x_sol = model.getAttr('x', xvars)
        tours = build_tours_in_sol(K, n, x_sol, xvars.keys())
        
        # Executar heurística lagrangiana para obter um limitante superior
        heuristic_sol = lagrangian_heuristic(dist, tours)
        ub = {'cost': heuristic_sol[0], 'tours': heuristic_sol[1]}

        # Atualizar melhor limitante superior, se necessário
        if ub['cost'] < best_ub['cost']:
            best_ub = ub

        # CRITÉRIOS DE PARADA:
        # - Optimalidade
        # - Fim das iterações
        opt_gap = (best_ub['cost'] - model.objVal) / best_ub['cost']
        if opt_gap < 10e-6 or i == (max_iterations - 1):
            break

        # Recuperar subgradiente
        sg_sol = model.getAttr('x', sgvars)

        # Denominador do passo é a soma dos quadrados dos valores do
        # subgradiente
        square_subgrad_sum = sum(
            [sg_sol[i,j]**2 for i in range(n) for j in range(i)]
        )

        # Atualizar multiplicadores, dando um passo em direção ao 
        # subgradiente com o intuito de maximizar o limitante inferior
        # retornado pela relaxação
        step = pi * (best_ub['cost'] - model.objVal) / square_subgrad_sum
        u = {
            (i,j): max(0.0, u[i,j] + step * sg_sol[i,j]) 
            for i in range(n) for j in range(i)
        }

    # Retornar dicionário com melhores limitantes encontrados e tempo de
    # execução total do método
    return {
        'best_lb': {'cost': model.objVal, 'tours': tours},
        'best_ub': best_ub,
        'runtime': runtime,
    }
