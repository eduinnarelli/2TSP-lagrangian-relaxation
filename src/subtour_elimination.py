'''
Nesse módulo consta a implementação da callback que adiciona dinamicamente 
cortes relativos às restrições de eliminação de subciclo ao modelo do K-TSP.
'''
from itertools import combinations
from tqdm import tqdm
import gurobipy as gp
from gurobipy import GRB
from utils import get_edges_in_tour, shortest_cycle

def subtour_elimination(model, where):
    '''
    Callback que, para uma solução ótima do K-TSP relaxado, verifica se essa
    solução viola restrições de eliminação de subciclo e, se sim, adiciona
    essas restrições ao modelo, que será re-otimizado.

    Args:
        model: o modelo associado a callback.
        where: indica da onde no processo de otimização a callback foi chamada.
    '''

    if where == GRB.Callback.MIPSOL:
        # Analisar cada rota t
        for t in range(model._K):

            # Criar lista de arestas na rota t selecionadas na solução
            x_sol = model.cbGetSolution(model._xvars)
            edges_in_tour = get_edges_in_tour(t, x_sol, model._xvars.keys())

            # Encontrar menor ciclo e verificar se viola restrição, i.e., se
            # não percorre todos os vértices, formando um subciclo
            cycle = shortest_cycle(model._n, edges_in_tour)
            if len(cycle) < model._n:

                # Adicionar restrições de eliminação de subciclo, para cada par 
                # de vértices do subciclo encontrado
                model.cbLazy(
                    gp.quicksum(
                        model._xvars[i, j, t]
                        for i, j in combinations(cycle, 2)
                    ) <= len(cycle)-1
                )

