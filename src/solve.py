'''
Nesse módulo consta uma função que modela e resolve o K-TSP (generalização do 
2-TSP) e a relaxação lagrangiana do 2-TSP, bem como chamadas a essa função
para análise do 2-TSP. 
'''

import sys
import argparse
import math
import pickle
import random
from itertools import combinations
from tqdm import tqdm
import gurobipy as gp
from gurobipy import GRB
from subtour_elimination import subtour_elimination
from subgradient import subgradient
from utils import build_tours_in_sol, shortest_cycle, print_solution

def k_tsp(K, n, dist, relaxed=False):
    '''
    Função que define e resolve o modelo exato ou relaxado para o K-TSP, dada uma 
    determinada instância. Aqui, K-TSP generaliza o TSP e o 2-TSP para qualquer K, 
    o que evita a implementação de modelos diferentes.

    Args:
        K: nº de caixeiros viajantes.
        n: nº de vértices do grafo.
        dist: dicionário de custo das arestas (i,j), i >= j.
        relaxed: booleano que indica se será resolvido o modelo original ou a
            relaxação lagrangiana.

    Returns:
        Dicionário da solução, contendo a solução ótima se resolvido o problema
        original ou os melhores limitantes se resolvida a relaxação lagrangiana.
    '''

    # Inicializar ambiente
    env = gp.Env(empty = True)
    env.setParam('OutputFlag', 0)
    env.start()

    # Inicializar modelo
    model = gp.Model(name = str(K) + '-tsp', env = env)

    # Adaptar o dicionário de distâncias de acordo com a quantidade de 
    # caixeiros
    distK = {
        (i, j, k):  dist[i, j] 
                    for i in range(n) for j in range(i) for k in range(K)
    }

    # Criar variáveis
    xvars = model.addVars(distK.keys(), obj=distK, vtype=GRB.BINARY, name='x')
    for i, j, k in xvars.keys():
        xvars[j, i, k] = xvars[i, j, k]  # grafo não-orientado

    # Restrições de grau 2, p/ cada rota k
    model.addConstrs(
        (xvars.sum(i, '*', k) == 2 for i in range(n) for k in range(K)), 
        name='deg-2'
    )

    # Salvar alguns atributos no modelo para acessá-los facilmente na callback
    model._n = n
    model._K = K
    model._xvars = xvars

    # Indicar limite de tempo da otimização e callback a ser chamada após a
    # solução ótima do modelo relaxado ser encontrada
    model.Params.lazyConstraints = 1
    model.Params.timeLimit = 1800.0

    # Restrições de disjunção entre arestas de diferentes rotas são incluídas no 
    # modelo do problema original...
    if not relaxed:

        # Incluir restrições e otimizar
        model.addConstrs(
            (xvars.sum(i, j, '*') <= 1 for i in range(n) for j in range(i)), 
            name='disj'
        )
        model.optimize(subtour_elimination)

        # Recuperar solução
        x_sol = model.getAttr('x', xvars)
        tours = build_tours_in_sol(K, n, x_sol, xvars.keys())

        # Retornar dicionário com solução ótima e tempo de execução
        return {
            'opt': {'cost': model.objVal, 'tours': tours},
            'runtime': model.Runtime,
        }

    # ... e dualizadas na Relaxação Lagrangiana
    else:

        # Criar variáveis para o subgradiente
        sgvars = model.addVars(dist.keys(), vtype=GRB.BINARY, name='sg')
        for i, j in sgvars.keys():
            sgvars[j, i] = sgvars[i, j]  # grafo não-orientado

        # As novas variáveis são associadas às restrições dualizadas, o que 
        # facilita na manipulação e extração desses valores
        model.addConstrs(
            (
                sgvars[i,j] == - 1 + xvars.sum(i, j, '*') 
                for i in range(n) for j in range(i)
            ),
            name='dualized'
        )

        # Resolver método do subgradiente
        return subgradient(model, sgvars, dist)

# O usuário indica pela linha do comando se ele deseja resolver o
# 2-TSP de forma exata ou relaxada
parser = argparse.ArgumentParser()
parser.add_argument('--relaxed', default=False, action='store_true')
args = parser.parse_args()
relaxed = vars(args)['relaxed']

# Carregar instâncias salvas em 'fixed_instances.pkl'
with open("instances/fixed_instances.pkl", "rb") as fp:
    instances = pickle.load(fp)

dash = '===================='

# Salvar output em um txt
output_name = 'rel_output.txt' if relaxed else 'opt_output.txt'
print(f"Soluções serão salvas em '{output_name}'")
sys.stdout = open('outputs/' + output_name, 'w')

for instance in tqdm(instances):
    n = instance['n']
    dist = instance['dist']

    # Resolver 2-TSP de forma exata ou relaxada
    sol_type = 'RELAXADA' if relaxed else 'EXATA'
    print(f'\n{dash} SOLUÇÃO {sol_type} DO 2-TSP PARA N = {n} {dash}\n')
    sol = k_tsp(2, n, dist, relaxed=relaxed)

    if relaxed:
        # Imprimir limitantes da relaxação lagrangiana
        print('Melhor limitante inferior encontrado:')
        print_solution(2, sol['best_lb']['tours'], sol['best_lb']['cost'])
        print('Melhor limitante superior encontrado:')
        print_solution(2, sol['best_ub']['tours'], sol['best_ub']['cost'])

    else:
        # Imprimir solução ótima
        print_solution(2, sol['opt']['tours'], sol['opt']['cost'])

    print(f"Tempo de execução: {sol['runtime']}s")

sys.stdout.close()
