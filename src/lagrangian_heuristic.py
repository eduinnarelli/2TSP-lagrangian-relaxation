'''
Nesse módulo consta a implementação da heurística lagrangiana para o K-TSP.
'''

import math
from disjoint_set import DisjointSet
from itertools import permutations

def tours_cost(dist, tours):
    '''
    Calcula a distância de k tours.

    Args:
        dist: Dicionário de custo das arestas (i,j), i >= j.
        tours: Lista contendo k tours.

    Returns:
        Distancia total dos tours.
    '''

    sum = 0
    for tour in tours:
        edges = zip(tour,tour[1:] + [tour[0]])
        for i, j in edges:
            sum += dist[max(i,j),min(i,j)]

    return sum

def lagrangian_heuristic(dist, tours, n):
    '''
    Dados k tours com possível interseção nas arestas, remove as arestas
    repetidas e as substitui de forma gulosa por arestas de distâncias baixas.

    Args:
        dist: Dicionário de custo das arestas (i,j), i >= j.
        tours: Lista contendo k tours.

    Returns:
        Tupla com custo dos novos tours e lista dos tours corrigidos.
    '''

    # Arestas ordenadas pelas distâncias
    sorted_edges = [k for k, _ in sorted(dist.items(), key=lambda item: item[1])]

    # Função auxiliar que corrige o segundo tour, retorna os tours corrigidos
    # e o custo total.
    def lagrangian_heuristic_assymmetric(tours):
        tours_new = []

        # O primeiro tour é fixo
        tour_fix = tours[0]
        tours_new.append(tour_fix.copy())
        # Arestas de tour_fix
        edges_fix = set(zip(tour_fix,tour_fix[1:] + [tour_fix[0]]))

        # Conserta cada um dos tours restantes
        pos = 1
        sorted_edges_local = sorted_edges
        while pos < len(tours):
            tour_wrong = tours[pos]

            # Estruturas para representar os fragmentos corretos do tour_wrong:
            # Vizinho de cada vértice (conectados por arestas que não estejam 
            # fixas)
            neighbors = {}
            # Conjunto para indicar que um vértice não tem as duas arestas
            orphans = set()
            # Estrutura de conjunto disjunto para representar a que fragmento 
            # um vértice pertencem
            ds = DisjointSet()

            for i in range(n):
                ds.find(i)
                neighbors[i] = []

            edges = list(zip(tour_wrong,tour_wrong[1:] + [tour_wrong[0]]))
            for (i, j) in edges: # Para cada aresta do tour_wrong
                # Se a aresta é inválida
                if (i, j) in edges_fix or (j, i) in edges_fix:
                    # i e j são órfãos
                    orphans.add(j)
                    orphans.add(i)
                else:
                    # adiciona aresta entre i e j
                    neighbors[i].append(j)
                    neighbors[j].append(i)
                    # i e j pertencem ao mesmo fragmento
                    ds.union(i,j)

            # Adiciona arestas até limpar o conjunto orphans,
            # não permite ciclos menores que n
            for (i, j) in sorted_edges:

                # Acabaram os órfãos
                if len(orphans) == 0:
                    break

                # Se (i,j) é uma candidata para adicionar
                if not (i, j) in edges_fix and not (j, i) in edges_fix and\
                   i in orphans and j in orphans:

                    # Se não vai formar um ciclo de tamanho menor que n
                    if not ds.connected(i,j) or len(orphans) <= 2:
                        # adiciona aresta entre i e j
                        neighbors[i].append(j)
                        neighbors[j].append(i)
                        ds.union(i,j)
                        # i e j não são mais órfãos
                        # a menos que só tenham um vizinho
                        if len(neighbors[i]) > 1:
                            orphans.discard(i)
                        if len(neighbors[j]) > 1:
                            orphans.discard(j)

            # Existem casos onde não podemos fechar o ciclo com o ordem 
            # de arestas dadas, esses casos são raros com K=2, mas quanto maior
            # o K maior a chance disso ocorrer. Nesses casos mudamos a ordem
            # e tentamos de novo.
            if len(orphans) > 0:
                sorted_edges_local = sorted_edges_local[1:] + [sorted_edges_local[0]]

                # Caso não seja possível formar o ciclo (K muito grande com
                # poucos vértices).
                if sorted_edges_local[0] ==\
                    max(dist.items(), key=lambda item: item[1])[0]:
                    raise Exception('Impossível formar K ciclos disjuntos')
                continue
            else:
                sorted_edges_local = sorted_edges
                pos += 1

            # Retorna o dicionário para uma lista
            tour_corrected = [0]
            u = 0
            v = neighbors[0][0]
            while v != 0:
                tour_corrected.append(v)
                if neighbors[v][0] == u:
                    u = v
                    v = neighbors[v][1]
                else:
                    u = v
                    v = neighbors[v][0]

            # Fixa o novo tour
            edges_fix = edges_fix.union(set(zip(tour_corrected,tour_corrected[1:] + [tour_corrected[0]])))

            tours_new.append(tour_corrected)

        return [tours_cost(dist, tours_new), tours_new]

    # A ordem em que processamos os tours afeta o resultado,
    # portanto testamos todas as ordens e escolhemos a melhor
    tours_best = []
    cost_best = math.inf
    for tours in permutations(tours):
        (cost_new,tours_new) = lagrangian_heuristic_assymmetric(tours)
        if cost_new < cost_best:
            tours_best = tours_new
            cost_best = cost_new

    return (cost_best,tours_best)
