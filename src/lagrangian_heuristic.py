'''
Nesse módulo consta a implementação da heurística lagrangiana para o K-TSP.
'''

import math
from itertools import permutations
from heapq import *

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
            sum += get_dist(dist,i,j)

    return sum

def get_dist(dist, i, j):
    return dist[max(i,j),min(i,j)]

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

    # Recupera as arestas de um tour
    def to_edges(tour):
        return zip(tour,tour[1:] + [tour[0]])

    # Função auxiliar que corrige o segundo tour, retorna os tours corrigidos
    # e o custo total.
    def lagrangian_heuristic_assymmetric(tours):
        tours_new = []

        # O primeiro tour é fixo
        tour_fix = tours[0]
        tours_new.append(tour_fix.copy())
        # Arestas de tour_fix
        edges_fix = set(to_edges(tour_fix))

        # Conserta cada um dos tours restantes
        for tour_wrong in tours[1:]:

            # Estruturas para representar os fragmentos corretos do tour_wrong:
            # Vizinho de cada vértice (conectados por arestas que não estejam 
            # fixas)
            neighbors = {}
            # Conjunto de arestas invalidas
            invalid_edges = set()
            # Conjunto de arestas validas
            valid_edges = set()
            # Pares de arestas candidatas para correção
            candidates = []

            for i in range(n):
                neighbors[i] = []

            # Para adicionar arestas
            def add_edge(i,j):
                # adiciona aresta entre i e j
                neighbors[i].append(j)
                neighbors[j].append(i)
                # Se i e j são validas
                if valid(i,j):
                    valid_edges.add((i,j))
                else:
                    invalid_edges.add((i,j))

            # Para gerar canditatos
            def find_canditates(i,j):
                if not valid(i,j):
                    for (k,l) in invalid_edges:
                        if not i in [k,l] and not j in [k,l]:
                            heappush(candidates,
                                     (get_dist(dist,i,k) +
                                      get_dist(dist,j,l),(i,j),(k,l),False))
                            heappush(candidates,
                                     (get_dist(dist,i,l) +
                                      get_dist(dist,j,k),(i,j),(k,l),True))

            # Para remover arestas
            def remove_edge(i,j):
                # adiciona aresta entre i e j
                neighbors[i].remove(j)
                neighbors[j].remove(i)
                # Se i e j são validas
                if valid(i,j):
                    valid_edges.discard((i,j))
                else:
                    invalid_edges.discard((i,j))

            # Para encontar o final de um caminho
            def end_of_path(i):
                u = i
                v = neighbors[i][0]
                while len(neighbors[v]) > 1:
                    if neighbors[v][0] == u:
                        u = v
                        v = neighbors[v][1]
                    else:
                        u = v
                        v = neighbors[v][0]
                return v

            # Para verificar se (i,j) não esta fixo
            def valid(i,j):
                return (not (i, j) in edges_fix and not (j, i) in edges_fix)

            edges = list(to_edges(tour_wrong))
            for (i, j) in edges: # Para cada aresta do tour_wrong
                add_edge(i,j)
                find_canditates(i,j)

            # Adiciona arestas até limpar o conjunto arestas invalidas,
            # não permite ciclos menores que n

            # Remove pares de arestas invalidas não consecutivas,
            # da preferencia para adicionar arestas com custos menores
            skip = set() # Candidatos deixados para depois
            while len(candidates) > 0:

                (_, (i,j), (k,l), flip) = heappop(candidates)

                # Uma das arestas já foi substituida
                if not (i,j) in invalid_edges or not (k,l) in invalid_edges:
                    continue

                # Trocar as arestas forma um ciclo
                remove_edge(i,j)
                remove_edge(k,l)
                if (end_of_path(i) == l) == flip:
                    if (min((i,j),(k,l)), max((i,j),(k,l))) in skip:
                        flip = not flip
                    else:
                        skip.add((min((i,j),(k,l)), max((i,j),(k,l))))
                        add_edge(i,j)
                        add_edge(k,l)
                        continue

                # Substitue ambas as arestas
                # (pelo menos uma das novas é valida)
                if flip:
                    add_edge(i,l)
                    find_canditates(i,l)
                    add_edge(j,k)
                    find_canditates(j,k)
                else:
                    add_edge(i,k)
                    find_canditates(i,k)
                    add_edge(j,l)
                    find_canditates(j,l)

            # Se não possue mais arestas invalidas não consecutivas...
            while len(invalid_edges) > 0:
                (i,j) = next(iter(invalid_edges))
                # ...encontra uma valida para remover
                for (k,l) in iter(valid_edges):
                    if not i in [k,l] and not j in [k,l]:
                        remove_edge(i,j)
                        remove_edge(k,l)
                        if valid(i,k) and valid(j,l) and\
                                end_of_path(i) == l:
                            add_edge(i,k)
                            add_edge(j,l)
                            break
                        if valid(i,l) and valid(j,k) and\
                                end_of_path(i) == k:
                            add_edge(i,l)
                            add_edge(j,k)
                            break
                        add_edge(i,j)
                        add_edge(k,l)

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

        assert len(tour_corrected) == n
        assert set(to_edges(tour_corrected)).isdisjoint(set(to_edges(tour_fix)))
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
