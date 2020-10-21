'''
Nesse módulo consta a implementação da heurística lagrangiana para o K-TSP.
'''

import math
from disjoint_set import DisjointSet
from itertools import permutations
from utils import next_permutation

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
            # Tamanho do fragmentos
            frag_sizes = {}
            # Estrutura de conjunto disjunto para representar a que fragmento 
            # um vértice pertencem
            ds = DisjointSet()

            for i in range(n):
                ds.find(i)
                frag_sizes[i] = 1
                neighbors[i] = []
                orphans.add(i)

            # Para adicionar arestas
            def add_edge(i,j):
                # adiciona aresta entre i e j
                neighbors[i].append(j)
                neighbors[j].append(i)
                # i e j pertencem ao mesmo fragmento
                ds.union(i,j)
                frag_sizes[ds.find(i)] = frag_sizes[ds.find(i)] + frag_sizes[ds.find(j)]
                # i e j não são mais órfãos
                # a menos que só tenham um vizinho
                if len(neighbors[i]) > 1:
                    orphans.discard(i)
                if len(neighbors[j]) > 1:
                    orphans.discard(j)

            # Para verificar se (i,j) não esta fixo
            def valid(i,j):
                return (not (i, j) in edges_fix and not (j, i) in edges_fix)

            # Para verificar se adicionar uma aresta não abriga a adição de uma
            # aresta invalida
            def force_invalid(i,j):
               new_size = frag_sizes[ds.find(i)] + frag_sizes[ds.find(j)]

               if new_size >= n - 1:
                   orphans.discard(i)
                   orphans.discard(j)
                   l = list(orphans)
                   if len(l) == 1:
                       if not valid(i,l[0]) or not valid(j,l[0]):
                           return True
                   elif not valid(l[0], l[1]):
                       return True
                   orphans.add(i)
                   orphans.add(j)

               return False

            edges = list(zip(tour_wrong,tour_wrong[1:] + [tour_wrong[0]]))
            for (i, j) in edges: # Para cada aresta do tour_wrong
                # Se a aresta é válida
                if valid(i,j):
                    add_edge(i,j)

            # Adiciona arestas até limpar o conjunto orphans,
            # não permite ciclos menores que n
            for (i, j) in sorted_edges:

                # Acabaram os órfãos
                if len(orphans) == 0:
                    break

                # Se (i,j) conecta órfãos
                if i in orphans and j in orphans:
                    # Se (i,j) é uma aresta válida, não vai formar um ciclo de
                    # tamanho menor que n e não vai obrigar a adição de uma
                    # aresta invalida
                    if valid and ( len(orphans) <= 2 or
                          not ds.connected(i,j) and not force_invalid(i,j)):
                        add_edge(i,j)

            # Para k>2, existem mais casos onde não podemos fechar o ciclo com
            # a ordem de arestas dadas. Nesses casos mudamos a ordem e tentamos
            # de novo.
            if len(orphans) > 0:
                assert len(tours) > 2
                sorted_edges_local = sorted_edges.copy()
                if not next_permutation(sorted_edges_local, key=lambda item: item[1]):
                    # Caso não seja possível formar o ciclo (K muito grande com
                    # poucos vértices).
                    raise Exception('Impossível formar K ciclos disjuntos')
                continue
            else:
                pos += 1
                if pos < len(tours):
                    sorted_edges_local = sorted_edges.copy()

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
