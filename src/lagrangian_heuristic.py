'''
Nesse módulo consta a implementação da heurística lagrangiana para o 2-TSP.

TODO: um próximo passo seria generalizar essa heurística para o K-TSP, uma vez 
que implementamos um modelo genérico e atualmente a heurística limita o uso da
relaxação lagrangiana do modelo para K = 2.
'''

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

def lagrangian_heuristic(dist, tours):
    '''
    Dados dois tours com possível interseção nas arestas, remove as arestas
    repetidas e as substitui de forma gulosa por arestas de distâncias baixas.

    Args:
        dist: Dicionário de custo das arestas (i,j), i >= j.
        tours: Lista contendo dois tours.
        
    Returns:
        Tupla com custo dos novos tours e lista dos tours corrigidos.
    '''

    # Arestas ordenadas pelas distâncias
    sorted_edges = [k for k, _ in sorted(dist.items(), key=lambda item: item[1])]

    # Para debugar
    def print_frangments(fragments):
        for v, (u,o) in sorted(fragments.items()):
            print(f"{v} -> {u} ({o}), ",end='')
        print("\n")

    # Função auxiliar que corrige o segundo tour, retorna os tours corrigidos
    # e o custo total.
    def lagrangian_heuristic_assymmetric(tour_fix,tour_wrong):
        tour_corrected = []

        # Arestas de tour_fix
        edges_fix = set(zip(tour_fix,tour_fix[1:] + [tour_fix[0]]))

        # Convertendo representação do tour_wrong para um dicionário. Cada 
        # vértice aponta para o próximo

        ######################################################################
        # Criamos um dicionário para representar os fragmentos corretos do
        # tour_wrong cada vértice aponta para o próximo (se a aresta já existe
        # no tour_fix o vértice aponta para None) e para o primeiro vértice do
        # fragmento ao qual ele pertence. Também construímos um dicionário
        # orphans para indicar os vértices que não são mais apontados por
        # ninguém, os valores correspondem ao último vértice dos fragmentos
        # ao qual esses vértices pertencem.
        fragments = {}
        orphans = {}
        edges = list(zip(tour_wrong,tour_wrong[1:] + [tour_wrong[0]]))

        # Encontra primeiro órfão
        pos = 0
        for (i, j) in edges:
            if (i, j) in edges_fix or (j, i) in edges_fix:
                orphans[j] = j # por enquanto o fragmento termina em j
                first_orphan = j
                break
            pos += 1
        else:
            # Se não achou um órfão, os ciclos já são disjuntos
            tours_new = (tour_fix.copy(), tour_wrong.copy())
            return [tours_cost(dist, tours_new), tours_new]

        orphan = first_orphan
        for i, j in edges[pos+1:] + edges[:pos+1]:
            first = orphan # o fragmento começa no órfão

            if (i, j) in edges_fix or (j, i) in edges_fix: # aresta inválida
                next_vert = None # final do fragmento
                orphan = j # novo órfão
                # Se ainda não sabemos onde o fragmento termina
                if orphan != first_orphan:
                    orphans[j] = j # por enquanto o novo fragmento termina em j
            else: # aresta válida
                next_vert = j
                orphans[orphan] = j # por enquanto o fragmento termina em j

            fragments[i] = [next_vert, first]
        ######################################################################

        # Adiciona arestas até limpar o conjunto orphans, não permite ciclos
        # menores que n
        for i, j in sorted_edges:
            if len(orphans) == 1: # só existe uma possibilidade de aresta
                beg,end = list(orphans.items())[0]
                fragments[end][0] = beg
                break

            # Se (i,j) é uma candidata para adicionar
            if not (i, j) in edges_fix and not (j, i) in edges_fix:

                # Se (i,j) fecha um fragmento
                if fragments[i][0] == None and j in orphans:
                    def try_add_edge(i,j):
                        # Se não vai formar um ciclo de tamanho menor que n
                        if fragments[i][1] != j:
                            fragments[i][0] = j # Adiciona a aresta
                            # Junta os fragmentos
                            fragments[orphans[j]][1] = fragments[i][1]
                            orphans[fragments[i][1]] = orphans[j]
                            orphans.pop(j) # j não é mais órfão
                    try_add_edge(i,j)

                # Grafo é não direcionado, precisamos ver (j,i) também
                elif fragments[j][0] == None and i in orphans:
                    try_add_edge(j,i)

        # Retorna o dicionário para uma lista
        tour_corrected = [0]
        v = fragments[0][0]
        while v != 0:
            tour_corrected.append(v)
            v = fragments[v][0]

        tours_new = (tour_fix.copy(), tour_corrected)
        return [tours_cost(dist, tours_new), tours_new]

    # A ordem em que processamos os tours afeta o resultado,
    # portanto testamos todas as ordens e escolhemos a melhor
    (cost01,tours01) = lagrangian_heuristic_assymmetric(tours[0],tours[1])
    (cost10,tours10) = lagrangian_heuristic_assymmetric(tours[1],tours[0])
    if cost01 < cost10:
        return (cost01,tours01)
    else:
        return (cost10,tours10)
