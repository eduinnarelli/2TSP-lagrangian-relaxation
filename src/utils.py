'''
Funções auxiliares compartilhadas por alguns métodos.
'''
import gurobipy as gp

def shortest_cycle(n, edges):
    '''
    Função que constrói o menor ciclo de um conjunto de arestas, em termos
    do número de vértices no ciclo.

    Args:
        n: nº de vértices.
        edges: lista de tuplas de arestas.
        
    Returns:
        Lista de vértices no menor ciclo, cada um conectado com o anterior
        e próximo da lista.
    '''

    unvisited = list(range(n))

    # Tamanho inicial tem um nó a mais, p/ forçar atualização
    cycle = range(n + 1) 

    while unvisited:  # 'True' enquanto pilha for não-vazia
        thiscycle = []
        neighbors = unvisited

        # Construir ciclo
        while neighbors:
            current = neighbors[0]
            thiscycle.append(current)
            unvisited.remove(current)
            # 'select(current, '*')' retorna todos vizinhos de 'current'
            neighbors = [
                j 
                for i, j in edges.select(current, '*')
                if j in unvisited
            ]

        # Atualizar menor ciclo, se preciso
        if len(cycle) > len(thiscycle):
            cycle = thiscycle

    return cycle

def get_edges_in_tour(tour_id, x_sol, all_edges):
    '''
    Função que, dado o conjunto de variáveis na solução, retorna as arestas
    (i, j) que estão presentes em uma rota.

    Args:
        tour_id: identificador da rota.
        x_sol: valores de 'x' (variáveis que indicam presença das arestas) na 
            solução.
        all_edges: lista de arestas (i,j) no grafo de entrada.

    Returns:
        Lista de arestas (i,j) na rota 'tour_id'.
    '''

    return gp.tuplelist(
        (i, j) 
        for i, j, k in all_edges
        if x_sol[i, j, k] > 0.5 and k == tour_id
    )       


def build_tours_in_sol(K, n, x_sol, all_edges):
    '''
    Função que constrói as rotas correspondentes aos valores das variáveis
    'x' na solução.

    Args:
        K: nº de caixeiros viajantes.
        n: nº de vértices do grafo.
        x_sol: valores de 'x' (variáveis que indicam presença das arestas) na 
            solução.
        all_edges: lista de arestas (i,j) no grafo de entrada.
    
    Returns:
        Lista de rotas na solução.
    '''

    tours = []
    for t in range(K):
        edges_in_tour = get_edges_in_tour(t, x_sol, all_edges) 
        tours.append(shortest_cycle(n, edges_in_tour))

    return tours

def print_solution(K, tours, cost):
    '''
    Função que imprime uma solução no stdout.

    Args:
        K: nº de caixeiros viajantes.
        tours: lista de rotas.
        cost: custo da solução.
    '''

    for t in range(K):
        print(f'Rota {t}: {tours[t]}')
    print(f'Custo: {cost}\n')

def next_permutation(l, key):
    '''
    Encontra a próxima permutação dos elementos de l em ordem lexicografica.
    Modifica a lista de entrada

    Args:
        l: lista para ser modificada
        key: função para obter a chave usada na comparação

    Returns:
        True se encontrou a próxima permutação, False caso já é a ultima
    '''

    if len(l) <= 1:
        return False

    for i in range(len(l) - 2, -1, -1):
        j = i + 1
        if key(l[i]) < key(l[j]):
            k = len(l) - 1
            while not (key(l[i]) < key(l[k])):
                k -= 1
            l[i], l[k] = l[k], l[i]
            l[j:] = l[-1:j-1:-1]
            return True
    return False
