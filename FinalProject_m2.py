# 2022년도 2학기 데이터마이닝 - FinalProject_m2
# 이름 : 이지윤, 정유라, 김송우
# 학번 : 2021247052, 2021247010, 2020253093

import sys
import time
from itertools import combinations
import networkx as nx
import numpy as np
from F1score import get_F1score


# 0. get data ============
def get_data():
    edges = []
    genes = {}
    graph = nx.Graph()
    with open('assignment5_input.txt',
              'r') as file:
        for line in file:
            [temp_gene1, temp_gene2] = line.split()
            if temp_gene1 not in genes:
                genes[temp_gene1] = set()
            if temp_gene2 not in genes:
                genes[temp_gene2] = set()

            genes[temp_gene1].add(temp_gene2)
            genes[temp_gene2].add(temp_gene1)
            edges.append(sorted([temp_gene1, temp_gene2]))
            graph.add_edge(temp_gene1, temp_gene2)
    return edges, genes, graph


def get_ground_truth():
    res = []
    with open('complex_merged.txt', 'r') as file:
        # with open(sys.argv[1], 'r') as file:
        for line in file:
            temp_list = line.strip().split(' ')
            tmp = []
            for i in range(len(temp_list)):
                tmp.append(temp_list[i])
            tmp.sort()
            res.append(tmp)
    return res


#######
def set_result_maxi(cliques) :
    res = []
    for clique in cliques:
        for genes in clique:
            res.append(list(genes))
    return res

def selective_joining(graph):
    cliques = [{i, j} for i, j in graph.edges() if i != j]
    k = 2
    while cliques:
        yield k, cliques

        cliques_new = set()
        for u, v in combinations(cliques, 2):
            w = u ^ v
            if len(w) == 2 and graph.has_edge(*w):
                cliques_new.add(tuple(u | v))

        cliques = [set(item) for item in set(frozenset(item) for item in cliques_new)]
        k += 1

def get_maximal_cliques(graph):
    res = []
    cliques_all = []
    for k, cliques in selective_joining(graph):
        cliques_all.append(cliques)
    res = set_result_maxi(cliques_all)
    return res

def delete_duplicate(cliques, graph):
    result = cliques[::-1]
    i = 0
    j = 0
    while True :
        n = len(result)-1
        try:
            clique_last = result[i]
        except:
            break
        result_next = [clique_temp for clique_temp in result if len(set(clique_last) & set(clique_temp)) == 0] # min(len(set(clique_last)), len(clique_temp))
        result_next.insert(i, clique_last)

        if n == len(result_next)-1:
            j += 1
            i += 1
        else:
            i += 1
            result = result_next.copy()

    # flatten : sum(result,[])
    return result

def get_initial_maximal_cliques(graph):
    res = []
    res = get_maximal_cliques(graph)
    res = delete_duplicate(res, graph)
    # with open('maximal_cliques_song.txt',
    #           'r') as file:
    #     for line in file:
    #         res.append(line.split())
    return res


# =========================

# 0. set data =============

def set_result(clusters, score):
    with open(f'./finalproj_output_m2.txt',
              'w') as file:
        for cluster in clusters:
            file.write(f'{len(cluster)} :')
            for gene in cluster:
                file.write(f' {gene} ')
            file.write(f'\n')
        file.write(f'\n')
        file.write(f'{score}')


# =========================


# 1. needed function ===========

def get_jaccard_coefs(edges, genes):
    jaccard_coefs = []
    for u, v in edges:
        nu = genes[u]
        nv = genes[v]

        union = nu | nv
        intersection = nu & nv
        jaccard_coef = len(intersection) / len(union)
        jaccard_coefs.append(jaccard_coef)
    return jaccard_coefs


def get_part(graph, genes):  # get genes from key
    result = {}
    edges = []
    for genes_key in graph:
        neighbors = genes[genes_key]
        result[genes_key] = set()
        for neighbor in neighbors:
            edges.append([genes_key, neighbor])
            result[genes_key].add(neighbor)
    edges = [list(edge) for edge in set(frozenset(edge) for edge in edges)]
    return result

def get_max_jaccrd_coef(jaccard_coefs,edges, genes) :

    # 2. 클러스터 내에 가장 큰 jccrd_cf 값 삭제
    coef_value_max = max(jaccard_coefs)
    coefs_max_indices = [idx for idx, jaccard_coef in enumerate(jaccard_coefs) if
                         jaccard_coef == coef_value_max]

    coefs_max_index = coefs_max_indices[0]
    return edges[coefs_max_index], coefs_max_index
# ================================



def get_cluster_of_gene(gene,clusters) :
    cluster_idx = [idx for idx, cluster in enumerate(clusters) if gene in cluster]
    if len(cluster_idx)  == 0 :
        return -1
    return cluster_idx[0] # cluster에서 제일 큰 것에 속해있는 걸로 반환

#### 중복제거 해야함
def recur1(genes, gene, bfr_jcd_coef, entire_genes, coefs, clusters):

    cluster = [cluster for cluster in clusters if gene in cluster]
    if len(cluster) != 0:
        neighbors = [n for n in genes[gene] if (n not in entire_genes) ] # 이전에 포함된 것들 뺴고 추가해보자
    else:
        neighbors = [n for n in genes[gene] if (n not in entire_genes)]  # 이전에 포함된 것들 뺴고 추가해보자
    if len(neighbors) == 0 :
        return entire_genes, sum(coefs) / len(entire_genes)
    edges_sub = [[n, gene] for n in neighbors]
    jaccard_coefs = get_jaccard_coefs(edges_sub, genes)

    # 이전의 jcd가 같거나 크다면
    if max(jaccard_coefs) <= bfr_jcd_coef:
        return entire_genes, sum(coefs)/len(entire_genes)
    jcd_coef_max = max(jaccard_coefs)
    max_jcd_indice = [idx for idx, jaccard_coef in enumerate(jaccard_coefs) if jaccard_coef == jcd_coef_max]
    neighbors_jcd_max = []
    for max_jcd_idx in max_jcd_indice :
        edges_sub[max_jcd_idx].remove(gene)
        neighbors_jcd_max.append(edges_sub[max_jcd_idx][0])

    val_max = 0
    entire_res = []
    for neighbor_jcd_max in neighbors_jcd_max:
        flag = 0
        entire_genes_sub = entire_genes.copy()
        # 지금까지 묶어온 트리들을 보면, gene-gene_recur-gene_now
        for entire_gene in entire_genes :
            cluster = [cluster for cluster in clusters if entire_gene in cluster]
            if len(cluster) != 0 and neighbor_jcd_max in cluster[0] :
                flag = 1
                continue

        # 만약 어떤 cluster에 속해 있으면 추가로 들어가지 않고 가장 큰 값을 계산하기를 멈춘다.
        if flag == 1 :
            entire_res = entire_genes.copy()
            val_max = sum(coefs)/len(entire_genes)
            continue
        entire_genes_sub.extend(neighbors_jcd_max)
        coefs_sub = coefs.copy()
        coefs_sub.append(max(jaccard_coefs))

        entire_new, val = recur1(genes, neighbor_jcd_max, max(jaccard_coefs), entire_genes_sub, coefs_sub, clusters)
        if val_max <= val :
            val_max = val
            entire_res = entire_new.copy()
    return entire_res, val_max


def merge_clusters(entire_res, clusters):
    cluster_indice = []
    for gene in entire_res:
        cluster_idx_temp = [idx for idx, cluster in enumerate(clusters) if gene in cluster]
        if len(cluster_idx_temp) == 0:
            cluster_indice.append(gene)
            continue
        cluster_idx_temp = cluster_idx_temp[0]
        cluster_indice.extend(clusters[cluster_idx_temp])
        del clusters[cluster_idx_temp]

    cluster_indice = list(set(cluster_indice))
    clusters.append(cluster_indice)
    return clusters


def merging_by_initial_clique(edges, genes, graph, clusters) :
    clustered_genes = list(set(sum(clusters, []))) # 중복 flatten node
    # cluster 내에 모든 gene 중에서 maximal clique에 포함되지 않은 gene들을 선택한다.
    not_clustered_genes = [gene for gene in genes.keys() if gene not in clustered_genes]

    for gene in not_clustered_genes :
        cluster_idx_of_main_gene = get_cluster_of_gene(gene, clusters)
        # 한개짜리가 어딘가 cluster에 속해있는 경우 그냥 넘어간다.
        if cluster_idx_of_main_gene != -1 :
            continue

        neighbors = genes[gene]
        # 한개짜리랑 연결된 모든 edge를 얻어옴
        edges_sub = []
        for n in neighbors :
            edges_sub.append([n, gene])
        # 주변에 연결된 edge의 jaccard 값을 얻어옴
        jaccard_coefs = get_jaccard_coefs(edges_sub, genes)
        # jaccard값 중 가장 큰 index들을 얻어온다.
        max_jcd_indice = [idx for idx,jaccard_coef in enumerate(jaccard_coefs) if jaccard_coef == max(jaccard_coefs)]

        # jaccard 값이 가장 큰 neighbor node들을 얻어온다.
        edges_sub_jcd_max  = [edges_sub[max_jcd_index] for max_jcd_index in max_jcd_indice]
        neighbors_jcd_max = []
        for edge_sub_jcd_max in edges_sub_jcd_max :
            edge_sub_jcd_max.remove(gene)
            neighbors_jcd_max.extend(edge_sub_jcd_max)
        # flatten
        # neighbors_jcd_max = sum(neighbors_jcd_max, [])

        # 만약 jaccard 값이 같은 거가 여러 개가 있다면 jaccard 값 평균이 큰 친구를 클러스터에 추가한다.
        entire_res = []
        val_max = 0
        for neighbor_jcd_max in neighbors_jcd_max :
            entire_new, val = recur1(genes, neighbor_jcd_max, max(jaccard_coefs), [gene, neighbor_jcd_max], [max(jaccard_coefs)], clusters)
            if val_max <= val :
                val_max = val
                entire_res = entire_new.copy()

        # 2. 각각 해당하는 클러스터들이 있다면 다 합쳐서 하나로 만들어버린다.
        clusters = merge_clusters(entire_res, clusters)

    return clusters



def main():
    edges, genes, graph = get_data()  # type(data_input) = nx.graph
    data_grnd_trth = get_ground_truth()
    clusters = get_initial_maximal_cliques(graph)
    clusters = clusters[::-1]
    clusters = merging_by_initial_clique(edges, genes, graph, clusters)
    clusters = sorted(clusters, key=len, reverse=True)

    score = get_F1score(clusters, data_grnd_trth)
    print(f'result : {score}')
    set_result(clusters, score)


if __name__ == '__main__':
    main()
