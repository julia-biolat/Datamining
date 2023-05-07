# 2022년도 2학기 데이터마이닝 - FinalProject_F1Score
# 이름 : 이지윤, 정유라, 김송우
# 학번 : 2021247052, 2021247010, 2020253093

import sys
import time
from itertools import combinations
import networkx as nx
import numpy as np
from F1score import get_F1score


# 0. get data ============
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

def get_assignment5(str):
    res = []
    with open(str, 'r') as file:  # 중복 제거 안 한 거
        for line in file:
            temp_list = line.strip().split()
            tmp = []
            for i in range(len(temp_list)):
                tmp.append(temp_list[i])
            tmp.sort()
            res.append(tmp)
    return res

def get_data(str):
    res = []
    with open(str, 'r') as file:
        for line in file:
            temp_list = line.strip().split()
            tmp = []
            for i in range(2,len(temp_list)):
                tmp.append(temp_list[i])
            tmp.sort()
            res.append(tmp)
    return res
# =========================

# 0. set data =============


# =========================


def main():
    # edges, genes, graph = get_data()  # type(data_input) = nx.graph
    data_grnd_trth = get_ground_truth()
    data_assignment5 = get_assignment5('assignment5_output.txt')
    data_assignment6 = get_data('assignment6_output.txt')
    data_final2_1 = get_data('finalproj_output_m1.txt')
    data_final2_2 = get_data('finalproj_output_m2.txt')
    data_final2_3 = get_data('finalproj_output_m3.txt')

    score_assignmet5 = get_F1score(data_assignment5, data_grnd_trth)
    score_assignmet6 = get_F1score(data_assignment6, data_grnd_trth)
    score_finalproj2_1 = get_F1score(data_final2_1, data_grnd_trth)
    score_finalproj2_2 = get_F1score(data_final2_2, data_grnd_trth)
    score_finalproj2_3 = get_F1score(data_final2_3, data_grnd_trth)

    print(f'assignment5(중복제거 안한 거) : {score_assignmet5}')
    print(f'assignment6 : {score_assignmet6}')
    print(f'finalproj m1 : {score_finalproj2_1}')
    print(f'finalproj m2 : {score_finalproj2_2}')
    print(f'finalproj m3 : {score_finalproj2_3}')

if __name__ == '__main__':
    main()

