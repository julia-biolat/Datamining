# 2022년도 2학기 데이터마이닝 - FinalProject_F1score
# 이름 : 이지윤, 정유라, 김송우
# 학번 : 2021247052, 2021247010, 2020253093

import numpy

def get_F1score(data, gr):
    f1 = []
    for i in range(len(data)):
        data_n = set(data[i])
        tmp = []
        for j in range(len(gr)):
            gr_n = set(gr[j])
            tp = data_n.intersection(gr_n)
            if len(tp) > 0:
                # print(i, j)
                tmp.append((2 * len(tp)) / (len(gr_n) + len(data_n)))
        tmp.sort()
        if len(tmp) >= 1:
            f1.append(tmp[len(tmp) - 1])
    mean = numpy.mean(f1)
    return mean

