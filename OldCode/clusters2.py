import pandas as pd
import numpy as np
import matplotlib.pyplot as matpl
import random
from quickselect import floyd_rivest as fr
from collections import Counter
import math

# read the csv file and remove undesired columns
def extract_csv(file_name):
    data_indexed = pd.read_csv(file_name)
    if 'Unnamed' in data_indexed.columns[0]:
        data = data_indexed.iloc[:, 4:]
    else:
        data = data_indexed.iloc[:, 3:]
    return data_indexed, data

# extract the NPs of interest to us at the moment
# function is aggregate function to perform,
# result is conditional used to drop columns
# data_frame is the pandas data frame to operate on
def extract_data(function, result, data_frame):
    res = data_frame.agg(function)
    for i in res.index:
        if res[i] == result:
            res = res.drop(i)
    return res

# input is a pandas data frame
# modified: denominator of jaccard index is:
# min{A, B}, where A = m_11+m_10 for column A and B is same for column B
def jaccard_similarity(data):
    size = data.columns.size
    jaccard_matrix = np.ndarray(shape=(size, size), dtype=float)

    index_1 = 0

    for _, column_A in data.iteritems():
        index_2 = 0
        for __, column_B in data.iteritems():
            A = 0
            B = 0
            m_11 = 0 

            for i in range(column_A.values.size):

                if column_A.values[i] == 1 and column_B.values[i] == 1:
                    m_11 += 1
                    A += 1
                    B += 1

                if column_A.values[i] == 0 and column_B.values[i] == 1:
                    # m_01 += 1
                    B += 1

                if column_A.values[i] == 1 and column_B.values[i] == 0:
                    A += 1
                    # m_10 += 1

            index = m_11 / min([A, B])
            jaccard_matrix[index_1][index_2] = index
            index_2 += 1
            
        index_1 += 1
    return jaccard_matrix


def jaccard_distance(j_matrix):
    # go through jaccard similarity and subtract each value from 1.
    # save that value in the new matrix.  return the matrix
    size = j_matrix.shape[0]
    jaccard_dis = np.ndarray(shape=(size, size), dtype=float)

    for i in range(size):
        for j in range(size):
            jaccard_dis[i][j] = 1 - j_matrix[i][j]

    return jaccard_dis


def assign_clusters(center_list, dists):
    # print("center list assignment: ", center_list)

    cluster_list = []
    for i in center_list:
        cluster_list.append([])

    for i in range(len(dists)):
        temp_list = [0,0,0]
        for j in range(len(center_list)):
            if i == center_list[j][0]:
                temp_list[j] = -100
                # print("guaranteeing center: ", i, " at ", j)
            else:
                temp_list[j] = dists[i][center_list[j][0]]
        index = np.argmin(temp_list)
        cluster_list[index].append(i)
    return cluster_list


def counter_cosine_similarity(c1, c2):
    terms = set(c1).union(c2)
    dotprod = sum(c1.get(k, 0) * c2.get(k, 0) for k in terms)
    magA = math.sqrt(sum(c1.get(k, 0)**2 for k in terms))
    magB = math.sqrt(sum(c2.get(k, 0)**2 for k in terms))
    return dotprod / (magA * magB)


def main():
    indexed, data = extract_csv('./HIST1.csv')
    hist1_data = extract_data('sum', 0, data)

    columns = hist1_data.index.intersection(data.columns)
    hist1_data_frame = data.loc[:, columns]

    # jaccard index
    jaccard_sim = jaccard_similarity(hist1_data_frame)
    # print(jaccard_sim)

    jaccard_dis = jaccard_distance(jaccard_sim)
    # print(jaccard_dis)

    # get random numbers as centers for the first iteration
    num_clusters = 3
    rand_nums = np.ndarray(shape=(num_clusters,1), dtype=int)
    for i in range(len(rand_nums)):
        rand_nums[i][0] = 0
    for i in range(len(rand_nums)):
        center_i = random.randint(0, 162)
        while center_i in rand_nums[:,0]:
            center_i = random.randint(0, 162)
        rand_nums[i][0] = center_i

    # start making the clusters
    cluster_list = assign_clusters(rand_nums, jaccard_dis)
    print(cluster_list)
    print("==============================================================================")

    oldest_iteration = cluster_list
    previous_iteration = cluster_list

    while True:
        # finding the new centers to get better clusters
        new_centers_wrt_clusters = np.ndarray(shape=(num_clusters,1), dtype=int)
        for i in range(len(new_centers_wrt_clusters)):
            new_centers_wrt_clusters[i][0] = 0

        # finding average dissimilarity within clusters
        distance_list = []
        for i in cluster_list:
            distance_list.append([])
        for k in range(len(cluster_list)):
            for i in range(len(cluster_list[k])):
                dis_sum = 0
                for j in range(len(cluster_list[k])):
                    dis_sum += jaccard_dis[i][j]
                distance_list[k].append(dis_sum)

        # finding the new centers and assigning them to a matrix
        for i in range(len(distance_list)):
            median = fr.nth_largest(distance_list[i], len(distance_list[i])//2)
            index = [x for x,e in enumerate(distance_list[i]) if e == median][0]
            new_centers_wrt_clusters[i][0] = index

        # print("new centers wrt clusters: ", new_centers_wrt_clusters)

        new_centers = np.ndarray(shape=(num_clusters, 1), dtype=int)
        for i in range(len(new_centers)):
            new_centers[i][0] = 0

        for i in range(len(new_centers)):
            index = new_centers_wrt_clusters[i][0]
            new_centers[i][0] = cluster_list[i][index]

        print("new centers wrt population: ", new_centers)

        cluster_list = assign_clusters(new_centers, jaccard_dis)
        # print(cluster_list)
        print("==============================================================================")

        # sentinel condition
        sim_list = []
        sim_list_2 = []
        for i in cluster_list:
            sim_list.append([])
            # sim_list_2.append([])

        for i in range(len(cluster_list)):
            sim_list[i].append(Counter(cluster_list[i]))
            sim_list[i].append(Counter(previous_iteration[i]))
        
        for i in range(len(sim_list)):
            sim_list_2.append(counter_cosine_similarity(sim_list[i][0], sim_list[i][1]))

        sentinel_list = []
        for i in cluster_list:
            sentinel_list.append(True)
            
        for i in range(len(sim_list_2)):
            if sim_list_2[i] <= 0.7:
                print("the bastard: ", i, " ", sim_list_2[i])
                sentinel_list[i] = False

        if Counter(sentinel_list)[True] >= 2:
            print("DONE")
            break
        else:
            previous_iteration = cluster_list
        

        # if cluster_list == previous_iteration or cluster_list == oldest_iteration:
        #     print("DONE")
        #     break
        # else:
        #     oldest_iteration = previous_iteration
        #     previous_iteration = cluster_list

    return
main()

