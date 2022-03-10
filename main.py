import cluster as cl
import pandas as pd
import numpy as np
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


def counter_cosine_similarity(c1, c2):
    terms = set(c1).union(c2)
    dotprod = sum(c1.get(k, 0) * c2.get(k, 0) for k in terms)
    magA = math.sqrt(sum(c1.get(k, 0)**2 for k in terms))
    magB = math.sqrt(sum(c2.get(k, 0)**2 for k in terms))
    return dotprod / (magA * magB)

def main():
    indexed, data = extract_csv('./CSVs/HIST1.csv')
    hist1_data = extract_data('sum', 0, data)
    hist1_data_indexed = extract_data('sum', 0, indexed)

    columns = hist1_data.index.intersection(data.columns)
    columns_indexed = hist1_data_indexed.index.intersection(indexed.columns)
    hist1_data_frame = data.loc[:, columns]

    # jaccard index
    jaccard_sim = jaccard_similarity(hist1_data_frame)
    # print(jaccard_sim)

    jaccard_dis = jaccard_distance(jaccard_sim)
    # print(jaccard_dis)

    global_sim_list = []
    global_cluster_list = []
    global_center_list = []

    for outer_index in range(5):
        cluster = cl.Clusters(3, jaccard_dis)

        cluster.set_random_centers()
        cluster_list = cluster.assign_clusters()
        cluster.set_cluster_list(cluster_list)

        previous_iteration = cluster_list
        while True:
            cluster.set_new_centers()
            cluster_list = cluster.assign_clusters()
            cluster.set_cluster_list(cluster_list)

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
                    # print("the bastard: ", i, " ", sim_list_2[i])
                    sentinel_list[i] = False
                
            if Counter(sentinel_list)[True] >= 2:
                global_sim_list.append(sim_list_2)
                global_cluster_list.append(cluster_list)
                # global_center_list.append(new_centers)
                print("DONE")
                break
            else:
                previous_iteration = cluster_list
                
        print(global_cluster_list)
            


main()        
