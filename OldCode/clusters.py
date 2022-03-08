import pandas as pd
import numpy as np
import matplotlib.pyplot as matpl
import random
from quickselect import floyd_rivest as fr


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
            # m_10 = 0
            # m_01 = 0

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

            # index = m_11 / (m_11 + m_01 + m_10)
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


# pseudo random choice of nuclear profiles.  Guarantees that if one NP has already
# been selected, it will not be selected again
def choose_rand_profiles(already_chosen, start_index, end_index):
    new_num = random.randint(start_index, end_index)
    while new_num in already_chosen:
        new_num = random.randint(start_index, end_index)
    already_chosen.append(new_num)
    return already_chosen


def assign_clusters(distance_matrix, cluster_list):
    clusters = []
    for i in cluster_list:
        clusters.append([])

    for i in range(len(distance_matrix)):
        distance_list = []
        for j in cluster_list:
            distance_list.append(distance_matrix[i][j[1]])
            # print(i, ": ", distance_matrix[i][j[1]])

        temp = min(distance_list)
        index = [x for x,y in enumerate(distance_list) if y == temp][0]
        # print(index, " ", distance_list)

        flag = True
        for j in range(len(cluster_list)):
            if cluster_list[j][1] == i:
                clusters[j].append(i)
                flag = False

        if flag:
            clusters[index].append(i)
    return clusters



def main():
    indexed, data = extract_csv('./HIST1.csv')
    hist1_data = extract_data('sum', 0, data)

    columns = hist1_data.index.intersection(data.columns)
    hist1_data_frame = data.loc[:, columns]

    # jaccard index
    jaccard_sim = jaccard_similarity(hist1_data_frame)
    # print(jaccard_sim)

    # for i in jaccard:
    #     print(np.isnan(i))

    jaccard_dis = jaccard_distance(jaccard_sim)
    # print(jaccard_dis)

    # matpl.imshow(jaccard_sim, cmap='hot', interpolation='none')
    # matpl.colorbar(matpl.pcolor(jaccard_sim))
    # matpl.title("Jaccard Similarity Index")
    # matpl.show()

    # matpl.imshow(jaccard_dis, cmap='hot', interpolation='nearest')
    # matpl.colorbar(matpl.pcolor(jaccard_dis))
    # matpl.title("Jaccard Distance Index")
    # matpl.show()


    # randomly select 3 NPs to be initial clusters
    rand_nums = []
    for i in range(3):
        choose_rand_profiles(rand_nums, 0, 162)

    # store cluster number with cluster position
    cluster_list = []
    for i in range(len(rand_nums)):
        cluster_list.append((i, rand_nums[i]))
    print(cluster_list)

    # create the clusters
    clusters = assign_clusters(jaccard_dis, cluster_list)
    for i in clusters:
        print(i)

    # improve the clusters
    previous_iteration = clusters
    previous_medoids = cluster_list

    while True:

        distances = []
        for i in clusters:
            distances.append([])

        """
        for every cluster:
          sum all the distances for each pair of elements in cluster
          find minimum sum => new center for that cluster
        """

        for i in range(len(clusters)):
            for j in range(len(clusters[i])):
                dist_sum = 0
                for k in range(len(clusters[i])):
                    dist_sum += jaccard_dis[j][k]
                distances[i].append(dist_sum)

        



        new_centers = []
        for i in range(len(distances)):
            # print("LENGTH: ", len(distances[i]))
            center = fr.nth_largest(distances[i], len(distances[i])//2)
            index = distances[i].index(center)
            new_centers.append((i, index))

        print("new: ", new_centers)
        clusters = assign_clusters(jaccard_dis, new_centers)
                    

        # for sentinel condition:
        #   check percent difference of clusters
        #     sum number of same elements in current and previous iterations,
        #     -> divide by avg len of clusters
        # new sentinel condition
        # lens = [0,0,0]
        # sim = [0,0,0]
        # final = 0
        
        # for i in range(len(clusters)):
        #     lens[i] = (len(clusters[i]) + len(previous_iteration[i])) / 2
            
        # for i in range(len(clusters)):
        #     for j in range(len(clusters[i])):
        #         for k in range(len(previous_iteration)):
        #             for l in range(len(previous_iteration[k])):
        #                 if clusters[i][j] == previous_iteration[k][l]:
        #                     sim[i] += 1
                        

        # for i in range(len(sim)):
        #     sim[i] = sim[i] / lens[i]

        # for i in range(len(sim)):
        #     final += sim[i]
        # final /= len(clusters)

        # print(final)

        # if final >= 90:
        #     break

        # # sentinel condition
        # if clusters == previous_iteration:
        #     break
        # previous_iteration = clusters


    return


main()













        # find the average distance of each element in its cluster from other elements in its cluster
        # """
        # strategy: 
        # for each cluster: 
        #   for each element x:
        #     sum up distance from x to all other elements in cluster, divide by size of cluster => y
        #     store y in list
        # """
        # average_distances_clustered = []
        # for i in clusters:
        #     average_distances_clustered.append([])

        # for i in range(len(clusters)):
        #     dist_sum = 0
        #     for j in range(len(clusters[i])):
        #         dist_sum += jaccard_dis[j][i]
        #     dist_sum /= len(clusters[i])
        #     average_distances_clustered[i].append(dist_sum)
                

        
        # find the median of each cluster
        # begin by clustering the distances
        # average_distances_clustered = []
        # for i in clusters:
        #     average_distances_clustered.append([])

        # for i in range(len(clusters)):
        #     for j in clusters[i]:
        #         average_distances_clustered[i].append(average_distances[j])
