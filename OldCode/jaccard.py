import pandas as pd
import numpy as np
import matplotlib.pyplot as matpl


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
def jaccard_similarity(data):
    size = data.columns.size
    jaccard_matrix = np.ndarray(shape=(size, size), dtype=float)

    index_1 = 0

    for _, column_A in data.iteritems():
        index_2 = 0

        for __, column_B in data.iteritems():
            m_11 = 0
            m_10 = 0
            m_01 = 0

            for i in range(column_A.values.size):

                if column_A.values[i] == 1 and column_B.values[i] == 1:
                    m_11 += 1

                if column_A.values[i] == 0 and column_B.values[i] == 1:
                    m_01 += 1

                if column_A.values[i] == 1 and column_B.values[i] == 0:
                    m_10 += 1

            index = m_11 / (m_11 + m_01 + m_10)
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

    matpl.imshow(jaccard_sim, cmap='hot', interpolation='none')
    matpl.colorbar(matpl.pcolor(jaccard_sim))
    matpl.title("Jaccard Similarity Index")
    matpl.show()

    matpl.imshow(jaccard_dis, cmap='hot', interpolation='nearest')
    matpl.colorbar(matpl.pcolor(jaccard_dis))
    matpl.title("Jaccard Distance Index")
    matpl.show()

    return

main()















# need to reorganize structure of class directory and code

# make matrix that shows similarity for each pair of nuclear profiles
# size -> 163x163 (NPs that see Hist1)
# number of genomic windows shared

# Jaccard index:
# For each pair of NPs:
#   m_11 = Number of times both saw a given window
#   m_10 = Number of times NP_A saw a window
#   m_01 = Number of time NP_B saw a window
#   Jaccard: m_11 / (m_10 + m_01 + m_11)

# for i in matrix columns:
#     let rest = 163 - i
#     for j in matrix rows:
#         if i[j] == 1:
#             m_10 += 1
#         if i[rest] == 1:
#             m_01 += 1
#         if i[j] == 1 and i[rest] == 1:
#             m_11 += 1
#     j_index = m_11 / (m_10 + m_01 + m_11)
#     similarity_matrix[i][rest] = j_index
