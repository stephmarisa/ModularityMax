# parameters:
#     p_G = networkx graph
#     p_k = number of clusters to be found
# returns: A list of k sets, where each set contains the nodes in one cluster
# clusters are determined via modularity-maximization
def modularity_cluster(p_G, p_k):
    
    # create container for sets representing clusters
    list_sets = []
    
    # initially add one set of all nodes to list_sets
    list_sets.append(set(p_G.nodes()))
    
    # form 2^x clusters, where 2^x is the largest power of 2 such that 2^x <= k
    for i in range(int(math.log(p_k, 2))):
        # repeatedly bisect all sets of nodes in list_sets
        list_sets = modularity_bisect_list(p_G, list_sets)
    
    # if k - 2^x > 1, then form additional clusters one at a time:
    if p_k - pow(2, int(math.log(p_k, 2))) > 0:
        # create a container to hold possible partitions
        list_lists = []
        # for as many additional clusters as needed:
        for j in range(p_k - pow(2, int(math.log(p_k, 2)))):
            # for each cluster in the original 2^x clusters
            for a_set in list_sets:
                a_list = list_sets.copy()
                # remove this cluster
                a_list.remove(a_set)
                # bisect the cluster that was removed
                (set_1, set_2) = modularity_bisect(p_G.subgraph(a_set))
                # add the two clusters formed by the bisection back into list_sets
                a_list.append(set_1)
                a_list.append(set_2)
                # add this possible partition to the list of possible partitions
                list_lists.append(a_list)
        # create a container to hold the modularity scores of each possible partition
        mod_scores = []
        # calculate the modularity score of each possible partition
        for a_list in list_lists:
            mod_scores.append(nx.algorithms.community.quality.performance(p_G, a_list))
        # return the partition corresponding to the highest modularity score
        return list_lists[mod_scores.index(max(mod_scores))]
    
    # return the list of sets corresponding to clusters (a.k.a. return the computed partition)
    return list_sets


# purpose: use modularity-maximization to bisect a subgraph
# parameters: p_G = networkx subgraph
# returns: A tuple of sets, where each set contains the nodes in one cluster of the bisected subgraph
def modularity_bisect(p_G):
    
    # compute the modularity maxtrix of the subgraph
    mod_matrix = nx.modularity_matrix(p_G)
    
    # compute the leading eigenvector of the modularity matrix
    eig = np.linalg.eigh(mod_matrix)
    leading_eigval_index = eig[0].argmax()
    leading_eigvect = eig[1][:, leading_eigval_index]
    
    # assign each node to one of two sets, based on the sign of the corresponding element in the leading eigenvector
    nodes_list = list(p_G.nodes())
    set_1 = set()
    set_2 = set()
    for i in range(len(leading_eigvect)):
        if leading_eigvect[i] < 0:
            set_1.add(nodes_list[i])
        else:
            set_2.add(nodes_list[i])
            
    # return these two sets of nodes
    return (set_1, set_2)


# purpose: use modularity-maximization to bisect a list of subgraphs
# parameters: 
#     p_G = original full networkx graph
#     p_list_sets = list of sets of nodes
# returns: A list of sets that is twice as long as p_list_sets; contains results of bisecting each set in p_list_sets
def modularity_bisect_list(p_G, p_list_sets):
    # create a container to hold the bisections of each set
    new_list_sets = []
    # bisect each set
    for a_set in p_list_sets:
        a_subgraph = p_G.subgraph(a_set)
        (set_1, set_2) = modularity_bisect(a_subgraph)
        new_list_sets.append(set_1)
        new_list_sets.append(set_2)
    return new_list_sets