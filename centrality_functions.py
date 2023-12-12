def get_betweenness_centrality_log10(model_graph_df):
    # calculate betweenness centrality
    model_betweenness_centrality = nx.betweenness_centrality(model_graph_df, weight='weight', normalized=False)
    
    # strip whitespace from keys and log-transform the values
    log10_betweenness_centrality = {key.strip(): np.log10(value + 1) for key, value in model_betweenness_centrality.items()}
   
    return log10_betweenness_centrality


def get_bridging_centrality_log10(model_graph_df):
    # calculate betweenness centrality
    betweenness = nx.betweenness_centrality(model_graph_df, weight='weight', normalized=False)

    # calculate the bridging coefficient for each node
    bridging_coefficient = {}
    for node in model_graph_df.nodes():
        
        # for directed graphs, consider successors and predecessors
        if model_graph_df.is_directed():
            successors = list(model_graph_df.successors(node))
            predecessors = list(model_graph_df.predecessors(node))
            degree_sum = sum(model_graph_df.out_degree(successor) for successor in successors) + \
                         sum(model_graph_df.in_degree(predecessor) for predecessor in predecessors)
        else:
            # for undirected graphs, consider neighbors
            neighbors = list(model_graph_df.neighbors(node))
            degree_sum = sum(model_graph_df.degree(neighbor) for neighbor in neighbors)

        bridging_coefficient[node] = 1 / degree_sum if degree_sum > 0 else 0

    # calculate bridging centrality
    bridging_centrality = {node: betweenness[node] * bridging_coefficient[node] for node in model_graph_df.nodes()}

    # strip whitespace from keys and log-transform the values
    log10_bridging_centrality = {key.strip(): np.log10(value + 1) for key, value in bridging_centrality.items()}

    return log10_bridging_centrality

def get_clustering_coefficient(model_graph_df):
    model_clustering_coefficient = nx.clustering(model_graph_df)

    # strip whitespace from keys and log-transform the values
    clustering_coefficient = {key.strip(): value for key, value in model_clustering_coefficient.items()}

    return clustering_coefficient

def get_degree(model_graph_df, log=False, directional=False):
    model_degree = dict(model_graph_df.degree())
    model_degree = {key.strip(): value for key, value in model_degree.items()}
    
    if (directional):
        c=0
    
    # strip whitespace from keys and log-transform the values
    if (log):
        model_degree = {key.strip(): np.log10(value + 1) for key, value in model_degree.items()}
    
    return model_degree