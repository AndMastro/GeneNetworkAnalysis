import networkx as nx
import sys
import csv
from tqdm import tqdm

if __name__ == "__main__":
    

    arg_len = len(sys.argv)

    if arg_len < 4:
        print("ERROR: not enough arguments. Specify file, centrality measure and [top, all].")
        sys.exit()

    if sys.argv[3] != "top" and sys.argv[3] != "all":
        print("ERROR: choose \"top\" to get only the central genes or \"all\" to get all the genes in the network.")
        sys.exit()

    
    DATA_FILE = sys.argv[1]
    centrality = sys.argv[2]

    top = sys.argv[3] == "top"

    edges = []

    with open(DATA_FILE, "r") as dataFile:
        for row in dataFile:
            line = row.strip().split("\t")
            edges.append((line[2], line[3]))
    
    netName = DATA_FILE
    G = nx.Graph(name = DATA_FILE)
    G.add_edges_from(edges)

    print(nx.info(G))

    if centrality == "degree":

        avg_degree = sum(dict(G.degree()).values())/G.number_of_nodes()
        print("Average " + centrality + ": " + str(avg_degree))

        deg = G.degree()

        deg_dict = {}
        for gene in deg:
            deg_dict[gene[0]] = gene[1]

        # Create graph usong only central genes (deg > avg_degree)


        central_genes = []
        central_measures = []

        if top:
            for pair in deg:
                if pair[1] > avg_degree:
                    central_genes.append(pair[0])
                    central_measures.append(pair[1])
        else:
            for pair in deg:
                central_genes.append(pair[0])
                central_measures.append(pair[1])

        central_genes = [x for _,x in sorted(zip(central_measures,central_genes), reverse=True)]
        central_edges = []

        for gene1, gene2 in edges:
            if gene1 in central_genes and gene2 in central_genes:
                central_edges.append((gene1, gene2))

        central_G = nx.Graph(name = "Central Genes Network")

        central_G.add_edges_from(central_edges)

        print(nx.info(central_G))
        
        SAVE_PATH = None
        if top:
            SAVE_PATH = DATA_FILE[:-4] + "_" + centrality +"_CENTRAL_GENES.txt"
        else:
            SAVE_PATH = DATA_FILE[:-4] + "_" + centrality + ".txt"

        with open(SAVE_PATH, "w+") as saveFile:
            saveFile.write("Avg " + centrality +": " + str(avg_degree) + "\n\n")
            saveFile.write("Gene\tDegree\n")
            for gene in central_genes:
                saveFile.write(gene + "\t" + str(deg_dict[gene]) + "\n")
        

        '''deg_list = sorted(list(G.degree()), reverse = True, key = lambda x:x[1])

        central_deg_list = []

        for gene in deg_list:
            if gene[0] in genes_list:
                central_deg_list.append(gene)
                
        print(sorted(central_deg_list, reverse = True, key = lambda x:x[1]))'''
        
    else:
        centrality_measure = None
        avg_centrality = None

        if centrality == "betweenness":

            centrality_measure = nx.betweenness_centrality(G)
            avg_centrality = sum(centrality_measure.values())/G.number_of_nodes()

        elif centrality == "closeness":

            centrality_measure = nx.closeness_centrality(G)
            avg_centrality = sum(centrality_measure.values())/G.number_of_nodes()

        elif centrality == "eigenvector":
            centrality_measure = nx.eigenvector_centrality(G)
            avg_centrality = sum(centrality_measure.values())/G.number_of_nodes()

        elif centrality == "eccentricity":
            centrality_measure = nx.eccentricity(G)
            avg_centrality = sum(centrality_measure.values())/G.number_of_nodes()

        elif centrality == "diameter":
            diameter = nx.diameter(G)
            print("\nNetwork diameter (max eccentricity): " + str(diameter))
            sys.exit()

        elif centrality == "average_distance":
            distance = nx.average_shortest_path_length(G)
            print("\nNetwork average distance (avg shortest path): " + str(distance))
            sys.exit()
        
        elif centrality == "radiality": #slighly differnet values from Centiscape (raking is the same). Difference is probably due to apporx. Check on documentation later
            diameter = nx.diameter(G)
            nodes = G.nodes()
            numNodes = G.number_of_nodes()
            centrality_measure = {}

            for gene in nodes:
                paths = nx.shortest_path_length(G, source = gene)
                for target in paths:
                    paths[target] = (diameter + 1) - paths[target]
                centrality_measure[gene] = sum(list(paths.values())) / (numNodes-1)
            
            avg_centrality = sum(centrality_measure.values())/numNodes
        
        elif centrality == "stress": # DIFFERNET VALUES FROM CENTISCAPE! WHY? TO CHECK. First ranked genes are the same, tho.

            all_paths = dict(nx.all_pairs_shortest_path(G))
            nodes = G.nodes()
            numNodes = G.number_of_nodes()
            centrality_measure = {}

            for gene in nodes:
                centrality_measure[gene] = 0
                for gene1 in all_paths:
                    for gene2 in all_paths:
                        if gene in all_paths[gene1][gene2] and gene != gene1 and gene != gene2:
                            centrality_measure[gene] += 1

            avg_centrality = sum(centrality_measure.values())/numNodes

        elif centrality == "centroid_value":
            
            centrality_measure = {}

            nodes = G.nodes()

            for gene1 in tqdm(nodes):
                gene1_paths = nx.shortest_path_length(G, source=gene1)
                gene1_closer = []
                for gene2 in nodes:
                    if gene1 != gene2:
                        gene2_paths = nx.shortest_path_length(G, source=gene2)
                        closer_1 = 0
                        closer_2 = 0
                        for gene_target in nodes:
                            if gene_target != gene1 and gene_target != gene2:
                                if gene_target not in gene1_paths and gene_target in gene2_paths:
                                    closer_2 += 1
                                elif gene_target not in gene2_paths and gene_target in gene1_paths:
                                    closer_1 += 1
                                elif gene_target not in gene2_paths and gene_target not in gene1_paths:
                                    pass
                                elif gene1_paths[gene_target] < gene2_paths[gene_target]:
                                    closer_1 += 1
                                elif gene1_paths[gene_target] > gene2_paths[gene_target]:
                                    closer_2 += 1
                                else:
                                    pass
                        gene1_closer.append(closer_1-closer_2)
                centrality_measure[gene1] = min(gene1_closer)

            avg_centrality = sum(centrality_measure.values())/G.number_of_nodes()

            

        elif centrality == "bridging": #VALUES AND RAKING COMPLETELY DIFFERENT FROM CENTISCAPE.
            
            centrality_measure = {}

            btw = nx.betweenness_centrality(G)

            BC = {}

            nodes = G.nodes()

            for gene in nodes:
                deg_gene = G.degree(gene)
                neighbours = G.neighbors(gene) #it returns an iterator! Once iterated, your are done (can your reset it?).
                #print(list(neighbors))
                deg_neig = []
                sigmas = []
                neighbours_list = []
                for n in neighbours:
                    neighbours_list.append(n)
                    deg_neig.append(G.degree(n))

                #print(neighbours_list)
                
                for n in neighbours_list:
                    edges_outgoing_neigborhood = 0
                    neigh_n = G.neighbors(n)
                    
                    for n_n in neigh_n:
                        if n_n != gene and n_n not in neighbours_list:
                            edges_outgoing_neigborhood +=1 #we count the nodes, not edges but it is the same computation
                    sigmas.append(edges_outgoing_neigborhood)

                #print(sigmas)
                   
                SIGMA = []
                for sig, d in zip(sigmas, deg_neig):
                    if d == 1:
                        SIGMA.append(0)
                    else:
                        SIGMA.append(sig/(d-1))

                BC[gene] = (1/deg_gene) * sum(SIGMA)

                #sys.exit()

            for gene in nodes:
                centrality_measure[gene] = btw[gene]*BC[gene]
            #print(centrality_measure["KDF1"])


            avg_centrality = sum(centrality_measure.values())/G.number_of_nodes()
            

        elif centrality == "edge_betweenness": #values are different from Centiscape but the ranking is the same

            centrality_measure = nx.edge_betweenness_centrality(G, normalized = False)
            avg_centrality = sum(centrality_measure.values())/G.number_of_edges()

            print("Average " + centrality + ": " + str(avg_centrality))

            central_edges = []
            central_measures = []

            if top:
                for edge in centrality_measure:
                    if centrality_measure[edge] > avg_centrality:
                        central_edges.append(edge)
                        central_measures.append(centrality_measure[edge])
            else:
                for edge in centrality_measure:
                    central_edges.append(edge)
                    central_measures.append(centrality_measure[edge])
            
            central_edges = [x for _,x in sorted(zip(central_measures,central_edges), reverse=True)]
            
            central_G = nx.Graph(name = "Central Genes Network")

            central_G.add_edges_from(central_edges)

            print(nx.info(central_G))

            central_genes = central_G.nodes()
            
            SAVE_PATH = None
            if top:
                SAVE_PATH = DATA_FILE[:-4] + "_" + centrality +"_CENTRAL_GENES.txt"
            else:
                SAVE_PATH = DATA_FILE[:-4] + "_" + centrality +".txt"
        

            with open(SAVE_PATH, "w+") as saveFile:
                saveFile.write("Avg " + centrality +": " + str(avg_centrality) + "\n\n")
                saveFile.write("Gene"+ "\n")
                for gene in central_genes:
                    saveFile.write(gene + "\n")

                saveFile.write("\nGene1\t"+ "Gene2\t" + centrality +"\n")
                for edge in central_edges:
                    saveFile.write(edge[0] + "\t" + edge[1] + "\t" + str(centrality_measure[edge]) + "\n")

            sys.exit()
            
        else:
            print("Centrality measure not recognized: choose from [degree, betweenness, closeness, eigenvector, eccentricity, diameter, average_distance, radiality, stress, centroid_value, bridging, edge_betweenness]")
            sys.exit()

        print("Average " + centrality + ": " + str(avg_centrality))

        # Create graph usong only central genes (deg > avg_degree)


        central_genes = []
        central_measures = []

        if top:
            for gene in centrality_measure:
                if centrality_measure[gene] > avg_centrality:
                    central_genes.append(gene)
                    central_measures.append(centrality_measure[gene])
        else:
            for gene in centrality_measure:
                central_genes.append(gene)
                central_measures.append(centrality_measure[gene])

        central_genes = [x for _,x in sorted(zip(central_measures,central_genes), reverse=True)]
        central_edges = []
        
        for gene1, gene2 in edges:
            if gene1 in central_genes and gene2 in central_genes:
                central_edges.append((gene1, gene2))
        central_G = nx.Graph(name = "Central Genes Network")

        central_G.add_nodes_from(central_genes)
        central_G.add_edges_from(central_edges)

        print(nx.info(central_G))

        central_G.nodes()

        SAVE_PATH = None
        if top:
            SAVE_PATH = DATA_FILE[:-4] + "_" + centrality +"_CENTRAL_GENES.txt"
        else:
            SAVE_PATH = DATA_FILE[:-4] + "_" + centrality +".txt"
        

        with open(SAVE_PATH, "w+") as saveFile:
            saveFile.write("Avg " + centrality +": " + str(avg_centrality) + "\n\n")
            saveFile.write("Gene\t"+ centrality +"\n")
            for gene in central_genes:
                saveFile.write(gene + "\t" + str(centrality_measure[gene]) + "\n")

        
    # To add: take in input interaction file and centrality measure to use. Then, return a file with central genes: GENE_NAME DEG_ORIGINAL DEG_CENTRAL_NETWORK + avg degree in both networks (may be useful). To that for now, then, decide if consider edge weight.
