import networkx as nx
import sys
import csv
from tqdm import tqdm
from os import listdir
from os.path import isfile, join
import numpy as np
from sklearn.preprocessing import normalize

if __name__ == "__main__":
    arg_len = len(sys.argv)
    central_genes = []

    if arg_len < 4:
        print("ERROR: not enough arguments. Specify itneraction file, centrality measures folder and type of analysis [count, rank, score]")
        sys.exit()

    if sys.argv[3] != "count" and sys.argv[3] != "rank" and sys.argv[3] != "score":
        print("Select and aggregation method from [count, rank, score]")
        sys.exit()

    mode = sys.argv[3]

    if mode == "count":
        LIST_FILE = sys.argv[2]
        
        fileList = [f for f in listdir(LIST_FILE) if isfile(join(LIST_FILE, f))]
        print(fileList)
        edges = []

        for fi in fileList:


            with open(LIST_FILE + "/" + fi, "r") as dataFile:
                g_list = []
                dataFile.readline()
                dataFile.readline()
                dataFile.readline()
                for row in dataFile:
                    
                    line = row.strip().split("\t")
                    
                    g_list.append(line[0])
            
            central_genes.append(g_list)

        #print(central_genes)
        DATA_FILE = sys.argv[1]

        edges = []

        with open(DATA_FILE, "r") as dataFile:
            for row in dataFile:
                line = row.strip().split("\t")
                edges.append((line[2], line[3]))
        
        netName = DATA_FILE
        G = nx.Graph(name = DATA_FILE)
        G.add_edges_from(edges)

        print(nx.info(G))

        all_genes = G.nodes()
        top_genes = {}

        for gene in all_genes:
            count = 0
            for centrality in central_genes:
                if gene in centrality:
                    count += 1
            top_genes[gene] = count
        
        avg_count = sum(list(top_genes.values()))/len(top_genes)
        #print(top_genes)
        top_gene_names = list(top_genes.keys())
        top_genes_scores = list(top_genes.values())

        top_gene_names_sorted = [x for _,x in sorted(zip(top_genes_scores,top_gene_names), reverse=True)]

        with open(DATA_FILE[:-4] + "_top_genes.txt", "w+") as saveFile:
            saveFile.write("Avg score: " + str(avg_count) + "\n\n")
            saveFile.write("Gene\tScore\n")
            for gene in top_gene_names_sorted:
                if top_genes[gene] > avg_count:
                    saveFile.write(gene +"\t" + str(top_genes[gene]) + "\n")

    elif mode == "score":
        LIST_FILE = sys.argv[2]
        
        genes_scores = {}
        fileList = [f for f in listdir(LIST_FILE) if isfile(join(LIST_FILE, f))]
        print(fileList)
        edges = []

        for fi in fileList:

            g_list = []
            v_list = []

            with open(LIST_FILE + "/" + fi, "r") as dataFile:
                g_list = []
                dataFile.readline()
                dataFile.readline()
                dataFile.readline()
                for row in dataFile:
                    
                    line = row.strip().split("\t")
                    g_list.append(line[0])
                    v_list.append(float(line[1]))
            
            v_list = (v_list / np.linalg.norm(v_list)).tolist()
            #v_list = normalize(np.array(v_list).reshape(-1, 1)).tolist()

            for gene in g_list:
                if gene not in genes_scores:
                    score_index = g_list.index(gene)
                    genes_scores[gene] = [v_list[score_index]]
                else:
                    genes_scores[gene].append(v_list[score_index])
                

        #print(central_genes)
        DATA_FILE = sys.argv[1]

        edges = []

        with open(DATA_FILE, "r") as dataFile:
            for row in dataFile:
                line = row.strip().split("\t")
                edges.append((line[2], line[3]))
        
        netName = DATA_FILE
        G = nx.Graph(name = DATA_FILE)
        G.add_edges_from(edges)

        print(nx.info(G))

        for gene in genes_scores:
            genes_scores[gene] = sum(genes_scores[gene]) / len(genes_scores[gene])
        
        avg_score = sum(list(genes_scores.values()))/len(genes_scores)
        #print(top_genes)
        top_gene_names = list(genes_scores.keys())
        top_genes_scores = list(genes_scores.values())

        top_gene_names_sorted = [x for _,x in sorted(zip(top_genes_scores,top_gene_names), reverse=True)]

        with open(DATA_FILE[:-4] + "_top_genes_normalized_scores.txt", "w+") as saveFile:
            saveFile.write("Avg score: " + str(avg_score) + "\n\n")
            saveFile.write("Gene\tScore\n")
            for gene in top_gene_names_sorted:
                if genes_scores[gene] > avg_score:
                    saveFile.write(gene +"\t" + str(genes_scores[gene]) + "\n")


        #print(nx.info(G))