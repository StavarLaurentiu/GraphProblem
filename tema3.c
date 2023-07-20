#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_NAME 51

//// Structures List ////

typedef struct l_node {
    int num, cost;
    struct l_node *next;
} l_node;

typedef struct List {
    int dimension;
    struct l_node *head;
} List;

//// Structures List ////

//// Structures Graf ////

typedef struct g_node {
    int num, depth;
    char name[MAX_NAME];
} g_node;

typedef struct Graf {
    int dimension;
    struct g_node **nodes;
    struct List **link;
} Graf;

//// Structures Graf ////

List *init_list();
int is_empty_list(List *list);
void add_list(List *list, int num, int cost);
void del_list(List *list, int num);
void free_list(List *list);

int get_index(g_node **nodes, char *name, int num_nodes);
Graf *init_graph(int num_nodes);
void dfs(Graf *graf, int current_node, int *visited);
int connex_components(Graf *graf, int num_nodes);
int min_index(int *v, int n, int *visited);
int min_index_float(float *v, int n, int *visited);
int prims(Graf *graf, int *visited, int num_nodes);
int cmpfunction(const void *a, const void *b);
void dijktra(Graf *graf, int source_index, int dest_index, int *prev);
void print_road(FILE *fid, Graf *graf, int source_index, int dest_index,
                int *prev);
int road_cost(Graf *graf, int source_index, int dest_index, int *prev);
int my_min(int a, int b);
int weight(Graf *graf, int source_index, int dest_index, int *prev);
void free_graph(Graf *graf);

//// List Metods ////

List *init_list() {
    List *list = (List *)malloc(1 * sizeof(List));
    list->head = NULL;
    list->dimension = 0;

    return list;
}

int is_empty_list(List *list) { return list->dimension == 0; }

void add_list(List *list, int num, int cost) {
    l_node *new_node = (l_node *)malloc(1 * sizeof(l_node));
    new_node->num = num;
    new_node->cost = cost;
    new_node->next = list->head;
    list->head = new_node;

    (list->dimension)++;
}

void del_list(List *list, int num) {
    if (is_empty_list(list)) return;

    l_node *current_node = list->head;
    l_node *prev_node = NULL;
    while (current_node != NULL) {
        if (current_node->num == num) break;
        prev_node = current_node;
        current_node = current_node->next;
    }

    // The node which has to be erased is not in the list
    if (current_node == NULL) return;

    // If the list has only one element
    if (prev_node == NULL) {
        list->head = list->head->next;
    } else {
        prev_node->next = current_node->next;
    }

    // Free the memory for the erased node
    free(current_node);

    (list->dimension)--;
}

void free_list(List *list) {
    l_node *current_node = list->head;

    while (current_node != NULL) {
        l_node *to_free = current_node;
        current_node = current_node->next;
        free(to_free);
    }

    free(list);
}

//// List Metods ////

// Return the index of a node given by it's name
int get_index(g_node **nodes, char *name, int num_nodes) {
    int i;
    for (i = 1; i <= num_nodes; i++) {
        if (strcmp(nodes[i]->name, name) == 0) return i;

        if (strcmp(nodes[i]->name, "-") == 0) {
            strcpy(nodes[i]->name, name);
            return i;
        }
    }

    return -1;
}

Graf *init_graph(int num_nodes) {
    Graf *graf = (Graf *)malloc(1 * sizeof(Graf));
    graf->dimension = num_nodes;

    graf->link = (List **)malloc((num_nodes + 1) * sizeof(List *));

    int i;
    for (i = 1; i <= num_nodes; i++) {
        graf->link[i] = init_list();
    }

    // Declare a vector of nodes(keeps the name and the index of the node)
    graf->nodes = (g_node **)malloc((num_nodes + 1) * sizeof(g_node *));
    for (i = 1; i <= num_nodes; i++) {
        graf->nodes[i] = (g_node *)malloc(1 * sizeof(g_node));
        graf->nodes[i]->num = i;
        strcpy(graf->nodes[i]->name, "-");
    }

    return graf;
}

void dfs(Graf *graf, int current_node, int *visited) {
    // Now the node is visited
    visited[current_node] = 1;

    // Visit all the neighbors of the current node
    l_node *neighbor = graf->link[current_node]->head;
    int i;
    for (i = 0; i < graf->link[current_node]->dimension; i++) {
        if (visited[neighbor->num] == 0) dfs(graf, neighbor->num, visited);
        neighbor = neighbor->next;
    }
}

// Return the number of connex components in the graph -- USES DFS
int connex_components(Graf *graf, int num_nodes) {
    // Declare a vector which determines if a node has been visited yet
    int *visited = (int *)malloc((num_nodes + 1) * sizeof(int));
    int i;
    for (i = 1; i <= num_nodes; i++) visited[i] = 0;

    // Count the number of connex components
    int connex = 0;
    for (i = 1; i <= num_nodes; i++) {
        if (visited[i] == 0) {
            dfs(graf, i, visited);
            connex++;
        }
    }

    // Free visited vector
    free(visited);

    return connex;
}

// Return the index of the minimum value in the array, cosidering only unvisited
// positions FOR INT *
int min_index(int *v, int n, int *visited) {
    int i, min_value = INT_MAX, min_index = -1;
    for (i = 1; i <= n; i++) {
        if (v[i] < min_value && visited[i] == 0) {
            min_value = v[i];
            min_index = i;
        }
    }

    return min_index;
}

// Return the index of the minimum value in the array, cosidering only unvisited
// positions FOR FLOAT *
int min_index_float(float *v, int n, int *visited) {
    int i, min_value = INT_MAX, min_index = -1;
    for (i = 1; i <= n; i++) {
        if (v[i] < min_value && visited[i] == 0) {
            min_value = v[i];
            min_index = i;
        }
    }

    return min_index;
}

// Computes the minimum cost of a spanning tree for a given graph
int prims(Graf *graf, int *visited, int num_nodes) {
    // Initiate the key vector
    int *key = (int *)malloc((num_nodes + 1) * sizeof(int));
    int i;
    for (i = 1; i <= num_nodes; i++) key[i] = INT_MAX;

    // Find the node from where to start the algorithm(the first one which is
    // not visited)
    int node_index;
    for (node_index = 1; node_index <= num_nodes; node_index++) {
        if (visited[node_index] == 0) break;
    }

    // The cost to the starting node is 0
    key[node_index] = 0;
    while (1) {
        // Calculate the node with the minimum cost which does not make a cicle
        node_index = min_index(key, num_nodes, visited);

        // If there are no more nodes with less than infinite cost break
        if (node_index == -1) break;

        // Now the node is visited
        visited[node_index] = 1;

        // Update the key array considering new costs
        l_node *neighbor = graf->link[node_index]->head;
        for (i = 1; i <= graf->link[node_index]->dimension; i++) {
            if (visited[neighbor->num] == 0 &&
                neighbor->cost < key[neighbor->num]) {
                key[neighbor->num] = neighbor->cost;
            }

            neighbor = neighbor->next;
        }
    }

    // Compute the cost as the sum of all elements in key that are smaller than
    // infinite
    int cost = 0;
    for (i = 1; i <= num_nodes; i++) {
        if (key[i] < INT_MAX) cost += key[i];
    }

    // Free key array
    free(key);

    return cost;
}

// Used for qsort
int cmpfunction(const void *a, const void *b) {
    return (*(int *)a - *(int *)b);
}

// Dijktra Algorithm
void dijktra(Graf *graf, int source_index, int dest_index, int *prev) {
    // Alloc and initiate arrays we will need
    int i;
    int *visited = (int *)malloc((graf->dimension + 1) * sizeof(int));
    for (i = 1; i <= graf->dimension; i++) visited[i] = 0;
    float *score = (float *)malloc((graf->dimension + 1) * sizeof(float));
    for (i = 1; i <= graf->dimension; i++) score[i] = INT_MAX;

    // The score to reach the source is 0
    score[source_index] = 0;

    // While destination is not visited
    while (visited[dest_index] == 0) {
        // Get the index of the minimum score value
        int current_index = min_index_float(score, graf->dimension, visited);

        // Now the node is visited
        visited[current_index] = 1;

        // For every neighbor update the score if it is the case
        l_node *node = graf->link[current_index]->head;
        while (node != NULL) {
            int next_index = node->num;
            float computed_score =
                (float)node->cost / graf->nodes[next_index]->depth;

            if (score[current_index] + computed_score < score[next_index]) {
                score[next_index] = score[current_index] + computed_score;
                prev[next_index] = current_index;
            }

            node = node->next;
        }
    }

    // Free used arrays
    free(visited);
    free(score);
}

// Print the road from source_index to dest_index using "prev" array
void print_road(FILE *fid, Graf *graf, int source_index, int dest_index,
                int *prev) {
    if (source_index == dest_index) {
        fprintf(fid, "%s", graf->nodes[source_index]->name);
    } else {
        print_road(fid, graf, prev[source_index], dest_index, prev);
        fprintf(fid, " %s", graf->nodes[source_index]->name);
    }
}

// Return the cost of the road from source_index to dest_index using "prev"
// array
int road_cost(Graf *graf, int source_index, int dest_index, int *prev) {
    if (source_index == dest_index) return 0;

    // Find cost of the link between prev[source_index] -> source_index
    l_node *current_node = graf->link[prev[source_index]]->head;
    int current_cost = INT_MAX;
    while (current_node != NULL) {
        if (current_node->num == source_index &&
            current_node->cost < current_cost) {
            current_cost = current_node->cost;
        }

        current_node = current_node->next;
    }

    // Add the cost of the link and continue in recursivity
    return current_cost + road_cost(graf, prev[source_index], dest_index, prev);
}

int my_min(int a, int b) {
    if (a < b) return a;

    return b;
}

// Return the max weight of the boat
int weight(Graf *graf, int source_index, int dest_index, int *prev) {
    if (prev[source_index] == dest_index)
        return graf->nodes[source_index]->depth;

    return my_min(graf->nodes[source_index]->depth,
                  weight(graf, prev[source_index], dest_index, prev));
}

void free_graph(Graf *graf) {
    int i;
    for (i = 1; i <= graf->dimension; i++) {
        free_list(graf->link[i]);
    }
    free(graf->link);

    for (i = 1; i <= graf->dimension; i++) {
        free(graf->nodes[i]);
    }
    free(graf->nodes);

    free(graf);
}

int main(int argc, char *argv[]) {
    // Open input/output files
    FILE *fid_in = fopen("tema3.in", "r");
    FILE *fid_out = fopen("tema3.out", "w");

    // Read input
    int num_nodes, num_links, i;
    fscanf(fid_in, "%d %d", &num_nodes, &num_links);

    // Declare and initiate the graph
    Graf *graf = init_graph(num_nodes);

    // Check for user input
    if (argc != 2) {
        printf("Invalid number of arguments\n");
        printf("Usage: ./tema3 [1 2]\n");
        goto final;
    }

    if (strcmp(argv[1], "1") == 0) {
        // Read all the links
        for (i = 0; i < num_links; i++) {
            int cost;
            char node1_name[MAX_NAME], node2_name[MAX_NAME];
            fscanf(fid_in, "%s %s %d", node1_name, node2_name, &cost);

            // Get the indexes of the two given nodes
            int node1_num = get_index(graf->nodes, node1_name, num_nodes);
            int node2_num = get_index(graf->nodes, node2_name, num_nodes);

            // Add links in graph
            add_list(graf->link[node1_num], node2_num, cost);
            add_list(graf->link[node2_num], node1_num, cost);
        }

        // Print the number of connex components
        int connex = connex_components(graf, num_nodes);
        fprintf(fid_out, "%d\n", connex);

        // Binary array - v[i] = 0 -> not visited; v[i] = 1 -> visited
        int *visited = (int *)malloc((num_nodes + 1) * sizeof(int));
        for (i = 1; i <= num_nodes; i++) visited[i] = 0;

        // Print the minimum cost of a spanning tree -- USES PRIMS ALGORITHM
        int *cost = (int *)malloc(connex * sizeof(int));
        for (i = 0; i < connex; i++) cost[i] = prims(graf, visited, num_nodes);
        qsort(cost, connex, sizeof(int), cmpfunction);
        for (i = 0; i < connex; i++) fprintf(fid_out, "%d\n", cost[i]);

        // Free cost and visited arrays
        free(visited);
        free(cost);
    } else if (strcmp(argv[1], "2") == 0) {
        int island_index, boat_index;

        // Read all the links
        for (i = 0; i < num_links; i++) {
            int cost;
            char node1_name[MAX_NAME], node2_name[MAX_NAME];
            fscanf(fid_in, "%s %s %d", node1_name, node2_name, &cost);

            // Get the indexes of the two given nodes
            int node1_num = get_index(graf->nodes, node1_name, num_nodes);
            int node2_num = get_index(graf->nodes, node2_name, num_nodes);

            // Add links in graph
            add_list(graf->link[node1_num], node2_num, cost);
        }

        // Read node depths
        for (i = 1; i <= num_nodes; i++) {
            int depth;
            char node_name[MAX_NAME];
            fscanf(fid_in, "%s %d", node_name, &depth);

            // Get the indexes of the given node
            int node_num = get_index(graf->nodes, node_name, num_nodes);

            // Verify if the read node is the island
            if (strcmp(node_name, "Insula") == 0) island_index = node_num;

            // Verify if the read node is the boat
            if (strcmp(node_name, "Corabie") == 0) boat_index = node_num;

            // Add links in graph
            graf->nodes[node_num]->depth = depth;
        }

        // Read the weight of the treasure
        int treasure_weight = 0;
        fscanf(fid_in, "%d", &treasure_weight);

        // Binary array - v[i] = 0 -> not visited; v[i] = 1 -> visited
        int *visited = (int *)malloc((num_nodes + 1) * sizeof(int));
        for (i = 1; i <= num_nodes; i++) visited[i] = 0;

        // Verify if the crew can reach the island from the boat
        dfs(graf, boat_index, visited);
        if (visited[island_index] == 0) {
            free(visited);
            fprintf(fid_out, "Echipajul nu poate ajunge la insula\n");
            goto final;
        }

        for (i = 1; i <= num_nodes; i++) visited[i] = 0;

        // Verify if the crew can reach the boat from the island
        dfs(graf, island_index, visited);
        if (visited[boat_index] == 0) {
            free(visited);
            fprintf(
                fid_out,
                "Echipajul nu poate transporta comoara inapoi la corabie\n");
            goto final;
        }

        // Use Dijktra to find the best path source_index -> dest_index
        int *prev = (int *)malloc((num_nodes + 1) * sizeof(int));
        for (i = 1; i <= num_nodes; i++) prev[i] = -1;

        dijktra(graf, island_index, boat_index, prev);

        // Print the determined road
        print_road(fid_out, graf, boat_index, island_index, prev);
        fprintf(fid_out, "\n");

        // Print the cost of the road
        fprintf(fid_out, "%d\n",
                road_cost(graf, boat_index, island_index, prev));

        // Print the max weight of the boat
        int max_weight = weight(graf, prev[boat_index], island_index, prev);
        fprintf(fid_out, "%d\n", max_weight);

        // Print the number of roads
        fprintf(fid_out, "%d\n", treasure_weight / max_weight);

        // Free used arrays
        free(visited);
        free(prev);
    } else {
        printf("Usage: ./tema3 [1 2]\n");
    }

final:

    // Free the memory allocated for the graph
    free_graph(graf);

    // Close the files
    fclose(fid_in);
    fclose(fid_out);

    return 0;
}