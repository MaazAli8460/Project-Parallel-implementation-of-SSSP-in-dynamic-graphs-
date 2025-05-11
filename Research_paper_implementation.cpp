#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <algorithm>
#include <omp.h>
#include <fstream>
#include <sstream>
#include <random>
#include <set>
#include <chrono>
#include <iomanip>
#include <mutex>
#include <atomic>
#include <memory>

// Data structure for a weighted edge
struct Edge {
    int source;
    int destination;
    double weight;
    
    Edge(int s, int d, double w) : source(s), destination(d), weight(w) {}
    
    bool operator<(const Edge& other) const {
        if (source != other.source)
            return source < other.source;
        if (destination != other.destination)
            return destination < other.destination;
        return weight < other.weight;
    }
    
    bool operator==(const Edge& other) const {
        return source == other.source && destination == other.destination;
    }
};

// Data structure for the graph
class Graph {
private:
    int numVertices;
    std::vector<std::vector<std::pair<int, double>>> adjacencyList;
    std::set<Edge> edgeSet;
    mutable std::mutex edgeSetMutex; // Mutable for const methods
    mutable std::vector<std::unique_ptr<std::mutex>> adjMutexes; // Mutable for const methods

public:
    Graph(int n) : numVertices(n) {
        adjacencyList.resize(n);
        adjMutexes.resize(n);
        for (int i = 0; i < n; ++i) {
            adjMutexes[i] = std::make_unique<std::mutex>();
        }
    }

    // Custom copy constructor
    Graph(const Graph& other) : numVertices(other.numVertices) {
        adjacencyList.resize(numVertices);
        adjMutexes.resize(numVertices);
        for (int i = 0; i < numVertices; ++i) {
            adjMutexes[i] = std::make_unique<std::mutex>();
            adjacencyList[i] = other.adjacencyList[i]; // Deep copy adjacency list
        }
        {
            std::lock_guard<std::mutex> lock(other.edgeSetMutex);
            edgeSet = other.edgeSet; // Deep copy edge set
        }
    }

    void addEdge(int source, int dest, double weight) {
        if (source < 0 || source >= numVertices || dest < 0 || dest >= numVertices) {
            return; // Invalid vertices
        }

        // Check if edge already exists
        if (hasEdge(source, dest)) {
            return;
        }
        
        {
            std::lock_guard<std::mutex> lock1(*adjMutexes[source]);
            std::lock_guard<std::mutex> lock2(*adjMutexes[dest]);
            adjacencyList[source].push_back(std::make_pair(dest, weight));
            adjacencyList[dest].push_back(std::make_pair(source, weight));
        }
        
        {
            std::lock_guard<std::mutex> lock(edgeSetMutex);
            edgeSet.insert(Edge(source, dest, weight));
            edgeSet.insert(Edge(dest, source, weight));
        }
    }

    void removeEdge(int source, int dest) {
        if (source < 0 || source >= numVertices || dest < 0 || dest >= numVertices) {
            return; // Invalid vertices
        }

        {
            std::lock_guard<std::mutex> lock1(*adjMutexes[source]);
            std::lock_guard<std::mutex> lock2(*adjMutexes[dest]);
            adjacencyList[source].erase(
                std::remove_if(
                    adjacencyList[source].begin(),
                    adjacencyList[source].end(),
                    [dest](const std::pair<int, double>& p) { return p.first == dest; }
                ),
                adjacencyList[source].end()
            );
            adjacencyList[dest].erase(
                std::remove_if(
                    adjacencyList[dest].begin(),
                    adjacencyList[dest].end(),
                    [source](const std::pair<int, double>& p) { return p.first == source; }
                ),
                adjacencyList[dest].end()
            );
        }
        
        {
            std::lock_guard<std::mutex> lock(edgeSetMutex);
            edgeSet.erase(Edge(source, dest, 0));
            edgeSet.erase(Edge(dest, source, 0));
        }
    }

    bool hasEdge(int source, int dest) const {
        if (source < 0 || source >= numVertices || dest < 0 || dest >= numVertices) {
            return false;
        }
        std::lock_guard<std::mutex> lock(*adjMutexes[source]);
        for (const auto& edge : adjacencyList[source]) {
            if (edge.first == dest) {
                return true;
            }
        }
        return false;
    }

    const std::vector<std::pair<int, double>>& getNeighbors(int vertex) const {
        if (vertex < 0 || vertex >= numVertices) {
            static std::vector<std::pair<int, double>> empty;
            return empty;
        }
        return adjacencyList[vertex];
    }

    int getNumVertices() const {
        return numVertices;
    }
    
    Edge getRandomEdge(std::mt19937& gen) const {
        std::lock_guard<std::mutex> lock(edgeSetMutex);
        if (edgeSet.empty()) {
            return Edge(-1, -1, 0);
        }
        
        std::uniform_int_distribution<int> dist(0, edgeSet.size() - 1);
        int randomIndex = dist(gen);
        
        auto it = edgeSet.begin();
        std::advance(it, randomIndex);
        
        return *it;
    }
    
    std::vector<Edge> getAllEdges() const {
        std::vector<Edge> allEdges;
        std::lock_guard<std::mutex> lock(edgeSetMutex);
        for (int i = 0; i < numVertices; ++i) {
            std::lock_guard<std::mutex> adjLock(*adjMutexes[i]);
            for (const auto& neighbor : adjacencyList[i]) {
                if (i < neighbor.first) {
                    allEdges.push_back(Edge(i, neighbor.first, neighbor.second));
                }
            }
        }
        return allEdges;
    }
};

// Data structure for SSSP tree
class SSSPTree {
public:
    std::vector<int> parent;
    std::vector<double> distance;
    std::vector<bool> affectedDel;
    std::vector<bool> affected;
    int source;

    SSSPTree(int numVertices, int src) : source(src) {
        parent.resize(numVertices, -1);
        distance.resize(numVertices, std::numeric_limits<double>::infinity());
        affectedDel.resize(numVertices, false);
        affected.resize(numVertices, false);
        distance[source] = 0;
    }

    bool isEdgeInTree(int u, int v) const {
        return (parent[v] == u || parent[u] == v);
    }
    
    void computeSSSP(const Graph& graph) {
        std::fill(distance.begin(), distance.end(), std::numeric_limits<double>::infinity());
        std::fill(parent.begin(), parent.end(), -1);
        distance[source] = 0;
        
        std::priority_queue<std::pair<double, int>, 
                           std::vector<std::pair<double, int>>, 
                           std::greater<std::pair<double, int>>> pq;
        
        pq.push(std::make_pair(0, source));
        
        while (!pq.empty()) {
            double dist = pq.top().first;
            int u = pq.top().second;
            pq.pop();
            
            if (dist > distance[u]) {
                continue;
            }
            
            for (const auto& neighbor : graph.getNeighbors(u)) {
                int

 v = neighbor.first;
                double weight = neighbor.second;
                
                if (distance[v] > distance[u] + weight) {
                    distance[v] = distance[u] + weight;
                    parent[v] = u;
                    pq.push(std::make_pair(distance[v], v));
                }
            }
        }
    }
};

void processChangedEdges(const Graph& originalGraph, Graph& updatedGraph, 
                        SSSPTree& tree, const std::vector<Edge>& deletedEdges, 
                        const std::vector<Edge>& insertedEdges) {
    std::fill(tree.affectedDel.begin(), tree.affectedDel.end(), false);
    std::fill(tree.affected.begin(), tree.affected.end(), false);

    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < deletedEdges.size(); ++i) {
        const Edge& edge = deletedEdges[i];
        int u = edge.source;
        int v = edge.destination;
        
        if (tree.isEdgeInTree(u, v)) {
            int affectedVertex = (tree.distance[u] > tree.distance[v]) ? u : v;
            
            #pragma omp critical
            {
                tree.distance[affectedVertex] = std::numeric_limits<double>::infinity();
                tree.affectedDel[affectedVertex] = true;
                tree.affected[affectedVertex] = true;
                tree.parent[affectedVertex] = -1;
            }
        }
        
        updatedGraph.removeEdge(u, v);
    }

    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < insertedEdges.size(); ++i) {
        const Edge& edge = insertedEdges[i];
        int u = edge.source;
        int v = edge.destination;
        double weight = edge.weight;
        
        int x = (tree.distance[u] > tree.distance[v]) ? v : u;
        int y = (tree.distance[u] > tree.distance[v]) ? u : v;
        
        if (tree.distance[y] > tree.distance[x] + weight) {
            #pragma omp critical
            {
                tree.distance[y] = tree.distance[x] + weight;
                tree.parent[y] = x;
                tree.affected[y] = true;
            }
        }
        
        updatedGraph.addEdge(u, v, weight);
    }
}

void updateAffectedVerticesSynchronous(const Graph& graph, SSSPTree& tree) {
    // Step 1: Update subtrees affected by deletion
    std::atomic<bool> hasAffectedDelVertices(true);
    
    while (hasAffectedDelVertices) {
        hasAffectedDelVertices = false;
        std::vector<int> currentAffectedDelVertices;
        
        for (int v = 0; v < graph.getNumVertices(); ++v) {
            if (tree.affectedDel[v]) {
                currentAffectedDelVertices.push_back(v);
                tree.affectedDel[v] = false;
            }
        }
        
        #pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < currentAffectedDelVertices.size(); ++i) {
            int v = currentAffectedDelVertices[i];
            
            for (int c = 0; c < graph.getNumVertices(); ++c) {
                if (tree.parent[c] == v) {
                    #pragma omp critical
                    {
                        tree.distance[c] = std::numeric_limits<double>::infinity();
                        tree.affectedDel[c] = true;
                        tree.affected[c] = true;
                        hasAffectedDelVertices = true;
                    }
                }
            }
        }
    }
    
    // Step 2: Update distances
    std::atomic<bool> hasAffectedVertices(true);
    
    while (hasAffectedVertices) {
        hasAffectedVertices = false;
        std::vector<int> currentAffectedVertices;
        
        for (int v = 0; v < graph.getNumVertices(); ++v) {
            if (tree.affected[v]) {
                currentAffectedVertices.push_back(v);
                tree.affected[v] = false;
            }
        }
        
        #pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < currentAffectedVertices.size(); ++i) {
            int v = currentAffectedVertices[i];
            
            for (const auto& neighborPair : graph.getNeighbors(v)) {
                int n = neighborPair.first;
                double weight = neighborPair.second;
                
                bool updated = false;
                if (tree.distance[n] > tree.distance[v] + weight) {
                    #pragma omp critical
                    {
                        if (tree.distance[n] > tree.distance[v] + weight) {
                            tree.distance[n] = tree.distance[v] + weight;
                            tree.parent[n] = v;
                            tree.affected[n] = true;
                            updated = true;
                        }
                    }
                }
                
                if (tree.distance[v] > tree.distance[n] + weight) {
                    #pragma omp critical
                    {
                        if (tree.distance[v] > tree.distance[n] + weight) {
                            tree.distance[v] = tree.distance[n] + weight;
                            tree.parent[v] = n;
                            tree.affected[v] = true;
                            updated = true;
                        }
                    }
                }
                
                if (updated) {
                    hasAffectedVertices = true;
                }
            }
        }
    }
}

// Placeholder asynchronous update (simplified, with relaxed synchronization)
void updateAffectedVerticesAsynchronous(const Graph& graph, SSSPTree& tree, int asynchronyLevel) {
    // Process affected verticesAGONwith limited synchronization
    std::atomic<bool> hasAffectedVertices(true);
    
    while (hasAffectedVertices) {
        hasAffectedVertices = false;
        std::vector<int> currentAffectedVertices;
        
        for (int v = 0; v < graph.getNumVertices(); ++v) {
            if (tree.affected[v] || tree.affectedDel[v]) {
                currentAffectedVertices.push_back(v);
                tree.affected[v] = false;
                tree.affectedDel[v] = false;
            }
        }
        
        // Shuffle vertices to simulate asynchrony
        std::random_device rd;
        std::mt19937 gen(rd());
        std::shuffle(currentAffectedVertices.begin(), currentAffectedVertices.end(), gen);
        
        #pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < currentAffectedVertices.size(); ++i) {
            int v = currentAffectedVertices[i];
            
            for (const auto& neighborPair : graph.getNeighbors(v)) {
                int n = neighborPair.first;
                double weight = neighborPair.second;
                
                bool updated = false;
                // Relaxed update (less frequent critical sections based on asynchronyLevel)
                if ((rand() % asynchronyLevel) == 0) {
                    #pragma omp critical
                    {
                        if (tree.distance[n] > tree.distance[v] + weight) {
                            tree.distance[n] = tree.distance[v] + weight;
                            tree.parent[n] = v;
                            tree.affected[n] = true;
                            updated = true;
                        }
                        if (tree.distance[v] > tree.distance[n] + weight) {
                            tree.distance[v] = tree.distance[n] + weight;
                            tree.parent[v] = n;
                            tree.affected[v] = true;
                            updated = true;
                        }
                    }
                } else {
                    // Non-atomic update (simulates asynchrony, may cause validation failures)
                    if (tree.distance[n] > tree.distance[v] + weight) {
                        tree.distance[n] = tree.distance[v] + weight;
                        tree.parent[n] = v;
                        tree.affected[n] = true;
                        updated = true;
                    }
                    if (tree.distance[v] > tree.distance[n] + weight) {
                        tree.distance[v] = tree.distance[n] + weight;
                        tree.parent[v] = n;
                        tree.affected[v] = true;
                        updated = true;
                    }
                }
                
                if (updated) {
                    hasAffectedVertices = true;
                }
            }
        }
    }
}

void generateRandomChanges(const Graph& graph, int numChanges, int insertionPercentage,
                         std::vector<Edge>& deletedEdges, std::vector<Edge>& insertedEdges) {
    if (insertionPercentage < 0 || insertionPercentage > 100) {
        std::cerr << "Error: insertionPercentage must be between 0 and 100" << std::endl;
        return;
    }
    
    int numInsertions = (numChanges * insertionPercentage) / 100;
    int numDeletions = numChanges - numInsertions;
    
    std::random_device rd;
    std::mt19937 gen(rd());
    
    std::vector<Edge> allEdges = graph.getAllEdges();
    numDeletions = std::min(numDeletions, (int)allEdges.size());
    
    std::shuffle(allEdges.begin(), allEdges.end(), gen);
    
    for (int i = 0; i < numDeletions; ++i) {
        deletedEdges.push_back(allEdges[i]);
    }
    
    std::uniform_int_distribution<int> vertexDist(0, graph.getNumVertices() - 1);
    std::uniform_real_distribution<double> weightDist(1.0, 10.0);
    std::set<std::pair<int, int>> insertedEdgeSet;
    
    while (insertedEdges.size() < numInsertions) {
        int src = vertexDist(gen);
        int dest = vertexDist(gen);
        double weight = weightDist(gen);
        
        if (src == dest || graph.hasEdge(src, dest) || 
            insertedEdgeSet.count(std::make_pair(src, dest)) ||
            insertedEdgeSet.count(std::make_pair(dest, src))) {
            continue;
        }
        
        insertedEdges.push_back(Edge(src, dest, weight));
        insertedEdgeSet.insert(std::make_pair(src, dest));
        insertedEdgeSet.insert(std::make_pair(dest, src));
    }
}

bool validateSSSP(const SSSPTree& computedTree, const SSSPTree& referenceTree) {
    for (size_t i = 0; i < computedTree.distance.size(); ++i) {
        if (std::abs(computedTree.distance[i] - referenceTree.distance[i]) > 1e-9) {
            std::cout << "Validation failed for vertex " << i 
                     << ": Computed distance = " << computedTree.distance[i]
                     << ", Reference distance = " << referenceTree.distance[i] << std::endl;
            return false;
        }
    }
    return true;
}

void updateSSSP(const Graph& originalGraph, int source, 
               const std::vector<Edge>& deletedEdges, 
               const std::vector<Edge>& insertedEdges, 
               bool useAsync, int asynchronyLevel) {
    Graph updatedGraph = originalGraph;
    SSSPTree tree(originalGraph.getNumVertices(), source);
    
    tree.computeSSSP(originalGraph);
    
    auto startTime = std::chrono::high_resolution_clock::now();
    
    auto processStartTime = std::chrono::high_resolution_clock::now();
    processChangedEdges(originalGraph, updatedGraph, tree, deletedEdges, insertedEdges);
    auto processEndTime = std::chrono::high_resolution_clock::now();
    double processingTime = std::chrono::duration<double>(processEndTime - processStartTime).count();
    
    auto updateStartTime = std::chrono::high_resolution_clock::now();
    if (useAsync) {
        updateAffectedVerticesAsynchronous(updatedGraph, tree, asynchronyLevel);
    } else {
        updateAffectedVerticesSynchronous(updatedGraph, tree);
    }
    auto updateEndTime = std::chrono::high_resolution_clock::now();
    double updateTime = std::chrono::duration<double>(updateEndTime - updateStartTime).count();
    
    auto endTime = std::chrono::high_resolution_clock::now();
    double totalTime = std::chrono::duration<double>(endTime - startTime).count();
    
    SSSPTree validationTree(updatedGraph.getNumVertices(), source);
    validationTree.computeSSSP(updatedGraph);
    bool isValid = validateSSSP(tree, validationTree);
    
    std::cout << "Timing Results (" << (useAsync ? "Asynchronous" : "Synchronous") << "):" << std::endl;
    std::cout << "--------------------------------------------" << std::endl;
    std::cout << "Processing changed edges time: " << std::fixed << std::setprecision(6) << processingTime << " seconds" << std::endl;
    std::cout << "Updating affected vertices time: " << std::fixed << std::setprecision(6) << updateTime << " seconds" << std::endl;
    std::cout << "Total update time: " << std::fixed << std::setprecision(6) << totalTime << " seconds" << std::endl;
    std::cout << "SSSP Validation: " << (isValid ? "Passed" : "Failed") << std::endl;
    
    std::cout << "\nStatistics:" << std::endl;
    std::cout << "---------------" << std::endl;
    std::cout << "Number of vertices: " << originalGraph.getNumVertices() << std::endl;
    std::cout << "Number of deletions: " << deletedEdges.size() << std::endl;
    std::cout << "Number of insertions: " << insertedEdges.size() << std::endl;
    if (useAsync) {
        std::cout << "Level of asynchrony: " << asynchronyLevel << std::endl;
    }
}

int main(int argc, char* argv[]) {
    if (argc != 4 && argc != 5) {
        std::cerr << "Usage: " << argv[0] << " <data.txt> <numOfChanges> <insertionPercentage> [asynchronyLevel]" << std::endl;
        return 1;
    }
    
    std::string dataFilename = argv[1];
    int numChanges = std::stoi(argv[2]);
    int insertionPercentage = std::stoi(argv[3]);
    bool useAsync = (argc == 5);
    int asynchronyLevel = useAsync ? std::stoi(argv[4]) : 1;
    
    if (insertionPercentage < 0 || insertionPercentage > 100) {
        std::cerr << "Error: insertionPercentage must be between 0 and 100" << std::endl;
        return 1;
    }
    if (numChanges < 0) {
        std::cerr << "Error: numOfChanges must be non-negative" << std::endl;
        return 1;
    }
    if (useAsync && asynchronyLevel <= 0) {
        std::cerr << "Error: asynchronyLevel must be positive" << std::endl;
        return 1;
    }
    
    std::ifstream dataFile(dataFilename);
    if (!dataFile.is_open()) {
        std::cerr << "Error: Cannot open file " << dataFilename << std::endl;
        return 1;
    }
    
    int maxVertex = -1;
    int src, dest;
    double weight;
    std::string line;
    
    while (std::getline(dataFile, line)) {
        std::istringstream iss(line);
        if (!(iss >> src >> dest >> weight)) {
            continue;
        }
        maxVertex = std::max(maxVertex, std::max(src, dest));
    }
    
    int numVertices = maxVertex + 1;
    
    dataFile.clear();
    dataFile.seekg(0, std::ios::beg);
    
    Graph graph(numVertices);
    
    while (std::getline(dataFile, line)) {
        std::istringstream iss(line);
        if (!(iss >> src >> dest >> weight)) {
            continue;
        }
        graph.addEdge(src, dest, weight);
    }
    
    std::vector<Edge> deletedEdges, insertedEdges;
    generateRandomChanges(graph, numChanges, insertionPercentage, deletedEdges, insertedEdges);
    
    // Run both synchronous and asynchronous updates if asynchronyLevel is provided
    std::cout << "\n=== incisionsUpdate ===\n" << std::endl;
    updateSSSP(graph, 0, deletedEdges, insertedEdges, false, asynchronyLevel);
    
    if (useAsync) {
        std::cout << "\n=== Asynchronous Update ===\n" << std::endl;
        updateSSSP(graph, 0, deletedEdges, insertedEdges, true, asynchronyLevel);
    }
    
    return 0;
}
