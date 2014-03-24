import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Random;
import java.util.Stack;

public class dijkstra {

	/*
	 * Class VNode is a class described vertices in a graph.
	 */
	class VNode {
		private int key;
		private ENode next;
		
		public VNode(int key) {
			this.key = key;
			next = null;
		}
		
		public int getKey() {
			return key;
		}
		
		public ENode getNext() {
			return next;
		}
		
		public void setNext(ENode eNode) {
			this.next = eNode;
		}
	}
	
	/*
	 * Class ENode is a class described the nodes on edges in a graph.
	 */
	class ENode {
		private int key;
		private ENode next;
		private int cost;
		
		public ENode(int key, int cost) {
			this.key = key;
			this.cost = cost;
			this.next = null;
		}
		
		public int getKey() {
			return key;
		}
		
		public ENode getNext() {
			return next;
		}
		
		public int getCost() {
			return cost;
		}
		
		public void setNext(ENode eNode) {
			this.next = eNode;
		}
	}
	
	/*
	 * Class Graph is a class described a undirected graph using adjacent list.
	 */
	class Graph {
		private List<VNode> vertices;
		private int verticesNum;
		private int edgeNum;
		
		public Graph(int verticesNum, int edgeNum) {
			vertices = new ArrayList<VNode>();
			this.verticesNum = verticesNum;
			this.edgeNum = edgeNum;
		}
		
		public int getVertexNum() {
			return verticesNum;
		}
		
		public int getEdgeNum() {
			return edgeNum;
		}
		
		public VNode getVertex(int key) {
			return vertices.get(key);
		}
		
		/*
		 * isDuplicate(int, int) method is to judge if the edge exists in the graph
		 * The two parameters are the start node and the end node
		 * It will return boolean value
		 */
		public boolean isDuplicateEdge(int vNodeKey, int eNodeKey) {
			if (vertices.get(vNodeKey).getNext() == null) {
				return false;
			} else {
				ENode nodePointer = vertices.get(vNodeKey).getNext();
				while(nodePointer != null) {
					if (eNodeKey == nodePointer.getKey()) {
						return true;
					} else {
						nodePointer = nodePointer.getNext();
					}
				}
				return false;
			}
		}
		
		/*
		 * addEdge(int, int, int) method is to add edge between two nodes in a graph
		 * The three parameters are the start node ,the end node and the cost
		 */
		public void addEdge(int vNodeKey, int eNodeKey, int cost) {
			ENode eNode = new ENode(eNodeKey, cost);
			if (vertices.get(vNodeKey).getNext() == null) {
				vertices.get(vNodeKey).setNext(eNode);
			} else {
				ENode nodePointer = vertices.get(vNodeKey).getNext();
				while (nodePointer.getNext() != null) {
					nodePointer = nodePointer.getNext();
				}
				nodePointer.setNext(eNode);
			}
			eNode = new ENode(vNodeKey, cost);
			if (vertices.get(eNodeKey).getNext() == null) {
				vertices.get(eNodeKey).setNext(eNode);
			} else {
				ENode nodePointer = vertices.get(eNodeKey).getNext();
				while (nodePointer.getNext() != null) {
					nodePointer = nodePointer.getNext();
				}
				nodePointer.setNext(eNode);
			}
		}
		
		/*
		 * The isConnected Method is to use DFS to judge whether the graph is connected or not.
		 * It will return boolean value.
		 */
		public boolean isConnected() {
			int visitedNodeNum = 0;
			int vNodeIndex = 0;
			Stack<Integer> verticesStack = new Stack<Integer>();
			HashMap<Integer, Boolean> visitedNode = new HashMap<Integer, Boolean>();
			verticesStack.push(0);
			visitedNode.put(0, true);
			visitedNodeNum++;
			VNode vNodePointer;
			while(visitedNodeNum < verticesNum && !verticesStack.isEmpty()) {
				boolean allVisited = false;
				vNodePointer = vertices.get(vNodeIndex);
				if (vNodePointer.getNext() != null) {
					ENode eNodePointer = vNodePointer.getNext();
					while (visitedNode.get(eNodePointer.getKey()) != null) {
						eNodePointer = eNodePointer.getNext();
						if (eNodePointer == null) {
							allVisited = true;
							vNodeIndex = verticesStack.pop();
							break;
						}
					}
					if (!allVisited) {
						if (visitedNode.get(vNodeIndex) == true) {
							verticesStack.push(vNodeIndex);
						} else {
							visitedNode.put(vNodeIndex, true);
						}
						visitedNode.put(eNodePointer.getKey(), false);
						verticesStack.push(eNodePointer.getKey());
						visitedNodeNum++;
						vNodeIndex = eNodePointer.getKey();
					}
				} else {
					vNodeIndex = verticesStack.pop();
				}
			}
			if (visitedNodeNum == verticesNum) {
				return true;
			} else {
				return false;
			}
		}
		
		/*
		 * printGraph method is to print out the graph as an adjacent list
		 */
		public void printGraph() {
			for (int n = 0; n < verticesNum; n++) {
				System.out.print(vertices.get(n).getKey() + " ");
				if (vertices.get(n).getNext() != null) {
					ENode nodePointer = vertices.get(n).getNext();
					System.out.print(nodePointer.getKey() + " ");
					while (nodePointer.getNext() != null) {
						nodePointer = nodePointer.getNext();
						System.out.print(nodePointer.getKey() + " ");
					}
				}
				System.out.print("\n");
			}
		}
	}
	
	/*
	 * The VertexInf class is designed for the table when finding the shortest paths.
	 * It contains key, distance, previous vertex and known status of the vertex.
	 */
	class VertexInf {
		
		private int key;
		private boolean known;
		private int distance;
		private int preVertexKey;
		
		public VertexInf(int key, int distance) {
			this.key = key;
			this.known = false;
			this.distance = distance;
			this.preVertexKey = 0;
		}
		
		public int getKey() {
			return key;
		}
		
		public boolean isKnown() {
			return known;
		}
		
		public int getDistance() {
			return distance;
		}
		
		public int preVertex() {
			return preVertexKey;
		}
		
		public void setKnown() {
			this.known = true;
		}
		public void setDistance(int distance) {
			this.distance = distance;
		}
		
		public void setPreVertexKey(int key) {
			this.preVertexKey = key;
		}
	}
	
	/*
	 * The FibonacciHeapNode is a class described nodes of a Fibonacci heap node.
	 */
	
	class FibonacciHeapNode<T> {
		
		T data;
		FibonacciHeapNode<T> left;
		FibonacciHeapNode<T> right;
		FibonacciHeapNode<T> parent;
		FibonacciHeapNode<T> child;
		boolean childCut;
		int degree;
		int key;
		
		public FibonacciHeapNode(T data, int key) {
			this.right = this;
			this.left = this;
			this.data = data;
			this.key = key;
		}
		
		public int getKey() {
			return key;
		}
		
		public T getData() {
			return data;
		}
		
		/*
		 * cascadingCut(FibonacciHeapNode<T>) method is to cut nodes with
		 * childCut = true from their sibling lists and inserted into top-level
		 * list.
		 * The parameter is the minimum node of this heap.
		 */
		public void cascadingCut(FibonacciHeapNode<T> minNode) {
			FibonacciHeapNode<T> parentNode = parent;
            if (parentNode != null) {
                if (childCut) {
                	parentNode.cut(this, minNode);
                	parentNode.cascadingCut(minNode);
                } else {
                    childCut = true;
                }
            }
        }
		
		/*
		 * cut(FibonacciHeapNode<T>, FibonacciHeapNode<T>) method is to cut a node from
		 * the parent, make it as a root node and decrease the degree
		 * The two parameters are the node which is going to be cut and the minimum node
		 * of this heap.
		 */
		public void cut(FibonacciHeapNode<T> node, FibonacciHeapNode<T> minNode) {
            node.left.right = node.right;
            node.right.left = node.left;
            degree--;
            
            if (degree == 0) {
                child = null;
            } else if (child == node) {
                child = node.right;
            }
            
            node.right = minNode;
            node.left = minNode.left;
            minNode.left = node;
            node.left.right = node;
            node.parent = null;
            node.childCut = false;
        }
		
		/*
		 * link(FibonacciHeapNode<T>) method is to link the node to the other node.
		 * The parameter is the parent node which will be linked to.
		 */
		public void link(FibonacciHeapNode<T> parentNode) {
            left.right = right;
            right.left = left;

            this.parent = parentNode;
            if (parentNode.child == null) {
                parentNode.child = this;
                right = this;
                left = this;
            } else {
                left = parentNode.child;
                right = parentNode.child.right;
                parentNode.child.right = this;
                right.left = this;
            }

            parentNode.degree++;

            childCut = false;
        }
	}
	
	/*
	 * The FibonacciHeap is a class described nodes of a Fibonacci heap.
	 */
	class FibonacciHeap<T> {
		
		private FibonacciHeapNode<T> minNode;
		private int nodesNum;
		private double constant = 1.0 / Math.log((1.0 + Math.sqrt(5.0)) / 2.0);
		
		public FibonacciHeap() {
		}
		
		public boolean isEmpty() {
			return minNode == null;
		}
		
		public void clear() {
			minNode = null;
			nodesNum = 0;
		}
		
		/*
		 * insert(T, int) method is to insert a new node into a heap.
		 * The two parameters are the element and key to be inserted.
		 * It will return the node.
		 */
		public FibonacciHeapNode<T> insert(T data, int key) {
	        FibonacciHeapNode<T> node = new FibonacciHeapNode<T>(data, key);

	        if (minNode != null) {
	            node.right = minNode;
	            node.left = minNode.left;
	            minNode.left = node;
	            node.left.right = node;
	            if (key < minNode.key) {
	            	minNode = node;
	            }
	        } else {
	        	minNode = node;
	        }
	        nodesNum++;
	        return node;
	    }
		
		/*
		 * decreaseKey(FibonacciHeapNode<T>, int) method is to decrease the key
		 * of the node.
		 * The two parameters are the node and its new key.
		 */
		private void decreaseKey(FibonacciHeapNode<T> node, int key) {
	        if (key > node.key) {
	            throw new IllegalArgumentException("The new key value is bigger than the old one!");
	        }
	        node.key = key;
	        FibonacciHeapNode<T> parentNode = node.parent;
	        if (parentNode != null && key < parentNode.key) {
	        	parentNode.cut(node, minNode);
	        	parentNode.cascadingCut(minNode);
	        }
	        if (key < minNode.key) {
	            minNode = node;
	        }
	    }
		
		/*
		 * removeMin() method is to remove the minimum node from the heap
		 * It will return the element in this heap.
		 */
		public T removeMin() {
	        FibonacciHeapNode<T> node = minNode;
	        if (node == null) {
	            return null;
	        }
	        if (node.child != null) {
	            node.child.parent = null;
	            for (FibonacciHeapNode<T> childNode = node.child.right; childNode != node.child; childNode = childNode.right) {
	                childNode.parent = null;
	            }
	            
	            FibonacciHeapNode<T> minLeftPointer = minNode.left;
	            FibonacciHeapNode<T> childLeftPointer = node.child.left;
	            minNode.left = childLeftPointer;
	            childLeftPointer.right = minNode;
	            node.child.left = minLeftPointer;
	            minLeftPointer.right = node.child;
	        }
	       
	        node.left.right = node.right;
	        node.right.left = node.left;
	        if (node == node.right) {
	            minNode = null;
	        } else {
	            minNode = node.right;
	            consolidate();
	        }
	        nodesNum--;
	        return node.data;
	    }
		
		/*
		 * delete(FibonacciHeapNode<T>) method is to delete a node from a heap.
		 * The parameter is the node to be deleted.
		 */
		public void delete(FibonacciHeapNode<T> node) {
			decreaseKey(node, -9999999);
			removeMin();
		}
		
		public FibonacciHeapNode<T> getMinNode() {
			return minNode;
		}
		
		/*
		 * consolidate() method is to join those trees whose degree are equal
		 * until there are no more trees of equal degree in the heap.
		 */
		private void consolidate() {
			
			int size = ((int) Math.floor(Math.log(nodesNum) * constant)) + 1;
			ArrayList<FibonacciHeapNode<T>> nodeList = new ArrayList<FibonacciHeapNode<T>>(size);
			for (int i = 0; i < size; i++) {
				nodeList.add(null);
			}
			
			FibonacciHeapNode<T> startNode = minNode;
			FibonacciHeapNode<T> currentNodePointer = minNode;
	        do {
	        	FibonacciHeapNode<T> node = currentNodePointer;
	        	FibonacciHeapNode<T> nodeRightPointer = currentNodePointer.right;
	            int deg = node.degree;
	            while (nodeList.get(deg) != null) {
	            	FibonacciHeapNode<T> existedNodePointer = nodeList.get(deg);
	                if (node.key > existedNodePointer.key) {
	                	FibonacciHeapNode<T> tempPointer = existedNodePointer;
	                	existedNodePointer = node;
	                	node = tempPointer;
	                }
	                if (existedNodePointer == startNode) {
	                    startNode = startNode.right;
	                }
	                if (existedNodePointer == nodeRightPointer) {
	                	nodeRightPointer = nodeRightPointer.right;
	                }
	                existedNodePointer.link(node);
	                nodeList.set(deg, null);
	                deg++;
	            }

	            nodeList.set(deg, node);
	            currentNodePointer = nodeRightPointer;
	        } while (currentNodePointer != startNode);

	        minNode = startNode;
	        
	        for (FibonacciHeapNode<T> node : nodeList) {
	            if (node != null && node.key < minNode.key) {
	                minNode = node;
	            }
	        }
	    }
	}
	
	/*
	 * randomGraphGenerate(int, int) method is to generate a random graph
	 * according to the amount of vertices and edges.
	 * Two parameters are the amount of vertices and edges.
	 * It will return the generated graph.
	 */
	private Graph randomGraphGenerate(int verticesNum, int edgesNum) {
		Graph randomGraph = new Graph(verticesNum, edgesNum);
		for (int key = 0; key < verticesNum; key++) {
			VNode vertex = new VNode(key);
			randomGraph.vertices.add(vertex);
		}
		if (edgesNum != verticesNum * (verticesNum - 1) / 2 && (double) edgesNum / (verticesNum * (verticesNum - 1) / 2) > 0.005) {
			for (int n = 0; n < edgesNum; n++) {
				int cost = random(1000) + 1;
				int vNodeKey = random(verticesNum);
				int eNodeKey = random(verticesNum);
				while (eNodeKey == vNodeKey) {
					eNodeKey = random(verticesNum);
				}
				while (randomGraph.isDuplicateEdge(vNodeKey, eNodeKey)) {
					vNodeKey = random(verticesNum);
					eNodeKey = random(verticesNum);
					while (eNodeKey == vNodeKey) {
						eNodeKey = random(verticesNum);
					}
				}
				randomGraph.addEdge(vNodeKey, eNodeKey, cost);
			}
		} else if (edgesNum == verticesNum * (verticesNum - 1) / 2) {
			for (int vNodeKey = 0; vNodeKey < verticesNum; vNodeKey++) {
				for (int eNodeKey = vNodeKey + 1; eNodeKey < verticesNum; eNodeKey++) {
					int cost = random(1000) + 1;
					randomGraph.addEdge(vNodeKey, eNodeKey, cost);
				}
			}
		} else {
			for (int vNodeKey = 0; vNodeKey < verticesNum - 1; vNodeKey++) {
				int cost = random(1000) + 1;
				randomGraph.addEdge(vNodeKey, vNodeKey + 1, cost);
			}
			for (int n = verticesNum - 2; n < edgesNum; n++) {
				int cost = random(1000) + 1;
				int vNodeKey = random(verticesNum);
				int eNodeKey = random(verticesNum);
				while (eNodeKey == vNodeKey) {
					eNodeKey = random(verticesNum);
				}
				while (randomGraph.isDuplicateEdge(vNodeKey, eNodeKey)) {
					vNodeKey = random(verticesNum);
					eNodeKey = random(verticesNum);
					while (eNodeKey == vNodeKey) {
						eNodeKey = random(verticesNum);
					}
				}
				randomGraph.addEdge(vNodeKey, eNodeKey, cost);
			}
		}
		return randomGraph;
	}
	
	/*
	 * random(int) method is to generate integer within the upper limit.
	 * The parameter is the maximum number that the method can generate.
	 * It will return the random number.
	 */
	private int random(int max) {
		int randomNum = new Random().nextInt(max);
		return randomNum;
	}
	
	/*
	 * simpleScheme(Graph, int) method is to use dijikstra's method find the
	 * shortest distance from the source node to other nodes with simple data
	 * structure array.
	 * The three parameters are the graph, the source node and whether to
	 * output the distance or not.
	 */
	private void simpleScheme(Graph graph, int sourceKey, boolean output) {
		ArrayList<VertexInf> routeTable = new ArrayList<VertexInf>();
		ArrayList<VertexInf> vertexQueue = new ArrayList<VertexInf>();
		int verticesNum = graph.getVertexNum();
		for (int key = 0; key < verticesNum; key++) {
			int distance = (key == sourceKey) ? 0 : 9999999;
			VertexInf vertex = new VertexInf(key, distance);
			routeTable.add(vertex);
			vertexQueue.add(vertex);
		}
		
		while (!vertexQueue.isEmpty()) {
			VertexInf vi = vertexQueue.get(0);
			int minIndex = 0;
			for (int index = 1; index < vertexQueue.size(); index++) {
				if (vi.getDistance() > vertexQueue.get(index).getDistance()) {
					vi = vertexQueue.get(index);
					minIndex = index;
				}
			}
			int vNodeKey = vi.getKey();
			vertexQueue.remove(minIndex);
			VertexInf vNodeInf = routeTable.get(vNodeKey);
			vNodeInf.setKnown();
			ENode eNodePointer = graph.getVertex(vNodeKey).getNext();
			while (eNodePointer != null) {
				int eNodeKey = eNodePointer.getKey();
				int cost = eNodePointer.getCost();
				VertexInf eNodeInf = routeTable.get(eNodeKey);
				int vNodeDistance = vNodeInf.getDistance();
				int eNodeDistance = eNodeInf.getDistance();
				if (!eNodeInf.isKnown() && (vNodeDistance + cost < eNodeDistance)) {
					eNodeInf.setPreVertexKey(vNodeKey);
					eNodeInf.setDistance(vNodeDistance + cost);
				}
				eNodePointer = eNodePointer.getNext();
			}
		}
		
		if (output) {
			for (int m = 0; m < routeTable.size(); m++) {
				System.out.println(routeTable.get(m).getDistance() + " ");
			}
		}
	}
	
	/*
	 * fHeapScheme(Graph, int) method is to use dijikstra's method find the
	 * shortest distance from the source node to other nodes with Fibonacci
	 * Heap.
	 * The three parameters are the graph, the source node and whether to
	 * output the distance or not.
	 */
	private void fHeapScheme(Graph graph, int sourceKey, boolean output) {
		ArrayList<VertexInf> routeTable = new ArrayList<VertexInf>();
		FibonacciHeap<VertexInf> vertexQueue = new FibonacciHeap<VertexInf>();
		HashMap<VertexInf, FibonacciHeapNode<VertexInf>> vertexTable = new HashMap<VertexInf, FibonacciHeapNode<VertexInf>>();
		int verticesNum = graph.getVertexNum();
		for (int key = 0; key < verticesNum; key++) {
			int distance = (key == sourceKey) ? 0 : 9999999;
			VertexInf vertex = new VertexInf(key, distance);
			routeTable.add(vertex);
			vertexTable.put(vertex, vertexQueue.insert(vertex, vertex.getDistance()));
		}
		
		while (!vertexQueue.isEmpty()) {
			VertexInf vNode = vertexQueue.removeMin();
			
			int vNodeKey = vNode.getKey();
			VertexInf vNodeInf = routeTable.get(vNodeKey);
			vNodeInf.setKnown();
			ENode eNodePointer = graph.getVertex(vNodeKey).getNext();
			while (eNodePointer != null) {
				int eNodeKey = eNodePointer.getKey();
				int cost = eNodePointer.getCost();
				VertexInf eNodeInf = routeTable.get(eNodeKey);
				int vNodeDistance = vNodeInf.getDistance();
				int eNodeDistance = eNodeInf.getDistance();
				if (!eNodeInf.isKnown() && (vNodeDistance + cost < eNodeDistance)) {
					eNodeInf.setPreVertexKey(vNodeKey);
					eNodeInf.setDistance(vNodeDistance + cost);
					vertexQueue.decreaseKey(vertexTable.get(eNodeInf), vNodeDistance + cost);
				}
				eNodePointer = eNodePointer.getNext();
			}
		}
		
		if (output) {
			for (int m = 0; m < routeTable.size(); m++) {
				System.out.println(routeTable.get(m).getDistance() + " ");
			}
		}
	}
	
	/*
	 * Class GraphAndSource is used to describe the information of
	 * the user input graph.
	 */
	class GraphAndSource {
		Graph graph;
		int sourceKey;
		
		public GraphAndSource(Graph graph, int sourceKey) {
			this.graph = graph;
			this.sourceKey = sourceKey;
		}
	}
	
	/*
	 * readGraphFromFile(String) method is to read the file which contains
	 * the information of a graph.
	 * The parameter is the file path.
	 * It will return GraphAndSource object.
	 */
	private GraphAndSource readGraphFromFile(String filePath) throws IOException{
		InputStream is = new FileInputStream(filePath);
		String line;
		BufferedReader br = new BufferedReader(new InputStreamReader(is));
		line = br.readLine();
		int sourceKey = Integer.parseInt(line);
		line = br.readLine();
		int verticesNum = Integer.parseInt(line.substring(0, line.indexOf(" ")));
		int edgesNum = Integer.parseInt(line.substring(line.indexOf(" ") + 1, line.length()));
		
		Graph graphFromFile = new Graph(verticesNum, edgesNum);
		for (int key = 0; key < verticesNum; key++) {
			VNode vertex = new VNode(key);
			graphFromFile.vertices.add(vertex);
		}
		line = br.readLine();
		while (line != null) {
			int vNodeKey = Integer.parseInt(line.substring(0, line.indexOf(" ")));
			int eNodeKey = Integer.parseInt(line.substring(line.indexOf(" ") + 1, line.indexOf(" ", line.indexOf(" ") + 1)));
			int cost = Integer.parseInt(line.substring(line.indexOf(" ", line.indexOf(" ") + 1) + 1, line.length()));
			graphFromFile.addEdge(vNodeKey, eNodeKey, cost);
			line = br.readLine();
		}
		
		br.close();
		is.close();
		
		return new GraphAndSource(graphFromFile, sourceKey);
	}
	
	public static void main(String[] args) {
		dijkstra d = new dijkstra();
		if (args.length == 0) {
			System.out.println("Error: No arguments");
			System.exit(0);
		}
		if (args[0].equals("-r")) {
			int verticesNum = Integer.parseInt(args[1]);
			double density = Double.parseDouble(args[2]) / 100.00;
			int sourceKey = Integer.parseInt(args[3]);
			int edgeNum = (int) Math.ceil(density * ((verticesNum * (verticesNum - 1.00)) / 2.00));
			Graph randomGraph = d.randomGraphGenerate(verticesNum, edgeNum);
			while (density != 1.00 && !randomGraph.isConnected()) {
				randomGraph = d.randomGraphGenerate(verticesNum, edgeNum);
			}
//			randomGraph.printGraph();
			int simpleTime = 0, fHeapTime = 0;
			for (int count = 0; count < 10; count++) {
				long startTime = System.currentTimeMillis();
				d.simpleScheme(randomGraph, sourceKey, false);
				long endTime = System.currentTimeMillis();
				simpleTime += (endTime - startTime);
				startTime = System.currentTimeMillis();
				d.fHeapScheme(randomGraph, sourceKey, false);
				endTime = System.currentTimeMillis();
				fHeapTime += (endTime - startTime);
			}
			System.out.println("Simple Scheme Avg Runing Time: " + (simpleTime / 10) + "ms");
			System.out.println("FHeap Scheme Avg Runing Time: " + (fHeapTime / 10) + "ms");
//			long startTime = System.currentTimeMillis();
//			d.simpleScheme(randomGraph, sourceKey, false);
//			long endTime = System.currentTimeMillis();
//			System.out.println("Simple Scheme Runing Time: " + (endTime - startTime) + "ms");
//			startTime = System.currentTimeMillis();
//			d.fHeapScheme(randomGraph, sourceKey, false);
//			endTime = System.currentTimeMillis();
//			System.out.println("FHeap Scheme Runing Time: " + (endTime - startTime) + "ms");
		} else if (args[0].equals("-s")) {
			String filePath = args[1];
			try {
				GraphAndSource gs = d.readGraphFromFile(filePath);
				Graph graph = gs.graph;
				int sourceKey = gs.sourceKey;
				d.simpleScheme(graph, sourceKey, true);
			} catch (IOException e) {
				System.out.println("Can't read the file: " + filePath);
			}
			
		} else if (args[0].equals("-f")) {
			String filePath = args[1];
			try {
				GraphAndSource gs = d.readGraphFromFile(filePath);
				Graph graph = gs.graph;
				int sourceKey = gs.sourceKey;
				d.fHeapScheme(graph, sourceKey, true);
			} catch (IOException e) {
				System.out.println("Can't read the file: " + filePath);
			}
		} else {
			System.out.println("Error: The mode " + args[0] + " does not exist");
			System.exit(0);
		}
		
		
	}
}
