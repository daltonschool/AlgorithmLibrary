import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;
import java.util.Stack;

public class Algorithms {
	
	final String FILE_PATH = "";
	
	public Algorithms() {
//		long start = System.currentTimeMillis();
//		
//		System.out.println(stirling(11, 5));	
//		
//		long finish = System.currentTimeMillis();
//		System.out.println(finish-start);
		
		//file reading:
//		List<String> data = new ArrayList<String>();
//		try {
//			BufferedReader br = new BufferedReader(new FileReader("FILE_PATH"));
//			String line;
//			while((line = br.readLine()) != null) {
//				data.add(line);
//			}
//		} catch (Exception e) {
//			e.printStackTrace();
//		}
		
//		p(permutations("ABCD"));
//		AllPaths.test();
		
//		DataWrapper d1 = new DataWrapper("Hi", 4);
//		DataWrapper d2 = new DataWrapper("Small guy", 1);
//		DataWrapper d3 = new DataWrapper("Big guy", 100);
//		
//		DataWrapper[] d = new DataWrapper[] {d1, d2, d3};
//		Arrays.sort(d, new Comparator<DataWrapper>() {
//			public int compare(DataWrapper a, DataWrapper b) {
//				return a.sortBy - b.sortBy;
//			}
//		});
//		
//		for(DataWrapper dw : d) {
//			System.out.println(dw.val);
//		}
		
		for(int i : primeFactors(111777)) {
			System.out.println(i);
		}
		
	}
	
	public ArrayList<String> permutations(String word) {
		ArrayList<String> res = new ArrayList<String>();
		permute("", word, res);
		return res;
	}
	
	public void permute(String prefix, String suffix, ArrayList<String> found) {
		if(suffix.length() == 0) {
			found.add(prefix);
			return;
		}
		
		for(int i = 0; i < suffix.length(); i++) {
			permute(prefix + suffix.charAt(i), suffix.substring(0, i) + suffix.substring(i+1), found);
		}
	}

	public boolean[] primeSieve(int size) {
		boolean[] isPrime = new boolean[size+1];
		
		// initially assume all integers >= 2 are prime
        for (int i = 2; i <= size; i++) {
            isPrime[i] = true;
        }
		
        // mark non-prime positions <= size as false
        for (int i = 2; i*i <= size; i++) {

            // if i is prime, then mark multiples of i as non-prime
            // suffices to consider multiples i, i+1, ..., N/i
            if (isPrime[i]) {
                for (int j = i; i*j <= size; j++) {
                    isPrime[i*j] = false;
                }
            }
        }
        
		return isPrime;
	}
	
	public List<Integer> primeFactors(int n) {
		List<Integer> factors = new ArrayList<Integer>();
		
		for(int i = 2; i*i <= n; i++) {
			while(n%i == 0) {
				factors.add(i);
				n /= i;
			}
		}
		
		if(n>1) factors.add(n);
		
		return factors;
	}
	
	public int gcd(int a, int b) {
		if(a==0) return b;
		return gcd(b%a, a);
	}
	
	/**
	 * Stirling numbers of the second kind denote the number of ways to
	 * partition a set of n objects into k non-empty subsets.
	 * 
	 * This implementation is very inefficient, so if you need to use
	 * large stirling numbers, use the iterative definition instead.
	 * 
	 * Or you could memoize it.
	 */
	public int stirling(int n, int k) {
		if(n==0 && k==0) return 1;
		else if(n==0 || k == 0) return 0;
		return k*stirling(n-1, k) + stirling(n-1, k-1);
		
		//iterative:
		//stirling(n, k) = 1/(k!) * sum from j=0 to k of (-1^j * kCj * (k-j)^n)
	}
	
//	public int[][] getAdjacentSpaces(int[][] someArray, int row, int col) { //does not include diagonal spaces
//		
//	}
	
	public int[][] getTouchingSpaces(int[][] someArray, int row, int col) { //includes diagonally touching spaces
		int[][] res = new int[8][2]; //max potential surrounding chars is 8, each char has 2 coords
		int counter = 0;
		
		//avoid IndexOutOfBounds:
		int[] loopVals = new int[] {-1, 1, -1, 1}; // { startingI, finalI, startingJ, finalJ }
		if(row==0) loopVals[0] = 0;
		if(row==someArray.length-1) loopVals[1] = 0;
		if(col==0) loopVals[2] = 0;
		if(col==someArray[0].length-1) loopVals[3] = 0;

		//find all touching spaces:
		for(int i = loopVals[0]; i <= loopVals[1]; i++) {
			for(int j = loopVals[2]; j <= loopVals[3]; j++) {
				if(i==0&&j==0) continue; //do not check the element itself
					res[counter][0] = row+i;
					res[counter][1]	= col+j;
					counter++;
			}
		}
		
		int[][] finalRes = new int[counter][2]; //this array will contain no null elements
		for(int i = 0; i < counter; i++) {
			finalRes[i][0] = res[i][0];
			finalRes[i][1] = res[i][1];
		}
		
		return finalRes;
	}
	
	public int sumDigits(int n) {
		String s = Integer.toString(n);
		int res = 0;
		for(int i = 0; i < s.length(); i++) {
			res += Integer.parseInt(Character.toString(s.charAt(i)));
		}
		return res;
	}
	
	public static int ctoi(char c)
	{
		return Integer.parseInt(Character.toString(c));
	}
	public static char itoc(int i)
	{
		return String.valueOf(i).charAt(0);
	}
	public static void p(int i)
	{
		System.out.println(i);
	}
	public static void p(double i)
	{
		System.out.println(i);
	}
	public static void p(char i)
	{
		System.out.println(i);
	}
	public static void p(String i)
	{
		System.out.println(i);
	}
	
	public static void p(int[] i)
	{
		for(int j = 0; j<i.length; j++)
		{
		System.out.println(i);
		}
	}
	public static void p(char[] i)
	{
		for(int j = 0; j<i.length; j++)
		{
		System.out.println(i);
		}
	}
	public static void p(double[] i)
	{
		for(int j = 0; j<i.length; j++)
		{
		System.out.println(i);
		}
	}
	public static void p(String[] i)
	{
		for(int j = 0; j<i.length; j++)
		{
		System.out.println(i);
		}
	}
	public static void p(int[][] x)
	{
		for(int i = 0; i<x.length; i++)
		{
			for(int j = 0; j<x[0].length; j++)
				System.out.print(x[i][j] + " ");
		}
		System.out.println();
	}
	
	
	public static void p(char[][] x)
	{
		for(int i = 0; i<x.length; i++)
		{
			for(int j = 0; j<x[0].length; j++)
				System.out.print(x[i][j] + " ");
		}
		System.out.println();
	}
	
	public static void p(double[][] x)
	{
		for(int i = 0; i<x.length; i++)
		{
			for(int j = 0; j<x[0].length; j++)
				System.out.print(x[i][j] + " ");
		}
		System.out.println();
	}
	
	public static void p(String[][] x)
	{
		for(int i = 0; i<x.length; i++)
		{
			for(int j = 0; j<x[0].length; j++)
				System.out.print(x[i][j] + " ");
		}
		System.out.println();
	}
	
	public static void p(ArrayList al) {
		for(Object item : al) {
			System.out.println(item);
		}
	}
	
	public static void main(String[] args) {
		new Algorithms();
	}
}

class Node {
	String val;
	
	public Node(String val) {
		this.val = val;
	}
	
	public boolean equals(Node other) {
		return this.hashCode() == other.hashCode();
	}
	
	public int hashCode() {
		return val.hashCode();
	}
}

class Graph {
	
	HashMap<Node, ArrayList<Node>> adj;
	
	public Graph() {
		adj = new HashMap<Node, ArrayList<Node>>();
	}
	
	public void addEdge(Node n1, Node n2) {
		if(adj.containsKey(n1))
			adj.get(n1).add(n2);
		else {
			ArrayList<Node> al = new ArrayList<Node>();
			al.add(n2);
			adj.put(n1, al);
		}
		
		if(adj.containsKey(n2))
			adj.get(n2).add(n1);
		else {
			ArrayList<Node> al = new ArrayList<Node>();
			al.add(n1);
			adj.put(n2, al);
		}
	}
	
	public ArrayList<Node> adjacentTo(Node n) {
		if(adj.get(n) == null) {
			adj.put(n, new ArrayList<Node>());
		}
		return adj.get(n);
	}
	
	public void print() {
		for(Entry<Node, ArrayList<Node>> entry : adj.entrySet()) {
			System.out.print(entry.getKey().val + ": ");
			for(Node n : entry.getValue()) {
				System.out.print(n.val + " ");
			}
			System.out.println();
		}
	}
}

class AllPaths {

    private Stack<Node> path  = new Stack<Node>();   	 // the current path
    private Set<Node> onPath  = new HashSet<Node>();     // the set of vertices on the path
    public ArrayList<ArrayList<Node>> foundPaths = new ArrayList<ArrayList<Node>>();

    public AllPaths(Graph G, Node s, Node t) { //find all paths from S to T
        enumerate(G, s, t);
    }

    // use DFS
    private void enumerate(Graph G, Node v, Node t) {

        // add node v to current path from s
        path.push(v);
        onPath.add(v);

        // found path from s to t - currently prints in reverse order because of stack        	
        if (v.equals(t)) {
        	ArrayList<Node> pathCopy = new ArrayList<Node>();
            for(Node n : path)
            	pathCopy.add(n);
             
        	foundPaths.add(pathCopy);
        }

        // consider all neighbors that would continue path with repeating a node
        else {
            for (Node w : G.adjacentTo(v)) {
                if (!onPath.contains(w)) enumerate(G, w, t);
            }
        }

        // done exploring from v, so remove from path
        path.pop();
        onPath.remove(v);
    }

    public static void test() {
        Graph graph = new Graph();
        
        Node A = new Node("A");
        Node B = new Node("B");
        Node C = new Node("C");
        Node D = new Node("D");
        Node E = new Node("E");
        Node F = new Node("F");
        Node G = new Node("G");
        
        graph.addEdge(A, B);
        graph.addEdge(A, C);
        graph.addEdge(C, D);
        graph.addEdge(D, E);
        graph.addEdge(C, F);
        graph.addEdge(B, F);
        graph.addEdge(F, D);
        graph.addEdge(D, G);
        graph.addEdge(E, G);
        graph.print();
        
        AllPaths allpaths = new AllPaths(graph, A, D);
    	System.out.println("Hi");
        for(ArrayList<Node> path : allpaths.foundPaths) {
        	for(Node n : path) {
        		System.out.print(n.val + " ");
        	}
        	System.out.println();
        }
    }

}

class Line {
	public double x1;
	public double x2;
	public double y1;
	public double y2;
	public double slope;
	public double yIntercept;
	
	public Line(double x1, double y1, double x2, double y2) {
		this.x1 = x1;
		this.y1 = y1;
		this.x2 = x2;
		this.y2 = y2;
		
		slope = (y2 - y1) / (x2 - x1);
		yIntercept = y1 - slope * x1;
	}
	
	public boolean intersects(Line other) {
		return true;
	}
	
	public double[] intersectPoint(Line other) {
		assert this.slope != other.slope;
		
		double x = (other.yIntercept - this.yIntercept) / (this.slope - other.slope);
		double y = slope * x + yIntercept;
		return new double [] {x, y};
	}
}

class DataWrapper {
	public String val;
	public int sortBy;
	
	public DataWrapper(String val, int sortBy) {
		this.val = val;
		this.sortBy = sortBy;
	}
}

//class Graph {
//	public ArrayList<GraphNode> nodes;
//	
//	public Graph(ArrayList<GraphNode> nodes) {
//		for(GraphNode n : nodes) {
//			this.nodes.add(n);
//		}
//	}
//	
//	public boolean areConnected(GraphNode n1, GraphNode n2) {
//		boolean res = n1.checkConnectedDFS(n2);
//		this.cleanGraph(); //make sure everything has visited == false again
//		return res;
//	}
//	
//	public ArrayList<GraphNode> getPathWithDFS(GraphNode from, GraphNode to) {
//		ArrayList<GraphNode> path = new ArrayList<GraphNode>();
//		from.findPath(path, to);
//		return path;
//	}
//	
//	public void cleanGraph() {
//		for(GraphNode n : nodes) {
//			n.visited = false;
//		}
//	}
//}
//
//class GraphNode {
//	public boolean visited;
//	public int id;
//	public ArrayList<GraphNode> links;
//	
//	public GraphNode(int id) {
//		this.visited = false;
//		this.id = id;
//		this.links = new ArrayList<GraphNode>();
//	}
//	
//	public boolean checkConnectedDFS(GraphNode other) {
//		if(this.id == other.id)
//			return true;
//		else {
//			this.visited = true;
//			boolean res = false;
//			for(GraphNode n : links) {
//				if(!n.visited) res = res || n.checkConnectedDFS(other);
//			}
//			return res;
//		}
//	}
//	
//	public ArrayList<GraphNode> findPath(ArrayList<GraphNode> pathSoFar, GraphNode to) {
//		this.visited = true;
//		pathSoFar.add(this);
//		if(this.id == to.id) {
//			return pathSoFar;
//		}
//		else {
//			for(GraphNode n : links) {
//				if(n.visited) continue;
//				
//			}
//		}
//	}
//}
