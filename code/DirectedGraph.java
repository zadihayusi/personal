// DirectedGraph.java

class DirectedGraph {
	protected final int N;
	protected final double[][] Cost;
	protected static final double INF = Double.MAX_VALUE/5.0;

	public static double infinity() {
		return INF;
	}

	public int getNumOfVertex() {
		return this.N;
	}
	public int getNumOfEdge() {
		int counter = 0;
		for(int i = 0; i < this.N; i++)
			for(int j = 0; j < this.N; j++)
				if(i != j && this.getEdgeCost(i, j) < INF)
					counter++;
		return counter;
	}
	public double getEdgeCost(int i, int j) {
		if(i < 0 || i >= this.N || j < 0 || j >= this.N)
			return INF;
		return this.Cost[i][j];
	}
	public void setEdgeCost(int i, int j, double cost) {
		if(i < 0 || i >= this.N || j < 0 || j >= this.N || i == j || cost < 0.0)
			return;
		this.Cost[i][j] = cost;
	}
	public void removeEdge(int i, int j) {
		this.setEdgeCost(i, j, INF);
	}

	public double[] dijkstra(int start) {
		final double[] C = new double[this.N];
		for(int i = 0; i < this.N; i++)
			C[i] = this.getEdgeCost(start, i);
		final boolean[] U = new boolean[this.N];
		for(int i = 0; i < this.N; i++)
			U[i] = true;
		for(int i = 0; i < this.N; i++) {
			double min = INF + 1.0;
			int argmin = 0;
			for(int j = 0; j < this.N; j++)
				if(U[j] && C[j] < min) {
					min = C[j];
					argmin = j;
				}
			U[argmin] = false;
			for(int j = 0; j < this.N; j++)
				if(U[j])
					C[j] = Math.min(C[j], C[argmin] + this.getEdgeCost(argmin, j));
		}
		System.out.println("Dijkstra algorithm finished.");
		return C;
	}
	public double[][] floyd() {
		final double[][] C = new double[this.N][this.N];
		for(int i = 0; i < this.N; i++)
			for(int j = 0; j < this.N; j++)
				C[i][j] = this.getEdgeCost(i, j);
		for(int i = 0; i < this.N; i++)
			for(int j = 0; j < this.N; j++)
				for(int k = 0; k < this.N; k++)
					C[j][k] = Math.min(C[j][k], C[j][i] + C[i][k]);
		System.out.println("Floyd algorithm finished.");
		return C;
	}

	protected DirectedGraph(int num_of_vertex) {
		this.N = num_of_vertex;
		this.Cost = new double[this.N][this.N];
		for(int i = 0; i < this.N; i++)
			for(int j = 0; j < this.N; j++)
				if(i == j)
					this.Cost[i][j] = 0.0;
				else
					this.Cost[i][j] = INF;
	}
	protected DirectedGraph(final double[][] cost) {
		this.N = cost.length;
		this.Cost = new double[this.N][this.N];
		for(int i = 0; i < this.N; i++)
			for(int j = 0; j < this.N; j++)
				if(i == j)
					this.Cost[i][j] = 0.0;
				else if(cost[i][j] < 0.0 || cost[i][j] > INF)
					this.Cost[i][j] = INF;
				else
					this.Cost[i][j] = cost[i][j];
	}
}
