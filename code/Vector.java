// Vector.java

class Vector {
	protected final int N;
	protected final double[] values;
	protected static final double INF = Double.MAX_VALUE/2.0;

	public void set(int i, double value) {
		this.values[i] = value;
	}
	public double get(int i) {
		return this.values[i];
	}
	public int dimension() {
		return this.N;
	}

	public Vector scale(double scale) {
		final Vector P = new Vector(this.dimension());
		for(int i = 0; i < this.dimension(); i++)
			P.set(i, this.get(i)*scale);
		return P;
	}
	public double norm() {
		double tmp = 0.0;
		for(int i = 0; i < this.dimension(); i++)
			tmp += this.get(i)*this.get(i);
		return Math.sqrt(tmp);
	}
	public Vector normalize() {
		final double norm = this.norm();
		if(norm == 0.0) {
			System.err.println("OPERATION INFEASIBLE");
			return null;
		}
		return this.scale(1.0/norm);
	}
	public Matrix convertToMatrix() {
		final Matrix matrix = new Matrix(this.dimension(), 1);
		for(int i = 0; i < this.dimension(); i++)
			matrix.set(i, 0, this.get(i));
		return matrix;
	}
	public Vector copy() {
		return this.scale(1.0);
	}
	public void display() {
		for(int i = 0; i < this.dimension(); i++)
			System.out.println(this.get(i));
	}

	public static Vector add(final Vector v1, final Vector v2) {
		final int N = v1.dimension();
		if(N != v2.dimension()) {
			System.err.println("OPERATION INFEASIBLE");
			return null;
		}
		final Vector S = new Vector(N);
		for(int i = 0; i < N; i++)
			S.set(i, v1.get(i) + v2.get(i));
		return S;
	}
	public static Vector sub(final Vector v1, final Vector v2) {
		final int N = v1.dimension();
		if(N != v2.dimension()) {
			System.err.println("OPERATION INFEASIBLE");
			return null;
		}
		final Vector D = new Vector(N);
		for(int i = 0; i < N; i++)
			D.set(i, v1.get(i) - v2.get(i));
		return D;
	}
	public static double getInnerProduct(final Vector v1, final Vector v2) {
		final int N = v1.dimension();
		if(N != v2.dimension()) {
			System.err.println("OPERATION INFEASIBLE");
			return -INF;
		}
		double tmp = 0.0;
		for(int i = 0; i < N; i++)
			tmp += v1.get(i)*v2.get(i);
		return tmp;
	}
	public static Vector multiply(final Matrix matrix, final Vector vector) {
		if(matrix.number_of_cols() != vector.dimension()) {
			System.err.println("OPERATION INFEASIBLE");
			return null;
		}
		final Matrix convert = vector.convertToMatrix();
		final Matrix product = Matrix.multiply(matrix, convert);
		return product.convertToVector();
	}
	public static double distance(final Vector v1, final Vector v2) {
		return Vector.sub(v1, v2).norm();
	}
	public static boolean compare(final Vector v1, final Vector v2) {
		final int N = v1.dimension();
		if(N != v2.dimension())
			return false;
		for(int i = 0; i < N; i++)
			if(v1.get(i) != v2.get(i))
				return false;
		return true;
	}
	public static Vector[] orthonormalize(final Vector[] v) {
		final int N = v[0].dimension();
		final Vector[] u = new Vector[v.length];
		for(int i = 0; i < v.length; i++) {
			if(N != v[i].dimension()) {
				System.err.println("OPERATION INFEASIBLE");
				return null;
			}
			u[i] = v[i].copy();
			for(int j = 0; j < i; j++) {
				final double tmp1 = Vector.getInnerProduct(u[j], v[i]);
				final double tmp2 = Vector.getInnerProduct(u[j], u[j]);
				u[i] = Vector.sub(u[i], u[j].scale(tmp1/tmp2));
			}
		}
		for(int i = 0; i < u.length; i++)
			u[i] = u[i].normalize();
		return u;
	}

	protected Vector(int N) {
		this.N = N;
		this.values = new double[this.N];
	}
	protected Vector(double[] values) {
		this.N = values.length;
		this.values = new double[this.N];
		for(int i = 0; i < this.N; i++)
			this.values[i] = values[i];
	}
}
