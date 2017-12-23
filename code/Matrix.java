// Matrix.java

import java.util.SplittableRandom;

class Matrix {
	protected final int NUMBER_OF_ROWS;
	protected final int NUMBER_OF_COLS;
	protected double[][] VALUES;

	protected static final double INF = Double.MAX_VALUE;
	protected static final double DELTA = Math.pow(10.0, -8.0);

	public void set(int i, int j, double value) {
		this.VALUES[i][j] = value;
	}
	public double get(int i, int j) {
		return this.VALUES[i][j];
	}
	public int number_of_rows() {
		return this.NUMBER_OF_ROWS;
	}
	public int number_of_cols() {
		return this.NUMBER_OF_COLS;
	}

	public Matrix scale(double scale) {
		final Matrix P = new Matrix(this.NUMBER_OF_ROWS, this.NUMBER_OF_COLS);
		for(int i = 0; i < this.NUMBER_OF_ROWS; i++)
			for(int j = 0; j < this.NUMBER_OF_COLS; j++)
				P.set(i, j, scale*this.get(i, j));
		return P;
	}
	public Matrix transpose() {
		final Matrix T = new Matrix(this.NUMBER_OF_COLS, this.NUMBER_OF_ROWS);
		for(int i = 0; i < this.NUMBER_OF_ROWS; i++)
			for(int j = 0; j < this.NUMBER_OF_COLS; j++)
				T.set(j, i, this.get(i, j));
		return T;
	}
	public double getTrace() {
		if(this.NUMBER_OF_ROWS != this.NUMBER_OF_COLS) {
			System.out.println("NOT A SQUARE MATRIX");
			return INF;
		}
		double trace = 0.0;
		for(int i = 0; i < this.NUMBER_OF_ROWS; i++)
			trace += this.get(i, i);
		return trace;
	}
	public double getDeterminant() {
		final int N = this.NUMBER_OF_ROWS;
		if(N != this.NUMBER_OF_COLS) {
			System.out.println("NOT A SQUARE MATRIX");
			return INF;
		}
		final Matrix A = this.copy();
		boolean flag = true;
		double tmp;
		for(int i = 0; i < N; i++) {
			if(A.get(i, i) == 0.0)
				for(int j = i + 1; j < N; j++)
					if(A.get(j, i) != 0.0) {
						for(int k = 0; k < N; k++) {
							tmp = A.get(i, k);
							A.set(i, k, A.get(j, k));
							A.set(j, k, tmp);
						}
						flag = !flag;
						break;
					}
			if(A.get(i, i) == 0.0)
				return 0.0;
			for(int j = i + 1; j < N; j++) {
				tmp = A.get(j, i)/A.get(i, i);
				for(int k = 0; k < N; k++)
					A.set(j, k, A.get(j, k) - A.get(i, k)*tmp);
			}
		}
		tmp = 1.0;
		for(int i = 0; i < N; i++)
			tmp = tmp*A.get(i, i);
		if(flag)
			return tmp;
		return -tmp;
	}
	public Matrix inverse() {
		return this.rowReduction(identity(this.NUMBER_OF_ROWS));
	}
	public Matrix rowReduction(final Matrix matrix) {
		final int N = this.NUMBER_OF_ROWS;
		if(N != this.NUMBER_OF_COLS || N != matrix.NUMBER_OF_ROWS) {
			System.out.println("NOT A SQUARE MATRIX");
			return null;
		}
		final Matrix A = this.copy();
		final Matrix B = matrix.copy();
		double tmp;
		for(int i = 0; i < N; i++) {
			if(A.get(i, i) == 0.0)
				for(int j = i + 1; j < N; j++)
					if(A.get(j, i) != 0.0) {
						for(int k = 0; k < N; k++) {
							tmp = A.get(i, k);
							A.set(i, k, A.get(j, k));
							A.set(j, k, tmp);
						}
						for(int k = 0; k < B.NUMBER_OF_COLS; k++) {
							tmp = B.get(i, k);
							B.set(i, k, B.get(j, k));
							B.set(j, k, tmp);
						}
						break;
					}
			if(A.get(i, i) == 0.0) {
				System.out.println("IS THIS A REGULAR MATRIX?");
				return null;
			}
			tmp = 1.0/A.get(i, i);
			for(int j = 0; j < N; j++)
				A.set(i, j, A.get(i, j)*tmp);
			for(int j = 0; j < B.NUMBER_OF_COLS; j++)
				B.set(i, j, B.get(i, j)*tmp);
			for(int j = 0; j < N; j++) {
				if(i == j)
					continue;
				tmp = A.get(j, i);
				for(int k = 0; k < N; k++)
					A.set(j, k, A.get(j, k) - A.get(i, k)*tmp);
				for(int k = 0; k < B.NUMBER_OF_COLS; k++)
					B.set(j, k, B.get(j, k) - B.get(i, k)*tmp);
			}
		}
		return B;
	}
	public Matrix[] getLU() {
		final int N = this.NUMBER_OF_ROWS;
		if(N != this.NUMBER_OF_COLS) {
			System.out.println("NOT A SQUARE MATRIX");
			return null;
		}
		final Matrix A = this.copy();
		final Matrix L = new Matrix(N, N);
		final Matrix U = new Matrix(N, N);
		for(int i = 0; i < N; i++) {
			L.set(i, i, A.get(i, i));
			U.set(i, i, 1.0);
			for(int j = i + 1; j < N; j++) {
				L.set(j, i, A.get(j, i));
				U.set(i, j, A.get(i, j)/L.get(i, i));
			}
			for(int j = i + 1; j < N; j++)
				for(int k = i + 1; k < N; k++)
					A.set(j, k, A.get(j, k) - L.get(j, i)*U.get(i, k));
		}
		return (new Matrix[]{L, U});
	}
	public double[] getEigenvalues(int iteration) {
		if(this.NUMBER_OF_ROWS != this.NUMBER_OF_COLS) {
			System.out.println("NOT A SQUARE MATRIX");
			return null;
		}
		Matrix A = this.copy();
		for(int i = 0; i < iteration; i++) {
			final Matrix[] LU = A.getLU();
			A = multiply(LU[1], LU[0]);
		}
		final double[] E = new double[this.NUMBER_OF_ROWS];
		for(int i = 0; i < E.length; i++)
			E[i] = A.get(i, i);
		return E;
	}
	public double[] getEigenvalues() {
		return this.getEigenvalues(10000);
	}
	public Vector powerMethod() {
		if(this.NUMBER_OF_ROWS != this.NUMBER_OF_COLS) {
			System.out.println("NOT A SQUARE MATRIX");
			return null;
		}
		final int N = this.NUMBER_OF_ROWS;
		Vector x = new Vector(N);
		final SplittableRandom rnd = new SplittableRandom(System.currentTimeMillis());
		for(int i = 0; i < N; i++)
			x.set(i, rnd.nextDouble(0.0, 1.0));
		x = x.normalize();
		boolean convergence = false;
		while(!convergence) {
			final Vector tmp = Vector.multiply(this, x).normalize();
			convergence = (Vector.sub(x, tmp).norm() < DELTA);
			x = tmp;
		}
		return x;
	}
	public Vector inversePowerMethod() {
		final Matrix inverse = this.inverse();
		if(inverse == null)
			return null;
		return inverse.powerMethod();
	}
	public Vector inverseIteration(double value) {
		final int N = this.NUMBER_OF_ROWS;
		if(N != this.NUMBER_OF_COLS) {
			System.out.println("NOT A SQUARE MATRIX");
			return null;
		}
		final Matrix I = Matrix.identity(N);
		final Matrix m = Matrix.subtract(this, I.scale(value)).inverse();
		if(m == null)
			return null;
		return m.powerMethod();
	}
	public Vector getEigenvector() {
		return this.powerMethod();
	}
	public Vector[] getEigenvectors(int k) {
		if(this.NUMBER_OF_ROWS != this.NUMBER_OF_COLS) {
			System.out.println("NOT A SQUARE MATRIX");
			return null;
		}
		final int N = this.NUMBER_OF_ROWS;
		Vector[] x = new Vector[k];
		final SplittableRandom rnd = new SplittableRandom(System.currentTimeMillis());
		for(int i = 0; i < k; i++) {
			x[i] = new Vector(N);
			for(int j = 0; j < N; j++)
				x[i].set(j, rnd.nextDouble(0.0, 1.0));
			x[i] = x[i].normalize();
		}
		x = Vector.orthonormalize(x);
		boolean convergence = false;
		while(!convergence) {
			Vector[] tmp = new Vector[k];
			for(int i = 0; i < k; i++)
				tmp[i] = Vector.multiply(this, x[i]);
			tmp = Vector.orthonormalize(tmp);
			convergence = true;
			for(int i = 0; i < k; i++)
				if(Vector.sub(x[i], tmp[i]).norm() >= DELTA)
					convergence = false;
			x = tmp;
		}
		return x;
	}
	public Vector[] getEigenvectors() {
		return this.getEigenvectors(this.NUMBER_OF_ROWS);
	}
	public Vector convertToVector() {
		if(this.NUMBER_OF_COLS != 1) {
			System.out.println("CANNOT CONVERT TO VECTOR");
			return null;
		}
		final Vector vector = new Vector(this.NUMBER_OF_ROWS);
		for(int i = 0; i < this.NUMBER_OF_ROWS; i++)
			vector.set(i, this.get(i, 0));
		return vector;
	}
	public Matrix copy() {
		final Matrix C = new Matrix(this.NUMBER_OF_ROWS, this.NUMBER_OF_COLS);
		for(int i = 0; i < this.NUMBER_OF_ROWS; i++)
			for(int j = 0; j < this.NUMBER_OF_COLS; j++)
				C.set(i, j, this.get(i, j));
		return C;
	}
	public void display(final String DelimitSignal) {
		for(int i = 0; i < this.NUMBER_OF_ROWS; i++)
			for(int j = 0; j < this.NUMBER_OF_COLS; j++) {
				System.out.print((float)this.get(i, j));
				if(j < this.NUMBER_OF_COLS - 1)
					System.out.print(DelimitSignal);
				else
					System.out.print("\n");
			}
	}
	public void display() {
		this.display("\t");
	}

	public static Matrix identity(int n) {
		final Matrix I = new Matrix(n, n);
		for(int i = 0; i < n; i++)
			I.set(i, i, 1.0);
		return I;
	}
	public static Matrix diagonal(double[] values) {
		final int N = values.length;
		final Matrix D = new Matrix(N, N);
		for(int i = 0; i < N; i++)
			D.set(i, i, values[i]);
		return D;
	}
	public static Matrix add(final Matrix matrix1, final Matrix matrix2) {
		if(matrix1.NUMBER_OF_ROWS != matrix2.NUMBER_OF_ROWS ||
				matrix1.NUMBER_OF_COLS != matrix2.NUMBER_OF_COLS) {
			System.out.println("OPERATION INFEASIBLE");
			return null;
		}
		final Matrix S = new Matrix(matrix1.NUMBER_OF_ROWS, matrix1.NUMBER_OF_COLS);
		for(int i = 0; i < matrix1.NUMBER_OF_ROWS; i++)
			for(int j = 0; j < matrix1.NUMBER_OF_COLS; j++)
				S.set(i, j, matrix1.get(i, j) + matrix2.get(i, j));
		return S;
	}
	public static Matrix subtract(final Matrix matrix1, final Matrix matrix2) {
		if(matrix1.NUMBER_OF_ROWS != matrix2.NUMBER_OF_ROWS ||
				matrix1.NUMBER_OF_COLS != matrix2.NUMBER_OF_COLS) {
			System.out.println("OPERATION INFEASIBLE");
			return null;
		}
		final Matrix D = new Matrix(matrix1.NUMBER_OF_ROWS, matrix1.NUMBER_OF_COLS);
		for(int i = 0; i < matrix1.NUMBER_OF_ROWS; i++)
			for(int j = 0; j < matrix1.NUMBER_OF_COLS; j++)
				D.set(i, j, matrix1.get(i, j) - matrix2.get(i, j));
		return D;
	}
	public static Matrix multiply(final Matrix matrix1, final Matrix matrix2) {
		if(matrix1.NUMBER_OF_COLS != matrix2.NUMBER_OF_ROWS) {
			System.out.println("OPERATION INFEASIBLE");
			return null;
		}
		final Matrix P = new Matrix(matrix1.NUMBER_OF_ROWS, matrix2.NUMBER_OF_COLS);
		for(int i = 0; i < matrix1.NUMBER_OF_ROWS; i++)
			for(int j = 0; j < matrix2.NUMBER_OF_COLS; j++) {
				double sum = 0.0;
				for(int k = 0; k < matrix1.NUMBER_OF_COLS; k++)
					sum += matrix1.get(i, k)*matrix2.get(k, j);
				P.set(i, j, sum);
			}
		return P;
	}
	public static boolean compare(final Matrix matrix1, final Matrix matrix2) {
		if(matrix1.NUMBER_OF_ROWS != matrix2.NUMBER_OF_ROWS ||
				matrix1.NUMBER_OF_COLS != matrix2.NUMBER_OF_COLS)
			return false;
		for(int i = 0; i < matrix1.NUMBER_OF_ROWS; i++)
			for(int j = 0; j < matrix1.NUMBER_OF_COLS; j++)
				if(matrix1.get(i, j) != matrix2.get(i, j))
					return false;
		return true;
	}

	protected Matrix(double[][] values) {
		this.NUMBER_OF_ROWS = values.length;
		this.NUMBER_OF_COLS = values[0].length;
		this.VALUES = new double[this.NUMBER_OF_ROWS][this.NUMBER_OF_COLS];
		for(int i = 0; i < this.NUMBER_OF_ROWS; i++)
			for(int j = 0; j < this.NUMBER_OF_COLS; j++)
				this.VALUES[i][j] = values[i][j];
	}
	protected Matrix(int number_of_rows, int number_of_cols) {
		this.NUMBER_OF_ROWS = number_of_rows;
		this.NUMBER_OF_COLS = number_of_cols;
		this.VALUES = new double[this.NUMBER_OF_ROWS][this.NUMBER_OF_COLS];
	}
}
