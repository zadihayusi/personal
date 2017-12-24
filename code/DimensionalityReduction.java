// DimensionalityReduction.java

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.LinkedList;

import javafx.embed.swing.SwingFXUtils;
import javafx.scene.Group;
import javafx.scene.Scene;
import javafx.scene.canvas.Canvas;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.image.WritableImage;
import javafx.scene.input.KeyCode;
import javafx.scene.input.KeyEvent;
import javafx.scene.paint.Color;
import javafx.stage.DirectoryChooser;
import javafx.stage.Stage;

import javax.imageio.ImageIO;

class DimensionalityReduction {
	protected static final int DEFAULT_NUM_OF_COMPONENTS = 2;
	protected static final double INF = Double.MAX_VALUE/2.0;
	protected static final int WIDTH = 1200;
	protected static final int HEIGHT = 800;

	protected static Vector[] readData(final String fname) throws IOException {
		final File file = new File(fname);
		final LinkedList<Vector> data = new LinkedList<Vector>();
		try {
			final FileReader fr = new FileReader(file);
			final BufferedReader br = new BufferedReader(fr);
			String line = br.readLine();
			final int D = line.split("\t").length;
			while(line != null) {
				final String[] split = line.split("\t");
				if(split.length != D) {
					line = br.readLine();
					continue;
				}
				final Vector v = new Vector(D);
				for(int i = 0; i < D; i++)
					v.set(i, Double.parseDouble(split[i]));
				data.add(v);
				line = br.readLine();
			}
			br.close();
			fr.close();
		} catch(IOException e) {
			System.err.println(e.toString());
			return null;
		}
		return data.toArray(new Vector[0]);
	}
	protected static Vector[] center(final Vector[] data) {
		final int N = data.length;
		final int D = data[0].dimension();
		for(int i = 0; i < N; i++)
			if(data[i].dimension() != D) {
				System.err.println("DATA INVALID");
				return null;
			}
		final Vector ave = average(data);
		final Vector[] center_data = new Vector[N];
		for(int i = 0; i < N; i++)
			center_data[i] = Vector.sub(data[i], ave);
		return center_data;
	}
	protected static Vector[] standardize(final Vector[] data) {
		final int N = data.length;
		final int D = data[0].dimension();
		for(int i = 0; i < N; i++)
			if(data[i].dimension() != D) {
				System.err.println("DATA INVALID");
				return null;
			}
		final Vector[] center_data = center(data);
		final Vector tmp = new Vector(D);
		for(int i = 0; i < N; i++)
			for(int j = 0; j < D; j++)
				tmp.set(j, tmp.get(j) + center_data[i].get(j)*center_data[i].get(j));
		final Vector var = tmp.scale(1.0/N);
		final Vector[] std_data = new Vector[N];
		for(int i = 0; i < N; i++) {
			std_data[i] = new Vector(D);
			for(int j = 0; j < D; j++)
				std_data[i].set(j, center_data[i].get(j)/Math.sqrt(var.get(j)));
		}
		return std_data;
	}
	protected static Vector average(final Vector[] data) {
		final int N = data.length;
		final int D = data[0].dimension();
		for(int i = 0; i < N; i++)
			if(data[i].dimension() != D) {
				System.err.println("DATA INVALID");
				return null;
			}
		Vector sum = new Vector(D);
		for(int i = 0; i < N; i++)
			sum = Vector.add(sum, data[i]);
		return sum.scale(1.0/N);
	}
	protected static Matrix covariance(final Vector[] data) {
		final int N = data.length;
		final int D = data[0].dimension();
		for(int i = 0; i < N; i++)
			if(data[i].dimension() != D) {
				System.err.println("DATA INVALID");
				return null;
			}
		final Vector ave = average(data);
		Matrix tmp = new Matrix(D, D);
		for(int i = 0; i < N; i++) {
			final Matrix dif = Vector.sub(data[i], ave).convertToMatrix();
			tmp = Matrix.add(tmp, Matrix.multiply(dif, dif.transpose()));
		}
		return tmp.scale(1.0/N);
	}
	protected static Matrix getEuclideanDistance(final Vector[] data) {
		final int N = data.length;
		final Matrix EuclideanDistance = new Matrix(N, N);
		for(int i = 0; i < N; i++)
			for(int j = i + 1; j < N; j++) {
				final double d = Vector.distance(data[i], data[j]);
				EuclideanDistance.set(i, j, d);
				EuclideanDistance.set(j, i, d);
			}
		return EuclideanDistance;
	}

	protected static void plot(final Vector[] data) throws IOException {
		final int N = data.length;
		final int D = 2;
		for(int i = 0; i < N; i++)
			if(data[i].dimension() != D)
				return;
		final Stage stage = new Stage();
		stage.setTitle("Dimensionnality Reduction");
		stage.centerOnScreen();
		stage.setResizable(false);
		final Group root = new Group();
		final Canvas canvas = new Canvas(WIDTH, HEIGHT);
		final GraphicsContext gc = canvas.getGraphicsContext2D();
		draw(data, gc);
		root.getChildren().add(canvas);
		final Scene scene = new Scene(root, WIDTH, HEIGHT);
		canvas.setOnKeyPressed((final KeyEvent event) -> {
			final WritableImage wimg = scene.snapshot(null);
			if(event.getCode() == KeyCode.ENTER) {
				final DirectoryChooser dc = new DirectoryChooser();
				final File directory = dc.showDialog(stage);
				final String fname = directory.getPath().toString() + "/Plot" +
					String.valueOf(System.currentTimeMillis()%10000) + ".png";
				final File figure = new File(fname);
				try {
					ImageIO.write(SwingFXUtils.fromFXImage(wimg, null), "png", figure);
				} catch(IOException e) {
					System.err.println(e.toString());
				}
			}
		});
		canvas.setFocusTraversable(true);
		stage.setScene(scene);
		stage.show();
	}
	protected static void draw(final Vector[] data, final GraphicsContext gc) {
		final int N = data.length;
		final int D = 2;
		for(int i = 0; i < N; i++)
			if(data[i].dimension() != D)
				return;
		double maxX = -INF;
		double minX = INF;
		double maxY = -INF;
		double minY = INF;
		for(int i = 0; i < N; i++) {
			if(data[i].get(0) > maxX)
				maxX = data[i].get(0);
			if(data[i].get(0) < minX)
				minX = data[i].get(0);
			if(data[i].get(1) > maxY)
				maxY = data[i].get(1);
			if(data[i].get(1) < minY)
				minY = data[i].get(1);
		}
		gc.setStroke(Color.BLACK);
		gc.setLineWidth(1.0);
		gc.strokeLine(WIDTH/10, HEIGHT/10, 9*WIDTH/10, HEIGHT/10);
		gc.strokeLine(9*WIDTH/10, HEIGHT/10, 9*WIDTH/10, 9*HEIGHT/10);
		gc.strokeLine(9*WIDTH/10, 9*HEIGHT/10, WIDTH/10, 9*HEIGHT/10);
		gc.strokeLine(WIDTH/10, 9*HEIGHT/10, WIDTH/10, HEIGHT/10);
		for(int i = 0; i < 10; i++) {
			gc.strokeText(String.valueOf((float)(minX + i*(maxX - minX)/10.0)),
				WIDTH/10 + i*8*WIDTH/100, 9*HEIGHT/10 + 10);
			gc.strokeText(String.valueOf((float)(minY + i*(maxY - minY)/10.0)),
				WIDTH/10 - 50, 9*HEIGHT/10 - i*8*HEIGHT/100);
		}
		final int radius = 4;
		gc.setFill(Color.BLUE);
		for(int i = 0; i < N; i++) {
			final double x = data[i].get(0);
			final double y = data[i].get(1);
			final int x_fix = (int)(WIDTH*(0.10 + 0.80*(x - minX)/(maxX - minX)));
			final int y_fix = (int)(HEIGHT*(0.90 - 0.80*(y - minY)/(maxY - minY)));
			gc.fillOval(x_fix - radius, y_fix - radius, 2*radius, 2*radius);
		}
	}

	public static Vector[] runPCA(final String fname, int num_of_components) throws IOException {
		return runPCA(readData(fname), num_of_components);
	}
	public static Vector[] runPCA(final String fname) throws IOException {
		return runPCA(fname, DEFAULT_NUM_OF_COMPONENTS);
	}
	public static Vector[] runPCA(final Vector[] data, int num_of_components) throws IOException {
		final int N = data.length;
		final int D = data[0].dimension();
		for(int i = 0; i < N; i++)
			if(data[i].dimension() != D) {
				System.err.println("DATA INVALID");
				return null;
			}
		final Vector[] std_data = standardize(data);
		final Matrix cov = covariance(std_data);
		final Vector[] eigenvectors = cov.getEigenvectors(num_of_components);
		final Matrix V = new Matrix(D, num_of_components);
		for(int i = 0; i < D; i++)
			for(int j = 0; j < num_of_components; j++)
				V.set(i, j, eigenvectors[j].get(i));
		final Vector[] PCA = new Vector[N];
		for(int i = 0; i < N; i++)
			PCA[i] = Vector.multiply(V.transpose(), std_data[i]);
		plot(PCA);
		return PCA;
	}
	public static Vector[] runPCA(final Vector[] data) throws IOException {
		return runPCA(data, DEFAULT_NUM_OF_COMPONENTS);
	}

	protected static Vector[] cMDS(final Matrix distance, int num_of_components) throws IOException {
		final int N = distance.number_of_rows();
		if(N != distance.number_of_cols()) {
			System.err.println("DATA INVALID");
			return null;
		}
		final Matrix e = new Matrix(N, 1);
		for(int i = 0; i < N; i++)
			e.set(i, 0, 1.0);
		final Matrix I = Matrix.identity(N);
		final Matrix P = Matrix.subtract(I, Matrix.multiply(e, e.transpose()).scale(1.0/N));
		final Matrix G = Matrix.multiply(P, Matrix.multiply(distance, P)).scale(-0.50);
		final Vector[] eigenvectors = G.getEigenvectors(num_of_components);
		final Matrix sqrt_lambda = new Matrix(num_of_components, num_of_components);
		for(int i = 0; i < num_of_components; i++) {
			double lambda = Vector.multiply(G, eigenvectors[i]).get(0)/eigenvectors[i].get(0);
			sqrt_lambda.set(i, i, Math.sqrt(lambda));
		}
		final Matrix V = new Matrix(N, num_of_components);
		for(int i = 0; i < num_of_components; i++)
			for(int j = 0; j < N; j++)
				V.set(j, i, eigenvectors[i].get(j));
		final Matrix X = Matrix.multiply(sqrt_lambda, V.transpose());
		final Vector[] MDS = new Vector[N];
		for(int i = 0; i < N; i++) {
			MDS[i] = new Vector(num_of_components);
			for(int j = 0; j < num_of_components; j++)
				MDS[i].set(j, X.get(j, i));
		}
		plot(MDS);
		return MDS;
	}
	public static Vector[] runMDS(final String fname, int num_of_components) throws IOException {
		return runMDS(readData(fname), num_of_components);
	}
	public static Vector[] runMDS(final String fname) throws IOException {
		return runMDS(fname, DEFAULT_NUM_OF_COMPONENTS);
	}
	public static Vector[] runMDS(final Vector[] data, int num_of_components) throws IOException {
		final int N = data.length;
		final int D = data[0].dimension();
		for(int i = 0; i < N; i++)
			if(data[i].dimension() != D) {
				System.err.println("DATA INVALID");
				return null;
			}
		final Vector[] std_data = standardize(data);
		return cMDS(getEuclideanDistance(std_data), num_of_components);
	}
	public static Vector[] runMDS(final Vector[] data) throws IOException {
		return runMDS(data, DEFAULT_NUM_OF_COMPONENTS);
	}

	public static Vector[] isomap(final String fname, int num_of_components) throws IOException {
		return isomap(readData(fname), num_of_components);
	}
	public static Vector[] isomap(final String fname) throws IOException {
		return isomap(fname, DEFAULT_NUM_OF_COMPONENTS);
	}
	public static Vector[] isomap(final Vector[] data, int num_of_components) throws IOException {
		final int N = data.length;
		final int D = data[0].dimension();
		for(int i = 0; i < N; i++)
			if(data[i].dimension() != D) {
				System.err.println("DATA INVALID");
				return null;
			}
		final Vector[] std_data = standardize(data);
		final Matrix EuclideanDistance = getEuclideanDistance(std_data);
		final DirectedGraph k_neighbor = new DirectedGraph(N);
		final int k = 20;
		for(int i = 0; i < N; i++) {
			final boolean[] U = new boolean[N];
			for(int j = 0; j < N; j++)
				U[i] = true;
			U[i] = false;
			for(int j = 0; j < k; j++) {
				double min = INF + 1.0;
				int argmin = -1;
				for(int l = 0; l < N; l++)
					if(U[l] && EuclideanDistance.get(i, l) < min) {
						min = EuclideanDistance.get(i, l);
						argmin = l;
					}
				if(argmin >= 0) {
					k_neighbor.setEdgeCost(i, argmin, min);
					k_neighbor.setEdgeCost(argmin, i, min);
					U[argmin] = false;
				}
			}
		}
		Matrix GeodesicDistance = new Matrix(k_neighbor.floyd());
		boolean isConnected = true;
		for(int i = 0; i < N; i++)
			for(int j = 0; j < N; j++)
				if(GeodesicDistance.get(i, j) >= DirectedGraph.infinity())
					isConnected = false;
		while(!isConnected) {
			for(int i = 0; i < N; i++) {
				double min = INF + 1.0;
				int argmin1 = -1;
				int argmin2 = -1;
				for(int j = 0; j < N; j++)
					for(int l = j + 1; l < N; l++)
						if(GeodesicDistance.get(j, l) >= DirectedGraph.infinity() &&
								EuclideanDistance.get(j, l) < min) {
							min = EuclideanDistance.get(j, l);
							argmin1 = j;
							argmin2 = l;
						}
				if(argmin1 >= 0 && argmin2 >= 0) {
					k_neighbor.setEdgeCost(argmin1, argmin2, min);
					k_neighbor.setEdgeCost(argmin2, argmin1, min);
					GeodesicDistance.set(argmin1, argmin2, 0.0);
				}
			}
			GeodesicDistance = new Matrix(k_neighbor.floyd());
			isConnected = true;
			for(int i = 0; i < N; i++)
				for(int j = 0; j < N; j++)
					if(GeodesicDistance.get(i, j) >= DirectedGraph.infinity())
						isConnected = false;
		}
		return cMDS(GeodesicDistance, num_of_components);
	}
	public static Vector[] isomap(final Vector[] data) throws IOException {
		return isomap(data, DEFAULT_NUM_OF_COMPONENTS);
	}
}
