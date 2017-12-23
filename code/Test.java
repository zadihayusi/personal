// Test.java

import java.util.SplittableRandom;

import javafx.application.Application;
import javafx.stage.Stage;

public class Test extends Application {
	@Override
	public void start(final Stage stage) throws Exception {
		final SplittableRandom rnd = new SplittableRandom(System.currentTimeMillis());
		final int N = 1000;
		final int D = 3;
		final double r = 0.010;
		final Vector[] data = new Vector[N];
		for(int i = 0; i < N; i++) {
			data[i] = new Vector(D);
			final double theta = i*3.0*Math.PI/N;
			data[i].set(0, (1.0 + r*rnd.nextDouble(-1.0, 1.0))*Math.cos(theta)*2.0*Math.PI/(theta + 2.0*Math.PI));
			data[i].set(1, (1.0 + r*rnd.nextDouble(-1.0, 1.0))*Math.sin(theta)*2.0*Math.PI/(theta + 2.0*Math.PI));
			for(int j = 2; j < D; j++)
				data[i].set(j, rnd.nextDouble());
		}
		DimensionalityReduction.isomap(data);
	}
}
