/**
 * Implementation of preconditioning to allow faster optimisation of a multi-dimensional cost function.
 * This is a wrapper around an existing CostFunction. It transforms the coordinate system using the provided forward and inverse
 * matrices. The transformed coordinate system should have an improved condition number compared to the enclosed CostFunction.
 *
 * @author Matthew Wakeling
 */
public class PreconditionedCostFunction extends CostFunction
{
	private CostFunction cf;
	private double[][] forward, reverse;
	private double[] storedState;

	/**
	 * Creates the new object.
	 *
	 * @param cf the enclosed CostFunction
	 * @param forward the matrix used to convert from coordinates in this object into the coordinates used by the enclosing CostFunction
	 * @param reverse the matrix used to convert from coordinates in the enclosing CostFunction into the coordinates used by this object
	 */
	public PreconditionedCostFunction(CostFunction cf, double[][] forward, double[][] reverse) {
		this.cf = cf;
		this.forward = forward;
		this.reverse = reverse;
		this.storedState = multiply(reverse, cf.getState());
	}

	/**
	 * Creates the new object, calculating the forward and reverse matrices using the second derivative of the cost function.
	 *
	 * @param cf the enclosed CostFunction
	 */
	public PreconditionedCostFunction(CostFunction cf) {
		this.cf = cf;
		double[][] secondDerivative = cf.getSecondDerivative();
		//for (int i = 0; i < secondDerivative.length; i++) {
		//	System.err.println("SecondDerivative: " + CostFunction.printArray(secondDerivative[i]));
		//}
		Jama.Matrix A = new Jama.Matrix(secondDerivative);
		Jama.EigenvalueDecomposition decomp = new Jama.EigenvalueDecomposition(A);
		double[] eigenvalues = decomp.getRealEigenvalues();
		//double lowestEigenvalue = Double.MAX_VALUE;
		//double highestEigenvalue = -Double.MAX_VALUE;
		//for (int i = 0; i < eigenvalues.length; i++) {
		//	//System.err.println("Eigenvalue\t" + repeat + "\t" + sampleNo + "\t" + i + "\t" + eigenvalues[i]);
		//	lowestEigenvalue = Math.min(eigenvalues[i], lowestEigenvalue);
		//	highestEigenvalue = Math.max(eigenvalues[i], highestEigenvalue);
		//}
		//System.err.println("Condition number: " + (highestEigenvalue / lowestEigenvalue) + "\tIteration limit: " + (ITERATION_LIMITS[repeat] * contributions.length));
		double[][] inverseEigenArray = new double[secondDerivative.length][];
		//for (int i = 0; i < secondDerivative.length; i++) {
		//	inverseEigenArray[i] = new double[secondDerivative.length];
		//	inverseEigenArray[i][i] = 1.0 / eigenvalues[i];
		//}
		//Jama.Matrix invD = new Jama.Matrix(inverseEigenArray);
		Jama.Matrix V = decomp.getV();
		//Jama.Matrix invA = V.times(invD.times(V.transpose()));
		//double[][] invAArray = invA.getArray();
		//for (int i = 0; i < secondDerivative.length; i++) {
		//	for (int o = 0; o < secondDerivative.length; o++) {
		//		System.err.println("InvSecond\t" + repeat + "\t" + sampleNo + "\t" + i + "\t" + o + "\t" + invAArray[i][o]);
		//	}
		//}
		for (int i = 0; i < secondDerivative.length; i++) {
			inverseEigenArray[i] = new double[secondDerivative.length];
			inverseEigenArray[i][i] = 1.0 / Math.sqrt(eigenvalues[i]);
		}
		Jama.Matrix invsrD = new Jama.Matrix(inverseEigenArray);
		Jama.Matrix P = V.times(invsrD.times(V.transpose()));
		for (int i = 0; i < secondDerivative.length; i++) {
			inverseEigenArray[i] = new double[secondDerivative.length];
			inverseEigenArray[i][i] = Math.sqrt(eigenvalues[i]);
		}
		Jama.Matrix srD = new Jama.Matrix(inverseEigenArray);
		Jama.Matrix invP = V.times(srD.times(V.transpose()));
		boolean hasNaN = false;
		this.forward = P.getArray();
		for (int i = 0; i < secondDerivative.length; i++) {
			//System.err.println("Forward: " + CostFunction.printArray(this.forward[i]));
			for (int o = 0; o < secondDerivative.length; o++) {
				if (Double.isNaN(this.forward[i][o])) {
					hasNaN = true;
				}
			}
		}
		this.reverse = invP.getArray();
		for (int i = 0; i < secondDerivative.length; i++) {
			//System.err.println("Reverse: " + CostFunction.printArray(this.reverse[i]));
			for (int o = 0; o < secondDerivative.length; o++) {
				if (Double.isNaN(this.reverse[i][o])) {
					hasNaN = true;
				}
			}
		}
		if (hasNaN) {
			// Eigenvector decomposition has failed somehow. Fall back on the identity matrix.
			for (int i = 0; i < secondDerivative.length; i++) {
				for (int o = 0; o < secondDerivative.length; o++) {
					this.forward[i][o] = (i == o) ? 1.0 : 0.0;
					this.reverse[i][o] = (i == o) ? 1.0 : 0.0;
				}
			}
		}
		this.storedState = multiply(reverse, cf.getState());
	}

	private double[] multiply(double[][] matrix, double[] array) {
		double[] retval = new double[array.length];
		for (int i = 0; i < array.length; i++) {
			for (int o = 0; o < array.length; o++) {
				retval[i] += array[o] * matrix[i][o];
			}
		}
		return retval;
	}

	public double[] getState() {
		return storedState;
		//return multiply(reverse, cf.getState());
	}

	public double getCost() {
		return cf.getCost();
	}

	public double[] getGradient(double cost) {
		double[] grad = cf.getGradient(cost);
		return multiply(forward, grad);
	}

	public double[][] getSecondDerivative() {
		throw new RuntimeException("Unimplemented");
	}

	public void setState(double[] state) {
		double[] encState = multiply(forward, state);
		cf.setState(encState);
		storedState = state;
	}
}
