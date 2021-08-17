import java.util.Formatter;

/**
 * Implementation of a cost function to allow optimisation of a multi-dimensional cost function. Default optimiser uses conjugate gradients.
 *
 * @author Matthew Wakeling
 */
public abstract class CostFunction
{
	public abstract double[] getState();
	public abstract double getCost();
	public abstract double[] getGradient(double cost);
	public abstract double[][] getSecondDerivative();
	public abstract void setState(double[] state);

	public static final double SMALL_STEP = 0.000000000000001;
	public static final double GOLD2 = (Math.sqrt(5.0) - 1.0) / 2.0;
	public static final double GOLD1 = 1.0 - GOLD2;
	public static final double GOLD_LIMIT = 1.05;

	public int optimise(int maxIters) {
		double[] state = getState();
		//if (state.length > 2) {
		//	System.err.println("Optimising " + state.length + " values");
		//}
		double step = 0.1;
		int iterations = 0;
		double[] previousDirection = new double[state.length];
		double biggestStep = 0.0;
		int lastResetIteration = 0;
		while (iterations++ < maxIters) {
			//System.err.println("State: " + printArray(state));
			//System.err.println("Step: " + step);
			double cost = getCost();
			double origCost = cost;
			//System.err.println("Cost at beginning of iteration: " + origCost);
			//setState(state);
			//System.err.println("Cost again: " + getCost());
			if (((iterations % 10000) == 0) && (state.length > 2)) {
			//if (state.length > 2) {
				System.err.println(iterations + "\t" + cost + "\t" + printArray(state));
			}
			// Every 30 iterations, reset the conjugate gradients.
			if (iterations - lastResetIteration >= 30) {
				//System.err.println("Resetting conjugate gradients after 30 iterations");
				for (int i = 0; i < previousDirection.length; i++) {
					previousDirection[i] = 0.0;
				}
				biggestStep = 0.0;
				lastResetIteration = iterations;
			}
			double[] gradient = getGradient(cost);
			//System.err.println("Grad:\t" + printArray(gradient));
			double gradientLength = 0.0;
			double beta = 0.0;
			double betaDiv = 0.0;
			for (int i = 0; i < state.length; i++) {
				gradient[i] = -gradient[i];
				beta += gradient[i] * (gradient[i] - previousDirection[i]);
				betaDiv += previousDirection[i] * previousDirection[i];
			}
			beta = betaDiv == 0.0 ? 0.0 : (beta / betaDiv);
			//if (state.length > 2) {
			//	System.err.println("Beta: " + beta + ", div: " + betaDiv);
			//}
			//beta = 0.0;
			for (int i = 0; i < state.length; i++) {
				gradient[i] = gradient[i] + beta * previousDirection[i];
				gradientLength += gradient[i] * gradient[i];
			}
			gradientLength = Math.sqrt(gradientLength);
			if (gradientLength <= 0.0) {
				//if (state.length > 2) {
					//System.err.println("Zero gradient. Used " + iterations + " iterations, cost: " + cost);
				//}
				return iterations;
			}
			//System.err.println("Cost: " + cost);
			for (int i = 0; i < state.length; i++) {
				previousDirection[i] = gradient[i];
				gradient[i] = gradient[i] / gradientLength;
			}
			//System.err.println("Grad:\t" + printArray(gradient));
			// Now do a line search in the direction of negative gradient.
			// Algorithm:
			// Start with step from previous search. If the new cost is NaN, then we need to reduce the step.
			// If the new cost is less than old cost, then the step might be too small. Record the lowest cost seen, and increase the step until the new cost is greater than the lowest cost seen.
			// Once we have done these two, then the minimum must be bracketed by the old step and the new step. Do golden section search.
			double[] newState = new double[state.length];
			for (int i = 0; i < state.length; i++) {
				newState[i] = state[i] + step * gradient[i];
			}
			//System.err.println("Gradient: " + printArray(gradient));
			//System.out.println("New state (before line search): " + printArray(newState));
			setState(newState);
			double newCost = getCost();
			//System.err.println("New cost (before line search): " + newCost);
			//if ((newCost >= Double.MAX_VALUE) || Double.isNaN(newCost)) {
			//	System.err.println("Invalid cost with initial step " + step);
			//}
			if (Double.isNaN(newCost) || (newCost > cost)) {
				while ((Double.isNaN(newCost) || (newCost > cost)) && (step > SMALL_STEP)) {
					step = step / 4.0;
					//System.err.println("Decrease step: New step: " + step);
					for (int i = 0; i < state.length; i++) {
						newState[i] = state[i] + step * gradient[i];
					}
					//System.err.println("Decrease step: New state: " + printArray(newState));
					setState(newState);
					newCost = getCost();
					//System.err.println("Decrease step: New cost: " + newCost);
					//if ((newCost >= Double.MAX_VALUE) || Double.isNaN(newCost)) {
					//	System.err.println("Invalid cost while decreasing step to " + step);
					//}
					if (step < SMALL_STEP) {
						if (beta == 0.0) {
							// Terminate optimisation, because we have a small step size, after we had a small step size in the previous iteration and reset the conjugate gradients.
							if (newCost < cost) {
								cost = newCost;
							} else {
								step = 0.0;
								cost = origCost;
								setState(state);
							}
							//if (state.length > 2) {
								//System.err.println("Short line search. Used " + iterations + " iterations, step: " + step + ", cost: " + cost + ", cost: " + getCost());
								//System.err.println("State:\t" + printArray(state));
							//}
							return iterations;
						}
						// Reset the conjugate gradients and try once more. Likely we will terminate optimisation in the next round.
						//System.err.println("Resetting conjugate gradients due to short line search");
						for (int i = 0; i < previousDirection.length; i++) {
							previousDirection[i] = 0.0;
						}
						biggestStep = 0.0;
						lastResetIteration = iterations;
					}
				}
				step = step * 4.0;
				for (int i = 0; i < state.length; i++) {
					newState[i] = state[i] + step * gradient[i];
				}
				setState(newState);
				newCost = getCost();
			} else if (newCost < cost) {
				while ((newCost < cost) && (step < 10.0)) {
					cost = newCost;
					step = step * 4.0;
					//System.err.println("Increase step: New step: " + step);
					for (int i = 0; i < state.length; i++) {
						newState[i] = state[i] + step * gradient[i];
					}
					//System.err.println("Increase step: New state: " + printArray(newState));
					setState(newState);
					newCost = getCost();
					//System.err.println("Increase step: New cost: " + newCost);
				}
			}
			// Now start golden section method.
			// Hold four steps - a, b, c, and d. a and d bracket the minimum, while b and c are GOLD1 and GOLD2 along the line between a and d. If b < c, then the minimum must be between a and c, otherwise between b and d.
			double stepA = 0.0;
			double stepD = step;
			double stepB = stepA + GOLD1 * (stepD - stepA);
			double stepC = stepA + GOLD2 * (stepD - stepA);
			double costA = origCost;
			double costD = newCost;
			for (int i = 0; i < state.length; i++) {
				newState[i] = state[i] + stepB * gradient[i];
			}
			setState(newState);
			double costB = getCost();
			for (int i = 0; i < state.length; i++) {
				newState[i] = state[i] + stepC * gradient[i];
			}
			setState(newState);
			double costC = getCost();
			while (stepD > stepA * GOLD_LIMIT + SMALL_STEP) {
				//System.err.println("Gold\t" + stepA + "\t" + costA + "\t" + stepB + "\t" + costB + "\t" + stepC + "\t" + costC + "\t" + stepD + "\t" + costD);
				//if ((costA >= Double.MAX_VALUE) || Double.isNaN(costA) || (costB >= Double.MAX_VALUE) || Double.isNaN(costB) || (costC >= Double.MAX_VALUE) || Double.isNaN(costC) || (costD >= Double.MAX_VALUE) || Double.isNaN(costD)) {
				//	System.err.println("Invalid cost while doing golden section:\t" + stepA + "\t" + stepB + "\t" + stepC + "\t" + stepD + "\t" + costA + "\t" + costB + "\t" + costC + "\t" + costD);
				//}
				if ((costB <= costC) || Double.isNaN(costC)) {
					stepD = stepC;
					costD = costC;
					stepC = stepB;
					costC = costB;
					stepB = stepA + GOLD1 * (stepD - stepA);
					for (int i = 0; i < state.length; i++) {
						newState[i] = state[i] + stepB * gradient[i];
					}
					setState(newState);
					costB = getCost();
				} else {
					stepA = stepB;
					costA = costB;
					stepB = stepC;
					costB = costC;
					stepC = stepA + GOLD2 * (stepD - stepA);
					for (int i = 0; i < state.length; i++) {
						newState[i] = state[i] + stepC * gradient[i];
					}
					setState(newState);
					costC = getCost();
				}
			}
			// Now, we have narrowed down to a small enough region that it should be parabolic. So, we repeatedly do:
			// Find three points A, B, and C, where cost(B)<cost(A) and cost(B)<cost(C)
			// Find the parabola described by the quadratic equation y = ax^2 + bx + c, which goes through all three points.
			// Find the minimum of that, which is at -b/2a.
			// While the distance from the minimum to point B is greater than SMALL_STEP, set up a new set of three points, with the minimum as point B.
			if (costC < costB) {
				// The three points are B, C, and D. Shift them.
				stepA = stepB;
				stepB = stepC;
				stepC = stepD;
				costA = costB;
				costB = costC;
				costC = costD;
			}
/*			for (int o = 0; o < 2; o++) {
				System.err.println("Starting parabolic solution with " + stepA + "\t" + costA + "\t" + stepB + "\t" + costB + "\t" + stepC + "\t" + costC);
			//while (Math.max(stepB - stepA, stepC - stepB) > SMALL_STEP) {
				// Solve quadratic. cost(x) = ax^2 + bx + c
				// First point is at x = stepA, y = costA.
				// costA = a(stepA^2) + b(stepA) + c
				// costB = a(stepB^2) + b(stepB) + c
				// costA - costB = a(stepA^2 - stepB^2) + b(stepA - stepB)
				// costA - costC = a(stepA^2 - stepC^2) + b(stepA - stepC)
				// (costA - costB)(stepA - stepC) - (costA - costC)(stepA - stepB) = a((stepA^2 - stepB^2)(stepA - stepC) - (stepA^2 - stepC^2)(stepA - stepB))
				// -costAstepC - costBstepA + costBstepC + costAstepB + costCstepA - costCstepB = a(-stepA^2stepC - stepB^2stepA + stepB^2stepC + stepA^2stepB + stepC^2stepA - stepC^2stepB)
				double a = (-costA * stepC - costB + stepA + costB * stepC + costA * stepB + costC * stepA - costC * stepB) / (-stepA * stepA * stepC - stepB * stepB * stepA + stepB * stepB * stepC + stepA * stepA * stepB + stepC * stepC * stepA - stepC * stepC * stepB);
				double b = (costB - costA + a * (stepA * stepA - stepB * stepB)) / (stepB - stepA);
				//double c = costA - a * stepA * stepA - b * stepA;
				double minStep = -b / 2.0 / a;
				// minStep = ((costA - costB 
				for (int i = 0; i < state.length; i++) {
					newState[i] = state[i] + minStep * gradient[i];
				}
				setState(newState);
				double minCost = getCost();
				System.err.println("Parabolic solution at " + minStep + " cost " + minCost);
				if (minCost < costB) {
					if (minStep < stepB) {
						stepC = stepB;
						costC = costB;
					} else {
						stepA = stepB;
						costA = costB;
					}
					stepB = minStep;
					costB = minCost;
				} else {
					System.err.println("Not using parabolic solution.");
					for (int p = 0; p <= 20; p++) {
						double q = stepA + (stepC - stepA) * p / 20.0;
						for (int i = 0; i < state.length; i++) {
							newState[i] = state[i] + q * gradient[i];
						}
						setState(newState);
						double r = getCost();
						double[] gradient2 = getGradient(r);
						double dotGradient = 0.0;
						for (int i = 0; i < state.length; i++) {
							dotGradient += gradient2[i] * gradient[i];
						}
						System.err.println(q + "\t" + r + "\t" + dotGradient);
					}
					o++;
				}
			}
			//step = costB < costC ? stepB : stepC;
			//cost = costB < costC ? costB : costC;
*/
			// Use a secant method, using the A and B values. We are looking for the place where the gradient is 0. This assumes that we have an accurate gradient.
			// The best estimate will be held in stepB.
			double gradA = 0.0;
			for (int i = 0; i < state.length; i++) {
				newState[i] = state[i] + stepA * gradient[i];
			}
			setState(newState);
			double[] gradientA = getGradient(costA);
			double gradB = 0.0;
			for (int i = 0; i < state.length; i++) {
				newState[i] = state[i] + stepB * gradient[i];
			}
			setState(newState);
			costB = getCost();
			double[] gradientB = getGradient(costB);
			for (int i = 0; i < state.length; i++) {
				gradA += gradientA[i] * gradient[i];
				gradB += gradientB[i] * gradient[i];
			}
			for (int secantStep = 0; secantStep < 2; secantStep++) {
				stepC = stepB - gradB * (stepB - stepA) / (gradB - gradA);
				for (int i = 0; i < state.length; i++) {
					newState[i] = state[i] + stepC * gradient[i];
				}
				setState(newState);
				costC = getCost();
				if (costC < costB) {
					double[] gradientC = getGradient(costC);
					double gradC = 0.0;
					for (int i = 0; i < state.length; i++) {
						gradC += gradientC[i] * gradient[i];
					}
					//System.err.println("Used secant method from " + stepA + "\t" + costA + "\t" + gradA + "\t" + stepB + "\t" + costB + "\t" + gradB + "\tfound " + stepC);
					stepA = stepB;
					stepB = stepC;
					costA = costB;
					costB = costC;
					gradA = gradB;
					gradB = gradC;
				} else {
					for (int i = 0; i < state.length; i++) {
						newState[i] = state[i] + stepB * gradient[i];
					}
					setState(newState);
					secantStep += 100;
				}
			}
			step = stepB;
			cost = costB;
			//step = stepA;
			//cost = costA;
			for (int i = 0; i < state.length; i++) {
				newState[i] = state[i] + step * gradient[i];
			}
			setState(newState);
			//System.err.println("End of optimiser loop");
			state = newState;
			if ((cost >= Double.MAX_VALUE) || Double.isNaN(cost)) {
				System.err.println("Invalid cost with step " + step);
			}
			//System.err.println("Iteration " + iterations + "\tStep: " + step + "\tcost: " + cost);
			//System.err.println("State:\t" + printArray(state));
			if (step < SMALL_STEP) {
				if (beta == 0.0) {
					// Terminate optimisation, because we have a small step size, after we had a small step size in the previous iteration and reset the conjugate gradients.
					if (newCost < cost) {
						cost = newCost;
					} else {
						setState(state);
					}
					//if (state.length > 2) {
						//System.err.println("Short line search. Used " + iterations + " iterations, cost: " + newCost);
					//}
					return iterations;
				}
				// Reset the conjugate gradients and try once more. Likely we will terminate optimisation in the next round.
				//System.err.println("Resetting conjugate gradients due to small step");
				for (int i = 0; i < previousDirection.length; i++) {
					previousDirection[i] = 0.0;
				}
				biggestStep = 0.0;
				lastResetIteration = iterations;
			} else if (step < biggestStep * 0.03) {
				//System.err.println("Resetting conjugate gradients due to comparatively small step");
				for (int i = 0; i < previousDirection.length; i++) {
					previousDirection[i] = 0.0;
				}
				biggestStep = 0.0;
				lastResetIteration = iterations;
			} else {
				biggestStep = Math.max(step, biggestStep);
			}
		}
		//if (state.length > 2) {
		//	System.err.println("Too many iterations. Used " + iterations + " iterations, cost: " + getCost());
		//}
		return iterations;
	}

	public static String printArray(double[] array) {
		//StringBuilder retval = new StringBuilder("" + array[0]);
		StringBuilder retval = new StringBuilder();
		Formatter f = new Formatter(retval);
		f.format("% .4f", array[0]);
		for (int i = 1; i < array.length; i++) {
			f.format("\t% .4f", array[i]);
			//retval.append("\t" + array[i]);
		}
		return retval.toString();
	}



}
