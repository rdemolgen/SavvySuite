import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.*;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Stack;
import java.util.regex.Pattern;

/**
 * Reads a VCF file, along with a sample name from the file and a contamination proportion, and tries to calculate which other sample in the file the contamination came from.
 *
 * @author Matthew Wakeling
 */
public class SavvyContaminationFinder
{
	public static final Pattern INDEL = Pattern.compile(".*[ACGT][ACGT].*");
	public static final double[] COST_LIMITS = new double[] {100.0, 100.0, 100.0, 5.0, 3.5, 3.2, 3.0, 10000.0};
	public static final int[] ITERATION_LIMITS = new int[] {1, 1, 2, 3, 3, 3, 3, 3};
	public static final double ERROR_RATE = 0.002;

	@SuppressWarnings("deprecation") public static void main(String[] args) throws Exception {
		String vcfFile = args[0];
		String sampleName = args[1];
		int maxDepth = (args.length > 2 ? Integer.parseInt(args[2]) : Integer.MAX_VALUE);
		int threadCount = (args.length > 3 ? Integer.parseInt(args[3]) : 1);
		VCFFileReader reader = new VCFFileReader(new File(vcfFile));
		List<String> sampleNames = reader.getFileHeader().getGenotypeSamples();
		int sampleCount = sampleNames.size();
		int minSampleNo = 0;
		int maxSampleNo = sampleCount - 1;
		if (!"all".equals(sampleName)) {
			for (int i = 0; i < sampleCount; i++) {
				if (sampleName.equals(sampleNames.get(i))) {
					minSampleNo = i;
					maxSampleNo = i;
				}
			}
		}
		List<List<Variant>> sampleVariants = new ArrayList<List<Variant>>();
		for (int i = 0; i < sampleCount; i++) {
			sampleVariants.add(new ArrayList<Variant>());
		}
		List<Variant> variants = new ArrayList<Variant>();
		boolean isFullVcf = true;
		for (VariantContext context : reader) {
			if ((!"X".equals(context.getChr())) && (!"Y".equals(context.getChr())) && (!context.getChr().startsWith("G")) && context.getFilters().isEmpty() && (!INDEL.matcher("" + context.getAlleles()).matches())) {
				boolean biallelic = true;
				for (int sampleNo = minSampleNo; sampleNo <= maxSampleNo; sampleNo++) {
					sampleName = sampleNames.get(sampleNo);
					Genotype g = context.getGenotype(sampleName);
					int[] ad = g.getAD();
					if ((ad == null) || (ad.length != 2)) {
						biallelic = false;
					}
				}
				if (biallelic) {
					int totalDepth = 0;
					for (int i = 0; i < sampleCount; i++) {
						Genotype g2 = context.getGenotype(sampleNames.get(i));
						int[] ad2 = g2.getAD();
						if (ad2 == null) {
							ad2 = new int[2];
							ad2[0] = 30;
							ad2[1] = 0;
						}
						totalDepth += ad2[0] + ad2[1];
					}
					if ((totalDepth >= 20 * sampleCount) && (totalDepth <= ((long) maxDepth) * sampleCount)) {
						int[] genotypes = new int[sampleCount];
						int[] refs = new int[sampleCount];
						int[] alts = new int[sampleCount];
						boolean[] types = new boolean[3];
						types[0] = false;
						types[1] = false;
						types[2] = false;
						for (int i = 0; i < sampleCount; i++) {
							Genotype g2 = context.getGenotype(sampleNames.get(i));
							int[] ad2 = g2.getAD();
							if (ad2 == null) {
								ad2 = new int[2];
								ad2[0] = 30;
								ad2[1] = 0;
							}
							double frac2 = (1.0 * ad2[1]) / (ad2[0] + ad2[1]);
							genotypes[i] = frac2 < 0.25 ? 0 : (frac2 < 0.75 ? 1 : 2);
							refs[i] = ad2[0];
							alts[i] = ad2[1];
							types[genotypes[i]] = true;
						}
						if (types[0] && types[1] && types[2]) {
							Variant variant = new Variant(genotypes, refs, alts);
							variants.add(variant);
							for (int sampleNo = minSampleNo; sampleNo <= maxSampleNo; sampleNo++) {
								sampleName = sampleNames.get(sampleNo);
								Genotype g = context.getGenotype(sampleName);
								int[] ad = g.getAD();
								double frac = (ad[1] * 1.0) / (ad[0] + ad[1]);
								if (frac >= 0.25) {
									sampleVariants.get(sampleNo).add(variant);
								}
							}
						}
					}
				}
			}
		}
		double[][] contributions = new double[sampleCount][];
		CostFunction[] costFunction = new CostFunction[sampleCount];
		int[] oldVariantCount = new int[sampleCount];
		for (int sampleNo = minSampleNo; sampleNo <= maxSampleNo; sampleNo++) {
			sampleName = sampleNames.get(sampleNo);
			System.err.println(sampleName + ": Found " + sampleVariants.get(sampleNo).size() + " informative variants");
			// Initial conditions for state - choose a state that will hopefully optimise down fairly quickly.
			contributions[sampleNo] = new double[sampleCount];
			for (int i = 0; i < sampleCount; i++) {
				contributions[sampleNo][i] = 0.000001;
			}
			contributions[sampleNo][sampleNo] = 1.0 - 0.000001 * (sampleCount - 1);
			oldVariantCount[sampleNo] = sampleVariants.get(sampleNo).size();
		}

		double[] lowestCosts = new double[sampleCount];
		for (int sampleNo = minSampleNo; sampleNo <= maxSampleNo; sampleNo++) {
			lowestCosts[sampleNo] = Double.MAX_VALUE;
		}
		// Fuzzy matching set to record what adjustments we have previously seen. We record the hash of the arrangement only.
		HashSet<Integer> previousAdjusted = new HashSet<Integer>();
		int[] hashCoefficients = new int[sampleCount * 7];
		java.util.Random rand = new java.util.Random();
		for (int i = 0; i < sampleCount * 7; i++) {
			hashCoefficients[i] = rand.nextInt();
		}
		int round = 0;
		for (int repeat = 0; (repeat < COST_LIMITS.length) && (round < 3000); repeat++) {
			round++;
			boolean adjustedHasChanged = false;
			int adjustedHash = 0;
			for (int sampleNo = minSampleNo; sampleNo <= maxSampleNo; sampleNo++) {
				sampleName = sampleNames.get(sampleNo);
				int highestContamSample = -1;
				double highestContam = -1.0;
				for (int i = 0; i < sampleCount; i++) {
					if (i != sampleNo) {
						if (contributions[sampleNo][i] > highestContam) {
							highestContam = contributions[sampleNo][i];
							highestContamSample = i;
						}
					}
				}
				System.err.println(sampleName + ": Most likely contaminant so far is " + (highestContamSample >= 0 ? sampleNames.get(highestContamSample) : "?") + " with contribution " + highestContam);

				//List<Double> costs = new ArrayList<Double>();
				int[] adjustedVariants = new int[7];
				if (repeat < COST_LIMITS.length - 1) {
					for (Variant v : variants) {
						int changeType = v.adjustGenotype(sampleNo, contributions[sampleNo], true);
						adjustedVariants[changeType]++;
					}
				
					Iterator<Variant> iter = sampleVariants.get(sampleNo).iterator();
					while (iter.hasNext()) {
						Variant v = iter.next();
						double vCost = v.getCost(contributions[sampleNo], sampleNo);
						//costs.add(vCost);
						if (vCost > COST_LIMITS[repeat]) {
							iter.remove();
							adjustedHasChanged = true;
						}
					}
				}
				//Collections.sort(costs);
				//for (int i = 0; i < costs.size(); i++) {
				//	System.err.println(i + "\t" + costs.get(i));
				//}
				int totalAdjusted = 0;
				for (int i = 1; i < 7; i++) {
					totalAdjusted += adjustedVariants[i];
					long hashPart = hashCoefficients[sampleNo * 7 + i] * adjustedVariants[i];
					adjustedHash = adjustedHash ^ ((int) (hashPart & 0xFFFFFFFF)) ^ ((int) ((hashPart >> 32) & 0xFFFFFFFF));
				}
				System.err.println(sampleName + ": Adjusted genotype for " + totalAdjusted + " variants and removed " + (oldVariantCount[sampleNo] - sampleVariants.get(sampleNo).size()) + " so far, now have " + sampleVariants.get(sampleNo).size() + " informative variants");
				System.err.println(sampleName + ": Adjustments:  " + adjustedVariants[1] + " 0/0->0/1   " + adjustedVariants[2] + " 0/0->1/1   " + adjustedVariants[3] + " 0/1->0/0   " + adjustedVariants[4] + " 0/1->1/1   " + adjustedVariants[5] + " 1/1->0/0   " + adjustedVariants[6] + " 1/1->0/1");

				Map<String, int[]> mergeVariantsMap = new HashMap<String, int[]>();
				Map<String, int[]> mergeVariantsPattern = new HashMap<String, int[]>();
				for (Variant v : sampleVariants.get(sampleNo)) {
					String pattern = v.getPattern();
					int[] depths = mergeVariantsMap.get(pattern);
					if (depths == null) {
						depths = new int[] {v.ref[sampleNo], v.alt[sampleNo]};
						mergeVariantsMap.put(pattern, depths);
						mergeVariantsPattern.put(pattern, v.samples);
					} else {
						depths[0] += v.ref[sampleNo];
						depths[1] += v.alt[sampleNo];
					}
				}
				List<Variant> mergedVariants = new ArrayList<Variant>();
				for (String pattern : mergeVariantsMap.keySet()) {
					int[] sampleGenotypes = mergeVariantsPattern.get(pattern);
					int[] sampleDepths = mergeVariantsMap.get(pattern);
					int[] refs = new int[sampleCount];
					int[] alts = new int[sampleCount];
					refs[sampleNo] = sampleDepths[0];
					alts[sampleNo] = sampleDepths[1];
					mergedVariants.add(new Variant(sampleGenotypes, refs, alts));
				}
				//System.err.println(sampleName + ": Merged " + sampleVariants.get(sampleNo).size() + " variants into " + mergedVariants.size());
				costFunction[sampleNo] = new Optimiser(contributions[sampleNo], mergedVariants, sampleNo);
				//System.err.println(sampleName + ": Cost: " + costFunction[sampleNo].getCost());
			}
			//System.err.println("Adjusted hash = " + adjustedHash + ", previous = " + previousAdjusted);
			adjustedHasChanged = !previousAdjusted.contains(adjustedHash);
			previousAdjusted.add(adjustedHash);
			if (adjustedHasChanged) {
				for (int sampleNo = minSampleNo; sampleNo <= maxSampleNo; sampleNo++) {
					lowestCosts[sampleNo] = Double.MAX_VALUE;
				}
			}

			System.err.println("Stage " + repeat + ", round " + round + ", iteration limit: " + (ITERATION_LIMITS[repeat] * contributions.length));
			Stack<Job> jobs = new Stack<Job>();
			for (int sampleNo = minSampleNo; sampleNo <= maxSampleNo; sampleNo++) {
				sampleName = sampleNames.get(sampleNo);
				Job job = new Job(sampleName, costFunction[sampleNo], ITERATION_LIMITS[repeat] * contributions.length);
				jobs.push(job);
			}
			List<Runner> runners = new ArrayList<Runner>();
			for (int i = 0; i < threadCount; i++) {
				Runner runner = new Runner(jobs);
				runners.add(runner);
				Thread t = new Thread(runner);
				t.start();
			}
			int maxIterations = 0;
			for (Runner runner : runners) {
				runner.waitForFinish();
				maxIterations = Math.max(maxIterations, runner.iterations);
			}

			boolean shouldRepeat = false;
			for (int sampleNo = minSampleNo; sampleNo <= maxSampleNo; sampleNo++) {
				double cost = costFunction[sampleNo].getCost();
				if (cost < lowestCosts[sampleNo]) {
					shouldRepeat = true;
					lowestCosts[sampleNo] = cost;
				}
			}
			//System.err.println("maxIterations: " + maxIterations);
			if (((repeat == 1) && adjustedHasChanged) || ((repeat == ITERATION_LIMITS.length - 1) && (maxIterations * 2 > ITERATION_LIMITS[repeat] * contributions.length))) {
				// Repeat the first optimisation until we get a decent solution.
				repeat--;
			} else {
				for (int sampleNo = minSampleNo; sampleNo <= maxSampleNo; sampleNo++) {
					lowestCosts[sampleNo] = Double.MAX_VALUE;
				}
			}
		}

		for (int sampleNo = minSampleNo; sampleNo <= maxSampleNo; sampleNo++) {
			sampleName = sampleNames.get(sampleNo);
			for (int i = 0; i < sampleCount; i++) {
				String otherSampleName = sampleNames.get(i);
				System.out.println(sampleName + "\t" + otherSampleName + "\t" + contributions[sampleNo][i]);
			}
		}
	}

	public static void testGradient(CostFunction cf) {
		double[] state = cf.getState();
		double cost = cf.getCost();
		double[] gradient = cf.getGradient(cost);
		for (int i = 0; i < gradient.length; i++) {
			double[] gradient2 = new double[4];
			double finite = 0.01;
			for (int o = 0; o < 4; o++) {
				double[] state2 = new double[state.length];
				for (int p = 0; p < state.length; p++) {
					state2[p] = state[p];
				}
				state2[i] += finite;
				cf.setState(state2);
				gradient2[o] = (cf.getCost() - cost) / finite;
				finite = finite / 10.0;
			}
			System.err.println("Gradient test: " + i + "\t" + gradient[i] + "\t" + gradient2[0] + "\t" + gradient2[1] + "\t" + gradient2[2] + "\t" + gradient2[3]);
		}
		cf.setState(state);
	}

	public static String formatArray(double[] array) {
		StringBuilder retval = new StringBuilder("[");
		for (int i = 0; i < array.length; i++) {
			if (i > 0) {
				retval.append(", ");
			}
			retval.append("" + array[i]);
		}
		retval.append("]");
		return retval.toString();
	}

	public static class Runner implements Runnable
	{
		Stack<Job> jobs;
		boolean finished = false;
		RuntimeException e = null;
		int iterations = 0;

		public Runner(Stack<Job> jobs) {
			this.jobs = jobs;
		}

		public void run() {
			while (true) {
				Job nextJob = null;
				synchronized(jobs) {
					if (jobs.isEmpty()) {
						synchronized(this) {
							finished = true;
							notifyAll();
						}
						return;
					}
					nextJob = jobs.pop();
				}
				try {
					nextJob.run();
					iterations = Math.max(iterations, nextJob.iterations);
				} catch (Exception ex) {
					synchronized(this) {
						e = new RuntimeException("Exception while running job", ex);
						finished = true;
						notifyAll();
					}
					return;
				}
			}
		}

		public synchronized void waitForFinish() {
			while (!finished) {
				try {
					wait();
				} catch (InterruptedException e) {
				}
			}
			if (e != null) {
				throw e;
			}
		}
	}

	public static class Job
	{
		String sampleName;
		CostFunction costFunction;
		int iterationLimit;
		int iterations;

		public Job(String sampleName, CostFunction costFunction, int iterationLimit) {
			this.sampleName = sampleName;
			this.costFunction = costFunction;
			this.iterationLimit = iterationLimit;
		}

		public void run() {
			long startTime = System.currentTimeMillis();
			PreconditionedCostFunction pcf = new PreconditionedCostFunction(costFunction);
			iterations = pcf.optimise(iterationLimit);
			double cost = pcf.getCost();
			long timeTaken = System.currentTimeMillis() - startTime;
			System.err.println(sampleName + ": Optimised down to cost " + cost + " in " + (timeTaken / 1000) + "." + ((timeTaken / 100) % 10) + ((timeTaken / 10) % 10) + (timeTaken % 10) + "s (" + iterations + " iterations)");
		}
	}

	public static class Optimiser extends CostFunction
	{
		private double[] contributions;
		private List<Variant> variants;
		private double[] state;
		int sampleNo;

		public Optimiser(double[] contributions, List<Variant> variants, int sampleNo) {
			this.contributions = contributions;
			this.variants = variants;
			this.sampleNo = sampleNo;
			this.state = new double[contributions.length - 1];
			for (int i = 0; i < state.length; i++) {
				int s = i >= sampleNo ? i + 1 : i;
				state[i] = contributions[s];
			}
		}

		public double[] getState() {
			//System.err.println("getState() = " + formatArray(state));
			double[] retval = new double[state.length];
			for (int i = 0; i < state.length; i++) {
				retval[i] = state[i];
			}
			return retval;
		}

		public double getCost() {
			//long start = System.currentTimeMillis();
			double cost = 0.0;
			for (Variant v : variants) {
				cost += v.getCost(contributions, sampleNo);
			}
			for (int i = 0; i < contributions.length; i++) {
				if (contributions[i] < 0.0) {
					cost += contributions[i] * contributions[i] * 10000000000.0;
				}
			}
			//System.err.println("cost(" + formatArray(contributions) + ") = " + cost);
			//if (Double.isNaN(cost)) {
			//	System.err.println("NaN cost with contributions " + CostFunction.printArray(contributions));
			//	cost = Double.MAX_VALUE;
			//}
			//System.err.println("cost() took " + (System.currentTimeMillis() - start) + "ms");
			return cost;
		}

		public double[] getGradient(double cost) {
			//long start = System.currentTimeMillis();
			double[] retval = new double[state.length];
			for (Variant v : variants) {
				double[] newContGrad = v.getGradient(contributions, sampleNo);
				for (int i = 0; i < state.length; i++) {
					int s = i >= sampleNo ? i + 1 : i;
					if (Double.isNaN(newContGrad[s])) {
						System.err.println("NaN while calculating gradient using variant " + v);
						System.err.println("Contributions: " + formatArray(contributions));
						System.err.println("newContGrad: " + formatArray(newContGrad));
						System.exit(1);
					}
					retval[i] += newContGrad[s] - newContGrad[sampleNo];
				}
			}
			double lastContGrad = contributions[sampleNo] < 0.0 ? contributions[sampleNo] * 20000000000.0 : 0.0;
			for (int i = 0; i < state.length; i++) {
				int s = i >= sampleNo ? i + 1 : i;
				if (contributions[s] < 0.0) {
					retval[i] += contributions[s] * 20000000000.0;
				}
				retval[i] -= lastContGrad;
			}
			//System.err.println("gradient(" + formatArray(contributions) + ") = " + formatArray(retval));
			//System.err.println("gradient() took " + (System.currentTimeMillis() - start) + "ms");
			return retval;
		}

		public double[][] getSecondDerivative() {
			double[][] sum = new double[contributions.length][];
			for (int i = 0; i < contributions.length; i++) {
				sum[i] = new double[contributions.length];
				if (contributions[i] < 0.0) {
					sum[i][i] = 20000000000.0;
				}
			}
			for (Variant v : variants) {
				double[][] vDer = v.getSecondDerivative(contributions, sampleNo);
				for (int i = 0; i < contributions.length; i++) {
					for (int o = 0; o < contributions.length; o++) {
						sum[i][o] += vDer[i][o];
					}
				}
			}
			double[][] retval = new double[state.length][];
			for (int i = 0; i < state.length; i++) {
				int si = i >= sampleNo ? i + 1 : i;
				retval[i] = new double[state.length];
				for (int o = 0; o < state.length; o++) {
					int so = o >= sampleNo ? o + 1 : o;
					retval[i][o] = sum[si][so] - sum[si][sampleNo] - sum[sampleNo][so] + sum[sampleNo][sampleNo];
				}
			}
			return retval;
		}

		public void setState(double[] state) {
			//System.err.println("setState(" + formatArray(state) + ")");
			//Exception e = new Exception("");
			//e.fillInStackTrace();
			//e.printStackTrace(System.err);
			for (int i = 0; i < this.state.length; i++) {
				this.state[i] = state[i];
			}
			contributions[sampleNo] = 1.0;
			for (int i = 0; i < this.state.length; i++) {
				int s = i >= sampleNo ? i + 1 : i;
				contributions[s] = this.state[i];
				contributions[sampleNo] -= this.state[i];
			}
			//System.err.println("setState(" + CostFunction.printArray(state) + ")");
			//System.err.println("Setting contributions to " + CostFunction.printArray(contributions));
			//System.err.println("setState() Cost: " + getCost());
		}
	}

	public static class Variant
	{
		private int[] samples, originalSamples, ref, alt;

		public Variant(int[] samples, int[] ref, int[] alt) {
			this.samples = samples;
			this.originalSamples = new int[samples.length];
			for (int i = 0; i < samples.length; i++) {
				this.originalSamples[i] = samples[i];
			}
			this.ref = ref;
			this.alt = alt;
		}

		public static final double PARABOLA_FRAC = 0.02;

		public double getCost(double[] contributions, int sampleNo) {
			double fraction = 0.0;
			for (int i = 0; i < samples.length; i++) {
				fraction += 0.5 * samples[i] * contributions[i];
			}
			double readFrac = (1.0 * alt[sampleNo]) / (ref[sampleNo] + alt[sampleNo]);
			return fractionCost(fraction, readFrac, ref[sampleNo], alt[sampleNo]);
		}

		public static double fractionCost(double fraction, double readFrac, int ref, int alt) {
			fraction = (1.0 - ERROR_RATE * 2.0) * fraction + ERROR_RATE;
			if (fraction < (readFrac == 0.0 ? -Double.MAX_VALUE : readFrac * PARABOLA_FRAC)) {
				return parabolaCost(fraction, readFrac * PARABOLA_FRAC, ref, alt);
			} else if (fraction > (readFrac == 1.0 ? Double.MAX_VALUE : 1.0 - (1.0 - readFrac) * PARABOLA_FRAC)) {
				return parabolaCost(fraction, 1.0 - (1.0 - readFrac) * PARABOLA_FRAC, ref, alt);
			} else {
				return -((alt > 0 ? alt * (Math.log(fraction) - Math.log((1.0 * alt) / (alt + ref))) : 0.0)
						+ (ref > 0 ? ref * (Math.log(1.0 - fraction) - Math.log((1.0 * ref) / (alt + ref))) : 0.0));
			}
		}

		public static double parabolaCost(double fraction, double switchFrac, int ref, int alt) {
			double readFrac = (1.0 * alt) / (ref + alt);
			double val = -((alt > 0 ? alt * (Math.log(switchFrac) - Math.log(readFrac)) : 0.0)
					+ (ref > 0 ? ref * (Math.log(1.0 - switchFrac) - Math.log(1.0 - readFrac)) : 0.0));
			double firstD = (ref > 0 ? ref / (1.0 - switchFrac) : 0.0) - (alt > 0 ? alt / switchFrac : 0.0);
			double secondD = (ref > 0 ? ref / (1.0 - switchFrac) / (1.0 - switchFrac) : 0.0) + (alt > 0 ? alt / switchFrac / switchFrac : 0.0);
			double a = secondD / 2.0;
			double b = firstD - secondD * switchFrac;
			double c = val - a * switchFrac * switchFrac - b * switchFrac;
			return a * fraction * fraction + b * fraction + c;
		}

		public double[] getGradient(double[] contributions, int sampleNo) {
			double fraction = 0.0;
			for (int i = 0; i < samples.length; i++) {
				fraction += 0.5 * samples[i] * contributions[i];
			}
			double readFrac = (1.0 * alt[sampleNo]) / (ref[sampleNo] + alt[sampleNo]);
			fraction = (1.0 - ERROR_RATE * 2.0) * fraction + ERROR_RATE;
			double baseGrad = ((ref[sampleNo] > 0 ? ref[sampleNo] / (1.0 - fraction) : 0.0) - (alt[sampleNo] > 0 ? alt[sampleNo] / fraction : 0.0)) * (1.0 - ERROR_RATE * 2.0);
			if (fraction < (readFrac == 0.0 ? -Double.MAX_VALUE : readFrac * PARABOLA_FRAC)) {
				baseGrad = parabolaGradient(fraction, readFrac * PARABOLA_FRAC, ref[sampleNo], alt[sampleNo]) * (1.0 - ERROR_RATE * 2.0);
			} else if (fraction > (readFrac == 1.0 ? Double.MAX_VALUE : 1.0 - (1.0 - readFrac) * PARABOLA_FRAC)) {
				baseGrad = parabolaGradient(fraction, 1.0 - (1.0 - readFrac) * PARABOLA_FRAC, ref[sampleNo], alt[sampleNo]) * (1.0 - ERROR_RATE * 2.0);
			}
			if (Double.isNaN(baseGrad)) {
				System.err.println("NaN gradient. reads = " + ref[sampleNo] + ":" + alt[sampleNo] + ", fraction = " + fraction);
			}
			double[] retval = new double[samples.length];
			for (int i = 0; i < samples.length; i++) {
				retval[i] = 0.5 * samples[i] * baseGrad;
			}
			return retval;
		}

		public static double parabolaGradient(double fraction, double switchFrac, int ref, int alt) {
			double readFrac = (1.0 * alt) / (ref + alt);
			double firstD = (ref > 0 ? ref / (1.0 - switchFrac) : 0.0) - (alt > 0 ? alt / switchFrac : 0.0);
			double secondD = (ref > 0 ? ref / (1.0 - switchFrac) / (1.0 - switchFrac) : 0.0) + (alt > 0 ? alt / switchFrac / switchFrac : 0.0);
			double b = firstD - secondD * switchFrac;
			return secondD * fraction + b;
		}

		public double[][] getSecondDerivative(double[] contributions, int sampleNo) {
			double fraction = 0.0;
			for (int i = 0; i < samples.length; i++) {
				fraction += 0.5 * samples[i] * contributions[i];
			}
			double readFrac = (1.0 * alt[sampleNo]) / (ref[sampleNo] + alt[sampleNo]);
			double valMult = (1.0 - ERROR_RATE * 2.0) * (1.0 - ERROR_RATE * 2.0);
			fraction = (1.0 - ERROR_RATE * 2.0) * fraction + ERROR_RATE;
			double baseVal = ((ref[sampleNo] > 0 ? ref[sampleNo] / (1.0 - fraction) / (1.0 - fraction) : 0.0) + (alt[sampleNo] > 0 ? alt[sampleNo] / fraction / fraction : 0.0)) * valMult;
			if (fraction < (readFrac == 0.0 ? -Double.MAX_VALUE : readFrac * PARABOLA_FRAC)) {
				baseVal = parabolaSecond(readFrac * PARABOLA_FRAC, ref[sampleNo], alt[sampleNo]) * valMult;
			} else if (fraction > (readFrac == 1.0 ? Double.MAX_VALUE : 1.0 - (1.0 - readFrac) * PARABOLA_FRAC)) {
				baseVal = parabolaSecond(1.0 - (1.0 - readFrac) * PARABOLA_FRAC, ref[sampleNo], alt[sampleNo]) * valMult;
			}
			double[][] retval = new double[samples.length][];
			for (int i = 0; i < samples.length; i++) {
				retval[i] = new double[samples.length];
				for (int o = 0; o < samples.length; o++) {
					retval[i][o] = 0.25 * samples[i] * samples[o] * baseVal;
				}
			}
			return retval;
		}

		public static double parabolaSecond(double switchFrac, int ref, int alt) {
			return (ref > 0 ? ref / (1.0 - switchFrac) / (1.0 - switchFrac) : 0.0) + (alt > 0 ? alt / switchFrac / switchFrac : 0.0);
		}

		public int adjustGenotype(int sampleNo, double[] contributions, boolean toHet) {
			// Adjust the genotype of the sample we are checking, because the genotype may have been called wrong because of the contamination.
			double fraction = 0.0;
			for (int i = 0; i < samples.length; i++) {
				if (i != sampleNo) {
					fraction += 0.5 * samples[i] * contributions[i];
				}
			}
			int newGenotype = -1;
			double bestCost = Double.MAX_VALUE;
			double readFrac = (1.0 * alt[sampleNo]) / (ref[sampleNo] + alt[sampleNo]);
			for (int i = 0; i <= 2; i++) {
				double varFraction = fraction + 0.5 * i * contributions[sampleNo];
				double newCost = fractionCost(varFraction, readFrac, ref[sampleNo], alt[sampleNo]);
				//System.err.println("Cost for sample " + sampleNo + " genotype " + i + " is " + newCost);
				if (newCost < bestCost) {
					newGenotype = i;
					bestCost = newCost;
				}
			}
			//System.err.println("Best genotype is " + newGenotype);
			double presentFraction = (1.0 * alt[sampleNo]) / (ref[sampleNo] + alt[sampleNo]);
			double fractionForSample = (presentFraction - fraction) / contributions[sampleNo];
			//int newGenotype = fractionForSample < 0.07 ? 0 : (fractionForSample <= 0.93 ? 1 : 2);
			if (toHet || (newGenotype != 1) || (originalSamples[sampleNo] == 1)) {
				//int oldGenotype = samples[sampleNo];
				//if (oldGenotype != newGenotype) {
					//double cost = getCost(contributions, sampleNo);
					samples[sampleNo] = newGenotype;
					//double newCost = getCost(contributions, sampleNo);
					//System.err.println("Adjusting genotype from " + oldGenotype + " (cost " + cost + ") to " + newGenotype + " (cost " + newCost + "). fraction = " + fraction + ", presentFraction = " + presentFraction + " (" + ref[sampleNo] + ":" + alt[sampleNo] + "), contribution = " + contributions[sampleNo] + ", fractionForSample = " + fractionForSample);
				//}
			}
			if (originalSamples[sampleNo] == samples[sampleNo]) {
				return 0;
			}
			if (originalSamples[sampleNo] == 0) {
				if (samples[sampleNo] == 1) {
					return 1;
				} else {
					return 2;
				}
			} else if (originalSamples[sampleNo] == 1) {
				if (samples[sampleNo] == 0) {
					return 3;
				} else {
					return 4;
				}
			} else {
				if (samples[sampleNo] == 0) {
					return 5;
				} else {
					return 6;
				}
			}
		}

		public String getPattern() {
			char[] characters = new char[samples.length];
			for (int i = 0; i < samples.length; i++) {
				characters[i] = (char) ('0' + samples[i]);
			}
			return new String(characters);
		}

		public String toString() {
			String retval = "[";
			for (int i = 0; i < samples.length; i++) {
				retval += samples[i] + "(" + ref[i] + ":" + alt[i] + ")";
			}
			return retval + "]";
		}

		public int[] getGenotypes() {
			return samples;
		}
	}
}
