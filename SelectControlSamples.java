import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

/**
 * Selects the best set of samples to use as control samples for SavvyCNV for calling CNVs.
 *
 * @author Matthew Wakeling
 */
public class SelectControlSamples
{
	public static void main(String[] args) throws Exception {
		// Two modes of operation:
		// 1. Create a summary of lots of samples.
		// 2. Select a set of controls for a set of cases from the summary.
		//
		// The selection should be done by finding the first N singular vectors, multiplying by the singular values, and finding the X samples with shortest distance to the sample.
		// SVD takes a long time with lots of samples, so we should create the summary by doing SVD on a random subset of the potential control samples, then finding the singular vectors for the rest of the samples using the dot product.
		// So, the summary has to contain:
		// 1. The first N singular vectors for genomic location, so that weights per-sample can be generated for new samples.
		// 2. The singular vectors for samples, created by dot product (which should be the same as the singular vector multiplied by the singular value), to compare the new samples to.
		
		List<String> samples = new ArrayList<String>();
		int divider = 1000000;
		int minReads = 20;
		int subsetCount = 50;
		String limitChromosome = null;
		String summaryFile = null;
		boolean cross = false;
		boolean stats = false;
		boolean svs = false;
		for (int i = 0; i < args.length; i++) {
			if ("-d".equals(args[i])) {
				i++;
				divider = Integer.parseInt(args[i]);
			} else if ("-minReads".equals(args[i])) {
				i++;
				minReads = Integer.parseInt(args[i]);
			} else if ("-subset".equals(args[i])) {
				i++;
				subsetCount = Integer.parseInt(args[i]);
			} else if ("-chr".equals(args[i])) {
				i++;
				limitChromosome = args[i];
			} else if ("-summary".equals(args[i])) {
				i++;
				summaryFile = args[i];
			} else if ("-cross".equals(args[i])) {
				cross = true;
			} else if ("-stats".equals(args[i])) {
				stats = true;
			} else if ("-svs".equals(args[i])) {
				svs = true;
			} else {
				samples.add(args[i]);
			}
		}
		if (summaryFile == null) {
			if (subsetCount > samples.size()) {
				subsetCount = samples.size();
			}
			System.err.println("Processing " + samples.size() + " samples - using subset of " + subsetCount + " samples to build SVD");
			System.err.println("Using divider of " + divider);
			System.err.println("Informative genome chunks have an average of " + minReads + " reads or more");
			// Mode 1 - create a summary. We have to select a random subset of samples, and run SVD on them.
			List<String> shuffledSamples = new ArrayList<String>(samples);
			Collections.shuffle(shuffledSamples);
			long[] totals = new long[subsetCount];
			Map<String, int[][]> arraysMap = new TreeMap<String, int[][]>();
			boolean includeDoubleClip = false;
			boolean includeSecondary = false;
			int minMq = -1;
			for (int i = 0; i < subsetCount; i++) {
				System.err.println("Reading " + i + " " + shuffledSamples.get(i));
				CoverageBinner in = new CoverageBinner(shuffledSamples.get(i));
				if (i == 0) {
					includeDoubleClip = in.getDoubleClipIncluded();
					includeSecondary = in.getSecondaryIncluded();
					minMq = in.getMinMq();
				} else if ((includeDoubleClip != in.getDoubleClipIncluded()) || (includeSecondary != in.getSecondaryIncluded()) || (minMq != in.getMinMq())) {
					System.err.println("Cannot analyse these samples because the CoverageBinner files were generated with different options so the data is not comparable.");
					System.err.println("File " + samples.get(0) + " has double-clipped reads " + (includeDoubleClip ? "included" : "excluded") + ", secondary/supplementary reads " + (includeSecondary ? "included" : "excluded") + ", and minimum MQ " + minMq);
					System.err.println("File " + samples.get(i) + " has double-clipped reads " + (in.getDoubleClipIncluded() ? "included" : "excluded") + ", secondary/supplementary reads " + (in.getSecondaryIncluded() ? "included" : "excluded") + ", and minimum MQ " + in.getMinMq());
					System.err.println("Note that the defaults options for CoverageBinner have changed - to create files compatible with the old version, please add the -includeDoubleClip option. Alternatively, for higher quality results, regenerate all your CoverageBinner files with the latest CoverageBinner with default options.");
					System.exit(1);
				}
				LinkedHashMap<String, int[]> vectors = in.getVectors(divider);
				for (Map.Entry<String, int[]> entry : vectors.entrySet()) {
					String chr = entry.getKey();
					if ((limitChromosome == null) || limitChromosome.equals(chr)) {
						int[][] array = arraysMap.get(chr);
						if (array == null) {
							array = new int[samples.size()][];
							arraysMap.put(chr, array);
						}
						array[i] = entry.getValue();
						for (int o = 0; o < array[i].length; o++) {
							totals[i] += array[i][o];
						}
					}
				}
			}
			System.err.println("Normalising");
			int totalPositions = 0;
			for (String chr : arraysMap.keySet()) {
				int[][] array = arraysMap.get(chr);
				int longestSample = 0;
				for (int i = 0; i < subsetCount; i++) {
					if (array[i] != null) {
						longestSample = (longestSample > array[i].length ? longestSample : array[i].length);
					} else {
						array[i] = new int[0];
					}
				}
				totalPositions += longestSample;
			}
			long total = 0L;
			for (int i = 0; i < subsetCount; i++) {
				total += totals[i];
			}
			System.err.println("Average reads in each bucket: " + ((total * 1.0) / subsetCount / totalPositions));
			// Work out which chunks of the genome we want to do SVD on.
			List<String> chunkChromosomes = new ArrayList<String>();
			List<Integer> chunkStarts = new ArrayList<Integer>();
			double[] scaleArray = new double[subsetCount];
			List<Long> totalReadsArray = new ArrayList<Long>();
			for (String chr : arraysMap.keySet()) {
				int[][] array = arraysMap.get(chr);
				boolean repeat = true;
				int start = 0;
				while (repeat) {
					repeat = false;
					long totalReads = 0;
					for (int i = 0; i < subsetCount; i++) {
						if (array[i].length > start) {
							totalReads += array[i][start];
							repeat = true;
						}
					}
					if ((totalReads * 50 * totalPositions > total) && (totalReads > subsetCount * minReads)) {
					//if (totalReads > 0) {
						chunkChromosomes.add(chr);
						chunkStarts.add(start);
						totalReadsArray.add(totalReads);
						for (int i = 0; i < subsetCount; i++) {
							if (array[i].length > start) {
								scaleArray[i] += (1.0 * array[i][start]) * subsetCount / totalReads;
							}
						}
					}
					start++;
				}
			}
			System.err.println("Number of genome chunks: " + chunkChromosomes.size());
			double scale = 0.0;
			for (int i = 0; i < subsetCount; i++) {
				scale += scaleArray[i];
			}
			for (int i = 0; i < subsetCount; i++) {
				scaleArray[i] = scaleArray[i] * subsetCount / scale;
			}
			double[][] aArray = new double[chunkChromosomes.size()][subsetCount];
			double[] sumArray = new double[chunkChromosomes.size()];
			for (int o = 0; o < chunkChromosomes.size(); o++) {
				String chr = chunkChromosomes.get(o);
				int start = chunkStarts.get(o);
				int[][] array = arraysMap.get(chr);
				double sum = 0.0;
				for (int i = 0; i < subsetCount; i++) {
					sum += (array[i].length > start ? array[i][start] / scaleArray[i] : 0.0);
				}
				sum = sum / subsetCount;
				sumArray[o] = sum;
				for (int i = 0; i < subsetCount; i++) {
					aArray[o][i] = Math.log((array[i].length > start ? array[i][start] / scaleArray[i] / sum : 0.0) + 0.01);
				}
			}
			// Don't need this data any more, so allow it to be GCed.
			arraysMap = null;
			Jama.Matrix A = new Jama.Matrix(aArray);
			if (subsetCount > chunkChromosomes.size()) {
				A = A.transpose();
			}
			long time1 = System.currentTimeMillis();
			Jama.SingularValueDecomposition decomp = new Jama.SingularValueDecomposition(A);
			long time2 = System.currentTimeMillis();
			System.err.println("Performed SVD in " + (time2 - time1) + "ms");

			time1 = System.currentTimeMillis();
			// Find inverse of U.
			Jama.SingularValueDecomposition invDecomp = new Jama.SingularValueDecomposition(decomp.getU());
			Jama.Matrix invS = invDecomp.getS();
			for (int i = 0; i < subsetCount; i++) {
				invS.getArray()[i][i] = 1.0 / invS.getArray()[i][i];
			}
			double[][] uInverse = invDecomp.getV().times(invS).times(invDecomp.getU().transpose()).getArray();
			//double[][] uInverse2 = decomp.getU().inverse().getArray();
			System.err.println("Performed inverse in " + (System.currentTimeMillis() - time1) + "ms");
			//for (int svNo = 0; svNo < subsetCount; svNo++) {
			//	double x = 0.0;
			//	for (int o = 0; o < uInverse[0].length; o++) {
			//		x += aArray[o][0] * uInverse[svNo][o];
			//	}
			//	System.out.println(svNo + "\t" + s[svNo][svNo] + "\t" + v[0][svNo] + "\t" + (s[svNo][svNo] * v[0][svNo]) + "\t" + x);
			//}
			// Now load each file in individually, and calculate SV for each one.
			if (cross) {
				double[][] vectors = new double[samples.size()][];
				for (int i = 0; i < samples.size(); i++) {
					vectors[i] = calculateSampleVector(subsetCount, chunkChromosomes, chunkStarts, totalReadsArray, scale, sumArray, uInverse, samples.get(i), limitChromosome, divider);
				}
				for (int i = 0; i < samples.size(); i++) {
					for (int o = 0; o < samples.size(); o++) {
						double sdistance = 0.0;
						for (int p = 0; p < subsetCount; p++) {
							double dist = vectors[i][p] - vectors[o][p];
							sdistance += dist * dist;
						}
						System.out.println(i + "\t" + o + "\t" + Math.sqrt(sdistance));
					}
				}
			} else if (svs) {
				for (int i = 0; i < samples.size(); i++) {
					double[] vector = calculateSampleVector(subsetCount, chunkChromosomes, chunkStarts, totalReadsArray, scale, sumArray, uInverse, samples.get(i), limitChromosome, divider);
					System.out.print(i + "\t" + samples.get(i));
					for (int svNo = 0; svNo < subsetCount; svNo++) {
						System.out.print("\t" + vector[svNo]);
					}
					System.out.println("");
				}
			} else {
				ObjectOutputStream out = new ObjectOutputStream(System.out);
				out.writeInt(subsetCount);
				out.writeObject(chunkChromosomes);
				out.writeObject(chunkStarts);
				out.writeObject(totalReadsArray);
				out.writeDouble(scale);
				out.writeObject(sumArray);
				out.writeObject(uInverse);
				out.writeObject(limitChromosome);
				out.writeInt(divider);
				out.writeObject(samples);
				for (int i = 0; i < samples.size(); i++) {
					double[] vector = calculateSampleVector(subsetCount, chunkChromosomes, chunkStarts, totalReadsArray, scale, sumArray, uInverse, samples.get(i), limitChromosome, divider);
//					System.out.print(i + "\t" + samples.get(i));
//					for (int svNo = 0; svNo < subsetCount; svNo++) {
//						System.out.print("\t" + vector[svNo]);
//					}
//					System.out.println("");
					out.writeObject(vector);
				}
				out.flush();
				out.close();
			}
		} else {
			ObjectInputStream in = new ObjectInputStream(new FileInputStream(summaryFile));
			int summarySubsetCount = in.readInt();
			@SuppressWarnings("unchecked") List<String> chunkChromosomes = (List<String>) in.readObject();
			@SuppressWarnings("unchecked") List<Integer> chunkStarts = (List<Integer>) in.readObject();
			@SuppressWarnings("unchecked") List<Long> totalReadsArray = (List<Long>) in.readObject();
			double scale = in.readDouble();
			double[] sumArray = (double[]) in.readObject();
			double[][] uInverse = (double[][]) in.readObject();
			limitChromosome = (String) in.readObject();
			divider = in.readInt();
			@SuppressWarnings("unchecked") List<String> summarySamples = (List<String>) in.readObject();
			double[][] vectors = new double[summarySamples.size()][];
			for (int i = 0; i < summarySamples.size(); i++) {
				vectors[i] = (double[]) in.readObject();
			}
			in.close();
			if (cross) {
				if (!samples.isEmpty()) {
					System.err.println("-cross was specified, ignoring samples specified in the command line");
				}
				for (int i = 0; i < summarySamples.size(); i++) {
					for (int o = 0; o < summarySamples.size(); o++) {
						double sdistance = 0.0;
						for (int p = 0; p < subsetCount; p++) {
							double dist = vectors[i][p] - vectors[o][p];
							sdistance += dist * dist;
						}
						System.out.println(i + "\t" + o + "\t" + Math.sqrt(sdistance));
					}
				}
			} else if (svs) {
				if (!samples.isEmpty()) {
					System.err.println("-svs was specified, ignoring samples specified in the command line");
				}
				for (int i = 0; i < summarySamples.size(); i++) {
					System.out.print(i + "\t" + summarySamples.get(i));
					for (int svNo = 0; svNo < subsetCount; svNo++) {
						System.out.print("\t" + vectors[i][svNo]);
					}
					System.out.println("");
				}
			} else {
				System.err.println("Processing " + samples.size() + " samples - finding " + subsetCount + " best matching samples out of pool of " + summarySamples.size());
				System.err.println("Using divider of " + divider);
				System.err.println("Informative genome chunks have an average of " + minReads + " reads or more");
				double[] sampleDistances = new double[summarySamples.size()];
				for (int i = 0; i < samples.size(); i++) {
					System.err.println("Reading " + i + " " + samples.get(i));
					double[] vector = calculateSampleVector(summarySubsetCount, chunkChromosomes, chunkStarts, totalReadsArray, scale, sumArray, uInverse, samples.get(i), limitChromosome, divider);
					for (int o = 0; o < summarySamples.size(); o++) {
						double sdistance = 0.0;
						for (int p = 0; p < summarySubsetCount; p++) {
							double dist = vector[p] - vectors[o][p];
							sdistance += dist * dist;
						}
						if (sdistance == 0.0) {
							// Sample is identical to the argument. Exclude it.
							sampleDistances[o] += Double.POSITIVE_INFINITY;
						} else {
							sampleDistances[o] += Math.sqrt(sdistance);
						}
					}
				}
				List<SortedSample> sortedSamples = new ArrayList<SortedSample>();
				for (int o = 0; o < summarySamples.size(); o++) {
					if ((!samples.contains(summarySamples.get(o))) && (sampleDistances[o] < Double.POSITIVE_INFINITY)) {
						sortedSamples.add(new SortedSample(summarySamples.get(o), sampleDistances[o]));
					}
				}
				Collections.sort(sortedSamples);
				int o = 0;
				int emitted = 0;
				while ((emitted < subsetCount) && (o < sortedSamples.size())) {
					SortedSample s = sortedSamples.get(o);
					if ((new File(s.getName())).canRead()) {
						if (stats) {
							System.out.println(s.getName() + "\t" + s.getDistance());
						} else {
							System.out.println(s.getName());
						}
						emitted++;
					} else {
						System.err.println("File " + s.getName() + " (score " + s.getDistance() + ") has gone missing");
					}
					o++;
				}
			}
		}
	}

	public static double[] calculateSampleVector(int subsetCount, List<String> chunkChromosomes, List<Integer> chunkStarts, List<Long> totalReadsArray, double scale, double[] sumArray, double[][] uInverse, String fileName, String limitChromosome, int divider) throws IOException, ClassNotFoundException {
		CoverageBinner in = new CoverageBinner(fileName);
		LinkedHashMap<String, int[]> vectors = in.getVectors(divider);
		double sampleScale = 0.0;
		for (int o = 0; o < chunkChromosomes.size(); o++) {
			String chr = chunkChromosomes.get(o);
			int start = chunkStarts.get(o);
			int[] array = vectors.get(chr);
			if ((array != null) && (array.length > start)) {
				sampleScale += (1.0 * array[start]) * subsetCount / totalReadsArray.get(o);
			}
		}
		sampleScale = sampleScale * subsetCount / scale;
		double[] sampleAArray = new double[chunkChromosomes.size()];
		for (int o = 0; o < chunkChromosomes.size(); o++) {
			String chr = chunkChromosomes.get(o);
			int start = chunkStarts.get(o);
			int[] array = vectors.get(chr);
			sampleAArray[o] = Math.log((((array != null) && (array.length > start)) ? array[start] / sampleScale / sumArray[o] : 0.0) + 0.01);
		}
		double[] retval = new double[subsetCount];
		for (int svNo = 0; svNo < subsetCount; svNo++) {
			double x = 0.0;
			for (int o = 0; o < uInverse[0].length; o++) {
				x += sampleAArray[o] * uInverse[svNo][o];
			}
			retval[svNo] = x;
		}
		return retval;
	}

	public static class SortedSample implements Comparable<SortedSample>
	{
		private String name;
		private double distance;

		public SortedSample(String name, double distance) {
			this.name = name;
			this.distance = distance;
		}

		public String getName() {
			return name;
		}

		public double getDistance() {
			return distance;
		}

		public int compareTo(SortedSample s) {
			if (s.distance > distance) {
				return -1;
			} else if (s.distance < distance) {
				return 1;
			} else {
				return name.compareTo(s.name);
			}
		}
	}
}
