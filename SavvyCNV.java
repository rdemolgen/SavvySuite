import java.io.BufferedOutputStream;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.ObjectInputStream;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

/**
 * Reads a set of CoverageBinner files, and uses a Hidden Markov Model to find copy number variants.
 * This will use the off-target reads.
 * This is useful for detecting large CNVs or trisomy.
 * The data is normalised using singular value decomposition.
 *
 * @author Matthew Wakeling
 */
public class SavvyCNV
{
	public static void main(String[] args) throws Exception {
		List<String> samples = new ArrayList<String>();
		Set<String> caseSamples = new HashSet<String>();
		int divider = 1000000;
		double cutoff = 0.25;
		double cutoffV = 5.25;
		boolean graph = false;
		boolean allGraphs = false;
		boolean writeData = false;
		String cytoBands = null;
		double transitionProb = 0.00001;
		double minProb = 0.00000000001;
		boolean dump = false;
		boolean errorModel = false;
		boolean intervals = false;
		boolean mosaic = false;
		String graphSize = "16in, 9in";
		double fontScale = 1.0;
		int svsBlanked = 5;
		int minReads = 20;
		int errorType = 0; // 0 means multiply, 1 means add, 2 means poisson.
		String limitChromosome = null;
		boolean sampleIsCase = true;
		for (int i = 0; i < args.length; i++) {
			if ("-d".equals(args[i])) {
				i++;
				divider = Integer.parseInt(args[i]);
			} else if ("-cutoff".equals(args[i])) {
				i++;
				cutoff = Double.parseDouble(args[i]);
			} else if ("-trans".equals(args[i])) {
				i++;
				transitionProb = Double.parseDouble(args[i]);
			} else if ("-g".equals(args[i])) {
				graph = true;
			} else if ("-a".equals(args[i])) {
				allGraphs = true;
				graph = true;
			} else if ("-cytoBands".equals(args[i])) {
				i++;
				cytoBands = args[i];
			} else if ("-minProb".equals(args[i])) {
				i++;
				minProb = exp(-Double.parseDouble(args[i]));
			} else if ("-dump".equals(args[i])) {
				dump = true;
			} else if ("-errorModel".equals(args[i])) {
				errorModel = true;
			} else if ("-intervals".equals(args[i])) {
				intervals = true;
			} else if ("-mosaic".equals(args[i])) {
				mosaic = true;
				System.err.println("Using mosaic mode");
			} else if ("-size".equals(args[i])) {
				i++;
				graphSize = args[i];
			} else if ("-fontscale".equals(args[i])) {
				i++;
				fontScale = Double.parseDouble(args[i]);
			} else if ("-sv".equals(args[i])) {
				i++;
				svsBlanked = Integer.parseInt(args[i]);
			} else if ("-minReads".equals(args[i])) {
				i++;
				minReads = Integer.parseInt(args[i]);
			} else if ("-addError".equals(args[i])) {
				errorType = 1;
			} else if ("-poisson".equals(args[i])) {
				errorType = 2;
			} else if ("-data".equals(args[i])) {
				writeData = true;
			} else if ("-chr".equals(args[i])) {
				i++;
				limitChromosome = args[i];
			} else if ("-case".equals(args[i])) {
				sampleIsCase = true;
			} else if ("-control".equals(args[i])) {
				sampleIsCase = false;
			} else {
				samples.add(args[i]);
				if (sampleIsCase) {
					caseSamples.add(args[i]);
				}
			}
		}
		double logTransProb = log(transitionProb);
		//System.err.println("Transition probability " + logTransProb);
		System.err.println("Processing " + samples.size() + " samples");
		System.err.println("Using divider of " + divider);
		System.err.println("Using noise cutoff of " + cutoff);
		System.err.println("Using transition probability of " + transitionProb + " (phred " + logTransProb + ")");
		System.err.println("Blanking " + svsBlanked + " singular vectors");
		System.err.println("Informative genome chunks have an average of " + minReads + " reads or more");
		if (errorType != 0) {
			System.err.println("Using " + (errorType == 1 ? "additive" : "poisson") + " error model");
		}
		if (svsBlanked >= samples.size()) {
			System.err.println(svsBlanked + " singular vectors being removed, but only " + samples.size() + " samples - try reducing the number of singular vectors with the -sv option, or increase the number of samples.");
			System.err.println("Cannot process samples.");
			System.exit(1);
		}
		if (svsBlanked * 2 > samples.size()) {
			System.err.println(svsBlanked + " singular vectors being removed, but only " + samples.size() + " samples - try reducing the number of singular vectors with the -sv option, or increase the number of samples.");
		}

		long[] totals = new long[samples.size()];
		Map<String, long[][]> arraysMap = new TreeMap<String, long[][]>();
		for (int i = 0; i < samples.size(); i++) {
			System.err.println("Reading " + i + " " + samples.get(i));
			ObjectInputStream in = new ObjectInputStream(new FileInputStream(samples.get(i)));
			@SuppressWarnings("unchecked") LinkedHashMap<String, int[]> vectors = (LinkedHashMap<String, int[]>) in.readObject();
			for (Map.Entry<String, int[]> entry : vectors.entrySet()) {
				String chr = entry.getKey();
				if ((limitChromosome == null) || limitChromosome.equals(chr)) {
					int[] toAdd = entry.getValue();
					long[][] array = arraysMap.get(chr);
					if (array == null) {
						array = new long[samples.size()][];
						arraysMap.put(chr, array);
					}
					if (array[i] == null) {
						array[i] = new long[0];
					}
					for (int o = 0; o < toAdd.length; o++) {
						int pos = (o * 200) / divider;
						if (toAdd[o] > 0) {
							if (array[i].length <= pos) {
								long[] newArray = new long[pos + pos / 2 + 1];
								for (int p = 0; p < array[i].length; p++) {
									newArray[p] = array[i][p];
								}
								array[i] = newArray;
							}
							array[i][pos] += toAdd[o];
							totals[i] += toAdd[o];
						}
					}
				}
			}
			in.close();
		}
		
		System.err.println("Normalising");
		int totalPositions = 0;
		for (String chr : arraysMap.keySet()) {
			long[][] array = arraysMap.get(chr);
			int longestSample = 0;
			for (int i = 0; i < samples.size(); i++) {
				if (array[i] != null) {
					longestSample = (longestSample > array[i].length ? longestSample : array[i].length);
				} else {
					array[i] = new long[0];
				}
			}
			totalPositions += longestSample;
		}
		long total = 0L;
		for (int i = 0; i < samples.size(); i++) {
			total += totals[i];
		}
		System.err.println("Average reads in each bucket: " + ((total * 1.0) / samples.size() / totalPositions));
		PrintStream out = new PrintStream(new BufferedOutputStream(System.out));
		// Work out which chunks of the genome we want to do SVD on.
		List<String> chunkChromosomes = new ArrayList<String>();
		List<Integer> chunkStarts = new ArrayList<Integer>();
		double[] scaleArray = new double[samples.size()];
		for (String chr : arraysMap.keySet()) {
			long[][] array = arraysMap.get(chr);
			boolean repeat = true;
			int start = 0;
			while (repeat) {
				repeat = false;
				long totalReads = 0;
				for (int i = 0; i < samples.size(); i++) {
					if (array[i].length > start) {
						totalReads += array[i][start];
						repeat = true;
					}
				}
				if ((totalReads * 50 * totalPositions > total) && (totalReads > samples.size() * minReads)) {
				//if (totalReads > 0) {
					chunkChromosomes.add(chr);
					chunkStarts.add(start);
					for (int i = 0; i < samples.size(); i++) {
						if (array[i].length > start) {
							scaleArray[i] += (1.0 * array[i][start]) * samples.size() / totalReads;
						}
					}
				}
				start++;
			}
		}
		System.err.println("Number of genome chunks: " + chunkChromosomes.size());
		if (intervals) {
			for (int o = 0; o < chunkChromosomes.size(); o++) {
				String chr = chunkChromosomes.get(o);
				int start = chunkStarts.get(o);
				System.out.println(chr + "\t" + (start * divider) + "\t" + ((start + 1) * divider));
			}
			System.exit(0);
		}
		double scale = 0.0;
		for (int i = 0; i < samples.size(); i++) {
			scale += scaleArray[i];
		}
		for (int i = 0; i < samples.size(); i++) {
			scaleArray[i] = scaleArray[i] * samples.size() / scale;
		}
		double[][] aArray = new double[chunkChromosomes.size()][samples.size()];
		double[][] eArray = new double[chunkChromosomes.size()][samples.size()];
		for (int o = 0; o < chunkChromosomes.size(); o++) {
			String chr = chunkChromosomes.get(o);
			int start = chunkStarts.get(o);
			long[][] array = arraysMap.get(chr);
			double sum = 0.0;
			for (int i = 0; i < samples.size(); i++) {
				sum += (array[i].length > start ? array[i][start] / scaleArray[i] : 0.0);
			}
			sum = sum / samples.size();
			for (int i = 0; i < samples.size(); i++) {
				aArray[o][i] = Math.log((array[i].length > start ? array[i][start] / scaleArray[i] / sum : 0.0) + 0.01);
				eArray[o][i] = scaleArray[i] * sum;
			}
		}
		Jama.Matrix A = new Jama.Matrix(aArray);
		if (samples.size() > chunkChromosomes.size()) {
			A = A.transpose();
		}
		long time1 = System.currentTimeMillis();
		Jama.SingularValueDecomposition decomp = new Jama.SingularValueDecomposition(A);
		long time2 = System.currentTimeMillis();
		System.err.println("Performed SVD in " + (time2 - time1) + "ms");
		Jama.Matrix S = decomp.getS();
		Process gnuplot = null;
		SplitPrintStream pipe = null;
		if (graph) {
			ProcessBuilder pb = new ProcessBuilder("gnuplot");
			//pb.redirectError(ProcessBuilder.Redirect.to(new File("/dev/null")));
			//pb.redirectOutput(ProcessBuilder.Redirect.to(new File("/dev/null")));
			pb.redirectError(ProcessBuilder.Redirect.INHERIT);
			pb.redirectOutput(ProcessBuilder.Redirect.INHERIT);
			gnuplot = pb.start();
			PrintStream pipe1 = new PrintStream(gnuplot.getOutputStream());
			String tempFileName = samples.get(0) + ".tempFile" + Math.random();
			PrintStream pipe2 = new PrintStream(new FileOutputStream(tempFileName + ".plot"));
			List<PrintStream> streams = new ArrayList<PrintStream>();
			streams.add(pipe1);
			streams.add(pipe2);
			pipe = new SplitPrintStream(streams);
			pipe.println("set terminal pdf noenhanced" + (fontScale == 1.0 ? "" : " fontscale " + fontScale) + " size " + graphSize);
			if (samples.size() <= chunkChromosomes.size()) {
				PrintStream svFile = new PrintStream(new FileOutputStream(tempFileName + ".SVs.csv"));
				for (int i = 0; i < samples.size(); i++) {
					svFile.println((i + 1) + "\t" + S.getArray()[i][i]);
				}
				svFile.flush();
				svFile.close();
				pipe.println("set output \"" + tempFileName + ".SVs.pdf\"");
				pipe.println("plot '" + tempFileName + ".SVs.csv' using 1:2:($1 - 0.4):($1 + 0.4):(0):2 with boxxy fillstyle solid notitle");
				pipe.flush();
			}
			pipe.println("yo(a) = (a > 13 && a <= 22 ? 27 - a : (a == 23 || a == 24 ? 14 : a))");
			pipe.println("xo(a) = (a == 22 || a == 21 ? 200000000 : (a == 20 || a == 19 ? 180000000 : (a == 24 ? 170000000 : (a == 16 || a == 17 || a == 18 ? 160000000 : (a == 15 ? 150000000 : (a == 14 ? 140000000 : 0))))))");
			pipe.println("set label \"Y\" at 166000000, 14 right");
			pipe.println("set label \"14\" at 136000000, 13 right");
			pipe.println("set label \"15\" at 146000000, 12 right");
			pipe.println("set label \"16\" at 156000000, 11 right");
			pipe.println("set label \"17\" at 156000000, 10 right");
			pipe.println("set label \"18\" at 156000000, 9 right");
			pipe.println("set label \"19\" at 176000000, 8 right");
			pipe.println("set label \"20\" at 176000000, 7 right");
			pipe.println("set label \"21\" at 196000000, 6 right");
			pipe.println("set label \"22\" at 196000000, 5 right");
			for (int i = 14; i <= 24; i++) {
				if (i == 23) {
					i++;
				}
				pipe.println("set arrow from xo(" + i + ") - 2000000, yo(" + i + ") to xo(" + i + "), yo(" + i + ") nohead");
				pipe.println("set arrow from xo(" + i + "), yo(" + i + ") - 0.5 to xo(" + i + "), yo(" + i + ") + 0.5 nohead front");
			}
			for (int i = 14; i <= 20; i++) {
				if ((i == 17) || (i == 19)) {
					i++;
				}
				pipe.println("set arrow from xo(" + i + "), yo(" + i + ") - 0.5 to xo(" + i + " + 1), yo(" + i + ") - 0.5 nohead");
			}
			pipe.println("set arrow from xo(14), yo(24) - 0.5 to xo(24), yo(24) - 0.5 nohead");
			pipe.println("set key above");
			pipe.println("set grid");
			pipe.println("set style fill solid border -1");
			pipe.println("set ytics 1");
			pipe.println("set xtics out nomirror");
			pipe.println("set ytics out nomirror");
			pipe.println("unset colorbox");
			pipe.println("min(a, b) = (a > b ? b : a)");
			pipe.println("set ytics (\"1\" 1, \"2\" 2, \"3\" 3, \"4\" 4, \"5\" 5, \"6\" 6, \"7\" 7, \"8\" 8, \"9\" 9, \"10\" 10, \"11\" 11, \"12\" 12, \"13\" 13, \"X\" 14)");
			pipe.println("set xtics (\"0\" 0, \"10M\" 10000000, \"20M\" 20000000, \"30M\" 30000000, \"40M\" 40000000, \"50M\" 50000000, \"60M\" 60000000, \"70M\" 70000000, \"80M\" 80000000, \"90M\" 90000000, \"100M\" 100000000, \"110M\" 110000000, \"120M\" 120000000, \"130M\" 130000000, \"140M\" 140000000, \"150M\" 150000000, \"160M\" 160000000, \"170M\" 170000000, \"180M\" 180000000, \"190M\" 190000000, \"200M\" 200000000, \"210M\" 210000000, \"220M\" 220000000, \"230M\" 230000000, \"240M\" 240000000, \"250M\" 250000000)");
			if (samples.size() <= chunkChromosomes.size()) {
				for (int i = 0; i < Math.min(svsBlanked * 2, samples.size()); i++) {
					PrintStream svFile = new PrintStream(new FileOutputStream(tempFileName + ".SV_" + i + "_chunk.csv"));
					boolean needBlankLine = false;
					for (String chr : arraysMap.keySet()) {
						if (needBlankLine) {
							svFile.println("");
							needBlankLine = false;
						}
						int lastStart = -1;
						for (int o = 0; o < chunkChromosomes.size(); o++) {
							String chunkChr = chunkChromosomes.get(o);
							if (chr.equals(chunkChr)) {
								int start = chunkStarts.get(o);
								if (start > lastStart + 1) {
									if (needBlankLine) {
										svFile.println("");
										needBlankLine = false;
									}
								}
								double val = decomp.getU().getArray()[o][i];
								svFile.println(chr + "\t" + (start * divider) + "\t" + val);
								svFile.println(chr + "\t" + (start * divider + divider) + "\t" + val);
								needBlankLine = true;
								lastStart = start;
							}
						}
					}
					svFile.flush();
					svFile.close();
					svFile = new PrintStream(new FileOutputStream(tempFileName + ".SV_" + i + "_sample.csv"));
					for (int o = 0; o < samples.size(); o++) {
						svFile.println(o + "\t" + decomp.getV().getArray()[o][i]);
					}
					svFile.flush();
					svFile.close();
					pipe.println("set output \"" + tempFileName + ".SV_" + i + ".pdf");
					pipe.println("plot [* to *] [-2.5 to 14.55] " + (cytoBands == null ? "" : "'< sed -e \"s/chr//;s/^X/23/;s/^Y/24/;s/gneg/0/;/acen/d;s/gvar/0/;s/stalk/0/;s/gpos//\" <" + cytoBands + " | grep -v \"_\"' using ($2 + xo($1)):(yo($1)):($2 + xo($1)):($3 + xo($1)):(yo($1) - 0.4):(yo($1) + 0.4):((65536 + 256 + 1) * int(220 - $5/4)) with boxxy fillstyle solid noborder linecolor rgb variable notitle, ") + "'< sed -e \"s/chr//;s/^X/23/;s/^Y/24/\" <" + tempFileName + ".SV_" + i + "_chunk.csv' using ($2 + xo($1)):(yo($1) + $3 * 4.0) with lines linewidth 2 linecolor rgb \"#000000\" title \"Singular vector " + (i + 1) + "\", '" + tempFileName + ".SV_" + i + "_sample.csv' using 1:($2 * 1.5 - 1.0):($1 - 0.4):($1 + 0.4):(-1):($2 * 1.5 - 1.0) with boxxy axes x2y1 fillstyle solid border -1 notitle");
				}
			}
			pipe.flush();
		}
		for (int i = 0; i < svsBlanked; i++) {
			S.getArray()[i][i] = 0.0;
		}
		Jama.Matrix Aprime = decomp.getU().times(S).times(decomp.getV().transpose());
		//Jama.Matrix Aprime = A;
		if (samples.size() > chunkChromosomes.size()) {
			Aprime = Aprime.transpose();
		}
		double[][] aPrimeArray = Aprime.getArray();
		if (dump) {
			for (int o = 0; o < chunkChromosomes.size(); o++) {
				String chunkChr = chunkChromosomes.get(o);
				int start = chunkStarts.get(o);
				System.out.print(chunkChr + "\t" + start);
				for (int i = 0; i < samples.size(); i++) {
					System.out.print("\t" + (Math.exp(aPrimeArray[o][i]) - 0.01));
				}
				System.out.println("");
			}
			System.exit(0);
		}
			
		double scoreSum[] = new double[samples.size()];
		double scoreSsum[] = new double[samples.size()];
		int scoreCount[] = new int[samples.size()];
		double posVariance[] = new double[aPrimeArray.length];
		double totalPosVariance = 0.0;
		int validChunks = 0;
		for (int o = 0; o < aPrimeArray.length; o++) {
			// Calculate variance of the normalised corrected read depth at this position. Exclude some of the outliers (say 3% on each side).
			double[] sortedDepths = new double[samples.size()];
			for (int i = 0; i < samples.size(); i++) {
				double val = Math.exp(aPrimeArray[o][i]) - 0.01;
				sortedDepths[i] = val;
			}
			Arrays.sort(sortedDepths);
			int count = 0;
			double sum = 0.0;
			double ssum = 0.0;
			for (int i = samples.size() / 30; i < samples.size() - (samples.size() / 30); i++) {
				double val = sortedDepths[i];
				sum += val;
				ssum += val * val;
				count++;
			}
			posVariance[o] = (ssum - (sum * sum / count)) / count;
			if (posVariance[o] < cutoff * cutoff) {
				for (int i = 0; i < samples.size(); i++) {
					double val = Math.exp(aPrimeArray[o][i]) - 1.01;
					scoreSum[i] += val;
					scoreSsum[i] += val * val;
					scoreCount[i]++;
				}
				totalPosVariance += posVariance[o];
				validChunks++;
			}
		}
		System.err.println("Number of low-noise genome chunks: " + validChunks);
		totalPosVariance = totalPosVariance / validChunks;
		System.err.println("Average noise: " + Math.sqrt(totalPosVariance));
		for (int sampleNo = 0; sampleNo < samples.size(); sampleNo++) {
			if (caseSamples.contains(samples.get(sampleNo))) {
				String tempFileName = samples.get(sampleNo) + ".tempFile" + Math.random();
				PrintWriter depthFile = null;
				PrintWriter cnvFile = null;
				if (graph) {
					depthFile = new PrintWriter(new BufferedWriter(new FileWriter(tempFileName + ".readDepth")), false);
					cnvFile = new PrintWriter(new BufferedWriter(new FileWriter(tempFileName + ".cnvs")), false);
				}
				PrintWriter dataFile = null;
				if (writeData) {
					dataFile = new PrintWriter(new BufferedWriter(new FileWriter(samples.get(sampleNo) + "." + divider + ".data")), false);
				}
				double sampleVariance = (scoreSsum[sampleNo] - (scoreSum[sampleNo] * scoreSum[sampleNo] / scoreCount[sampleNo])) / scoreCount[sampleNo];
				// sampleVariance now contains the variance for this sample, including all the CNVs.
				// We need to do an initial CNV detection, remove these from consideration, and recalculate.
				// Also, posVariance[i] contains the variance for position number i.
				// Run viterbi to get initial CNV calls
				List<State> cnvs = viterbi(arraysMap, false, null, null, posVariance, aPrimeArray, divider, logTransProb, minProb, sampleNo, sampleVariance, chunkChromosomes, chunkStarts, cutoffV, mosaic, totalPosVariance, aArray, errorType, eArray, null);
				// Now re-calculate sampleVariance, without those parts that have a detected CNV.
				double sum = 0.0;
				double ssum = 0.0;
				int count = 0;
				for (int o = 0; o < aPrimeArray.length; o++) {
					double val = Math.exp(aPrimeArray[o][sampleNo]) - 1.01;
					if (posVariance[o] < cutoff * cutoff) {
						String chr = chunkChromosomes.get(o);
						int start = chunkStarts.get(o) * divider;
						boolean notCovered = true;
						for (State cnv : cnvs) {
							if (chr.equals(cnv.getChr()) && (start >= cnv.getStart()) && (start < cnv.getEnd())) {
								notCovered = false;
							}
						}
						if (notCovered) {
							sum += val;
							ssum += val * val;
							count++;
						}
					}
				}
				double newSampleVariance = (ssum - (sum * sum / count)) / count;
				if (errorModel) {
					for (int o = 0; o < aPrimeArray.length; o++) {
						double val = Math.exp(aPrimeArray[o][sampleNo]) - 1.01;
						double eVal = eArray[o][sampleNo];
						if (posVariance[o] < cutoff * cutoff) {
							String chr = chunkChromosomes.get(o);
							int start = chunkStarts.get(o) * divider;
							boolean notCovered = true;
							for (State cnv : cnvs) {
								if (chr.equals(cnv.getChr()) && (start >= cnv.getStart()) && (start < cnv.getEnd())) {
									notCovered = false;
								}
							}
							if (notCovered) {
								System.out.println(Math.sqrt(posVariance[o] * newSampleVariance / totalPosVariance) + "\t" + Math.sqrt(posVariance[o] + newSampleVariance) + "\t" + (1.0 / Math.sqrt(eVal)) + "\t" + val);
							}
						}
					}
				} else {
					cnvs = viterbi(arraysMap, graph, depthFile, cnvFile, posVariance, aPrimeArray, divider, logTransProb, minProb, sampleNo, newSampleVariance, chunkChromosomes, chunkStarts, cutoffV, mosaic, totalPosVariance, aArray, errorType, eArray, dataFile);
					int delCount = 0;
					int dupCount = 0;
					for (State evalState : cnvs) {
						int blockLength = (int) ((evalState.getEnd() - evalState.getStart()) / divider);
						if (evalState.getState() == 1) {
							out.println(evalState.getChr() + "\t" + evalState.getStart() + "\t" + evalState.getEnd() + "\tDeletion\t" + evalState.getCount() + "\t" + blockLength + "\t" + evalState.getProb() + "\t" + (evalState.getProb() / blockLength) + "\t" + evalState.getProportion() + "\t" + samples.get(sampleNo));
							delCount++;
						} else if (evalState.getState() == 3) {
							out.println(evalState.getChr() + "\t" + evalState.getStart() + "\t" + evalState.getEnd() + "\tDuplication\t" + evalState.getCount() + "\t" + blockLength + "\t" + evalState.getProb() + "\t" + (evalState.getProb() / blockLength) + "\t" + evalState.getProportion() + "\t" + samples.get(sampleNo));
							dupCount++;
						}
					}
					System.err.println(Math.sqrt(sampleVariance) + "\t" + Math.sqrt(newSampleVariance) + "\t" + delCount + "\t" + dupCount + "\t" + totals[sampleNo] + "\t" + samples.get(sampleNo));
					if (writeData) {
						dataFile.flush();
						dataFile.close();
						dataFile = null;
					}
					if (graph) {
						depthFile.flush();
						depthFile.close();
						depthFile = null;
						cnvFile.println("200\t0\t0\t0\t1");
						cnvFile.println("200\t0\t0\t0\t3");
						cnvFile.flush();
						cnvFile.close();
						cnvFile = null;
						if ((delCount > 0) || (dupCount > 0) || allGraphs) {
							pipe.println("set output \"" + samples.get(sampleNo) + "." + divider + ".cnvs.pdf");
							pipe.println("plot [* to *] [0.5 to 14.55] " + (cytoBands == null ? "" : "'< sed -e \"s/chr//;s/^X/23/;s/^Y/24/;s/gneg/0/;/acen/d;s/gvar/0/;s/stalk/0/;s/gpos//\" <" + cytoBands + " | grep -v \"_\"' using ($2 + xo($1)):(yo($1)):($2 + xo($1)):($3 + xo($1)):(yo($1) - 0.4):(yo($1) + 0.4):((65536 + 256 + 1) * int(220 - $5/4)) with boxxy fillstyle solid noborder linecolor rgb variable notitle, ") + "'< sed -e \"s/chr//;s/^X/23/;s/^Y/24/\" <" + tempFileName + ".cnvs | grep \"1$\"' using ($2 + xo($1)):(yo($1)):($2 + xo($1)):($3 + xo($1)):(yo($1) + 0.4 - $4 * 0.6):(yo($1) + 0.4) with boxxy fillstyle solid noborder linecolor rgb \"#FFAAAA\" title \"" + delCount + " Deletion" + (delCount == 1 ? "" : "s") + "\", '< sed -e \"s/chr//;s/^X/23/;s/^Y/24/\" <" + tempFileName + ".cnvs | grep \"3$\"' using ($2 + xo($1)):(yo($1)):($2 + xo($1)):($3 + xo($1)):(yo($1) - 0.4):(yo($1) - 0.4 + $4 * 0.6) with boxxy fillstyle solid noborder linecolor rgb \"#8888FF\" title \"" + dupCount + " Duplication" + (dupCount == 1 ? "" : "s") + "\", '< sed -e \"s/chr//;s/^X/23/;s/^Y/24/\" <" + tempFileName + ".readDepth' using ($2 + xo($1)):(yo($1) - min($4, 1.25) * 0.4):(yo($1) + min($4, 1.25) * 0.4) with filledcurves linecolor rgb \"#33AA33\" notitle, '' using ($2 + xo($1)):(yo($1) + (min($3, 2.25) - 1.0) * 0.4) with lines linewidth 2 linecolor rgb \"#000000\" title \"" + samples.get(sampleNo) + "\"");
							pipe.flush();
						}
					}
				}
			}
		}
		out.flush();
		if (graph) {
			pipe.close();
			gnuplot.waitFor();
		}
	}

	public static List<State> viterbi(Map<String, long[][]> arraysMap, boolean graph, PrintWriter depthFile, PrintWriter cnvFile, double[] posVariance, double[][] aPrimeArray, int divider, double logTransProb, double minProb, int sampleNo, double sampleVariance, List<String> chunkChromosomes, List<Integer> chunkStarts, double cutoffV, boolean mosaic, double totalPosVariance, double[][] aArray, int errorType, double[][] eArray, PrintWriter dataFile) {
		boolean needBlankLine = false;
		double[] dArray = new double[aPrimeArray.length];
		double[] beforeArray = new double[aArray.length];
		double[] dEArray = new double[eArray.length];
		for (int i = 0; i < aPrimeArray.length; i++) {
			dArray[i] = Math.exp(aPrimeArray[i][sampleNo]) - 0.01;
			beforeArray[i] = Math.exp(aArray[i][sampleNo]) - 0.01;
			dEArray[i] = eArray[i][sampleNo];
		}
		List<State> retval = new ArrayList<State>();
		for (String chr : arraysMap.keySet()) {
			long time1 = System.currentTimeMillis();
			// dArray[o] contains the normalised corrected read depth for position chunkStarts.get(o) (divided by divider) for this chromosome.
			// Now do the Viterbi algorithm on each sample for this chromosome.
			// We need to hold onto three separate paths, corresponding to deletion, normal, and duplication.
			// The three paths have a common start of null, leading up to a last state.
			// There is a transition penalty for switching between states.
			// State probability is calculated using the stddev calculated by the sqrt of the sum of the sample and position variances.
			// cdf(x) is the cumulative density function for the Gaussian distribution from -inf to x with a stddev calculated above.
			// The deletion state probability is (1.0 - cdf(x - 0.5)) - divide the stddev by sqrt(2) because stddev of a poisson is proportional to sqrt(E)
			// The insertion state probability is cdf(x - 1.5)
			// The normal state probability is min(cdf(x - 1.0), 1.0 - cdf(x - 1.0)) - multiply the sttdev by sqrt(2)
			State delState = null;
			State norState = null;
			State dupState = null;
			double delProb = logTransProb;
			double norProb = 0.0;
			double dupProb = logTransProb;
			int lastStart = -1000;
			for (int o = 0; o < chunkChromosomes.size(); o++) {
				if (chr.equals(chunkChromosomes.get(o))) {
					int start = chunkStarts.get(o);
					if ((start > lastStart + 1) && needBlankLine) {
						if (graph) {
							depthFile.println("");
						}
						if (dataFile != null) {
							dataFile.println("");
						}
						needBlankLine = false;
					}
					double stddev = Math.sqrt(posVariance[o] * sampleVariance / totalPosVariance);
					if (errorType == 1) {
						// Additive error model
						stddev = Math.sqrt(posVariance[o] + sampleVariance);
					} else if (errorType == 2) {
						// Poisson error model
						stddev = 1.0 / Math.sqrt(dEArray[o]);
					}
					double val = dArray[o];
					double newDelProb = logProbDel(val, stddev, minProb, mosaic);
					double newDupProb = logProbDup(val, stddev, minProb, mosaic);
					if (graph) {
						depthFile.println(chr + "\t" + (start * divider) + "\t" + val + "\t" + stddev + "\t" + beforeArray[o]);
						depthFile.println(chr + "\t" + (start * divider + divider) + "\t" + val + "\t" + stddev + "\t" + beforeArray[o]);
					}
					if (dataFile != null) {
						dataFile.println(chr + "\t" + (start * divider) + "\t" + (start * divider + divider) + "\t" + val + "\t" + stddev + "\t" + beforeArray[o] + "\t" + newDelProb + "\t" + newDupProb);
					}
					needBlankLine = true;
					if (posVariance[o] < cutoffV * cutoffV) {
						double newNorProb = 0.0;
						double fromDel = delProb + newDelProb;
						double fromNor = norProb + logTransProb + newDelProb;
						double fromDup = dupProb + logTransProb + newDelProb;
						double nextDelProb;
						State nextDelState;
						if ((fromDel > fromNor) && (fromDel > fromDup)) {
							// The new delState should be taken from the old delState.
							nextDelState = new State(chr, start * divider, start * divider + divider, 1, delState, newDelProb - newNorProb, val);
							nextDelProb = fromDel;
						} else if (fromNor > fromDup) {
							nextDelState = new State(chr, start * divider, start * divider + divider, 1, norState, newDelProb - newNorProb, val);
							nextDelProb = fromNor;
						} else {
							nextDelState = new State(chr, start * divider, start * divider + divider, 1, dupState, newDelProb - newNorProb, val);
							nextDelProb = fromDup;
						}
						fromDel = delProb + logTransProb + newNorProb;
						fromNor = norProb + newNorProb;
						fromDup = dupProb + logTransProb + newNorProb;
						double nextNorProb;
						State nextNorState;
						if ((fromDel > fromNor) && (fromDel > fromDup)) {
							// The new norState should be taken from the old delState.
							nextNorState = new State(chr, start * divider, start * divider + divider, 2, delState, 0.0, val);
							nextNorProb = fromDel;
						} else if (fromNor > fromDup) {
							nextNorState = new State(chr, start * divider, start * divider + divider, 2, norState, 0.0, val);
							nextNorProb = fromNor;
						} else {
							nextNorState = new State(chr, start * divider, start * divider + divider, 2, dupState, 0.0, val);
							nextNorProb = fromDup;
						}
						fromDel = delProb + logTransProb + newDupProb;
						fromNor = norProb + logTransProb + newDupProb;
						fromDup = dupProb + newDupProb;
						double nextDupProb;
						State nextDupState;
						if ((fromDel > fromNor) && (fromDel > fromDup)) {
							// The new dupState should be taken from the old delState.
							nextDupState = new State(chr, start * divider, start * divider + divider, 3, delState, newDupProb - newNorProb, val);
							nextDupProb = fromDel;
						} else if (fromNor > fromDup) {
							nextDupState = new State(chr, start * divider, start * divider + divider, 3, norState, newDupProb - newNorProb, val);
							nextDupProb = fromNor;
						} else {
							nextDupState = new State(chr, start * divider, start * divider + divider, 3, dupState, newDupProb - newNorProb, val);
							nextDupProb = fromDup;
						}
						delState = nextDelState;
						delProb = nextDelProb;
						norState = nextNorState;
						norProb = nextNorProb;
						dupState = nextDupState;
						dupProb = nextDupProb;
					}
					lastStart = start;
				}
			}
			delProb = delProb + logTransProb;
			dupProb = dupProb + logTransProb;
			//System.err.println("Probabilities\t" + samples.get(sampleNo) + "\t" + chr + "\t" + delProb + "\t" + norProb + "\t" + dupProb);
			State evalState = (delProb > norProb) && (delProb > dupProb) ? delState : (norProb > dupProb ? norState : dupState);
			ArrayList<State> states = new ArrayList<State>();
			while (evalState != null) {
				states.add(evalState);
				evalState = evalState.getPrevious();
			}
			for (int p = states.size() - 1; p >= 0; p--) {
				evalState = states.get(p);
				int blockLength = (int) ((evalState.getEnd() - evalState.getStart()) / divider);
				if ((evalState.getState() == 1) || (evalState.getState() == 3)) {
					if (graph) {
						cnvFile.println(chr + "\t" + evalState.getStart() + "\t" + evalState.getEnd() + "\t" + ((evalState.getCount() * 1.0) / blockLength) + "\t" + evalState.getState());
					}
					retval.add(evalState);
				}
			}
		}
		return retval;
	}

	public static double logProbDel(double x, double stddev, double minProb, boolean mosaic) {
		if (mosaic) {
			return log(minProb + 1.0 / (minProb + exp(cdf(x - 1.0, stddev) + 8.0)));
		}
		if (stddev > 1.0 - Math.sqrt(0.5)) {
			x += stddev + Math.sqrt(0.5) - 1.0;
		}
		return log(minProb + 1.0 / (minProb + exp(cdf(x - 1.0, stddev) - cdf(0.5 - x, stddev / Math.sqrt(2.0)))));
	}

	public static double logProbDup(double x, double stddev, double minProb, boolean mosaic) {
		if (mosaic) {
			return log(minProb + 1.0 / (minProb + exp(cdf(1.0 - x, stddev) + 8.0)));
		}
		if (stddev > Math.sqrt(0.5) - 0.5) {
			x -= stddev - Math.sqrt(0.5) + 0.5;
		}
		return log(minProb + 1.0 / (minProb + exp(cdf(1.0 - x, stddev) - cdf(x - 1.5, stddev * Math.sqrt(2.0)))));
	}

	public static double cdf(double x, double stddev) {
		x = x / stddev;
		if (x < -6.0) {
			return 10.0 * (-2.0 * x * x / Math.PI - Math.sqrt(2.0)) / Math.log(10.0);
		} else {
			return log(0.5 * (1 + Math.signum(x) * Math.sqrt(1 - Math.exp(-2.0 * x * x / Math.PI))));
		}
	}

	public static double exp(double x) {
		return Math.exp(x * Math.log(10.0) / 10.0);
	}

	public static double log(double x) {
		return 10.0 * Math.log(x) / Math.log(10.0);
	}

	public static class State
	{
		private State previous;
		private String chr;
		private int start, end, state, count;
		private double prob, proportion;

		public State(String chr, int start, int end, int state, State previous, double prob, double proportion) {
			this.chr = chr;
			this.state = state;
			this.previous = previous;
			this.end = end;
			this.start = start;
			this.count = 1;
			this.prob = prob;
			this.proportion = proportion;
			if (previous != null) {
				if (state == previous.state) {
					this.start = previous.start;
					this.prob = prob + previous.prob;
					this.proportion = proportion + previous.proportion;
					this.count = previous.count + 1;
					this.previous = previous.previous;
					//System.err.println("New state(" + state + "), prob " + this.prob + "\t" + prob);
				}
			}
		}

		public State getPrevious() {
			return previous;
		}

		public String getChr() {
			return chr;
		}

		public int getState() {
			return state;
		}

		public int getStart() {
			return start;
		}

		public int getEnd() {
			return end;
		}

		public double getProb() {
			return prob;
		}

		public double getProportion() {
			return proportion / count;
		}

		public int getCount() {
			return count;
		}

		public String toString() {
			return start + "-" + end + " " + state + " -> " + previous;
		}
	}

	public static class SplitPrintStream
	{
		private List<PrintStream> streams;

		public SplitPrintStream(List<PrintStream> streams) {
			this.streams = streams;
		}

		public void println(String line) {
			for (PrintStream stream : streams) {
				stream.println(line);
			}
		}

		public void flush() {
			for (PrintStream stream : streams) {
				stream.flush();
			}
		}

		public void close() {
			for (PrintStream stream : streams) {
				stream.close();
			}
		}
	}
}
