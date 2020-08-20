import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;
import java.util.regex.Pattern;

/**
 * Reads a set of SavvyCNV data files, and performs joint calling of CNVs, assuming that all samples are in the same family.
 *
 * @author Matthew Wakeling
 */
public class SavvyCNVJointCaller
{
	public static final Pattern TAB = Pattern.compile("\t");

	public static void main(String[] args) throws Exception
	{
		List<String> samples = new ArrayList<String>();
		double transitionProb = 0.00001;
		double minProb = 0.00000000001;
		boolean mosaic = false;
		for (int i = 0; i < args.length; i++) {
			if ("-trans".equals(args[i])) {
				i++;
				transitionProb = Double.parseDouble(args[i]);
			} else if ("-minProb".equals(args[i])) {
				i++;
				minProb = SavvyCNV.exp(-Double.parseDouble(args[i]));
			} else if ("-mosaic".equals(args[i])) {
				mosaic = true;
				System.err.println("Using mosaic mode");
			} else {
				samples.add(args[i]);
			}
		}
		double logTransProb = SavvyCNV.log(transitionProb);
		//System.err.println("Transition probability " + logTransProb);
		System.err.println("Processing " + samples.size() + " samples");
		System.err.println("Using transition probability of " + transitionProb + " (phred " + logTransProb + ")");

		TreeSet<Interval<Map<String, ReadDepth>>> data = new TreeSet<Interval<Map<String, ReadDepth>>>();
		for (String sample : samples) {
			loadData(data, sample);
		}
		viterbi(data, samples, logTransProb, minProb, mosaic);
	}

	public static void loadData(TreeSet<Interval<Map<String, ReadDepth>>> data, String sample) throws IOException {
		BufferedReader in = new BufferedReader(new FileReader(sample));
		String line = in.readLine();
		while (line != null) {
			String[] split = TAB.split(line);
			if (split.length >= 5) {
				String chr = split[0];
				int start = Integer.parseInt(split[1]);
				int end = Integer.parseInt(split[2]);
				double val = Double.parseDouble(split[3]);
				double stddev = Double.parseDouble(split[4]);
				Interval<Map<String, ReadDepth>> search = new Interval<Map<String, ReadDepth>>(chr, start, end, null);
				Interval<Map<String, ReadDepth>> result = data.floor(search);
				if ((result != null) && chr.equals(result.getChromosome()) && (start == result.getStart())) {
					if (end != result.getEnd()) {
						throw new RuntimeException("Divider of multiple samples should be identical.");
					}
					result.getData().put(sample, new ReadDepth(val, stddev));
				} else {
					Map<String, ReadDepth> map = new HashMap<String, ReadDepth>();
					map.put(sample, new ReadDepth(val, stddev));
					data.add(new Interval<Map<String, ReadDepth>>(chr, start, end, map));
				}
			}
			line = in.readLine();
		}
	}

	public static void viterbi(TreeSet<Interval<Map<String, ReadDepth>>> data, List<String> samples, double logTransProb, double minProb, boolean mosaic) {
		// We are doing a viterbi algorithm, deciding between 3^samples states, corresponding to all possible combinations of
		// each sample being loss, normal, or gain. Each transition has the same penalty, regardless of how many samples it changes
		// the status of.
		// The states are numbered, according to sum of 3^sampleno * (0 for normal, 1 for loss, 2 for gain).
		int sampleCount = samples.size();
		int stateCount = 1;
		for (int i = 0; i < sampleCount; i++) {
			stateCount = stateCount * 3;
		}
		State[] states = new State[stateCount];
		double[] probabilities = new double[stateCount];
		String lastChromosome = null;
		for (Interval<Map<String, ReadDepth>> datum : data) {
			String chr = datum.getChromosome();
			if (!chr.equals(lastChromosome)) {
				if (lastChromosome != null) {
					printCalls(lastChromosome, states, probabilities, logTransProb, samples);
					//System.out.println(states[0]);
				}
				probabilities[0] = 0.0;
				states[0] = null;
				for (int i = 1; i < stateCount; i++) {
					probabilities[i] = logTransProb;
					states[i] = null;
				}
			}
			int start = datum.getStart();
			int end = datum.getEnd();
			Map<String, ReadDepth> depthMap = datum.getData();
			ReadDepth[] depths = new ReadDepth[sampleCount];
			double meanStddev = 0.0;
			int validCount = 0;
			//System.out.print(chr + "\t" + start + "\t" + end);
			for (int i = 0; i < sampleCount; i++) {
				depths[i] = depthMap.get(samples.get(i));
				if (depths[i] != null) {
					meanStddev += depths[i].getStddev();
					validCount++;
			//		System.out.print("\t" + depths[i].getVal() + "\t" + depths[i].getStddev());
			//	} else {
			//		System.out.print("\t-\t-");
				}
			}
			meanStddev = meanStddev / validCount;
			//System.out.print("\t" + meanStddev);
			double[] delProb = new double[sampleCount];
			double[] dupProb = new double[sampleCount];
			double[] newVal = new double[sampleCount];
			for (int i = 0; i < sampleCount; i++) {
				if (depths[i] != null) {
					delProb[i] = SavvyCNV.logProbDel(depths[i].getVal(), depths[i].getStddev(), minProb, mosaic);
					dupProb[i] = SavvyCNV.logProbDup(depths[i].getVal(), depths[i].getStddev(), minProb, mosaic);
					newVal[i] = depths[i].getVal();
		//			System.out.print("\t" + delProb[i] + "\t" + dupProb[i]);
				} else {
					delProb[i] = 0.0;
					dupProb[i] = 0.0;
					newVal[i] = 1.0;
		//			System.out.print("\t-\t-");
				}
			}
			double[] probabilityChanges = new double[stateCount];
			for (int i = 0; i < stateCount; i++) {
				int[] decodedState = decodeState(i, sampleCount);
				for (int sample = 0; sample < sampleCount; sample++) {
					probabilityChanges[i] += decodedState[sample] == 0 ? 0.0 : (decodedState[sample] == 1 ? delProb[sample] : dupProb[sample]);
				}
		//		System.out.print("\t" + probabilityChanges[i]);
			}
			State[] nextStates = new State[stateCount];
			double[] nextProbabilities = new double[stateCount];
			for (int i = 0; i < stateCount; i++) {
				double bestProb = -Double.MAX_VALUE;
				int bestFromState = -1;
				for (int fromState = 0; fromState < stateCount; fromState++) {
					double prob = probabilities[fromState] + probabilityChanges[i];
					if (i != fromState) {
						prob += logTransProb;
					}
					if (prob > bestProb) {
						bestProb = prob;
						bestFromState = fromState;
					}
				}
				nextStates[i] = new State(chr, start, end, i, states[bestFromState], delProb, dupProb, newVal);
				nextProbabilities[i] = bestProb;
			}
			states = nextStates;
			probabilities = nextProbabilities;
			//System.out.println("");
			lastChromosome = chr;
		}
		if (lastChromosome != null) {
			printCalls(lastChromosome, states, probabilities, logTransProb, samples);
			//System.out.println(states[0]);
		}
	}

	public static void printCalls(String chr, State[] states, double[] probabilities, double logTransProb, List<String> samples) {
		int sampleCount = samples.size();
		List<State> reverseStates = new ArrayList<State>();
		int bestState = -1;
		double bestProb = -Double.MAX_VALUE;
		for (int i = 0; i < states.length; i++) {
			double prob = probabilities[i] + (i == 0 ? 0.0 : logTransProb);
			if (prob > bestProb) {
				bestProb = prob;
				bestState = i;
			}
		}
		{
			State state = states[bestState];
			while (state != null) {
				reverseStates.add(state);
				state = state.getPrevious();
			}
		}

		for (int i = reverseStates.size() - 1; i >= 0; i--) {
			State state = reverseStates.get(i);
			if (state.getState() != 0) {
				int[] decodedState = decodeState(state.getState(), sampleCount);
				for (int sample = 0; sample < sampleCount; sample++) {
					System.out.println(chr + "\t" + state.getStart() + "\t" + state.getEnd() + "\t" + (decodedState[sample] == 0 ? "Normal" : (decodedState[sample] == 1 ? "Deletion" : "Duplication")) + "\t" + state.getCount()[sample] + "\t" + state.getDelProb()[sample] + "\t" + state.getDupProb()[sample] + "\t" + (decodedState[sample] == 0 ? "-" : (decodedState[sample] == 1 ? "" + state.getDelProb()[sample] / state.getCount()[sample] : state.getDupProb()[sample] / state.getCount()[sample])) + "\t" + state.getProportion()[sample] + "\t" + samples.get(sample));
				}
			}
		}
	}

	public static int[] decodeState(int state, int sampleCount) {
		int[] retval = new int[sampleCount];
		for (int i = 0; i < sampleCount; i++) {
			retval[i] = state % 3;
			state = state / 3;
		}
		return retval;
	}

	public static class State
	{
		private State previous;
		private String chr;
		private int start, end, state;
		private int[] count;
		private double[] delProb, dupProb, proportion;

		public State(String chr, int start, int end, int state, State previous, double[] delProb, double[] dupProb, double[] proportion) {
			this.chr = chr;
			this.state = state;
			this.end = end;
			if ((previous != null) && (state == previous.state)) {
				this.start = previous.start;
				double[] newDelProb = new double[delProb.length];
				double[] newDupProb = new double[delProb.length];
				double[] newProportion = new double[delProb.length];
				this.count = new int[delProb.length];
				for (int i = 0; i < delProb.length; i++) {
					if ((delProb[i] != 0.0) || (dupProb[i] != 0.0)) {
						newDelProb[i] = delProb[i] + previous.delProb[i];
						newDupProb[i] = dupProb[i] + previous.dupProb[i];
						newProportion[i] = proportion[i] + previous.proportion[i];
						this.count[i] = previous.count[i] + 1;
					} else {
						newDelProb[i] = previous.delProb[i];
						newDupProb[i] = previous.dupProb[i];
						newProportion[i] = previous.proportion[i];
						this.count[i] = previous.count[i];
					}
				}
				this.delProb = newDelProb;
				this.dupProb = newDupProb;
				this.proportion = newProportion;
				this.previous = previous.previous;
			} else {
				this.start = start;
				this.previous = previous;
				this.count = new int[delProb.length];
				this.proportion = new double[delProb.length];
				for (int i = 0; i < count.length; i++) {
					if ((delProb[i] != 0.0) || (dupProb[i] != 0.0)) {
						this.count[i] = 1;
						this.proportion[i] = proportion[i];
					} else {
						this.count[i] = 0;
						this.proportion[i] = 0.0;
					}
				}
				this.delProb = delProb;
				this.dupProb = dupProb;
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

		public double[] getDelProb() {
			return delProb;
		}

		public double[] getDupProb() {
			return dupProb;
		}

		public double[] getProportion() {
			double[] retval = new double[proportion.length];
			for (int i = 0; i < proportion.length; i++) {
				retval[i] = proportion[i] / count[i];
			}
			return retval;
		}

		public int[] getCount() {
			return count;
		}

		public String toString() {
			return start + "-" + end + " " + state + " -> " + previous;
		}
	}


	public static class ReadDepth
	{
		private double val, stddev;

		public ReadDepth(double val, double stddev) {
			this.val = val;
			this.stddev = stddev;
		}

		public double getVal() {
			return val;
		}

		public double getStddev() {
			return stddev;
		}

		public String toString() {
			return val + "+-" + stddev;
		}
	}
}
