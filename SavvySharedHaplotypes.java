import htsjdk.samtools.*;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.*;
import java.io.BufferedInputStream;
import java.io.EOFException;
import java.io.IOException;
import java.io.File;
import java.io.FileInputStream;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipException;

public class SavvySharedHaplotypes
{
	public static final int WIDTH = 100000;
	public static final int READ_LENGTH = 200;

	/**
	 * Process a sample to detect homozygous regions or shared haplotypes. The first argument is a VCF file containing common variants for lots of WGS samples. The second argument is a BAM file containing the sample to be analysed. The optional third argument is the BAM file of the other sample.
	 *
	 * @author Matthew Wakeling
	 */
	public static void main(String[] args) throws IOException, ClassNotFoundException {
		int variants = 0;
		ObjectInputStream vcf = null;
		try {
			vcf = new ObjectInputStream(new BufferedInputStream(new GZIPInputStream(new FileInputStream(args[0]))));
		} catch (ZipException e) {
			vcf = new ObjectInputStream(new BufferedInputStream(new FileInputStream(args[0])));
		}
		BamReader bamReader = new BamReader(args[1]);
		BamReader bamReader2 = null;
		double negativeMultiplier = 20.0; // This is how much less frequent the discordant read pairs must be before calling a discordant region.
		if (args.length > 2) {
			bamReader2 = new BamReader(args[2]);
			negativeMultiplier = 10.0;
		}
		TreeMap<Integer, VariantArray> storedVariants = null;
		TreeMap<Integer, Boolean> storedBases = null;
		TreeMap<Integer, VariantArray> storedVariants2 = null;
		TreeMap<Integer, Boolean> storedBases2 = null;
		String currentChr = "";
		Viterbi viterbi = null;
		boolean hasMoreVariants = true;
		VariantArray vArray = null;
		try {
			vArray = (VariantArray) vcf.readObject();
		} catch (EOFException e) {
			hasMoreVariants = false;
		}
		while (hasMoreVariants) {
			String contextChr = vArray.getChr();
			if (!(currentChr.equals(contextChr))) {
				if (viterbi != null) {
					viterbi.finish();
				}
				storedVariants = new TreeMap<Integer, VariantArray>();
				currentChr = contextChr;
				storedBases = new TreeMap<Integer, Boolean>();
				viterbi = new Viterbi(currentChr, 40.0);
				if (bamReader2 != null) {
					storedVariants2 = new TreeMap<Integer, VariantArray>();
					storedBases2 = new TreeMap<Integer, Boolean>();
				}
			}
			int[] bases = bamReader.getCounts(contextChr, vArray.getPosition());
			int[] bases2 = bases;
			if (bamReader2 != null) {
				bases2 = bamReader2.getCounts(contextChr, vArray.getPosition());
			}
			int wildCount = bases[baseIndex(vArray.getWild())];
			int varCount = bases[baseIndex(vArray.getVar())];
			if (wildCount + varCount > 0) {
				//System.out.println(contextChr + ":" + vArray.getPosition() + "\t" + bases[0] + "\t" + bases[1] + "\t" + bases[2] + "\t" + bases[3] + "\t" + context.getAlleles() + "\t" + wildCount + "\t" + varCount);
				if ((wildCount == 0) || (varCount == 0)) {
					if ((bamReader2 == null) && (wildCount + varCount > 2)) {
						viterbi.addSignal(new ViterbiSignal(vArray.getPosition(), 1.1));
						System.err.println(currentChr + "\t" + vArray.getPosition() + "\t1.1");
					}
					createSignal(wildCount, varCount, storedBases, storedVariants, vArray, viterbi, negativeMultiplier, currentChr);
				} else if (bamReader2 == null) {
					viterbi.addSignal(new ViterbiSignal(vArray.getPosition(), -6.0));
					System.err.println(currentChr + "\t" + vArray.getPosition() + "\t-1.1");
				}
			}
			wildCount = bases2[baseIndex(vArray.getWild())];
			varCount = bases2[baseIndex(vArray.getVar())];
			if ((wildCount + varCount > 0) && ((wildCount == 0) || (varCount == 0))) {
				if (bamReader2 != null) {
					createSignal(wildCount, varCount, storedBases2, storedVariants2, vArray, viterbi, negativeMultiplier, currentChr);
				}
				boolean thisBase = varCount > 0;
				storedBases.put(vArray.getPosition(), thisBase);
				storedVariants.put(vArray.getPosition(), vArray);
			}
			if (bamReader2 != null) {
				wildCount = bases[baseIndex(vArray.getWild())];
				varCount = bases[baseIndex(vArray.getVar())];
				if ((wildCount + varCount > 0) && ((wildCount == 0) || (varCount == 0))) {
					boolean thisBase = varCount > 0;
					storedBases2.put(vArray.getPosition(), thisBase);
					storedVariants2.put(vArray.getPosition(), vArray);
				}
			}
			try {
				vArray = (VariantArray) vcf.readObject();
			} catch (EOFException e) {
				hasMoreVariants = false;
			}
		}
		if (viterbi != null) {
			viterbi.finish();
		}
	}

	public static void createSignal(int wildCount, int varCount, TreeMap<Integer, Boolean> storedBases, TreeMap<Integer, VariantArray> storedVariants, VariantArray vArray, Viterbi viterbi, double negativeMultiplier, String currentChr) {
		boolean thisBase = varCount > 0;
		Iterator<Map.Entry<Integer, Boolean>> storedBaseIter = storedBases.entrySet().iterator();
		while (storedBaseIter.hasNext()) {
			Map.Entry<Integer, Boolean> storedBase = storedBaseIter.next();
			if (storedBase.getKey() < vArray.getPosition() - WIDTH) {
				storedBaseIter.remove();
			}
		}
		Iterator<Map.Entry<Integer, VariantArray>> storedIter = storedVariants.entrySet().iterator();
		while (storedIter.hasNext()) {
			Map.Entry<Integer, VariantArray> storedEntry = storedIter.next();
			if (storedEntry.getKey() < vArray.getPosition() - WIDTH) {
				storedIter.remove();
			} else if (storedEntry.getKey() < vArray.getPosition() - READ_LENGTH) {
				boolean otherBase = storedBases.get(storedEntry.getKey());
				double rsquared = vArray.getRsquared(storedEntry.getValue());
				double rsquaredMinus = (thisBase == otherBase) ? rsquared : -rsquared;
				if (rsquaredMinus < -0.8) {
					viterbi.addSignal(new ViterbiSignal((vArray.getPosition() + storedEntry.getKey()) / 2, (rsquaredMinus + 0.8) * 5.0 * negativeMultiplier));
					System.err.println(currentChr + "\t" + ((vArray.getPosition() + storedEntry.getKey()) / 2) + "\t" + ((rsquaredMinus + 0.8) * 5.0) + "\t" + vArray.getPosition() + "\t" + storedEntry.getKey());
				} else if (rsquaredMinus > 0.8) {
					viterbi.addSignal(new ViterbiSignal((vArray.getPosition() + storedEntry.getKey()) / 2, (rsquaredMinus - 0.8) * 5.0));
					System.err.println(currentChr + "\t" + ((vArray.getPosition() + storedEntry.getKey()) / 2) + "\t" + ((rsquaredMinus - 0.8) * 5.0) + "\t" + vArray.getPosition() + "\t" + storedEntry.getKey());
				}
			}
		}
	}

	public static class BamReader
	{
		public static final Pattern SIMPLE_CIGAR = Pattern.compile("[0-9]*M");

		private SAMRecordIterator iter;
		private SamReader in;
		private String currentChromosome;
		private LinkedList<SAMRecord> records;
		private SAMRecord lastRecord = null;

		public BamReader(String fileName) throws IOException {
			SamReaderFactory factory = SamReaderFactory.makeDefault();
			in = factory.open(new File(fileName));
		}

		public int[] getCounts(String chr, int pos) throws IOException {
			if (!chr.equals(currentChromosome)) {
				if (iter != null) {
					iter.close();
				}
				currentChromosome = chr;
				records = new LinkedList<SAMRecord>();
				iter = in.queryOverlapping(chr, 0, 0);
				lastRecord = null;
			}
			int[] retval = new int[5];
			while (iter.hasNext() && ((lastRecord == null) || (lastRecord.getAlignmentStart() <= pos))) {
				lastRecord = iter.next();
				if ((lastRecord.getAlignmentEnd() >= pos) && SIMPLE_CIGAR.matcher(lastRecord.getCigarString()).matches()) {
					records.add(lastRecord);
				}
			}
			Iterator<SAMRecord> recordIter = records.iterator();
			while (recordIter.hasNext()) {
				SAMRecord record = recordIter.next();
				if (record.getAlignmentEnd() < pos) {
					recordIter.remove();
				} else {
					int relative = pos - record.getAlignmentStart();
					String bases = record.getReadString();
					if ((relative >= 0) && (relative < bases.length())) {
						char base = bases.charAt(relative);
						retval[baseIndex(base)]++;
					}
				}
			}
			return retval;
		}
	}

	public static int baseIndex(char base) {
		int retval = 4;
		switch (base) {
			case 'A':
				retval = 0;
				break;
			case 'C':
				retval = 1;
				break;
			case 'G':
				retval = 2;
				break;
			case 'T':
				retval = 3;
				break;
		}
		return retval;
	}

	public static class Viterbi
	{
		// Viterbi algorithm. We need to hold onto two independent paths, corresponding to two end states, which are normal and homozygous.
		// There is a transition penalty for switching betweeen normal and homozygous, but not from chromosome start or to chromosome end.
		private String chromosome;
		private TreeSet<ViterbiSignal> signals = new TreeSet<ViterbiSignal>();
		private ViterbiState hetState, homState;
		private double hetProb, homProb;
		private double transition;


		public Viterbi(String chromosome, double transition) {
			this.chromosome = chromosome;
			this.hetProb = 0.0;
			this.homProb = 0.0;
			this.transition = transition;
		}

		public void addSignal(ViterbiSignal signal) {
			signals.add(signal);
			while (signals.first().getPos() < signal.getPos() - WIDTH) {
				ViterbiSignal first = signals.first();
				signals.remove(first);
				processSignal(first);
			}
		}

		public void finish() {
			for (ViterbiSignal signal : signals) {
				processSignal(signal);
			}
			ArrayList<ViterbiState> states = new ArrayList<ViterbiState>();
			ViterbiState evalState = hetProb > homProb ? hetState : homState;
			while (evalState != null) {
				states.add(evalState);
				evalState = evalState.getPrevious();
			}
			for (int i = states.size() - 1; i >= 0; i--) {
				evalState = states.get(i);
				if (evalState.getState() == 2) {
					System.out.println(chromosome + "\t" + evalState.getStart() + "\t" + evalState.getEnd());
				}
			}
		}

		public void finishWithIntervals(TreeSet<Interval<Set<String>>> intervals) {
			for (ViterbiSignal signal : signals) {
				processSignal(signal);
			}
			ArrayList<ViterbiState> states = new ArrayList<ViterbiState>();
			ViterbiState evalState = hetProb > homProb ? hetState : homState;
			Set<String> data = new HashSet<String>();
			data.add("n");
			while (evalState != null) {
				states.add(evalState);
				if (evalState.getState() == 2) {
					intervals.add(new Interval<Set<String>>(chromosome, evalState.getStart(), evalState.getEnd(), data));
				}
				evalState = evalState.getPrevious();
			}
		}

		public void processSignal(ViterbiSignal signal) {
			double hetToHet = hetProb - signal.getSignal();
			double homToHet = homProb - transition - signal.getSignal();
			double newHetProb = Math.max(hetToHet, homToHet);
			double hetToHom = hetProb + signal.getSignal() - transition;
			double homToHom = homProb + signal.getSignal();
			double newHomProb = Math.max(hetToHom, homToHom);
			ViterbiState newHetState = new ViterbiState(signal.getPos(), 1, hetToHet > homToHet ? hetState : homState);
			homState = new ViterbiState(signal.getPos(), 2, hetToHom > homToHom ? hetState : homState);
			hetState = newHetState;
			hetProb = newHetProb;
			homProb = newHomProb;
		}
	}

	public static class ViterbiSignal implements Comparable<ViterbiSignal>
	{
		private static int seqStatic = 0;
		private int pos, seq;
		private double signal;

		public ViterbiSignal(int pos, double signal) {
			this.pos = pos;
			this.seq = seqStatic++;
			this.signal = signal;
		}

		public int getPos() {
			return pos;
		}

		public double getSignal() {
			return signal;
		}

		public int compareTo(ViterbiSignal v) {
			if (pos > v.pos) {
				return 1;
			} else if (pos < v.pos) {
				return -1;
			} else if (seq > v.seq) {
				return 1;
			} else if (seq < v.seq) {
				return -1;
			}
			return 0;
		}
	}

	public static class ViterbiState
	{
		private ViterbiState previous;
		private int start, end, state;

		public ViterbiState(int pos, int state, ViterbiState previous) {
			this.state = state;
			this.previous = previous;
			this.end = pos;
			this.start = pos;
			if (previous != null) {
				if (state == previous.state) {
					start = previous.start;
					this.previous = previous.previous;
				}
			}
		}

		public ViterbiState getPrevious() {
			return previous;
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

		public String toString() {
			return start + "-" + end + " " + state + " -> " + previous;
		}
	}
}
