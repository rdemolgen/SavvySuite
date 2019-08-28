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
import java.util.Collections;
import java.util.HashMap;
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

public class SavvyHomozygosity
{
	public static final int WIDTH = 100000;
	public static final int READ_LENGTH = 200;
	public static final Map<String, String> addChrMap = new HashMap<String, String>();
	public static final Map<String, String> delChrMap = new HashMap<String, String>();
	static {
		for (int i = 1; i < 23; i++) {
			addChrMap.put("" + i, "chr" + i);
			delChrMap.put("chr" + i, "" + i);
		}
		addChrMap.put("X", "chrX");
		delChrMap.put("chrX", "X");
		addChrMap.put("Y", "chrY");
		delChrMap.put("chrY", "Y");
		addChrMap.put("MT", "chrM");
		delChrMap.put("chrM", "MT");
	}


	/*
	 * Method - We want to find homozygous regions using single off-target reads from the BAM file, using linkage disequilibrium.
	 * The files are big, so we need to ensure performance is sensible.
	 * Procedure:
	 *
	 * 1. Load a variant from the VCF. Check to see if there is a read covering it from the BAM. If not, then discard it.
	 * 2. Otherwise, iterate through all variants stored to the left.
	 *   2.1. If the variant is more than WIDTH to the left of the current variant, then discard it.
	 *   2.2. Otherwise, calculate linkage from VCF variants, produce score from reads.
	 *
	 * So, we need to keep track of both variants from the VCF and "genotypes" from the BAM file, between the current variant and WIDTH bases to the left in the same chromosome.
	 */

	/**
	 * Process a sample to detect homozygous regions. The first argument is a VCF file containing common variants for lots of WGS samples. The second argument is a BAM file containing the sample to be analysed.
	 *
	 * @author Matthew Wakeling
	 */
	public static void main(String[] args) throws IOException, ClassNotFoundException {
		int variants = 0;
		double startTransition = 80.0;
		double negativeMult = 20.0;
		// Variants are read from an ObjectInputStream as VariantArray objects.
		ObjectInputStream vcf = null;
		try {
			vcf = new ObjectInputStream(new BufferedInputStream(new GZIPInputStream(new FileInputStream(args[0]))));
		} catch (ZipException e) {
			vcf = new ObjectInputStream(new BufferedInputStream(new FileInputStream(args[0])));
		}
		BamReader bamReader = new BamReader(args[1]);
		boolean outputPoints = false;
		boolean gotTranslation = false;
		Map<String, String> translation = null;
		for (int i = 2; i < args.length; i++) {
			if ("-points".equals(args[i])) {
				outputPoints = true;
			} else if ("-trans".equals(args[i])) {
				i++;
				startTransition = Double.parseDouble(args[i]);
			} else if ("-neg".equals(args[i])) {
				i++;
				negativeMult = Double.parseDouble(args[i]);
			}
		}
		TreeMap<Integer, VariantArray> storedVariants = null;
		TreeMap<Integer, Boolean> storedBases = null;
		String currentChr = "";
		MultiViterbi viterbi = null;
		boolean hasMoreVariants = true;
		VariantArray vArray = null;
		try {
			vArray = (VariantArray) vcf.readObject();
			boolean varChr = vArray.getChr().startsWith("chr");
			boolean bamChr = false;
			try {
				bamReader.getCounts("1", 20000);
			} catch (Exception e) {
				bamChr = true;
			}
			if (varChr && (!bamChr)) {
				translation = delChrMap;
			} else if (bamChr && (!varChr)) {
				translation = addChrMap;
			}
		} catch (EOFException e) {
			hasMoreVariants = false;
		}
		while (hasMoreVariants) {
			String contextChr = vArray.getChr();
			if (translation != null) {
				contextChr = translation.get(contextChr);
			}
			if (contextChr == null) {
				if (!gotTranslation) {
					System.err.println("Variants for some chromosomes are ignored - there is no translation between naming schemes. Chromosome that triggered this message is \"" + vArray.getChr() + "\"");
					gotTranslation = true;
				}
			} else {
				if (!(currentChr.equals(contextChr))) {
					if (viterbi != null) {
						viterbi.finish();
					}
					storedVariants = new TreeMap<Integer, VariantArray>();
					currentChr = contextChr;
					storedBases = new TreeMap<Integer, Boolean>();
					viterbi = new MultiViterbi(currentChr, startTransition);
				}
				int[] bases = bamReader.getCounts(contextChr, vArray.getPosition());
				int wildCount = bases[baseIndex(vArray.getWild())];
				int varCount = bases[baseIndex(vArray.getVar())];
				if (wildCount + varCount > 0) {
					//System.out.println(vArray.getChr() + ":" + vArray.getStart() + "\t" + bases[0] + "\t" + bases[1] + "\t" + bases[2] + "\t" + bases[3] + "\t" + vArray.getWild() + ">" + vArray.getVar() + "\t" + wildCount + "\t" + varCount);
					if ((wildCount == 0) || (varCount == 0)) {
						if (wildCount + varCount > 2) {
							viterbi.addSignal(new ViterbiSignal(vArray.getPosition(), 1.1));
							if (outputPoints) {
								System.err.println(currentChr + "\t" + vArray.getPosition() + "\t" + vArray.getPosition() + "\t1.1");
							}
						}
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
									viterbi.addSignal(new ViterbiSignal((vArray.getPosition() + storedEntry.getKey()) / 2, (rsquaredMinus + 0.8) * 5.0 * negativeMult));
									if (outputPoints) {
										System.err.println(currentChr + "\t" + vArray.getPosition() + "\t" + storedEntry.getKey() + "\t" + ((rsquaredMinus + 0.8) * 5.0));
									}
								} else if (rsquaredMinus > 0.8) {
									viterbi.addSignal(new ViterbiSignal((vArray.getPosition() + storedEntry.getKey()) / 2, (rsquaredMinus - 0.8) * 5.0));
									if (outputPoints) {
										System.err.println(currentChr + "\t" + vArray.getPosition() + "\t" + storedEntry.getKey() + "\t" + ((rsquaredMinus - 0.8) * 5.0));
									}
								}
							}
						}
						storedVariants.put(vArray.getPosition(), vArray);
						storedBases.put(vArray.getPosition(), thisBase);
					} else {
						viterbi.addSignal(new ViterbiSignal(vArray.getPosition(), -6.0));
						if (outputPoints) {
							System.err.println(currentChr + "\t" + vArray.getPosition() + "\t" + vArray.getPosition() + "\t-1.1");
						}
					}
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

	public static class MultiViterbi
	{
		// Run the Viterbi algorithm with multiple transition penalties, to merge together large homozygous regions.
		private List<Viterbi> viterbis;

		public MultiViterbi(String chromosome, double startTransition) {
			viterbis = new ArrayList<Viterbi>();
			for (double transition = startTransition; transition < startTransition * 200.0; transition = transition * 1.5) {
				viterbis.add(new Viterbi(chromosome, transition));
			}
		}

		public void addSignal(ViterbiSignal signal) {
			for (Viterbi viterbi : viterbis) {
				viterbi.addSignal(signal);
			}
		}

		public void finish() {
			ArrayList<Interval<Boolean>> intervals = new ArrayList<Interval<Boolean>>();
			for (Viterbi viterbi : viterbis) {
				//System.err.println("Finishing Viterbi(" + viterbi.getTransition() + ")");
				TreeSet<Interval<Set<String>>> union = new TreeSet<Interval<Set<String>>>();
				viterbi.finishWithIntervals(union);
				for (Interval<Set<String>> interval : union) {
					//System.err.println("Viterbi(" + viterbi.getTransition() + "): Copying interval " + interval);
					intervals.add(new Interval<Boolean>(interval.getChromosome(), interval.getStart(), 0, Boolean.FALSE));
					intervals.add(new Interval<Boolean>(interval.getChromosome(), interval.getEnd(), 0, Boolean.TRUE));
				}
			}
			Collections.sort(intervals);
			int popularity = 0;
			for (Interval<Boolean> interval : intervals) {
				if (Boolean.FALSE == interval.getData()) {
					if (popularity == 0) {
						System.out.print(interval.getChromosome() + "\t" + interval.getStart() + "\t");
					}
					popularity++;
				} else {
					if (popularity == 1) {
						System.out.println(interval.getStart());
					}
					popularity--;
				}
			}
		}
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
			//System.err.println("Starting Viterbi(" + chromosome + ", " + transition + ")");
			this.chromosome = chromosome;
			this.hetProb = 0.0;
			this.homProb = 0.0;
			this.transition = transition;
		}

		public double getTransition() {
			return transition;
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
