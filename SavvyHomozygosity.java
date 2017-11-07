import htsjdk.samtools.*;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.*;
import java.io.IOException;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Pattern;

public class SavvyHomozygosity
{
	public static final int WIDTH = 100000;
	public static final Pattern INDEL = Pattern.compile(".*[ACGT][ACGT].*");

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
	public static void main(String[] args) throws IOException {
		int variants = 0;
		VCFFileReader vcf = new VCFFileReader(new File(args[0]));
		BamReader bamReader = new BamReader(args[1]);
		List<String> vcfSampleNames = vcf.getFileHeader().getGenotypeSamples();
		TreeMap<Integer, VariantContext> storedVariants = null;
		TreeMap<Integer, Boolean> storedBases = null;
		String currentChr = "";
		MultiViterbi viterbi = null;
		for (VariantContext context : vcf) {
			@SuppressWarnings("deprecation") String contextChr = context.getChr();
			if (!(currentChr.equals(contextChr))) {
				if (viterbi != null) {
					viterbi.finish();
				}
				storedVariants = new TreeMap<Integer, VariantContext>();
				currentChr = contextChr;
				storedBases = new TreeMap<Integer, Boolean>();
				viterbi = new MultiViterbi(currentChr);
			}
			String alleles = "" + context.getAlleles();
			boolean indelMatch = INDEL.matcher(alleles).matches();
			if (!indelMatch) {
				int ac = 0;
				int an = 0;
				for (String sample : vcfSampleNames) {
					Genotype contextGt = context.getGenotype(sample);
					if (contextGt.isHomVar()) {
						ac += 2;
					} else if (contextGt.isHet()) {
						ac++;
					}
					an += 2;
				}
				if (ac * 4 > an) {
					int[] bases = bamReader.getCounts(contextChr, context.getStart());
					int wildCount = bases[baseIndex(context.getAlleles().get(0).toString().charAt(0))];
					int varCount = bases[baseIndex(context.getAlleles().get(1).toString().charAt(0))];
					if (wildCount + varCount > 0) {
						//System.out.println(context.getChr() + ":" + context.getStart() + "\t" + bases[0] + "\t" + bases[1] + "\t" + bases[2] + "\t" + bases[3] + "\t" + context.getAlleles() + "\t" + wildCount + "\t" + varCount);
						if ((wildCount == 0) || (varCount == 0)) {
							if (wildCount + varCount > 2) {
								viterbi.addSignal(new ViterbiSignal(context.getStart(), 1.1));
								//System.err.println(currentChr + "\t" + context.getStart() + "\t1.1");
							}
							boolean thisBase = varCount > 0;
							Iterator<Map.Entry<Integer, Boolean>> storedBaseIter = storedBases.entrySet().iterator();
							while (storedBaseIter.hasNext()) {
								Map.Entry<Integer, Boolean> storedBase = storedBaseIter.next();
								if (storedBase.getKey() < context.getStart() - WIDTH) {
									storedBaseIter.remove();
								}
							}
							Iterator<Map.Entry<Integer, VariantContext>> storedIter = storedVariants.entrySet().iterator();
							while (storedIter.hasNext()) {
								Map.Entry<Integer, VariantContext> storedEntry = storedIter.next();
								if (storedEntry.getKey() < context.getStart() - WIDTH) {
									storedIter.remove();
								} else {
									boolean otherBase = storedBases.get(storedEntry.getKey());
									int[] grid = new int[9];
									for (String sample : vcfSampleNames) {
										Genotype contextGt = context.getGenotype(sample);
										Genotype storedGt = storedEntry.getValue().getGenotype(sample);
										if (contextGt.isHomVar()) {
											if (storedGt.isHomVar()) {
												grid[8]++;
											} else if (storedGt.isHet()) {
												grid[7]++;
											} else {
												grid[6]++;
											}
										} else if (contextGt.isHet()) {
											if (storedGt.isHomVar()) {
												grid[5]++;
											} else if (storedGt.isHet()) {
												grid[4]++;
											} else {
												grid[3]++;
											}
										} else {
											if (storedGt.isHomVar()) {
												grid[2]++;
											} else if (storedGt.isHet()) {
												grid[1]++;
											} else {
												grid[0]++;
											}
										}
									}
									int[] grid2 = new int[4];
									// grid is the 3x3 grid of genotypes. We now need to convert that into a 2x2 grid, separating the diploid haplotypes if possible. The only case where we cannot work it out is grid[4], where both variants are 0/1.
									// grid2[0] is the ref/ref case. grid[0] (0/0 - 0/0) provides two, while grid[1] (0/0 - 0/1) and grid[3] (0/1 - 0/0) provide one each.
									grid2[0] = 2 * grid[0] + grid[1] + grid[3];
									// grid2[1] is the ref/var case. grid[1] (0/0 - 0/1) provides one, grid[2] (0/0 - 1/1) provides two, and grid[5] (0/1 - 1/1) provides one.
									grid2[1] = grid[1] + 2 * grid[2] + grid[5];
									// grid2[2] is the var/ref case. grid[3] (0/1 - 0/0) provides one, grid[6] (1/1 - 0/0) provides two, and grid[7] (1/1 - 0/1) provides one.
									grid2[2] = grid[3] + 2 * grid[6] + grid[7];
									// grid2[3] is the var/var case. grid[5] (0/1 - 1/1) provides one, grid[7] (1/1 - 0/1) provides one, and grid[8] (1/1 - 1/1) provides two.
									grid2[3] = grid[5] + grid[7] + 2 * grid[8];
									// We can now calculate dprime
									int total = grid2[0] + grid2[1] + grid2[2] + grid2[3];
									double pab = (grid2[3] * 1.0) / total;
									double pa = ((grid2[2] + grid2[3]) * 1.0) / total;
									double pb = ((grid2[1] + grid2[3]) * 1.0) / total;
									double d = pab - pa * pb;
									double dmin = d < 0.0 ? Math.max(-pa * pb, -(1.0 - pa) * (1.0 - pb)) : Math.min(pa * (1.0 - pb), pb * (1.0 - pa));
									double dprime = d / dmin;
									//System.out.println(context.getChr() + "\t" + context.getStart() + "\t" + storedEntry.getKey() + "\t" + (context.getStart() - storedEntry.getKey()) + "\t" + grid2[0] + "\t" + grid2[1] + "\t" + grid2[2] + "\t" + grid2[3] + "\t" + d + "\t" + dmin + "\t" + dprime);
									if ((pa > 0.0) && (pb > 0.0) && (pa < 1.0) && (pb < 1.0)) {
										double rsquared = (d / Math.sqrt(pa * (1.0 - pa) * pb * (1.0 - pb)));
										//System.out.println(context.getChr() + "\t" + context.getStart() + "\t" + storedEntry.getKey() + "\t" + (context.getStart() - storedEntry.getKey()) + "\t" + grid2[0] + "\t" + grid2[1] + "\t" + grid2[2] + "\t" + grid2[3] + "\t" + d + "\t" + dmin + "\t" + dprime + "\t" + rsquared + "\t[" + grid[0] + "," + grid[1] + "," + grid[2] + "," + grid[3] + "," + grid[4] + "," + grid[5] + "," + grid[6] + "," + grid[7] + "," + grid[8] + "]\t" + ((thisBase == otherBase) ? rsquared : -rsquared));
										double rsquaredMinus = (thisBase == otherBase) ? rsquared : -rsquared;
										if (rsquaredMinus < -0.8) {
											viterbi.addSignal(new ViterbiSignal((context.getStart() + storedEntry.getKey()) / 2, (rsquaredMinus + 0.8) * 100.0));
											//System.err.println(currentChr + "\t" + ((context.getStart() + storedEntry.getKey()) / 2) + "\t" + ((rsquaredMinus + 0.8) * 5.0));
										} else if (rsquaredMinus > 0.8) {
											viterbi.addSignal(new ViterbiSignal((context.getStart() + storedEntry.getKey()) / 2, (rsquaredMinus - 0.8) * 5.0));
											//System.err.println(currentChr + "\t" + ((context.getStart() + storedEntry.getKey()) / 2) + "\t" + ((rsquaredMinus - 0.8) * 5.0));
										}
									}
								}
							}
							storedVariants.put(context.getStart(), context);
							storedBases.put(context.getStart(), thisBase);
						} else {
							viterbi.addSignal(new ViterbiSignal(context.getStart(), -6.0));
							//System.err.println(currentChr + "\t" + context.getStart() + "\t-1.1");
						}
					}
				}
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

		public MultiViterbi(String chromosome) {
			viterbis = new ArrayList<Viterbi>();
			for (double transition = 80.0; transition < 10000.0; transition = transition * 1.5) {
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
