import htsjdk.samtools.*;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.*;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.OutputStream;
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
import java.util.zip.GZIPOutputStream;

public class PrepareLinkageData
{
	public static final int WIDTH = 100000;
	public static final Pattern INDEL = Pattern.compile(".*[ACGT][ACGT].*");
	/*
	 * Method - We want to find homozygous regions using single off-target reads from the BAM file, using linkage disequilibrium.
	 * The files are big, so we need to ensure performance is sensible.
	 * This program takes a VCF file containing linkage information, removes variants that do not have linkage, are not biallelic, or are not SNVs, and writes a binary-format file containing only the relevant data.
	 * This improves the performance of SavvyHomozygosity and SavvySharedHaplotypes, as a lot of the time spent is loading the linkage data.
	 * Procedure:
	 *
	 * 1. Load a variant from the VCF. Check to see if it is biallelic, is an SNV, and has an appropriate frequency. Discard if not.
	 * 2. Otherwise, iterate through all variants stored to the left.
	 *   2.1. If the variant is more than WIDTH to the left of the current variant, then discard it.
	 *   2.2. Otherwise, calculate linkage from VCF variants, produce score from reads.
	 */

	/**
	 * Process a VCF file containing whole genome genotypes for a large number of samples. The first argument is a VCF file containing common variants for lots of WGS samples.
	 *
	 * @author Matthew Wakeling
	 */
	public static void main(String[] args) throws IOException {
		VCFFileReader vcf = new VCFFileReader(new File(args[0]));
		boolean allVariants = (args.length > 1 ? "-all".equals(args[1]) : false);
		List<String> vcfSampleNames = vcf.getFileHeader().getGenotypeSamples();
		TreeMap<Integer, VariantArray> storedVariants = null;
		String currentChr = "";
		Output output = new Output(System.out);
		for (VariantContext context : vcf) {
			@SuppressWarnings("deprecation") String contextChr = context.getChr();
			if (!(currentChr.equals(contextChr))) {
				storedVariants = new TreeMap<Integer, VariantArray>();
				currentChr = contextChr;
			}
			String alleles = "" + context.getAlleles();
			boolean indelMatch = INDEL.matcher(alleles).matches();
			if ((!indelMatch) && (context.getAlleles().size() == 2)) {
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
				if ((ac * 8 > an) && (ac * 8 < an * 7)) {
					VariantArray vArray = new VariantArray(vcfSampleNames, context);
					if (allVariants) {
						output.add(vArray);
					} else {
						Iterator<Map.Entry<Integer, VariantArray>> storedIter = storedVariants.entrySet().iterator();
						while (storedIter.hasNext()) {
							Map.Entry<Integer, VariantArray> storedEntry = storedIter.next();
							if (storedEntry.getKey() < context.getStart() - WIDTH) {
								storedIter.remove();
							} else {
								double rsquared = vArray.getRsquared(storedEntry.getValue());
								if ((rsquared < -0.8) || (rsquared > 0.8)) {
									output.add(vArray);
									output.add(storedEntry.getValue());
								}
							}
						}
						storedVariants.put(context.getStart(), vArray);
					}
				}
			}
		}
		output.finish();
		output.close();
	}

	public static class Output
	{
		private TreeMap<Integer, VariantArray> variantsInFlight = new TreeMap<Integer, VariantArray>();
		private ObjectOutputStream out;
		private String chr = null;
		private static int count = 0;

		public Output(OutputStream o) throws IOException {
			this.out = new ObjectOutputStream(new BufferedOutputStream(new GZIPOutputStream(o)));
		}

		public void add(VariantArray context) throws IOException {
			@SuppressWarnings("deprecation") String chromosome = context.getChr();
			if (!chromosome.equals(chr)) {
				finish();
				chr = chromosome;
			}
			if (!variantsInFlight.containsKey(context.getPosition())) {
				variantsInFlight.put(context.getPosition(), context);
			}
			while (variantsInFlight.firstKey() < context.getPosition() - WIDTH) {
				outputContext(variantsInFlight.remove(variantsInFlight.firstKey()));
			}
		}

		public void finish() throws IOException {
			for (Map.Entry<Integer, VariantArray> entry : variantsInFlight.entrySet()) {
				outputContext(entry.getValue());
			}
			variantsInFlight.clear();
			out.flush();
		}

		public void outputContext(VariantArray context) throws IOException {
			if ((++count) % 10000 == 0) {
				System.err.println("Written " + count + " variants, now on " + context.getChr() + ":" + context.getPosition());
			}
			out.writeObject(context);
			out.reset();
		}

		public void close() throws IOException {
			out.flush();
			out.close();
		}
	}
}
