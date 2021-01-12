import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.*;
import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.EOFException;
import java.io.IOException;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.ObjectInputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipException;

/**
 * Reads a list of CNVs to validate (such as produced by GenomeStrip2 LCNV pipeline), and calculates the average allele balance parameters for each one, from the VCF files provided.
 * The first argument should be a file containing CNV calls. The first four columns should be chromosome, start, end, and sample name.
 * The second argument should be a genotype file created by PrepareLinkageData.
 * The remaining arguments should be vcf file names.
 * The original contents of the file will be output, with four extra columns showing number of variants in the CNV, and the popularity of 25%, 33%, and 50% allele ratio variants.
 * With a heterozygous deletion, all popularities will be low. For a heterozygous duplication, the 33% value will be highest. For a homozygous duplication, the 25% ratio will be highest. For a homozygous deletion, the number of variants will be low.
 *
 * @author Matthew Wakeling
 */
public class ValidateGenomeCnvsUsingAlleleBalance
{
	public static final Pattern TAB = Pattern.compile("\t");

	public static void main(String[] args) throws IOException, InterruptedException, EOFException, ClassNotFoundException {
		int argNo = 0;
		int chromColumn = 0;
		int startColumn = 1;
		int endColumn = 2;
		int sampleNameColumn = 3;
		if ("-columns".equals(args[argNo])) {
			argNo++;
			chromColumn = Integer.parseInt(args[argNo]) - 1;
			argNo++;
			startColumn = Integer.parseInt(args[argNo]) - 1;
			argNo++;
			endColumn = Integer.parseInt(args[argNo]) - 1;
			argNo++;
			sampleNameColumn = Integer.parseInt(args[argNo]) - 1;
			argNo++;
		}
		Map<String, List<Cnv>> cnvs = new HashMap<String, List<Cnv>>();
		BufferedReader in = new BufferedReader(new FileReader(args[argNo++]));
		String line = in.readLine();
		System.out.println(line + "\tVariants\t25%\t33%\t50%");
		line = in.readLine();
		while (line != null) {
			try {
				String[] split = TAB.split(line);
				List<Cnv> list = cnvs.get(split[sampleNameColumn]);
				if (list == null) {
					list = new ArrayList<Cnv>();
					cnvs.put(split[sampleNameColumn], list);
				}
				list.add(new Cnv(split[chromColumn], Integer.parseInt(split[startColumn]), Integer.parseInt(split[endColumn]), line));
			} catch (NumberFormatException e) {
			}
			line = in.readLine();
		}
		in.close();
		String referenceName = args[argNo++];
		Map<String, Map<Integer, VariantArray>> ref = null;
		HardyWeinbergLimit hardy = null;
		if (!("none".equals(referenceName))) {
			ref = loadReference(referenceName);
			hardy = new HardyWeinbergLimit(ref.values().iterator().next().values().iterator().next().getSamples().length, 0.0001);
		}
		for (; argNo < args.length; argNo++) {
			String vcfFile = args[argNo];
			System.err.println("Loading " + vcfFile);
			VCFFileReader reader = new VCFFileReader(new File(vcfFile));
			List<String> sampleNames = reader.getFileHeader().getGenotypeSamples();
			for (String sampleName : sampleNames) {
				List<Cnv> cnvList = cnvs.get(sampleName);
				if (cnvList == null) {
					cnvList = cnvs.get("all");
				}
				if (cnvList != null) {
					for (Cnv cnv : cnvList) {
						double sum25 = 0.0;
						double sum33 = 0.0;
						double sum50 = 0.0;
						double count = 0.0;
						int variants = 0;
						int hetVariants = 0;
						int totalReads = 0;
						Iterator<VariantContext> iterator = reader.query(cnv.getChr(), cnv.getStart(), cnv.getEnd());
						List<int[]> ads = new ArrayList<int[]>();
						while (iterator.hasNext()) {
							VariantContext context = iterator.next();
							boolean hardyPass = true;
							if (ref != null) {
								hardyPass = false;
								Map<Integer, VariantArray> chrRef = ref.get(context.getChr());
								if (chrRef != null) {
									VariantArray vArray = chrRef.get(context.getStart());
									if (vArray != null) {
										int[] rcounts = new int[3];
										for (byte b : vArray.getSamples()) {
											rcounts[b]++;
										}
										if ("X".equals(context.getChr()) || hardy.pass(rcounts[0], rcounts[1], rcounts[2])) {
											hardyPass = true;
										}
									}
								}
							}
							if (hardyPass) {
								Genotype g = context.getGenotype(sampleName);
								int[] ad = g.getAD();
								if ((ad != null) && (!GenotypeType.HOM_REF.equals(g.getType())) && (ad.length == 2) && (ad[0] + ad[1] > 0)) {
									variants++;
									if (GenotypeType.HET.equals(g.getType())) {
										ads.add(ad);
										hetVariants++;
										int high = ad[0] + ad[1];
										totalReads += high;
										count += high * high;
										sum25 += (calculateScore(ad[0], ad[1], 0.25) + calculateScore(ad[0], ad[1], 0.75)) * high * high;
										sum33 += (calculateScore(ad[0], ad[1], 0.3333) + calculateScore(ad[0], ad[1], 0.6667)) * high * high;
										sum50 += 2.0 * calculateScore(ad[0], ad[1], 0.5) * high * high;
									}
								}
							}
						}
						double low = 0.0;
						double high = 0.5;
						double step = 0.05;
						double firstPeak = 0.0;
						double peakHeight = 0.0;
						while (step > 0.00001) {
							peakHeight = 0.0;
							boolean fallen = false;
							for (double i = low; (i < high) && (!fallen); i += step) {
								double h = 0.0;
								for (int[] ad : ads) {
									h += (calculateScore(ad[0], ad[1], i) + calculateScore(ad[0], ad[1], 1.0 - i)) * (ad[0] + ad[1]) * (ad[0] + ad[1]);
								}
								if (h >= peakHeight) {
									peakHeight = h;
									firstPeak = i;
//								} else {
//									fallen = true;
								}
							}
							low = Math.max(0.0, firstPeak - step);
							high = Math.min(0.5, firstPeak + step);
							step = step / 8;
						}
						System.out.println(cnv.getLine() + "\t" + sampleName + "\t" + variants + "\t" + hetVariants + "\t" + totalReads + "\t" + (sum25 / count) + "\t" + (sum33 / count) + "\t" + (sum50 / count) + "\t" + ((hetVariants * 1.0) / variants) + "\t" + firstPeak + "\t" + (peakHeight / count));
					}
				}
			}
			reader.close();
		}
	}

	public static double calculateScore(int ad0, int ad1, double frac) {
		double mfrac = 1.0 - frac;
		int low = ad0 > ad1 ? ad1 : ad0;
		int high = ad0 + ad1;
		double mult = 1.0;
		int div1 = 2;
		int div2 = low;
		int div3 = high - low;
		for (int i = high - low + 1; i <= high + 1; i++) {
			mult = mult * i;
			while ((mult > 1000000.0) && (div1 <= low)) {
				mult = mult / div1;
				div1++;
			}
			while ((mult > 1000000.0) && (div2 > 0)) {
				mult = mult * frac;
				div2--;
			}
			while ((mult > 1000000.0) && (div3 > 0)) {
				mult = mult * mfrac;
				div3--;
			}
		}
		while (div1 <= low) {
			mult = mult / div1;
			div1++;
		}
		mult = mult * Math.pow(frac, div2);
		mult = mult * Math.pow(mfrac, div3);
		return mult;
	}

	public static class Cnv
	{
		private String chr, line;
		private int start, end;

		public Cnv(String chr, int start, int end, String line) {
			this.chr = chr;
			this.start = start;
			this.end = end;
			this.line = line;
		}

		public String getChr() {
			return chr;
		}

		public int getStart() {
			return start;
		}

		public int getEnd() {
			return end;
		}

		public String getLine() {
			return line;
		}
	}

	public static Map<String, Map<Integer, VariantArray>> loadReference(String fileName) throws IOException, ClassNotFoundException {
		ObjectInputStream vcf = null;
		try {
			vcf = new ObjectInputStream(new BufferedInputStream(new GZIPInputStream(new FileInputStream(fileName))));
		} catch (ZipException e) {
			vcf = new ObjectInputStream(new BufferedInputStream(new FileInputStream(fileName)));
		}
		HashMap<String, Map<Integer, VariantArray>> ref = new HashMap<String, Map<Integer, VariantArray>>();
		boolean hasMoreVariants = true;
		VariantArray vArray = null;
		try {
			vArray = (VariantArray) vcf.readObject();
		} catch (EOFException e) {
			hasMoreVariants = false;
		}
		while (hasMoreVariants) {
			Map<Integer, VariantArray> chrRef = ref.get(vArray.getChr());
			if (chrRef == null) {
				chrRef = new HashMap<Integer, VariantArray>();
				ref.put(vArray.getChr(), chrRef);
			}
			chrRef.put(vArray.getPosition(), vArray);
			try {
				vArray = (VariantArray) vcf.readObject();
			} catch (EOFException e) {
				hasMoreVariants = false;
			}
		}
		return ref;
	}

	public static class HardyWeinbergLimit
	{
		private int samples;
		private int[] lowLimit;
		private int[] highLimit;

		public HardyWeinbergLimit(int samples, double pValue) {
			this.samples = samples;
			this.lowLimit = new int[samples + 1];
			this.highLimit = new int[samples + 1];
			for (int homVar = 0; homVar <= samples; homVar++) {
				while (FishersExactTest.hardyS(samples - homVar - lowLimit[homVar], lowLimit[homVar], homVar) < pValue) {
					lowLimit[homVar]++;
				}
				highLimit[homVar] = samples - homVar;
				while (FishersExactTest.hardyS(samples - homVar - highLimit[homVar], highLimit[homVar], homVar) < pValue) {
					highLimit[homVar]--;
				}
				//System.out.println(samples + "\t" + homVar + "\t" + lowLimit[homVar] + "\t" + highLimit[homVar]);
			}
		}

		public boolean pass(int homRef, int het, int homVar) {
			return (lowLimit[homVar] <= het) && (highLimit[homVar] >= het);
		}
	}
}
