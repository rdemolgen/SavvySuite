import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.util.zip.ZipException;
import htsjdk.samtools.*;

/**
 * Reads a BAM file and counts the number of reads in each 200-bp bin across the entire genome
 *
 * @author Matthew Wakeling
 */
public class CoverageBinner
{
	public static final int DEFAULT_MINIMUM_MQ = 30;
	public static final int DEFAULT_BIN_SIZE = 200;

	public static final Pattern DOUBLE_CLIP = Pattern.compile("[0-9]*[SH].*[SH]");

	public static void main(String[] args) throws Exception {
		if ("info".equals(args[0])) {
			CoverageBinner in = new CoverageBinner(args[1]);
			for (String key : in.metadata.keySet()) {
				System.out.println(key + ": " + in.metadata.get(key));
			}
			System.exit(0);
		} else if ("dump".equals(args[0])) {
			CoverageBinner in = new CoverageBinner(args[1]);
			int binSize = Integer.parseInt(args[3]);
			Map<String, int[]> vectors = in.getVectors(binSize);
			int[] vector = vectors.get(args[2]);
			if (vector == null) {
				System.err.println("No chromosome " + args[2] + " in file");
				System.exit(1);
			}
			for (int i = 0; i < vector.length; i++) {
				System.out.println(i + "\t" + (i * binSize) + "\t" + vector[i]);
			}
			System.exit(0);
		}
		String inputFileName = null;
		String referenceFileName = null;
		int binSize = DEFAULT_BIN_SIZE;
		int minMq = DEFAULT_MINIMUM_MQ;
		boolean includeDoubleClip = false;
		boolean includeSecondary = true;
		String sampleName = null;
		for (int argNo = 0; argNo < args.length; argNo++) {
			if ("-d".equals(args[argNo])) {
				argNo++;
				binSize = Integer.parseInt(args[argNo]);
			} else if ("-mmq".equals(args[argNo])) {
				argNo++;
				minMq = Integer.parseInt(args[argNo]);
			} else if ("-includeDoubleClip".equals(args[argNo])) {
				includeDoubleClip = true;
			} else if ("-excludeSecondary".equals(args[argNo])) {
				includeSecondary = false;
			} else if ("-R".equals(args[argNo])) {
				argNo++;
				referenceFileName = args[argNo];
			} else if ("-s".equals(args[argNo])) {
				argNo++;
				sampleName = args[argNo];
			} else if (inputFileName == null) {
				inputFileName = args[argNo];
			} else {
				System.err.println("Could not understand option \"" + args[argNo] + "\"");
				System.err.println("Usage: java CoverageBinner [options] input.sam/bam/cram >output.coverageBinner");
				System.err.println("Options:");
				System.err.println("   -d n                Change the bin size from the default of 200");
				System.err.println("   -mmq n              Change the minimum mapping quality from the default of 30");
				System.err.println("   -R ref.fasta        Location of the reference genome FASTA file (only required for CRAM files)");
				System.err.println("   -includeDoubleClip  Include double-clipped reads - these are excluded by default as they are often misaligned reads from contamination from a different species, for instance bacterial contamination of a saliva DNA sample.");
				System.err.println("   -excludeSecondary   Exclude reads that are marked as secondary or supplementary");
				System.err.println("   -s sampleName       Record the given sample name instead of the one in the BAM file");
				System.err.println("Usage: java CoverageBinner info input.coverageBinner");
				System.err.println("Usage: java CoverageBinner dump input.coverageBinner chromosome binSize");
				System.exit(1);
			}
		}
		SamReaderFactory factory = SamReaderFactory.makeDefault();
		factory.validationStringency(ValidationStringency.LENIENT);
		if (referenceFileName != null) {
			factory.referenceSequence(new File(referenceFileName));
		}
		SamReader in = factory.open(new File(inputFileName));
		SAMFileHeader header = in.getFileHeader();
		if ((sampleName == null) && (header != null)) {
			List<SAMReadGroupRecord> readGroups = header.getReadGroups();
			if ((readGroups != null) && (!readGroups.isEmpty())) {
				sampleName = readGroups.get(0).getSample();
			}
		}
		SAMRecordIterator iter = in.iterator();
		LinkedHashMap<String, int[]> vectors = new LinkedHashMap<String, int[]>();
		int progress = 0;
		int resizeCount = 0;
		while (iter.hasNext()) {
			SAMRecord record = iter.next();
			String cigar = record.getCigarString();
			Matcher m = DOUBLE_CLIP.matcher(cigar);
			if ((record.getMappingQuality() >= minMq) && ((!m.matches()) || includeDoubleClip) && ((!record.isSecondaryOrSupplementary()) || includeSecondary)) {
				String chr = record.getReferenceName();
				int start = record.getAlignmentStart();
				int end = record.getAlignmentEnd();
				int middle = start + ((end - start) / 2);
				int index = middle / binSize;
				int[] vector = vectors.get(chr);
				if (vector == null) {
					vector = new int[index + 1000];
					vectors.put(chr, vector);
				}
				if (index >= vector.length) {
					int newSize = Math.max((int) (vector.length * 1.5), index + 1000);
					int[] newVector = new int[newSize];
					for (int i = 0; i < vector.length; i++) {
						newVector[i] = vector[i];
					}
					vector = newVector;
					vectors.put(chr, vector);
					resizeCount++;
				}
				vector[index]++;
				if ((++progress) % 10000000 == 0) {
					System.err.println("Processed " + progress + " reads, resized " + resizeCount + " times, now on " + chr + ":" + start);
				}
			}
		}
		LinkedHashMap<String, byte[]> byteArrays = new LinkedHashMap<String, byte[]>();
		for (Map.Entry<String, int[]> entry : vectors.entrySet()) {
			int[] vector = entry.getValue();
			int size = 0;
			for (int i = 0; i < vector.length; i++) {
				if (vector[i] != 0) {
					size = i;
				}
			}
			size++;
			byte[] byteArray = new byte[(size + 1) * 5];
			byteArray[0] = (byte) (size >> 24);
			byteArray[1] = (byte) (size >> 16);
			byteArray[2] = (byte) (size >> 8);
			byteArray[3] = (byte) size;
			int pos = 4;
			for (int i = 0; i < size; i++) {
				int v = vector[i];
				if (v < 128) {
					byteArray[pos++] = (byte) v;
				} else if (v < 16384) {
					byteArray[pos++] = (byte) (128 | (v >> 8));
					byteArray[pos++] = (byte) v;
				} else if (v < 2097152) {
					byteArray[pos++] = (byte) (192 | (v >> 16));
					byteArray[pos++] = (byte) (v >> 8);
					byteArray[pos++] = (byte) v;
				} else {
					byteArray[pos++] = (byte) 224;
					byteArray[pos++] = (byte) (v >> 24);
					byteArray[pos++] = (byte) (v >> 16);
					byteArray[pos++] = (byte) (v >> 8);
					byteArray[pos++] = (byte) v;
				}
			}
			byte[] newByteArray = new byte[pos];
			for (int i = 0; i < pos; i++) {
				newByteArray[i] = byteArray[i];
			}
			byteArrays.put(entry.getKey(), newByteArray);
		}
		ObjectOutputStream out = new ObjectOutputStream(new GZIPOutputStream(System.out));
		Map<String, Object> metadata = new LinkedHashMap<String, Object>();
		metadata.put("d", binSize);
		metadata.put("mmq", minMq);
		metadata.put("doubleclip", includeDoubleClip);
		metadata.put("secondary", includeSecondary);
		metadata.put("sourceFile", inputFileName);
		metadata.put("sample", sampleName);
		out.writeObject(metadata);
		out.writeObject(byteArrays);
		out.flush();
		out.close();
	}

	private LinkedHashMap<String, Object> metadata;
	private LinkedHashMap<String, int[]> vectors;

	@SuppressWarnings("unchecked")
	public CoverageBinner(String fileName) throws IOException {
		ObjectInputStream in = null;
		try {
			in = new ObjectInputStream(new BufferedInputStream(new GZIPInputStream(new FileInputStream(fileName))));
			try {
				metadata = (LinkedHashMap<String, Object>) in.readObject();
				metadata.put("format", "new");
				Map<String, byte[]> byteArrays = (Map<String, byte[]>) in.readObject();
				vectors = new LinkedHashMap<String, int[]>();
				for (String chr : byteArrays.keySet()) {
					byte[] byteArray = byteArrays.get(chr);
					int size = ((((int) byteArray[0]) & 255) << 24) | ((((int) byteArray[1]) & 255) << 16) | ((((int) byteArray[2]) & 255) << 8) | (((int) byteArray[3]) & 255);
					int[] vector = new int[size];
					int pos = 4;
					for (int i = 0; i < size; i++) {
						int v = ((int) byteArray[pos++]) & 255;
						if (v < 128) {
							vector[i] = v;
						} else if (v < 192) {
							int b2 = ((int) byteArray[pos++]) & 255;
							vector[i] = ((v & 63) << 8) | b2;
						} else if (v < 224) {
							int b2 = ((int) byteArray[pos++]) & 255;
							int b3 = ((int) byteArray[pos++]) & 255;
							vector[i] = ((v & 31) << 16) | (b2 << 8) | b3;
						} else if (v == 224) {
							int b2 = ((int) byteArray[pos++]) & 255;
							int b3 = ((int) byteArray[pos++]) & 255;
							int b4 = ((int) byteArray[pos++]) & 255;
							int b5 = ((int) byteArray[pos++]) & 255;
							vector[i] = (b2 << 24) | (b3 << 16) | (b4 << 8) | b5;
						} else {
							throw new IOException("Unexpected first byte " + v);
						}
					}
					if (pos != byteArray.length) {
						throw new IOException("Decoded " + pos + " bytes but data is " + byteArray.length + " bytes long");
					}
					vectors.put(chr, vector);
				}
			} catch (ClassNotFoundException e) {
				throw new IOException(e);
			} finally {
				if (in != null) {
					in.close();
				}
			}
		} catch (ZipException e) {
			in = new ObjectInputStream(new BufferedInputStream(new FileInputStream(fileName)));
			try {
				vectors = (LinkedHashMap<String, int[]>) in.readObject();
			} catch (ClassNotFoundException e2) {
				throw new IOException(e2);
			} finally {
				if (in != null) {
					in.close();
				}
			}
			metadata = new LinkedHashMap<String, Object>();
			metadata.put("d", 200);
			metadata.put("mmq", 30);
			metadata.put("doubleclip", true);
			metadata.put("secondary", true);
			metadata.put("sourceFile", fileName);
			metadata.put("sample", fileName);
			metadata.put("format", "old");
		}
		metadata.put("coverageFile", fileName);
	}

	public int getBinSize() {
		return (Integer) metadata.get("d");
	}

	public int getMinMq() {
		return (Integer) metadata.get("mmq");
	}

	public boolean getDoubleClipIncluded() {
		return (Boolean) metadata.get("doubleclip");
	}

	public boolean getSecondaryIncluded() {
		return (Boolean) metadata.get("secondary");
	}

	public String getSourceFile() {
		return (String) metadata.get("sourceFile");
	}

	public String getSampleName() {
		return (String) metadata.get("sample");
	}

	public boolean isNewFormat() {
		return "new".equals(metadata.get("format"));
	}

	public LinkedHashMap<String, int[]> getVectors(int binSize) throws IOException {
		int fileBinSize = getBinSize();
		if (binSize % fileBinSize != 0) {
			throw new IOException("Cannot use a bin size of " + binSize + " because the file " + metadata.get("coverageFile") + " has a bin size of " + fileBinSize + ". Requested bin size must be a multiple of " + fileBinSize + ".");
		}
		if (binSize == fileBinSize) {
			return vectors;
		}
		LinkedHashMap<String, int[]> retval = new LinkedHashMap<String, int[]>();
		int multiplier = binSize / fileBinSize;
		for (String chr : vectors.keySet()) {
			int[] vector = vectors.get(chr);
			int[] newVector = new int[(vector.length + multiplier - 1) / multiplier];
			for (int i = 0; i < vector.length; i++) {
				newVector[i / multiplier] += vector[i];
			}
			retval.put(chr, newVector);
		}
		return retval;
	}
}
