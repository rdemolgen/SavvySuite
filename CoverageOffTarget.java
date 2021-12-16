import java.io.FileInputStream;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Reads a load of CoverageBinner files, calculates the proportion of the genome that is covered by a mean read depth greater than a set amount, then calculates for each file what proportion of the reads are in that part of the genome.
 * Yes, that involves loading each file twice.
 *
 * @author Matthew Wakeling
 */
public class CoverageOffTarget
{
	public static final double DEFAULT_THRESHOLD = 5.0;
	public static final int TOTAL_CHUNKS = 14323926;

	public static void main(String[] args) throws Exception {
		double threshold = DEFAULT_THRESHOLD;
		int readLength = -1;
		Map<String, long[]> totals = new LinkedHashMap<String, long[]>();
		ArrayList<String> inputFiles = new ArrayList<String>();
		for (int argNo = 0; argNo < args.length; argNo++) {
			if ("-readLength".equals(args[argNo])) {
				argNo++;
				readLength = Integer.parseInt(args[argNo]);
			} else if ("-threshold".equals(args[argNo])) {
				argNo++;
				threshold = Double.parseDouble(args[argNo]);
			} else {
				String arg = args[argNo];
				inputFiles.add(arg);
				System.err.println("Loading file " + arg);
				ObjectInputStream in = new ObjectInputStream(new FileInputStream(arg));
				LinkedHashMap<String, int[]> vectors = (LinkedHashMap<String, int[]>) in.readObject();
				for (Map.Entry<String, int[]> entry : vectors.entrySet()) {
					int[] toAdd = entry.getValue();
					long[] total = totals.get(entry.getKey());
					if (total == null) {
						total = new long[(int) (toAdd.length * 1.1)];
						totals.put(entry.getKey(), total);
					} else {
						if (toAdd.length > total.length) {
							long[] swap = new long[(int) (toAdd.length * 1.1)];
							for (int i = 0; i < total.length; i++) {
								swap[i] = total[i];
							}
							total = swap;
							totals.put(entry.getKey(), swap);
						}
					}
					for (int i = 0; i < toAdd.length; i++) {
						total[i] += toAdd[i];
					}
				}
				in.close();
			}
		}
		// Loaded all the files, and stored the total. We can now use this to determine which parts of the genome are targeted.
		int chunksCovered = 0;
		for (Map.Entry<String, long[]> entry : totals.entrySet()) {
			long[] totalArray = entry.getValue();
			for (int i = 0; i < totalArray.length; i++) {
				if (totalArray[i] > threshold * inputFiles.size()) {
					chunksCovered++;
				}
			}
		}
		// Assume that the human genome has about 2864785223 (from GRCh37) mappable bases. Chunk size is 200. Therefore total number of chunks is 14323926
		int chunksOffTarget = TOTAL_CHUNKS - chunksCovered;
		System.out.println("Chunks covered >X" + threshold + ": " + chunksCovered + "\tProportion of the genome covered >X" + threshold + ": " + ((1.0 * chunksCovered) / TOTAL_CHUNKS));
		if (readLength > -1) {
			System.out.println("Sample\ttotal_reads\treads_in_chunks_>X" + threshold + "\tread_proportion_in_>X" + threshold + "\toff_target_reads\toff_target_read_depth");
		} else {
			System.out.println("Sample\ttotal_reads\treads_in_chunks_>X" + threshold + "\tread_proportion_in_>X" + threshold + "\toff_target_reads");
		}
		for (String arg : inputFiles) {
			long totalReads = 0;
			long targetedReads = 0;
			ObjectInputStream in = new ObjectInputStream(new FileInputStream(arg));
			LinkedHashMap<String, int[]> vectors = (LinkedHashMap<String, int[]>) in.readObject();
			for (Map.Entry<String, int[]> entry : vectors.entrySet()) {
				int[] toAdd = entry.getValue();
				long[] total = totals.get(entry.getKey());
				for (int i = 0; i < toAdd.length; i++) {
					if (total[i] > threshold * inputFiles.size()) {
						targetedReads += toAdd[i];
					}
					totalReads += toAdd[i];
				}
			}
			System.out.println(arg + "\t" + totalReads + "\t" + targetedReads + "\t" + ((1.0 * targetedReads) / totalReads) + "\t" + (totalReads - targetedReads) + "\t" + (readLength > -1 ? "\t" + ((totalReads - targetedReads) * readLength / 200.0 / chunksOffTarget): ""));
		}
	}
}