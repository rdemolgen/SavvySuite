import java.io.File;
import java.io.ObjectOutputStream;
import java.util.LinkedHashMap;
import htsjdk.samtools.*;

/**
 * Reads a BAM file and counts the number of reads in each 200-bp bin across the entire genome
 *
 * @author Matthew Wakeling
 */
public class CoverageBinner
{
	public static final int MINIMUM_MQ = 30;
	public static final int BIN_SIZE = 200;

	public static void main(String[] args) throws Exception {
		// Single argument is the name of the BAM file. The vector of integers will be written to the stdout.
		SamReaderFactory factory = SamReaderFactory.makeDefault();
		factory.validationStringency(ValidationStringency.LENIENT);
		SamReader in = factory.open(new File(args[0]));
		SAMRecordIterator iter = in.iterator();
		LinkedHashMap<String, int[]> vectors = new LinkedHashMap<String, int[]>();
		int progress = 0;
		int resizeCount = 0;
		while (iter.hasNext()) {
			SAMRecord record = iter.next();
			if (record.getMappingQuality() >= MINIMUM_MQ) {
				String chr = record.getReferenceName();
				int start = record.getAlignmentStart();
				int end = record.getAlignmentEnd();
				int middle = start + ((end - start) / 2);
				int index = middle / BIN_SIZE;
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
		ObjectOutputStream out = new ObjectOutputStream(System.out);
		out.writeObject(vectors);
		out.flush();
		out.close();
	}
}
