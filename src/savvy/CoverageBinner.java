package savvy;
import java.io.ObjectOutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;

import htsjdk.samtools.*;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;

/**
 * Reads a BAM file and counts the number of reads in each 200-bp bin across the entire genome
 *
 * @author Matthew Wakeling
 * 
 * change 2019-08-02 : Pierre Lindenbaum @yokofakun, moving to AbstractApplication, add support for cram
 */
public class CoverageBinner extends AbstractApplication
{
	protected final static Log LOG=Log.getInstance(CoverageBinner.class);

	public static final int DEFAULT_MINIMUM_MQ = 30;
	public static final int BIN_SIZE = 200;

	public static void main(final String[] args) throws Exception {
		new CoverageBinner().instanceMainWithExit(args);
		}
	@Override
	protected Log getLogger() {
		return LOG;
		}
	
	@Override
	public int doWork(List<String> args) throws Exception {
		int min_mapq = DEFAULT_MINIMUM_MQ;
		Path refPath = null;
		Path outPath = null;
		int optind=0;
		while(optind<args.size()) {
			if(args.get(optind).equals("-R") || args.get(optind).equals("--reference")) {
				++optind;
				refPath = Paths.get(args.get(optind));
				}
			else if(args.get(optind).equals("-q") || args.get(optind).equals("--mapq")) {
				++optind;
				min_mapq = Integer.parseInt(args.get(optind));
				}
			else if(args.get(optind).equals("-o") || args.get(optind).equals("--out")) {
				++optind;
				outPath = Paths.get(args.get(optind));
				}
			else if(args.get(optind).equals("-h") || args.get(optind).equals("--help")) {
				System.err.println(" -R|--reference <file> optional reference sequence");
				System.err.println(" -q|--mapq min-mapq <int> ["+DEFAULT_MINIMUM_MQ+"]");
				System.err.println(" -o|--out output <file> file or stdout");
				}
			else if(args.get(optind).equals("--"))
				{
				optind++;
				break;
				}
			else if(args.get(optind).startsWith("-"))
				{
				getLogger().error("unknown option "+args.get(optind));
				return -1;
				}
			else 
				{
				break;
				}
			++optind;
			}
		
		args = args.subList(optind, args.size());
		final Path inputBam = Paths.get(oneAndOnlyOneFile(args));
		
		// Single argument is the name of the BAM file. The vector of integers will be written to the stdout.
		final SamReaderFactory factory = SamReaderFactory.makeDefault();
		if(refPath!=null) factory.referenceSource(new ReferenceSource(refPath));
 		factory.validationStringency(ValidationStringency.LENIENT);
		try(SamReader in = factory.open(inputBam)) {
			final SAMFileHeader header = in.getFileHeader();
			final SAMSequenceDictionary dict = header.getSequenceDictionary();
			if(dict==null || dict.isEmpty()) {
				getLogger().error("no dictionary found in  "+inputBam);
				return -1;
				}
			final List<int[]> tid2vectors = new ArrayList<int[]>();
			
			ProgressLogger progress = new ProgressLogger(getLogger(),10_000_000);
			try(SAMRecordIterator iter = in.iterator()) {
				//LinkedHashMap<String, int[]> vectors = new LinkedHashMap<String, int[]>();
				while (iter.hasNext()) {
					final SAMRecord record = iter.next();
					if(record.getReadUnmappedFlag()) continue;
					if (record.getMappingQuality() < min_mapq) continue;
				
					final int tid = record.getReferenceIndex();
					if(tid <0 ) continue;
					final int start = record.getAlignmentStart();
					final int end = record.getAlignmentEnd();
					final int middle = start + ((end - start) / 2);
					int index = middle / BIN_SIZE;
					
					while(tid2vectors.size()<= tid) tid2vectors.add(null);
					int[] vector = tid2vectors.get(tid);
					if (vector == null) {
						vector = new int[index + 1000];
						tid2vectors.set(tid, vector);
					}
					if (index >= vector.length) {
						final int newSize = Math.max((int) (vector.length * 1.5), index + 1000);
						final int[] newVector = new int[newSize];
						System.arraycopy(vector, 0, newVector, 0, vector.length);
						vector = newVector;
						tid2vectors.set(tid, vector);
						}
					vector[index]++;
					progress.record(record);
					}
				}// end of iterator
			
			final LinkedHashMap<String,int[]> hash = new LinkedHashMap<>(tid2vectors.size());
			for(int i=0;i< dict.size() && i< tid2vectors.size();i++) {
				if(tid2vectors.get(i)==null) continue;
				hash.put(dict.getSequence(i).getSequenceName(), tid2vectors.get(i));
				}
			final ObjectOutputStream out;
			if(outPath==null) {
				out = new ObjectOutputStream(System.out);
				}
			else
				{
				out =  new ObjectOutputStream(Files.newOutputStream(outPath));
				}
			out.writeObject(hash);
			out.flush();
			out.close();	
			}// end of samreader
		 return 0;		
		}
}
