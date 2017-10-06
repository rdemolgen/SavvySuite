import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.Set;
import java.util.TreeSet;

/**
 * Object representing an interval on a chromosome.
 *
 * @author Matthew Wakeling
 */
public class Interval<T> implements Comparable<Interval<?>>
{
	private String chromosome;
	private int start, end;
	private T data;

	public Interval(String chromosome, int start, int end, T data) {
		this.chromosome = chromosome;
		this.start = start;
		this.end = end;
		this.data = data;
	}

	public String getChromosome() {
		return chromosome;
	}

	public int getStart() {
		return start;
	}

	public int getEnd() {
		return end;
	}

	public T getData() {
		return data;
	}

	public void setData(T data) {
		this.data = data;
	}

	public int compareTo(Interval<?> i) {
		if (!chromosome.equals(i.chromosome)) {
			int num = Integer.MAX_VALUE;
			int inum = Integer.MAX_VALUE;
			try {
				num = Integer.parseInt(chromosome);
			} catch (NumberFormatException e) {
				// Not a number
			}
			try {
				inum = Integer.parseInt(i.chromosome);
			} catch (NumberFormatException e) {
				// Not a number
			}
			if (num > inum) {
				return 1;
			} else if (num < inum) {
				return -1;
			}
			return chromosome.compareTo(i.chromosome);
		} else if (start > i.start) {
			return 1;
		} else if (start < i.start) {
			return -1;
		} else {
			return 0;
		}
	}

	public String toString() {
		return chromosome + ":" + start + "-" + end + " " + data;
	}

	public static void addInterval(TreeSet<Interval<Set<String>>> intervals, Interval<Set<String>> interval) {
		//System.err.println("Adding " + interval);
		Interval<Set<String>> before = intervals.lower(interval);
		if ((before == null) && (!intervals.isEmpty())) {
			before = intervals.first();
		}
		String chr = interval.getChromosome();
		int start = interval.getStart();
		int end = interval.getEnd();
		Interval<Set<String>> endInterval = new Interval<Set<String>>(chr, end, 0, null);
		while ((before != null) && (before.compareTo(endInterval) <= 0)) {
			if (before.getChromosome().equals(chr) && (before.getEnd() >= start)) {
				intervals.remove(before);
				if (start > before.getStart()) {
					intervals.add(new Interval<Set<String>>(chr, before.getStart(), start - 1, before.getData()));
				}
				if (start < before.getStart()) {
					intervals.add(new Interval<Set<String>>(chr, start, before.getStart() - 1, interval.getData()));
					start = before.getStart();
				}
				Set<String> strings = new HashSet<String>();
				strings.addAll(interval.getData());
				strings.addAll(before.getData());
				if (end >= before.getEnd()) {
					intervals.add(new Interval<Set<String>>(chr, start, before.getEnd(), strings));
					start = before.getEnd() + 1;
				} else {
					intervals.add(new Interval<Set<String>>(chr, start, end, strings));
					intervals.add(new Interval<Set<String>>(chr, end + 1, before.getEnd(), before.getData()));
					start = end + 1;
				}
			}
			before = intervals.higher(before);
		}
		if (start <= end) {
			intervals.add(new Interval<Set<String>>(chr, start, end, interval.getData()));
		}
		//System.err.println("Intervals: " + intervals);
	}

	public static void loadTarget(TreeSet<Interval<Set<String>>> union, String targetFile, String targetName) throws IOException {
		//System.err.println("Loading target intervals \"" + targetFile + "\"");
		BufferedReader in = new BufferedReader(new FileReader(targetFile));
		String line = in.readLine();
		while (line != null) {
			try {
				if (!line.startsWith("@")) {
					int tabPos = line.indexOf('\t');
					if (tabPos > -1) {
						String chr = line.substring(0, tabPos);
						line = line.substring(tabPos + 1);
						tabPos = line.indexOf('\t');
						if (tabPos > -1) {
							int start = Integer.parseInt(line.substring(0, tabPos));
							line = line.substring(tabPos + 1);
							tabPos = line.indexOf('\t');
							if (tabPos > -1) {
								int end = Integer.parseInt(line.substring(0, tabPos));
								Set<String> names = new HashSet<String>();
								names.add("zzz" + targetName + "\t");
								addInterval(union, new Interval<Set<String>>(chr, start, end, names));
							} else {
								int end = Integer.parseInt(line);
								Set<String> names = new HashSet<String>();
								names.add("zzz" + targetName + "\t");
								addInterval(union, new Interval<Set<String>>(chr, start, end, names));
							}
						}
					}
				}
			} catch (NumberFormatException e) {
			}
			line = in.readLine();
		}
		in.close();
	}

	public static void writeBed(TreeSet<Interval<Set<String>>> intervals, PrintStream out) {
		for (Interval i : intervals) {
			System.out.println(i.getChromosome() + "\t" + (i.getStart() - 1) + "\t" + i.getEnd() + "\t+\t");
		}
	}

	public static void writeIntervalList(TreeSet<Interval<Set<String>>> intervals, PrintStream out) {
		out.println("@SQ\tSN:1\tLN:249250621\t\t");
		out.println("@SQ\tSN:2\tLN:243199373\t\t");
		out.println("@SQ\tSN:3\tLN:198022430\t\t");
		out.println("@SQ\tSN:4\tLN:191154276\t\t");
		out.println("@SQ\tSN:5\tLN:180915260\t\t");
		out.println("@SQ\tSN:6\tLN:171115067\t\t");
		out.println("@SQ\tSN:7\tLN:159138663\t\t");
		out.println("@SQ\tSN:8\tLN:146364022\t\t");
		out.println("@SQ\tSN:9\tLN:141213431\t\t");
		out.println("@SQ\tSN:10\tLN:135534747\t\t");
		out.println("@SQ\tSN:11\tLN:135006516\t\t");
		out.println("@SQ\tSN:12\tLN:133851895\t\t");
		out.println("@SQ\tSN:13\tLN:115169878\t\t");
		out.println("@SQ\tSN:14\tLN:107349540\t\t");
		out.println("@SQ\tSN:15\tLN:102531392\t\t");
		out.println("@SQ\tSN:16\tLN:90354753\t\t");
		out.println("@SQ\tSN:17\tLN:81195210\t\t");
		out.println("@SQ\tSN:18\tLN:78077248\t\t");
		out.println("@SQ\tSN:19\tLN:59128983\t\t");
		out.println("@SQ\tSN:20\tLN:63025520\t\t");
		out.println("@SQ\tSN:21\tLN:48129895\t\t");
		out.println("@SQ\tSN:22\tLN:51304566\t\t");
		out.println("@SQ\tSN:X\tLN:155270560\t\t");
		out.println("@SQ\tSN:Y\tLN:59373566\t\t");
		out.println("@SQ\tSN:MT\tLN:16569\t\t");
		out.println("@SQ\tSN:GL000207.1\tLN:4262\t\t");
		out.println("@SQ\tSN:GL000226.1\tLN:15008\t\t");
		out.println("@SQ\tSN:GL000229.1\tLN:19913\t\t");
		out.println("@SQ\tSN:GL000231.1\tLN:27386\t\t");
		out.println("@SQ\tSN:GL000210.1\tLN:27682\t\t");
		out.println("@SQ\tSN:GL000239.1\tLN:33824\t\t");
		out.println("@SQ\tSN:GL000235.1\tLN:34474\t\t");
		out.println("@SQ\tSN:GL000201.1\tLN:36148\t\t");
		out.println("@SQ\tSN:GL000247.1\tLN:36422\t\t");
		out.println("@SQ\tSN:GL000245.1\tLN:36651\t\t");
		out.println("@SQ\tSN:GL000197.1\tLN:37175\t\t");
		out.println("@SQ\tSN:GL000203.1\tLN:37498\t\t");
		out.println("@SQ\tSN:GL000246.1\tLN:38154\t\t");
		out.println("@SQ\tSN:GL000249.1\tLN:38502\t\t");
		out.println("@SQ\tSN:GL000196.1\tLN:38914\t\t");
		out.println("@SQ\tSN:GL000248.1\tLN:39786\t\t");
		out.println("@SQ\tSN:GL000244.1\tLN:39929\t\t");
		out.println("@SQ\tSN:GL000238.1\tLN:39939\t\t");
		out.println("@SQ\tSN:GL000202.1\tLN:40103\t\t");
		out.println("@SQ\tSN:GL000234.1\tLN:40531\t\t");
		out.println("@SQ\tSN:GL000232.1\tLN:40652\t\t");
		out.println("@SQ\tSN:GL000206.1\tLN:41001\t\t");
		out.println("@SQ\tSN:GL000240.1\tLN:41933\t\t");
		out.println("@SQ\tSN:GL000236.1\tLN:41934\t\t");
		out.println("@SQ\tSN:GL000241.1\tLN:42152\t\t");
		out.println("@SQ\tSN:GL000243.1\tLN:43341\t\t");
		out.println("@SQ\tSN:GL000242.1\tLN:43523\t\t");
		out.println("@SQ\tSN:GL000230.1\tLN:43691\t\t");
		out.println("@SQ\tSN:GL000237.1\tLN:45867\t\t");
		out.println("@SQ\tSN:GL000233.1\tLN:45941\t\t");
		out.println("@SQ\tSN:GL000204.1\tLN:81310\t\t");
		out.println("@SQ\tSN:GL000198.1\tLN:90085\t\t");
		out.println("@SQ\tSN:GL000208.1\tLN:92689\t\t");
		out.println("@SQ\tSN:GL000191.1\tLN:106433\t\t");
		out.println("@SQ\tSN:GL000227.1\tLN:128374\t\t");
		out.println("@SQ\tSN:GL000228.1\tLN:129120\t\t");
		out.println("@SQ\tSN:GL000214.1\tLN:137718\t\t");
		out.println("@SQ\tSN:GL000221.1\tLN:155397\t\t");
		out.println("@SQ\tSN:GL000209.1\tLN:159169\t\t");
		out.println("@SQ\tSN:GL000218.1\tLN:161147\t\t");
		out.println("@SQ\tSN:GL000220.1\tLN:161802\t\t");
		out.println("@SQ\tSN:GL000213.1\tLN:164239\t\t");
		out.println("@SQ\tSN:GL000211.1\tLN:166566\t\t");
		out.println("@SQ\tSN:GL000199.1\tLN:169874\t\t");
		out.println("@SQ\tSN:GL000217.1\tLN:172149\t\t");
		out.println("@SQ\tSN:GL000216.1\tLN:172294\t\t");
		out.println("@SQ\tSN:GL000215.1\tLN:172545\t\t");
		out.println("@SQ\tSN:GL000205.1\tLN:174588\t\t");
		out.println("@SQ\tSN:GL000219.1\tLN:179198\t\t");
		out.println("@SQ\tSN:GL000224.1\tLN:179693\t\t");
		out.println("@SQ\tSN:GL000223.1\tLN:180455\t\t");
		out.println("@SQ\tSN:GL000195.1\tLN:182896\t\t");
		out.println("@SQ\tSN:GL000212.1\tLN:186858\t\t");
		out.println("@SQ\tSN:GL000222.1\tLN:186861\t\t");
		out.println("@SQ\tSN:GL000200.1\tLN:187035\t\t");
		out.println("@SQ\tSN:GL000193.1\tLN:189789\t\t");
		out.println("@SQ\tSN:GL000194.1\tLN:191469\t\t");
		out.println("@SQ\tSN:GL000225.1\tLN:211173\t\t");
		out.println("@SQ\tSN:GL000192.1\tLN:547496\t\t");
		out.println("@RG\tID:Output_by_Java_Program\tPL:ILLUMINA\tLB:hg19");
		for (Interval i : intervals) {
			System.out.println(i.getChromosome() + "\t" + i.getStart() + "\t" + i.getEnd() + "\t+\t");
		}
	}
}
