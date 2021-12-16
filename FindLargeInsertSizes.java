import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeSet;

import htsjdk.samtools.*;

/**
 * Reads a BAM file and identifies all the read pairs that have an insert size greater than 1,500.
 *
 * @author Matthew Wakeling
 */
public class FindLargeInsertSizes
{
	public static final int MAX_INSERT = 1000;
	public static final int WINDOW = 600;
	public static final int MINIMUM_MQ = 20;

	public static void main(String[] args) {
		String probandFileName = null;
		ArrayList<String> otherFileNames = new ArrayList<String>();
		String analyseChr = null;
		int analyseStart = 0;
		int analyseEnd = 0;
		boolean stddevs = false;
		int maxInsert = MAX_INSERT;
		int minMq = MINIMUM_MQ;
		for (int i = 0; i < args.length; i++) {
			if ("-limit".equals(args[i])) {
				i++;
				analyseChr = args[i];
				i++;
				analyseStart = Integer.parseInt(args[i]);
				i++;
				analyseEnd = Integer.parseInt(args[i]);
			} else if ("-stddev".equals(args[i])) {
				stddevs = true;
			} else if ("-maxInsert".equals(args[i])) {
				i++;
				maxInsert = Integer.parseInt(args[i]);
			} else if ("-minMq".equals(args[i])) {
				i++;
				minMq = Integer.parseInt(args[i]);
			} else {
				if (probandFileName == null) {
					probandFileName = args[i];
				}
				otherFileNames.add(args[i]);
			}
		}
		SamReaderFactory factory = SamReaderFactory.makeDefault();
		SamReader in = factory.open(new File(probandFileName));
		SAMRecordIterator iter = analyseChr == null ? in.iterator() : in.queryOverlapping(analyseChr, analyseStart, analyseEnd);
		TreeSet<Insert> current = new TreeSet<Insert>();
		ArrayList<SamReader> otherReaders = new ArrayList<SamReader>();
		for (String otherFileName : otherFileNames) {
			otherReaders.add(factory.open(new File(otherFileName)));
		}
		int chrNo = 1;
		// This is used to sort the reads in the same chromosome order as the BAM file.
		Map<String, Integer> chrNoMap = new HashMap<String, Integer>();
		while (iter.hasNext()) {
			SAMRecord record = iter.next();
			if (((Math.abs(record.getInferredInsertSize()) >= maxInsert) || (!record.getReferenceName().equals(record.getMateReferenceName())) || (record.getReadNegativeStrandFlag() == record.getMateNegativeStrandFlag())) && (record.getMappingQuality() >= minMq) && checkCigar(record.getReadNegativeStrandFlag(), record.getCigarString()) && checkCigar(record.getMateNegativeStrandFlag(), record.getStringAttribute("MC"))) {
				if (!chrNoMap.containsKey(record.getReferenceName())) {
					chrNoMap.put(record.getReferenceName(), chrNo++);
				}
				//System.out.println(record.getReferenceName() + ":" + record.getAlignmentStart() + (record.getReadNegativeStrandFlag() ? " R\t" : " F\t") + record.getMateReferenceName() + ":" + record.getMateAlignmentStart() + (record.getMateNegativeStrandFlag() ? " R" : " F"));
				Insert newInsert = new Insert(record.getReferenceName(), record.getMateReferenceName(), record.getAlignmentStart(), record.getAlignmentStart() + cigarLength(record.getCigarString()) - 1, record.getReadNegativeStrandFlag(), record.getMateAlignmentStart(), record.getMateAlignmentStart() + cigarLength(record.getStringAttribute("MC")) - 1, record.getMateNegativeStrandFlag(), chrNoMap);
				current.add(newInsert);
				Insert lastInsert = current.first();
				while ((lastInsert.getStart1() < newInsert.getStart1() - WINDOW) || (!lastInsert.getChr().equals(newInsert.getChr()))) {
					current.remove(lastInsert);
					if (!lastInsert.isProcessed()) {
						//System.out.println(lastInsert.getChr() + ":" + lastInsert.getStart1() + (lastInsert.getNeg1() ? " R\t" : " F\t") + lastInsert.getChr() + ":" + lastInsert.getStart2() + (lastInsert.getNeg2() ? " R" : " F") + " inspected");
						long start1Total = lastInsert.getStart1();
						double start1SqTotal = (1.0 * lastInsert.getStart1()) * (1.0 * lastInsert.getStart1());
						long start2Total = lastInsert.getStart2();
						double start2SqTotal = (1.0 * lastInsert.getStart2()) * (1.0 * lastInsert.getStart2());
						int startCount = 1;
						int start1Min = Integer.MAX_VALUE;
						int start1Max = Integer.MIN_VALUE;
						int start2Min = Integer.MAX_VALUE;
						int start2Max = Integer.MIN_VALUE;
						for (Insert insert : current) {
							if ((!insert.isProcessed()) && insert.getChr().equals(lastInsert.getChr()) && insert.getOtherChr().equals(lastInsert.getOtherChr()) && (insert.getNeg1() == lastInsert.getNeg1()) && (insert.getNeg2() == lastInsert.getNeg2()) && (Math.abs(insert.getStart1() - lastInsert.getStart1()) < WINDOW) && (Math.abs(insert.getStart2() - lastInsert.getStart2()) < WINDOW)) {
								//System.out.println(insert.getChr() + ":" + insert.getStart1() + (insert.getNeg1() ? " R\t" : " F\t") + insert.getChr() + ":" + insert.getStart2() + (insert.getNeg2() ? " R" : " F") + " processed");
								start1Total += insert.getStart1();
								start1SqTotal += (1.0 * insert.getStart1()) * (1.0 * insert.getStart1());
								start2Total += insert.getStart2();
								start2SqTotal += (1.0 * insert.getStart2()) * (1.0 * insert.getStart2());
								start1Min = Math.min(start1Min, insert.getStart1());
								start1Max = Math.max(start1Max, insert.getEnd1());
								start2Min = Math.min(start2Min, insert.getStart2());
								start2Max = Math.max(start2Max, insert.getEnd2());
								startCount++;
								insert.setProcessed();
							}
						}
						// Require at least three reads.
						// Compare to number of weird reads in this locality, to avoid processing areas where there are read pairs pointing everywhere.
						if ((startCount > 3) && (startCount * 5 + 1 > current.size())) {
							String desc = "";
							boolean neg1 = lastInsert.getNeg1();
							boolean neg2 = lastInsert.getNeg2();
							if (!lastInsert.getChr().equals(lastInsert.getOtherChr())) {
								desc = "Translocation";
							} else if (!neg1) {
								if (neg2) {
									desc = start1Total < start2Total ? "Deletion_left_edge" : "Duplication_right_edge";
								} else {
									desc = start1Total < start2Total ? "Inversion_outer_left" : "Inversion_inner_right";
								}
							} else {
								if (neg2) {
									desc = start1Total < start2Total ? "Inversion_inner_left" : "Inversion_outer_right";
								} else {
									desc = start1Total > start2Total ? "Deletion_right_edge" : "Duplication_left_edge";
								}
							}
							int pos1 = neg1 ? start1Min : start1Max;
							int pos2 = neg2 ? start2Min : start2Max;
							double stddev1 = Math.sqrt((start1SqTotal - ((1.0 * start1Total) * (1.0 * start1Total) / startCount)) / startCount);
							double stddev2 = Math.sqrt((start2SqTotal - ((1.0 * start2Total) * (1.0 * start2Total) / startCount)) / startCount);
							// Check other samples.
							StringBuilder otherText = new StringBuilder();
							for (int i = 0; i < otherReaders.size(); i++) {
								long otherStart1Total = 0L;
								double otherStart1SqTotal = 0.0;
								long otherStart2Total = 0L;
								double otherStart2SqTotal = 0.0;
								int otherStartCount = 0;
								int otherTotalCount = 0;
								SAMRecordIterator otherIter = otherReaders.get(i).queryOverlapping(lastInsert.getChr(), pos1 - WINDOW, pos1 + WINDOW);
								while (otherIter.hasNext()) {
									SAMRecord otherRecord = otherIter.next();
									if (lastInsert.getChr().equals(otherRecord.getReferenceName()) && (otherRecord.getMappingQuality() >= minMq) && (otherRecord.getReadNegativeStrandFlag() == lastInsert.getNeg1()) && (Math.abs(otherRecord.getAlignmentStart() - pos1) < WINDOW) && checkCigar(otherRecord.getReadNegativeStrandFlag(), otherRecord.getCigarString()) && checkCigar(otherRecord.getMateNegativeStrandFlag(), otherRecord.getStringAttribute("MC"))) {
										otherTotalCount++;
										if (lastInsert.getChr().equals(otherRecord.getReferenceName()) && lastInsert.getOtherChr().equals(otherRecord.getMateReferenceName()) && (otherRecord.getMappingQuality() >= minMq) && (otherRecord.getReadNegativeStrandFlag() == lastInsert.getNeg1()) && (otherRecord.getMateNegativeStrandFlag() == lastInsert.getNeg2()) && (Math.abs(otherRecord.getAlignmentStart() - pos1) < WINDOW) && (Math.abs(otherRecord.getMateAlignmentStart() - pos2) < WINDOW)) {
											int start1 = otherRecord.getAlignmentStart();
											int start2 = otherRecord.getMateAlignmentStart();
											otherStart1Total += start1;
											otherStart1SqTotal += (1.0 * start1) * (1.0 * start1);
											otherStart2Total += start2;
											otherStart2SqTotal += (1.0 * start2) * (1.0 * start2);
											otherStartCount++;
										}
									}
								}
								otherIter.close();
								double otherStddev1 = otherStartCount == 0 ? 0.0 : Math.sqrt((otherStart1SqTotal - ((1.0 * otherStart1Total) * (1.0 * otherStart1Total) / otherStartCount)) / otherStartCount);
								double otherStddev2 = otherStartCount == 0 ? 0.0 : Math.sqrt((otherStart2SqTotal - ((1.0 * otherStart2Total) * (1.0 * otherStart2Total) / otherStartCount)) / otherStartCount);
								if (stddevs) {
									otherText.append("\t" + (otherStartCount == 0 ? "?" : lastInsert.getChr() + ":" + (otherStart1Total / otherStartCount)) + "\t+-" + otherStddev1 + "\t" + (otherStartCount == 0 ? "?" : lastInsert.getOtherChr() + ":" + (otherStart2Total / otherStartCount)) + "\t+-" + otherStddev2 + "\t" + otherStartCount + "\t" + otherTotalCount);
								} else {
									otherText.append("\t" + otherStartCount + "\t" + otherTotalCount);
								}
							}
							if (stddevs) {
								System.out.println(lastInsert.getChr() + ":" + pos1 + "\t+-" + stddev1 + (lastInsert.getNeg1() ? "\tR\t" : "\tF\t") + lastInsert.getOtherChr() + ":" + pos2 + "\t+-" + stddev2 + (lastInsert.getNeg2() ? "\tR\t" : "\tF\t") + startCount + "\t" + (lastInsert.getChr().equals(lastInsert.getOtherChr()) ? "" + (pos2 - pos1) : "?") + "\t" + desc + otherText.toString());
							} else {
								System.out.println(lastInsert.getChr() + ":" + pos1 + (lastInsert.getNeg1() ? "\tR\t" : "\tF\t") + lastInsert.getOtherChr() + ":" + pos2 + (lastInsert.getNeg2() ? "\tR\t" : "\tF\t") + startCount + "\t" + (lastInsert.getChr().equals(lastInsert.getOtherChr()) ? "" + (pos2 - pos1) : "?") + "\t" + desc + otherText.toString());
							}
						}
					}
					lastInsert = current.first();
				}
			}
		}
	}

	public static boolean checkCigar(boolean reverse, String cigar) {
		if (cigar == null) {
			return false;
		}
		if (reverse) {
			// The read is pointing leftwards, so the left-hand side can be clipped, but we don't want the right-hand side to be clipped.
			// If the right-hand side were clipped, that would indicate that we are reading an inclusion that has been mis-mapped here.
			return !((cigar.charAt(cigar.length() - 1) == 'H') || (cigar.charAt(cigar.length() - 1) == 'S'));
		} else {
			// We don't want the left-hand side to be clipped.
			int charPos = 0;
			while ((cigar.charAt(charPos) >= '0') && (cigar.charAt(charPos) <= '9')) {
				charPos++;
			}
			return !((cigar.charAt(charPos) == 'H') || (cigar.charAt(charPos) <= '9'));
		}
	}

	public static int cigarLength(String cigar) {
		int retval = 0;
		for (int i = 0; i < cigar.length(); i++) {
			int number = 0;
			while ((cigar.charAt(i) >= '0') && (cigar.charAt(i) <= '9')) {
				number = number * 10 + (cigar.charAt(i) - '0');
				i++;
			}
			switch (cigar.charAt(i)) {
				case 'M':
				case 'D':
					retval += number;
					break;
				case 'S':
				case 'H':
				case 'I':
					break;
				default:
					throw new RuntimeException("Unrecognised cigar symbol " + cigar.charAt(i) + " in " + cigar);
			}
		}
		return retval;
	}

	public static class Insert implements Comparable<Insert>
	{
		public static int newSeq = 0;

		String chr, otherChr;
		private int start1, start2, end1, end2, seq;
		private boolean neg1, neg2, processed;
		private Map<String, Integer> chrNoMap;

		public Insert(String chr, String otherChr, int start1, int end1, boolean neg1, int start2, int end2, boolean neg2, Map<String, Integer> chrNoMap) {
			this.chr = chr;
			this.otherChr = otherChr;
			this.start1 = start1;
			this.end1 = end1;
			this.neg1 = neg1;
			this.start2 = start2;
			this.end2 = end2;
			this.neg2 = neg2;
			this.seq = newSeq++;
			this.processed = false;
			this.chrNoMap = chrNoMap;
		}

		public String getChr() {
			return chr;
		}

		public String getOtherChr() {
			return otherChr;
		}

		public int getStart1() {
			return start1;
		}

		public int getEnd1() {
			return end1;
		}

		public boolean getNeg1() {
			return neg1;
		}

		public int getStart2() {
			return start2;
		}

		public int getEnd2() {
			return end2;
		}

		public boolean getNeg2() {
			return neg2;
		}

		public boolean isProcessed() {
			return processed;
		}

		public void setProcessed() {
			processed = true;
		}

		public int compareTo(Insert i) {
			int chrNo = chrNoMap.get(chr);
			int iChrNo = chrNoMap.get(i.chr);
			if (chrNo < iChrNo) {
				return -1;
			} else if (chrNo > iChrNo) {
				return 1;
			} else if (start1 < i.start1) {
				return -1;
			} else if (start1 > i.start1) {
				return 1;
			} else if (seq < i.seq) {
				return -1;
			} else if (seq > i.seq) {
				return 1;
			} else {
				return 0;
			}
		}
	}
}
