import htsjdk.variant.variantcontext.*;
import java.io.Serializable;
import java.util.List;

/**
 * Class describing a variant with an array of genotypes.
 *
 * @author Matthew Wakeling
 */
public class VariantArray implements Serializable
{
	private String chromosome;
	private int position;
	private char wild, var;
	private byte[] samples;

	@SuppressWarnings("deprecation") public VariantArray(List<String> vcfSampleNames, VariantContext context) {
		chromosome = context.getChr();
		position = context.getStart();
		wild = context.getAlleles().get(0).toString().charAt(0);
		var = context.getAlleles().get(1).toString().charAt(0);
		samples = new byte[vcfSampleNames.size()];
		for (int i = 0; i < samples.length; i++) {
			Genotype gt = context.getGenotype(vcfSampleNames.get(i));
			samples[i] = (byte) (gt.isHomVar() ? 2 : (gt.isHet() ? 1 : 0));
		}
	}

	public String getChr() {
		return chromosome;
	}

	public int getPosition() {
		return position;
	}

	public char getWild() {
		return wild;
	}

	public char getVar() {
		return var;
	}

	public byte[] getSamples() {
		return samples;
	}

	public double getRsquared(VariantArray v) {
		int[] grid = new int[9];
		for (int i = 0; i < samples.length; i++) {
			grid[((int) samples[i]) * 3 + ((int) v.samples[i])]++;
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
		if ((pa > 0.0) && (pb > 0.0) && (pa < 1.0) && (pb < 1.0)) {
			return (d / Math.sqrt(pa * (1.0 - pa) * pb * (1.0 - pb)));
		}
		return 0.0;
	}

	public String toString() {
		StringBuilder retval = new StringBuilder();
		retval.append(chromosome + ":" + position + " " + wild + ">" + var + " ");
		int[] counts = new int[3];
		for (byte b : samples) {
			counts[b]++;
			retval.append("" + b);
		}
		for (int i = 0; i < 3; i++) {
			retval.append("\t" + counts[i]);
		}
		return retval.toString();
	}
}
