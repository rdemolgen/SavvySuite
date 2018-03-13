import java.util.Stack;

/**
 * Methods to calculate the Fisher's Exact Test probability values from a 2x2 contingency table.
 *
 * @author Matthew Wakeling
 */
public class FishersExactTest
{
	public static void main(String[] args) {
		//System.out.println(fisher(1, 3, 3, 2));
		/*int s = 100;
		int a = 34;
		for (int homVar = 0; homVar <= a / 2; homVar++) {
			int het = a - homVar * 2;
			int homRef = s - homVar - het;
			System.out.println("hardyP(" + homRef + ", " + het + ", " + homVar + ")\t" + hardyP(homRef, het, homVar) + "\t" + hardyS(homRef, het, homVar));
		}*/
		int s = 252;
		for (int het = 0; het <= s; het++) {
			for (int homRef = 0; homRef <= s; homRef++) {
				if (het + homRef <= s) {
					System.out.println(het + "\t" + homRef + "\t" + hardyS(homRef, het, s - het - homRef));
				} else {
					System.out.println(het + "\t" + homRef + "\t0");
				}
			}
		}
		//System.out.println(hardyS(0, 0, 252));
	}

	/**
	 * Calculate the Fisher's Exact Test probabilities, as an array of double. The first element is the left 1-tail. The second is the right 1-tail. The third is the 2-tail probability.
	 */
	public static double[] fisher(int a, int b, int c, int d) {
		int ab = a + b;
		int ac = a + c;
		int bd = b + d;
		int cd = c + d;
		int minAD = Math.min(a, d);
		int minBC = Math.min(b, c);
		double pExact = hypergeometric(a, b, c, d);
		if (minAD < minBC) {
			double pLower = 0.0;
			for (int i = 1; i <= minAD; i++) {
				pLower += hypergeometric(a - i, b + i, c + i, d - i);
			}
			double pHigher = 0.0;
			int i = minBC;
			boolean probLower = true;
			while (probLower && (i > 0)) {
				double p = hypergeometric(a + i, b - i, c - i, d + i);
				if (p < pExact) {
					pHigher += p;
				} else {
					probLower = false;
				}
				i--;
			}
			//System.out.println("Left: " + (pLower + pExact) + ", Right: " + (1.0 - pLower) + ", 2-tail: " + (pLower + pExact + pHigher));
			return new double[] {(pLower + pExact), (1.0 - pLower), (pLower + pExact + pHigher)};
		} else {
			double pHigher = 0.0;
			for (int i = 1; i < minBC; i++) {
				pHigher += hypergeometric(a + i, b - i, c - i, d + i);
			}
			double pLower = 0.0;
			int i = minAD;
			boolean probLower = true;
			while (probLower && (i > 0)) {
				double p = hypergeometric(a - i, b + i, c + i, d - i);
				if (p < pExact) {
					pLower += p;
				} else {
					probLower = false;
				}
				i--;
			}
			//System.out.println("Left: " + (1.0 - pHigher) + ", Right: " + (pHigher + pExact) + ", 2-tail: " + (pLower + pExact + pHigher));
			return new double[] {(1.0 - pHigher), (pHigher + pExact), (pLower + pExact + pHigher)};
		}
	}

	/**
	 * Calculate the probability of the hypergeometric distribution.
	 */
	public static double hypergeometric(int a, int b, int c, int d) {
		int ab = a + b;
		int ac = a + c;
		int bd = b + d;
		int cd = c + d;
		int abcd = a + b + c + d;
		// Result is fact(ab) * fact(ac) * fact(bd) * fact(cd) / fact(a) / fact(b) / fact(c) / fact(d) / fact(abcd)
		// Or alternatively, (fact(ab) / fact(a)) * (fact(bd) / fact(b)) * (fact(ac) / fact(c)) * (fact(cd) / fact(d)) / fact(abcd)
		Stack<Integer> mult = new Stack<Integer>();
		Stack<Integer> div = new Stack<Integer>();
		for (int i = a + 1; i <= ab; i++) {
			mult.push(i);
		}
		for (int i = b + 1; i <= bd; i++) {
			mult.push(i);
		}
		for (int i = c + 1; i <= ac; i++) {
			mult.push(i);
		}
		for (int i = d + 1; i <= cd; i++) {
			mult.push(i);
		}
		for (int i = 2; i <= abcd; i++) {
			div.push(i);
		}
		//System.out.println(mult + " / " + div);
		double retval = 1.0;
		while (!(mult.empty() && div.empty())) {
			if (((retval > 1.0) && (!div.empty())) || mult.empty()) {
				retval = retval / div.pop();
			} else {
				retval = retval * mult.pop();
			}
		}
		//System.err.println("hyp(" + a + ", " + b + ", " + c + ", " + d + ") = " + retval);
		return retval;
	}

	/**
	 * Calculate the two-tailed Hardy-Weinberg significance, given the number of homRef, het, and homVar samples.
	 */
	public static double hardyS(int homRef, int het, int homVar) {
		int s = homRef + het + homVar;
		int a = het + 2 * homVar;
		double prob = hardyP(homRef, het, homVar);
		double totalProb = prob;
		int homVarCount = Math.max(0, a - s);
		double posProb = 0.0;
		while ((homVarCount < homVar) && (posProb < prob)) {
			int hetCount = a - homVarCount * 2;
			int homRefCount = s - homVarCount - hetCount;
			posProb = hardyP(homRefCount, hetCount, homVarCount);
			if (posProb < prob) {
				totalProb += posProb;
			}
			homVarCount++;
		}
		homVarCount = a / 2;
		posProb = 0.0;
		while ((homVarCount > homVar) && (posProb < prob)) {
			int hetCount = a - homVarCount * 2;
			int homRefCount = s - homVarCount - hetCount;
			posProb = hardyP(homRefCount, hetCount, homVarCount);
			if (posProb < prob) {
				totalProb += posProb;
			}
			homVarCount--;
		}
		return totalProb;
	}

	private static double[][][] hardyPCache = null; // Dimensions are samples, homVar, het

	/**
	 * Calculate the Hardy-Weinberg probability, given the number of homRef, het, and homVar samples, given the total allele freq.
	 */
	public static double hardyP(int homRef, int het, int homVar) {
		// Cache result, because it can take a while to calculate.
		int s = homRef + het + homVar;
		if (hardyPCache == null) {
			hardyPCache = new double[s + 1][][];
		}
		if (hardyPCache.length < s + 1) {
			double[][][] newHardyPCache = new double[s + 1][][];
			for (int i = 0; i < hardyPCache.length; i++) {
				newHardyPCache[i] = hardyPCache[i];
			}
			hardyPCache = newHardyPCache;
		}
		if (hardyPCache[s] == null) {
			hardyPCache[s] = new double[s + 1][s + 1];
			for (int i = 0; i <= s; i++) {
				for (int o = 0; o <= s; o++) {
					hardyPCache[s][i][o] = -1.0;
				}
			}
		}
		if (hardyPCache[s][homVar][het] > -1.0) {
			return hardyPCache[s][homVar][het];
		}
		int ab = homRef * 2 + het;
		int bc = het + 2 * homVar;
		int n = homRef + het + homVar;
		Stack<Integer> mult = new Stack<Integer>();
		Stack<Integer> div = new Stack<Integer>();
		for (int i = 2; i <= n; i++) {
			mult.add(i);
		}
		for (int i = 2; i <= homRef; i++) {
			div.add(i);
		}
		for (int i = 2; i <= het; i++) {
			div.add(i);
		}
		for (int i = 2; i <= homVar; i++) {
			div.add(i);
		}
		choose(div, mult, 2 * n, bc);
		//System.err.println("mult: " + mult + ", div: " + div);
		for (int i = 0; i < het; i++) {
			mult.add(2);
		}
		double retval = 1.0;
		while (!(mult.empty() && div.empty())) {
			if (((retval > 1.0) && (!div.empty())) || mult.empty()) {
				retval = retval / div.pop();
			} else {
				retval = retval * mult.pop();
			}
		}
		//System.err.println("hardyP(" + homRef + ", " + het + ", " + homVar + ") = " + retval);
		hardyPCache[s][homVar][het] = retval;
		return retval;
	}

	/**
	 * Add elements to the multiply and divide stack for n choose k.
	 */
	public static void choose(Stack<Integer> mult, Stack<Integer> div, int n, int k) {
		//System.err.println("choose(" + n + ", " + k + ")");
		for (int i = k + 1; i<= n; i++) {
			//System.err.println("Adding " + i + " to mult");
			mult.add(i);
		}
		for (int i = 2; i <= n - k; i++) {
			//System.err.println("Adding " + i + " to div");
			div.add(i);
		}
	}
}
