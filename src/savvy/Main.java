package savvy;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import htsjdk.samtools.util.Log;
/* 2019-08-02 Pierre Lindenbaum central point of Savvy */
public class Main extends AbstractApplication {
	protected final static Log LOG=Log.getInstance(Main.class);
	
	private static List<Class<? extends AbstractApplication>> subprogams=
			Arrays.asList(
					CoverageBinner.class,
					SavvyCNV.class,
					PrepareLinkageData.class,
					SavvyHomozygosity.class,
					SavvySharedHaplotypes.class,
					SavvyVcfHomozygosity.class
					);
	
	@Override
	protected Log getLogger() {
		return LOG;
		}
	@Override
	public int doWork(final List<String> args) throws Exception {
		
		if(args.isEmpty()) {
			getLogger().error("Program name missing: Must be one of :\n" + 
					subprogams.stream().map(C->C.getSimpleName()).
					collect(Collectors.joining("\n")));
			return -1;
			}
		final Class<? extends AbstractApplication> clazz = subprogams.stream().
				filter(C->C.getSimpleName().equalsIgnoreCase(args.get(0).replaceAll("[\\- ]", ""))).
				findFirst().
				orElse(null);
		
				
		if( clazz==null)
			{
			getLogger().error("Unknown sub-program "+args.get(0)+". Available are:"
					+ subprogams.stream().map(C->C.getSimpleName()).collect(Collectors.joining(" ")));
			return -1;
			}
		else
			{
			final AbstractApplication subProg = clazz.newInstance();
			return subProg.doWork(args.subList(1, args.size()));
			}
		}
	
	public static void main(final String[] args) {
		new Main().instanceMainWithExit(args);
	}
}
