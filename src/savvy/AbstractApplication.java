package savvy;

import java.util.Arrays;
import java.util.List;

import htsjdk.samtools.util.Log;


public abstract class AbstractApplication {

public abstract int doWork(final List<String> args) throws Exception;

protected abstract Log getLogger();

protected String oneAndOnlyOneFile(final List<String> args)	 {
	if(args.size()!=1) {
		throw new IllegalArgumentException("Expected one and only one argument but got "+args.size()+" "+String.join(" ",args));
		}
	return args.get(0);
	}

public int instanceMain(final String args[]) {
	try {
		return doWork(Arrays.asList(args));
		}
	catch(final Throwable err) {
		getLogger().error(err);
		return -1;
		}
	}
	
public void instanceMainWithExit(final String args[]) {
	int ret = instanceMain(args);
	System.exit(ret);
	}
}
