package picard.analysis;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.PositionalArguments;
import picard.cmdline.programgroups.Metrics;

import java.io.File;
import java.io.FileReader;
import java.util.List;

/**
 * Compare two metrics files.
 */
@CommandLineProgramProperties(
        usage = CompareMetrics.USAGE_SUMMARY + CompareMetrics.USAGE_DETAIL,
        usageShort = CompareMetrics.USAGE_SUMMARY,
        programGroup = Metrics.class
)
public class CompareMetrics extends CommandLineProgram {
    static final String USAGE_SUMMARY =  "Compares two metrics files";
    static final String USAGE_DETAIL = "Compares the headers of the two input metrics files that have the same structure," +
            " but come from different underlying metric classes.  Outputs can be either equal" +
            " or not equal. <br /> "  +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar CompareMetrics \\<br />" +
            "     -Mymetricfile1.txt \\<br />" +
            "     -Mymetricfile2.txt" +
            "</pre>" +
            "<hr />";

    @PositionalArguments(minElements = 2, maxElements = 2)
    public List<File> metricsFiles;

    private static final Log log = Log.getInstance(CompareMetrics.class);

    @Override
    protected int doWork() {
        IOUtil.assertFilesAreReadable(metricsFiles);
        final MetricsFile<?, ?> metricsA = new MetricsFile();
        final MetricsFile<?, ?> metricsB = new MetricsFile();
        try {
            metricsA.read(new FileReader(metricsFiles.get(0)));
            metricsB.read(new FileReader(metricsFiles.get(1)));
            final boolean areEqual = metricsA.areMetricsEqual(metricsB) && metricsA.areHistogramsEqual(metricsB);
            final String status = areEqual ? "EQUAL" : "NOT EQUAL";
            log.info("Files " + metricsFiles.get(0) + " and " + metricsFiles.get(1) + "are " + status);
        } catch (final Exception e) {
            throw new PicardException(e.getMessage());
        }
        return 0;
    }

    public static void main(String[] argv) {
        new CompareMetrics().instanceMainWithExit(argv);
    }
}
