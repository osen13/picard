package picard.analysis;

import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.programgroups.Metrics;

/**
 * Computes a number of metrics that are useful for evaluating coverage and performance of whole genome sequencing
 * experiments, same implementation as CollectWgsMetrics, with different defaults: lacks baseQ and mappingQ filters
 * and has much higher coverage cap.
 *
 * @author farjoun
 */
@CommandLineProgramProperties(
        usage = CollectRawWgsMetrics.USAGE_SUMMARY + CollectRawWgsMetrics.USAGE_DETAILS,
        usageShort = CollectRawWgsMetrics.USAGE_SUMMARY,
        programGroup = Metrics.class
)
public class CollectRawWgsMetrics extends CollectWgsMetrics{
    static final String USAGE_SUMMARY = "Writes whole genome sequencing-related metrics for a SAM or BAM file";
    static final String USAGE_DETAILS = "Useful metrics for evaluating coverage and performance of whole " +
            "genome sequencing experiments.  Although similar to the CollectWgsMetrics tool," +
            " the default thresholds for CollectRawWgsMetrics are less stringent.  For example, the CollectRawWgsMetrics" +
            " has lower base and mapping quality score thresholds as well a higher coverage cap than the CollectWgsMetrics" +
            " tool.<br /><br />" +
            "Histogram output is optional and for a given run, displays two separate outputs on the y-axis while using a single set" +
            " of values for the x-axis.  Specifically, the first column in the histogram table (x-axis) is labeled \"coverage\" and represents" +
            " different possible coverage depths.  However, it is also represents the range of values for the base quality scores and thus should probably be" +
            " labeled \"sequence depth and base quality scores\".  " +
            "The second and third columns (y-axes) correspond to the numbers of bases at a specific sequence depth" +
            " \"count\" and the numbers of bases at a particular base quality score \"baseq_count\" respectively." +
            "<br />" +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar CollectRawWgsMetrics \\<br />" +
            "     -I=MyBAM.bam \\<br />" +
            "     -O=RAWwgsmetrics.txt \\<br />" +
            "     -R=ReferenceSeq.fasta \\<br />" +
            "     -INCLUDE_BQ_HISTOGRAM=true" +
            "</pre>" +
            "<hr />" +
            "For detailed explanations of the output metrics, please see: " +
            "http://broadinstitute.github.io/picard/picard-metric-definitions.html#CollectWgsMetrics.WgsMetrics" +
            "<hr />";
    @Option(shortName="MQ", doc="Minimum mapping quality for a read to contribute coverage.")
    public int MINIMUM_MAPPING_QUALITY = 0;

    @Option(shortName="Q", doc="Minimum base quality for a base to contribute coverage.")
    public int MINIMUM_BASE_QUALITY = 3;

    @Option(shortName="CAP", doc="Treat bases with coverage exceeding this value as if they had coverage at this value.")
    public int COVERAGE_CAP = 100000;

    // rename the class so that in the metric file it is annotated differently.
    public static class RawWgsMetrics extends WgsMetrics {}

    @Override
    protected WgsMetrics generateWgsMetrics() {
        return new RawWgsMetrics();
    }

}
