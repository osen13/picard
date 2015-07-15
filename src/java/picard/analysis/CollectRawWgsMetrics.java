package picard.analysis;

import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.programgroups.Metrics;

/**
 * Computes a number of metrics that are useful for evaluating coverage and performance of whole genome sequencing experiments, same implementation as CollectWgsMetrics, with different defaults: lacks baseQ and mappingQ filters and has much higher coverage cap.
 *
 * @author farjoun
 */
@CommandLineProgramProperties(
        usage = "Computes a number of metrics that are useful for evaluating coverage and performance of whole " +
                "genome sequencing experiments using raw sequence data.  Although similar to the CollectWgsMetrics tool," +
                " the defaults thresholds are different e.g. the CollectRawWgsMetrics has lower base and mapping quality" +
                " score thresholds as well a higher coverage cap than the CollectWgsMetrics tool.<br /><br />" +
                "Histogram output is optional and displays the numbers of reads for each particular depth as well as the " +
                "mean value for the summed base-quality scores for the reads." +
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
        "<hr />",
        usageShort = "Writes whole genome sequencing-related metrics for a SAM or BAM file",
        programGroup = Metrics.class
)
public class CollectRawWgsMetrics extends CollectWgsMetrics{

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
