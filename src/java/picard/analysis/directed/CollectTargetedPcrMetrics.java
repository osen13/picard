package picard.analysis.directed;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.IntervalList;
import picard.analysis.MetricAccumulationLevel;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.programgroups.Metrics;

import java.io.File;
import java.util.List;
import java.util.Set;

/**
 * Collect metric information for target pcr metrics runs.  See CollectTargetedMetrics and TargetPcrMetricsCollector for
 * more information
 */
@CommandLineProgramProperties(
        usage = CollectTargetedPcrMetrics.USAGE_SUMMARY + CollectTargetedPcrMetrics.USAGE_DETAIL,
        usageShort = CollectTargetedPcrMetrics.USAGE_SUMMARY,
        programGroup = Metrics.class
)
public class CollectTargetedPcrMetrics extends CollectTargetedMetrics<TargetedPcrMetrics, TargetedPcrMetricsCollector> {
    static final String USAGE_SUMMARY = "Produces targeted PCR-related metrics for a given SAM/BAM file.  ";
    static final String USAGE_DETAIL = "<br /><br />This tool calculates a set of metrics to TSCA sequencing from an aligned SAM or" +
            "BAM file.  TruSeq Custom Amplicon (TSCA) is an Illumina kit that allows researchers to sequence targeted " +
            "regions of interest in a genome.  Although there are many applications for using this targeted " +
            "approach, it is frequently used for genomic regions of high or low GC content.  For additional information, please see:" +
            "<br /><br />" +
            "www.illumina.com/content/dam/illumina-marketing/documents/products/datasheets/datasheet_truseq_custom_amplicon.pdf." +
            "<br /><br />" +
            "If a reference sequence is provided, AT/GC dropout metrics will be calculated and the PER_TARGET_COVERAGE" +
            " option can be used to output GC content and mean coverage information for each target.  " +
            "The AT/GC dropout metrics indicate the degree of inadequate coverage of a " +
            "particular region based on its AT or GC content." +
            "" +
            "<h4>Usage Example</h4>" +
            "<pre>" +
            "java -jar picard.jar CollectTargetedPcrMetrics \\<br /> " +
            "     -I=input.bam \\<br /> " +
            "     -O=PCRMetrics.bam \\<br /> " +
            "     -R=ReferenceSequence.fasta \\<br /> " +
            "     -AMPLICON_INTERVALS=Ampliconintervallist.vcf \\<br /> " +
            "     -TARGET_INTERVALS=Targetedintervallist.vcf " +
            "</pre>" +
            "" +
            "For explanations of the output metrics, see " +
            "http://broadinstitute.github.io/picard/picard-metric-definitions.html#TargetedPcrMetrics" +
            "<hr />"
            ;
    @Option(shortName = "AI", doc = "An interval list file that contains the locations of the baits used.")
    public File AMPLICON_INTERVALS;

    @Option(shortName = "N", doc = "Custom amplicon set name. If not provided it is inferred from the filename of the AMPLICON_INTERVALS intervals.", optional = true)
    public String CUSTOM_AMPLICON_SET_NAME;

    /**
     * @return AMPLICON_INTERVALS
     */
    @Override
    protected IntervalList getProbeIntervals() {
        return IntervalList.fromFile(AMPLICON_INTERVALS);
    }

    /**
     * @return CUSTOM_AMPLICON_SET_NAME
     */
    @Override
    protected String getProbeSetName() {
        return CUSTOM_AMPLICON_SET_NAME != null ? CUSTOM_AMPLICON_SET_NAME : CollectTargetedMetrics.renderProbeNameFromFile(AMPLICON_INTERVALS);
    }

    /** Stock main method. */
    public static void main(final String[] argv) {
        System.exit(new CollectTargetedPcrMetrics().instanceMain(argv));
    }

    @Override
    protected TargetedPcrMetricsCollector makeCollector(final Set<MetricAccumulationLevel> accumulationLevels,
                                                        final List<SAMReadGroupRecord> samRgRecords,
                                                        final ReferenceSequenceFile refFile,
                                                        final File perTargetCoverage,
                                                        final IntervalList targetIntervals,
                                                        final IntervalList probeIntervals,
                                                        final String probeSetName) {
        return new TargetedPcrMetricsCollector(accumulationLevels, samRgRecords, refFile, perTargetCoverage, targetIntervals, probeIntervals, probeSetName);
    }
}
