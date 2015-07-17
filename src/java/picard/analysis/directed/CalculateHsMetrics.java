/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package picard.analysis.directed;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.StringUtil;
import picard.analysis.MetricAccumulationLevel;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.programgroups.Metrics;

import java.io.File;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

/**
 * Calculates a set of HS metrics from a sam or bam file.  See HsMetricsCollector and CollectTargetedMetrics for more details.
 *
 * @author Tim Fennell
 */

@CommandLineProgramProperties(
        usage = CalculateHsMetrics.USAGE_SUMMARY + CalculateHsMetrics.USAGE_DETAILS,
        usageShort = CalculateHsMetrics.USAGE_SUMMARY,
        programGroup = Metrics.class
)
public class CalculateHsMetrics extends CollectTargetedMetrics<HsMetrics, HsMetricCollector> {
    static final String USAGE_SUMMARY = "Calculates a set of Hybrid Selection-specific metrics from an aligned SAM or BAM file "  ;
    static final String USAGE_DETAILS = "Calculates hybrid Selection-specific metrics from aligned SAM or BAM files. <br /><br />  " +
            "Hybrid selection (HS) enables targeted sequencing analysis via the capture of specified genomic DNA " +
            "sequences (doi:10.1038/nbt.1523).  It is commonly used to characterize exon sequences from genomic DNA " +
            "or filter out bacterial DNA sequences from clinical samples.<br /><br /> " +
            "The technique captures unique regions of genomic DNA (targets) from a biological samples via the hybridization of synthetically" +
            " produced RNAs (baits) via solution hybridization.  The baits are synthesized with biotinylated nucleotides" +
            " to facilitate capture of bait:target hybrids on streptavidin beads." +
            "Targets captured by these baits are amplified, sequenced, and processed for variant calling." +
            "" +
            "Sequences are aligned to the reference sequence and the alignment " +
            "efficacies are characterized and output.  Bait intervals are the sequences used for baits and " +
            "target intervals are the sequences of interest in the sample genomic DNA. "  +
            "" +
            "The tool produces both AT_DROPOUT and GC_DROPOUT metrics, which indicate the percentage of " +
            "reads dropped from an analysis due to the inability to map to the reference as result of excessively " +
            "AT-rich or GC-rich regions respectfully. <br /><br />" +
            "The PER_TARGET_COVERAGE option can be used to output GC content and mean sequence depth information for every target." +
            "<br /><br />" +
            "<h4>Usage Example:</h4>"+
            "<pre>" +
            "java -jar picard.jar CalculateHsMetrics \\<br />" +
            "     -I=Input.bam \\<br />" +
            "     -O=HSMetrics.txt \\<br />" +
            "     -R=ReferenceSequence.fasta \\<br />" +
            "     -BAIT_INTERVALS=Baitintervallist.vcf \\<br />" +
            "     -TARGET_INTERVALS=Targetintervallist.vcf" +
            "</pre> "   +
            "<hr />";
    @Option(shortName = "BI", doc = "An interval list file that contains the locations of the baits used.", minElements=1)
    public List<File> BAIT_INTERVALS;

    @Option(shortName = "N", doc = "Bait set name. If not provided it is inferred from the filename of the bait intervals.", optional = true)
    public String BAIT_SET_NAME;

    @Override
    protected IntervalList getProbeIntervals() {
        for (final File file : BAIT_INTERVALS) IOUtil.assertFileIsReadable(file);
        return IntervalList.fromFiles(BAIT_INTERVALS);
    }

    @Override
    protected String getProbeSetName() {
        if (BAIT_SET_NAME != null) {
            return BAIT_SET_NAME;
        } else {
            final SortedSet<String> baitSetNames = new TreeSet<String>();
            for (final File file : BAIT_INTERVALS) {
                baitSetNames.add(CollectTargetedMetrics.renderProbeNameFromFile(file));
            }
            return StringUtil.join(".", baitSetNames);
        }
    }

    /** Stock main method. */
    public static void main(final String[] argv) {
        System.exit(new CalculateHsMetrics().instanceMain(argv));
    }

    @Override
    protected HsMetricCollector makeCollector(final Set<MetricAccumulationLevel> accumulationLevels,
                                              final List<SAMReadGroupRecord> samRgRecords,
                                              final ReferenceSequenceFile refFile,
                                              final File perTargetCoverage,
                                              final IntervalList targetIntervals,
                                              final IntervalList probeIntervals,
                                              final String probeSetName) {
        return new HsMetricCollector(accumulationLevels, samRgRecords, refFile, perTargetCoverage, targetIntervals, probeIntervals, probeSetName);
    }
}
