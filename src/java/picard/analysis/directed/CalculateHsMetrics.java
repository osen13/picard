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
    static final String USAGE_SUMMARY = "Calculates hybrid Selection-specific metrics from aligned SAM or BAM files.  ";
    static final String USAGE_DETAILS = "Hybrid selection (HS) enables targeted sequencing analysis via the capture of specified genomic DNA " +
            "sequences (doi:10.1038/nbt.1523).  It is commonly used to characterize exon sequences from genomic DNA " +
            "or filter out bacterial DNA sequences from clinical samples.  <br /><br /> " +
            "" +
            "The technique involves the capture of unique regions of genomic DNA (targets) using synthetic RNA (or DNA) baits." +
            "   The baits are synthesized with biotinylated nucleotides to facilitate capture of bait:target hybrids on" +
            " streptavidin beads.  The captured target sequences are amplified, sequenced, and processed for variant calling <br /><br />." +
            "" +
            "This tool requires an aligned SAM or BAM file as well as bait and target interval files.  " +
            "The bait and target files are VCF formatted.  For information on VCF files, please see:" +
            "<br /><br />www.1000genomes.org/wiki/analysis/variant%20call%20format/vcf-variant-call-format-version-41" +
            "" +
            "<br /><br />Tool provides multiple outputs including:" +
            "<li>The numbers of aligned reads that pass the sequencing vendor's quality filter (PF) " +
            "<li>The numbers of reads and their constituent bases that map to the reference sequence " +
            "<li>The per read and per base coverage depth for both bait and target intervals<br /><br />" +
            "The tool also calculates hybrid selection penalty values; defined as the cost incurred to get 80% of target bases" +
            " to a specific amount of depth coverage.  For example, if a design consists of 10 megabases of target and" +
            " 10X coverage depth is required, then the amount of sequencing needed is the PF_ALIGNED_BASES = 10^7 * 10 * HS_PENALTY_10X.<br /><br />" +
            "" +
            "Both AT_DROPOUT and GC_DROPOUT metrics are indicated.  These metrics indicate the percentage of " +
            "reads dropped from an analysis due to the inability to map to the reference as result of excessively " +
            "AT- or GC-rich regions respectfully.<br /><br />" +

            "The PER_TARGET_COVERAGE option can be used to output GC content and mean sequence depth information for every target interval." +
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

            "Additional information on output metrics can be found at:<br /><br />" +
            "http://broadinstitute.github.io/picard/picard-metric-definitions.html#HsMetrics" +
            "<hr />"
            ;
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
