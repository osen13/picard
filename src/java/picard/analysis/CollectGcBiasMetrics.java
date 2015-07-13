/*
 * The MIT License
 *
 * Copyright (c) 2015 The Broad Institute
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

package picard.analysis;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.StringUtil;
import picard.PicardException;
import picard.analysis.directed.GcBiasMetricsCollector;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.programgroups.Metrics;
import picard.metrics.GcBiasMetrics;
import picard.util.RExecutor;

import java.io.File;
import java.text.NumberFormat;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Tool to collect information about GC bias in the reads in a given BAM file. Computes
 * the number of windows (of size specified by WINDOW_SIZE) in the genome at each GC%
 * and counts the number of read starts in each GC bin.  What is output and plotted is
 * the "normalized coverage" in each bin - i.e. the number of reads per window normalized
 * to the average number of reads per window across the whole genome.
 *
 * @author Tim Fennell
 * edited by Kylee Bergin
 */
@CommandLineProgramProperties(
        usage = CollectGcBiasMetrics.USAGE_SUMMARY + CollectGcBiasMetrics.USAGE_DETAILS,
        usageShort = CollectGcBiasMetrics.USAGE_SUMMARY,
        programGroup = Metrics.class
)
public class CollectGcBiasMetrics extends SinglePassSamProgram {
    static final String USAGE_SUMMARY = "Collects information regarding GC bias, from a SAM/BAM input file";
    static final String USAGE_DETAILS = "Tool that collects information about G + C content, enabling the user to" +
            " determine how GC-content affects the sequencing and alignment of their data.<br /> <br />" +
            "" +
            "Briefly, G + C content of DNA can affect base-calling accuracy due in part, to the hydrogen bonding" +
            " differences between G/C and A/T pairs.<br /><br />" +
            "" +
            "For each run, the corresponding reference sequence is divided into bins or windows based on the percentage of G + C" +
            " content ranging from 0 - 100%.  " +
            "This produces a distribution of which, each bin is a defined length of sequence, the default " +
            "value is set at 100 bases.  The distribution plots the frequency of bins with a specific GC content.  " +
            "The mean of the distribution will vary among organisms.  However, human DNA has a mean GC-content of 40%, suggesting " +
            "preponderance of AT-rich regions.  <br /><br />" +

            "GC-bias is calculated and output as both summary (optional) and detailed metrics (required).  " +
            "The GcBiasSummaryMetrics provides high-level metrics that capture run-specific bias information including" +
            " WINDOW_SIZE, ALIGNED_READS, TOTAL_CLUSTERS, AT_DROPOUT, and GC_DROPOUT.  While WINDOW_SIZE refers to the" +
            "numbers of bases used for the distribution (above), the ALIGNED_READS and" +
            " TOTAL_CLUSTERS are the total number of aligned reads and the total number of reads (after filtering) " +
            "produced in a run.  In addition, the tool produces both AT_DROPOUT and GC_DROPOUT metrics, which indicate the percentage of" +
            " reads dropped from an analysis due to the inability to map to the reference as result of excessively" +
            " AT-rich or GC-rich regions respectfully. <br /><br />" +
            "" +
            "GcBiasDetailedMetrics produces both a chart (pdf) and a table of data.  These data include GC percentages " +
            "for each bin (GC), the numbers of windows corresponding to each bin (WINDOWS), the numbers of reads that start within a bin (READ_STARTS), " +
            "and the mean base quality of the reads that correspond to a specific GC-content distribution window (MEAN_BASE_QUALITY)." +
            "  In addition, NORMALIZED_COVERAGE is a relative measure of sequence coverage by the reads at a particular GC-content." +
            "  The percentage of \"coverage\" or depth in a GC bin is calculated by dividing the number of reads of a particular GC content, " +
            "by the mean number of reads of all GC bins.  A number of 1 represents mean coverage, a number less than " +
            "one represents lower than mean coverage (e.g. 0.5 means half as much coverage as average) while a " +
            "number greater than one represents higher than mean coverage (e.g. 3.1 means this GC bin has 3.1 times" +
            " more reads per window than average).  Tool also plots mean base-quality scores of the reads within each" +
            " GC-content bin, enabling the user to determine how base-quality scores vary with GC-content.<br />"+

            "<h4>Usage Example:</h4>"+
            "<pre>" +
            "java -jar picard.jar CollectGcBiasMetrics \\<br />"+
            "     -I=Input.bam \\<br />"+
            "     -O=GCBiasMetrics.txt \\<br />"+
            "     -CHART=GCBiasMetrics.pdf \\<br />"+
            "     -R=ReferenceSequence.fasta"+
            "</pre>"+
            "For detailed explanations of the output metrics, please see:" +
            "https://broadinstitute.github.io/picard/picard-metric-definitions.html#GcBiasMetrics" +
            "<hr />"
    ;
    /** The location of the R script to do the plotting. */
    private static final File R_SCRIPT = new File("src/scripts/picard/analysis/","gcBias.R");

    // Usage and parameters

    @Option(shortName = "CHART", doc = "The PDF file to render the chart to.")
    public File CHART_OUTPUT;

    @Option(shortName = "S", doc = "The text file to write summary metrics to.", optional = true)
    public File SUMMARY_OUTPUT;

    @Option(doc = "The size of windows on the genome that are used to bin reads.")
    public int WINDOW_SIZE = 100;

    @Option(doc = "For summary metrics, exclude GC windows that include less than this fraction of the genome.")
    public double MINIMUM_GENOME_FRACTION = 0.00001;

    @Option(shortName = "BS", doc = "Whether the SAM or BAM file consists of bisulfite sequenced reads.")
    public boolean IS_BISULFITE_SEQUENCED = false;

    @Option(shortName = "LEVEL", doc = "The level(s) at which to accumulate metrics.")
    public Set<MetricAccumulationLevel> METRIC_ACCUMULATION_LEVEL = CollectionUtil.makeSet(MetricAccumulationLevel.ALL_READS);

    // Calculates GcBiasMetrics for all METRIC_ACCUMULATION_LEVELs provided
    private GcBiasMetricsCollector multiCollector;

    //windowSize is the size of the scanning window that goes over the reference
    private final int windowSize = WINDOW_SIZE;
    final int[] windowsByGc = new int[WINDOWS];

    // Histograms to track the number of windows at each GC, and the number of read starts
    // at windows of each GC. Need 101 to get from 0-100.
    private static final int WINDOWS = 101;

    //Hash map of gc[] with reference name as key
    private final Map<String, byte[]> gcByRef = new HashMap<String, byte[]>();

    ////////////////////////////////////////////////////////////////////////////
    // Stock main method
    ////////////////////////////////////////////////////////////////////////////
    public static void main(final String[] args) {
        System.exit(new CollectGcBiasMetrics().instanceMain(args));
    }

    /////////////////////////////////////////////////////////////////////////////
    // Setup calculates gc[] for the reference. Must be done at startup to avoid
    // missing reference sequences in the case of small files that may
    // not have reads aligning to every reference sequence
    /////////////////////////////////////////////////////////////////////////////
    @Override
    protected void setup(final SAMFileHeader header, final File samFile) {
        IOUtil.assertFileIsWritable(CHART_OUTPUT);

        if (SUMMARY_OUTPUT != null) IOUtil.assertFileIsWritable(SUMMARY_OUTPUT);

        IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);

        final ReferenceSequenceFile refFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(REFERENCE_SEQUENCE);
        ReferenceSequence ref;

        while ((ref = refFile.nextSequence()) != null) {
            final byte[] refBases = ref.getBases();
            final String refName = ref.getName();
            StringUtil.toUpperCase(refBases);
            final int refLength = refBases.length;
            final int lastWindowStart = refLength - windowSize;
            final byte[] gc = calculateAllGcs(refBases, windowsByGc, lastWindowStart);
            gcByRef.put(refName, gc);
        }
        //Delegate actual collection to GcBiasMetricCollector
        multiCollector = new GcBiasMetricsCollector(METRIC_ACCUMULATION_LEVEL, gcByRef, windowsByGc, header.getReadGroups(), windowSize, IS_BISULFITE_SEQUENCED);
    }

    ////////////////////////////////////////////////////////////////////////////
    // MultiCollector acceptRead
    ////////////////////////////////////////////////////////////////////////////
    @Override
    protected void acceptRead(final SAMRecord rec, final ReferenceSequence ref) {
        multiCollector.acceptRecord(rec, ref);
    }

    /////////////////////////////////////////////////////////////////////////////
    // Write out all levels of normalized coverage metrics to a file
    /////////////////////////////////////////////////////////////////////////////
    @Override
    protected void finish() {
        multiCollector.finish();
        final MetricsFile<GcBiasMetrics, Integer> file = getMetricsFile();
        final MetricsFile<GcBiasDetailMetrics, ?> detailMetricsFile = getMetricsFile();
        final MetricsFile<GcBiasSummaryMetrics, ?> summaryMetricsFile = getMetricsFile();
        multiCollector.addAllLevelsToFile(file);
        final List<GcBiasMetrics> gcBiasMetricsList = file.getMetrics();
        for(final GcBiasMetrics gcbm : gcBiasMetricsList){
            final List<GcBiasDetailMetrics> gcDetailList = gcbm.DETAILS.getMetrics();
            for(final GcBiasDetailMetrics d : gcDetailList) {
                detailMetricsFile.addMetric(d);
            }
            summaryMetricsFile.addMetric(gcbm.SUMMARY);
        }
        detailMetricsFile.write(OUTPUT);
        summaryMetricsFile.write(SUMMARY_OUTPUT);

        final NumberFormat fmt = NumberFormat.getIntegerInstance();
        fmt.setGroupingUsed(true);
        RExecutor.executeFromFile(R_SCRIPT,
                OUTPUT.getAbsolutePath(),
                SUMMARY_OUTPUT.getAbsolutePath(),
                CHART_OUTPUT.getAbsolutePath(),
                String.valueOf(WINDOW_SIZE));
    }

    /////////////////////////////////////////////////////////////////////////////
    // Calculcate all the GC values for all windows
    /////////////////////////////////////////////////////////////////////////////
    private byte[] calculateAllGcs(final byte[] refBases, final int[] windowsByGc, final int lastWindowStart) {
        final CalculateGcState state = new CalculateGcState();
        final int refLength = refBases.length;
        final byte[] gc = new byte[refLength + 1];
        for (int i = 1; i < lastWindowStart; ++i) {
            final int windowEnd = i + windowSize;
            final int windowGc = calculateGc(refBases, i, windowEnd, state);
            gc[i] = (byte) windowGc;
            if (windowGc != -1) windowsByGc[windowGc]++;
        }
        return gc;
    }

    /////////////////////////////////////////////////////////////////////////////
    // Calculates GC as a number from 0 to 100 in the specified window.
    // If the window includes more than five no-calls then -1 is returned.
    /////////////////////////////////////////////////////////////////////////////
    private int calculateGc(final byte[] bases, final int startIndex, final int endIndex, final CalculateGcState state) {
        if (state.init) {
            state.init = false;
            state.gcCount = 0;
            state.nCount = 0;
            for (int i = startIndex; i < endIndex; ++i) {
                final byte base = bases[i];
                if (base == 'G' || base == 'C') ++state.gcCount;
                else if (base == 'N') ++state.nCount;
            }
        } else {
            final byte newBase = bases[endIndex - 1];
            if (newBase == 'G' || newBase == 'C') ++state.gcCount;
            else if (newBase == 'N') ++state.nCount;

            if (state.priorBase == 'G' || state.priorBase == 'C') --state.gcCount;
            else if (state.priorBase == 'N') --state.nCount;
        }
        state.priorBase = bases[startIndex];
        if (state.nCount > 4) return -1;
        else return (state.gcCount * 100) / (endIndex - startIndex);
    }

    /////////////////////////////////////////////////////////////////////////////
    // Keeps track of current GC calculation state
    /////////////////////////////////////////////////////////////////////////////
    class CalculateGcState {
        boolean init = true;
        int nCount;
        int gcCount;
        byte priorBase;
    }
}


