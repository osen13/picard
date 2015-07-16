package picard.analysis.artifacts;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.filter.AggregateFilter;
import htsjdk.samtools.filter.AlignedFilter;
import htsjdk.samtools.filter.DuplicateReadFilter;
import htsjdk.samtools.filter.FailsVendorReadQualityFilter;
import htsjdk.samtools.filter.InsertSizeFilter;
import htsjdk.samtools.filter.MappingQualityFilter;
import htsjdk.samtools.filter.NotPrimaryAlignmentFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.IntervalListReferenceSequenceMask;
import htsjdk.samtools.util.StringUtil;
import picard.PicardException;
import picard.analysis.SinglePassSamProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.programgroups.Metrics;
import picard.util.DbSnpBitSetUtil;
import picard.analysis.artifacts.SequencingArtifactMetrics.*;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import static htsjdk.samtools.util.CodeUtil.getOrElse;

/**
 * Quantify substitution errors caused by mismatched base pairings during various
 * stages of sample / library prep.
 *
 * We measure two distinct error types - artifacts that are introduced before
 * the addition of the read1/read2 adapters ("pre adapter") and those that are
 * introduced after target selection ("bait bias"). For each of these, we provide
 * summary metrics as well as detail metrics broken down by reference context
 * (the ref bases surrounding the substitution event).
 *
 * For a deeper explanation, see Costello et al. 2013:
 * http://www.ncbi.nlm.nih.gov/pubmed/23303777
 *
 * @author mattsooknah
 *
 */
@CommandLineProgramProperties(
        usage = CollectSequencingArtifactMetrics.USAGE_SUMMARY + CollectSequencingArtifactMetrics.USAGE_DETAILS,
        usageShort = CollectSequencingArtifactMetrics.USAGE_SUMMARY,
        programGroup = Metrics.class
)
public class CollectSequencingArtifactMetrics extends SinglePassSamProgram {
    static final String USAGE_SUMMARY = "Collect metrics to quantify single-base sequencing artifacts.";
    static final String USAGE_DETAILS = "This tool examines two sources of sequencing errors; bait-bias and pre-adapter artifacts resulting " +
            "from hybrid selection protocols.  The hybrid selection platform enables selection of specific sequences from a pool" +
            " of genomic DNA for targeted sequencing analyses via pull-down assays (doi:10.1038/nbt.1523).  " +
            "Typical applications include the selection of exome sequences or pathogen-specific sequences in complex biological samples.  " +
            "Briefly, baits are RNA (or sometimes DNA) molecules synthesized with biotinylated nucleotides.  " +
            "The biotinylated nucleotides are ligands for streptavidin enabling enabling RNA:DNA hybrids to be captured in solution.  " +
            "The hybridization targets are sheared genomic DNA fragments, which have been \"polished\" with synthetic" +
            " adapters to facilitate PCR cloning downstream." +
            "  Hybridization of the baits with the denatured targets is followed by selective capture of the" +
            " RNA:DNA “hybrids” using streptavidin-coated beads via pull-down assays or columns.<br /><br />" +
            "" +
            "" +
            "Systematic errors, ultimately leading to sequence bias and incorrect variant calls, can arise at" +
            " several steps.  For example, the shearing of target genomic DNA leading to oxidation of an amine" +
            " of guanine at position 8 \"8-oxoguanine\" (doi:10.1093/nar/gks1443).  This is considered a pre-adapter " +
            "artifact since it can arise as" +
            " a result of DNA shearing and prior to the ligation of the PCR adapters.  These artifacts occur on the" +
            " original template strand, before the addition of adapters, so they correlate with read number" +
            " orientation in a specific way.  For example, the well-known \"8-oxoguanine\" artifact occurs when a G" +
            " on the template strand is oxidized, giving it an affinity for binding to A rather than the usual C." +
            " Thus, PCR will introduce apparent G>T substitutions in read 1 and C>A in read 2. In the resulting" +
            " alignments, a given G>T or C>A observation could either be: 1. a true mutation 2. an 8-oxoguanine artifact" +
            " 3. some other kind of artifact.  On average, we assume that 1 and 3 will not display this read" +
            " number / orientation bias, so their contributions will cancel out in the calculation. <br /><br />" +
            "" +
            "A second type of bias is known as a single bait bias or a reference bias artifact. These artifacts" +
            " occur during or after the target selection step, and correlate with substitution rates that are" +
            " \"biased\", or higher for sites having one base on the reference/positive strand relative to" +
            " sites having the complementary base on that strand.  For example, a G>T artifact during the target" +
            " selection step might result in a higher (G>T)/(C>A) substitution rate at sites with a G on the" +
            " positive strand (and C on the negative), relative to sites with the flip (C positive)/(G negative)." +
            " This is known as the \"G-Ref\" artifact. <br /><br />" +
            "" +
            "This tool produces four files; summary and detail metrics files for both pre-adapter and " +
            "bait-bias artifacts.  The detail metrics show the error rates for each type of base" +
            " substitution within every possible triplet base configuration.  Error rates associated with these" +
            " substitutions are Phred-scaled and provided as quality scores, the lower the value, the more" +
            " likely it is that an alternate base call is due to an artifact.  The summary metrics provide likelihood " +
            "information on the \"worst-case\" errors. <br />" +
            "" +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar CollectSequencingArtifactMetrics \\<br />" +
            "     -I=Input.bam \\<br />" +
            "     -O=Artifactmetrics.txt \\<br />" +
            "     -R=HumanSequence.fasta" +
            "</pre>" +
            "" +
            "For additional information, please see" +
            "http://broadinstitute.github.io/picard/picard-metric-definitions.html#SequencingArtifactMetrics." +
            "PreAdapterDetailMetrics <br />" +
            "http://broadinstitute.github.io/picard/picard-metric-definitions.html#SequencingArtifactMetrics." +
            "PreAdapterSummaryMetrics <br />" +
            "http://broadinstitute.github.io/picard/picard-metric-definitions.html#SequencingArtifactMetrics." +
            "BaitBiasDetailMetrics <br />" +
            "http://broadinstitute.github.io/picard/picard-metric-definitions.html#SequencingArtifactMetrics." +
            "BaitBiasSummaryMetrics <br />" +
            "<hr />" ;
    @Option(doc = "An optional list of intervals to restrict analysis to.", optional = true)
    public File INTERVALS;

    @Option(doc = "VCF format dbSNP file, used to exclude regions around known polymorphisms from analysis.", optional = true)
    public File DB_SNP;

    @Option(shortName = "Q", doc = "The minimum base quality score for a base to be included in analysis.")
    public int MINIMUM_QUALITY_SCORE = 20;

    @Option(shortName = "MQ", doc = "The minimum mapping quality score for a base to be included in analysis.")
    public int MINIMUM_MAPPING_QUALITY = 30;

    @Option(shortName = "MIN_INS", doc = "The minimum insert size for a read to be included in analysis.")
    public int MINIMUM_INSERT_SIZE = 60;

    @Option(shortName = "MAX_INS", doc = "The maximum insert size for a read to be included in analysis. Set to 0 to have no maximum.")
    public int MAXIMUM_INSERT_SIZE = 600;

    @Option(shortName = "UNPAIRED", doc = "Include unpaired reads. If set to true then all paired reads will be included as well - " +
            "MINIMUM_INSERT_SIZE and MAXIMUM_INSERT_SIZE will be ignored.")
    public boolean INCLUDE_UNPAIRED = false;

    @Option(shortName = "TANDEM", doc = "Set to true if mate pairs are being sequenced from the same strand, " +
            "i.e. they're expected to face the same direction.")
    public boolean TANDEM_READS = false;

    @Option(doc = "When available, use original quality scores for filtering.")
    public boolean USE_OQ = true;

    @Option(doc = "The number of context bases to include on each side of the assayed base.")
    public int CONTEXT_SIZE = 1;

    @Option(doc = "If specified, only print results for these contexts in the detail metrics output. " +
                  "However, the summary metrics output will still take all contexts into consideration.")
    public Set<String> CONTEXTS_TO_PRINT = new HashSet<String>();

    private static final String UNKNOWN_LIBRARY = "UnknownLibrary";
    private static final String UNKNOWN_SAMPLE = "UnknownSample";

    private File preAdapterSummaryOut;
    private File preAdapterDetailsOut;
    private File baitBiasSummaryOut;
    private File baitBiasDetailsOut;

    private IntervalListReferenceSequenceMask intervalMask;
    private DbSnpBitSetUtil dbSnpMask;
    private SamRecordFilter recordFilter;

    private final Set<String> samples = new HashSet<String>();
    private final Set<String> libraries = new HashSet<String>();
    private final Map<String, ArtifactCounter> artifactCounters = new HashMap<String, ArtifactCounter>();

    public static void main(final String[] args) {
        new CollectSequencingArtifactMetrics().instanceMainWithExit(args);
    }

    @Override
    protected String[] customCommandLineValidation() {
        final List<String> messages = new ArrayList<String>();

        final int contextFullLength = 2 * CONTEXT_SIZE + 1;
        if (CONTEXT_SIZE < 0) messages.add("CONTEXT_SIZE cannot be negative");
        for (final String context : CONTEXTS_TO_PRINT) {
            if (context.length() != contextFullLength) {
                messages.add("Context " + context + " is not the length implied by CONTEXT_SIZE: " + contextFullLength);
            }
        }

        if (MINIMUM_INSERT_SIZE < 0) messages.add("MINIMUM_INSERT_SIZE cannot be negative");
        if (MAXIMUM_INSERT_SIZE < 0) messages.add("MAXIMUM_INSERT_SIZE cannot be negative");
        if (MAXIMUM_INSERT_SIZE > 0 && MAXIMUM_INSERT_SIZE < MINIMUM_INSERT_SIZE) {
            messages.add("MAXIMUM_INSERT_SIZE cannot be less than MINIMUM_INSERT_SIZE unless set to 0");
        }

        return messages.isEmpty() ? null : messages.toArray(new String[messages.size()]);
    }

    @Override
    protected void setup(final SAMFileHeader header, final File samFile) {
        preAdapterSummaryOut = new File(OUTPUT + SequencingArtifactMetrics.PRE_ADAPTER_SUMMARY_EXT);
        preAdapterDetailsOut = new File(OUTPUT + SequencingArtifactMetrics.PRE_ADAPTER_DETAILS_EXT);
        baitBiasSummaryOut = new File(OUTPUT + SequencingArtifactMetrics.BAIT_BIAS_SUMMARY_EXT);
        baitBiasDetailsOut = new File(OUTPUT + SequencingArtifactMetrics.BAIT_BIAS_DETAILS_EXT);

        IOUtil.assertFileIsWritable(preAdapterSummaryOut);
        IOUtil.assertFileIsWritable(preAdapterDetailsOut);
        IOUtil.assertFileIsWritable(baitBiasSummaryOut);
        IOUtil.assertFileIsWritable(baitBiasDetailsOut);

        for (final SAMReadGroupRecord rec : header.getReadGroups()) {
            samples.add(getOrElse(rec.getSample(), UNKNOWN_SAMPLE));
            libraries.add(getOrElse(rec.getLibrary(), UNKNOWN_LIBRARY));
        }

        if (INTERVALS != null) {
            IOUtil.assertFileIsReadable(INTERVALS);
            intervalMask = new IntervalListReferenceSequenceMask(IntervalList.fromFile(INTERVALS).uniqued());
        }

        if (DB_SNP != null) {
            IOUtil.assertFileIsReadable(DB_SNP);
            dbSnpMask = new DbSnpBitSetUtil(DB_SNP, header.getSequenceDictionary());
        }

        // set record-level filters
        final List<SamRecordFilter> filters = new ArrayList<SamRecordFilter>();
        filters.add(new FailsVendorReadQualityFilter());
        filters.add(new NotPrimaryAlignmentFilter());
        filters.add(new DuplicateReadFilter());
        filters.add(new AlignedFilter(true)); // discard unmapped reads
        filters.add(new MappingQualityFilter(MINIMUM_MAPPING_QUALITY));
        if (!INCLUDE_UNPAIRED) {
            final int effectiveMaxInsertSize = (MAXIMUM_INSERT_SIZE == 0) ? Integer.MAX_VALUE : MAXIMUM_INSERT_SIZE;
            filters.add(new InsertSizeFilter(MINIMUM_INSERT_SIZE, effectiveMaxInsertSize));
        }
        recordFilter = new AggregateFilter(filters);

        // set up the artifact counters
        final String sampleAlias = StringUtil.join(",", new ArrayList<String>(samples));
        for (final String library : libraries) {
            artifactCounters.put(library, new ArtifactCounter(sampleAlias, library, CONTEXT_SIZE, TANDEM_READS));
        }
    }

    @Override
    protected void acceptRead(final SAMRecord rec, final ReferenceSequence ref) {
        // see if the whole read should be skipped
        if (recordFilter.filterOut(rec)) return;

        // check read group + library
        final String library = (rec.getReadGroup() == null) ? UNKNOWN_LIBRARY : getOrElse(rec.getReadGroup().getLibrary(), UNKNOWN_LIBRARY);
        if (!libraries.contains(library)) {
            // should never happen if SAM is valid
            throw new PicardException("Record contains library that is missing from header: " + library);
        }

        // iterate over aligned positions
        for (final AlignmentBlock block : rec.getAlignmentBlocks()) {
            for (int offset = 0; offset < block.getLength(); offset++) {
                // remember, these are 1-based!
                final int readPos = block.getReadStart() + offset;
                final int refPos = block.getReferenceStart() + offset;

                /**
                 * Skip regions outside of intervals.
                 *
                 * NB: IntervalListReferenceSequenceMask.get() has side-effects which assume
                 * that successive ReferenceSequence's passed to this method will be in-order
                 * (e.g. it will break if you call acceptRead() with chr1, then chr2, then chr1
                 * again). So this only works if the underlying iteration is coordinate-sorted.
                 */
                if (intervalMask != null && !intervalMask.get(ref.getContigIndex(), refPos)) continue;

                // skip dbSNP sites
                if (dbSnpMask != null && dbSnpMask.isDbSnpSite(ref.getName(), refPos)) continue;

                // skip the ends of the reference
                final int contextStartIndex = refPos - CONTEXT_SIZE - 1;
                final int contextFullLength = 2 * CONTEXT_SIZE + 1;
                if (contextStartIndex < 0 || contextStartIndex + contextFullLength > ref.length()) continue;

                // skip contexts with N bases
                final String context = StringUtil.bytesToString(ref.getBases(), contextStartIndex, contextFullLength).toUpperCase();
                if (context.contains("N")) continue;

                // skip low BQ sites
                if (failsBaseQualityCutoff(readPos, rec)) continue;

                // skip N bases in read
                final char readBase = Character.toUpperCase((char) rec.getReadBases()[readPos - 1]);
                if (readBase == 'N') continue;

                // count the base!
                artifactCounters.get(library).countRecord(context, readBase, rec);
            }
        }
    }

    @Override
    protected void finish() {
        final MetricsFile<PreAdapterSummaryMetrics, Integer> preAdapterSummaryMetricsFile = getMetricsFile();
        final MetricsFile<PreAdapterDetailMetrics, Integer> preAdapterDetailMetricsFile = getMetricsFile();
        final MetricsFile<BaitBiasSummaryMetrics, Integer> baitBiasSummaryMetricsFile = getMetricsFile();
        final MetricsFile<BaitBiasDetailMetrics, Integer> baitBiasDetailMetricsFile = getMetricsFile();

        for (final ArtifactCounter counter : artifactCounters.values()) {
            // build metrics
            counter.finish();

            // write metrics
            preAdapterSummaryMetricsFile.addAllMetrics(counter.getPreAdapterSummaryMetrics());
            baitBiasSummaryMetricsFile.addAllMetrics(counter.getBaitBiasSummaryMetrics());

            for (final PreAdapterDetailMetrics preAdapterDetailMetrics : counter.getPreAdapterDetailMetrics()) {
                if (CONTEXTS_TO_PRINT.size() == 0 || CONTEXTS_TO_PRINT.contains(preAdapterDetailMetrics.CONTEXT)) {
                    preAdapterDetailMetricsFile.addMetric(preAdapterDetailMetrics);
                }
            }
            for (final BaitBiasDetailMetrics baitBiasDetailMetrics : counter.getBaitBiasDetailMetrics()) {
                if (CONTEXTS_TO_PRINT.size() == 0 || CONTEXTS_TO_PRINT.contains(baitBiasDetailMetrics.CONTEXT)) {
                    baitBiasDetailMetricsFile.addMetric(baitBiasDetailMetrics);
                }
            }

        }

        preAdapterDetailMetricsFile.write(preAdapterDetailsOut);
        preAdapterSummaryMetricsFile.write(preAdapterSummaryOut);
        baitBiasDetailMetricsFile.write(baitBiasDetailsOut);
        baitBiasSummaryMetricsFile.write(baitBiasSummaryOut);
    }

    @Override
    protected boolean usesNoRefReads() { return false; }

    /**
     * Check if this read base fails the base quality cutoff.
     */
    private boolean failsBaseQualityCutoff(final int oneIndexedPos, final SAMRecord rec) {
        final byte qual;
        if (USE_OQ && rec.getOriginalBaseQualities() != null) {
            qual = rec.getOriginalBaseQualities()[oneIndexedPos - 1];
        } else {
            qual = rec.getBaseQualities()[oneIndexedPos - 1];
        }
        return (qual < MINIMUM_QUALITY_SCORE);
    }
}
