package org.theseed.proteins.cluster;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.theseed.counters.PairCounter;
import org.theseed.genome.Genome;
import org.theseed.genome.GenomeDirectory;
import org.theseed.proteins.Role;
import org.theseed.proteins.RoleMap;

/**
 * This is the primary class for processing role couplings.  It reads genomes from a
 * directory, runs them through the coupler, and then saves the coupler to a file.
 * It accepts the following command-line parameters.
 *
 * -t			minimum togetherness threshold for two features to be considered
 * 				coupled (default 0.80)
 * -m			minimum number of occurrences for a pair to be considered coupled
 * 				(default 10)
 * -v			write progress messages to STDERR
 *
 * --create		create a new coupling table in the specified coupler file; otherwise,
 * 				if the coupler file exists it will be loaded before the new genome
 * 				directory is processed
 *
 * --compare	compare the coupling results to the coupler in a specified file;
 * 				for every coupling in our results file that exceeds the threshold,
 * 				we want to insure it exceeds the threshold in the specified file
 *
 * If "--create" is specified, the following options are required; otherwise they
 * are ignored.
 *
 * -g			the maximum distance allowed for two features to be considered
 * 				neighbors (default 500)
 * -R			name of a file containing the useful roles; the file is tab-delimited,
 * 				each record containing a role ID and a role name
 *
 * If "--compare" is specified, the following options are also used.
 *
 * -n			minimum number of role occurrences for a pair to be considered in
 * 				error (default 20)
 * -u			minimum togetherness threshold for two features to be considered
 * 				coupled in the comparison (default 0.70)
 *
 * The positional parameters are the name of the coupler file and the name of a
 * directory containing the genomes to process.  If "--create" is specified, the
 * coupler file will be used for output only; otherwise, it will be read in to
 * initialize the coupler and the new genomes will be incorporated into it before
 * the coupler is written back out.
 *
 * @author Bruce Parrello
 *
 */
public class RoleCouplingProcessor {

    // FIELDS
    /** current coupler database */
    private RoleCoupleCounter coupler;
    /** role map for create mode */
    private RoleMap roleSet;

    // COMMAND LINE

    /** help option */
    @Option(name="-h", aliases={"--help"}, help=true)
    private boolean help;

    /** TRUE if we want progress messages */
    @Option(name="-v", aliases= {"--verbose", "--debug"}, usage="display progress on STDERR")
    private boolean debug;

    /** maximum distance between neighboring features */
    @Option(name="-g", aliases={"--gap"}, metaVar="500",
            usage="maximum distance between neighboring features (create only)")
    private int gap;

    /** file containing the useful roles */
    @Option(name="-R", aliases={"--roles", "--roleFile"}, metaVar="roleFile",
            usage="file of useful roles (create only)")
    private File roleFile;

    /** minimum togetherness threshold for output report */
    @Option(name="-t", aliases={"--minScore", "--minStrength"}, metaVar="0.80",
            usage="minimum fraction of times coupled features found together")
    private double togetherness;

    /** minimum togetherness threshold for comparison table */
    @Option(name="-u", aliases={"--compareScore", "--compareStrength"}, metaVar="0.70",
            usage="minimum fraction of times coupled features found together in comparison")
    private double otherTogetherness;

    /** minimum number of occurrences for output report */
    @Option(name="-m", aliases={"--minCount", "--minFreq"}, metaVar="10",
            usage="minimum number of times coupled features found together")
    private int minCount;

    /** minimum number of occurrences for output report */
    @Option(name="-n", aliases={"--compareCount", "--comareFreq"}, metaVar="20",
            usage="minimum number of times coupled features found for valid comparison")
    private int otherCount;

    /** creation / reuse flag */
    @Option(name="--create", usage="create new coupler file")
    private boolean createMode;

    /** compare option */
    @Option(name="--compare", usage="compare results to another coupler file")
    private File compareFile;

    /** coupler file */
    @Argument(index=0, metaVar="coupler.ser", usage="name of coupler file",
            required=true)
    private File couplerFile;

    /** input genome directories */
    @Argument(index=1, metaVar="genomeDir1 genomeDir2 ...", multiValued=true,
            usage="directory of input genomes")
    private List<File> genomeDirs;

    /**
     * Parse command-line options to specify the parameters of this object.
     *
     * @param args	an array of the command-line parameters and options
     *
     * @return TRUE if successful, FALSE if the parameters are invalid
     */
    public boolean parseCommand(String[] args) {
        boolean retVal = false;
        // Set the defaults.
        this.gap = 500;
        this.debug = false;
        this.roleFile = null;
        this.compareFile = null;
        this.togetherness = 0.80;
        this.minCount = 10;
        this.otherCount = 20;
        this.otherTogetherness = 0.70;
        this.createMode = false;
        this.genomeDirs = new ArrayList<File>();
        // Parse the command line.
        CmdLineParser parser = new CmdLineParser(this);
        try {
            parser.parseArgument(args);
            if (this.help) {
                parser.printUsage(System.err);
            } else {
                // We need to do a sanity check here.
                if (! this.createMode) {
                    // If we are not creating, then the coupler file must exist.
                    if (! couplerFile.exists()) {
                        throw new FileNotFoundException("Coupler file must exist unless --create is specified.");
                    }
                } else {
                    // If we are creating, then get the role file.
                    if (this.roleFile == null) {
                        throw new IllegalArgumentException("Role file required in create mode.");
                    } else {
                        this.roleSet = RoleMap.load(this.roleFile);
                    }
                }
                // Insure the genome directories are valid.
                for (File genomeDir : genomeDirs) {
                    if (! genomeDir.isDirectory()) {
                        throw new FileNotFoundException(genomeDir.getPath() + " is not a valid directory.");
                    }
                }
                // Insure the compare file exists.
                if (this.compareFile != null && ! this.compareFile.exists()) {
                    throw new FileNotFoundException(this.compareFile.getPath() + " is not found.");
                }
                // We made it this far, we can run the application.
                retVal = true;
            }
        } catch (CmdLineException e) {
            System.err.println(e.getMessage());
            // For parameter errors, we display the command usage.
            parser.printUsage(System.err);
        } catch (IOException e) {
            System.err.println(e.getMessage());
        }
        return retVal;
    }

    /**
     * Excute the operation indicated by the parameters.
     */
    public void run() {
        try {
            // Start by loading the coupler.  We create a blank if we are in create mode; otherwise we
            // read it.
            if (this.createMode) {
                if (debug) System.err.println("Initializing new coupling counter.");
                this.coupler = new RoleCoupleCounter(this.roleSet, this.gap);
            } else {
                if (debug) System.err.println("Loading coupling counter from " + this.couplerFile.getPath() + ".");
                this.coupler = RoleCoupleCounter.load(this.couplerFile);
            }
            // Loop through the genome directories, adding their data to the coupling counts.
            for (File dirFile : this.genomeDirs) {
                if (debug) System.err.println("Processing genome directory " + dirFile.getPath() + ".");
                GenomeDirectory genomeDir = new GenomeDirectory(dirFile.getPath());
                for (Genome genome : genomeDir) {
                    if (debug) System.err.println("Parsing genome " + genome + ".");
                    this.coupler.countCouplings(genome);
                }
            }
            // Save the coupling data.
            if (debug) System.err.println("Saving coupling data to " + this.couplerFile.getPath() + ".");
            this.coupler.save(this.couplerFile);
            // If we have a comparison file, read it in here.  Note we will also track the number of
            // compare failures.
            int failureCount = 0;
            int totalCount = 0;
            RoleCoupleCounter comparator = null;
            if (this.compareFile != null) {
                if (debug) System.err.println("Loading comparator from " + this.compareFile + ".");
                comparator = RoleCoupleCounter.load(this.compareFile);
            }
            // Write the data that meets the thresholds.
            if (debug) System.err.println("Writing output.");
            // Compute the headers based on whether or not there is a comparison.
            if (comparator != null) {
                System.out.println("role_id1\trole_id2\tfraction\tcount\totherFrac\totherCount\totherFound\tfailure");
            } else {
                System.out.println("role_id1\trole_id2\tfraction\tcount");
            }
            List<PairCounter<Role>.Count> goodPairs = this.coupler.getPairCounts(this.togetherness, this.minCount);
            for (PairCounter<Role>.Count goodPair : goodPairs) {
                System.out.format("%s\t%s\t%4.2g\t%d", goodPair.getKey1().getId(), goodPair.getKey2().getId(),
                        goodPair.togetherness(), goodPair.getCount());
                totalCount++;
                if (comparator == null) {
                    // No comparison.  End the line.
                    System.out.println();
                } else {
                    // Comparison required.  Get the count.  Note we ignore the result if neither
                    // role appeared in the sample.
                    int otherPairCount = comparator.getCount(goodPair.getKey1(), goodPair.getKey2());
                    int appearances = comparator.getCount(goodPair.getKey1()) +
                            comparator.getCount(goodPair.getKey2());
                    if (appearances > 0) {
                        double otherFrac = comparator.getTogetherness(goodPair.getKey1(), goodPair.getKey2());
                        int otherRoleCount = appearances - otherPairCount;
                        String errorFlag = "";
                        if (otherFrac < this.otherTogetherness && appearances >= this.otherCount) {
                            errorFlag = "Y";
                            failureCount++;
                        }
                        System.out.format("\t%4.2g\t%d\t%d\t%s%n", otherFrac, otherPairCount, otherRoleCount,
                                errorFlag);
                    } else {
                        System.out.println();
                    }

                }
            }
            if (debug && comparator != null) {
                System.err.println("Failure count for comparison is " + failureCount + ".");
            }
            if (debug) System.err.println(totalCount + " couplings found.");
        } catch (IOException e) {
            // Percolate the error.
            System.err.println("Error processing command: " + e.getMessage());
        }
    }

}
