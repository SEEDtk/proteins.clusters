package org.theseed.proteins.cluster;

import java.io.File;
import java.util.List;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.theseed.genomes.GenomeDirectory;

/**
 * This is the primary class for processing role couplings.  It reads genomes from a
 * directory, runs them through the coupler, and then saves the coupler to a file.
 * It accepts the following command-line parameters.
 *
 * -t			minimum togetherness threshold for two features to be considered
 * 				coupled (default 0.80)
 * -m			minimum number of occurrences for a pair to be considered coupled
 * 				(default 10)
 *
 * --create		create a new coupling table in the specified coupler file; otherwise,
 * 				if the coupler file exists it will be loaded before the new genome
 * 				directory is processed
 *
 * If "--create" is specified, the following options are required; otherwise they
 * are ignored.
 *
 * -g			the maximum distance allowed for two features to be considered
 * 				neighbors (default 500)
 * -R			name of a file containing the useful roles; the file is tab-delimited,
 * 				each record containing a role ID and a role name
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
    /** current directory of input genomes */
    private GenomeDirectory inputDir;
    /** current coupler database */
    private RoleCoupleCounter coupler;

    // COMMAND LINE

    /** help option */
    @Option(name="-h", aliases={"--help"}, help=true)
    private boolean help;

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

    /** minimum number of occurrences for output report */
    @Option(name="-m", aliases={"--minCount", "--minFreq"}, metaVar="10",
            usage="minimum number of times coupled features found together")
    private int minCount;

    /** creation / reuse flag */
    @Option(name="--create", usage="create new coupler file")
    private boolean createMode;

    /** coupler file */
    @Argument(index=0, metaVar="coupler.ser", usage="name of coupler file",
            required=true)
    private File couplerFile;

    /** input genome directories */
    @Argument(index=1, metaVar="genomeDir1 genomeDir2 ...", multiValued=true,
            usage="directory of input genomes")
    private List<File> genomeDirs;


}
