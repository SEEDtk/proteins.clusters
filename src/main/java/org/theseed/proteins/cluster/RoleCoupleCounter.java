/**
 *
 */
package org.theseed.proteins.cluster;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Scanner;
import java.util.regex.Pattern;

import org.theseed.counters.CountMap;
import org.theseed.counters.PairCounter;
import org.theseed.genome.Contig;
import org.theseed.genome.Feature;
import org.theseed.genome.FeatureList;
import org.theseed.genome.Genome;
import org.theseed.proteins.Role;
import org.theseed.proteins.RoleMap;

/**
 * This class counts functionally-coupled roles in genomes.  After a sufficient number of genomes have
 * been counted, it can be used to produce a list of useful functional couplings.
 *
 * @author Bruce Parrello
 *
 */
public class RoleCoupleCounter {

    // FIELDS
    /** counting map for coupled roles */
    PairCounter<Role> roleCounts;
    /** table of interesting roles */
    RoleMap usefulRoles;
    /** maximum gap between neighbors */
    int gap;

    /** delimiter pattern for the scanner */
    public static final Pattern DELIM_PATTERN = Pattern.compile("\\t|\\r?\\n");

    /**
     * Create a role-coupling counter for the specified roles with the specified feature gap.
     *
     * @param goodRoles	role map containing the acceptable roles
     * @param gap		maximum distance between neighboring roles
     */
    public RoleCoupleCounter(RoleMap goodRoles, int gap) {
        this.roleCounts = new PairCounter<Role>();
        this.usefulRoles = goodRoles;
        this.gap = gap;
    }

    /**
     * Create a role-coupling counter with no predefined roles and the specified feature gap.
     *
     * @param gap	maximum distance between neighboring roles
     */
    protected RoleCoupleCounter(int gap) {
        this.roleCounts = new PairCounter<Role>();
        this.usefulRoles = new RoleMap();
        this.gap = gap;
    }

    /**
     * Store a useful role with a specified count.
     *
     * @param roleId	ID of the role
     * @param roleName	name of the role
     * @param count		number of occurrences of the role
     */
    protected void addRole(String roleId, String roleName, int count) {
        Role newRole = new Role(roleId, roleName);
        this.usefulRoles.register(newRole);
        this.roleCounts.recordOccurrences(newRole, count);
    }

    /**
     * Store an occurrence count for a pair of roles.
     *
     * @param roleId1	ID of the first role
     * @param roleId2	ID of the second role
     * @param count		number of occurrences together
     */
    protected void addPair(String roleId1, String roleId2, int count) {
        Role role1 = this.getRoleById(roleId1);
        Role role2 = this.getRoleById(roleId2);
        this.roleCounts.recordPairings(role1, role2, count);
    }

    /**
     * @return the role with the specified ID
     *
     * Unlike getRole, this throws an error if the role is not found.
     *
     * @param roleId	ID of the desired role
     */
    private Role getRoleById(String roleId) {
        Role retVal = this.getRole(roleId);
        if (retVal == null)
            throw new IllegalArgumentException("No role found with ID " + roleId + ".");
        return retVal;
    }

    /**
     * Count the couplings in the specified genome.
     *
     * @param genome	genome of interest
     */
    public void countCouplings(Genome genome) {
        for (Contig contig : genome.getContigs()) {
            FeatureList contigFeatures = genome.getContigFeatures(contig.getId());
            // Use the position object to move through the features in the contig.
            FeatureList.Position position = contigFeatures.new Position();
            while (position.hasNext()) {
                Feature current = position.next();
                // Get this feature's roles.
                Collection<Role> currentRoles = current.getUsefulRoles(this.usefulRoles);
                if (currentRoles.size() > 0) {
                    // We have useful roles.  Find the other features in the area.
                    Collection<Feature> neighbors = position.within(this.gap);
                    // Extract the other feature roles.
                    ArrayList<Role> neighborRoles = new ArrayList<Role>(neighbors.size() + 4);
                    for (Feature feat : neighbors) {
                        neighborRoles.addAll(feat.getUsefulRoles(this.usefulRoles));
                    }
                    // Now we loop through our roles, counting the neighbor roles.
                    for (Role role : currentRoles) {
                        roleCounts.recordOccurrence(role, neighborRoles);
                    }
                }
            }
        }
    }

    /**
     * @return a collection of counts for the roles that appeared, sorted from most frequent to least
     */
    public List<CountMap<Role>.Count> getRoleCounts() {
        return this.roleCounts.sortedItemCounts();
    }

    /**
     * @return the number of times a role appeared
     *
     * @param role	role of interest
     */
    public int getCount(Role role) {
        return this.roleCounts.getCount(role);
    }

    /**
     * @return the number of times a pair of roles appeared together
     *
     * @param role1	first role of interest
     * @param role2	second role of interest
     */
    public int getCount(Role role1, Role role2) {
        return this.roleCounts.getCount(role1, role2);
    }

    /**
     * @return a list of counts for the role pairs that appeared, sorted from most frequent to least
     */
    public List<PairCounter<Role>.Count> getPairCounts() {
        return this.roleCounts.sortedCounts(0.0, 0);
    }

    /**
     * @return a list of counts for the role pairs with a togetherness and count higher than
     * 		   the specified amounts, sorted from most frequent to least frequent
     *
     * @param minTogether	minimum acceptable togetherness fraction
     * @param minCount		minimum acceptable occurrence count
     */
    public List<PairCounter<Role>.Count> getPairCounts(double minTogether, int minCount) {
        return this.roleCounts.sortedCounts(minTogether, minCount);
    }

    /**
     * @return the strength of the specified coupling, that is, the fraction of times the two roles
     * 		   occur together
     *
     * @param role1	first role of interest
     * @param role2 second role of interest
     */
    public double getTogetherness(Role role1, Role role2) {
        PairCounter<Role>.Count counter = this.roleCounts.getPairCount(role1, role2);
        double retVal = 0.0;
        if (counter != null) {
            retVal = counter.togetherness();
        }
        return retVal;
    }

    /**
     * @return the role object for the role with the specified ID, or NULL if it does not exist
     */
    public Role getRole(String roleId) {
        return this.usefulRoles.getItem(roleId);
    }

    /**
     * Save this role-coupling counter to the specified file.  The counter data is saved as a text
     * file to permit easy manipulation in other languages. The basic file format consists of a
     * header record with the gap, a line of column labels, then one record per useful role (count, ID, name),
     * a line of more column labels, then one record per role pair (count, ID1, ID2, togetherness).
     *
     * @param outFile	output file
     *
     * @throws IOException
     */
    public void save(File outFile) throws IOException {
        PrintWriter writer = new PrintWriter(outFile);
        writer.format("%d\tRole-Coupling Database%n", this.gap);
        writer.println("count\trole_id\trole_name%n");
        Collection<Role> roles = this.usefulRoles.objectValues();
        for (Role role : roles) {
            writer.format("%d\t%s\t%s%n", this.getCount(role), role.getId(), role.getName());
        }
        writer.println("role1_id\trole2_id\tcount\ttogetherness");
        List<PairCounter<Role>.Count> sortedCounts = this.getPairCounts();
        for (PairCounter<Role>.Count count : sortedCounts) {
            writer.format("%d\t%s\t%s\t%4.2g%n", count.getCount(), count.getKey1().getId(),
                    count.getKey2().getId(), count.togetherness());
        }
        writer.close();
    }

    /**
     * Load a role-coupling counter from the specified file.
     *
     * @param inFile	input file
     *
     * @return a new role-coupling counter read from the file
     *
     * @throws IOException
     */
    public static RoleCoupleCounter load(File inFile) throws IOException {
        Scanner reader = new Scanner(inFile);
        reader.useDelimiter(DELIM_PATTERN);
        // Read the gap number and skip the rest of the heading line
        int gap = reader.nextInt();
        RoleCoupleCounter retVal = new RoleCoupleCounter(gap);
        reader.nextLine();
        // Skip the role table header.
        reader.nextLine();
        // Loop through the role table.
        while (reader.hasNextInt()) {
            int count = reader.nextInt();
            String roleId = reader.next();
            String roleDesc = reader.next();
            retVal.addRole(roleId, roleDesc, count);
        }
        // Skip the end-of-line token on the last line.
        reader.nextLine();
        // Skip the pair table header.
        reader.nextLine();
        // Loop through the pair table.
        while (reader.hasNextInt()) {
            int count = reader.nextInt();
            String roleId1 = reader.next();
            String roleId2 = reader.next();
            retVal.addPair(roleId1, roleId2, count);
            // Skip the togetherness number.
            reader.nextLine();
        }
        reader.close();
        return retVal;
    }

}
