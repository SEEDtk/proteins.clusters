/**
 *
 */
package org.theseed.proteins.cluster;

import java.util.ArrayList;
import java.util.Collection;

import org.theseed.counters.CountMap;
import org.theseed.counters.PairCounter;
import org.theseed.genomes.Contig;
import org.theseed.genomes.Feature;
import org.theseed.genomes.FeatureList;
import org.theseed.genomes.Genome;
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
    public Collection<CountMap<Role>.Count> getRoleCounts() {
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
     * @return a collection of counts for the role pairs that appeared, sorted from most frequent to least
     */
    public Collection<PairCounter<Role>.Count> getPairCounts() {
        return this.roleCounts.sortedCounts();
    }

    /**
     * @return the role object for the role with the specified ID, or NULL if it does not exist
     */
    public Role getRole(String roleId) {
        return this.usefulRoles.get(roleId);
    }

}
