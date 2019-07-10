package org.theseed.proteins.cluster;

import junit.framework.Test;

import junit.framework.TestCase;
import junit.framework.TestSuite;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;

import org.theseed.counters.CountMap;
import org.theseed.counters.PairCounter;
import org.theseed.genomes.Contig;
import org.theseed.genomes.Feature;
import org.theseed.genomes.Genome;
import org.theseed.proteins.Role;
import org.theseed.proteins.RoleMap;
/**
 * Unit test for simple App.
 */
public class AppTest
    extends TestCase
{
    /**
     * Create the test case
     *
     * @param testName name of the test case
     */
    public AppTest( String testName )
    {
        super( testName );
    }

    /**
     * @return the suite of tests being tested
     */
    public static Test suite()
    {
        return new TestSuite( AppTest.class );
    }


    /**
     * Test of coupling with a fake genome
     * @throws IOException
     */
    public void testCouples() throws IOException
    {
        Genome fakeGenome = new Genome("12345.6", "Bacillus praestrigiae Narnia", "Bacteria", 11);
        fakeGenome.addContig(new Contig("con1", "agct", 11));
        fakeGenome.addFeature(new Feature("fig|12345.6.peg.1",  "Role 1", "con1", "+",  100,  300));
        fakeGenome.addFeature(new Feature("fig|12345.6.peg.2",  "Role 2", "con1", "-",  100,  400));
        fakeGenome.addFeature(new Feature("fig|12345.6.peg.3",  "Role 3", "con1", "+",  200,  500));
        fakeGenome.addFeature(new Feature("fig|12345.6.peg.4",  "Role 4", "con1", "-", 1000, 1200));
        fakeGenome.addFeature(new Feature("fig|12345.6.peg.5",  "Role 5", "con1", "+", 1010, 1300));
        fakeGenome.addFeature(new Feature("fig|12345.6.peg.6",  "Role 6 / Role 1", "con1", "-", 3300, 4000));
        fakeGenome.addFeature(new Feature("fig|12345.6.peg.7",  "Role 2 # comment", "con1", "-", 5000, 5100));
        fakeGenome.addFeature(new Feature("fig|12345.6.peg.8",  "Role 3 @ Role X", "con1", "+", 5150, 5200));
        fakeGenome.addFeature(new Feature("fig|12345.6.peg.9",  "Role 1", "con1", "+", 5250, 5400));
        fakeGenome.addFeature(new Feature("fig|12345.6.peg.10", "Role 2", "con1", "-", 5401, 5450));
        RoleMap goodRoles = new RoleMap();
        goodRoles.register("Role 1", "Role 2", "Role 3", "Role 4", "Role 5", "Role 6");
        RoleCoupleCounter testCounter = new RoleCoupleCounter(goodRoles, 100);
        testCounter.countCouplings(fakeGenome);
        ArrayList<String> roleIds = new ArrayList<String>(Arrays.asList("invalid", "Role1n1", "Role2n1", "Role3n1", "Role4n1",
                "Role5n1", "Role6n1"));
        Role[] roles = roleIds.stream().map(x -> testCounter.getRole(x)).toArray(Role[]::new);
        assertEquals("Wrong role returned by getRole", "Role1n1", roles[1].getId());
        assertEquals("r1 count wrong", 3, testCounter.getCount(roles[1]));
        assertEquals("r2 count wrong", 3, testCounter.getCount(roles[2]));
        assertEquals("r3 count wrong", 2, testCounter.getCount(roles[3]));
        assertEquals("r4 count wrong", 1, testCounter.getCount(roles[4]));
        assertEquals("r5 count wrong", 1, testCounter.getCount(roles[5]));
        assertEquals("r6 count wrong", 1, testCounter.getCount(roles[6]));
        assertEquals("r1/2 count wrong", 2, testCounter.getCount(roles[1], roles[2]));
        assertEquals("r1/3 count wrong", 2, testCounter.getCount(roles[1], roles[3]));
        assertEquals("r1/6 count wrong", 0, testCounter.getCount(roles[1], roles[6]));
        assertEquals("r2/3 count wrong", 2, testCounter.getCount(roles[2], roles[3]));
        assertEquals("r4/5 count wrong", 1, testCounter.getCount(roles[4], roles[5]));
        Collection<CountMap<Role>.Count> rolesCounted = testCounter.getRoleCounts();
        assertThat("Wrong list of roles counted.",
                rolesCounted.stream().map(x -> x.getKey().getId()).collect(Collectors.toList()),
                containsInAnyOrder("Role1n1", "Role2n1", "Role3n1", "Role4n1",
                        "Role5n1", "Role6n1"));
        int prev = Integer.MAX_VALUE;
        for (CountMap<Role>.Count counter : rolesCounted) {
            assertThat("Counts out of order.", prev, greaterThanOrEqualTo(counter.getCount()));
            assertEquals("Counts do not match.", testCounter.getCount(counter.getKey()), counter.getCount());
            prev = counter.getCount();
        }
        prev = Integer.MAX_VALUE;
        List<PairCounter<Role>.Count> pairsCounted = testCounter.getPairCounts();
        assertEquals("Wrong number of pairs found.", 4, pairsCounted.size());
        for (PairCounter<Role>.Count counter : pairsCounted) {
            assertThat("Counts out of order.", prev, greaterThanOrEqualTo(counter.getCount()));
            assertEquals("Counts do not match.", testCounter.getCount(counter.getKey1(), counter.getKey2()),
                    counter.getCount());
            prev = counter.getCount();
        }
        // Count this genome again.
        testCounter.countCouplings(fakeGenome);
        assertEquals("r1 count wrong on pass 2", 6, testCounter.getCount(roles[1]));
        assertEquals("r2 count wrong on pass 2", 6, testCounter.getCount(roles[2]));
        assertEquals("r3 count wrong on pass 2", 4, testCounter.getCount(roles[3]));
        assertEquals("r4 count wrong on pass 2", 2, testCounter.getCount(roles[4]));
        assertEquals("r5 count wrong on pass 2", 2, testCounter.getCount(roles[5]));
        assertEquals("r6 count wrong on pass 2", 2, testCounter.getCount(roles[6]));
        assertEquals("r1/2 count wrong on pass 2", 4, testCounter.getCount(roles[1], roles[2]));
        assertEquals("r1/3 count wrong on pass 2", 4, testCounter.getCount(roles[1], roles[3]));
        assertEquals("r1/6 count wrong on pass 2", 0, testCounter.getCount(roles[1], roles[6]));
        assertEquals("r2/3 count wrong on pass 2", 4, testCounter.getCount(roles[2], roles[3]));
        assertEquals("r4/5 count wrong on pass 2", 2, testCounter.getCount(roles[4], roles[5]));
        // Rebuilt the genome.
        fakeGenome = new Genome("12345.6", "Bacillus praestrigiae Narnia", "Bacteria", 11);
        fakeGenome.addContig(new Contig("con1", "agct", 11));
        fakeGenome.addFeature(new Feature("fig|12345.6.peg.1",  "Role Y", "con1", "+",  100,  300));
        fakeGenome.addFeature(new Feature("fig|12345.6.peg.2",  "Role A @ Role 1 / Role 2", "con1", "-",  100,  400));
        fakeGenome.addFeature(new Feature("fig|12345.6.peg.3",  "Role 3", "con1", "+",  200,  500));
        fakeGenome.addFeature(new Feature("fig|12345.6.peg.8",  "Role 3 @ Role X", "con1", "+", 5150, 5200));
        fakeGenome.addFeature(new Feature("fig|12345.6.peg.9",  "Role 1", "con1", "+", 5250, 5400));
        fakeGenome.addFeature(new Feature("fig|12345.6.peg.10", "Role 2", "con1", "-", 5401, 5450));
        testCounter.countCouplings(fakeGenome);
        assertEquals("r1 count wrong on pass 3", 8, testCounter.getCount(roles[1]));
        assertEquals("r2 count wrong on pass 3", 8, testCounter.getCount(roles[2]));
        assertEquals("r3 count wrong on pass 3", 6, testCounter.getCount(roles[3]));
        assertEquals("r4 count wrong on pass 3", 2, testCounter.getCount(roles[4]));
        assertEquals("r5 count wrong on pass 3", 2, testCounter.getCount(roles[5]));
        assertEquals("r6 count wrong on pass 3", 2, testCounter.getCount(roles[6]));
        assertEquals("r1/2 count wrong on pass 3", 5, testCounter.getCount(roles[1], roles[2]));
        assertEquals("r1/3 count wrong on pass 3", 6, testCounter.getCount(roles[1], roles[3]));
        assertEquals("r1/6 count wrong on pass 3", 0, testCounter.getCount(roles[1], roles[6]));
        assertEquals("r2/3 count wrong on pass 3", 5, testCounter.getCount(roles[2], roles[3]));
        assertEquals("r4/5 count wrong on pass 3", 2, testCounter.getCount(roles[4], roles[5]));
        Role roleX = testCounter.getRole("RoleX");
        assertEquals("Role X was counted.", 0, testCounter.getCount(roleX));
        assertEquals("r2/3 togetherness wrong.", 0.556, testCounter.getTogetherness(roles[2], roles[3]), 0.001);
        // Test the pair counts with threshold.
        pairsCounted = testCounter.getPairCounts(0.5, 3);
        assertEquals("Wrong number of pairs found.", 2, pairsCounted.size());
        prev = Integer.MAX_VALUE;
        for (PairCounter<Role>.Count counter : pairsCounted) {
            assertThat("Counts out of order.", prev, greaterThanOrEqualTo(counter.getCount()));
            assertEquals("Counts do not match.", testCounter.getCount(counter.getKey1(), counter.getKey2()),
                    counter.getCount());
            assertThat("Min count too low.", counter.getCount(), greaterThanOrEqualTo(2));
            assertThat("Togetherness too low.", counter.togetherness(), greaterThanOrEqualTo(0.5));
            prev = counter.getCount();
        }

        File saveFile = new File("src/test", "couples.ser");
        testCounter.save(saveFile);
        RoleCoupleCounter loadCounter = RoleCoupleCounter.load(saveFile);
        // Get the role counts from both tables.
        List<CountMap<Role>.Count> oldRoles = testCounter.getRoleCounts();
        List<CountMap<Role>.Count> newRoles = loadCounter.getRoleCounts();
        assertEquals("Role lists do not match after load", oldRoles.size(), newRoles.size());
        for (CountMap<Role>.Count oldCount : oldRoles) {
            int newCount = loadCounter.getCount(oldCount.getKey());
            assertEquals("Counts wrong after load for " + oldCount.getKey().getId(), oldCount.getCount(), newCount);
        }
        // Get the pair counts from both tables.
        List<PairCounter<Role>.Count> oldPairs = testCounter.getPairCounts();
        List<PairCounter<Role>.Count> newPairs = loadCounter.getPairCounts();
        assertEquals("Pair lists do not match after load", oldPairs.size(), newPairs.size());
        for (PairCounter<Role>.Count oldCount : oldPairs) {
            int newCount = loadCounter.getCount(oldCount.getKey1(), oldCount.getKey2());
            assertEquals("Counts wrong after load for " + oldCount.getKey1().getId() + "/" + oldCount.getKey2().getId(),
                    oldCount.getCount(), newCount);
        }
    }




}
