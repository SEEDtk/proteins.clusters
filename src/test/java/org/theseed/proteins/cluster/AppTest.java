package org.theseed.proteins.cluster;

import junit.framework.Test;

import junit.framework.TestCase;
import junit.framework.TestSuite;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
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
     */
    public void testCouples()
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
        Collection<PairCounter<Role>.Count> pairsCounted = testCounter.getPairCounts();
        assertEquals("Wrong number of pairs found.", 4, pairsCounted.size());
        for (PairCounter<Role>.Count counter : pairsCounted) {
            assertThat("Counts out of order.", prev, greaterThanOrEqualTo(counter.getCount()));
            assertEquals("Counts do not match.", testCounter.getCount(counter.getKey1(), counter.getKey2()),
                    counter.getCount());
            prev = counter.getCount();
        }
    }
}
