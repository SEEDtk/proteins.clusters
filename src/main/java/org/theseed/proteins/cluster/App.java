package org.theseed.proteins.cluster;

/**
 * Process genomes to find functionally-coupled roles.
 *
 */
public class App
{
    public static void main( String[] args )
    {
        RoleCouplingProcessor runObject = new RoleCouplingProcessor();
        boolean ok = runObject.parseCommand(args);
        if (ok) {
            runObject.run();
        }
    }
}
