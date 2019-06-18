/*
 * Created on Jul 10, 2007
 *
 */
package org.reactome.r3;

import java.io.IOException;

import org.junit.Test;
import org.reactome.r3.util.FileUtility;

/**
 * This class is used to analyze protein families.
 * @author guanming
 *
 */

public class ProteinFamilyAnalizer {
    
    /**
     * Extract human proteins from the whole Panther UniProt classification file.
     * @throws IOException
     */
    @Test
    public void extractHumanProteinsFromPanther() throws IOException {
        String dirName = "/Users/wgm/Documents/caBIG_R3/datasets/Panther/";
        String fileName = dirName + "PTHR6.1_uniprot8.9";
        String outFileName = dirName + "Human_PTHR6.1_uniprot8.9.txt";
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        FileUtility outFu = new FileUtility();
        outFu.setOutput(outFileName);
        String line = null;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (tokens[1].endsWith("_HUMAN"))
                outFu.printLine(line);
        }
        outFu.close();
        fu.close();;
    }
    
}
