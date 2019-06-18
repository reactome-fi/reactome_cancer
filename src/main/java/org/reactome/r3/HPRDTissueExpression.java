/*
 * Created on Apr 1, 2009
 *
 */
package org.reactome.r3;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.junit.Test;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.R3Constants;
//import org.reactome.weka.FeatureChecker;
//import org.reactome.weka.PositiveChecker;

/**
 * This class is used to handle tissue expression data set downloaded from HPRD.
 * @author wgm
 *
 */
public class HPRDTissueExpression {
    private final String TISSUE_EXP_FILE = R3Constants.HPRD_DIR + "TISSUE_EXPRESSIONS.txt";
    private FileUtility fu = new FileUtility();
    
    public HPRDTissueExpression() {
    }
    
    /**
     * This method is used to test tissue expressions.
     * @throws Exception
     */
    @Test
    public void testCoTissueExpFeature() throws Exception {
//        final Map<String, Set<String>> proteinToTissue = loadProteinToTissues();
//        ReactomeDataAnalyzer analyzer = new ReactomeDataAnalyzer();
//        Set<String> fis = analyzer.loadFIsFromFile();
//        FeatureChecker checker = new FeatureChecker();
//        PositiveChecker posChecker = new PositiveChecker() {
//            public boolean isPositive(String pair) {
//                return isTissueShared(pair, proteinToTissue);
//            }
//        };
//        checker.checkFeatureOddsRatio(posChecker);
    }
    
    private boolean isTissueShared(String pair,
                                   Map<String, Set<String>> proteinToTissues) {
        int index = pair.indexOf(" ");
        String id1 = pair.substring(0, index);
        Set<String> tissues1 = proteinToTissues.get(id1);
        if (tissues1 == null)
            return false;
        String id2 = pair.substring(index + 1);
        Set<String> tissues2 = proteinToTissues.get(id2);
        if (tissues2 == null)
            return false;
        for (String tissue : tissues1) {
            if (tissues2.contains(tissue))
                return true;
        }
        return false;
    }
    
    /**
     * This method is used to load protein to tissues map.
     * @return
     * @throws IOException
     */
    public Map<String, Set<String>> loadProteinToTissues() throws IOException {
        Map<String, String> hprdToUniprot = loadHPRDToUniProtMap();
        fu.setInput(TISSUE_EXP_FILE);
        String line = null;
        Map<String, Set<String>> proteinToTissues = new HashMap<String, Set<String>>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String hprd = tokens[0];
            String uniprot = hprdToUniprot.get(hprd);
            if (uniprot == null)
                continue; // Cannot use this expression
            Set<String> tissues = proteinToTissues.get(uniprot);
            if (tissues == null) {
                tissues = new HashSet<String>();
                proteinToTissues.put(uniprot, tissues);
            }
            tissues.add(tokens[3]);
        }
        return proteinToTissues;
    }
    
    /**
     * Check the downloaded tissue expression file.
     * @throws IOException
     */
    @Test
    public void checkTissueExpressionDataset() throws IOException {
        // Check how many tissues have been used
        fu.setInput(TISSUE_EXP_FILE);
        String line = null;
        Set<String> tissues = new HashSet<String>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            tissues.add(tokens[3]);
        }
        fu.close();
        System.out.println("Total tissues: " + tissues.size());
        for (String tissue : tissues)
            System.out.println(tissue);
    }
    
    private Map<String, String> loadHPRDToUniProtMap() throws IOException {
        // This file provides one-to-one mapping only.
        String mapFile = R3Constants.HPRD_DIR + "HPRD_ID_MAPPINGS.txt";
        Map<String, String> hprdToUniProt = new HashMap<String, String>();
        String line = null;
        fu.setInput(mapFile);
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (tokens[6].equals("-"))
                continue; // Escape it
            hprdToUniProt.put(tokens[0],
                              tokens[6]);
        }
        fu.close();
        return hprdToUniProt;
    }
    
    @Test
    public void testLoadHPRDToUniProtMap() throws IOException {
        Map<String, String> map = loadHPRDToUniProtMap();
        System.out.println("Total map size: " + map.size());
    }
    
}
