/*
 * Created on Mar 8, 2012
 *
 */
package org.reactome.r3;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import org.junit.Test;
import org.reactome.data.UniProtAnalyzer;
import org.reactome.fi.util.FeatureChecker;
import org.reactome.fi.util.PositiveChecker;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.FisherExact;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.R3Constants;

/**
 * This class is used to process radiation hybrid generated genetic interaction
 * file downloaded from Lin et al, Genome Research 20: 1122-1132 (2010).
 * @author gwu
 *
 */
public class RadiationHybridAnalyzer {
    private FileUtility fu = new FileUtility();
    private Set<String> rhInteractions = null;
    private final double FDR_THRESHOLD = 0.001;
    private final String RH_FILE_NAME = R3Constants.RH_DIR_NAME + "fully_combined_RH_network.txt";
    
    public RadiationHybridAnalyzer() {
    }
    
    @Test
    public void checkOverlappingViaGenes() throws Exception {
        Set<String> rhInteractions = loadRHInteractionsInGenes();
        Set<String> reactomeFIs = fu.loadInteractions(R3Constants.RESULT_DIR + "FIs_Reactome.txt");
        Map<String, Set<String>> geneToAcces = new UniProtAnalyzer().generateGeneNameToUniAccessMap(true);
        Map<String, Set<String>> accToGenes = InteractionUtilities.switchKeyValues(geneToAcces);
        Set<String> reactomeGeneFIs = new HashSet<String>();
        for (String fi : reactomeFIs) {
            String[] tokens = fi.split("\t");
            Set<String> genes1 = accToGenes.get(tokens[0]);
            Set<String> genes2 = accToGenes.get(tokens[1]);
            if (genes1 == null || genes2 == null)
                continue;
            for (String gene1 : genes1) {
                for (String gene2 : genes2) {
                    if (gene1.equals(gene2))
                        continue;
                    if (gene1.compareTo(gene2) < 0) 
                        reactomeGeneFIs.add(gene1 + "\t" + gene2);
                    else
                        reactomeGeneFIs.add(gene2 + "\t" + gene1);
                }
            }
        }
        // Do some filtering
        Set<String> rhGenes = InteractionUtilities.grepIDsFromInteractions(rhInteractions);
        Set<String> reactomeGenes = InteractionUtilities.grepIDsFromInteractions(reactomeGeneFIs);
        Set<String> sharedGenes = InteractionUtilities.getShared(rhGenes, reactomeGenes);
        rhInteractions = InteractionUtilities.getFIs(sharedGenes, rhInteractions);
        reactomeGeneFIs = InteractionUtilities.getFIs(sharedGenes, reactomeGeneFIs);
        Set<String> sharedGeneFIs = InteractionUtilities.getShared(rhInteractions, reactomeGeneFIs);
        int total = sharedGenes.size() * (sharedGenes.size() - 1) / 2;
        int f11 = sharedGeneFIs.size();
        int f12 = reactomeGeneFIs.size() - sharedGeneFIs.size();
        int f21 = rhInteractions.size() - sharedGeneFIs.size();
        int f22 = total - f11 - f12 - f21;
        System.out.println("Matrix: " + f11 + ", " + f12 + ", " + f21 + ", " + f22);
    }
    
    /**
     * Check overlapping between Reactome FIs and RH networks. This test is differnt from the 
     * method reported in the original paper: genes used in test have been filtered so that
     * both genes in the both network are shared.
     * @throws Exception
     */
    @Test
    public void checkOverlapping() throws Exception {
        String[] srcFileNames = new String[] {
//                R3Constants.RESULT_DIR + "FIs_Reactome.txt",
                R3Constants.RESULT_DIR + "FIs_KEGG.txt",
//                R3Constants.RESULT_DIR + "FIs_Pathway Interaction Database.txt",
//                R3Constants.RESULT_DIR + "FIs_BioCarta - Imported by PID.txt",
//                R3Constants.RESULT_DIR + "FIs_TRED.txt",
//                R3Constants.RESULT_DIR + "FIs_pantherdb.txt",
//                R3Constants.IREFINDEX_HUMAN_PPI_FILE
                R3Constants.LEE_GENE_EXP_FILE,
                R3Constants.PRIETO_GENE_EXP_FILE
        };
        for (String srcFileName : srcFileNames) {
            int index = srcFileName.lastIndexOf("/");
            System.out.println("\nSource file name: " + srcFileName.substring(index + 1));
            //        Set<String> reactomeFIs = new ReactomeAnalyzer().loadFIsFromFile();
            Set<String> reactomeFIs = fu.loadInteractions(srcFileName);
            // Just for test
//            reactomeFIs = new FeatureChecker().generateRandomPairs(reactomeFIs);
            //        reactomeFIs = replaceSpaceByTab(reactomeFIs);
            System.out.println("Total Reactome FIs: " + reactomeFIs.size());
            Set<String> reactomeProteins = InteractionUtilities.grepIDsFromInteractions(reactomeFIs);
            loadRHInteractions();
            Set<String> rhProteins = InteractionUtilities.grepIDsFromInteractions(rhInteractions);
            System.out.println("Total RH networks: " + rhInteractions.size());
            // Get shared proteins
            Set<String> sharedProteins = InteractionUtilities.getShared(reactomeProteins, rhProteins);
            // Need to prune the RH interactions so that only proteins have been in
            // the Reactome FIs are used as reported in the original paper.
            for (Iterator<String> it = rhInteractions.iterator(); it.hasNext();) {
                String rhFI = it.next();
                String[] tokens = rhFI.split("\t");
                if (sharedProteins.contains(tokens[0]) && sharedProteins.contains(tokens[1]))
                    continue;
                it.remove();
            }
            System.out.println("Total RH after filtering: " + rhInteractions.size());
            for (Iterator<String> it = reactomeFIs.iterator(); it.hasNext();) {
                String fi = it.next();
                String[] tokens = fi.split("\t");
                if (sharedProteins.contains(tokens[0]) && sharedProteins.contains(tokens[1]))
                    continue;
                it.remove();
            }
            //        Set<String> rhProteins = InteractionUtilities.grepIDsFromInteractions(rhInteractions);
            //        System.out.println("Total RH proteins: " + rhProteins.size());
            Set<String> shared = InteractionUtilities.getShared(reactomeFIs, rhInteractions);
            System.out.println("Total shared: " + shared.size());
            int totalPairs = sharedProteins.size() * (sharedProteins.size() - 1) / 2;
            int f11 = shared.size();
            int f12 = reactomeFIs.size() - f11;
            int f21 = rhInteractions.size() - f11;
            int f22 = totalPairs - f11 - f12 - f21;
            FisherExact fisherExact = new FisherExact(f22);
            System.out.println("Matrix: " + f11 + ", " + f12 + ", " + f21 + ", " + f22);
//            double pvalue = fisherExact.getP(f11, f12, f21, f22);
            //        System.out.println("Pvalue from Fisher exact: " + pvalue);
            //        pvalue = MathUtilities.calculateHypergeometricPValue(totalPairs, 
            //                                                             reactomeFIs.size(),
            //                                                             rhInteractions.size(), 
            //                                                             shared.size());
            //        System.out.println("Pvalue from hypergeometric: " + pvalue);
        }
    }

    private Set<String> replaceSpaceByTab(Set<String> reactomeFIs) {
        Set<String> copy = new HashSet<String>();
        for (String fi : reactomeFIs) {
            fi = fi.replace(" ", "\t");
            copy.add(fi);
        }
        reactomeFIs = copy;
        return reactomeFIs;
    }
    
    /**
     * Check odds ratio for this feature if it is used in NBC.
     * @throws Exception 
     */
    @Test
    public void testRHFeature() throws Exception {
        FeatureChecker checker = new FeatureChecker();
        PositiveChecker posChecker = new PositiveChecker() {
            public boolean isPositive(String pair) {
                try {
                    return checkIfInteracting(pair);
                }
                catch(IOException e) {
                    e.printStackTrace();
                }
                return false;
            }
        };
        checker.checkFeatureOddsRatio(posChecker);
    }
    
    private boolean checkIfInteracting(String pair) throws IOException {
        if (rhInteractions == null)
            loadRHInteractions();
        return rhInteractions.contains(pair);
    }
    
    private Set<String> loadRHInteractionsInGenes() throws IOException {
        Map<Double, Double> fdrToPvalue = loadFDRToPValue();
        Double pvalue = fdrToPvalue.get(FDR_THRESHOLD);
        fu.setInput(RH_FILE_NAME);
        String line = fu.readLine();
        Set<String> geneInteractions = new HashSet<String>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (pvalue < new Double(tokens[2]))
                continue; // This is inclusive
            String gene1 = tokens[0];
            String gene2 = tokens[1];
            if (gene1.equals(gene2))
                continue;
            if (gene1.compareTo(gene2) < 0)
                geneInteractions.add(gene1 + "\t" + gene2);
            else
                geneInteractions.add(gene2 + "\t" + gene1);
        }
        fu.close();
        System.out.println("Total RH interactions in UniProt accessions: " + geneInteractions.size());
        return geneInteractions;
    }
    
    private void loadRHInteractions() throws IOException {
        rhInteractions = new HashSet<String>();
        Map<Double, Double> fdrToPvalue = loadFDRToPValue();
        Double pvalue = fdrToPvalue.get(FDR_THRESHOLD);
        UniProtAnalyzer uniProtAnalyzer = new UniProtAnalyzer();
        Map<String, Set<String>> geneNameToUniAcces = uniProtAnalyzer.generateGeneNameToUniAccessMap(true);
        fu.setInput(RH_FILE_NAME);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (pvalue < new Double(tokens[2]))
                continue; // This is inclusive
            Set<String> uniAccess1 = geneNameToUniAcces.get(tokens[0]);
            Set<String> uniAccess2 = geneNameToUniAcces.get(tokens[1]);
            if (uniAccess1 == null || uniAccess2 == null)
                continue;
            for (String acc1 : uniAccess1) {
                for (String acc2 : uniAccess2) {
                    if (acc1.equals(acc2))
                        continue;
                    if (acc1.compareTo(acc2) < 0)
                        rhInteractions.add(acc1 + "\t" + acc2);
                    else
                        rhInteractions.add(acc2 + "\t" + acc1);
                }
            }
        }
        fu.close();
        System.out.println("Total RH interactions in UniProt accessions: " + rhInteractions.size());
//        int count = 0;
//        for (String i : rhInteractions) {
//            System.out.println(i);
//            count ++;
//            if (count == 100)
//                break;
//        }
    }
    
    @Test
    public void checkNumbers() throws IOException {
        String fileName;
        String line;
        Map<Double, Double> fdrToPvalue = loadFDRToPValue();
        fileName = R3Constants.RH_DIR_NAME + "fully_combined_RH_network.txt";
        for (Double fdr : fdrToPvalue.keySet()) {
            Double pvalue = fdrToPvalue.get(fdr);
            System.out.println(fdr + " -> " + pvalue);
            fu.setInput(fileName);
            line = fu.readLine();
            int count = 0;
            while ((line = fu.readLine()) != null) {
                String[] tokens = line.split("\t");
                if (pvalue >= new Double(tokens[2]))
                    count ++;
            }
            System.out.println("Total interactions: " + count);
            System.out.println();
            fu.close();
        }
    }

    private Map<Double, Double> loadFDRToPValue() throws IOException {
        Map<Double, Double> fdrToPvalue = new HashMap<Double, Double>();
        String fileName = R3Constants.RH_DIR_NAME + "fdr_thresholds_for_fully_combined_RH_network.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            if (line.trim().length() == 0)
                break;
            String[] tokens = line.split("\t");
            fdrToPvalue.put(new Double(tokens[0]),
                            new Double(tokens[1]));
        }
        fu.close();
        return fdrToPvalue;
    }
    
}
