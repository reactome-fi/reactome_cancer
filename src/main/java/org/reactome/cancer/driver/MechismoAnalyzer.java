/*
 * Created on Apr 10, 2017
 *
 */
package org.reactome.cancer.driver;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.junit.Test;
import org.reactome.r3.util.FileUtility;

/**
 * @author gwu
 *
 */
public class MechismoAnalyzer {
    private String dirName = "results/mechismo_04_17/";
    private String outputFileName = dirName + "mechismo_op.txt";
    private FileUtility fu = new FileUtility();
    
    /**
     * Default constructor.
     */
    public MechismoAnalyzer() {
    }
    
    @Test
    public void checkMechismoScores() throws IOException {
        String fileName = dirName + "SOS1_RAS.txt";
        
        fu.setInput(fileName);
        String line = null;
        Map<String, String> mutationToScore = new HashMap<>();
        Set<String> scores = new HashSet<>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            
//            for (int i = 0; i < tokens.length; i++) {
//                System.out.println(i + "\t" + tokens[i]);
//            }
//            if (true)
//                break;
            
            if (tokens.length < 15) {
//                System.out.println(line);
                continue;
            }
            // Make sure there are AA mutations
            if (tokens[4].length() > 0 && tokens[5].length() > 0) {
                mutationToScore.put(tokens[4] + "\t" + tokens[5], tokens[14]);
                scores.add(tokens[14]);
            }
        }
        System.out.println("Total mutations: " + mutationToScore.size()); // Should be 20 * 19 = 380 
        System.out.println("Total scores: " + new HashSet<>(mutationToScore.values()).size());
        System.out.println("Actual scores: " + scores.size());
        for (String mutation : mutationToScore.keySet())
            System.out.println(mutation + "\t" + mutationToScore.get(mutation));
        System.out.println("\nActual extra scores: ");
        scores.removeAll(mutationToScore.values());
        for (String score : scores)
            System.out.println(score);
    }
    
    @Test
    public void pickRows() throws IOException {
        String[] targetGenes = new String[] {
                "HRAS", "NRAS", "KRAS", "SOS1"
        };
        Set<String> geneSet = new HashSet<>();
        for (String gene : targetGenes)
            geneSet.add(gene);
        
        String fileName = outputFileName;
//        fileName = dirName + "HRAS_mechismo_op.tsv";
        fu.setInput(fileName);
        String line = null;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (tokens.length < 19)
                continue; // Cannot get interactors
//            System.out.println(tokens[0] + "\t" + tokens[18]);
            if (geneSet.contains(tokens[0]) && geneSet.contains(tokens[18]))
                System.out.println(line);
        }
        fu.close();
    }
    
}
