/*
 * Created on Jul 16, 2015
 *
 */
package org.reactome.cancerhallmark;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import org.junit.Test;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;

/**
 * This class is used to map Reactome pathways to cancer hallmarks.
 * @author gwu
 *
 */
public class PathwayAndHallmarksMapping {
    private final String DIR_NAME = "/Users/gwu/Dropbox/OHSU/CancerHallmarks/";
    
    /**
     * Default constructor.
     */
    public PathwayAndHallmarksMapping() {
    }
    
    @Test
    public void generateHallmarkToPathwaysMap() throws IOException {
        String fileName = DIR_NAME + "Pathways121514_Hallmarks_EnrcihedDriverGenes_071615.txt";
        Map<String, Set<String>> hallmarkToPathways = new HashMap<String, Set<String>>();
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String hallmarks = removeQuotation(tokens[2]);
            String pathway = tokens[1];
            tokens = hallmarks.split(",");
            for (String hallmark : tokens)
                InteractionUtilities.addElementToSet(hallmarkToPathways,
                                                     hallmark, 
                                                     pathway);
        }
        fu.close();
        String outFileName = DIR_NAME + "HallmarkToPathways_071615.txt";
        fu.setOutput(outFileName);
        fu.printLine("Hallmark\tNumber\tPathways");
        for (String hallmark : hallmarkToPathways.keySet()) {
            Set<String> pathways = hallmarkToPathways.get(hallmark);
            fu.printLine(hallmark + "\t" + 
                               pathways.size() + "\t" + 
                               InteractionUtilities.joinStringElements(", ", pathways));
        }
        fu.close();
    }
    
    private String removeQuotation(String token) {
        if (!token.startsWith("\""))
            return token;
        return token.substring(1, token.length() - 1);
    }
    
    /**
     * Merge a mapping file from reactome pathways to cancer hallmarks and another cancer drivers enrichment
     * result file together. The enrichment results were generatd from ReactomeFIViz.
     * @throws IOException
     */
    @Test
    public void mergeTwoFiles() throws IOException {
        String dirName = DIR_NAME;
        // Load pathways to hallmarks information
        String pathwayToHallmarksFile = dirName + "ReactomePathways121514.txt";
        FileUtility fu = new FileUtility();
        fu.setInput(pathwayToHallmarksFile);
        String line = fu.readLine();
        Map<String, String> p2hPathwayToLine = new HashMap<String, String>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String pathway = removeQuotation(tokens[1]);
            p2hPathwayToLine.put(pathway, line);
        }
        fu.close();
        // Load pathways to enrichment results
        String pathwayEnrichmentFile = dirName + "Annotations_CombinedDrivers_CGC_2Natures_071615.txt";
        fu.setInput(pathwayEnrichmentFile);
        line = fu.readLine();
        Map<String, String> pePathwayToLine = new HashMap<String, String>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            pePathwayToLine.put(tokens[0], line);
        }
        fu.close();
        // merge two maps together
        String outFileName = dirName + "Pathways121514_Hallmarks_EnrcihedDriverGenes_071615.txt";
        fu.setOutput(outFileName);
        StringBuilder builder = new StringBuilder();
        builder.append("DB_ID\tPathway\tHallmarks\tRatioOfProteinsInPathway\tNumberOfProteinsInPathway\tProteinsFromDrivers\tP-value\tFDR\tHitGenes");
        fu.printLine(builder.toString());
        builder.setLength(0);
        int counter = 0;
        for (String pathway : p2hPathwayToLine.keySet()) {
            String hallmarkLine = p2hPathwayToLine.get(pathway);
            String[] tokens = hallmarkLine.split("\t");
            // Get the first three tokens
            for (int i = 0; i < 3; i++)
                builder.append(tokens[i]).append("\t");
            String enrichmentLine = pePathwayToLine.get(pathway);
            if (enrichmentLine == null) {
                System.err.println(pathway + " doesn't have an enrichment line!");
                counter ++;
            }
            else {
                tokens = enrichmentLine.split("\t");
                // Get the second to the last tokens
                for (int i = 1; i < tokens.length; i++) {
                    // Remove "<" in FDR for each sorting
                    String token = tokens[i];
                    if (token.startsWith("<"))
                        token = token.substring(1);
                    builder.append(token).append("\t");
                }
            }
            // Remove the last tab
            builder.deleteCharAt(builder.length() - 1);
            fu.printLine(builder.toString());
            builder.setLength(0);
        }
        fu.close();
        System.out.println("Total not matched: " + counter);
    }
    
}
