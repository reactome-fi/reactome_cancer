/*
 * Created on Aug 4, 2009
 *
 */
package org.reactome.r3;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.junit.Test;
import org.reactome.r3.graph.BreadthFirstSearch;
import org.reactome.r3.graph.GraphAnalyzer;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.MathUtilities;
import org.reactome.r3.util.R3Constants;

/**
 * This method is used to analyze OMIM related data. Some of code is copied from Xin.
 * @author wgm
 *
 */
public class OMIMAnalyzer {
    private FileUtility fu = new FileUtility();
    
    public OMIMAnalyzer() {
    }
    
    /**
     * This method is used to check all FIs (pathways + predicted) for pathway genes only
     * to minimize the false negative problem.
     */
    @Test
    public void checkOMIMAssociationForAllFIsInPathways() throws Exception {
        String intFileName = R3Constants.GENE_FI_FILE_NAME;
        //String pathwayFileName = R3Constants.RESULT_DIR + "FIsInGene_Pathway_041709.txt";
        String pathwayFileName = R3Constants.RESULT_DIR + "FIsInGene_Predicted_041709.txt";
        Set<String> pathwayFIs = fu.loadInteractions(pathwayFileName);
        Set<String> pathwayFIGenes = InteractionUtilities.grepIDsFromInteractions(pathwayFIs);
        Set<String> allFIs = fu.loadInteractions(intFileName);
        Set<String> allFIGenes = InteractionUtilities.grepIDsFromInteractions(allFIs);
        Set<String> allFIsInPathwayGenes = InteractionUtilities.getFIs(pathwayFIGenes, allFIs);
        System.out.println("Total pathway FIs: " + pathwayFIs.size());
        System.out.println("Total pathway genes: " + pathwayFIGenes.size());
        System.out.println("Total FIs in pathway genes: " + allFIsInPathwayGenes.size());
        BreadthFirstSearch bfs = new BreadthFirstSearch();
        Map<String, Set<String>> geneToPartners = bfs.generateIdToPartnersMap(allFIsInPathwayGenes);
        final Map<String, Integer> geneToDegree = new HashMap<String, Integer>();
        for (String gene : geneToPartners.keySet()) {
            Set<String> partners = geneToPartners.get(gene);
            geneToDegree.put(gene, partners.size());
        }
        List<String> geneList = new ArrayList<String>(geneToDegree.keySet());
        // Sort proteins based on connection degrees
        Collections.sort(geneList, new Comparator<String>() {
           public int compare(String protein1, String protein2) { 
               Integer degree1 = geneToDegree.get(protein1);
               Integer degree2 = geneToDegree.get(protein2);
               return degree2 - degree1;
           }
        });
        List<String> hubs = new ArrayList<String>();
        List<String> nonHubs = new ArrayList<String>();
        Set<String> omimGenes = loadOMIMGenes();
        for (int i = 1; i < 11; i++) {
            double percentile = 0.01 * i;
            System.out.println("Hub percentile: " + percentile);
            getHubNonHubNodes(hubs, 
                              nonHubs,
                              geneList,
                              geneToDegree, 
                              percentile);
            reportOMIMResults(allFIGenes,
                              hubs, 
                              nonHubs,
                              omimGenes);
        }
    }
    
    /**
     * This method is used to calculate percentage of omim associations for hub and 
     * non hub proteins or genes.
     * @throws Exception
     */
    @Test
    public void checkOMIMAssociation() throws Exception {
        String intFileNames[] = new String[] {
//                R3Constants.GENE_FI_FILE_NAME,
//                R3Constants.RESULT_DIR + "FIsInGene_Pathway_041709.txt",
                R3Constants.RESULT_DIR + "FIsInGene_Predicted_041709.txt",
//                "results/v2/FI73InGene_061008.txt",
//                "results/v2/FI73InGene_Pathway_061008.txt",
//                "results/v2/FI73InGene_Predicated_061008.txt"
        };
        GraphAnalyzer graphAnalyzer = new GraphAnalyzer();
        for (String intFileName : intFileNames) {
            System.out.println("interaction file name: " + intFileName);
            Set<String> fis = fu.loadInteractions(intFileName);
            Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
            List<String> hubs = new ArrayList<String>();
            List<String> nonHubs = new ArrayList<String>();
            final Map<String, Integer> nodeToDegree = graphAnalyzer.generateProteinToDegree(intFileName);
            List<String> nodeList = new ArrayList<String>(nodeToDegree.keySet());
            // Sort proteins based on connection degrees
            Collections.sort(nodeList, new Comparator<String>() {
               public int compare(String protein1, String protein2) { 
                   Integer degree1 = nodeToDegree.get(protein1);
                   Integer degree2 = nodeToDegree.get(protein2);
                   return degree2 - degree1;
               }
            });
            Set<String> omimGenes = loadOMIMGenes();
            System.out.println("Total OMIM genes: " + omimGenes.size());
            for (int i = 1; i < 11; i++) {
                double percentage = 0.01 * i;
                System.out.println("Percentile: " + percentage);
                getHubNonHubNodes(hubs, 
                                  nonHubs,
                                  nodeList, 
                                  nodeToDegree,
                                  percentage);
                reportOMIMResults(fiGenes,
                                  hubs,
                                  nonHubs,
                                  omimGenes);
            }
        }
    }

    private void reportOMIMResults(Set<String> fiGenes,
                                   List<String> hubs,
                                   List<String> nonHubs,
                                   Set<String> omimGenes) {
        int hubSize = hubs.size();
        hubs.retainAll(omimGenes);
        int hubOmimSize = hubs.size();
        int nonHubSize = nonHubs.size();
        nonHubs.retainAll(omimGenes);
        int nonHubOmimSize = nonHubs.size(); 
        double zvalue = MathUtilities.calculateZValue(hubOmimSize, hubSize,
                                                      nonHubOmimSize, nonHubSize);
        double pvalue = MathUtilities.calTwoTailStandardNormalPvalue(zvalue);
        double hubRatio = (double) hubOmimSize / hubSize;
        double nonHubRatio = (double) nonHubOmimSize / nonHubSize;
        System.out.println("Hubs: " + hubRatio);
        System.out.println("Non-Hubs: " + nonHubRatio);
        System.out.println("pvalue: " + pvalue);
        fiGenes.retainAll(omimGenes);
        System.out.println("FI genes with OMIM associations: " + fiGenes.size());
        System.out.println();
    }
    
    private void getHubNonHubNodes(List<String> hubs,
                                   List<String> nonHubs,
                                   List<String> nodeList,
                                   Map<String, Integer> nodeToDegree,
                                   double percentage) {
        hubs.clear();
        nonHubs.clear();
        int hubSize = (int) (nodeToDegree.size() * percentage);
        for (int i = 0; i < nodeList.size(); i++) {
            if (i < hubSize)
                hubs.add(nodeList.get(i));
            else
                nonHubs.add(nodeList.get(i));
        }
    }
    
    /**
     * Use this method to load OMIM genes that are labeled as phenotype.
     * @return
     * @throws IOException
     */
    private Set<String> loadOMIMGenes() throws IOException {
        // Need to map gene names to gene ids
        Map<String, String> geneIdToName = new HashMap<String, String>();
        String fileName = R3Constants.NCBI_DIR + "Homo_sapiens.gene_info";
        fu.setInput(fileName);
        String line = null;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            geneIdToName.put(tokens[1], tokens[2]);
        }
        fu.close();
        fileName = R3Constants.NCBI_DIR + "mim2gene.txt";
        fu.setInput(fileName);
        line = fu.readLine();
        Set<String> omimGenes = new HashSet<String>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (tokens[2].equals("phenotype")) {
                String geneName = geneIdToName.get(tokens[1]);
                if (geneName != null)
                    omimGenes.add(geneName);
            }
        }
        fu.close();
        return omimGenes;
    }
    
}
