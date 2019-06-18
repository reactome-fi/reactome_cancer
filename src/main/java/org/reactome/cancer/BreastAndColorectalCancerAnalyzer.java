/*
 * Created on Mar 1, 2010
 *
 */
package org.reactome.cancer;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.junit.Test;
import org.reactome.annotate.AnnotationType;
import org.reactome.annotate.GeneSetAnnotation;
import org.reactome.annotate.PathwayBasedAnnotator;
import org.reactome.r3.graph.SpectralPartitionNetworkCluster;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.R3Constants;

public class BreastAndColorectalCancerAnalyzer extends CancerResequenceDataSetAnalyzer {
    private final String WORKING_DIR = CancerConstants.BREAST_AND_COLORECTAL_DIR;
    private String cancerType = "Breast";
//    private String cancerType = "Colorectal";
    
    public BreastAndColorectalCancerAnalyzer() {
    }
    
    public Map<String, Set<String>> loadSampleToMutatedGenes(String cancerType,
                                                             boolean useValidaton) throws IOException {
        Map<String, Set<String>> sampleToGenes = new HashMap<String, Set<String>>();
        String fileName = WORKING_DIR + "Mutations.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        Set<String> mutationTypes = new HashSet<String>();
        while ((line = fu.readLine()) != null) {
            if (line.trim().length() == 0)
                break;
            String[] tokens = line.split("\t");
            String type = tokens[3];
            if (!type.equals(cancerType))
                continue;
            String mutationType = tokens[8];
            if (mutationType.equals("UTR") || mutationType.equals("Synonymous"))
                continue;
            if (!useValidaton) {
                if (!tokens[4].equals("Discovery"))
                    continue;
            }
//            if (!tokens[4].equals("Validation"))
//                continue;
//            if (!tokens[4].equals("Discovery"))
//                continue;
            mutationTypes.add(mutationType);
            String gene = tokens[0];
            String tumor = tokens[2];
            InteractionUtilities.addElementToSet(sampleToGenes, tumor, gene);
        }
        fu.close();
        //System.out.println("Total mutation types: " + mutationTypes);
        return sampleToGenes;
    }
    
    @Test
    public void annotateClusters() throws Exception {
//        Map<String, Set<String>> sampleToGenes = loadSampleToSequenceAlteredGenes();
//        Set<String> allGenes = InteractionUtilities.grepAllGenes(sampleToGenes);
//        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
//        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
//        allGenes.retainAll(fiGenes);
//        System.out.println("All genes in FIs for randomization: " + allGenes.size());
        String clusterFileName = WORKING_DIR + cancerType + "CancerClusters030210.txt";
        NetworkClusterAnalyzer clusterAnalyzer = new NetworkClusterAnalyzer();
        clusterAnalyzer.annotateNetworkClusters(clusterFileName, 4);
//        clusterAnalyzer.annotateNetworkClusters(clusterFileName, 
//                                                4,
//                                                allGenes);
    }
    
    @Test
    public void annotateFIsInClusters() throws Exception {
        String fileName = WORKING_DIR + "FIsInClusterListFor" + cancerType + "Cancer030210.txt";
        Set<String> fis = fu.loadInteractions(fileName);
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        PathwayBasedAnnotator annotator = new PathwayBasedAnnotator();
        annotator.setPValueThreshold(0.01);
        annotator.setFDRThreshold(0.25);
        annotator.annotateGenesWithFDR(fiGenes, AnnotationType.Pathway);
    }
    
    @Test
    public void annotateHierarchicalClusters() throws Exception {
        String clusterFileName = WORKING_DIR + "ClusterListFor" + cancerType + "Cancer030210.txt";
        System.out.println("Cluster file: " + clusterFileName);
        Map<String, List<String>> headerToIds = loadClusters(clusterFileName);
        for (String tmp : headerToIds.keySet()) {
            System.out.println(tmp);
            List<String> names = headerToIds.get(tmp);
            annotateNames(names);
            System.out.println();
        }
    }
    
    /**
     * This method is used to check sample distribution in network clusters.
     * @throws Exception
     */
    @Test
    public void checkSamplesInNetworkClusters() throws Exception {
//        String clusterFileName = WORKING_DIR + cancerType + "CancerClusters030210.txt";
//        String clusterFileName = WORKING_DIR + cancerType + "CancerSpectralClusters032610.txt";
        String clusterFileName = WORKING_DIR + cancerType + "CancerClusters042810.txt";
        Map<String, Set<String>> sampleToAlteredGenes = loadSampleToMutatedGenes(cancerType,
                                                                                 true);
        //Map<String, Set<String>> sampleToAlteredGenes = loadSampleToMutations();
        System.out.println("Total samples: " + sampleToAlteredGenes.size());
        NetworkClusterAnalyzer clusterAnalyzer = new NetworkClusterAnalyzer();
        List<Integer> targetClusters = new ArrayList<Integer>();
        targetClusters.add(0);
        targetClusters.add(1);
        //targetClusters.add(2);
        //targetClusters.add(3);
        targetClusters.add(4);
        double percentage = clusterAnalyzer.checkSamplesInClusters(sampleToAlteredGenes, 
                                                                   clusterFileName,
                                                                   targetClusters);
        System.out.println("Samples in clusters 0, 1, 4: " + percentage);
        // Check for clusters 0, 1:
        targetClusters.clear();
        targetClusters.add(0);
        targetClusters.add(1);
        targetClusters.add(2);
        targetClusters.add(3);
        percentage = clusterAnalyzer.checkSamplesInClusters(sampleToAlteredGenes, 
                                                            clusterFileName,
                                                            targetClusters);
        System.out.println("Samples in clusters 0, 1, 2, 3: " + percentage);
        // Check for clusters 0, 2:
        targetClusters.clear();
        targetClusters.add(0);
        targetClusters.add(1);
        targetClusters.add(2);
        targetClusters.add(4);
        percentage = clusterAnalyzer.checkSamplesInClusters(sampleToAlteredGenes, 
                                                            clusterFileName,
                                                            targetClusters);
        System.out.println("Samples in clusters 0, 1, 2, 4: " + percentage);
        // Check for clusters 0, 2:
        targetClusters.clear();
        targetClusters.add(0);
        targetClusters.add(1);
        targetClusters.add(2);
        percentage = clusterAnalyzer.checkSamplesInClusters(sampleToAlteredGenes, 
                                                            clusterFileName,
                                                            targetClusters);
        System.out.println("Samples in clusters 0, 1, 2: " + percentage);
//        // Check for clusters 0, 2:
//        targetClusters.clear();
//        targetClusters.add(0);
//        targetClusters.add(4);
//        percentage = clusterAnalyzer.checkSamplesInClusters(sampleToAlteredGenes, 
//                                                            clusterFileName,
//                                                            targetClusters);
//        System.out.println("Samples in clusters 0, 4: " + percentage);
    }
    
    @Test
    public void permutationTestForSamplesInNetworkClustersOnEdgeBetweenness() throws Exception {
        String clusterFileName = WORKING_DIR + cancerType + "CancerClusters030210.txt";
        Map<String, Set<String>> sampleToAlteredGenes = loadSampleToMutatedGenes(cancerType,
                                                                                 false);

        System.out.println("Total samples: " + sampleToAlteredGenes.size());
        NetworkClusterAnalyzer clusterAnalyzer = new NetworkClusterAnalyzer();
        List<Set<String>> clusters = clusterAnalyzer.loadNetworkClusters(clusterFileName);
        List<Integer> targetClusters = new ArrayList<Integer>();
        targetClusters.add(0);
        targetClusters.add(1);
        //targetClusters.add(2);
        //targetClusters.add(3);
        targetClusters.add(4);
        //double target = 0.5454545454545454;
        //double target = 0.8181818181818182;
        List<Integer> targetClusterSizes = new ArrayList<Integer>();
        for (Integer index : targetClusters) {
            Set<String> cluster = clusters.get(index);
            targetClusterSizes.add(cluster.size());
            System.out.println(index + ": " + cluster.size());
        }
        clusterAnalyzer.permutationTestForSamplesInClustersOnBetweenness(sampleToAlteredGenes, 
                                                                         targetClusters,
                                                                         targetClusterSizes);
    }
    
    @Test
    public void permutationTestForSamplesInNetworkClusters() throws Exception {
        //String clusterFileName = WORKING_DIR + cancerType + "CancerClusters030210.txt";
        //String clusterFileName = WORKING_DIR + cancerType + "CancerClusters042810.txt";
        //System.out.println("Loading cluster file: " + clusterFileName);
        List<Set<String>> clusters = clusterAlteredGenesViaSpectral();
        
        Map<String, Set<String>> sampleToAlteredGenes = loadSampleToMutatedGenes(cancerType,
                                                                                 false);

        System.out.println("Total samples: " + sampleToAlteredGenes.size());
        NetworkClusterAnalyzer clusterAnalyzer = new NetworkClusterAnalyzer();
        List<Integer> targetClusters = new ArrayList<Integer>();
        targetClusters.add(0);
        targetClusters.add(1);
        //targetClusters.add(4);
//        double target = 0.8181818181818182;
        double target = 1.0;
        System.out.println("\nPermutation for samples in clusters " + targetClusters +":");
        clusterAnalyzer.permutationTestForSamplesInNetworkClusters(sampleToAlteredGenes,
                                                                   clusters,//clusterFileName,
                                                                   targetClusters, 
                                                                   target,
                                                                   false);
        
        targetClusters.clear();
        targetClusters.add(0);
        targetClusters.add(1);
        targetClusters.add(2);
        target = 0.8181818181818182;
//        double target = 1.0;
        System.out.println("\nPermutation for samples in clusters " + targetClusters +":");
        clusterAnalyzer.permutationTestForSamplesInNetworkClusters(sampleToAlteredGenes,
                                                                   clusters, //clusterFileName,
                                                                   targetClusters, 
                                                                   target,
                                                                   false);
    }
    
    @Test
    public void compareTwoNetworkClusters() throws IOException {
        NetworkClusterAnalyzer clusterAnalyzer = new NetworkClusterAnalyzer();
        String fileName = WORKING_DIR + cancerType + "CancerClusters042810.txt";
        //String fileName = TCGA_GBM_DIR + "ClustersInTCGAAlteredGenes042810.txt";
        List<Set<String>> scienceClusters = clusterAnalyzer.loadNetworkClusters(fileName);
        //fileName = TCGA_GBM_DIR + "ClustersInTCGAAlteredGenes071609.txt";
        fileName = WORKING_DIR + cancerType + "CancerClusters030210.txt";
        List<Set<String>> tcgaClusters = clusterAnalyzer.loadNetworkClusters(fileName);
        for (Set<String> cluster : tcgaClusters)
            System.out.println(cluster);
        // Check sharing
        System.out.println("TCGA_Cluster\tSize\tScience_Cluster\tSize\tShared\tShared_Percentage");
        // Just want to see the first five clusters
        for (int i = 0; i < tcgaClusters.size(); i++) {
            Set<String> cluster1 = tcgaClusters.get(i);
            for (int j = 0; j < scienceClusters.size(); j++) {
                Set<String> cluster2 = scienceClusters.get(j);
                int min = Math.min(cluster1.size(), cluster2.size());
                Set<String> copy = new HashSet<String>(cluster2);
                copy.retainAll(cluster1);
                double shared = copy.size() / (double) min;
                if (copy.size() == 0)
                    continue;
                System.out.println(i + "\t" + cluster1.size() + "\t" +
                                   j + "\t" + cluster2.size() + "\t" +
                                   copy.size() + "\t" + shared);
            }
        }
    }
    
    private List<Set<String>> clusterAlteredGenesViaSpectral() throws Exception {
        Map<String, Set<String>> sampleToGenes = loadSampleToMutatedGenes(cancerType, 
                                                                          false);
        System.out.println("Total samples in " + cancerType + " cancer: " + sampleToGenes.size());
        Set<String> genes = InteractionUtilities.grepAllGenes(sampleToGenes);
        System.out.println("Total genes: " + genes.size());
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> genesInFIs = InteractionUtilities.grepIDsFromInteractions(fis);
        genes.retainAll(genesInFIs);
        System.out.println("Total genes in the FI network: " + genes.size());
        Set<String> fisInGenes = InteractionUtilities.getFIs(genes, fis);
        Set<String> genesInFI = InteractionUtilities.grepIDsFromInteractions(fisInGenes);
        System.out.println("Gene in FI: " + genesInFI.size());
        NetworkClusterAnalyzer clusterAnalyzer = new NetworkClusterAnalyzer();
        SpectralPartitionNetworkCluster spectralClustering = new SpectralPartitionNetworkCluster();
        long time11 = System.currentTimeMillis();
        List<Set<String>> clusters = spectralClustering.cluster(fisInGenes);
        long time21 = System.currentTimeMillis();
        return clusters;
    }
    
    public void test() throws Exception {
        // Get a list of genes
        Set<String> geneSet = new HashSet<String>();
        // Spectral partitioning based network clustering
        SpectralPartitionNetworkCluster spectralClustering = new SpectralPartitionNetworkCluster();
        List<Set<String>> clusters = spectralClustering.clusterGenes(geneSet);
        // If you want to output the above clusters
        NetworkClusterAnalyzer clusterAnalyzer = new NetworkClusterAnalyzer();
        String outputFileName = ""; // Give it a name
        clusterAnalyzer.outputNetworkClusters(clusters, outputFileName);
        // Annotate network clusters
        // In order to use the following class, you have to set up an attribute for
        // AnnotationHelper as follows:
        // private String proteinNameToPathwayFile = R3Constants.RESULT_DIR + "ProteinNameToTopics072512.txt";
        PathwayBasedAnnotator annotator = new PathwayBasedAnnotator();
        // Go through each module and annotate it
        // You may do some filtering here (e.g. don't annotate small modules)
        for (Set<String> cluster : clusters) {
            if (cluster.size() < 3)
                continue; // Filter out modules having 2 genes only
            // Choose AnnotationType for different annotations (e.g. Pathway, GO BP, etc).
            List<GeneSetAnnotation> annotations = annotator.annotateGenesWithFDR(cluster, AnnotationType.Pathway);
            // Do whatever you want with annotations
        }
    }
    
    @Test
    public void clusterAlteredGenes() throws Exception {
        Map<String, Set<String>> sampleToGenes = loadSampleToMutatedGenes(cancerType, 
                                                                          false);
        System.out.println("Total samples in " + cancerType + " cancer: " + sampleToGenes.size());
        Set<String> genes = InteractionUtilities.grepAllGenes(sampleToGenes);
        System.out.println("Total genes: " + genes.size());
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> genesInFIs = InteractionUtilities.grepIDsFromInteractions(fis);
        genes.retainAll(genesInFIs);
        System.out.println("Total genes in the FI network: " + genes.size());
        Set<String> fisInGenes = InteractionUtilities.getFIs(genes, fis);
        Set<String> genesInFI = InteractionUtilities.grepIDsFromInteractions(fisInGenes);
        System.out.println("Gene in FI: " + genesInFI.size());
//        Set<String> fisInGenes = InteractionUtilities.getFIs(genes, fis);
//        String fiFileName = WORKING_DIR + cancerType + "FIsInMutatedGenes031510.txt";
//        fu.saveInteractions(fisInGenes,
//                            fiFileName);
        NetworkClusterAnalyzer clusterAnalyzer = new NetworkClusterAnalyzer();
//        long time1 = System.currentTimeMillis();
//        List<Set<String>> clusters = clusterAnalyzer.cluster(genes); 
//                                                             //0.105);
//        long time2 = System.currentTimeMillis();
//        System.out.println("Time for edge betweenness: " + (time2 - time1));
//        double modularity = clusterAnalyzer.calculateModularity(clusters, fisInGenes);
//        System.out.println("Modularity from edge betweenness: " + modularity);
        
        SpectralPartitionNetworkCluster spectralClustering = new SpectralPartitionNetworkCluster();
        //double modularity = spectralClustering.calculateModualarity(clusters, fisInGenes);
        //System.out.println("Modularity using new: " + modularity);
        long time11 = System.currentTimeMillis();
        List<Set<String>> clusters = spectralClustering.cluster(fisInGenes);
        long time21 = System.currentTimeMillis();
        System.out.println("\nTime for spectral: " + (time21 - time11));
        double modularity = clusterAnalyzer.calculateModularity(clusters, fisInGenes);
        System.out.println("Modularity from spectral: " + modularity);
        modularity = spectralClustering.calculateModualarity(clusters, fisInGenes);
        System.out.println("Modularity using new: " + modularity);
        
        //String clusterFileName = WORKING_DIR + cancerType + "CancerClusters030210.txt";
        //String clusterFileName = WORKING_DIR + cancerType + "CancerSpectralClusters032610.txt";
        String clusterFileName = WORKING_DIR + cancerType + "CancerClusters042810.txt";
        clusterAnalyzer.outputNetworkClusters(clusters,
                                              clusterFileName);
        for (int i = 0; i < clusters.size(); i++)
            System.out.println(i + "\t" + clusters.get(i).size());
    }
    
    @Test
    public void checkCoverage() throws IOException {
        Map<String, Set<String>> sampleToGenes = loadSampleToMutatedGenes(cancerType, 
                                                                          true);
        System.out.println("Total samples in Breast cancer: " + sampleToGenes.size());
        Set<String> genes = InteractionUtilities.grepAllGenes(sampleToGenes);
        System.out.println("Total genes: " + genes.size());
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_FILE_NAME);
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        int total = genes.size();
        genes.retainAll(fiGenes);
        System.out.println("Genes in FI: " + genes.size() + 
                           " (" +  (double) genes.size() / total + ")");
    }
    
    @Test
    public void testLoad() throws IOException {
        Map<String, Set<String>> sampleToGenes = loadSampleToMutatedGenes("Breast");
        System.out.println("Total samples in Breast cancer: " + sampleToGenes.size());
        List<String> samples = new ArrayList<String>(sampleToGenes.keySet());
        Collections.sort(samples);
        for (String sample : samples)
            System.out.println(sample);
        Set<String> genes = InteractionUtilities.grepAllGenes(sampleToGenes);
        System.out.println("Total genes: " + genes.size());
    }
    
    @Test
    public void hierarchicalClusterMutatedGenes() throws IOException {
        Map<String, Set<String>> sampleToGenes = loadSampleToMutatedGenes(cancerType, 
                                                                          true);
        Set<String> selectedGenes = selectGenesInSamples(2,
                                                         sampleToGenes);
        System.out.println("Total genes: " + selectedGenes.size());
        hierarchicalClusterGenes(new ArrayList<String>(selectedGenes));
    }
    
    @Test
    public void outputHitGenes() throws IOException {
        Map<String, Set<String>> sampleToGenes = loadSampleToMutatedGenes(cancerType, 
                                                                          true);
        Set<String> selectedGenes = selectGenesInSamples(2,
                                                         sampleToGenes);
        System.out.println("Total genes: " + selectedGenes.size());
        String fileName = WORKING_DIR + cancerType + "HitGenes030310.txt";
        fu.saveInteractions(selectedGenes, fileName);
    }
    
    @Test
    public void outputGeneToSampleCount() throws IOException {
        Map<String, Set<String>> sampleToGenes = loadSampleToMutatedGenes(cancerType, 
                                                                          true);
        Map<String, Set<String>> geneToSamples = InteractionUtilities.switchKeyValues(sampleToGenes);
        Set<String> selectedGenes = selectGenesInSamples(2,
                                                         sampleToGenes);
        System.out.println("Total genes: " + selectedGenes.size());
        String fileName = WORKING_DIR + cancerType + "GeneToSampleCount030310.txt";
        fu.setOutput(fileName);
        fu.printLine("Gene\tSampleCount");
        for (String gene : selectedGenes)
            fu.printLine(gene  + "\t" + geneToSamples.get(gene).size());
        fu.close();
    }
    
    @Test
    public void outputAllGeneToSampleCount() throws IOException {
        Map<String, Set<String>> sampleToGenes = loadSampleToMutatedGenes(cancerType, 
                                                                          false);
        Map<String, Set<String>> geneToSamples = InteractionUtilities.switchKeyValues(sampleToGenes);
        System.out.println("Total genes: " + geneToSamples.size());
        String fileName = WORKING_DIR + cancerType + "AllGeneToSampleCount031510.txt";
        fu.setOutput(fileName);
        fu.printLine("Gene\tSampleCount");
        for (String gene : geneToSamples.keySet())
            fu.printLine(gene  + "\t" + geneToSamples.get(gene).size());
        fu.close();
    }
    
    @Test
    public void calculateAveragePathInClusters() throws IOException {
        String clusterListFileName = WORKING_DIR + "ClusterListFor" + cancerType + "Cancer030210.txt";
        calculateAverageShortestPathInClusters(clusterListFileName);
    }

}
