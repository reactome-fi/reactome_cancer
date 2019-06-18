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
import org.reactome.annotate.PathwayBasedAnnotator;
import org.reactome.r3.graph.SpectralPartitionNetworkCluster;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.R3Constants;

/**
 * Many methods related to the Pancreatic cancers are grouped in this class.
 * @author wgm
 *
 */
public class PancreaticCancerAnalyzer extends CancerResequenceDataSetAnalyzer implements CancerConstants {
    public static final String PANCREATIC_DATA_DIR = PANCREATIC_DIR + "PancreaticDataset/";
    
    public PancreaticCancerAnalyzer() {
    }
    
    @Test
    public void annotateFIsInClusters() throws Exception {
        String fileName = PANCREATIC_DIR + "FIsInClusterList030210.txt";
        Set<String> fis = fu.loadInteractions(fileName);
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        PathwayBasedAnnotator annotator = new PathwayBasedAnnotator();
        annotator.setPValueThreshold(0.01);
        annotator.setFDRThreshold(0.25);
        annotator.annotateGenesWithFDR(fiGenes, AnnotationType.Pathway);
    }
    
    @Test
    public void outputHitGenes() throws IOException {
        Map<String, Set<String>> sampleToGenes = loadSampleToSequenceAlteredGenes();
        Set<String> selectedGenes = selectGenesInSamples(2,
                                                         sampleToGenes);
        System.out.println("Total genes: " + selectedGenes.size());
        String fileName = PANCREATIC_DIR + "HitGenes030310.txt";
        fu.saveInteractions(selectedGenes, fileName);
    }
    
    @Test
    public void outputGeneToSampleCount() throws IOException {
        Map<String, Set<String>> sampleToGenes = loadSampleToSequenceAlteredGenes();
        Map<String, Set<String>> geneToSamples = InteractionUtilities.switchKeyValues(sampleToGenes);
        Set<String> selectedGenes = selectGenesInSamples(2,
                                                         sampleToGenes);
        System.out.println("Total genes: " + selectedGenes.size());
        String fileName = PANCREATIC_DIR + "GeneToSampleCount030310.txt";
        fu.setOutput(fileName);
        fu.printLine("Gene\tSampleCount");
        for (String gene : selectedGenes)
            fu.printLine(gene  + "\t" + geneToSamples.get(gene).size());
        fu.close();
    }
    
    @Test
    public void outputAllGeneToSampleCount() throws IOException {
        Map<String, Set<String>> sampleToGenes = loadSampleToSequenceAlteredGenes();
        Map<String, Set<String>> geneToSamples = InteractionUtilities.switchKeyValues(sampleToGenes);
        System.out.println("Total genes: " + geneToSamples.size());
        String fileName = PANCREATIC_DIR + "AllGeneToSampleCount031510.txt";
        fu.setOutput(fileName);
        fu.printLine("Gene\tSampleCount");
        for (String gene : geneToSamples.keySet())
            fu.printLine(gene  + "\t" + geneToSamples.get(gene).size());
        fu.close();
    }
    
    @Test
    public void permutationTestForSamplesInNetworkClusters() throws Exception {
        String clusterFileName = PANCREATIC_DIR + "ClustersForAlteredGenes030110.txt";
        Map<String, Set<String>> sampleToAlteredGenes = loadSampleToSequenceAlteredGenes();
//        Map<String, Set<String>> sampleToAlteredGenes = loadSampleToMutations();
        System.out.println("Total samples: " + sampleToAlteredGenes.size());
        NetworkClusterAnalyzer clusterAnalyzer = new NetworkClusterAnalyzer();
        List<Integer> targetClusters = new ArrayList<Integer>();
        targetClusters.add(0);
        targetClusters.add(1);
        targetClusters.add(3);
        double target = 0.5416666666666666;
        System.out.println("\nPermutation for samples in clusters 0, 1, 3:");
        clusterAnalyzer.permutationTestForSamplesInNetworkClusters(sampleToAlteredGenes,
                                                                   clusterFileName,
                                                                   targetClusters, 
                                                                   target,
                                                                   false);
        // Check for clusters 0, 1:
        targetClusters.clear();
        targetClusters.add(0);
        targetClusters.add(1);
        target = 0.75;
        System.out.println("\nPermutation for samples in clusters 0, 1 for 0.75:");
        clusterAnalyzer.permutationTestForSamplesInNetworkClusters(sampleToAlteredGenes,
                                                                   clusterFileName,
                                                                   targetClusters, 
                                                                   target,
                                                                   false);
        // Check for clusters 0, 2:
        targetClusters.clear();
        targetClusters.add(0);
        targetClusters.add(2);
        target = 0.7083333333333334;
        System.out.println("\nPermutation for samples in clusters 0, 2 for 0.71:");
        clusterAnalyzer.permutationTestForSamplesInNetworkClusters(sampleToAlteredGenes,
                                                                   clusterFileName,
                                                                   targetClusters, 
                                                                   target,
                                                                   false);
    }
    
    /**
     * This method is used to check sample distribution in network clusters.
     * @throws Exception
     */
    @Test
    public void checkSamplesInNetworkClusters() throws Exception {
        String clusterFileName = PANCREATIC_DIR + "ClustersForAlteredGenes030110.txt";
        Map<String, Set<String>> sampleToAlteredGenes = loadSampleToSequenceAlteredGenes();
        //Map<String, Set<String>> sampleToAlteredGenes = loadSampleToMutations();
        System.out.println("Total samples: " + sampleToAlteredGenes.size());
        NetworkClusterAnalyzer clusterAnalyzer = new NetworkClusterAnalyzer();
        List<Integer> targetClusters = new ArrayList<Integer>();
        targetClusters.add(0);
        targetClusters.add(1);
        targetClusters.add(3);
        double percentage = clusterAnalyzer.checkSamplesInClusters(sampleToAlteredGenes, 
                                                                   clusterFileName,
                                                                   targetClusters);
        System.out.println("Samples in clusters 0, 1, 3: " + percentage);
        // Check for clusters 0, 1:
        targetClusters.clear();
        targetClusters.add(0);
        targetClusters.add(1);
        percentage = clusterAnalyzer.checkSamplesInClusters(sampleToAlteredGenes, 
                                                            clusterFileName,
                                                            targetClusters);
        System.out.println("Samples in clusters 0, 1: " + percentage);
        // Check for clusters 0, 2:
        targetClusters.clear();
        targetClusters.add(0);
        targetClusters.add(2);
        percentage = clusterAnalyzer.checkSamplesInClusters(sampleToAlteredGenes, 
                                                            clusterFileName,
                                                            targetClusters);
        System.out.println("Samples in clusters 0, 2: " + percentage);
    }
    
    
    @Test
    public void annotateClusters() throws Exception {
//        Map<String, Set<String>> sampleToGenes = loadSampleToSequenceAlteredGenes();
//        Set<String> allGenes = InteractionUtilities.grepAllGenes(sampleToGenes);
//        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
//        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
//        allGenes.retainAll(fiGenes);
//        System.out.println("All genes in FIs for randomization: " + allGenes.size());
        String clusterFileName = PANCREATIC_DIR + "ClustersForAlteredGenes030110.txt";
        NetworkClusterAnalyzer clusterAnalyzer = new NetworkClusterAnalyzer();
        clusterAnalyzer.annotateNetworkClusters(clusterFileName, 4);
//        clusterAnalyzer.annotateNetworkClusters(clusterFileName, 
//                                                4,
//                                                allGenes);
    }
    
    @Test
    public void clusterAlteredGenes() throws Exception {
        //        filterSampleToGenes(sampleToAlteredGenes, 2);
        //        //Map<String, Set<String>> sampleToAlteredGenes = loadSampleToMutations();
        Map<String, Set<String>> sampleToAlteredGenes = loadSampleToSequenceAlteredGenes();
        System.out.println("Total samples: " + sampleToAlteredGenes.size());
        Set<String> allGenes = new HashSet<String>();
        for (Set<String> set : sampleToAlteredGenes.values())
            allGenes.addAll(set);
        System.out.println("Total altered genes: " + allGenes.size());
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> fisInGenes = InteractionUtilities.getFIs(allGenes, fis);
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fisInGenes);
        System.out.println("Total genes in FIs: " + fiGenes.size());
//        if (true)
//            return;
//        String fiFileName = PANCREATIC_DATA_DIR + "FIsInAllAlteredGenes031510.txt";
//        fu.saveInteractions(fisInGenes,
//                            fiFileName);
        NetworkClusterAnalyzer cluster = new NetworkClusterAnalyzer();
        long time1 = System.currentTimeMillis();
        List<Set<String>> clusters = cluster.cluster(allGenes, 0.120);
        long time2 = System.currentTimeMillis();
        System.out.println("Time for edge betweenness: " + (time2 - time1));
        double modularity = cluster.calculateModularity(clusters, fisInGenes);
        System.out.println("Modularity from edge betweenness: " + modularity);
        
        SpectralPartitionNetworkCluster spectralClustering = new SpectralPartitionNetworkCluster();
        modularity = spectralClustering.calculateModualarity(clusters, fisInGenes);
        System.out.println("Modularity using new: " + modularity);
        long time11 = System.currentTimeMillis();
        clusters = spectralClustering.cluster(fisInGenes);
        long time21 = System.currentTimeMillis();
        System.out.println("\nTime for spectral: " + (time21 - time11));
        modularity = cluster.calculateModularity(clusters, fisInGenes);
        System.out.println("Modularity from spectral: " + modularity);
        modularity = spectralClustering.calculateModualarity(clusters, fisInGenes);
        System.out.println("Modularity using new: " + modularity);
        if (true)
            return;

        //System.out.println("Total clusters: " + clusters.size());
//        String clusterFileName = PANCREATIC_DIR + "ClustersForAlteredGenes030110.txt";
        String clusterFileName = PANCREATIC_DIR + "ClustersForAlteredGenesSpectral030110.txt";
        cluster.outputNetworkClusters(clusters, 
                                      clusterFileName);
        //List<Set<String>> clusters = cluster.loadNetworkClusters(clusterFileName);
        int index = 0;
        for (Set<String> set : clusters){
            System.out.println(index + ": " + set.size());
            index ++;
        }
//        //cluster.annotateNetworkClusters(clusters, 4);
    }
    
    public Map<String, Set<String>> loadSampleToSequenceAlteredGenes() throws IOException {
        Map<String, Set<String>> sampleToMutatedGenes = loadSampleToMutations();
        Map<String, Set<String>> sampleToCNVGenes = loadSampleToCNVGenes();
        Map<String, Set<String>> sampleToAlteredGenes = new HashMap<String, Set<String>>();
        for (String sample : sampleToCNVGenes.keySet()) {
            Set<String> cnvGenes = sampleToCNVGenes.get(sample);
            Set<String> mutatedGenes = sampleToMutatedGenes.get(sample.toLowerCase());
            Set<String> all = new HashSet<String>(mutatedGenes);
            all.addAll(cnvGenes);
            sampleToAlteredGenes.put(sample, all);
        }
        return sampleToAlteredGenes;
    }
    
    public Map<String, Set<String>> loadSampleToCNVGenes() throws IOException {
        String[] fileNames = new String[] {
                PANCREATIC_DATA_DIR + "HomoDeletion.txt",
                PANCREATIC_DATA_DIR + "Amplifications.txt"
        };
        return loadAllScienceSampleToCNVGenes(fileNames, 6, 8);
    }
    
    /**
     * This method is used to load sample to mutated genes map.
     * @return
     * @throws IOException
     */
    public Map<String, Set<String>> loadSampleToMutations() throws IOException {
        String fileName = PANCREATIC_DATA_DIR + "SomaticMutationsInDiscoveryScreen.txt";
        Map<String, Set<String>> sampleToGenes = new HashMap<String, Set<String>>();
        fu.setInput(fileName);
        String line = fu.readLine();
        Set<String> mutationTypes = new HashSet<String>();
        List<String> excludedTypes = new ArrayList<String>();
        excludedTypes.add("3'UTR");
        excludedTypes.add("5'UTR");
        excludedTypes.add("Synonymous");
        while ((line = fu.readLine()) != null) {
            if (line.length() == 0)
                break;
            String[] tokens = line.split("\t");
            String type = tokens[6];
            if (excludedTypes.contains(type))
                continue;
            mutationTypes.add(type);
            // Case should not be sensitive
            String sample = tokens[2].toLowerCase();
            String gene = tokens[0];
            InteractionUtilities.addElementToSet(sampleToGenes, sample, gene);
        }
        fu.close();
        System.out.println("Mutation types: " + mutationTypes);
        return sampleToGenes;
    }
    
    @Test
    public void testLoadSampleToAlteredGenes() throws IOException {
        Map<String, Set<String>> sampleToGenes = loadSampleToMutations();
        System.out.println("Total samples: " + sampleToGenes.size());
        List<String> samples = new ArrayList<String>(sampleToGenes.keySet());
        Collections.sort(samples);
        for (String sample : samples)
            System.out.println(sample);
        Set<String> totalGenes = new HashSet<String>();
        for (Set<String> set : sampleToGenes.values())
            totalGenes.addAll(set);
        System.out.println("Total mutated genes: " + totalGenes.size());
        Map<String, Set<String>> sampleToCNVGenes = loadSampleToCNVGenes();
        System.out.println("Total samples in CNVs: " + sampleToCNVGenes.size());
        for (Set<String> set : sampleToCNVGenes.values())
            totalGenes.addAll(set);
        System.out.println("Total mutated and CNV genes: " + totalGenes.size());
        Map<String, Set<String>> sampleToAlteredGenes = loadSampleToSequenceAlteredGenes();
        summarizeSampleToGenes(sampleToAlteredGenes);
    }
    
    @Test
    public void checkAlteredGenesCoverage() throws IOException {
        Map<String, Set<String>> sampleToAlteredGenes = loadSampleToSequenceAlteredGenes();
        for (String sample : sampleToAlteredGenes.keySet())
            System.out.println(sample + "\t" + sampleToAlteredGenes.get(sample).size());
        Set<String> allGenes = new HashSet<String>();
        for (Set<String> set : sampleToAlteredGenes.values())
            allGenes.addAll(set);
        System.out.println("Total genes: " + allGenes.size());
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        Set<String> shared = InteractionUtilities.getShared(fiGenes, allGenes);
        System.out.println("Total shared: " + shared.size() + " (" + shared.size() / (double) allGenes.size() + ")");
//        allGenes.removeAll(shared);
//        List<String> geneList = new ArrayList<String>(allGenes);
//        for (String gene : geneList)
//            System.out.println(gene);
        int total = allGenes.size();
        // Do a filtering
        filterSampleToGenes(sampleToAlteredGenes, 2);
        allGenes.clear();
        for (Set<String> set : sampleToAlteredGenes.values())
            allGenes.addAll(set);
        System.out.println("Total genes in 2 samples: " + allGenes.size());
        double percentage = (double) allGenes.size() / total;
        System.out.println("Percentage: " + percentage);
        // Check Science GBM genes
        Map<String, Set<String>> gbmSampleToAlteredGenes = getScienceGBMSampleToAlteredGenes();
        for (String sample : gbmSampleToAlteredGenes.keySet())
            System.out.println(sample + "\t" + gbmSampleToAlteredGenes.get(sample).size());
        allGenes.clear();
        for (Set<String> set : gbmSampleToAlteredGenes.values())
            allGenes.addAll(set);
        System.out.println("\nTotal genes in GBM: " + allGenes.size());
        total = allGenes.size();
        // Do a filtering
        filterSampleToGenes(gbmSampleToAlteredGenes, 2);
        allGenes.clear();
        for (Set<String> set : gbmSampleToAlteredGenes.values())
            allGenes.addAll(set);
        System.out.println("Total genes in GBM 2 samples: " + allGenes.size());
        percentage = (double) allGenes.size() / total;
        System.out.println("Percentage: " + percentage);
    }
    
    @Test
    public void hierarchicalClusterMutatedGenes() throws IOException {
        Map<String, Set<String>> sampleToGenes = loadSampleToSequenceAlteredGenes();
        Set<String> selectedGenes = selectGenesInSamples(2,
                                                         sampleToGenes);
        System.out.println("Total genes: " + selectedGenes.size());
        hierarchicalClusterGenes(new ArrayList<String>(selectedGenes));
    }
    
    @Test
    public void calculateAveragePathInClusters() throws IOException {
        String clusterListFileName = PANCREATIC_DIR + "ClusterList030210.txt";
        calculateAverageShortestPathInClusters(clusterListFileName);
    }
    
    @Test
    public void annotateHierarchicalClusters() throws Exception {
        String clusterFileName = PANCREATIC_DIR + "ClusterList030210.txt";
        System.out.println("Cluster file: " + clusterFileName);
        Map<String, List<String>> headerToIds = loadClusters(clusterFileName);
        for (String tmp : headerToIds.keySet()) {
            System.out.println(tmp);
            List<String> names = headerToIds.get(tmp);
            annotateNames(names);
            System.out.println();
        }
    }
}
