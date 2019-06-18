/*
 * Created on Jun 9, 2009
 *
 */
package org.reactome.cancer;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math.random.RandomData;
import org.apache.commons.math.random.RandomDataImpl;
import org.gk.util.StringUtils;
import org.junit.Test;
import org.reactome.annotate.AnnotationType;
import org.reactome.annotate.PathwayBasedAnnotator;
import org.reactome.data.GODataAnalyzer;
import org.reactome.data.GODataAnalyzerV2;
import org.reactome.r3.cluster.HierarchicalClusterNode;
import org.reactome.r3.graph.JungGraphUtilities;
import org.reactome.r3.graph.NetworkModularityCalculator;
import org.reactome.r3.graph.SpectralPartitionNetworkCluster;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.MathUtilities;
import org.reactome.r3.util.R3Constants;

import edu.uci.ics.jung.algorithms.cluster.EdgeBetweennessClusterer;
import edu.uci.ics.jung.algorithms.importance.BetweennessCentrality;
import edu.uci.ics.jung.algorithms.importance.Ranking;
import edu.uci.ics.jung.graph.Graph;

/**
 * This class is used to do network based clustering analysis.
 * @author wgm
 *
 */
public class NetworkClusterAnalyzer {
    private FileUtility fu = new FileUtility();
    private JungGraphUtilities gu = new JungGraphUtilities();
    
    public NetworkClusterAnalyzer() {
    }
    
    /**
     * This method is used to check overlapping between two network clustering results.
     * @param modules1
     * @param modules2
     * @param size
     * @param genes
     */
    public void checkNetworkModuleOverlapping(List<Set<String>> modules1,
                                              List<Set<String>> modules2,
                                              int size,
                                              double pvalueCutff,
                                              Set<String> genes) throws Exception {
        System.out.println("Module1\tSize1\tModule2\tSize2\tShared\tP-Value");
        for (int i = 0; i < modules1.size(); i++) {
            Set<String> module1 = modules1.get(i);
            if (module1.size() < size)
                continue;
            for (int j = 0; j < modules2.size(); j++) {
                Set<String> module2 = modules2.get(j);
                if (module2.size() < size)
                    continue;
                Set<String> shared = InteractionUtilities.getShared(module1, module2);
                double pvalue = MathUtilities.calculateHypergeometricPValue(genes.size(),
                                                                            module1.size(),
                                                                            module2.size(),
                                                                            shared.size());
                if (pvalue <= pvalueCutff)
                    System.out.println(i + "\t" + module1.size() + "\t" +
                                       j + "\t" + module2.size() + "\t" +
                                       shared.size() + "\t" + pvalue);
            }
        }
    }
    
    public void permutationTestBasedOnModularity(Map<String, Set<String>> sampleToAlteredGenes) throws Exception {
        // Will focus on FI genes only since only we know these genes
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        //Set<String> fis = fu.loadInteractions(R3Constants.RESULT_DIR + "FIsInGene_Pathway_041709.txt");
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        Set<String> totalGenes = new HashSet<String>();
        for (String sample : sampleToAlteredGenes.keySet()) {
            Set<String> genes = sampleToAlteredGenes.get(sample);
            genes.retainAll(fiGenes);
            //System.out.println(sample + ": " + genes.size());
        }
        for (Set<String> set : sampleToAlteredGenes.values())
            totalGenes.addAll(set);
        System.out.println("Total genes: " + totalGenes.size());
        Set<String> fisInGenes = InteractionUtilities.getFIs(totalGenes, fis);
        SpectralPartitionNetworkCluster clustering = new SpectralPartitionNetworkCluster();
        List<Set<String>> originalClusters = clustering.cluster(fisInGenes);
        NetworkModularityCalculator calculator = new NetworkModularityCalculator();
        double wModularity = calculator.calculateWeightedModularity(originalClusters, fisInGenes);
        double modularity = calculateModularity(originalClusters, fisInGenes);
        System.out.println("Original modularity: " + modularity + " (" + originalClusters.get(0).size() + ")");
        System.out.println("weighted: " + wModularity);
        // Print out the sizes
        for (Set<String> cluster : originalClusters) {
            modularity = calculator.calculateModularity(cluster, fisInGenes);
            System.out.println(cluster.size() + "\t" + modularity);
        }
        System.out.println();
        int permutationNumber = 100;
        List<Double> permutatedModularities = new ArrayList<Double>();
        List<Double> wPermutatedModularities = new ArrayList<Double>();
        List<Integer> biggestSizes = new ArrayList<Integer>();
        List<Integer> permutatedSizes = new ArrayList<Integer>();
        List<Double> individualModularities = new ArrayList<Double>();
        //Map<String, Set<String>> sampleToPermuGenes = new HashMap<String, Set<String>>();
        RandomData randomizer = new RandomDataImpl();
        for (int j = 0; j < permutationNumber; j++) {
            // Need to generate a same number of genes
//            for (String sample : sampleToAlteredGenes.keySet()) {
//                Set<String> genes = sampleToAlteredGenes.get(sample);
//                Set<String> permutated = null;
//                if (genes.size() == 0)
//                    permutated = new HashSet<String>();
//                else {
//                    permutated = MathUtilities.randomSampling(totalGenes, 
//                                                              genes.size(),
//                                                              randomizer);
//                }
//                sampleToPermuGenes.put(sample, permutated);
//            }
//            Set<String> allSampleGenes = InteractionUtilities.grepAllGenes(sampleToPermuGenes);
            Set<String> allSampleGenes = MathUtilities.randomSampling(fiGenes, 
                                                                      totalGenes.size());
            Set<String> fisInSampleGenes = InteractionUtilities.getFIs(allSampleGenes, fis);
//            List<Set<String>> permuClusters = cluster(allSampleGenes,
//                                                      0.10,
//                                                      fisInSampleGenes);
            List<Set<String>> permuClusters = clustering.cluster(fisInSampleGenes);
            double permModularity = calculator.calculateWeightedModularity(permuClusters, 
                                                         fisInSampleGenes);
            //sampleToPermuGenes.clear();
            System.out.println("Permutation " + j + ": " + permModularity + " (" + 
                               permuClusters.get(0).size() + ")");
            wPermutatedModularities.add(permModularity);
            permModularity = calculateModularity(permuClusters,
                                                 fisInSampleGenes);
            permutatedModularities.add(permModularity);
            for (Set<String> cluster : permuClusters) {
                permutatedSizes.add(cluster.size());
                modularity = calculator.calculateModularity(cluster, fisInSampleGenes);
                individualModularities.add(modularity);
            }
            biggestSizes.add(permuClusters.get(0).size());
        }
//        Collections.sort(permutatedModularities, new Comparator<Double>() {
//            public int compare(Double v1, Double v2) {
//                return v2.compareTo(v1);
//            }
//        });
//        Collections.sort(permutatedSizes, new Comparator<Integer>() {
//            public int compare(Integer size1, Integer size2) {
//                return size2.compareTo(size1);
//            }
//        });
//        Collections.sort(individualModularities, new Comparator<Double>() {
//            public int compare(Double v1, Double v2) {
//                return v2.compareTo(v1);
//            }
//        });
        System.out.println("\nModularities from permutations:");
        System.out.println("Modularity\tWeighted_Modularty\tSize_of_Largest_Module");
        for (int i = 0; i < permutatedModularities.size(); i++) {
            System.out.println(permutatedModularities.get(i) + "\t" + 
                               wPermutatedModularities.get(i) + "\t" + 
                               biggestSizes.get(i));
        }
        System.out.println("\nSizes from permutations:");
        for (Integer size : permutatedSizes)
            System.out.println(size);
        System.out.println("\nModularity for individual cluster:");
        System.out.println("Size\tModularity");
        for (int i = 0; i < individualModularities.size(); i++)
            System.out.println(permutatedSizes.get(i) + "\t" +
                               individualModularities.get(i));
    }
    
    public void permutationTestForSamplesInClustersOnBetweenness(Map<String, Set<String>> sampleToAlteredGenes,
                                                                 List<Integer> targetClusters,
                                                                 List<Integer> targetClusterSizes) throws Exception {
        // Will focus on FI genes only since only we know these genes
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        //Set<String> fis = fu.loadInteractions(R3Constants.RESULT_DIR + "FIsInGene_Pathway_041709.txt");
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        Set<String> totalGenes = new HashSet<String>();
        for (String sample : sampleToAlteredGenes.keySet()) {
            Set<String> genes = sampleToAlteredGenes.get(sample);
            genes.retainAll(fiGenes);
            //System.out.println(sample + ": " + genes.size());
        }
        for (Set<String> set : sampleToAlteredGenes.values())
            totalGenes.addAll(set);
        System.out.println("Total genes: " + totalGenes.size());
        int permutationNumber = 100;
        List<Double> permutPercentages = new ArrayList<Double>();
        Map<String, Set<String>> sampleToPermuGenes = new HashMap<String, Set<String>>();
        RandomData randomizer = new RandomDataImpl();
        for (int j = 0; j < permutationNumber; j++) {
            // Need to generate a same number of genes
            for (String sample : sampleToAlteredGenes.keySet()) {
                Set<String> genes = sampleToAlteredGenes.get(sample);
                Set<String> permutated = null;
                if (genes.size() == 0)
                    permutated = new HashSet<String>();
                else
                    permutated = MathUtilities.randomSampling(totalGenes, 
                                                              genes.size(),
                                                              randomizer);
                sampleToPermuGenes.put(sample, permutated);
            }
            Set<String> allSampleGenes = InteractionUtilities.grepAllGenes(sampleToPermuGenes);
            Set<String> fisInSampleGenes = InteractionUtilities.getFIs(allSampleGenes, fis);
            List<Set<String>> permuClusters = cluster(allSampleGenes,
                                                      0.190,
                                                      fisInSampleGenes);
            boolean isSizePermitted = true;
            for (int i = 0; i < targetClusters.size(); i++) {
                int index = targetClusters.get(i);
                int size = targetClusterSizes.get(i);
                int randomSize = permuClusters.get(index).size();
                if (randomSize < size) {
                    isSizePermitted = false;
                    break;
                }
            }
            // Check how many samples having genes in target clusters
            int counter = 0;
            if (isSizePermitted) {
                for (String sample : sampleToPermuGenes.keySet()) {
                    Set<String> genes = sampleToPermuGenes.get(sample);
                    // Want to get the list of clusters
                    List<Integer> clusterIds = new ArrayList<Integer>();
                    for (int i = 0; i < permuClusters.size(); i++) {
                        Set<String> set = permuClusters.get(i);
                        for (String gene : genes) {
                            if (set.contains(gene)) {
                                clusterIds.add(i);
                                break;
                            }
                        }
                    }
                    if (isSampleInClusters(clusterIds, 
                                           targetClusters))
                        counter ++;
                }
            }
            double percent = (double) counter / sampleToPermuGenes.size();
            permutPercentages.add(percent);
            sampleToPermuGenes.clear();
            System.out.println("Permutation " + j + ": " + percent);
        }
        Collections.sort(permutPercentages, new Comparator<Double>() {
            public int compare(Double v1, Double v2) {
                return v2.compareTo(v1);
            }
        });
        for (int i = 0; i < permutPercentages.size(); i++)
            System.out.println(i + ": " + permutPercentages.get(i));
        //        calculatePValue(target, 
        //                        permutationNumber, 
        //                        permutPercentages);
        //        String outFileName = TCGA_GBM_DIR + "ClustersInTCGAAlteredGenes071609Permutation10000_071609.txt";
        //        fu.setOutput(outFileName);
        //        for (Double d : permutPercentages)
        //            fu.printLine(d + "");
        //        fu.close();
    }
    
    public void permutationTestForSamplesInNetworkClusters(Map<String, Set<String>> sampleToAlteredGenes,
                                                           String clusterFileName,
                                                           List<Integer> targetClusters,
                                                           double target,
                                                           boolean permutateClusters) throws Exception {
        List<Set<String>> clusters = loadNetworkClusters(clusterFileName);
        permutationTestForSamplesInNetworkClusters(sampleToAlteredGenes,
                                                   clusters,
                                                   targetClusters, 
                                                   target,
                                                   permutateClusters);
    }

    public void permutationTestForSamplesInNetworkClusters(Map<String, Set<String>> sampleToAlteredGenes,
                                                           List<Set<String>> clusters,
                                                           List<Integer> targetClusters,
                                                           double target,
                                                           boolean permutateClusters)
            throws IOException {
        // Will focus on FI genes only since only we know these genes
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        //Set<String> fis = fu.loadInteractions(R3Constants.RESULT_DIR + "FIsInGene_Pathway_041709.txt");
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        Set<String> totalGenes = new HashSet<String>();
        for (String sample : sampleToAlteredGenes.keySet()) {
            Set<String> genes = sampleToAlteredGenes.get(sample);
            genes.retainAll(fiGenes);
            //System.out.println(sample + ": " + genes.size());
            System.out.print(genes.size() + ", ");
        }
        System.out.println();
        for (Set<String> set : sampleToAlteredGenes.values())
            totalGenes.addAll(set);
        System.out.println("Total genes: " + totalGenes.size());
        
//        int index = 0;
//        for (Set<String> cluster : clusters) {
//            System.out.println(index + ": " + cluster.size());
//            index ++;
//        }
        int permutationNumber = 10000;
        List<Double> permutPercentages = new ArrayList<Double>();
        Map<String, Set<String>> sampleToPermuGenes = new HashMap<String, Set<String>>();
        RandomData randomizer = new RandomDataImpl();
        if (permutateClusters) {
            for (int j = 0; j < permutationNumber; j++) {
                // Do a quick check using random clusters
                List<Set<String>> randomClusters = generateRandomClusters(clusters, 
                                                                          totalGenes, 
                                                                          randomizer);
                int counter1 = 0;
                for (String sample : sampleToAlteredGenes.keySet()) {
                    Set<String> genes = sampleToAlteredGenes.get(sample);
                    // Want to get the list of clusters
                    List<Integer> clusterIds = new ArrayList<Integer>();
                    for (int i = 0; i < randomClusters.size(); i++) {
                        Set<String> set = randomClusters.get(i);
                        for (String gene : genes) {
                            if (set.contains(gene)) {
                                clusterIds.add(i);
                                break;
                            }
                        }
                    }
                    if (isSampleInClusters(clusterIds, 
                                           targetClusters))
                        counter1 ++;
                }
                double percent = (double) counter1 / sampleToAlteredGenes.size();
                permutPercentages.add(percent);
            }
        }
        else {
//            // The following is testing code
//            sampleToAlteredGenes = new CancerResequenceDataSetAnalyzer().getScienceGBMSampleToAlteredGenes();
//            for (String sample : sampleToAlteredGenes.keySet()) {
//                Set<String> set = sampleToAlteredGenes.get(sample);
//                set.retainAll(fiGenes);
//            }
            
            for (int j = 0; j < permutationNumber; j++) {
                // Need to generate a same number of genes
                for (String sample : sampleToAlteredGenes.keySet()) {
                    Set<String> genes = sampleToAlteredGenes.get(sample);
                    Set<String> permutated = null;
                    if (genes.size() == 0)
                        permutated = new HashSet<String>();
                    else
                        permutated = MathUtilities.randomSampling(totalGenes, 
                                                                  genes.size(),
                                                                  randomizer);
                    sampleToPermuGenes.put(sample, permutated);
                }
                // Check how many samples having genes in target clusters
                int counter = 0;
                for (String sample : sampleToPermuGenes.keySet()) {
                    Set<String> genes = sampleToPermuGenes.get(sample);
                    // Want to get the list of clusters
                    List<Integer> clusterIds = new ArrayList<Integer>();
                    for (int i = 0; i < clusters.size(); i++) {
                        Set<String> set = clusters.get(i);
                        for (String gene : genes) {
                            if (set.contains(gene)) {
                                clusterIds.add(i);
                                break;
                            }
                        }
                    }
                    if (isSampleInClusters(clusterIds, 
                                           targetClusters))
                        counter ++;
                }
                double percent = (double) counter / sampleToPermuGenes.size();
                permutPercentages.add(percent);
                sampleToPermuGenes.clear();
            }
        }
        Collections.sort(permutPercentages, new Comparator<Double>() {
            public int compare(Double v1, Double v2) {
                return v2.compareTo(v1);
            }
        });
//        for (int i = 0; i < permutPercentages.size(); i++)
//            System.out.println(i + ": " + permutPercentages.get(i));
        calculatePValue(target, 
                        permutationNumber, 
                        permutPercentages);
//        String outFileName = TCGA_GBM_DIR + "ClustersInTCGAAlteredGenes071609Permutation10000_071609.txt";
//        fu.setOutput(outFileName);
//        for (Double d : permutPercentages)
//            fu.printLine(d + "");
//        fu.close();
    }
    
    private List<Set<String>> generateRandomClusters(List<Set<String>> clusters,
                                                     Set<String> totalGenes,
                                                     RandomData randomizer) {
        List<Set<String>> samples = new ArrayList<Set<String>>();
        Set<String> copy = new HashSet<String>(totalGenes);
        for (Set<String> cluster : clusters) {
            int size = cluster.size();
            Set<String> tmp = MathUtilities.randomSampling(copy, 
                                                           size,
                                                           randomizer);
            samples.add(tmp);
            copy.removeAll(tmp);
        }
        return samples;
    }
    
    private void calculatePValue(double target, 
                                 int permutationNumber,
                                 List<Double> permutPercentages) {
        int index = -1;
        for (int i = 0; i < permutPercentages.size(); i++) {
            double pert = permutPercentages.get(i);
            if (pert < target) {
                index = i;
                break;
            }
        }
        if (index == 0)
            System.out.println("pvalue < " + 1.0 / permutationNumber);
        else if (index == -1)
            System.out.println("pvalue = 1.0");
        else
            System.out.println("pvalue: " + (double) index / permutationNumber);
        double max = permutPercentages.get(0);
        System.out.println("Max value: " + max);
        System.out.println("Min value: " + permutPercentages.get(permutPercentages.size() - 1));
        System.out.println("Target value: " + target);
    }
    
    public double checkSamplesInClusters(Map<String, Set<String>> sampleToAlteredGenes,
                                         String clusterFileName,
                                         List<Integer> targetClusters) throws IOException {
        List<Set<String>> clusters = loadNetworkClusters(clusterFileName);
        System.out.println("Total clusters: " + clusters.size());
        // Check how many samples having 0 and 1 clusters
        int counter = 0;
        System.out.println("Sample\tClusters\tSharedGenes");
        List<String> posSamples = new ArrayList<String>();
        List<String> negSamples = new ArrayList<String>();
        for (String sample : sampleToAlteredGenes.keySet()) {
            Set<String> genes = sampleToAlteredGenes.get(sample);
            // Want to get the list of clusters
            List<Integer> clusterIds = new ArrayList<Integer>();
            List<String> sharedGenes = new ArrayList<String>();
            for (int i = 0; i < clusters.size(); i++) {
                Set<String> set = clusters.get(i);
                for (String gene : genes) {
                    if (set.contains(gene)) {
                        clusterIds.add(i);
                        break;
                    }
                }
                Set<String> shared = new HashSet<String>(set);
                shared.retainAll(genes);
                sharedGenes.add(StringUtils.join(":", new ArrayList<String>(shared)));
            }
            System.out.println(sample + "\t" + clusterIds + "\t" + sharedGenes);
            if (isSampleInClusters(clusterIds,
                                   targetClusters)) {
                counter ++;
                posSamples.add(sample);
            }
            else
                negSamples.add(sample);
        }
        double percent = (double) counter / sampleToAlteredGenes.size();
//        System.out.println("Total samples checked: " + sampleToAlteredGenes.size());
//        System.out.println("Total samples having 0 and 1: " + counter + "(" + percent + ")");
//        System.out.println("Total samples in cluster 0: " + cluster0);
//        System.out.println("Total samples in cluster 1: " + cluster1);
        return percent;
    }
    
    private boolean isSampleInClusters(List<Integer> clusterIds,
                                       List<Integer> targetClusters) {
        for (Integer target : targetClusters) {
            if (!clusterIds.contains(target))
                return false;
        }
        return true;
//        // Check for a sample has mutation in both cluster 0 and 1.
//        if (clusterIds.contains(0) && 
//            //clusterIds.contains(1) && 
//            clusterIds.contains(2))
//            return true;
////         Check for a sample has a mutation in either cluster 0 or 1.
////        if (clusterIds.contains(0) || clusterIds.contains(1))
////            return true;
////        if (clusterIds.contains(0))
////            return true;
//        return false;
    }
    
    /**
     * This method is used to label network clusters based on the highest centrality in 
     * each cluster.
     * @param clusters
     * @return
     */
    public List<String> labelNetworkClusters(List<Set<String>> clusters) throws IOException {
        return gu.labelNetworkClusters(clusters);
    }
    
    public Set<String> grepAllGenesInClusters(List<Set<String>> clusters) {
        Set<String> set = new HashSet<String>();
        for (Set<String> cluster : clusters)
            set.addAll(cluster);
        return set;
    }
    
    /**
     * This method is used to check network clustering results for a set of random genes
     * @throws Exception
     */
    @Test
    public void permutationTestForClusters() throws Exception {
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> genesInFIs = InteractionUtilities.grepIDsFromInteractions(fis);
        List<Double> ratios = new ArrayList<Double>();
        int sampleSize = 562;
        PathwayBasedAnnotator annotator = new PathwayBasedAnnotator();
        annotator.setFDRThreshold(0.25);
        annotator.setPValueThreshold(0.05);
        // Note: the whole FI network is used.
        GODataAnalyzerV2 goAnalyzer = new GODataAnalyzerV2();
        Map<String, Set<String>> proteinToTerms = goAnalyzer.loadProteinToGOCCTerms();
        // Need to convert protein in UniProt to gene names
//        Map<String, Set<String>> geneToGOTerms = annotator.convertProteinToGeneForGO(proteinToTerms);
//        // Convert GO Term ids to names
//        Map<String, Set<String>> geneToGONames = new GODataAnalyzer().convertIdsToNames(geneToGOTerms);
//        // Targets
//        String targetClusterFileName = CancerConstants.PANCREATIC_DIR + "ClustersForAlteredGenes030110.txt";
//        List<Set<String>> targetClusters = loadNetworkClusters(targetClusterFileName);
//        for (int i = 0; i < 5; i++) {
//            System.out.println("Target cluster " + i + ": " + targetClusters.get(i).size());
//        }
//        for (int i = 0; i < 100; i++) {
//            System.out.println("Permutation " + i);
//            Set<String> sample = MathUtilities.randomSampling(genesInFIs, sampleSize);
//            Set<String> fisInSample = InteractionUtilities.getFIs(sample, fis);
//            Set<String> connectectedGenes = InteractionUtilities.grepIDsFromInteractions(fisInSample);
//            double ratio = (double) connectectedGenes.size() / sampleSize;
//            ratios.add(ratio);
//            List<Set<String>> clusters = cluster(sample, 
//                                                 0.120,
//                                                 fis);
//            // Compare the first cluster
//            if (clusters.get(0).size() < targetClusters.get(0).size()) {
//                System.out.println("First cluster size is smaller: " + clusters.get(0).size());
//                continue;
//            }
//            int index = 0;
//            for (Set<String> cluster : clusters) {
//                if (cluster.size() >= 20) {
//                    System.out.println(index + ": " + cluster.size());
//                    //annotator.annotateGenesUsingGOWithFDR(cluster, "CC");
//                    annotator.annotateGeneSet(cluster,
//                                              geneToGONames);
//                }
//                index ++;
//            }
//            System.out.println();
//        }
//        Collections.sort(ratios);
//        for (Double d : ratios)
//            System.out.println(d);
    }
    
    public void annotateNetworkClusters(String fileName, 
                                        int clusterSizeCutoff) throws Exception {
        List<Set<String>> clusters = loadNetworkClusters(fileName);
        annotateNetworkClusters(clusters, clusterSizeCutoff);
    }
    
    public void annotateNetworkClusters(String fileName,
                                        int clusterSizeCutoff,
                                        Set<String> randomGenes) throws Exception {
        List<Set<String>> clusters = loadNetworkClusters(fileName);
        annotateNetworkClusters(clusters, 
                                clusterSizeCutoff,
                                randomGenes);
    }

    public void annotateNetworkClusters(List<Set<String>> clusters,
                                        int clusterSizeCutoff) throws Exception {
        annotateNetworkClusters(clusters,
                                clusterSizeCutoff, 
                                null);
    }
    
    public void annotateNetworkClusters(List<Set<String>> clusters,
                                        int clusterSizeCutoff,
                                        Set<String> randomGenes) throws Exception {
        PathwayBasedAnnotator annotator = new PathwayBasedAnnotator();
        annotator.setFDRThreshold(0.25);
        annotator.setPValueThreshold(0.05);
        annotator.setRandomGenes(randomGenes);
        int index = 0;
        for (Set<String> cluster : clusters) {
            if (cluster.size() < clusterSizeCutoff)
                continue;
            System.out.println("\nAnnotating cluster " + index + ": " + cluster.size());
            annotator.annotateGenesWithFDR(cluster, AnnotationType.Pathway);
//            annotator.annotateGenesUsingGOWithFDR(cluster, "CC");
            index ++;
        }
    }
    
    /**
     * Filter a list of network clusters based on a passed cluster size. After filtering,
     * clusters with sizes >= size are kept.
     * @param clusters
     * @param size
     */
    public void filterNetworkClusters(List<Set<String>> clusters,
                                      int size) {
        for (Iterator<Set<String>> it = clusters.iterator(); it.hasNext();) {
            Set<String> cluster = it.next();
            if (cluster.size() < size)
                it.remove();
        }
    }
    
    public List<Set<String>> loadNetworkClusters(String fileName) throws IOException {
        List<Set<String>> rtn = new ArrayList<Set<String>>();
        fu.setInput(fileName);
        String line = fu.readLine();
        String[] headTokens = line.split("\t");
        if (line.equals("Module (class=java.lang.Integer)")) {
            while ((line = fu.readLine()) != null) {
                int index = line.indexOf("=");
                String gene = line.substring(0, index);
                int cluster = new Integer(line.substring(index + 1));
                Set<String> set = null;
                if (cluster < rtn.size()) {
                    set = rtn.get(cluster);
                }
                else {
                    set = new HashSet<String>();
                    rtn.add(set);
                }
                set.add(gene);
            }
        }
        else if ((headTokens.length == 6 || headTokens.length == 5) &&
                 headTokens[0].equals("Module") &&
                 headTokens[headTokens.length - 1].equals("Node List")) { // Format copied from Reactome FI plug-in
            while ((line = fu.readLine()) != null) {
                String[] tokens = line.split("\t");
                tokens = tokens[headTokens.length - 1].split(",");
                Set<String> cluster = new HashSet<String>();
                for (String token : tokens)
                    cluster.add(token);
                rtn.add(cluster);
            }
        }
        else {
            // Based on MCL or other types
            String[] tokens = line.split("\t");
            Set<String> cluster = new HashSet<String>();
            for (String token : tokens)
                cluster.add(token);
            rtn.add(cluster);
            while ((line = fu.readLine()) != null) {
                tokens = line.split("\t");
                cluster = new HashSet<String>();
                for (String token : tokens)
                    cluster.add(token);
                rtn.add(cluster);
            }
        }
        fu.close();
        return rtn;
    }
    
    @Test
    public void clusterSubnetwork() throws Exception {
        String fiFileName = NatureGBMAnalyzer.TCGA_GBM_DIR + "FIsInHomoCNVMutatedGenes.txt";
        Set<String> fis = fu.loadInteractions(fiFileName);
        Set<String> genes = InteractionUtilities.grepIDsFromInteractions(fis);
        List<Set<String>> clusterList = cluster(genes);
        String outFileName = NatureGBMAnalyzer.TCGA_GBM_DIR + "ClustersInFIsInHomoCNVMutatedGenes070709.txt";
        outputNetworkClusters(clusterList, outFileName);
        double modularity = calculateModularity(clusterList, fis);
        System.out.println("Modularity: " + modularity);
    }
    
    @Test
    public void checkNetworkModuality() throws Exception {
        String fileName = NatureGBMAnalyzer.TCGA_GBM_DIR + "ClustersInFIsInHomoCNVMutatedGenes.txt";
        List<Set<String>> clusters = loadNetworkClusters(fileName);
        String fiFileName = NatureGBMAnalyzer.TCGA_GBM_DIR + "FIsInHomoCNVMutatedGenes.txt";
        Set<String> fis = fu.loadInteractions(fiFileName);
        double modularity = calculateModularity(clusters, fis);
        System.out.println("Modularity: " + modularity);
    }

    public void outputNetworkClusters(List<Set<String>> clusterList,
                                      String outFileName) throws IOException {
        fu.setOutput(outFileName);
        // Export as attribute file
        fu.printLine("Module (class=java.lang.Integer)");
        for (int i = 0; i < clusterList.size(); i++) {
            Set<String> cluster = clusterList.get(i);
            for (String gene : cluster)
                fu.printLine(gene + "=" + i);
        }
        fu.close();
    }
    
    /**
     * This method is used to do network clustering for a set of genes.
     * @param genes
     * @return
     * @throws Exception
     */
    public List<Set<String>> cluster(Collection<String> genes, double removeEdgeRatio) throws Exception {
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        return cluster(genes, removeEdgeRatio, fis);
    }

    public List<Set<String>> cluster(Collection<String> genes,
                                     double removeEdgeRatio, 
                                     Set<String> fis) {
        return gu.cluster(genes, removeEdgeRatio, fis);
    }
    
    /**
     * This is a very long process and should NOT run at a laptop or desktop computer.
     * @throws Exception
     */
    @Test
    public void checkEdgeBetweenness() throws Exception {
        Set<String> fis  = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        Graph<String, String> graph = gu.createJungGraph(fiGenes, fis);
        edu.uci.ics.jung.algorithms.scoring.BetweennessCentrality<String, String> betweenness = new edu.uci.ics.jung.algorithms.scoring.BetweennessCentrality<String, String>(graph);
        for (String edge : graph.getEdges()) {
            double score = betweenness.getEdgeScore(edge);
            System.out.println(edge + ": " + score);
        }
    }
    
    /**
     * This method tries to find an optimal numbers of removing edge by testing 100 ratios.
     * @param genes
     * @return
     * @throws Exception
     */
    public List<Set<String>> cluster(Collection<String> genes) throws Exception {
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        return cluster(genes, fis);
    }

    public List<Set<String>> cluster(Collection<String> genes, Set<String> fis) {
        return gu.cluster(genes, fis);
    }
    
    /**
     * Similar to generatePathwayClusters(), except the clustering methods are based on Jung
     * graph library.
     * @throws Exception
     */
    @Test
    public void generatePathwayClustersWithClustering() throws Exception {
//        Map<String, Set<String>> pathwayToGenes = new TopicAnalyzer().getTopicToNamesMap();
//        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_FILE_NAME);
//        Set<String> allGenes = InteractionUtilities.grepIDsFromInteractions(fis);
//        int c = 1;
//        String outFileName = R3Constants.RESULT_DIR + "Pathways_Reactome_28_plus_i_betweeness_clusters_v2.gmt";
//        fu.setOutput(outFileName);
//        for (String pathway : pathwayToGenes.keySet()) {
//            Set<String> genes = pathwayToGenes.get(pathway);
//            genes.retainAll(allGenes);
//            System.out.println(c + ": " + pathway + "(" + genes.size() + ")");
//            // Convert to FIs
//            Set<String> fisInGenes = InteractionUtilities.getFIs(genes, fis);
//            EdgeBetweennessClusterer<String, String> clusterer = new EdgeBetweennessClusterer<String, String>(fisInGenes.size() / 10);
//            Graph<String, String> graph = gu.createJungGraph(genes, fisInGenes);
//            Collection<Set<String>> clusters = clusterer.transform(graph);
//            //VoltageClusterer<String, String> clusterer = new VoltageClusterer<String, String>(graph, 10);
//            //Collection<Set<String>> clusters = clusterer.cluster(10);
//            List<Set<String>> clusterList = new ArrayList<Set<String>>(clusters);
//            Collections.sort(clusterList, new Comparator<Set<String>>() {
//                public int compare(Set<String> cluster1, Set<String> cluster2) {
//                    return cluster2.size() - cluster1.size();
//                }
//            });
//            int index = 0;
//            for (Set<String> cluster : clusterList) {
//                System.out.println("cluster " + index + ": " + cluster.size());
//                index ++;
//            }
//            // Output the largest container
//            StringBuilder builder = new StringBuilder();
//            for (int i = 0; i < clusterList.size(); i++) {
//                Set<String> cluster = clusterList.get(i);
//                if (cluster.size() < 10)
//                    continue; // Cutoff size: 10
//                builder.append(pathway).append(" ").append(i).append("\tNA");
//                for (String gene : cluster)
//                    builder.append("\t").append(gene);
//                fu.printLine(builder.toString());
//                builder.setLength(0);
//            }
//            BetweennessCentrality<String, String> centrality = new BetweennessCentrality<String, String>(graph, 
//                    true, true);
//            centrality.evaluate();
//            List<Ranking<?>> rankings = centrality.getRankings();
//            int index = 0;
//            double normalizer = (genes.size() - 1) * (genes.size() - 2) / 2.0;
//            for (Ranking<?> ranking : rankings) {
//                String object = (String) ranking.getRanked();
//                if (genes.contains(object)) {
//                    System.out.println(ranking.getRanked() + "\t" + ranking.rankScore + "\t" + ranking.rankScore / normalizer);
//                    index ++;
//                    if (index == 10)
//                        break;
//                }
//            }
//            c ++;
//            //break;
//        }
//        fu.close();
    }
    
    @Test
    public void checkPathwayCentrality() throws Exception {
//        Map<String, Set<String>> pathwayToGenes = new TopicAnalyzer().getTopicToNamesMap();
//        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
//        Set<String> allGenes = InteractionUtilities.grepIDsFromInteractions(fis);
//        int c = 1;
//        String outFileName = R3Constants.RESULT_DIR + "Pathways_Reactome_28_plus_i_betweeness_centrality.gmt";
//        fu.setOutput(outFileName);
//        for (String pathway : pathwayToGenes.keySet()) {
//            Set<String> genes = pathwayToGenes.get(pathway);
//            genes.retainAll(allGenes);
//            //fu.printLine("#" + pathway + "\t" + genes.size());
//            // Convert to LinLogLayout data structure
//            Set<String> fisInGenes = InteractionUtilities.getFIs(genes, fis);
//            Graph<String, String> graph = gu.createJungGraph(genes, fisInGenes);
//            BetweennessCentrality<String, String> centrality = new BetweennessCentrality<String, String>(graph,true, true);
//            centrality.evaluate();
//            List<Ranking<?>> rankings = centrality.getRankings();
//            double normalizer = (genes.size() - 1) * (genes.size() - 2) / 2.0;
//            StringBuilder builder = new StringBuilder();
//            builder.append(pathway).append("\tNA");
//            for (Ranking<?> ranking : rankings) {
//                String object = (String) ranking.getRanked();
//                if (genes.contains(object)) {
//                    //fu.printLine(ranking.getRanked() + "\t" + ranking.rankScore + "\t" + ranking.rankScore / normalizer);
//                    if (ranking.rankScore / normalizer >= 0.0001) {
//                        builder.append("\t").append(ranking.getRanked());
//                    }
//                }
//            }
//            fu.printLine(builder.toString());
//            c ++;
//            //break;
//        }
        fu.close();
    }
    
    /**
     * This method is used to calculate modularity based on paper: Newman, MEG & Girvan, M.
     * Physical Review E 69, 026113 (2004). Finding and evaluating community structure in networks.
     * In this paper, the modularity is defined as: Q = Sigma (e - a*a), which can be converted as
     * Q = Sigma (Eii / Et - (Ei / Et) * (Ei / Et)). For a much cleared definition, see Newman, MEJ 
     * PNAS 103: 8577-8582 (2006).
     * @param clusters
     * @param fis
     * @return
     */
    public double calculateModularity(Collection<Set<String>> clusters,
                                       Set<String> fis) {
        return new NetworkModularityCalculator().calculateModularity(clusters, fis);
    }
    
    /**
     * This method is used to generate a file that maps samples to a lists of booleans.
     * @param sampleToGenes
     * @param selectedGenes
     * @param outFileName
     * @throws IOException
     */
    public void generateSampleToSelectedGenesMatrix(Map<String, Set<String>> sampleToGenes,
                                                   Set<String> selectedGenes,
                                                   String outFileName) throws IOException {
        // Make selected genes as list for easy viewing
        List<String> geneList = new ArrayList<String>(selectedGenes);
        Collections.sort(geneList);
        Map<String, List<Boolean>> sampleToVector = new HashMap<String, List<Boolean>>();
        for (String sample : sampleToGenes.keySet()) {
            Set<String> genes = sampleToGenes.get(sample);
            List<Boolean> values = new ArrayList<Boolean>();
            for (String gene : geneList) {
                values.add(genes.contains(gene));
            }
            sampleToVector.put(sample, values);
        }
        outputSampleToBooleanList(outFileName, 
                                  sampleToVector,
                                  geneList);
    }
    
    public void generateSampleToModulesMatrix(Map<String, Set<String>> sampleToGenes,
                                              String clusterFileName,
                                              String outFileName,
                                              int sizeCutoff) throws Exception {
        Map<String, List<Boolean>> sampleToVector = generateSampleToNetworkModuleVectors(sampleToGenes, 
                                                                                         clusterFileName,
                                                                                         sizeCutoff);
        outputSampleToBooleanList(outFileName, 
                                  sampleToVector,
                                  null);
    }
    
    private void outputSampleToBooleanList(String outFileName,
                                           Map<String, List<Boolean>> sampleToVector,
                                           List<String> headNames)
            throws IOException {
        // Peel the size of list
        List<Boolean> list = sampleToVector.values().iterator().next();
        fu.setOutput(outFileName);
        StringBuilder builder = new StringBuilder();
        // Total cluster
        int totalCluster = list.size();
        builder.append("Sample");
        if (headNames == null) {
            for (int i = 0; i < totalCluster; i++)
                builder.append("\tModule").append(i);
        }
        else {
            for (String name : headNames)
                builder.append("\t").append(name);
        }
        fu.printLine(builder.toString());
        for (String sample : sampleToVector.keySet()) {
            builder.setLength(0);
            builder.append(sample);
            List<Boolean> vector = sampleToVector.get(sample);
            for (Boolean b : vector) {
                builder.append("\t");
                if (b)
                    builder.append(1);
                else
                    builder.append(0);
            }
            fu.printLine(builder.toString());
        }
        fu.close();
    }
    
    public void generateSampleToModulesMatrixInDegrees(Map<String, Set<String>> sampleToGenes,
                                                       String clusterFileName,
                                                       String outFileName,
                                                       int sizeCutoff) throws Exception {
        Map<String, List<Double>> sampleToVector = generateSampleToNetworkModuleVectorsInDegrees(sampleToGenes,
                                                                                                  clusterFileName,
                                                                                                  sizeCutoff);
        generateSampleToGeneSetDoubleMatrix(outFileName, 
                                            null,
                                            sampleToVector);
    }
    
    public void generateSampleToPathwayMatrixInDegrees(Map<String, Set<String>> sampleToGenes,
                                                       String outFileName,
                                                       String dbSymbol,
                                                       int sizeCutoff) throws Exception {
//        TopicAnalyzer topicAnalyzer = new TopicAnalyzer();
//        Map<String, Set<String>> pathwayToGenes = topicAnalyzer.getTopicToNamesMap();
//        Map<String, Set<String>> pathwayToFIs = topicAnalyzer.getTopicToFIInNameMap();
//        // Did a pre-filter
//        for (Iterator<String> it = pathwayToGenes.keySet().iterator(); it.hasNext();) {
//            String pathway = it.next();
//            Set<String> genes = pathwayToGenes.get(pathway);
//            if (genes.size() < sizeCutoff) {
//                it.remove();
//                continue;
//            }
//            if (dbSymbol != null && !pathway.endsWith(dbSymbol)) // Use KEGG only
//                it.remove();
//        }
//        System.out.println("Total pathways used: " + pathwayToGenes.size());
//        // Sort based on names
//        List<String> pathways = new ArrayList<String>(pathwayToGenes.keySet());
//        Collections.sort(pathways);
//        List<Set<String>> geneSets = new ArrayList<Set<String>>();
//        List<Set<String>> fisInGeneSets = new ArrayList<Set<String>>(); 
//        for (Iterator<String> it = pathways.iterator(); it.hasNext();) {
//            String pathway = it.next();
//            Set<String> geneset = pathwayToGenes.get(pathway);
//            Set<String> fis = pathwayToFIs.get(pathway);
//            if (geneset == null || fis == null) {
//                it.remove();  // Just remove it
//                continue;
//            }
//            geneSets.add(geneset);
//            fisInGeneSets.add(fis);
//        }
//        Map<String, List<Double>> sampleToVectors = generateSampleToGeneSetVectorsInDegrees(sampleToGenes, 
//                                                                                            sizeCutoff,
//                                                                                            geneSets,
//                                                                                            fisInGeneSets);
//        generateSampleToGeneSetDoubleMatrix(outFileName, 
//                                            pathways,
//                                            sampleToVectors);
    }

    private void generateSampleToGeneSetDoubleMatrix(String outFileName,
                                                     List<String> headerNames,
                                                     Map<String, List<Double>> sampleToVector)
            throws IOException {
        // Peel the size of list
        List<Double> list = sampleToVector.values().iterator().next();
        fu.setOutput(outFileName);
        StringBuilder builder = new StringBuilder();
        // Total cluster
        int totalCluster = list.size();
        builder.append("Sample");
        if (headerNames == null)
            for (int i = 0; i < totalCluster; i++)
                builder.append("\tModule").append(i);
        else 
            for (String header : headerNames)
                builder.append("\t").append(header);
        fu.printLine(builder.toString());
        for (String sample : sampleToVector.keySet()) {
            builder.setLength(0);
            builder.append(sample);
            List<Double> vector = sampleToVector.get(sample);
            for (Double b : vector) {
                builder.append("\t").append(b);
            }
            fu.printLine(builder.toString());
        }
        fu.close();
    }
    
    public void generateSampleToModulesMatrixInPValue(Map<String, Set<String>> sampleToGenes,
                                                      String clusterFileName,
                                                      String outFileName,
                                                      int sizeCutoff,
                                                      int totalGenes) throws Exception {
        Map<String, List<Double>> sampleToVector = generateSampleToNetworkModuleVectorsInPValue(sampleToGenes, 
                                                                                                 clusterFileName,
                                                                                                 sizeCutoff,
                                                                                                 totalGenes);
        generateSampleToGeneSetDoubleMatrix(outFileName,
                                            null,
                                            sampleToVector);
    }
    
    private Map<String, List<Boolean>> generateSampleToNetworkModuleVectors(Map<String, Set<String>> sampleToGenes,
                                                                            String clusterFileName,
                                                                            int sizeCutoff) throws IOException {
        List<Set<String>> clusters = loadNetworkClusters(clusterFileName);
        // Want to create a sample to a boolean vector map
        Map<String, List<Boolean>> sampleToVector = new HashMap<String, List<Boolean>>();
        for (String sample : sampleToGenes.keySet()) {
            Set<String> alteredGenes = sampleToGenes.get(sample);
            List<Boolean> vector = new ArrayList<Boolean>();
            for (Set<String> cluster : clusters) {
                if (cluster.size() < sizeCutoff)
                    continue;
                Set<String> shared = InteractionUtilities.getShared(cluster, alteredGenes);
                if (shared.size() > 0)
                    vector.add(Boolean.TRUE);
                else
                    vector.add(Boolean.FALSE);
                //vector.add(shared.size());
            }
            //System.out.println(sample + ": " + StringUtils.join(", ", vector));
            sampleToVector.put(sample, vector);
        }
        return sampleToVector;
    }
    
    private Map<String, List<Double>> generateSampleToNetworkModuleVectorsInDegrees(Map<String, Set<String>> sampleToGenes,
                                                                                     String clusterFileName,
                                                                                     int sizeCutoff) throws IOException {
        List<Set<String>> clusters = loadNetworkClusters(clusterFileName);
        return generateSampleToGeneSetVectorsInDegrees(sampleToGenes,
                                                       sizeCutoff, 
                                                       clusters,
                                                       null);
    }
    
    private Map<String, List<Double>> generateSampleToGeneSetVectorsInDegrees(Map<String, Set<String>> sampleToGenes,
                                                                              int sizeCutoff,
                                                                              List<Set<String>> clusters,
                                                                              List<Set<String>> fisInClusters) throws IOException {
        // Filter out clusters based on size
        for (java.util.Iterator<Set<String>> it = clusters.iterator(); it.hasNext();) {
            Set<String> cluster = it.next();
            if (cluster.size() < sizeCutoff)
                it.remove();
        }
        // Pre-calculate degrees for genes in clusters
        List<Map<String, Integer>> geneToDegreesInClusters = new ArrayList<Map<String,Integer>>();
        if (fisInClusters == null) {
            Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
            //Set<String> fis = fu.loadInteractions(R3Constants.RESULT_DIR + "FIsInGene_Pathway_041709.txt");
            for (Set<String> cluster : clusters) {
                Set<String> fisInCluster = InteractionUtilities.getFIs(cluster, fis);
                Map<String, Integer> geneToDegree = InteractionUtilities.generateProteinToDegree(fisInCluster);
                geneToDegreesInClusters.add(geneToDegree);
            }
        }
        else {
            for (int i = 0; i < clusters.size(); i++) {
                Set<String> fisInCluster = fisInClusters.get(i);
                Map<String, Integer> geneToDegree = InteractionUtilities.generateProteinToDegree(fisInCluster);
                geneToDegreesInClusters.add(geneToDegree);
            }
        }
        // Used to normalization
        List<Integer> totalDegrees = new ArrayList<Integer>();
        for (int i = 0; i < clusters.size(); i++) {
            int total = 0;
            Set<String> cluster = clusters.get(i);
            Map<String, Integer> geneToDegrees = geneToDegreesInClusters.get(i);
            for (String gene : geneToDegrees.keySet()) {
                total += geneToDegrees.get(gene);
            }
            totalDegrees.add(total);
        }
        // Want to create a sample to a boolean vector map
        Map<String, List<Double>> sampleToVector = new HashMap<String, List<Double>>();
        for (String sample : sampleToGenes.keySet()) {
            Set<String> alteredGenes = sampleToGenes.get(sample);
            List<Double> vector = new ArrayList<Double>();
            for (int i = 0; i < clusters.size(); i++) {
                Set<String> cluster = clusters.get(i);
                Set<String> shared = InteractionUtilities.getShared(cluster, alteredGenes);
                if (shared.size() == 0)
                    vector.add(0.0);
                else {
                    Map<String, Integer> geneToDegree = geneToDegreesInClusters.get(i);
                    int totalDegree = totalDegrees.get(i);
                    // Need to get the degree
                    int total = 0;
                    for (String gene : shared) {
                        Integer degree = geneToDegree.get(gene);
                        if (degree != null)
                            total += degree;
                    }
                    //vector.add((double)total);
                    vector.add((double) total / totalDegree);
                }
            }
            //System.out.println(sample + ": " + StringUtils.join(", ", vector));
            sampleToVector.put(sample, vector);
        }
        return sampleToVector;
    }
    
    private Map<String, List<Double>> generateSampleToNetworkModuleVectorsInPValue(Map<String, Set<String>> sampleToGenes,
                                                                                   String clusterFileName,
                                                                                   int sizeCutoff,
                                                                                   int totalGense) throws Exception {
        List<Set<String>> clusters = loadNetworkClusters(clusterFileName);
        // Want to create a sample to a boolean vector map
        Map<String, List<Double>> sampleToVector = new HashMap<String, List<Double>>();
        for (String sample : sampleToGenes.keySet()) {
            Set<String> alteredGenes = sampleToGenes.get(sample);
            List<Double> vector = new ArrayList<Double>();
            for (Set<String> cluster : clusters) {
                if (cluster.size() < sizeCutoff)
                    continue;
                Set<String> shared = InteractionUtilities.getShared(cluster, alteredGenes);
                Double pvalue = MathUtilities.calculateHypergeometricPValue(totalGense,
                                                                            alteredGenes.size(),
                                                                            cluster.size(), 
                                                                            shared.size());
                vector.add(pvalue);
            }
            sampleToVector.put(sample, vector);
        }
        return sampleToVector;
    }
    
    /**
     * Hierarchical cluster a set of samples based on network modules.
     * @throws Exception
     */
    public List<HierarchicalClusterNode> hierarchicalClusterSamplesOnNetworkModules(Map<String, Set<String>> sampleToGenes,
                                                                        String clusterFileName) throws Exception {
        Map<String, List<Boolean>> sampleToVector = generateSampleToNetworkModuleVectors(sampleToGenes,
                                                                                         clusterFileName,
                                                                                         0);
        // Create pair-wise distances among samples
        Map<String, Double> samplePairToDist = new HashMap<String, Double>();
        List<String> sampleList = new ArrayList<String>(sampleToGenes.keySet());
        for (int i = 0; i < sampleList.size() - 1; i++) {
            String sample1 = sampleList.get(i);
            List<Boolean> vector1 = sampleToVector.get(sample1);
            for (int j = i + 1; j < sampleList.size(); j++) {
                String sample2 = sampleList.get(j);
                List<Boolean> vector2 = sampleToVector.get(sample2);
                double dist = MathUtilities.calculateHammingDistance(vector1, vector2);
                //double dist = MathUtilities.calculateNetworkDistance(vector1,
                //                                                     vector2, 
                //                                                     clusters);
                samplePairToDist.put(sample1 + "\t" + sample2,
                                     dist);
            }
        }
        // Do a hierarchical clustering
        CancerResequenceDataSetAnalyzer analyzer = new CancerResequenceDataSetAnalyzer();
        List<HierarchicalClusterNode> sampleClusters = analyzer.hierarchicalCluster(sampleList, 
                                                                        samplePairToDist,
                                                                        false);
        return sampleClusters;
    }
    
}
